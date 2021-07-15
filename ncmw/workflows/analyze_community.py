import hydra
from omegaconf import DictConfig, OmegaConf

import logging
import socket
import importlib

import time

import random
import numpy as np
import os
import pandas as pd
import json


from utils import (
    get_models,
    get_result_path,
    SEPERATOR,
    save_model,
    get_model_paths,
    get_mediums,
)
from community import (
    BagOfReactionsModel,
    ShuttleCommunityModel,
    compute_pairwise_growth_relation_per_weight,
    compute_fair_weights,
    compute_community_summary,
    compute_dominant_weights,
)

import cobra


@hydra.main(config_path="../../data/hydra", config_name="config.yaml")
def run_community(cfg: DictConfig) -> None:
    log = logging.getLogger(__name__)
    log.info(OmegaConf.to_yaml(cfg))
    log.info(f"Hostname: {socket.gethostname()}")

    seed = cfg.seed
    random.seed(seed)
    np.random.seed(seed)
    log.info(f"Random seed: {seed}")

    start_time = time.time()
    name = cfg.name

    PATH_res = get_result_path(name)
    PATH = PATH_res + SEPERATOR + "community"
    log.info(f"Working directory: {PATH}")
    try:
        if not os.path.exists(PATH):
            os.mkdir(PATH)
        if not os.path.exists(PATH + SEPERATOR + "community_models"):
            os.mkdir(PATH + SEPERATOR + "community_models")
        if not os.path.exists(PATH + SEPERATOR + "medium"):
            os.mkdir(PATH + SEPERATOR + "medium")
        if not os.path.exists(PATH + SEPERATOR + "experiments"):
            os.mkdir(PATH + SEPERATOR + "experiments")
        log.info("Created folder structure successfully")
    except:
        pass

    log.info("Loading models from setup")
    try:
        models = get_models(
            "snm3_models", prefix=PATH_res + SEPERATOR + "setup" + SEPERATOR
        )
        num_models = len(models)
    except:
        raise ValueError(
            "We require snm3 models for the next steps, please run the setup!"
        )

    log.info("Generating community models")
    community_models = []
    if cfg.community.bag_of_reactions_model:
        log.info("Building BagOfReactions Community model")
        community_models += [BagOfReactionsModel(models)]
    if cfg.community.shuttle_reactions_model:
        log.info("Building shuttle reactions model")
        community_models += [ShuttleCommunityModel(models)]

    log.info(f"Saving models: {cfg.community.save_models}")
    if cfg.community.save_models:
        for m in community_models:
            path = PATH + SEPERATOR + "community_models" + SEPERATOR + type(m).__name__
            m.save(path)
            log.info(f"Saving community models: {path}")

    # Medias
    medium_prefix = PATH_res + SEPERATOR + "analysis" + SEPERATOR
    mediums = get_mediums("medium", medium_prefix)

    # Weights to consider
    weights = []
    for m in community_models:
        ws = []
        if cfg.community.weights.one_to_one:
            ws.append(np.ones(num_models))
        if cfg.community.weights.fair:
            ws.append(compute_fair_weights(m))
        if cfg.community.weights.dominant_weights:
            ws.extend(compute_dominant_weights(m))
        if cfg.community.weights.custom != []:
            for weight in cfg.community.weights:
                assert (
                    len(weight) == num_models
                ), "The custom weights must have dimension equal to the number of models in the community!"
                ws.append(weight)
        weights.append(ws)

    log.info(f"Compute community COMPM")
    COMPMS = []
    SNM3s = []
    for m in community_models:
        COMPM = dict()
        for medium in mediums:
            if "COMPM" in medium:
                COMPM = {**COMPM, **mediums[medium]}
        COMPMS.append(COMPM)
        # COMPMS
        path = (
            PATH
            + SEPERATOR
            + "medium"
            + SEPERATOR
            + "COMPM_"
            + type(m).__name__
            + ".json"
        )
        with open(path, "w+") as f:
            json.dump(COMPM, f)
        log.info(f"Saving COMPMs: {path}")
        # Also save default medium
        path = (
            PATH
            + SEPERATOR
            + "medium"
            + SEPERATOR
            + "SNM3_"
            + type(m).__name__
            + ".json"
        )
        with open(path, "w+") as f:
            json.dump(m.medium, f)
        SNM3s.append(m.medium)
        log.info(f"Saving COMPMs: {path}")

    # COOPM
    log.info(f"Compute community COOPMs")
    COOPMS_enforced = dict(
        zip([tuple(w) for w in weights[0]], [[] for _ in weights[0]])
    )
    for m, compm, weight in zip(community_models, COMPMS, weights):
        for w in weight:
            m.medium = compm
            m.weights = w
            MBR = m.slim_optimize()
            coopm = m.computeCOOPM(MBR)
            COOPMS_enforced[tuple(w)].append(coopm)
            path = (
                PATH
                + SEPERATOR
                + "medium"
                + SEPERATOR
                + "COOPM_survival"
                + str([np.round(i, 2) for i in w])
                + type(m).__name__
                + ".json"
            )
            with open(path, "w+") as f:
                json.dump(coopm, f)

    COOPMS_not_enforced = dict(
        zip([tuple(w) for w in weights[0]], [[] for _ in weights[0]])
    )

    for m, compm, weight in zip(community_models, COMPMS, weights):
        for w in weight:
            m.medium = compm
            m.weights = w
            MBR = m.slim_optimize()
            coopm = m.computeCOOPM(MBR, enforce_survival=False)
            COOPMS_not_enforced[tuple(w)].append(coopm)
            path = (
                PATH
                + SEPERATOR
                + "medium"
                + SEPERATOR
                + "COOPM_"
                + str([np.round(i, 2) for i in w])
                + type(m).__name__
                + ".json"
            )
            with open(path, "w+") as f:
                json.dump(coopm, f)

    log.info(f"Compute community growth summary on snm3")
    for i, m in enumerate(community_models):
        df = pd.DataFrame()
        # SNM3
        m.medium = SNM3s[i]
        df_snm3 = compute_community_summary(m, weights[i])
        df_snm3["Medium"] = "SNM3"
        df = df.append(df_snm3)

        # COMPM
        m.medium = COMPMS[i]
        df_compm = compute_community_summary(m, weights[i])
        df_compm["Medium"] = "COMPM"
        df = df.append(df_compm)

        # COOPM enforced
        for weight in weights[i]:
            coopm = COOPMS_enforced[tuple(weight)][i]
            m.medium = coopm
            df_coopm = compute_community_summary(m, [weight])
            df_coopm["Medium"] = "COOPM_enforced_survival"
            df = df.append(df_coopm)

            coopm = COOPMS_not_enforced[tuple(weight)][i]
            m.medium = coopm
            df_coopm = compute_community_summary(m, [weight])
            df_coopm["Medium"] = "COOPM"
            df = df.append(df_coopm)

        path = (
            PATH
            + SEPERATOR
            + "experiments"
            + SEPERATOR
            + type(m).__name__
            + "_growth_summary.csv"
        )
        df.to_csv(path)

    log.info("Compute pairwise growth relationships for different weights.")
    for m in community_models:
        print(compute_pairwise_growth_relation_per_weight(m, 0, 1))


if __name__ == "__main__":
    run_community()
