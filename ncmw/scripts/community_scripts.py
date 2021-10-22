import hydra
from omegaconf import DictConfig, OmegaConf

import logging
import socket

import time

import random
import numpy as np
import pandas as pd
from copy import deepcopy

import sys, os
import json

file_dir = os.path.dirname(os.path.dirname(__file__))
sys.path.append(file_dir)


from ncmw.utils.utils_io import (
    get_models,
    get_result_path,
    SEPERATOR,
    save_model,
    get_model_paths,
)

from ncmw.visualization import (
    plot_pairwise_growth_relation_per_weight,
    plot_community_interaction,
    plot_posterior_samples_for_observations,
    plot_community_uptake_graph,
    plot_species_interaction,
    plot_community_summary,
)

from ncmw.utils.utils_io import (
    get_models,
    get_result_path,
    SEPERATOR,
    save_model,
    get_model_paths,
    get_mediums,
)
from ncmw.community import (
    BagOfReactionsModel,
    ShuttleCommunityModel,
    compute_fair_weights,
    compute_community_summary,
    compute_dominant_weights,
)

from ncmw.utils import (
    get_models,
    get_result_path,
    SEPERATOR,
)


@hydra.main(config_path="../../data/hydra", config_name="config.yaml")
def run_community_hydra(cfg: DictConfig) -> None:
    run_community(cfg)


def run_community(cfg: DictConfig) -> None:
    log = logging.getLogger(__name__)
    log.setLevel(logging.INFO)
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
        if not os.path.exists(PATH_res):
            os.mkdir(PATH_res)
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

        path = (
            PATH
            + SEPERATOR
            + "community_models"
            + SEPERATOR
            + "BagOfReactionsModel.pkl"
        )
        if os.path.exists(path):
            log.info("Loading BagOfReactions Community model")
            m = BagOfReactionsModel.load(path)
            ids = [mod.id for mod in m.models]
            correct = True
            for model in models:
                correct = correct and (model.id in ids)
            if correct:
                community_models += [m]
            else:
                community_models += [BagOfReactionsModel(models)]
        else:
            log.info("Building BagOfReactions Community model")
            community_models += [BagOfReactionsModel(models)]
    if cfg.community.shuttle_reactions_model:
        path = (
            PATH
            + SEPERATOR
            + "community_models"
            + SEPERATOR
            + "ShuttleCommunityModel.pkl"
        )
        if os.path.exists(path):
            log.info("Loading shuttle reactions model")
            m = ShuttleCommunityModel.load(path)
            ids = [mod.id for mod in m.models]
            correct = True
            for model in models:
                correct = correct and (model.id in ids)
            if correct:
                community_models += [m]
            else:
                community_models += [ShuttleCommunityModel(models)]
        else:
            log.info("Building shuttle reactions model")
            community_models += [ShuttleCommunityModel(models)]

    log.info(f"Saving models: {cfg.community.save_models}")
    if cfg.community.save_models:
        for m in community_models:
            path = PATH + SEPERATOR + "community_models" + SEPERATOR + type(m).__name__
            m.save(path)
            log.info(f"Saving community models: {path}")

    log.info(f"Community summary on default medium: ")
    for m in community_models:
        if m._type == "compartmentalized":
            fig = plot_community_summary(m, cfg.visualization.names)
            fig.savefig(
                PATH + SEPERATOR + "experiments" + SEPERATOR + f"community_summary.pdf"
            )

    # Medias
    medium_prefix = PATH_res + SEPERATOR + "analysis" + SEPERATOR
    mediums = get_mediums("medium", medium_prefix)

    if cfg.community.pairwise_growth:
        log.info("Compute pairwise growth relationships for different weights.")
        for m in community_models:
            fig = plot_pairwise_growth_relation_per_weight(m, cfg.visualization.names)
            fig.savefig(
                PATH
                + SEPERATOR
                + "experiments"
                + SEPERATOR
                + f"pairwise_growth_relationship_{type(m).__name__}.pdf"
            )

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
    log.info(f"We consider following weights: {weights}")

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

    log.info(f"Compute community visualization")
    for m in community_models:
        if isinstance(m, ShuttleCommunityModel):
            m.weights = np.ones(len(m.weights))
            medium_compm = COMPMS[0]
            m.medium = medium_compm
            df = m.summary()
            fig2 = plot_community_interaction(m, df, names=cfg.visualization.names)
            fig2.savefig(
                PATH
                + SEPERATOR
                + "experiments"
                + SEPERATOR
                + type(m).__name__
                + "_compm_community_interaction.pdf"
            )
            fig2 = plot_community_uptake_graph(m, df, names=cfg.visualization.names)
            fig2.savefig(
                PATH
                + SEPERATOR
                + "experiments"
                + SEPERATOR
                + type(m).__name__
                + "_compm_uptake_flux.png"
            )

            fig2 = plot_species_interaction(m, df, names=cfg.visualization.names)
            fig2.savefig(
                PATH
                + SEPERATOR
                + "experiments"
                + SEPERATOR
                + type(m).__name__
                + "_compm_species_interaction.pdf"
            )

            medium_coopm = m.computeCOOPM(m.slim_optimize(), enforce_survival=True)
            m.medium = medium_coopm
            df = m.summary()
            fig1 = plot_community_interaction(m, df, names=cfg.visualization.names)
            fig1.savefig(
                PATH
                + SEPERATOR
                + "experiments"
                + SEPERATOR
                + type(m).__name__
                + "_coopm_community_interaction.pdf"
            )
            fig1 = plot_community_uptake_graph(m, df, names=cfg.visualization.names)
            fig1.savefig(
                PATH
                + SEPERATOR
                + "experiments"
                + SEPERATOR
                + type(m).__name__
                + "_coopm_uptake_flux.png",
                dpi=1000,
            )

            fig1 = plot_species_interaction(m, df, names=cfg.visualization.names)
            fig1.savefig(
                PATH
                + SEPERATOR
                + "experiments"
                + SEPERATOR
                + type(m).__name__
                + "_coopm_species_interaction.pdf"
            )

    # COOPM
    log.info(f"Compute community COOPMs enforced survival")
    COOPMS_enforced = dict()
    for m, compm, weight in zip(community_models, COMPMS, weights):
        for w in weight:
            m.medium = compm
            m.weights = w
            MBR = m.slim_optimize()
            if MBR is None:
                MBR = 0
            try:
                coopm = m.computeCOOPM(MBR)
            except:
                coopm = dict()
            if tuple(w) in COOPMS_enforced:
                COOPMS_enforced[tuple(w)].append(coopm)
            else:
                COOPMS_enforced[tuple(w)] = [coopm]
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

    COOPMS_not_enforced = dict()

    log.info(f"Compute community COOPMs enforced survival")
    for m, compm, weight in zip(community_models, COMPMS, weights):
        for w in weight:
            log.info(f"Computing coopm for community with weights {w}")
            m.medium = compm
            m.weights = w
            MBR = m.slim_optimize()
            if MBR is None:
                MBR = 0
            try:
                coopm = m.computeCOOPM(MBR, enforce_survival=False)
            except:
                coopm = dict()
            if tuple(w) in COOPMS_not_enforced:
                COOPMS_not_enforced[tuple(w)].append(coopm)
            else:
                COOPMS_not_enforced[tuple(w)] = [coopm]
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
            coopms = COOPMS_enforced[tuple(weight)]
            if len(coopms) > 1:
                coopm = coopms[i]
            else:
                coopm = coopms[0]
            m.medium = coopm
            df_coopm = compute_community_summary(m, [weight])
            df_coopm["Medium"] = "COOPM_enforced_survival"
            df = df.append(df_coopm)

            coopms = COOPMS_not_enforced[tuple(weight)]
            if len(coopms) > 1:
                coopm = coopms[i]
            else:
                coopm = coopms[0]
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

    log.info("Compute weight posteriors")
    for m in community_models:
        m.weights = np.ones(len(m.weights))
        medium_coopm = COOPMS_enforced[tuple(m.weights)]
        maximal_growths = [m.single_optimize(i) for i in range(len(m.models))]
        balanced_observation = np.ones(len(m.models)) * np.min(maximal_growths)
        observations = [balanced_observation]
        for i in range(len(m.models)):
            dominant_observation = np.ones(len(m.models)) * min(
                0.1, np.min(maximal_growths)
            )
            dominant_observation[i] = maximal_growths[i]
            observations.append(dominant_observation)

        figs = plot_posterior_samples_for_observations(m, observations)
        for i in range(len(figs)):
            figs[i].savefig(
                PATH
                + SEPERATOR
                + "experiments"
                + SEPERATOR
                + type(m).__name__
                + f"_coopm_weight_posterior_{i}.pdf"
            )

        medium_compm = COMPMS[0]
        m.medium = medium_compm
        maximal_growths = [m.single_optimize(i) for i in range(len(m.models))]
        balanced_observation = np.ones(len(m.models)) * np.min(maximal_growths)
        observations = [balanced_observation]
        for i in range(len(m.models)):
            dominant_observation = np.ones(len(m.models)) * min(
                0.1, np.min(maximal_growths)
            )
            dominant_observation[i] = maximal_growths[i]
            observations.append(dominant_observation)

        figs = plot_posterior_samples_for_observations(m, observations)
        for i in range(len(figs)):
            figs[i].savefig(
                PATH
                + SEPERATOR
                + "experiments"
                + SEPERATOR
                + type(m).__name__
                + f"_compm_weight_posterior_{i}.pdf"
            )

    end_time = time.time()
    runtime = end_time - start_time
    log.info(f"Finished Workflow in {runtime} seconds")
