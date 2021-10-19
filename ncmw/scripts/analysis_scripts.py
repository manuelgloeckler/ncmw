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
import glob
import json

import matplotlib.pyplot as plt

file_dir = os.path.dirname(os.path.dirname(__file__))
sys.path.append(file_dir)


from ncmw.utils.utils_io import (
    get_models,
    get_result_path,
    SEPERATOR,
)


from ncmw.analysis import (
    compute_fvas,
    jaccard_similarity_matrices,
    compute_COMPM,
    table_ex_transport,
    sekretion_uptake_fva,
    sekretion_uptake_fba,
    compute_uptake_sekretion_table,
)
from ncmw.visualization import (
    plot_full_fva,
    plot_medium_fva_range,
    jacard_index_similarity_heatmap,
    plot_scaled_medium_growth,
    uptake_sekretion_venn_diagrams,
    plot_growth_sensitivity,
)


@hydra.main(config_path="../../data/hydra", config_name="config.yaml")
def run_analysis_hydra(cfg: DictConfig) -> None:
    run_analysis(cfg)


def run_analysis(cfg: DictConfig) -> None:
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
    PATH = PATH_res + SEPERATOR + "analysis"
    log.info(f"Working directory: {PATH}")

    try:
        if not os.path.exists(PATH_res):
            os.mkdir(PATH_res)
        if not os.path.exists(PATH):
            os.mkdir(PATH)
        if not os.path.exists(PATH + SEPERATOR + "medium"):
            os.mkdir(PATH + SEPERATOR + "medium")
        if not os.path.exists(PATH + SEPERATOR + "growth"):
            os.mkdir(PATH + SEPERATOR + "growth")
        if not os.path.exists(PATH + SEPERATOR + "flux_analysis"):
            os.mkdir(PATH + SEPERATOR + "flux_analysis")
        if not os.path.exists(PATH + SEPERATOR + "similarity"):
            os.mkdir(PATH + SEPERATOR + "similarity")
        if not os.path.exists(PATH + SEPERATOR + "sekretion_uptake"):
            os.mkdir(PATH + SEPERATOR + "sekretion_uptake")
    except:
        pass

    if cfg.analysis.require_setup:
        model_PATH = PATH_res + SEPERATOR + "setup" + SEPERATOR + "snm3_models"
        assert os.path.exists(model_PATH), "We require the setup to run first."
        log.info("Loading models")
        models = get_models(
            "snm3_models", prefix=PATH_res + SEPERATOR + "setup" + SEPERATOR
        )
    else:
        log.info("Loading original models, WITHOUT setup")
        models = get_models(cfg.setup.models)

    # Check which results are already there
    models = np.array(models)
    # models_mask = np.zeros(len(models), dtype=np.int32)
    # for filepath in glob.iglob(PATH + SEPERATOR + "flux_analysis" + SEPERATOR + "*"):
    #     models_already_done = []
    #     for i in range(len(models)):
    #         if models[i].id in filepath:
    #             models_already_done.append(i)
    #     for m in list(set(models_already_done)):
    #         models_mask[m] = 1

    log.info("Generating fva results")
    # dfs = compute_fvas(models[models_mask], cfg.analysis.fva_fraction)
    dfs = compute_fvas(models, cfg.analysis.fva_fraction)
    for i, (model, df) in enumerate(zip(models, dfs)):
        # if models_mask[i] == 0:
        sol = model.optimize()
        df["flux"] = sol.fluxes
        df.to_csv(
            PATH + SEPERATOR + "flux_analysis" + SEPERATOR + "fva_" + model.id + ".csv"
        )
        log.info(df)
    #     elif models_mask[i] == 1:
    #         df = pd.read_csv(
    #             PATH
    #             + SEPERATOR
    #             + "flux_analysis"
    #             + SEPERATOR
    #             + "fva_"
    #             + model.id
    #             + ".csv",
    #             index_col=0,
    #         )
    #         log.info(df)
    #         dfs.insert(i, df)

    log.info("Computing COMPM media for modles")
    # mediums = compute_COMPM(
    #     [models[i] for i in range(len(models)) if models_mask[i] == 1],
    #     [dfs[i] for i in range(len(dfs)) if models_mask[i] == 1],
    # )
    mediums = compute_COMPM(models, dfs)
    for model, medium in zip(models, mediums):
        with open(
            PATH + SEPERATOR + "medium" + SEPERATOR + "COMPM_" + model.id + ".json", "w"
        ) as f:
            json.dump(medium, f)

    log.info("Plotting fva results")
    for model, df in zip(models, dfs):
        plot_full_fva(
            df,
            PATH
            + SEPERATOR
            + "flux_analysis"
            + SEPERATOR
            + "full_fva_plot"
            + model.id
            + ".pdf",
        )
        plot_medium_fva_range(
            model,
            PATH
            + SEPERATOR
            + "flux_analysis"
            + SEPERATOR
            + "medium_fva_plot"
            + model.id
            + ".pdf",
            cfg.analysis.fva_fraction,
        )

    log.info("Computing Uptake Sekretions + Transports")
    uptakes = []
    sekretions = []
    for i, model in enumerate(models):
        # Transport reactions
        if cfg.analysis.check_transport:
            transport_check = table_ex_transport(model)
            transport_check.to_csv(
                PATH
                + SEPERATOR
                + "sekretion_uptake"
                + SEPERATOR
                + model.id
                + "_transport_summary.csv"
            )

        # Sekretion uptakes
        if cfg.analysis.sekretion_uptake == "fva":
            uptake, sekretion = sekretion_uptake_fva(dfs[i])
        else:
            uptake, sekretion = sekretion_uptake_fba(model)
        log.info(
            f"Model: {model.id} has {len(uptake)} uptakes and {len(sekretions)} sekretions."
        )
        uptakes.append(uptake)
        sekretions.append(sekretion)

        # Plots for uptake sensitivity!
        fig = plot_growth_sensitivity(model, uptake)
        fig.savefig(
            PATH
            + SEPERATOR
            + "growth"
            + SEPERATOR
            + f"{model.id}_uptake_growth_sensitivity.pdf"
        )

    for i in range(len(models)):
        for j in range(i + 1, len(models)):
            log.info(f"Comparing flux interchange between {models[i]} and {models[j]}")
            uptake_sekretion_table = compute_uptake_sekretion_table(
                models[i].id,
                models[j].id,
                deepcopy(uptakes[i]),
                deepcopy(uptakes[j]),
                deepcopy(sekretions[i]),
                deepcopy(sekretions[j]),
            )
            uptake_sekretion_table.to_csv(
                PATH
                + SEPERATOR
                + "sekretion_uptake"
                + SEPERATOR
                + f"{models[i].id}_{models[j].id}_uptake_sekretion_summary.csv"
            )

    # if models_mask.sum() < len(models):
    if True:
        log.info(f"Plotting venn diagrams for uptake sekretion overlaps")

        fig = uptake_sekretion_venn_diagrams(
            models,
            uptakes,
            sekretions,
            names=cfg.visualization.names,
            cmap=cfg.visualization.cmap,
        )
        fig.savefig(
            PATH
            + SEPERATOR
            + "sekretion_uptake"
            + SEPERATOR
            + "uptake_sekretion_overlap_plot.pdf"
        )

    # if models_mask.sum() < len(models):
    if True:
        log.info("Computing Jacard Similarities")
        df_met, df_rec, df_ro = jaccard_similarity_matrices(models)
        df_met.to_csv(
            PATH
            + SEPERATOR
            + "similarity"
            + SEPERATOR
            + "jacard_similarities_metabolies.csv"
        )
        df_rec.to_csv(
            PATH
            + SEPERATOR
            + "similarity"
            + SEPERATOR
            + "jacard_similarities_reactions.csv"
        )
        df_ro.to_csv(
            PATH
            + SEPERATOR
            + "similarity"
            + SEPERATOR
            + "jacard_similarities_exchanges.csv"
        )
        fig = jacard_index_similarity_heatmap(
            df_met, df_rec, df_ro, names=cfg.visualization.names
        )
        fig.savefig(
            PATH + SEPERATOR + "similarity" + SEPERATOR + "similarity_summary.pdf"
        )

    # if models_mask.sum() < len(models):
    if True:
        log.info("Computing scaled medium growth plot")
        kwargs = cfg.visualization.scaled_medium_growth_plot
        fig = plot_scaled_medium_growth(
            models,
            kwargs.min_scale,
            kwargs.max_scale,
            kwargs.evaluations,
            names=cfg.visualization.names,
            cmap=cfg.visualization.cmap,
        )
        fig.savefig(
            PATH + SEPERATOR + "growth" + SEPERATOR + "scaled_medium_growth_plot.pdf"
        )

    end_time = time.time()
    log.info(f"Job finished in {end_time-start_time} seconds")
