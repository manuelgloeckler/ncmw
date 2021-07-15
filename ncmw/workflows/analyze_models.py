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


from utils import get_models, get_result_path, SEPERATOR, save_model, get_model_paths
from analysis import compute_fvas, jaccard_similarity_matrices, compute_COMPM
from visualization import (
    plot_full_fva,
    plot_medium_fva_range,
    jacard_index_similarity_heatmap,
    plot_scaled_medium_growth,
)

import cobra


@hydra.main(config_path="../../data/hydra", config_name="config.yaml")
def run_analysis(cfg: DictConfig) -> None:
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
    PATH = PATH_res + SEPERATOR + "analysis"
    log.info(f"Working directory: {PATH}")

    try:
        os.mkdir(PATH)
        os.mkdir(PATH)
        os.mkdir(PATH + SEPERATOR + "medium")
        os.mkdir(PATH + SEPERATOR + "growth")
        os.mkdir(PATH + SEPERATOR + "flux_analysis")
        os.mkdir(PATH + SEPERATOR + "similarity")
    except:
        pass

    log.info("Loading models")
    models = get_models(
        "snm3_models", prefix=PATH_res + SEPERATOR + "setup" + SEPERATOR
    )

    log.info("Generating fva results")
    dfs = compute_fvas(models, cfg.analysis.fva_fraction)
    for model, df in zip(models, dfs):
        sol = model.optimize()
        df["flux"] = sol.fluxes
        df.to_csv(
            PATH + SEPERATOR + "flux_analysis" + SEPERATOR + "fva_" + model.id + ".csv"
        )
        log.info(df)

    log.info("Computing COMPM media for modles")
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
    fig = jacard_index_similarity_heatmap(df_met, df_rec, df_ro)
    fig.savefig(PATH + SEPERATOR + "similarity" + SEPERATOR + "similarity_summary.pdf")

    log.info("Computing scaled medium growth plot")
    kwargs = cfg.visualization.scaled_medium_growth_plot
    fig = plot_scaled_medium_growth(
        models, kwargs.min_scale, kwargs.max_scale, kwargs.evaluations
    )
    fig.savefig(
        PATH + SEPERATOR + "growth" + SEPERATOR + "scaled_medium_growth_plot.pdf"
    )

    log.info("Computing Uptake Sekretions + Transports")


if __name__ == "__main__":
    run_analysis()
