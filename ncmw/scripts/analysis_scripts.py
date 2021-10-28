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
import pickle
import json

import matplotlib.pyplot as plt

from ncmw.visualization.similarity_visualization import expected_community_interaction

file_dir = os.path.dirname(os.path.dirname(__file__))
sys.path.append(file_dir)


from ncmw.utils.utils_io import (
    get_models,
    get_result_path,
    SEPERATOR,
    check_for_substring_in_folder,
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
    cross_feed_venn_diagrams,
)


@hydra.main(config_path="../../data/hydra", config_name="config.yaml")
def run_analysis_hydra(cfg: DictConfig) -> None:
    run_analysis(cfg)


def run_analysis(cfg: DictConfig) -> None:
    log = logging.getLogger(__name__)
    log.setLevel(logging.INFO)
    global_loger = logging.getLogger()
    global_loger.setLevel(logging.WARNING)
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
        if not os.path.exists(PATH + SEPERATOR + "secretion_uptake"):
            os.mkdir(PATH + SEPERATOR + "secretion_uptake")
    except:
        pass

    # Load configurations over previous runs
    if os.path.exists(PATH + SEPERATOR + ".configs"):
        with open(PATH + SEPERATOR + ".configs", "rb") as f:
            old_cfg = pickle.load(f)
    else:
        old_cfg = cfg

    with open(PATH + SEPERATOR + ".configs", "wb") as f:
        pickle.dump(cfg, f)

    if cfg.analysis.require_setup:
        model_PATH = PATH_res + SEPERATOR + "setup" + SEPERATOR + "snm3_models"
        assert os.path.exists(
            model_PATH
        ), "We require the setup to run first if you set analysis.require_setup=true!"
        log.info("Loading models")
        models = get_models(
            "snm3_models", prefix=PATH_res + SEPERATOR + "setup" + SEPERATOR
        )
    else:
        log.info("Loading original models, WITHOUT setup, as specified.")
        models = get_models(cfg.setup.models)

    done_dfs = []
    done_models = []
    if cfg.analysis.fva_fraction == old_cfg.analysis.fva_fraction:
        for model in models:
            if check_for_substring_in_folder(
                PATH + SEPERATOR + "flux_analysis", model.id, type="*.csv"
            ):
                log.info(f"Fva was already done for {model.id}")
                df = pd.read_csv(
                    PATH
                    + SEPERATOR
                    + "flux_analysis"
                    + SEPERATOR
                    + model.id
                    + "_fva_fba_results"
                    + ".csv",
                    index_col=0,
                )
                log.info(df)
                done_dfs.append(df)
                done_models.append(model)

    log.info("Generating fva results")
    remaining_models = [m for m in models if m not in done_models]
    remaining_dfs = compute_fvas(remaining_models, cfg.analysis.fva_fraction)
    for i, (model, df) in enumerate(zip(remaining_models, remaining_dfs)):
        sol = model.optimize()
        df["flux"] = sol.fluxes
        df.to_csv(
            PATH
            + SEPERATOR
            + "flux_analysis"
            + SEPERATOR
            + model.id
            + "_fva_fba_results"
            + ".csv"
        )
        log.info(df)
    models = done_models + remaining_models
    dfs = done_dfs + remaining_dfs

    log.info("Computing COMPM media for models")
    mediums = compute_COMPM(models, dfs)
    for model, medium in zip(models, mediums):
        with open(
            PATH + SEPERATOR + "medium" + SEPERATOR + model.id + "_COMPM" + ".json", "w"
        ) as f:
            json.dump(medium, f)

    log.info("Plotting fva results")
    for model, df in zip(models, dfs):
        if (
            cfg.analysis.fva_fraction == old_cfg.analysis.fva_fraction
            and check_for_substring_in_folder(
                PATH + SEPERATOR + "flux_analysis",
                model.id + "_full_fba_plot",
                type="*.pdf",
            )
        ):
            log.info("Plots are already computed")
            continue
        plot_full_fva(
            df,
            PATH
            + SEPERATOR
            + "flux_analysis"
            + SEPERATOR
            + model.id
            + "_full_fba_plot"
            + ".pdf",
        )
        if (
            cfg.analysis.fva_fraction == old_cfg.analysis.fva_fraction
            and check_for_substring_in_folder(
                PATH + SEPERATOR + "flux_analysis",
                model.id + "_medium_fva_plot",
                type="*.pdf",
            )
        ):
            log.info("Plots are already computed")
            continue
        plot_medium_fva_range(
            model,
            PATH
            + SEPERATOR
            + "flux_analysis"
            + SEPERATOR
            + model.id
            + "_medium_fva_plot"
            + ".pdf",
            cfg.analysis.fva_fraction,
        )

    log.info("Computing Uptake Secretions + Transports")
    uptakes = []
    sekretions = []
    for i, model in enumerate(models):
        # Transport reactions
        if cfg.analysis.check_transport:
            if check_for_substring_in_folder(
                PATH + SEPERATOR + "secretion_uptake", model.id + "_transport_summary"
            ):
                pass
            else:
                transport_check = table_ex_transport(model)
                transport_check.to_csv(
                    PATH
                    + SEPERATOR
                    + "secretion_uptake"
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
            f"Model: {model.id} has {len(uptake)} uptakes and {len(sekretions)} secretions."
        )
        uptakes.append(uptake)
        sekretions.append(sekretion)

        # Plots for uptake sensitivity!

        if (
            cfg.analysis.sekretion_uptake == old_cfg.analysis.sekretion_uptake
            and not check_for_substring_in_folder(
                PATH + SEPERATOR + "growth", model.id + "_uptake_growth_sensitivity"
            )
        ):
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
            if not check_for_substring_in_folder(
                PATH + SEPERATOR + "secretion_uptake",
                f"{models[i].id}_{models[j].id}_uptake_secretion_summary",
            ):
                log.info("Already computed secretion uptake")
                continue
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
                + "secretion_uptake"
                + SEPERATOR
                + f"{models[i].id}_{models[j].id}_uptake_secretion_summary.csv"
            )

    log.info(f"Plotting venn diagrams for uptake secretion overlaps")
    if not check_for_substring_in_folder(
        PATH + SEPERATOR + "secretion_uptake",
        f"uptake_secretion_overlap_plot",
    ):

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
            + "secretion_uptake"
            + SEPERATOR
            + "uptake_secretion_overlap_plot.pdf"
        )
    else:
        log.info("Already computed secretion uptake plots")

    log.info(f"Plotting expected community interactions")
    if not check_for_substring_in_folder(
        PATH + SEPERATOR + "secretion_uptake",
        f"expected_community_interactions",
    ):
        fig = expected_community_interaction(
            models,
            uptakes,
            sekretions,
            names=cfg.visualization.names,
        )
        fig.savefig(
            PATH
            + SEPERATOR
            + "secretion_uptake"
            + SEPERATOR
            + "expected_community_interactions.pdf"
        )
    else:
        log.info("Already computed expecate community interactions")

    if not check_for_substring_in_folder(
        PATH + SEPERATOR + "secretion_uptake",
        f"cross_feeding_plot",
    ):
        log.info("Computing cross feeding")
        fig = cross_feed_venn_diagrams(
            models,
            uptakes,
            sekretions,
            names=cfg.visualization.names,
            cmap=cfg.visualization.cmap,
        )
        fig.savefig(
            PATH + SEPERATOR + "secretion_uptake" + SEPERATOR + "cross_feeding_plot.pdf"
        )

    # if models_mask.sum() < len(models):
    if not check_for_substring_in_folder(
        PATH + SEPERATOR + "similarity",
        f"jacard_similarities",
    ):
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

    log.info("Computed scaled medium growth plot")
    if not check_for_substring_in_folder(
        PATH + SEPERATOR + "growth",
        f"scaled_medium_growth_plot",
    ):
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
    else:
        log.info("Already computed scaled medium growth plot")

    end_time = time.time()
    log.info(f"Job finished in {end_time-start_time} seconds")
