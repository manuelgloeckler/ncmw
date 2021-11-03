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
import json, pickle

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
    plot_posterior_samples_for_observation,
    plot_community_uptake_graph,
    plot_species_interaction,
    plot_community_summary,
    plot_weight_growth_pairplot,
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
    community_weight_posterior,
)

from ncmw.utils import (
    get_models,
    get_result_path,
    check_for_substring_in_folder,
    SEPERATOR,
)


@hydra.main(config_path="../../data/hydra", config_name="config.yaml")
def run_community_hydra(cfg: DictConfig) -> None:
    run_community(cfg)


def run_community(cfg: DictConfig) -> None:
    log = logging.getLogger(__name__)
    log.setLevel(logging.INFO)
    cobra_loger = logging.getLogger()
    cobra_loger.setLevel(logging.ERROR)
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
        if not os.path.exists(PATH + SEPERATOR + "weight_inference"):
            os.mkdir(PATH + SEPERATOR + "weight_inference")
        if not os.path.exists(PATH + SEPERATOR + "experiments_BagOfReactionsModel"):
            os.mkdir(PATH + SEPERATOR + "experiments_BagOfReactionsModel")
        if not os.path.exists(PATH + SEPERATOR + "experiments_ShuttleCommunityModel"):
            os.mkdir(PATH + SEPERATOR + "experiments_ShuttleCommunityModel")
        log.info("Created folder structure successfully")
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

    log.info("Loading models from setup")
    if cfg.community.base_models == "setup":
        try:
            models = get_models(
                "snm3_models", prefix=PATH_res + SEPERATOR + "setup" + SEPERATOR
            )
        except:
            raise ValueError(
                "We require snm3 models for the next steps, please run the setup! Alternative specifiy comunity.base_models with a path to the models you want to use"
            )
    else:
        try:
            models = get_models("", prefix=cfg.community.base_models)
        except:
            raise ValueError(f"No models found at {cfg.community.base_models}")

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
        if (
            check_for_substring_in_folder(
                PATH + SEPERATOR + "community_models", "BagOfReactionsModel.pkl"
            )
            and old_cfg.community.base_models == cfg.community.base_models
        ):
            log.info("Loading BagOfReactions community model")
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
        kwargs = cfg.community.shuttle_reaction_params
        if not isinstance(kwargs["shared_reactions"], list):
            kwargs["shared_reactions"] = None
            old_cfg.community.shuttle_reaction_params["shared_reactions"] = None
        if (
            check_for_substring_in_folder(
                PATH + SEPERATOR + "community_models", "ShuttleCommunityModel.pkl"
            )
            and old_cfg.community.base_models == cfg.community.base_models
            and all(
                [
                    old_cfg.community.shuttle_reaction_params[key]
                    == cfg.community.shuttle_reaction_params[key]
                    for key in cfg.community.shuttle_reaction_params
                ]
            )
        ):
            log.info("Loading shuttle reactions model")
            m = type(m).load(path)
            ids = [mod.id for mod in m.models]
            correct = True
            for model in models:
                correct = correct and (model.id in ids)
            if correct:
                community_models += [m]
            else:
                community_models += [ShuttleCommunityModel(models, **kwargs)]
        else:
            log.info("Building shuttle reactions model")
            community_models += [ShuttleCommunityModel(models, **kwargs)]

    log.info(f"Set correct weights: {cfg.community.main_weights}")
    for m in community_models:
        if cfg.community.main_weights == "ones":
            weights = np.ones(len(m.models))
        elif cfg.community.main_weights == "fair":
            weights = compute_fair_weights(m)
        else:
            weights = cfg.community.main_weights
            assert len(weights) == len(
                m.models
            ), "The custom main weights must be a iterable of floats for each member of the community!"
        m.weights = weights

    log.info(f"Saving models: {cfg.community.save_models}")
    if cfg.community.save_models:
        for m in community_models:
            path = PATH + SEPERATOR + "community_models" + SEPERATOR + type(m).__name__
            m.save(path)
            log.info(f"Saving community models: {path}")
    if cfg.community.save_as_sbml:
        for m in community_models:
            path = (
                PATH
                + SEPERATOR
                + "community_models"
                + SEPERATOR
                + type(m).__name__
                + ".xml"
            )
            m.save_as_sbml(path)
            log.info(f"Saving community models as sbml: {path}")

    # Medias
    all_medias = {}

    if cfg.community.medium["default"]:
        log.info("Computing default for community")
        # Also save default medium
        all_medias["default"] = m.medium
        path = PATH + SEPERATOR + "medium" + SEPERATOR + "DEFAULT" + ".json"
        with open(path, "w+") as f:
            json.dump(m.medium, f)
        log.info(f"Saving default: {path}")

    if cfg.community.medium["compm"] or cfg.community.medium["coopm"]:
        log.info("Computing COMPM for community")
        medium_prefix = PATH_res + SEPERATOR + "analysis" + SEPERATOR
        mediums = get_mediums("medium", medium_prefix)
        COMPM = dict()
        for medium in mediums.values():
            for key, val in medium.items():
                COMPM[key] = float(abs(val))
        all_medias["compm"] = COMPM
        if cfg.community.medium["compm"]:
            path = PATH + SEPERATOR + "medium" + SEPERATOR + "COMPM" + ".json"
            with open(path, "w+") as f:
                json.dump(COMPM, f)
            log.info(f"Saving COMPMs: {path}")

    for m in community_models:
        if cfg.community.medium["coopm"]:
            log.info("Computing COOPM medium")
            kwargs = cfg.community.coopm_params
            m.medium = all_medias["compm"]
            coopm = m.compute_COOPM(m.slim_optimize(), **kwargs)
            all_medias["coopm"] = coopm
        for medium_name, medium in all_medias.items():
            path_to_save = (
                PATH
                + SEPERATOR
                + "experiments_"
                + type(m).__name__
                + SEPERATOR
                + medium_name
            )
            if not os.path.exists(path_to_save):
                os.mkdir(path_to_save)
            log.info(f"Community FBA/FVA results on medium: {medium_name}")
            df_growth_summary = pd.DataFrame()

            m = type(m).load(
                PATH
                + SEPERATOR
                + "community_models"
                + SEPERATOR
                + f"{type(m).__name__}.pkl"
            )
            m.medium = medium

            log.info("Computing FVA")
            fraction_of_optimal = cfg.community.cooperative_tradeoff_params["alpha"]
            if cfg.community.conpute_community_fva:
                df_fva = m.flux_variability_analysis(
                    **cfg.community.community_fva_params
                )
                df_fva.columns = [
                    f"minimum (fraction of optimal: {fraction_of_optimal})",
                    f"maximum (fraction of optimal: {fraction_of_optimal})",
                ]
            else:
                df_fva = pd.DataFrame()
            # TODO after cooperative tradeoff or somthing this gets stuck in a loop for
            # some reason....................
            growth, single_growths, sol = m.optimize(**cfg.community.optimize_params)
            (
                growth_tradeoff,
                single_growths_tradeoff,
                sol_tradeoff,
            ) = m.cooperative_tradeoff(**cfg.community.cooperative_tradeoff_params)
            log.info(
                f"Achieved community growth:{growth}, with individual growth: {single_growths}"
            )
            log.info(
                f"Achieved community growth with cooperative tradeoff:{growth_tradeoff}, with individual growth: {single_growths_tradeoff}"
            )
            df_growth_summary["Weights"] = list(m.weights) + [None]
            df_growth_summary["FBA Growth"] = single_growths + [growth]
            df_growth_summary[
                "Cooperative tradeoff Growth"
            ] = single_growths_tradeoff + [growth_tradeoff]
            df_growth_summary.index = [m.id for m in models] + ["Community growth"]
            df_growth_summary.to_csv(path_to_save + SEPERATOR + f"growth_analysis.csv")

            df_fba = sol.to_frame()
            df_fba_cooperative_tradeoff = sol_tradeoff.to_frame()

            df_fva["FBA"] = df_fba["fluxes"]
            df_fva[
                f"Cooperative tradeoff (fraction of optimal: {fraction_of_optimal}"
            ] = df_fba_cooperative_tradeoff["fluxes"]
            log.info("Saving all flux analysis")
            df_fva.to_csv(path_to_save + SEPERATOR + f"flux_analysis.csv")

            if cfg.community.pairwise_growth:
                # we reload it here as it, check why need
                m = type(m).load(
                    PATH
                    + SEPERATOR
                    + "community_models"
                    + SEPERATOR
                    + f"{type(m).__name__}.pkl"
                )
                m.medium = medium
                log.info("Compute pairwise growth relationships for different weights.")
                fig = plot_pairwise_growth_relation_per_weight(
                    m,
                    cfg.visualization.names,
                    **cfg.community.pairwise_growth_params,
                )
                fig.savefig(
                    path_to_save + SEPERATOR + "pairwise_growth_relationship.pdf"
                )
            if m._type == "compartmentalized" and isinstance(m, ShuttleCommunityModel):
                m = type(m).load(
                    PATH
                    + SEPERATOR
                    + "community_models"
                    + SEPERATOR
                    + f"{type(m).__name__}.pkl"
                )
                m.medium = medium

                summary_1 = m.summary(**cfg.community.optimize_params)
                summary_1.to_csv(path_to_save + SEPERATOR + f"flux_summary.csv")
                fig = plot_community_summary(m, summary_1, cfg.visualization.names)
                fig.savefig(path_to_save + SEPERATOR + f"community_summary.pdf")

                summary_2 = m.summary(
                    cooperative_tradeoff=cfg.community.cooperative_tradeoff_params[
                        "alpha"
                    ]
                )
                summary_2.to_csv(path_to_save + SEPERATOR + f"flux_summary.csv")
                fig = plot_community_summary(m, summary_2, cfg.visualization.names)
                fig.savefig(path_to_save + SEPERATOR + f"community_summary.pdf")

                log.info(f"Compute community visualization")
                fig2 = plot_community_interaction(
                    m, summary_1, names=cfg.visualization.names
                )
                fig2.savefig(path_to_save + SEPERATOR + "community_exchange.pdf")
                fig2 = plot_community_uptake_graph(
                    m, summary_1, names=cfg.visualization.names
                )
                fig2.savefig(
                    path_to_save + SEPERATOR + "community_uptake.png", dpi=1000
                )

                fig2 = plot_species_interaction(
                    m, summary_1, names=cfg.visualization.names
                )
                fig2.savefig(path_to_save + SEPERATOR + "species_interaction.pdf")

                fig2 = plot_community_interaction(
                    m, summary_2, names=cfg.visualization.names
                )
                fig2.savefig(
                    path_to_save
                    + SEPERATOR
                    + "community_exchange_cooperative_tradeoff.pdf"
                )
                fig2 = plot_community_uptake_graph(
                    m, summary_2, names=cfg.visualization.names
                )
                fig2.savefig(
                    path_to_save
                    + SEPERATOR
                    + "community_uptake_cooperative_tradeoff.png",
                    dpi=1000,
                )

                fig2 = plot_species_interaction(
                    m, summary_2, names=cfg.visualization.names
                )
                fig2.savefig(
                    path_to_save
                    + SEPERATOR
                    + "species_interaction_cooperative_tradeoff.pdf"
                )

    log.info("Compute weight posteriors")
    if cfg.community.compute_infer_weights:
        for m in community_models:
            m.weights = np.ones(len(m.weights))
            m.medium = all_medias[cfg.community.infer_weights.medium]
            if cfg.community.infer_weights.observed_individual_biomass == "balanced":
                maximal_growths = [m.single_optimize(i) for i in range(len(m.models))]
                balanced_observation = np.ones(len(m.models)) * np.min(maximal_growths)
                observation = balanced_observation
            else:
                observation = cfg.community.infer_weights.observed_individual_biomass
                assert len(
                    cfg.community.infer_weights.observed_individual_biomass
                ) == len(
                    m.models
                ), "You have to specify for each community member an observed growth."
                for o in observation:
                    assert o >= 0, "The observed growth must be greater than zero!"

            alpha = None
            if cfg.community.infer_weights.competitive_tradeoff:
                alpha = cfg.community.infer_weights.competitive_tradeoff_alpha
            posterior, weights, growths = community_weight_posterior(
                m,
                num_simulations=cfg.community.infer_weights.simulations_for_different_weights,
                enforce_survival=cfg.community.infer_weights.enforce_survival,
                cooperative_tradeoff=alpha,
            )

            fig = plot_weight_growth_pairplot(
                m, weights, growths, names=cfg.visualization.names
            )
            fig.savefig(
                PATH
                + SEPERATOR
                + "weight_inference"
                + SEPERATOR
                + type(m).__name__
                + "weight_growth_pairplot.pdf"
            )

            fig = plot_posterior_samples_for_observation(m, posterior, observation)
            fig.savefig(
                PATH
                + SEPERATOR
                + "weight_inference"
                + SEPERATOR
                + type(m).__name__
                + "weight_posterior.pdf"
            )

    end_time = time.time()
    runtime = end_time - start_time
    log.info(f"Finished Workflow in {runtime} seconds")
