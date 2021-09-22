import hydra
from omegaconf import DictConfig, OmegaConf

import logging
import socket
import importlib

import time

import random
import numpy as np
import pandas as pd
from copy import deepcopy

import sys, os
import glob
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
from ncmw.setup_models import *

from ncmw.visualization import (
    plot_pairwise_growth_relation_per_weight,
    plot_community_interaction,
    plot_posterior_samples_for_observations,
    plot_community_uptake_graph,
    plot_species_interaction,
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
)

from ncmw.utils import (
    get_models,
    get_result_path,
    SEPERATOR,
    save_model,
    get_model_paths,
)

import cobra


def run_setup(cfg: DictConfig) -> None:
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
    PATH = PATH_res + SEPERATOR + "setup"

    try:
        if not os.path.exists(PATH_res):
            os.mkdir(PATH_res)
        if not os.path.exists(PATH):
            os.mkdir(PATH)
        if not os.path.exists(PATH + SEPERATOR + "snm3_models"):
            os.mkdir(PATH + SEPERATOR + "snm3_models")
        if not os.path.exists(PATH + SEPERATOR + "quality_report"):
            os.mkdir(PATH + SEPERATOR + "quality_report")
        if not os.path.exists(PATH + SEPERATOR + "gapfill"):
            os.mkdir(PATH + SEPERATOR + "gapfill")
    except:
        raise ValueError("Could not generate output directory")
    log.info(f"Generating result directory in {PATH}")

    models_folder = cfg.setup.models
    log.info(f"Loading all models in the folder DATA/{models_folder}")
    models = get_models(models_folder)
    for model in models:
        log.info("Id: " + model.id)
        log.info(f"Reactions: {len(model.reactions)}")
        log.info(f"Metabolites: {len(model.metabolites)}")

    # Check if models allready exists then stop
    for filepath in glob.iglob(PATH + SEPERATOR + "snm3_models" + SEPERATOR + "*.xml"):
        models_already_done = []
        for i in range(len(models)):
            if models[i].id in filepath:
                models_already_done.append(models[i])
        for m in models_already_done:
            log.info(f"Alread done {m.id}")
            models.remove(m)

    log.info(f"Improve models with FastCC: {cfg.setup.fastcc}")
    if cfg.setup.fastcc:
        consistent_models = []
        reports = []
        for model in models:
            cmodel, df = create_consistent_model(model)
            consistent_models.append(cmodel)
            reports.append(df)
            log.info(df)
            out_file = (
                PATH
                + SEPERATOR
                + "quality_report"
                + SEPERATOR
                + "fast_cc"
                + model.id
                + ".csv"
            )
            df.to_csv(out_file)

        models = consistent_models

    log.info(f"Set default configs {cfg.setup.configs} and medium {cfg.setup.medium}")
    files = []
    gapfill_dict = {"Id": [], "Additions": [], "Growth": []}
    for i, model_i in enumerate(models):
        if cfg.setup.set_bounds_and_medium:
            model = set_default_configs_and_snm3_medium(
                deepcopy(model_i), cfg.setup.configs, cfg.setup.medium
            )
        else:
            model = model_i

        growth = model.slim_optimize()
        log.info(f"Growth on SNM3: {growth}")
        if growth < cfg.eps:
            log.warning(
                f"The model {model.id} has no growth on selected Medium, we will try a automated gapfilling strategy"
            )

            if cfg.setup.gapfill == "model":
                # Todo add capability to different gapfilling models...

                model, gapfilled_reactions = gapfill_model(model, cfg.eps)
                log.info(
                    f"Following reactions were added to the model: {gapfilled_reactions}"
                )
                gapfill_dict["Additions"].append(gapfilled_reactions)
            elif cfg.setup.gapfill == "medium":
                model, extended_exchanges = gapfill_medium(model, cfg.eps)
                print(model)
                log.info(
                    f"Following metabolites were added to the medium: {extended_exchanges}"
                )
                gapfill_dict["Additions"].append(extended_exchanges)
            elif cfg.setup.gapfill == False:
                pass
            else:
                raise ValueError(
                    "We only support gapfilling strategies 'model' or 'medium'!"
                )
            growth = model.slim_optimize()
            gapfill_dict["Id"].append(model.id)
            gapfill_dict["Growth"].append(growth)
            log.info("You may have to check if the extensions are 'plausible'...")
            log.info(f"Gapfilling succeded, Growth: {growth}")

        if cfg.setup.set_bounds_and_medium:
            model.id += cfg.setup.medium.split(".")[0]
        file = PATH + SEPERATOR + "snm3_models" + SEPERATOR + model.id + ".xml"
        files.append(file)
        save_model(model, file)
    df = pd.DataFrame(gapfill_dict)
    df.to_csv(PATH + SEPERATOR + "gapfill" + SEPERATOR + "gapfill_report.csv")

    if cfg.setup.memote_evaluation:
        # Start two memote jobs in parallel
        log.info("Computing memote reports")
        i = 0
        for file, model in zip(files, models):
            out_file = (
                PATH
                + SEPERATOR
                + "quality_report"
                + SEPERATOR
                + "memote"
                + model.id
                + ".html"
            )
            p = score_memote(
                file, out_file, solver_timout=str(cfg.setup.memote_solver_time_out)
            )
            i += 1
            if (i + 1) % 5 == 0:
                p.wait()

    end_time = time.time()
    runtime = end_time - start_time
    log.info(f"Finished Workflow in {runtime} seconds")


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
            medium_coopm = m.computeCOOPM(m.slim_optimize())
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
                + "_coopm_uptake_flux.pdf"
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
                + "_compm_uptake_flux.pdf"
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

    assert os.path.exists(
        PATH_res + SEPERATOR + "setup" + SEPERATOR + "snm3_models"
    ), "We require the setup to run first."

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

    log.info("Computing Uptake Sekretions + Transports")
    uptakes = []
    sekretions = []
    for i, model in enumerate(models):
        # Transport reactions
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
        PATH + SEPERATOR + "sekretion_uptake"
        if cfg.analysis.sekretion_uptake == "fva":
            uptake, sekretion = sekretion_uptake_fva(dfs[i])
        else:
            uptake, sekretion = sekretion_uptake_fba(model)
        uptakes.append(uptake)
        sekretions.append(sekretion)

    for i in range(len(models)):
        for j in range(i + 1, len(models)):
            uptake_sekretion_table = compute_uptake_sekretion_table(
                models[i].id,
                models[j].id,
                uptakes[i],
                uptakes[j],
                sekretions[i],
                sekretions[j],
            )
            uptake_sekretion_table.to_csv(
                PATH
                + SEPERATOR
                + "sekretion_uptake"
                + SEPERATOR
                + f"{models[i].id}_{models[j].id}_uptake_sekretion_summary.csv"
            )

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
    fig.savefig(PATH + SEPERATOR + "similarity" + SEPERATOR + "similarity_summary.pdf")

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
