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

file_dir = os.path.dirname(os.path.dirname(__file__))
sys.path.append(file_dir)


from ncmw.utils.utils_io import (
    DATA_PATH,
    get_models,
    get_result_path,
    SEPERATOR,
    save_model,
    get_model_paths,
)
from ncmw.setup_models import *


@hydra.main(config_path="../../data/hydra", config_name="config.yaml")
def run_setup_hydra(cfg: DictConfig) -> None:
    run_setup(cfg)


def run_setup(cfg: DictConfig) -> None:
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
    except ValueError:
        raise ValueError("Could not generate output directory, maybe we do not have permission's to do so?")
        
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

    if cfg.setup.fastcc:
        log.info(f"Improve models with FastCC: {cfg.setup.fastcc}")
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

    files = []
    gapfill_dict = {"Id": [], "Additions": [], "Growth": []}
    for i, model_i in enumerate(models):
        if cfg.setup.set_bounds_and_medium:
            log.info(
                f"Set default configs {cfg.setup.configs} and medium {cfg.setup.medium} for {model_i.id}"
            )
            model = set_default_configs_and_snm3_medium(
                deepcopy(model_i), cfg.setup.configs, cfg.setup.medium
            )
        else:
            log.info(f"Keep model {model_i.id} as they are")
            model = model_i

        growth = model.slim_optimize()
        log.info(f"Growth on medium: {growth}")
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
        log.info("Computing memote reports, this can take several minutes...")
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

