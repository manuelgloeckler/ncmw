
import hydra
from omegaconf import DictConfig, OmegaConf 

import logging
import coloredlogs
import socket
import importlib

import time

import random 
import numpy as np
import os 
import pandas as pd


from utils import get_models, get_result_path, SEPERATOR, save_model, get_model_paths
from setup_models import *

import cobra


@hydra.main(config_path="../../data/hydra", config_name="config.yaml")
def run_setup(cfg : DictConfig) -> None:
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
    PATH = PATH_res + SEPERATOR +  "setup"

    try:
        os.mkdir(PATH_res)
        os.mkdir(PATH)
        os.mkdir(PATH + SEPERATOR + "snm3_models")
        os.mkdir(PATH + SEPERATOR + "quality_report")
    except:
        pass
    log.info(f"Generating result directory directory in {PATH}")


    models_folder = cfg.setup.models
    log.info(f"Loading all models in the folder DATA/{models_folder}")
    models = get_models(models_folder)
    for model in models:
        log.info("Id: " + model.id)
        log.info(f"Reactions: {len(model.reactions)}")
        log.info(f"Metabolites: {len(model.metabolites)}")

    log.info(f"Improve models with FastCC: {cfg.setup.fastcc}")
    if cfg.setup.fastcc:
        consistent_models = []
        reports = []
        for model in models:
            cmodel, df = create_consistent_model(model)
            cobra.io.write_sbml_model(cmodel, PATH + SEPERATOR + "snm3_models" + SEPERATOR + "consistent_" + cmodel.id)
            consistent_models.append(cmodel)
            reports.append(df)
            log.info(df)
            out_file = PATH + SEPERATOR + "quality_report" + SEPERATOR + "fast_cc" + model.id +".csv"
            df.to_csv(out_file)

        models = consistent_models

    
    
    log.info(f"Set default configs {cfg.setup.configs} and medium {cfg.setup.medium}")
    files = []
    for model in models:
        model = set_default_configs_and_snm3_medium(model, cfg.setup.configs, cfg.setup.medium)
        growth = model.slim_optimize()
        log.info(f"Growth on SNM3: {growth}")
        if growth < cfg.eps:
            log.warning(f"The model {model.id} has no growth on selected Medium, we will try a automated gapfilling strategy")

            if cfg.setup.gapfill == "model":
                # Todo add capability to different gapfilling models...
                model = gapfill_model(model, cfg.eps)
            elif cfg.setup.gapfill == "medium":
                model = gapfill_medium(model, cfg.eps)
            elif cfg.setup.gapfill==False:
                pass
            else:
                raise ValueError("We only support gapfilling strategies 'model' or 'medium'!")
            growth = model.slim_optimize()
            log.info(f"Gapfilling succeded, Growth: {growth}")
            
        model.id += "_snm3"
        file = PATH + SEPERATOR + "snm3_models" + SEPERATOR + model.id + ".xml"
        files.append(file)
        save_model(model,file)

        


    if cfg.setup.memote_evaluation:
        # Start two memote jobs in parallel
        log.info("Computing memote reports")
        i = 0
        for file, model in zip(files, models):
            out_file = PATH + SEPERATOR + "quality_report" + SEPERATOR + "memote" + model.id + ".html"
            p = score_memote(file, out_file, solver_timout=str(cfg.setup.memote_solver_time_out))
            i += 1
            if (i+1) % 2 == 0:
                p.wait()

    
    end_time = time.time()
    runtime = end_time-start_time
    log.info(f"Finished Workflow in {runtime} seconds")
    
    
    

if __name__ == "__main__":
    run_setup()





    
