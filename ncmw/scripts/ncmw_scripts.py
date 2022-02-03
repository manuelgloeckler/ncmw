from .setup_scripts import run_setup as run1
from .analysis_scripts import run_analysis as run2
from .community_scripts import run_community as run3

import hydra
from omegaconf import DictConfig, OmegaConf

import logging
import socket

import time
import random
import numpy as np


@hydra.main(config_path="../../data/hydra", config_name="config.yaml")
def run_ncmw(cfg: DictConfig) -> None:
    log = logging.getLogger(__name__)
    log.info(OmegaConf.to_yaml(cfg))
    log.info(f"Hostname: {socket.gethostname()}")

    seed = cfg.seed
    random.seed(seed)
    np.random.seed(seed)
    log.info(f"Random seed: {seed}")

    start_time = time.time()
    if cfg.run_setup:
        log.info("Running setup...")
        run1(cfg)

    if cfg.run_analysis:
        log.info("Running analysis...")
        run2(cfg)

    if cfg.run_community:
        log.info("Running community...")
        run3(cfg)

    end_time = time.time()
    runtime = end_time - start_time
    log.info(f"Finished Workflow in {runtime} seconds")
