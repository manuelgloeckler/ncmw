import os
import sys
import json

SEPERATOR = os.path.sep
PATH = (
    os.path.dirname(os.path.dirname(os.path.dirname(os.path.realpath(__file__))))
    + SEPERATOR
)
DATA_PATH = PATH + "data" + SEPERATOR
RESULT_PATH = PATH + "results" + SEPERATOR


def get_default_configs(configs: str = "default.json"):
    """Returns the config dictionary"""
    with open(DATA_PATH + "configs" + SEPERATOR + configs, "r") as f:
        configs_dict = json.load(f)
    return configs_dict


def get_default_medium(medium: str = "snm3.json"):
    """Returns the default medium"""
    with open(DATA_PATH + "medium" + SEPERATOR + medium, "r") as f:
        medium_dict = json.load(f)
    return medium_dict


def get_biomass_reaction(model):
    """Return the biomass reaction of a model"""
    objective_str = str(list(model.objective.variables)[0])
    for rec in model.reactions:
        if rec.id in objective_str:
            return rec
