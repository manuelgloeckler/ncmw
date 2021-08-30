import os
import sys
import json

from cobra.io import (
    read_sbml_model,
    load_json_model,
    load_yaml_model,
    load_matlab_model,
    write_sbml_model,
)

SEPERATOR = os.path.sep
PATH = (
    os.path.dirname(os.path.dirname(os.path.dirname(os.path.realpath(__file__))))
    + SEPERATOR
)
DATA_PATH = PATH + "data" + SEPERATOR
RESULT_PATH = PATH + "results" + SEPERATOR


def get_result_path(name: str):
    """Returns the path to the results"""
    return RESULT_PATH + name


def get_model_paths(folder: str):
    """ Gets all paths for the models."""
    directory = DATA_PATH + folder
    files = []
    for filename in os.listdir(directory):
        if (
            filename.endswith(".xml")
            or filename.endswith(".json")
            or filename.endswith(".yaml")
            or filename.endswith(".mat")
        ):
            files.append(directory + SEPERATOR + filename)
    return files


def save_model(model, path: str):
    """ Saves a model in sbml format. """
    write_sbml_model(model, path)


def get_mediums(folder, prefix=DATA_PATH):
    """ Returns the default medium """
    directory = prefix + folder
    mediums = dict()
    for filename in os.listdir(directory):
        with open(directory + SEPERATOR + filename, "r") as f:
            mediums[filename] = json.load(f)
    return mediums

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


def get_models(folder: str, prefix=DATA_PATH):
    """Returns all models in the directory"""
    directory = prefix + folder
    models = []
    for filename in os.listdir(directory):
        if filename.endswith(".xml"):
            model = read_sbml_model(directory + SEPERATOR + filename)
        elif filename.endswith(".json"):
            model = load_json_model(directory + SEPERATOR + filename)
        elif filename.endswith(".yaml"):
            model = load_yaml_model(directory + SEPERATOR + filename)
        elif filename.endswith(".mat"):
            model = load_matlab_model(directory + SEPERATOR + filename)
        else:
            raise ValueError("Unknown format")
        models.append(model)
    return models
