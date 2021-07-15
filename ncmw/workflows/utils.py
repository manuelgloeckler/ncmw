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
    write_sbml_model(model, path)


def get_mediums(folder, prefix=DATA_PATH):
    directory = prefix + folder
    mediums = dict()
    for filename in os.listdir(directory):
        with open(directory + SEPERATOR + filename, "r") as f:
            mediums[filename] = json.load(f)
    return mediums


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
            pass
        models.append(model)
    return models
