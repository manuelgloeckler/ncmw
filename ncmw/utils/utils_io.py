import os
import glob
import json

import pandas as pd
import numpy as np

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


def get_reference_network_weights():
    df = pd.read_csv(DATA_PATH + "configs/network_map.csv", sep=";")

    organisms = list(set(df["SOURCE"].tolist()))
    N = len(organisms)

    positive = np.zeros((N, N))
    negative = np.zeros((N, N))
    for i in range(N):
        for j in range(N):
            org1 = organisms[i]
            org2 = organisms[j]
            df_org = df[df["SOURCE"] == org1]
            org1_to_org2 = df_org["TARGET"] == org2
            directed = df_org[org1_to_org2]["DIRECTED"] == "TRUE"
            undirected = df_org[org1_to_org2]["DIRECTED"] == "FALSE"
            p = df_org[org1_to_org2]["INTERACTION"] == "p"
            n = df_org[org1_to_org2]["INTERACTION"] == "n"
            positive[i, j] += p.sum()
            positive[j, i] += p[undirected].sum()
            negative[i, j] += n.sum()
            negative[j, i] += n[undirected].sum()

    weights = positive - negative
    return organisms, weights


def get_model_paths(folder: str):
    """Gets all paths for the models."""
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
    """Saves a model in sbml format."""
    write_sbml_model(model, path)


def get_mediums(folder, prefix=DATA_PATH):
    """Returns the default medium"""
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


def check_for_substring_in_folder(PATH: str, substring: str, type: str = "*") -> bool:
    for filepath in glob.iglob(PATH + SEPERATOR + type):
        if substring in filepath:
            return True

    return False


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
        elif filename.startswith("."):
            pass
        else:
            raise ValueError("Unknown format")
        models.append(model)
    return models
