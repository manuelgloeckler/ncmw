import os
import sys
import json 

from cobra.io import read_sbml_model, load_json_model, load_yaml_model, load_matlab_model, write_sbml_model

SEPERATOR = os.path.sep
PATH = os.path.dirname(os.path.dirname(os.path.dirname(os.path.realpath(__file__)))) + SEPERATOR 
DATA_PATH = PATH + "data" + SEPERATOR
RESULT_PATH = PATH + "results" + SEPERATOR



def get_result_path(name:str):
    """ Returns the path to the results """
    return RESULT_PATH + name

def get_model_paths(folder:str):
    directory = DATA_PATH + folder
    files = []
    for filename in os.listdir(directory):
        if filename.endswith(".xml") or filename.endswith(".json") or filename.endswith(".yaml") or filename.endswith(".mat"):
            files.append(directory + SEPERATOR + filename)
    return files

def save_model(model, path:str):
    with open(path, "w") as f:
        write_sbml_model(model,f)



def get_models(folder:str):
    """ Returns all models in the directory """
    directory = DATA_PATH + folder
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
            raise NotImplementedError("We only suppport models that follow the sbml, json, yaml or matlab format.")
        models.append(model)
    return models