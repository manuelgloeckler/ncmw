import cobra 
import cobra.test
from cobra.core import Model
import pandas as pd 

from subprocess import call
import os
import json


def gapfill_model(model:Model, eps=1e-6, fill_model_base="ecoli"):
    growth = model.slim_optimize()
    if growth > eps:
        # Already has growth gapfilling is unnecessary
        return model
    if isinstance(fill_model_base, Model):
        fill_model = fill_model_base
    elif fill_model_base == "ecoli" or fill_model_base=="salmonella":
        test_model = cobra.test.create_test_model(fill_model_base)
        fill_model = cobra.Model("universal_reactions")
        fill_model.add_reactions(test_model.reactions)
    solution = cobra.flux_analysis.gapfill(model, fill_model, demand_reactions=False, iterations=1)[0]
    model.add_reactions(solution)

    assert model.slim_optimize() > eps, "The model still have no growth..."

    return model


def set_default_configs_and_snm3_medium(model:Model, configs:str="default.json", medium="snm3.json"):
    """This 

    Args:
        model (Model): A cobra model
        configs (str): File name of a config file in json format
        medium (str): File name of a medium file in json format 

    Returns:
        [type]: [description]
    """    
    SEPERATOR = os.path.sep
    PATH = os.path.dirname(os.path.dirname(os.path(__file__))) + SEPERATOR + "data" + SEPERATOR
    configs_dict = json.load(PATH + "configs" + SEPERATOR + configs)
    medium_dict = json.load(PATH + "medium" + SEPERATOR + medium)
    for key, val in configs_dict.items():
        key = key.split(".")
        if "reactions" in key[0]:
            for reaction in model.reactions:
                reaction_dic = reaction.__dict__
                reaction_dic["_" + key[1]] = val
        else:
            reaction = model.reactions.get_by_id(key[0])
            reaction_dic = reaction.__dict__
            reaction_dic["_" + key[1]] = val
    for key in model.medium:
        if key not in model.exchanges:
            del medium_dict[key]
    model.medium = medium_dict
    return model

def score_memote(model_path:str, report_path:str, solver_timout:int=120):
    """Generates a memote evaluation report for the model quality.

    NOTE: This typically requires a rather consistent model, otherwise it breaks 
    NOTE: This can take a while on a local computer, especially if the solver_timeout is high

    Args:
        model_path (str): Path to the model file: Typically SBML format
        report_path (str): Path to the model file: Typically SBML format
        solver_timout (int): Time in seconds until the solver timeouts.
    """    
    call(["memote", "snapshot", "--filename", report_path, "--solver-timeout", solver_timout, model_path])
    



def create_consistent_model(model:Model, verbose:bool=True):
    """This will create a more consistent model

    Args:
        model (cobra.core.Model): Cobra metabolic model class
        verbose (bool, optional): Print out results in console . Defaults to True.

    Returns:
        cobra.core.Model: Returns a more consistent cobra model.
        pd.DataFrame: Returns some statistics that may be imprved.
    """    
    blocked_reactions = cobra.flux_analysis.find_blocked_reactions(model)
    met_formulas = cobra.manipulation.check_metabolite_compartment_formula(model)
    mass_charge_balance = cobra.manipulation.check_mass_balance(model)

    
    if verbose:
        print("Nonconsistent model values:")
        print("Number of blocked reactions: ", blocked_reactions)
        print("Number of problematic metabolite formulas: ", met_formulas)
        print("Number of mass charge balance violations: ", mass_charge_balance)

    consistent_model = fastcc(model)
    consistent_model.id = model.id + "_consistent"
    blocked_reactions_consistent = cobra.flux_analysis.find_blocked_reactions(consistent_model)
    met_formulas_consistent = cobra.manipulation.check_metabolite_compartment_formula(consistent_model)
    mass_charge_balance_consistent = cobra.manipulation.check_mass_balance(consistent_model)

    if verbose:
        print("Consistent model values:")
        print("Number of blocked reactions: ", blocked_reactions_consistent)
        print("Number of problematic metabolite formulas: ", met_formulas_consistent)
        print("Number of mass charge balance violations: ", mass_charge_balance_consistent)

    data = {'Blocked reactions': [blocked_reactions,blocked_reactions_consistent],
            'Metabolite formula problems':[met_formulas,met_formulas_consistent], 
            'Mass charge balance violations':[mass_charge_balance,mass_charge_balance_consistent]}
    df = pd.DataFrmae(data, columns=[model.id, consistent_model.id])

    return consistent_model, df





def fastcc(model:Model):
    """FastCC algorithm to increase model quality by resolving conflicts and removing 
    unnecessary e.g. blocked pathways
    https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1003424

    Returns:
        cobra.core.Model: Returns a more consistent cobra model.
    """ 
    consistent_model = cobra.flux_analysis.fastcc(model)
    return consistent_model

