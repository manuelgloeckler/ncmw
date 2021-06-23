import cobra 
import cobra.test
from cobra.core import Model
import pandas as pd 

import subprocess

from .utils import get_default_configs, get_default_medium


def gapfill_model(model:Model, eps=1e-6, fill_model_base="ecoli"):
    """ Adds reactions to the model, such that it has growth in the given medium

    Args:
        model (Model): Cobra model
        eps ([type], optional): Minimum growth value. Defaults to 1e-6.
        fill_model_base (str, optional): The base set of reactions to consider . Defaults to "ecoli".

    Returns:
        (Model): Cobra model that has growth if algorithm succeeds
    """
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

def gapfill_medium(model:Model, eps=1e-6):
    raise NotImplementedError()


def set_default_configs_and_snm3_medium(model:Model, configs:str="default.json", medium="snm3.json"):
    """This 

    Args:
        model (Model): A cobra model
        configs (str): File name of a config file in json format
        medium (str): File name of a medium file in json format 

    Returns:
        (Model): Cobra model
    """    
    configs_dict = get_default_configs(configs)
    medium_dict = get_default_medium(medium)
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
    exchanges = [ex.id for ex in model.exchanges]
    model_medium = dict()
    for key in medium_dict:
        if key in exchanges:
            model_medium[key] = medium_dict[key]
    model.medium = model_medium
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
    try:
        p = subprocess.Popen(["memote", "report", "snapshot", "--filename", report_path, "--solver-timeout", solver_timout, model_path])
        return p
    except ValueError("It seems that the model cannot be evaluated or the installation of memote is broken..."):
        print("You may consider to score the model online: https://memote.io/")
    



def create_consistent_model(model:Model):
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

    consistent_model = fastcc(model)
    consistent_model.id = model.id + "_consistent"
    blocked_reactions_consistent = cobra.flux_analysis.find_blocked_reactions(consistent_model)
    met_formulas_consistent = cobra.manipulation.check_metabolite_compartment_formula(consistent_model)
    mass_charge_balance_consistent = cobra.manipulation.check_mass_balance(consistent_model)



    data = {'Blocked reactions': [len(blocked_reactions),len(blocked_reactions_consistent)],
            'Metabolite formula problems':[len(met_formulas),len(met_formulas_consistent)], 
            'Mass charge balance violations':[len(mass_charge_balance),len(mass_charge_balance_consistent)]}
    df = pd.DataFrame(data, index=["Original Model", "Consistent model"])

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

