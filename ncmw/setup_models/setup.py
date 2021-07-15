import cobra
import cobra.test
from cobra.core import Model
import pandas as pd

import subprocess

from .utils import get_default_configs, get_default_medium, get_biomass_reaction


def gapfill_model(model: Model, eps=1e-6, fill_model_base="ecoli"):
    """Adds reactions to the model, such that it has growth in the given medium

    Args:
        model (Model): Cobra model
        eps ([type], optional): Minimum growth value. Defaults to 1e-6.
        fill_model_base (str, optional): The base set of reactions to consider . Defaults to "ecoli".

    Returns:
        (Model): Cobra model that has growth if algorithm succeeds
        (List): List of reactions that were added
    """
    growth = model.slim_optimize()
    if growth > eps:
        # Already has growth gapfilling is unnecessary
        return model
    if isinstance(fill_model_base, Model):
        fill_model = fill_model_base
    elif fill_model_base == "ecoli" or fill_model_base == "salmonella":
        test_model = cobra.test.create_test_model(fill_model_base)
        fill_model = cobra.Model("universal_reactions")
        fill_model.add_reactions(test_model.reactions)
    solution = cobra.flux_analysis.gapfill(
        model, fill_model, demand_reactions=False, iterations=1
    )[0]
    model.add_reactions(solution)

    assert model.slim_optimize() > eps, "The model still have no growth..."

    return model, solution


def gapfill_medium(model: Model, eps=1e-6):
    """This will add the minimal set of exchange reactions such that the model
    has more than eps growth.

    Args:
        model (Model): Cobra model which has less than eps growth
        eps ([type], optional): Value for which we consider the model to have zero growth . Defaults to 1e-6.

    Returns:
        (Model): Cobra model with extended medium
        (List): List of extended metabolites
    """
    model_help = model.copy()
    # if model_help.slim_optimize() > eps:
    #     # Already feasible model.
    #     return model, []
    # We can gapfill any exchange reaction that currently is not in the medium
    gapfillable = set([ex.id for ex in model_help.exchanges]).difference(
        set(model.medium.keys())
    )
    print(f"There are {len(gapfillable)} many metabolites to fill the medium")

    biomass = get_biomass_reaction(model_help)
    # Binary variables: Theta_i
    # This is an indicator which is zero if the metabolite should be added to the medium
    thetas = []
    for i in range(len(gapfillable)):
        thetas.append(model_help.problem.Variable("theta_" + str(i), type="binary"))

    # Constraints for exchanges, which are turned of for theta_i = 1
    theta_constraints = []
    for i, id in enumerate(gapfillable):
        reaction = model_help.reactions.get_by_id(id)
        min_bound = -10
        reaction.lower_bound = min_bound
        cons = model_help.problem.Constraint(
            (reaction.flux_expression + min_bound * thetas[i]), lb=min_bound, ub=1000
        )
        theta_constraints.append(cons)

    # Constraints for growth rates, which must be at least 10% MBR
    constraint_growth = model_help.problem.Constraint(
        biomass.flux_expression, lb=eps, ub=1000
    )

    # Adding new variables and constraints.
    model_help.add_cons_vars(thetas)
    model_help.add_cons_vars(theta_constraints)
    model_help.add_cons_vars(constraint_growth)

    # Objevtive is maximising turned of exchanges, that is sum of theta_is
    objective = model_help.problem.Objective(sum(thetas), direction="max")
    model_help.objective = objective
    model_help.solver.update()

    # Model optimization
    sol = model_help.optimize()
    # Change medium and check if it worked
    new_exchanges = [
        ex.id
        for ex in model_help.exchanges
        if ex.flux < 0 and ex.id not in model.medium
    ]
    extended_medium = model.medium
    for id in new_exchanges:
        extended_medium[id] = 10
    model.medium = extended_medium
    # assert model.slim_optimize() > eps, "The medium extension failed for some reason..."

    return model, new_exchanges


def set_default_configs_and_snm3_medium(
    model: Model, configs: str = "default.json", medium="snm3.json"
):
    """This

    Args:
        model (Model): A cobra model
        configs (str): File name of a config file in json format
        medium (str): File name of a medium file in json format

    Returns:
        (Model): Cobra model
    """
    # Set bounds
    configs_dict = get_default_configs(configs)
    medium_dict = get_default_medium(medium)
    for key, val in configs_dict.items():
        key = key.split(".")
        if "reactions" in key[0]:
            for reaction in model.reactions:
                reaction_dic = reaction.__dict__
                if reaction_dic["_" + key[1]] != 0:
                    reaction_dic["_" + key[1]] = val
        else:
            reaction = model.reactions.get_by_id(key[0])
            reaction_dic = reaction.__dict__
            if reaction_dic["_" + key[1]] != 0:
                reaction_dic["_" + key[1]] = val

    # Set medium
    exchanges = [ex.id for ex in model.exchanges]
    model_medium = dict()
    for key in medium_dict:
        if key in exchanges:
            model_medium[key] = medium_dict[key]
    model.medium = model_medium
    return model


def score_memote(model_path: str, report_path: str, solver_timout: int = 120):
    """Generates a memote evaluation report for the model quality.

    NOTE: This typically requires a rather consistent model, otherwise it breaks
    NOTE: This can take a while on a local computer, especially if the solver_timeout is high

    Args:
        model_path (str): Path to the model file: Typically SBML format
        report_path (str): Path to the model file: Typically SBML format
        solver_timout (int): Time in seconds until the solver timeouts.
    """
    try:
        p = subprocess.Popen(
            [
                "memote",
                "report",
                "snapshot",
                "--filename",
                report_path,
                "--solver-timeout",
                solver_timout,
                model_path,
            ]
        )
        return p
    except ValueError(
        "It seems that the model cannot be evaluated or the installation of memote is broken..."
    ):
        print("You may consider to score the model online: https://memote.io/")


def create_consistent_model(model: Model):
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
    blocked_reactions_consistent = cobra.flux_analysis.find_blocked_reactions(
        consistent_model
    )
    met_formulas_consistent = cobra.manipulation.check_metabolite_compartment_formula(
        consistent_model
    )
    mass_charge_balance_consistent = cobra.manipulation.check_mass_balance(
        consistent_model
    )

    data = {
        "Blocked reactions": [
            len(blocked_reactions),
            len(blocked_reactions_consistent),
        ],
        "Metabolite formula problems": [
            len(met_formulas),
            len(met_formulas_consistent),
        ],
        "Mass charge balance violations": [
            len(mass_charge_balance),
            len(mass_charge_balance_consistent),
        ],
    }
    df = pd.DataFrame(data, index=["Original Model", "Consistent model"])

    return consistent_model, df


def fastcc(model: Model):
    """FastCC algorithm to increase model quality by resolving conflicts and removing
    unnecessary e.g. blocked pathways
    https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1003424

    Returns:
        cobra.core.Model: Returns a more consistent cobra model.
    """
    consistent_model = cobra.flux_analysis.fastcc(model)
    return consistent_model
