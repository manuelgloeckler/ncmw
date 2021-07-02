from pandas import DataFrame
import pandas as pd 
from cobra.medium import minimal_medium
from cobra.flux_analysis import flux_variability_analysis

from typing import Iterable, Callable

def compute_fvas(models:Iterable, fraction:float):
    dfs = []
    for model in models:
        fva = flux_variability_analysis(model, fraction_of_optimum=fraction)
        dfs.append(fva)
    return dfs

def compute_COMPM(models: list, fvas: list=None):
    """ Computes the COMPM medium, given all the fva results

    Args:
        fvas (list): [description]
        reduction (Callable, optional): [description]. Defaults to min.
    """
    if fvas is None:
        fvas = compute_fvas(models, 1.0)
    df = pd.concat(fvas)
    df = df.groupby(df.index).min()
    mediums = []
    for model in models:
        medium = model.medium
        for key in medium:
            if key in df.index:
                flux = df.loc[key].minimum
                medium[key] = -flux
        mediums.append(medium)

    for model, medium in zip(models, mediums):
        with model as m:
            max_growth = m.slim_optimize()
            m.medium = medium 
            growth = m.slim_optimize()
            assert max_growth == growth, "In the COMPM medium all community members must reach maximum growth rate, check the input!"

    return mediums

        

# def compute_COMPM_pairwise(fva_model1:DataFrame, fva_model2:DataFrame, reduction:Callable=min):
#     """ Computes the COMPM medium, given the fva results of two models.

#     Args:
#         fva_model1 (DataFrame): Dataframe of the first model fva result
#         fva_model2 (DataFrame): Dataframe of the second model fva result
#         reduction (DataFrame): Function f(X,Y) -> Z that reduces two fluxes to one
#                                typically the min for COMPM

#     Raises:
#         ValueError: Invalid value for an exchange reactions

#     Returns:
#         (dict): Medium for the first model.
#         (dict): Medium for the second model.
#     """
#     ex_ids1 = [idx for idx in fva_model1.index if "EX_" in idx]
#     ex_ids2 = [idx for idx in fva_model2.index if "EX_" in idx]
#     ex_union_ids = list(set(ex_ids1).union(set(ex_ids2)))
#     COMPM_1 = dict()
#     COMPM_2 = dict()
#     for ex in ex_union_ids:
#         if ex in fva_model2.index and ex in fva_model1.index:
#             # Take minimum of both minima if common
#             flux = reduction(fva_model2.loc[ex].minimum, fva_model1.loc[ex].minimum)
#             if flux < 0:
#                     COMPM_1[ex] = abs(flux)
#                     COMPM_2[ex] = abs(flux)
#         elif ex in fva_model2.index:
#             # Take just the minimum if not common
#             flux = fva_model2.loc[ex].minimum
#             if flux < 0:
#                     COMPM_1[ex] = abs(flux)
#         elif ex in fva_model1.index:
#             # Take just the minimum if not common
#             flux = fva_model1.loc[ex].minimum
#             if flux < 0:
#                     COMPM_2[ex] = abs(flux)
#         else:
#             raise ValueError("Unknown exchange reaction. Ur fva results seem to be cryptic...")
#     return COMPM_1, COMPM_2 

def minimal_medium(model1, model2,):
    pass