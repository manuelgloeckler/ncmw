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

