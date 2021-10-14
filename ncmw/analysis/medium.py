from pandas import DataFrame
import pandas as pd
from cobra.medium import minimal_medium
from cobra.flux_analysis import flux_variability_analysis

from typing import Iterable, Callable, List
import numpy as np


def compute_fvas(models: Iterable, fraction: float) -> List:
    """Compute the FVA result for all models

    Args:
        models: List of models
        fraction: Fraction of maximal biomass rate that must be achived

    Returns:
        list: List of DataFrames containing minimum and maximum fluxes

    """
    dfs = []
    for model in models:
        fva = flux_variability_analysis(model, fraction_of_optimum=fraction)
        dfs.append(fva)
    return dfs


def compute_COMPM(models: list, fvas: list = None):
    r"""Computes the COMPM medium, given all the fva results. The COMPM is defined as the
       medium, in which all models can achive their maximum biomass rate (MBR) if they
       are alone.

       Mathematically it requires the FVA results for each reaction contained in the
       medium! This contains the minimum flux required by the models to obtain
       MBR. We denote it by $FVA_{min; r}^{m}$ for any reaction $r$ of model $m$.

       For any reaction $r \in M$ within the medium $M$ we then define the COMPM as
       $$ COMPM = \{ \min_{m \in Models} FVA_{min;r}^{m} \}_{r \in M}$$

    Args:
        fvas (list): [description]
        reduction (Callable, optional): [description]. Defaults to min.
    """
    if fvas is None:
        fvas = compute_fvas(models, 1.0)
    print(fvas)
    df = pd.concat(list(fvas))
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
            assert np.isclose(
                max_growth, growth
            ), "In the COMPM medium all community members must reach maximum growth rate, check the input!"

    return mediums


def compute_exchange_medium(models, fvas, default=-10):
    if fvas is None:
        fvas = compute_fvas(models, 1.0)
    df = pd.concat(fvas)
    df_min = df.groupby(df.index).min()
    df_max = df.groupby(df.index).max()

    exchanges = set()
    for model in models:
        for ex in model.exchanges:
            exchanges.add(ex.id)
    # TODO This may is wrong...
    mediums = []
    for model in models:
        medium = dict()
        for ex in list(exchanges):
            if ex in df_min.index and ex in df_max.index:
                if df_min.loc[ex].minimum != 0 or df_min.loc[ex].maximum != 0:
                    if ex in [ex.id for ex in model.exchanges]:
                        medium[ex] = default
