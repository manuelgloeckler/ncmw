from pandas import DataFrame
import pandas as pd
from cobra.medium import minimal_medium
from cobra.flux_analysis import flux_variability_analysis

from typing import Iterable, Optional, List
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


def compute_COMPM(models: list, fvas: Optional[list] = None) -> List:
    r"""Computes the COMPM medium, given all the fva results. The COMPM is defined as
    the medium, in which all models can achive their maximum biomass rate (MBR) if they are alone.

    Mathematically it requires the FVA results for each reaction contained in the
    medium! This contains the minimum flux required by the models to obtain
    MBR. We denote it by :math:`\text{FVA}_{min; r}^{m_i}` for any reaction :math:`r` of model :math:`m_i`.

    For reactions :math:`r \in M` within the medium :math:`M` of :math:`K` models we
    define the COMPM as

    .. math:: COMPM = \left\{ \min_{k = 1...K} \text{FVA}_{min;r}^{m_k} \right\}_{r \in M}

    Args:
        models (list): List of cobra metabolit models
        fvas (list, optional): List of dataframes containing FVA results for the models.
                               If argument is "None", then it will be recomputed based on the given models.


    Returns:
        mediums (list): List of COMPM mediums for each model
    """
    if fvas is None:
        fvas = compute_fvas(models, 1.0)
    df = pd.concat(list(fvas))
    df = df.groupby(df.index).min()
    mediums = []
    for model in models:
        medium = model.medium
        for key in medium:
            if key in df.index:
                flux = df.loc[key].minimum
                if flux < 0:
                    medium[key] = -flux
                else:
                    medium[key] = 0
        mediums.append(medium)

    for model, medium in zip(models, mediums):
        with model as m:
            max_growth = m.slim_optimize()
            m.medium = medium
            growth = m.slim_optimize()
            assert (
                growth >= max_growth - 1e-10
            ), "In the COMPM medium all community members must reach maximum growth rate, check the input!"

    return mediums

