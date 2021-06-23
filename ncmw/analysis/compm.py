from pandas import DataFrame


def compute_COMPM(fva_model1:DataFrame, fva_model2:DataFrame):
    """ Computes the COMPM medium, given the fva results of two models.

    Args:
        fva_model1 (DataFrame): Dataframe of the first model fva result
        fva_model2 (DataFrame): Dataframe of the second model fva result

    Raises:
        ValueError: Invalid value for an exchange reactions

    Returns:
        (dict): Medium for the first model.
        (dict): Medium for the second model.
    """
    ex_ids1 = [idx for idx in fva_model1.index if "EX_" in idx]
    ex_ids2 = [idx for idx in fva_model2.index if "EX_" in idx]
    ex_union_ids = ex_ids1.union(ex_ids2)
    COMPM_1 = dict()
    COMPM_2 = dict()
    for ex in ex_union_ids:
        if ex in fva_model2.index and ex in fva_model1.index:
            # Take minimum of both minima if common
            flux = min(fva_model2.loc[ex].minimum, fva_model1.loc[ex].minimum)
            if flux < 0:
                    COMPM_1[ex] = abs(flux)
                    COMPM_2[ex] = abs(flux)
        elif ex in fva_model2.index:
            # Take just the minimum if not common
            flux = fva_model2.loc[ex].minimum
            if flux < 0:
                    COMPM_1[ex] = abs(flux)
        elif ex in fva_model1.index:
            # Take just the minimum if not common
            flux = fva_model1.loc[ex].minimum
            if flux < 0:
                    COMPM_2[ex] = abs(flux)
        else:
            raise ValueError("Unknown exchange reaction. Ur fva results seem to be cryptic...")
    return COMPM_1, COMPM_2 