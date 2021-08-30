from pandas.core.algorithms import isin
import pytest


from ncmw.utils import get_models
import numpy as np
from pandas import DataFrame

from ncmw.community.community_models import (
    BagOfReactionsModel,
    ShuttleCommunityModel,
    create_stoichiometry_matrix,
)

MODELS = get_models("models")
COMMUNITY_MODELS = [BagOfReactionsModel, ShuttleCommunityModel]


@pytest.mark.slow
@pytest.mark.parametrize("community", COMMUNITY_MODELS)
def test_community_models(community):
    model = community(MODELS)
    N = len(model.models)
    growth = model.slim_optimize()
    assert growth > 0
    growth2, single_growths = model.optimize()
    assert abs(growth2 - growth) < 1e-3
    # Test necessary functionality
    try:
        for i in range(len(MODELS)):
            growth_single = model.single_optimize(i)
            growth_single_ref = MODELS[i].slim_optimize()
            assert growth_single is not None
            assert growth_single_ref is not None
        assert abs(sum(single_growths) - growth2) < 1e-3
    except:
        assert False
    try:
        summary = model.summary()
        assert isinstance(summary, DataFrame)
        coopm = model.computeCOOPM(0.1)
        assert isinstance(coopm, dict)
        coopm = model.computeCOOPM(growth, enforce_survival=False)
        assert isinstance(coopm, dict)
        summary = model.compute_convex_combination(np.ones(N) / N, maxbiomass=0.1)
        assert np.isclose(growth, 0.1)
        assert isinstance(summary, DataFrame)
    except:
        assert False


@pytest.mark.parametrize("model", MODELS)
def test_create_stoichiometry_matrix(model):
    S, met, rec = create_stoichiometry_matrix(model)
    shape = S.shape
    # assert len(met) == len(model.metabolites)
    assert len(rec) == len(model.reactions)
    assert shape[1] == len(rec)
