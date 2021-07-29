import pytest


from ncmw.utils import get_models

from ncmw.community.community_models import (
    BagOfReactionsModel,
    ShuttleCommunityModel,
    create_stoichiometry_matrix,
)

MODELS = get_models("models")
COMMUNITY_MODELS = [BagOfReactionsModel, ShuttleCommunityModel]


@pytest.mark.parametrize("community", COMMUNITY_MODELS)
def test_community_models(community):
    model = community(MODELS)

    growth = model.slim_optimize()
    assert growth > 0
    growth2, single_growths = model.optimize()
    assert abs(growth2 - growth) < 1e-3
    for i in range(len(MODELS)):
        growth_single = model.single_optimize(i)
        growth_single_ref = MODELS[i].slim_optimize()
    assert abs(sum(single_growths) - growth2) < 1e-3
    summary = model.summary()
    assert len(summary) > 0
    coopm = model.computeCOOPM(growth)
    assert len(coopm) > 0
    coopm = model.computeCOOPM(growth, enforce_survival=False)
    assert len(coopm) > 0


@pytest.mark.parametrize("model", MODELS)
def test_create_stoichiometry_matrix(model):
    S, met, rec = create_stoichiometry_matrix(model)
    shape = S.shape
    # assert len(met) == len(model.metabolites)
    assert len(rec) == len(model.reactions)
    assert shape[1] == len(rec)
