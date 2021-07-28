import pytest

from ncmw.setup_models.utils import (
    get_default_medium,
    get_biomass_reaction,
    get_default_configs,
    get_default_medium,
)

from ncmw.workflows.utils import get_models

from community.community_models import BagOfReactionsModel, ShuttleCommunityModel

MODELS = get_models("models")
COMMUNITY_MODELS = [BagOfReactionsModel, ShuttleCommunityModel]


@pytest.mark.parametrize("community", COMMUNITY_MODELS)
def test_community_models(community):
    model = community(MODELS)

    growth = model.slim_optimize()
    assert growth > 0
    growth2, single_growths = model.optimize()
    assert growth2 == growth
    assert sum(single_growths) == growth2
    summary = model.summary()
