import pytest
import cobra

from ncmw.setup_models.utils import (
    get_default_medium,
    get_biomass_reaction,
    get_default_configs,
)


def test_default_configs():
    cfgs = get_default_configs()
    assert (
        cfgs.reactions.lower_bound < cfgs.reactions.upper_bound
    ), "The lower bound must be smaller than the upper one"
