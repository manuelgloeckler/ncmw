import pytest
import os

from ncmw.utils import (
    get_default_medium,
    get_biomass_reaction,
    get_default_configs,
    get_default_medium,
    check_for_substring_in_folder,
    get_model_paths,
)
from ncmw.setup_models.setup import (
    gapfill_medium,
    gapfill_model,
    set_default_configs_and_snm3_medium,
    create_consistent_model,
    score_memote,
)
from ncmw.utils import get_models

# Take standard model tests folder for testing
MODELS = get_models("models")


def test_default_configs():
    cfgs = get_default_configs()
    assert (
        cfgs["reactions.lower_bound"] < cfgs["reactions.upper_bound"]
    ), "The lower bound must be smaller than the upper one"


def test_get_model_paths():
    paths = get_model_paths("models")
    assert len(MODELS) == len(paths)


def test_get_default_medium():
    medium = get_default_medium()
    # Test default values
    for key, val in medium.items():
        if "EX_o2_e" in key:
            assert val == 20
        elif "EX_fe" in key:
            assert val == 0.1
        else:
            assert val == 10


def test_check_for_substring_in_folder():
    assert check_for_substring_in_folder(".", "test")
    assert not check_for_substring_in_folder(".", "test", ".txt")
    assert not check_for_substring_in_folder(".", "asdfasdfasdfasdf")


@pytest.mark.parametrize("model", MODELS)
def test_get_biomass_reaction(model):
    biomass_f = get_biomass_reaction(model)
    sol = model.optimize()
    assert sol.objective_value == biomass_f.flux


@pytest.mark.parametrize("model", MODELS)
def test_gap_fill_model(model):
    exchanges = [ex.id for ex in model.exchanges]
    model.medium = dict(
        [(key, val) for key, val in get_default_medium().items() if key in exchanges]
    )
    if model.slim_optimize() < 1e-10:
        model = gapfill_model(model)[0]
        assert model.slim_optimize() > 1e-10


@pytest.mark.parametrize("model", MODELS)
def test_gap_fill_medium(model):
    exchanges = [ex.id for ex in model.exchanges]
    model.medium = dict(
        [(key, val) for key, val in get_default_medium().items() if key in exchanges]
    )
    if model.slim_optimize() < 1e-10:
        model = gapfill_medium(model)[0]
        assert model.slim_optimize() > 1e-10


@pytest.mark.parametrize("model", MODELS)
def test_create_consistent_model(model):
    m, report = create_consistent_model(model)

    # More consistent model should be better or equal in all disciplinces
    data = report.to_numpy()
    assert all(data[0, :] >= data[1, :])
    assert m.slim_optimize() > 1e-10


@pytest.mark.parametrize("model_path", get_model_paths("models"))
def test_score_memote(model_path):
    test_output = (
        os.path.dirname(os.path.join(os.getcwd(), os.listdir(os.getcwd())[0]))
        + "/test.html"
    )
    p = score_memote(model_path, test_output, solver_timout="0")
    p.wait()
    assert os.path.exists(test_output)
    try:
        os.remove(test_output)
    except:
        pass


@pytest.mark.parametrize("model", MODELS)
def test_set_default_configs_and_snm3_medium(model):
    model = set_default_configs_and_snm3_medium(model)
    cfgs = get_default_configs()
    medium = get_default_medium()
    lb = cfgs["reactions.lower_bound"]
    ub = cfgs["reactions.upper_bound"]
    for rec in model.reactions:
        assert rec.lower_bound >= lb
        assert rec.upper_bound <= ub
    for key in model.medium:
        assert key in medium
        assert model.medium[key] == medium[key]
