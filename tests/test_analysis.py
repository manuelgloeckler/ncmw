import pytest

from ncmw.setup_models.utils import (
    get_default_medium,
    get_biomass_reaction,
    get_default_configs,
    get_default_medium,
)

from ncmw.workflows.utils import get_models

from ncmw.analysis.medium import compute_COMPM, compute_fvas
from ncmw.analysis.similarity import (
    jaccard_similarity,
    jaccard_similarity_matrices,
    resource_overlap,
)

MODELS = get_models("models")


@pytest.mark.slow
def test_compute_COMPM_and_compute_fvas():
    fvas = compute_fvas(MODELS, 1.0)
    compms = compute_COMPM(MODELS, fvas)

    for compm1 in compms:
        for compm2 in compms:
            for key in compm1:
                if key in compm2:
                    assert compm1[key] == compm2[key]
                assert compm1[key] <= 20
                assert compm1[key] >= 0


def test_jaccard_similarity_and_resourve_overlap():
    assert jaccard_similarity(MODELS[0], MODELS[0]) == (1.0, 1.0)
    assert jaccard_similarity(MODELS[0], MODELS[1]) != (0.0, 0.0)
    assert resource_overlap(MODELS[0], MODELS[0]) == 1.0
    assert resource_overlap(MODELS[0], MODELS[1]) != 0.0


def test_jaccard_similarity_matrices():
    df1, df2, df3 = jaccard_similarity_matrices(MODELS)
    for model1 in MODELS:
        for model2 in MODELS:
            if model1 == model2:
                assert df1[model1.id].loc[model2.id] == 1.0
            else:
                assert df1[model1.id].loc[model2.id] != 1.0
                assert df1[model1.id].loc[model2.id] >= 0
