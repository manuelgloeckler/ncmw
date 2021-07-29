import pytest
import numpy as np

from ncmw.setup_models.utils import (
    get_default_medium,
    get_biomass_reaction,
    get_default_configs,
    get_default_medium,
)

from ncmw.utils import get_models

from ncmw.analysis.medium import compute_COMPM, compute_fvas
from ncmw.analysis.similarity import (
    jaccard_similarity,
    jaccard_similarity_matrices,
    resource_overlap,
)
from ncmw.analysis.sekretion_uptake import sekretion_uptake_fba, sekretion_uptake_fva, compute_all_uptake_sekretion_tables

MODELS = get_models("models")
FVAS  = compute_fvas(MODELS, 1.0)

@pytest.mark.slow
def test_compute_COMPM_and_compute_fvas():
    fvas = FVAS
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
                assert df2 [model1.id].loc[model2.id] == 1.0
                assert df3[model1.id].loc[model2.id] == 1.0
            else:
                assert df1[model1.id].loc[model2.id] != 1.0
                assert df1[model1.id].loc[model2.id] >= 0

                assert df2[model1.id].loc[model2.id] != 1.0
                assert df2[model1.id].loc[model2.id] >= 0

                assert df3[model1.id].loc[model2.id] != 1.0
                assert df3[model1.id].loc[model2.id] >= 0

@pytest.mark.parametrize("model", MODELS)
def test_sekretion_uptake_fba(model):
    uptake,sekretion = sekretion_uptake_fba(model)
    assert len(uptake) >= 0
    assert len(sekretion) >= 0
    #model.optimize()
    for ex in uptake:
        if "EX_" in ex:
            r = model.reactions.get_by_id(ex)
            assert r.flux <= 0 
    for ex in sekretion:
        if "EX_" in ex:
            r = model.reactions.get_by_id(ex)
            assert r.flux >= 0 
@pytest.mark.parametrize("fva",FVAS)
def test_sekretion_uptake_fva(fva):
    uptake,sekretion = sekretion_uptake_fva(fva)
    assert len(uptake) >= 0
    assert len(sekretion) >= 0

def test_compute_all_uptake_sekretion_tables():
    dfs1 = compute_all_uptake_sekretion_tables(MODELS)
    dfs2 = compute_all_uptake_sekretion_tables(MODELS, FVAS)
    for df1,df2 in zip(dfs1,dfs2):
        assert len(df1) <= len(df2)
