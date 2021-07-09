
import cobra 
from cobra import Model
import re

import pandas as pd
import numpy as np

def transport_reactions(model:Model):
    """This function return a list of potential transport reactions

    NOTE: This is a cheap heuristic

    Args:
        model (Model): A cobra model

    Returns:
        list: List of names that potentially are transport reaction
    """
    compartment_name = ["_" + id for id in model.compartments.keys()]
    res = []
    for rec in model.reactions:
        for i,c1 in enumerate(compartment_name):
            for c2 in compartment_name[i+1:]:
                if c1 in rec.reaction and c2 in rec.reaction:
                    res.append(rec.id)
    return res


def table_ex_transport(model:Model):
    """This method heuristically checks if all exchange reaction has an associated transporter

    Args:
        model (Model): Cobra model

    Returns:
        DataFrame: Table of indicators (0 indicates abscence, 1 indicates presence)
    """
    compartments = [id for id in model.compartments.keys()]
    metabolites_ex = [key[3:-2] for key in model.medium]
    metabolites_comp = []
    transport_reaction = transport_reactions(model)
    for c in compartments:
        metabolites_comp.append([met for met in model.metabolites if c in met.compartment])
    df = dict(zip(metabolites_ex, [[0 for _ in range(len(compartments))] for _ in range(len(metabolites_ex))]))
    
    for met in metabolites_ex:
        met_id = re.compile(str(met) + "_.")
        hits = []
        for met_c in metabolites_comp:
            hits.append(list(filter(lambda x: re.match(met_id,x.id),met_c)))
        for i,hits_c in enumerate(hits):
            for hit in hits_c:
                for rec in [rec.id for rec in hit.reactions]:
                    if rec in transport_reaction:
                        df[met][i] = 1
    df = pd.DataFrame(df).T
    df.columns = compartments
    return df

def sekretion_uptake_fba(model:Model):
    """This gives the uptake and sekreation reaction in a FBA solution

    NOTE: This is not unique!

    Args:
        model (Model): A cobra model

    Returns:
        list: List of uptake reactions 
        list: List of sekretion reactions
    """
    summary = model.summary()
    uptake = [id for id in summary.uptake_flux.index if summary.uptake_flux.loc[id]["flux"] > 0]
    sekretion = [id for id in summary.secretion_flux.index if summary.secretion_flux.loc[id]["flux"] < 0]
    return uptake, sekretion 

def pad_dict_list(dict_list, padel):
    lmax = 0
    for lname in dict_list.keys():
        lmax = max(lmax, len(dict_list[lname]))
    for lname in dict_list.keys():
        ll = len(dict_list[lname])
        if  ll < lmax:
            dict_list[lname] += [padel] * (lmax - ll)
    return dict_list

def sekretion_uptake_fva(fva):
    """This computes the uptake and sekreation reaction using FVA, this is UNIQUE!

    Args:
        fva (DataFrame): Fva results

    Returns:
        list: List of uptake reactions.
        list: List of sekretion reactions.
    """
    ex_fva = fva.loc[fva.index.str.contains("EX_")]
    uptake = ex_fva[ex_fva["minimum"] < 0].index.tolist()
    sekretion = ex_fva[ex_fva["maximum"] > 0].index.tolist()
    return uptake, sekretion

def compute_uptake_sekretion_table(model_name1, model_name2, uptake1, uptake2, sekretion1, sekretion2):
    # Common up/sekretion from SA to DP
    sek2_up1 = []
    for sek in sekretion2:
        for up in uptake1:
            if str(sek) == str(up):
                sek2_up1.append(str(sek))
    sek1_up2 = []
    for sek in sekretion1:
        for up in uptake2:
            if str(sek) == str(up):
                sek1_up2.append(str(sek))
    df_dict = {f"{model_name1} Uptake":uptake1,f"{model_name1} Secretion":sekretion1, f"{model_name1} -> {model_name2}":sek1_up2, f"{model_name2} -> {model_name1}":sek2_up1, f"{model_name2} Secretion":sekretion2, f"{model_name2} Uptake":uptake2}
    df_dict = pad_dict_list(df_dict, "na")
    df = pd.DataFrame(df_dict) 
    return df

def compute_all_uptake_sekretion_tables(models:list, fvas = None):
    N = len(models)
    if fvas is not None:
        assert len(fvas) == len(models)
    dfs = []
    for i in range(N):
        for j in range(i+1,N):
            if fvas is None:
                df = compute_uptake_sekretion_table(models[i].id,models[j].id, *sekretion_uptake_fba(models[i]), *sekretion_uptake_fba(models[j]))
            else:
                df = compute_uptake_sekretion_table(models[i].id,models[j].id, *sekretion_uptake_fva(fvas[i]), *sekretion_uptake_fva(fvas[j]))
            dfs.append(df)
    return dfs

def compute_uptake_growth_relationship(ids, model, h=100, upper_bound=None):
    """ Computes the growth for multiple fixed flux values, while keeping the others variable."""
    medium = model.medium.copy()
    fluxes = []
    growths = []
    for up in ids:
        old_f = medium[up]
        flux = np.linspace(old_f,0 , h)
        if not (upper_bound ==None):
            flux = np.linspace(upper_bound,0,h)
        growth = []
        with model:
            for f in flux:
                medium[up] = f
                model.medium = medium 
                growth.append(model.slim_optimize())
            medium[up] = old_f
            fluxes.append(flux)
            growths.append(np.array(growth))
    return fluxes, growths
