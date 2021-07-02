from typing import Iterable
from cobra.core import Model
import pandas as pd

import numpy as np

def get_shared_metabolites_counts(model1:Model, model2:Model):
    """This method return the number of unique metabolites in both models .

    Args:
        model1 (Model): First cobra model
        model2 (Model): Second cobra model

    Returns:
        (int): Total number of metabolites in both models
        (int): Total number of shared metabolites
    """    
    met1 = set(map(lambda x:x.id, model1.metabolites))
    met2 = set(map(lambda x:x.id, model2.metabolites))
    met_ids = set(met1)
    met_ids = met_ids.union(set(met2))
    common_met_ids = met1.intersection(met2)

    return len(met_ids), len(common_met_ids)

def get_shared_reactions_counts(model1:Model, model2:Model):
    """ This computes the number of shared reactions

    Args:
        model1 (Model): First cobra model
        model2 (Model): Second cobra model
        prefix (str): Location to write file. Default is ""

    Returns:
        (int): Total number of reactions in both models
        (int): Total number of shared reactions
    """    

    rec1 = set(map(lambda x:x.id, model1.reactions))
    rec2 = set(map(lambda x:x.id, model2.reactions))
    rec_ids = set(rec1)
    rec_ids = rec_ids.union(set(rec2))
    common_rec_ids = rec1.intersection(rec2)

    return len(rec_ids), len(common_rec_ids)

def jaccard_similarity(model1:Model, model2:Model):
    """This returns the Jacard Similarity of both models with respect to the set 
    of metabolites and reactions

    Args:
        model1 (Model): First cobra model
        model2 (Model): Second cobra model

    Returns:
        (float): Jacard similarity of metabolite sets
        (float): Jacard similarity of reaction sets
    """
    total_mets, common_mets = get_shared_metabolites_counts(model1,model2)
    total_recs, common_recs = get_shared_reactions_counts(model1,model2)
    j_met = common_mets/total_mets
    j_rec = common_recs/total_recs
    return j_met, j_rec

def jaccard_similarity_matrices(models: Iterable):
    """The methods takes an Iterable of cobra models and returns a dictionary
    containing all pairwise jaccard similarities.

    Args:
        models (Iterable): List of cobra models

    Returns:
        (df): DataFrame of all jaccard similarities for metabolites indexed by the model ids.
        (df): DataFrame of all jaccard similarities for reaction indexed by the model ids.
        (df): DataFrame for resourece overlap indexed by the model ids.
    """    
    N = len(models)
    matrix_met = np.eye(N)
    matrix_rec = np.eye(N)
    matrix_ro = np.eye(N)
    index = [model.id for model in models]
    for i, model1 in enumerate(models):
        for j,model2 in zip(range(i+1,N),models[i+1:]):
            ji_met, ji_rec = jaccard_similarity(model1,model2)
            ji_ro = resource_overlap(model1,model2)
            matrix_met[i,j] = ji_met
            matrix_met[j,i] = ji_met
            matrix_rec[i,j] = ji_rec 
            matrix_rec[j,i] = ji_rec
            matrix_ro[i,j] = ji_ro 
            matrix_ro[j,i] = ji_ro

    df1 = pd.DataFrame(matrix_met, index=index, columns=index)
    df2 = pd.DataFrame(matrix_rec, index=index, columns=index)
    df3 = pd.DataFrame(matrix_ro, index=index, columns=index)

    return df1,df2,df3


def resource_overlap(model1: Model, model2:Model):
    """Coputes the resource overlap between two models

    Args:
        model1 (Model): Cobra model
        model2 (Model): Cobra model

    Returns:
        (float): Jacard index of resource overlap
    """
    in_ex1 = set([ex.id for ex in model1.exchanges if ex.lower_bound < 0])
    in_ex2 = set([ex.id for ex in model2.exchanges if ex.lower_bound < 0])

    common = in_ex1.intersection(in_ex2)
    union = in_ex1.union(in_ex2)

    return len(common)/len(union)


def write_out_common_metabolites(model1:Model, model2, prefix:str="common_reactions.csv"):
        # Check for correctness
    common_metabolits = [met for met in model1.metabolites if met in model2.metabolites]
    # Write csv
    df_dict = {"ID":[], "NAME":[], "FORMULA":[], "COMPARTMENT":[]}
    for met in common_metabolits:
        df_dict["ID"].append(met.id)
        df_dict["NAME"].append(met.name)
        df_dict["FORMULA"].append(met.formula)
        df_dict["COMPARTMENT"].append(met.compartment)
    df_common_met = pd.DataFrame(df_dict)
    df_common_met.to_csv(prefix)
    df_common_met.head()

def write_out_common_reactions(model1:Model, model2:Model, prefix:str="common_metabolites.csv"):
    common_reactions = [rec for rec in model1.reactions if rec in model2.reactions]
    # Write csv
    df_dict = {"ID":[], "NAME":[], "REACTION":[], "LOWER_BOUND":[], "UPPER_BOUND":[]}
    for rec in common_reactions:
        df_dict["ID"].append(rec.id)
        df_dict["NAME"].append(rec.name)
        df_dict["REACTION"].append(rec.reaction)
        df_dict["LOWER_BOUND"].append(rec.lower_bound)
        df_dict["UPPER_BOUND"].append(rec.upper_bound)

    def highlight_col(x):
        #copy df to new - original data are not changed
        df = x.copy()
        #mark exchange reactions by yellow color
        mask = x["ID"].str.contains("EX_")
        df.loc[mask, :] = 'background-color: yellow'
        df.loc[~mask,:] = 'background-color: white'
        return df   

    writer = pd.ExcelWriter(prefix)
    df_common_rec = pd.DataFrame(df_dict)
    df_common_rec.style.apply(highlight_col, axis=None).to_excel(writer)
    writer.save()
    writer.close()

