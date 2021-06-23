from typing import Iterable
from cobra.core import Model
import pandas as pd

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
    j_met = total_mets/common_mets
    j_rec = total_recs/common_recs
    return j_met, j_rec

def jaccard_similarity_matrix(models: Iterable):
    """The methods takes an Iterable of cobra models and returns a dictionary
    containing all pairwise jaccard similarities.

    Args:
        models (Iterable): List of cobra models

    Returns:
        (dict): Dictionary of all similarities indexed by the model ids!
    """    
    similarities = dict()
    for i, model1 in enumerate(models):
        for model2 in models[i:]:
            j_index = jaccard_similarity(model1,model2)
            similarities[(model1.id,model2.id)] = j_index
            similarities[(model2.id,model1.id)] = j_index
    return similarities


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

