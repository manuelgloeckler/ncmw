import numpy as np
import pandas as pd


def compute_pairwise_growth_relation_per_weight(
    community_model, idx1, idx2, num_alphas=1000
):
    weights = community_model.weights
    N = len(weights)
    weight_mask1 = np.zeros(N)
    weight_mask1[idx1] = 1.0
    weight_mask2 = np.zeros(N)
    weight_mask2[idx2] = 1.0
    alpha = np.linspace(0, 1, num_alphas).reshape(-1, 1).repeat(N, 1)
    alpha_weights = alpha * weight_mask1 + (1 - alpha) * weight_mask2

    growth1 = np.zeros(num_alphas)
    growth2 = np.zeros(num_alphas)
    for i in range(num_alphas):
        community_model.weights = alpha_weights[i]
        _, single_growths = community_model.optimize()
        growth1[i] = single_growths[idx1]
        growth2[i] = single_growths[idx2]
    return alpha.flatten(), growth1, growth2


def compute_fair_weights(community_model):
    """Compute fair weights, fair in the sense that organisms with high single growth
    has a low weights...



    Args:
        community_model: Cobra Model

    Returns:
        array: Numpy array of weights

    """
    N = len(community_model.models)
    weights = np.array([community_model.single_optimize(i) for i in range(N)])

    return 1 - weights / weights.sum()


def compute_dominant_weights(community_model, high_val=0.9):
    N = len(community_model.models)
    others_w = (1 - high_val) / (N - 1)
    weights = []
    for i in range(N):
        base_weight = np.ones(N) * others_w
        base_weight[i] = high_val
        weights.append(base_weight)
    return weights


def compute_community_summary(community_model, weights):
    df_total = pd.DataFrame()
    for i in range(len(weights)):
        community_model.weights = weights[i]
        sol = community_model.optimize()
        df = pd.DataFrame(
            dict(zip([m.id for m in community_model.models], np.round(sol[1], 3))),
            index=[i],
        )
        df["Total"] = np.round(sol[0], 3)
        df["Weight"] = str(community_model.weights)
        df_total = df_total.append(df)
    return df_total
