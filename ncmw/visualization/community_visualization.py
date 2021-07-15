import matplotlib.pyplot as plt
import numpy as np


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
