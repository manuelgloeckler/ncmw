import numpy as np
import pandas as pd

import networkx as nx

import torch
from sbi.inference import SNLE


def compute_community_interaction_graph(model, df):
    df_help = df[df.columns[:-1]]
    medium_col = np.zeros(len(df_help))
    for i, ex in enumerate(df_help.index):
        if ex in model.medium:
            medium_col[i] = model.medium[ex]
        else:
            medium_col[i] = 0.0
    # Set output to zero
    help_array = df_help.to_numpy()
    help_array[help_array >= 0] = 0
    # Non medium associated inputs -> interactions !
    species_interaction = (help_array.sum(1) + medium_col) < -1e-6
    species_interaction
    df = df[species_interaction]
    df = df.drop("Shuttle Reaction", 1)

    # Build interaction graph
    G = nx.DiGraph()
    # Build nodes
    for i, col in enumerate(df.columns):
        names = [n[3:-2] + f"_{i}" for n in df[df[col] != 0].index]
        G.add_nodes_from(names)

    for n in G.nodes:
        i = int(n[-1])
        G.nodes[n]["class"] = model.models[i].id.split("_")[0]
    return G, df


def community_weight_posterior(model):
    N = len(model.models)

    def simulator(thetas):
        xs = []
        for weight in thetas:
            model.weights = weight.numpy()
            growths = torch.tensor(model.optimize()[-1])
            xs.append(growths)
        return torch.vstack(xs).float()

    prior = torch.distributions.Dirichlet(torch.ones(N))
    prior.set_default_validate_args(False)

    thetas = prior.sample((N * 1000,))
    xs = simulator(thetas)

    inf = SNLE(prior)
    density_estimator = inf.append_simulations(thetas, xs).train()
    posterior = inf.build_posterior(
        density_estimator,
        sample_with="vi",
        vi_parameters={"flow": "spline_autoregressive", "bound": 15, "num_bins": 15},
    )
    return posterior


def compute_pairwise_growth_relation_per_weight(
    community_model, idx1, idx2, num_alphas=1000
):
    weights = community_model.weights
    N = len(weights)
    weight_mask1 = np.zeros(N)
    weight_mask1[idx1] = 1.0
    weight_mask2 = np.zeros(N)
    weight_mask2[idx2] = 1.0
    alpha = np.linspace(0, 1, num_alphas)
    alphas = alpha.reshape(-1, 1).repeat(N, 1)
    alpha_weights = alphas * weight_mask1 + (1 - alphas) * weight_mask2

    growth1 = np.zeros(num_alphas)
    growth2 = np.zeros(num_alphas)
    for i in range(num_alphas):
        community_model.weights = alpha_weights[i]
        _, single_growths = community_model.optimize()
        growth1[i] = single_growths[idx1]
        growth2[i] = single_growths[idx2]
    return alpha, growth1, growth2


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
