import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import matplotlib

def plot_scaled_medium_growth(
    models,
    min_scale=0.1,
    max_scale=100,
    evaluations=1000,
    names: dict = dict(),
    cmap: str =None,
):
    """Plots the growth for scaling the medium by a variable factor K.
    
    
    
    Args:
        models: Cobra models
        min_scale: Minimal scale
        max_scale: Maximal scale
        evaluations: Number of scales between min_scale and max_scale
        names: Names dictionary
        cmap: Color map for colors.
    
    Returns:
        fig: Plot
    
    """
    fig = plt.figure(figsize=(8, 5))
    cmap = matplotlib.cm.get_cmap(cmap)

    # Title and axes
    plt.title("Scaled medium growth behaviour")
    plt.xlabel("Factor")
    plt.ylabel("Growth")
    K = np.linspace(min_scale, max_scale, evaluations)
    for i, model in enumerate(models):
        growths = []
        for k in K:
            with model as model:
                medium = model.medium
                for key in medium:
                    medium[key] *= k
                model.medium = medium
                growths += [model.slim_optimize()]
        plt.plot(K, growths, color=cmap(i/len(models)))
    model_name = []
    for model in models:
        if model.id in names:
            model_name.append(names[model.id])
        else:
            model_name.append(model.id.split("_")[0])
    plt.legend([n for n in model_name])
    return fig
