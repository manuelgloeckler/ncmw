import matplotlib.pyplot as plt 
import pandas as pd
import numpy as np


def plot_scaled_medium_growth(models, min_scale=0.1, max_scale=100, evaluations=1000):
    fig = plt.figure(figsize=(8,5))
    # Title and axes
    plt.title("Scaled medium growth behaviour")
    plt.xlabel("Factor")
    plt.ylabel("Growth")
    K = np.linspace(min_scale, max_scale, evaluations)
    for i,model in enumerate(models):
        growths = []
        for k in K:
            with model as model:
                medium = model.medium
                for key in medium:
                    medium[key] *= k
                model.medium = medium
                growths += [model.slim_optimize()]
        plt.plot(K, growths)
    plt.legend([model.id for model in models])
    return fig
        
