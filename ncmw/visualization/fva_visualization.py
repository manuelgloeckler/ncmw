
import matplotlib as plt 
import pandas as pd

def plot_full_fva(fva:pd.DataFrame, label_only_at_flux:int=5):
    """Plots all exchange fluxes, but only labels the one with large flux values.

    Args:
        fva (pd.DataFrame): FVA dataset
        label_only_at_flux (int, optional): For which value the label should be shown. Defaults to 5.
    """    
    fig = plt.figure(figsize=(20,5))
    non_zero_flux = fva["flux"].loc[fva["flux"] != 0]
    plt.plot(non_zero_flux.index, non_zero_flux, ".--")
    high_vals = fva.index[fva["flux"] > label_only_at_flux].tolist() + fva.index[fva["flux"] < -label_only_at_flux].tolist()
    _ = plt.xticks(high_vals, rotation=90)
    plt.tight_layout()
    plt.ylabel("Flux")
    fig.savefig("fva_all_fluxes.pdf")

def plot_medium_fva_range(model, growth_threshold = 1.0):
    fva_ex = model.summary(fva=growth_threshold)
    fva_ex = fva_ex.to_frame()
    medium = model.medium

    x = []
    y_fluxes = []
    y_min_flux = []
    y_max_flux = []
    for _, data in fva_ex.iterrows():
        id = str(data.metabolite[:-2])
        f = data.factor*data.flux 
        f_min = data.factor*data.minimum 
        f_max = data.factor*data.maximum

        if id in list(medium.keys()) and not (f_min == 0 and f_max == 0):
            x.append(id)
            y_fluxes.append(f)
            y_min_flux.append(f_min)
            y_max_flux.append(f_max)
    
    fig = plt.figure(figsize=(10,5), dpi= 160)
    plt.vlines(x=x, ymin=y_min_flux, ymax=y_max_flux, alpha=0.8, color="C0")

    # Plot FBA flux value as black dot
    plt.scatter(x=x, y= y_fluxes, color="black",s=5)

    #Plot maximum/minimum value in green
    for y_min, x1, y_max in zip(y_min_flux, x, y_max_flux):
        t = plt.text(x1,y_min+0.1, round(y_min, 2), horizontalalignment='left', color="green",
                    verticalalignment='bottom', size=4)
        t = plt.text(x1,y_max-1.5, round(y_max, 2), horizontalalignment='left', color="red",
                    verticalalignment='bottom', size=4)

    # Decorations    
    plt.title('Exchange flux range', fontdict={'size':18})
    plt.xticks(x,rotation=90)
    plt.grid(linestyle='--', alpha=0.5)
    plt.ylabel("$Flux$")
    plt.show()
    fig.savefig("ex_flux.pdf")
