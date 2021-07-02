import matplotlib.pyplot as plt 
import seaborn as sns

def jacard_index_similarity_heatmap(ji_met, ji_rec, ji_ro, figsize=(10,5)):
    """This plots the similarity between the community members as heatmap.

    Args:
        ji_met (df): Jaccard index matrix metabolites
        ji_rec (df): Jaccard index matrix reactions
        ji_ro (df): Jaccard index matrix resource overlap
        figsize (tuple, optional): [description]. Defaults to (10,5).

    Returns:
        [fig]: Figure
    """
    fig, axes = plt.subplots(1,3, sharey=True, figsize=figsize)
    titles = ["Metabolites", "Reactions", "Exchanges (RO)"]
    cbar_ax = fig.add_axes([1.0, .4, .03, .4])
    for i,df in enumerate([ji_met,ji_rec,ji_ro]):
        axes[i].set_title(titles[i])
        sns.heatmap(df, annot=True,cmap="RdBu", vmin=0.0, vmax=1.0, ax=axes[i],cbar_ax=cbar_ax)
        axes[i].set_yticklabels(df.index,rotation=0)
        axes[i].set_xticklabels(df.index,rotation=90) 
    plt.tight_layout()
    return fig