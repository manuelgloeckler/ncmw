import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib_venn import venn2
import matplotlib


def jacard_index_similarity_heatmap(ji_met, ji_rec, ji_ro, names: dict = dict()):
    """This plots the similarity between the community members as heatmap.

    Args:
        ji_met (df): Jaccard index matrix metabolites
        ji_rec (df): Jaccard index matrix reactions
        ji_ro (df): Jaccard index matrix resource overlap
        figsize (tuple, optional): [description]. Defaults to (10,5).

    Returns:
        [fig]: Figure
    """
    N = len(ji_met)
    fig, axes = plt.subplots(1, 3, sharey=True, figsize=(N * 3 + 2, N + 2))
    titles = ["Metabolites", "Reactions", "Exchanges (RO)"]
    cbar_ax = fig.add_axes([1.0, 0.4, 0.03, 0.4])
    model_name = []
    for id in ji_met.index:
        if id in names:
            model_name.append(names[id])
        else:
            model_name.append(id)

    ji_met.index = model_name
    ji_rec.index = model_name
    ji_ro.index = model_name

    for i, df in enumerate([ji_met, ji_rec, ji_ro]):
        axes[i].set_title(titles[i])
        sns.heatmap(
            df, annot=True, cmap="RdBu", vmin=0.0, vmax=1.0, ax=axes[i], cbar_ax=cbar_ax
        )
        axes[i].set_yticklabels([str(id).split("_")[0] for id in df.index], rotation=0)
        axes[i].set_xticklabels([str(id).split("_")[0] for id in df.index], rotation=90)
    plt.tight_layout()
    return fig


def uptake_sekretion_venn_diagrams(
    models, uptakes, sekretions, names: dict = dict(), cmap=None
):
    N = len(models)
    cmap = matplotlib.cm.get_cmap(cmap)
    titles = []
    for model in models:
        if model.id in names:
            titles.append(names[model.id])
        else:
            titles.append(model.id.split("_")[0])
    uptakes = [set(up) for up in uptakes]
    sekretions = [set(sek) for sek in sekretions]

    assert len(uptakes) == N, "We need for each model a uptake and sekretion."
    assert len(sekretions) == N, "We need for each model a uptake and sekretion."

    fig, axes = plt.subplots(N, N, figsize=(2 * 5, N * 2))
    for i in range(N):
        for j in range(i + 1, N):
            axes[i, j].set_title("Uptake")
            axes[i, j].xaxis.set_visible(False)
            venn2(
                [uptakes[i], uptakes[j]],
                ax=axes[i, j],
                set_labels=[titles[i], titles[j]],
                set_colors=[cmap(i / N), cmap(j / N)],
            )

    for i in range(N):
        axes[i, i].axis("off")

    for j in range(N):
        for i in range(j + 1, N):
            axes[i, j].set_title("Sekretion")
            venn2(
                [sekretions[i], sekretions[j]],
                ax=axes[i, j],
                set_labels=[titles[i], titles[j]],
                set_colors=[cmap(i / N), cmap(j / N)],
            )

    for i in range(N):
        axes[-1, i].xaxis.set_visible(True)
        axes[-1, i].set_xlabel(models[i].id)

    fig.tight_layout()
    return fig
