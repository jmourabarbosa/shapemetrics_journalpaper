# %%
# -*- coding: utf-8 -*-
"""
@author: Amin
"""

import matplotlib.pyplot as plt
from sklearn import preprocessing
from matplotlib.patches import Rectangle

import collections
import numpy as np
import seaborn as sns
import pandas as pd

# %%
def plot_response_matrices(
        brain_regions, responses_dict, 
        acronyms_dict, eids_to_plot, dpi=300,
        save=False,file=None
    ):

    n_to_plot = len(eids_to_plot)
    fig, ax = plt.subplots(1, len(eids_to_plot), figsize=(3*n_to_plot, 6), dpi=dpi, sharey='all', sharex='all', squeeze=False)

    cmap = sns.color_palette('mako', as_cmap=True)
    scaler = preprocessing.MinMaxScaler()
    scale = lambda x: scaler.fit_transform(x)
    ax = ax.flatten()
    with sns.plotting_context('paper'):
        # sort columns of to_plot0 and to_plot1 by acronym 
        # ('VISa', 'CA1', 'DG', 'LP', 'PO')
        sorted_acronyms = ('VISa', 'CA1', 'DG', 'LP', 'PO')

        to_plot = collections.defaultdict(list)
        sorted_acronyms_dict = collections.defaultdict(list)
        for i in range(n_to_plot):
            eid_ = eids_to_plot[i]
            # sort columns by acronym
            for acronym in sorted_acronyms:
                idx = np.where(acronyms_dict[eid_] == acronym)[0]
                if len(idx) > 0:
                    to_plot[i].append(responses_dict[eid_].copy()[:, idx])
                    sorted_acronyms_dict[eid_].extend([acronym]*len(idx))
            to_plot[i] = scale(np.concatenate(to_plot[i], axis=1))

            for acr_idx, acronym in enumerate(sorted_acronyms_dict[eid_]):
                col = brain_regions.rgb[np.where(brain_regions.acronym == acronym)[0][0]]
                ax[i].add_patch(Rectangle((acr_idx, 0), -10, 10, color=col/255))

            ax[i].imshow(scale(to_plot[i]), aspect='auto', cmap=cmap)
            ax[i].set(xlabel='Unit', ylabel='Condition x Time bin', title=f'Experiment \n{eid_}',
                      )
            #make font smaller for title 
            ax[i].title.set_fontsize(8)
        fig.tight_layout
        sns.despine(bottom=True)
        

    if save:
        plt.savefig(file+'.png',format='png')
        plt.savefig(file+'.pdf',format='pdf')
        plt.close('all')
    else:
        plt.show()
        


# %%
def plot_var_explained(n_components,num_sessions,var_explained):
    with sns.plotting_context('paper'):
        fig, ax = plt.subplots(1, 1, figsize=(10, 5))
        ax.bar(
            range(num_sessions), var_explained.values(), 
            color='k'
        )
        sns.despine()
        fig.tight_layout()

    ax.set(
        title=f'Dim reduce session -> {n_components} components',
        xlabel='session',
        ylabel='variance explained'
    )


# %%
def plot_dist_mat(D,dist_mat,score_method):
    num_alpha = len(D)
    fig_distmat, ax = plt.subplots(
        1, len(D), figsize=(15, 3), sharex='all', sharey='all'
    )
    with sns.plotting_context('paper'):
        vmin, vmax = 0, np.max([np.max(v) for v in D.values()])

        cmap = sns.color_palette('mako', as_cmap=True)
        cmap.set_bad(color='k')
        for i, (alpha, dist_mat) in enumerate(D.items()):
            dist_mat_to_plot = dist_mat.copy()
            dist_mat_to_plot[np.diag_indices_from(dist_mat)] = np.nan
            sns.heatmap(
                dist_mat_to_plot, 
                square=True, cmap=cmap, cbar=False, vmax=vmax, ax=ax[i]
            )
            ax[i].set(title=f'alpha: {alpha:.2f}\nscore_method:{score_method}', xticks=[], yticks=[])
            ax[i].set(xlabel='session', ylabel='session')
        
        fig_distmat.tight_layout()


# %%
def plot_performance(
    filtered_eids,num_sessions,all_performance,
    psychometric,all_reaction_times,pca_2
):
    cols = sns.color_palette('colorblind', num_sessions)

    def plot_2d(pca_data, cols, ax=None):
        ax.scatter(pca_data[:, 0], pca_data[:, 1], c=cols, s=100)
        # add text of each area next to dot, offset by a little bit so it doesn't overlap
        ax.set(xlabel='PC 1', ylabel='PC 2',  title=f'Shape space, post-stim response| PCA({min_clusters}d) -> compute metric -> MDS(20d) -> PCA(2d)')

    df_plot = pd.DataFrame(
        {'eid': filtered_eids, 
        # 'performance': all_performance, 
        'performance': psychometric[:,0], 
        'reaction_time': all_reaction_times,
        'x': pca_2[:, 0], 
        'y': pca_2[:, 1], 
        'z': pca_2[:, 2]
        })

    # limit colors between the 4th min and max of all_performance
    cols = sns.color_palette('viridis', num_sessions)
    idx = np.argsort(all_performance)

    fig, ax= plt.subplots(1, 1, figsize=(10, 10))
    with sns.plotting_context():
        sns.scatterplot(data=df_plot, x='x', y='y', 
        hue='performance', palette='viridis',ax=ax)#, hue_norm=(.8, 1.))


# %%
# import matplotlib.pyplot as plt
# for i in range(len(Ds)):
#     plt.imshow(Ds[i])
#     plt.title(titles[i])
#     plt.show()

# # %%
# plt.scatter(
#     D_neural[~np.eye(D_neural.shape[0],dtype=bool)],
#     D_cc[~np.eye(D_cc.shape[0],dtype=bool)],
# )
# plt.show()

# import scipy
# print(scipy.stats.pearsonr(
#     D_neural[~np.eye(D_neural.shape[0],dtype=bool)],
#     D_cc[~np.eye(D_cc.shape[0],dtype=bool)],
# ))

# print(scipy.stats.spearmanr(
#     D_neural[~np.eye(D_neural.shape[0],dtype=bool)],
#     D_cc[~np.eye(D_cc.shape[0],dtype=bool)],
# ))

# %%
