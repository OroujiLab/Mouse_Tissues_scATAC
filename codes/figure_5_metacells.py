#!/usr/bin/env python

import os
import pandas as pd
import numpy as np

import scanpy as sc
import pyranges as pr
import warnings

import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns

import SEACells

data_dir = os.path.expanduser('/cluster/projects/epigenomics/Aminnn/scATAC/SEACells/Mouse_9Tissues/export/')


# Peaks data
from scipy.io import mmread
counts = mmread(data_dir + '/peak_counts/counts.mtx')


# Cell and peak information
cells = pd.read_csv(data_dir + '/peak_counts/cells.csv', index_col=0).iloc[:, 0]
peaks = pd.read_csv(data_dir + '/peak_counts/peaks.csv', index_col=0)
peaks.index = peaks['seqnames'] + ':' + peaks['start'].astype(str) + '-' + peaks['end'].astype(str)
#peaks.head(10)



ad = sc.AnnData(counts.T)
ad.obs_names = cells
ad.var_names = peaks.index
for col in peaks.columns:
  ad.var[col] = peaks[col]


ad.X = ad.X.tocsr()

ad.obsm['X_svd'] = pd.read_csv(data_dir + '/svd.csv', index_col=0).loc[ad.obs_names, : ].values

cell_meta = pd.read_csv(data_dir + '/cell_metadata.csv', index_col=0).loc[ad.obs_names, : ]
for col in cell_meta.columns:
  ad.obs[col] = cell_meta[col].values

# Gene scores
gene_scores = pd.read_csv(data_dir + '/gene_scores.csv', index_col=0).T


ad.obsm['GeneScores'] = gene_scores.loc[ad.obs_names, :].values
ad.uns['GeneScoresColums'] = gene_scores.columns.values

ad.write(data_dir + '/Mouse_All.h5ad')


# Filter to macrophages only
mac_types = ['Macrophages', 'Alveolar Macrophages', 'Microglial Cells']
ad_macro = ad[ad.obs['Clusters50'].isin(mac_types)].copy()

# Tissue composition heatmap
tissue_comp = ad_macro.obs.groupby(['SEACell', 'Sample']).size().unstack(fill_value=0)
tissue_comp_norm = tissue_comp.div(tissue_comp.sum(axis=1), axis=0)

# Sort tissues by macrophage count
mac_per_tissue = ad_macro.obs.groupby('Sample').size().sort_values(ascending=False)
tissue_order = mac_per_tissue.index.tolist()

# Sort columns by dominant tissue
dom_tissue = tissue_comp_norm[tissue_order].idxmax(axis=1)
col_order = []
for t in tissue_order:
    col_order.extend(dom_tissue[dom_tissue == t].index.tolist())

# Plot heatmap
plot_mat = tissue_comp_norm.loc[col_order, tissue_order].T
plt.figure(figsize=(12, 6))
sns.heatmap(plot_mat, cmap='YlOrRd', cbar_kws={'label': 'Fraction'})
plt.xlabel('Metacell', fontweight='bold')
plt.ylabel('Tissue', fontweight='bold')
plt.title('Macrophage Metacell Tissue Composition', fontweight='bold')
plt.tight_layout()
plt.savefig('macrophage_tissue_heatmap.pdf', dpi=300)
plt.show()

# Stacked bar plot
samples = ["Spleen", "Heart", "Pancreas", "Lung", "Small_Intestine", "Liver", "Colon", "Brain", "Kidney"]
colours = ["#F9B712", "#208A42", "#8A9FD1", "#FEE500", "#C06CAB", "#F47D2B", "#272E6A", "#D51F26", "#89288F"]

mac_counts = ad_macro.obs.groupby('Sample').size().reindex(samples, fill_value=0)

fig, ax = plt.subplots(figsize=(8, 1.5))
left = 0
for tissue, color in zip(samples, colours):
    width = mac_counts.loc[tissue]
    ax.barh(0, width, left=left, color=color, edgecolor='none')
    left += width

ax.set_xlim(0, left)
ax.set_yticks([])
ax.set_xlabel('Number of macrophage cells')
plt.tight_layout()
plt.savefig('macrophage_stacked_bar.pdf', dpi=300)
plt.show()


#FIBROBLAST

ad_fibro = ad[ad.obs['Clusters50'] == 'Fibroblasts'].copy()

tissue_comp_f = ad_fibro.obs.groupby(['SEACell', 'Sample']).size().unstack(fill_value=0)
tissue_comp_norm_f = tissue_comp_f.div(tissue_comp_f.sum(axis=1), axis=0)

fibro_per_tissue = ad_fibro.obs.groupby('Sample').size().sort_values(ascending=False)
tissue_order_f = fibro_per_tissue.index.tolist()

dom_tissue_f = tissue_comp_norm_f[tissue_order_f].idxmax(axis=1)
col_order_f = []
for t in tissue_order_f:
    col_order_f.extend(dom_tissue_f[dom_tissue_f == t].index.tolist())

plot_mat_f = tissue_comp_norm_f.loc[col_order_f, tissue_order_f].T
plt.figure(figsize=(12, 6))
sns.heatmap(plot_mat_f, cmap='YlOrRd', cbar_kws={'label': 'Fraction'})
plt.xlabel('Metacell', fontweight='bold')
plt.ylabel('Tissue', fontweight='bold')
plt.title('Fibroblast Metacell Tissue Composition', fontweight='bold')
plt.tight_layout()
plt.savefig('fibroblast_tissue_heatmap.pdf', dpi=300)
plt.show()

# Save metacell objects
ad_macro.write(data_dir + '/Metacells_Macro.h5ad')
ad_fibro.write(data_dir + '/Metacells_Fibro.h5ad')
