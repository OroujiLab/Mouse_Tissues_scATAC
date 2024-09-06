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
