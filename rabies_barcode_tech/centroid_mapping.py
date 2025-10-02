#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 1 09:00:00 2025
@author: ike
"""


import os
import math
import time
import numpy as np
import pandas as pd
import seaborn as sns
import multiprocess as mp
import matplotlib.pyplot as plt

from scipy import stats
from pathlib import Path

import scanpy as sc
import anndata as ad


pd.options.mode.chained_assignment = None

REF_DIR = Path(".../")


def corr2(a, b): #Build corr2 (Pearson correlation coefficient) matlab function
    a -= np.mean(a)
    b -= np.mean(b)

    r = (a*b).sum() / math.sqrt((a*a).sum() * (b*b).sum())
    return r


def centroid_mapping(referenceCounts,  # sparse reference matrix
                     referenceGenes,  # reference variable gene list
                     newCounts,  # sparse sample matrix
                     newGenes,  # sample variable gene list
                     referenceClusters,
                     # reference metadata column with cell type annotations
                     normSeqDepth=True,
                     # normalizes matrices by dividing by total counts per cell, multiplying by a scale factor, then taking the log of that value + 1
                     scaleFactor=10000,
                     # scale factor for normalization, scanpy and Seurat use 10,000
                     returnCmeans=False,
                     # whether to return matrix of correlation values for each cell x centroid combo
                     totalClusters=None,
                     # if you want to specify a different number of unique clusters
                     writeAssignmentDataframe=[],
                     # where to save cell type assignment dataframe to on local machine
                     writeCorrScoresDataframe=[]):  # where to save correlation coefficient matrix on local machine

    starttime = time.time()
    print(starttime)

    if normSeqDepth:
        # Genes counts for each cell are divided by the total UMIs for that cell, then multiplied by a scale factor.
        print('Normalizing sample matrix to sequencing depth per cell')
        newCounts.loc['total_counts'] = newCounts.sum(numeric_only=True,
                                                      axis=0)
        normalized = newCounts.divide(newCounts.loc['total_counts'], axis=1)
        normalized = normalized * scaleFactor
        # This is then natural-log transformed using log1p
        X = np.log1p(normalized)
    else:
        # Or process with just the raw matrix with log1p scaling
        print('No sequencing depth normalization')
        X = np.log1p(X)

    # Repeat the same normalization steps as above, but for reference dataset
    if normSeqDepth:
        print('Normalizing reference matrix to sequencing depth per cell')
        referenceCounts.loc['total_counts'] = referenceCounts.sum(
            numeric_only=True, axis=0)
        normalized = referenceCounts.divide(
            referenceCounts.loc['total_counts'], axis=1)
        normalized = normalized * scaleFactor
        T = np.log1p(normalized)
    else:
        print('No sequencing depth normalization')
        T = np.log1p(T)

    # Save the normalized copies of the matrices
    norm_ref = T.copy(deep=True)
    norm_sample = X.copy(deep=True)

    # Identify overlapping genes between reference dataset and query dataset
    gg = sorted(list(set(referenceGenes) & set(newGenes)))

    # Report back how many genes are in common+being used for query mapping
    print('Using a common set of ' + str(len(gg)) + ' genes.')
    print()

    # For both datasets, pull all rows corresponding to variable features
    T = T.loc[gg]
    X = X.loc[gg]

    # If clusters are already specified, set K to that number
    if totalClusters is not None:
        K = totalClusters
    # Otherwise, manually calculate the number of total clusters by adding 1 to total number of reference clusters
    else:
        types = referenceClusters.unique()
        types = pd.DataFrame(types)
        types.columns = ['type_assignment']
        types['cluster_number'] = types.index
        # Make dictionary to convert to integer cluster labels
        cell_type_dict = dict(
            zip(types['type_assignment'], types['cluster_number']))
        numerical_clusters = referenceClusters.replace(cell_type_dict)
        K = np.max(numerical_clusters) + 1

    # Assign clusters to cells in reference atlas
    T.columns = (ref_metadata['type_updated'])
    clusters = ref_metadata['type_updated'].unique()

    # Build centroids table by calculating the mean for each variable gene for all cells within a cluster
    centroids = pd.DataFrame(index=T.index)
    for celltype in clusters:
        corresponding_cells = T[celltype]
        means = corresponding_cells.mean(axis=1)
        means = pd.DataFrame(means)
        means.columns = [celltype]
        centroids = pd.concat([centroids, means], axis=1)

        # For each cell, calculate correlation coefficient across all variable genes for each centroid
    # print(centroids)
    input_cells = X.columns
    Cmeans = []
    for sample_cell in input_cells:
        individual_Cmeans = [sample_cell]
        for celltype in clusters:
            calculated = corr2(X[sample_cell], centroids[celltype])
            individual_Cmeans.append(calculated)
        Cmeans.append(individual_Cmeans)
    corr_scores = pd.DataFrame(Cmeans)
    corr_scores = corr_scores.set_index(0)
    corr_scores.columns = clusters

    # Assign clusters based on highest value concordance (excluding any NaNs)
    type_assignment = corr_scores.idxmax(axis=1)

    # Save cluster assignments and correlation coefficients to local machine
    type_assignment.to_csv(writeAssignmentDataframe)
    corr_scores.to_csv(writeCorrScoresDataframe)

    print('Done assigning cell types! :)')
    print()
    print('This took...')
    endtime = time.time()
    print(endtime - starttime)
    print('seconds!')

    # Return normalized matrices and type assignments
    if returnCmeans:
        return type_assignment, corr_scores, norm_ref, norm_sample
    else:
        return type_assignment, norm_ref, norm_sample


def __main__():
    # Load in reference variable genes, metadata, and sparse matrix
    ref_data = pd.read_csv(REF_DIR.joinpath('ref_matrix.csv'), index_col=0)
    ref_metadata = pd.read_csv(REF_DIR.joinpath('ref_metadata.csv'))
    ref_variable_genes = pd.read_csv(REF_DIR.joinpath('ref_var_genes.csv'))

    # Revert ref_data
    raw_ref_data = np.expm1(ref_data)
    print(raw_ref_data.shape)


