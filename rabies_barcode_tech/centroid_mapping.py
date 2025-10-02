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

from pathlib import Path


REF_DIR = Path(".../")


def corr2(
        a: np.ndarray,
        b: np.ndarray
):
    """Build corr2 (Pearson correlation coefficient) matlab function"""
    a = np.ravel(a).astype(float) - np.mean(a)
    b = np.ravel(b).astype(float) - np.mean(b)
    r = np.dot(a, b) / np.sqrt(np.dot(a, a) * np.dot(b, b))
    return r


def centroid_mapping(
        reference_counts,
        reference_genes,
        new_counts,
        new_genes,
        reference_clusters,
        norm_seq_depth: bool = True,
        scale_factor: int = 10_000,
        return_c_means: bool = False,
        total_clusters: int = None,
        write_assignment_df: list = None,
        write_corr_scores_df: list = None
):
    """
    Parameters
    ----------
    reference_counts :
        sparse reference matrix
    reference_genes :
        reference variable gene list
    new_counts :
        sparse sample matrix
    new_genes :
        sample variable gene list
    reference_clusters :
        reference metadata column with cell type annotations
    norm_seq_depth : bool, optional
        normalizes matrices by dividing by total counts per cell, multiplying
        by a scale factor, then taking the log of that value + 1
    scale_factor : int, optional
        scale factor for normalization, scanpy and Seurat use 10,000
    return_c_means : bool, optional
        whether to return matrix of correlation values for each cell x centroid
        combo
    total_clusters : int, optional
        if you want to specify a different number of unique clusters
    write_assignment_df : list, optional
        where to save cell type assignment dataframe to on local machine
    write_corr_scores_df : list, optional
        where to save correlation coefficient matrix on local machine

    Returns
    -------

    """
    start_time = time.time()
    print(start_time)
    write_assignment_df = (
        [] if write_assignment_df is None else write_assignment_df)
    write_corr_scores_df = (
        [] if write_corr_scores_df is None else write_corr_scores_df)

    """as written, assumes gene x cell matrix"""
    if norm_seq_depth:
        # Genes counts for each cell are divided by the total UMIs for that cell, then multiplied by a scale factor.
        print('Normalizing sample matrix to sequencing depth per cell')
        new_counts = new_counts.divide(
            new_counts.sum(numeric_only=True, axis=0), axis=1) * scale_factor

    X = np.log1p(new_counts)  # log1p normalization

    # Repeat the same normalization steps as above, but for reference dataset
    if norm_seq_depth:
        print('Normalizing reference matrix to sequencing depth per cell')
        reference_counts.loc['total_counts'] = reference_counts.sum(
            numeric_only=True, axis=0)
        normalized = reference_counts.divide(
            reference_counts.loc['total_counts'], axis=1)
        normalized = normalized * scale_factor
        T = np.log1p(normalized)
    else:
        print('No sequencing depth normalization')
        T = np.log1p(T)

    # Save the normalized copies of the matrices
    norm_ref = T.copy(deep=True)
    norm_sample = X.copy(deep=True)

    # Identify overlapping genes between reference dataset and query dataset
    gg = sorted(list(set(reference_genes) & set(new_genes)))

    # Report back how many genes are in common+being used for query mapping
    print('Using a common set of ' + str(len(gg)) + ' genes.')
    print()

    # For both datasets, pull all rows corresponding to variable features
    T = T.loc[gg]
    X = X.loc[gg]

    # If clusters are already specified, set K to that number
    if total_clusters is not None:
        K = total_clusters
    # Otherwise, manually calculate the number of total clusters by adding 1 to total number of reference clusters
    else:
        types = reference_clusters.unique()
        types = pd.DataFrame(types)
        types.columns = ['type_assignment']
        types['cluster_number'] = types.index
        # Make dictionary to convert to integer cluster labels
        cell_type_dict = dict(
            zip(types['type_assignment'], types['cluster_number']))
        numerical_clusters = reference_clusters.replace(cell_type_dict)
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
    type_assignment.to_csv(write_assignment_df)
    corr_scores.to_csv(write_corr_scores_df)

    print('Done assigning cell types! :)')
    print()
    print('This took...')
    endtime = time.time()
    print(endtime - start_time)
    print('seconds!')

    # Return normalized matrices and type assignments
    if return_c_means:
        return type_assignment, corr_scores, norm_ref, norm_sample

    return type_assignment, norm_ref, norm_sample


def calculate_embeddings(
        referenceCounts,  # reference matrix
        referenceGenes,  # reference variable genes list
        newCounts,  # query matrix
        newGenes,
        # if using, query matrix gene list--otherwise, can be reference variable genes list subset to genes measured in query dataset
        referenceAtlas,
        # coordinates representing the embeddings for each cell from the reference atlas
        scaleFactor=10000,  # normalization scale factor
        selectMedian=True,
        # whether to calculate the median or a weighted average for the query embeddings
        knn=10,
        # number of neighbors to consider when selecting median/weighted mean
        write_assignment_dataframe=[]
        # where to save assignment dataframe to on local machine
):

    starttime = time.time()
    print(starttime)

    # Genes counts for each cell are divided by the total UMIs for that cell, then multiplied by a scale factor.
    print('Normalizing sample matrix to sequencing depth per cell')
    newCounts.loc['total_counts'] = newCounts.sum(numeric_only=True,
                                                  axis=0)
    normalized = newCounts.divide(newCounts.loc['total_counts'], axis=1)
    normalized = normalized * scaleFactor
    # This is then natural-log transformed using log1p
    X = np.log1p(normalized)

    print('Normalizing reference matrix to sequencing depth per cell')
    referenceCounts.loc['total_counts'] = referenceCounts.sum(
        numeric_only=True, axis=0)
    normalized = referenceCounts.divide(
        referenceCounts.loc['total_counts'], axis=1)
    normalized = normalized * scaleFactor
    T = np.log1p(normalized)

    # Identify overlapping genes between reference dataset and query dataset (already did this in R)
    gg = sorted(list(set(referenceGenes) & set(newGenes)))

    # Report back how many genes are in common+being used for query mapping
    print('Using a common set of ' + str(len(gg)) + ' genes.')
    print()

    # For query dataset, pull all rows corresponding to variable features
    # If sparse matrix, send to dense
    T = T.loc[gg]
    X = X.loc[gg]

    ref_population = T.columns
    input_cells = X.columns
    assignmentPositions = pd.DataFrame()
    individual_assignment = []
    print('Beginning correlation calculations')
    if selectMedian == True:
        print('Assigning coordinates based on median')

    if selectMedian == False:
        print('Assigning coordinates based on weighted means')

    global build_correlations

    def build_correlations(sample_cell):
        global build_correlations
        # global assignmentPositions
        print('.', end='', flush=True)
        individual_corr = []
        for ref_cell in ref_population:
            calculated = corr2(X[sample_cell], T[ref_cell])
            individual_corr.append(calculated)
        ind = referenceAtlas[np.argpartition(individual_corr, -knn)][-knn:]

        if selectMedian == True:
            individual_assignment = np.median(ind, axis=0)

        if selectMedian == False:
            res = np.array(individual_corr)
            weights = res[np.argpartition(res, -knn)][-knn:]
            weights = np.transpose(weights)
            individual_assignment = np.average(ind, axis=0,
                                               weights=weights)

        temp = pd.DataFrame(individual_assignment)
        temp.columns = [sample_cell]
        return temp

    #
    pool = mp.Pool(32)
    print('Calling workers')  # Create a multiprocessing Pool
    assignmentPositions = pd.concat(
        pool.map(build_correlations, input_cells), axis=1)

    # Save output to local machine
    print(assignmentPositions)
    assignmentPositions = assignmentPositions.T
    assignmentPositions.to_csv(write_assignment_dataframe)

    print('Done finding UMAP coords! :)')
    print()
    print('This took...')
    endtime = time.time()
    print(endtime - starttime)
    print('seconds!')

    return assignmentPositions


def __main__():
    # Load in reference variable genes, metadata, and sparse matrix
    ref_data = pd.read_csv(REF_DIR.joinpath('ref_matrix.csv'), index_col=0)
    ref_metadata = pd.read_csv(REF_DIR.joinpath('ref_metadata.csv'))
    ref_variable_genes = pd.read_csv(REF_DIR.joinpath('ref_var_genes.csv'))

    # Revert ref_data
    raw_ref_data = np.expm1(ref_data)
    print(raw_ref_data.shape)


