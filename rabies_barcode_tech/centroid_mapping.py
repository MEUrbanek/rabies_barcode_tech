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
from rich.progress import track


def corr2(
        a: np.ndarray,
        b: np.ndarray
):
    """Build corr2 (Pearson correlation coefficient) matlab function"""
    a = np.ravel(a).astype(float) - np.mean(a)
    b = np.ravel(b).astype(float) - np.mean(b)
    r = np.dot(a, b) / np.sqrt(np.dot(a, a) * np.dot(b, b))
    return r


def corr_2_matrix(
        a: np.ndarray,
        b: np.ndarray
):
    """corr2 func extended to 2d input matrices.

    Parameters
    ----------
    a : np.ndarray
        2d array to correlate. Has shape (samples_a, observations)
    b : np.ndarray
        Another 2d array to correlate. Has shape (samples_b, observations)

    Returns
    -------
    corr_matrix
        2d array with shape (samples_a, samples_b). The value at [i, j] is the
        Pearson correlation between the i_th sample in array `a` and j_th
        sample in array `b`.

    """
    # Center columns (genes)
    mat_a = a - a.mean(axis=1, keepdims=True)  # cell x gene matrix
    mat_b = b - b.mean(axis=1, keepdims=True)  # cell x gene matrix

    # Compute numerator: X^T @ Y (dot product of centered columns)
    numerator = mat_a @ mat_b.T  # shape: (n_cells_a, n_cells_b)

    # Compute denominator: column-wise norms, shape = (n_cells1, n_cells2)
    denominator = np.outer(
        np.linalg.norm(mat_a, axis=1), np.linalg.norm(mat_b, axis=1))

    # Pearson correlation matrix
    corr_matrix = numerator / denominator
    return corr_matrix


def intersect_index(
        *args
):
    """Filter dataframes based on index intersection.

    Parameters
    ----------
    *args : pd.DataFrame
        Dataframes to filter.

    Returns
    -------
    dfs : list[pd.DataFrame]
        Dataframes filtered to intersection of index values.
    g : pd.Index
        Intersection of all indices for dataframes.
    """
    g = args[0].index
    for df in args[1:]:
        g = g.intersection(df.index)

    g = g.sort_index(ascending=True)
    dfs = [df.loc[g] for df in args]
    return dfs, g


def normalize_counts(
        count_matrix: pd.DataFrame,
        scalar: float = None,
        gene_x_cell: bool = True,
        log: bool = True
):
    """Normalize a count matrix to total UMI counts per cell.

    Parameters
    ----------
    count_matrix : pd.DataFrame
        One axis corresponds to identified genes. Other axis corresponds to
        sequenced cells.
    scalar : float, optional
        Scale factor with which to adjust normalized counts.
        Defaults to None, in which case counts are not scaled.
    gene_x_cell : bool, optional
        If True, `count_matrix` is a [gene, cell] matrix. If False,
        `count_matrix` is a [cell, gene] matrix.
    log : bool, optional
        If True, take logarithm of counts + 1

    Returns
    -------
    new_count_matrix : pd.DataFrame
        `count_matrix` normalized to total UMI counts per cell and scaled.

    Notes
    -----
    Genes counts for each cell are divided by the total UMIs for that cell,
    then multiplied by a scale factor.
    """
    counts = count_matrix.sum(numeric_only=True, axis=int(not gene_x_cell))
    new_counts = count_matrix.divide(counts, axis=int(gene_x_cell))
    new_counts = new_counts if scalar is None else new_counts * scalar
    new_counts = np.log1p(new_counts) if log else new_counts
    return new_counts


def centroid_mapping(
        ref_counts,
        new_counts,
        ref_metadata,
        write_assignment_df: Path = None,
        write_corr_scores_df: Path = None,
        step: int = 1_000
):
    """
    Parameters
    ----------
    ref_counts :
        sparse reference matrix. Normalized and filtered to relevant genes.
    new_counts :
        sparse sample matrix. Normalized and filtered to relevant genes.
    ref_metadata :
    write_assignment_df : Path
        where to save cell type assignment dataframe to on local machine
    write_corr_scores_df : Path
        where to save correlation coefficient matrix on local machine
    step : int, optional

    Returns
    -------
    type_assignment : pd.Series
    corr_Scores : pd.DataFrame
    """
    start_time = time.time()
    print(start_time)

    # filter to cells in count matrix and reference metadata
    ref_counts = ref_counts.loc[
        ref_counts.index.intersection(ref_metadata.index)]
    ref_counts['type_updated'] = ref_metadata.loc[
        ref_counts.index, 'type_updated']

    # build centroids table of mean for each variable gene within a cluster
    ref_counts = ref_counts.groupby("type_updated").mean()  # cluster x gene

    # For each cell, calculate correlations across all genes for each centroid
    corr_scores = []
    for c in track(range(0, new_counts.shape[0], step), description="corr..."):
        subset = new_counts.iloc[c:c+step]
        corr_scores += [pd.DataFrame(
            corr_2_matrix(subset.to_numpy(), ref_counts.to_numpy()),
            columns=ref_counts.index, index=subset.index)]

    corr_scores = pd.concat(corr_scores)
    if write_corr_scores_df is not None:
        corr_scores.to_csv(write_corr_scores_df)

    # Assign clusters based on highest value concordance (excluding any NaNs)
    type_assignment = corr_scores.idxmax(axis='columns')
    if write_assignment_df is not None:
        type_assignment.to_csv(write_assignment_df)

    print('Done assigning cell types! :)')
    print(f'This took... {time.time() - start_time} seconds')
    return type_assignment, corr_scores


def calculate_embeddings(
        reference_counts,
        reference_genes,
        new_counts,
        new_genes,
        reference_atlas,
        scale_factor: float = 10_000,
        select_median: bool = True,
        knn: int = 10,
        write_assignment_dataframe: Path = None
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
        if using, query matrix gene list--otherwise, can be reference variable
        genes list subset to genes measured in query dataset
    reference_atlas :
        coordinates representing the embeddings for each cell from the
        reference atlas
    scale_factor : int, optional
        scale factor for normalization, scanpy and Seurat use 10,000
    select_median : bool, optional
        whether to calculate the median or a weighted average for the query
        embeddings
    knn : int, optional
        number of neighbors to consider when selecting median/weighted mean
    write_assignment_dataframe : Path, optional
        where to save assignment dataframe to on local machine
    step : int, optional

    Returns
    -------

    """
    starttime = time.time()
    print(starttime)

    # # Identify overlapping genes between reference dataset and query dataset
    # gg = sorted(list(set(reference_genes).intersection(new_genes)))
    # print(f'Using a common set of {len(gg)} genes.')
    #
    # # For both datasets, pull all rows corresponding to variable features
    # new_counts = new_counts.loc[gg]
    # reference_counts = reference_counts.loc[gg]
    #
    # print('Normalizing sample matrix to sequencing depth per cell')
    # print('Normalizing reference matrix to sequencing depth per cell')
    #
    # # This is then natural-log transformed using log1p
    # X = np.log1p(normalize_counts(new_counts, scale_factor))
    # T = np.log1p(normalize_counts(reference_counts, scale_factor))


    ref_population = T.columns
    input_cells = X.columns
    assignmentPositions = pd.DataFrame()
    individual_assignment = []
    print('Beginning correlation calculations')
    if select_median == True:
        print('Assigning coordinates based on median')

    if select_median == False:
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
        ind = reference_atlas[np.argpartition(individual_corr, -knn)][-knn:]

        if select_median == True:
            individual_assignment = np.median(ind, axis=0)

        if select_median == False:
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


def __main__(
        ref_data_path: Path,
        ref_meta_path: Path,
        ref_gene_path: Path,
        scale_factor: float = 10_000,
        norm_seq_depth: bool = True,
        step: int = 1_000
):
    # extract relevant genes from reference file
    ref_variable_genes = pd.read_csv(ref_gene_path, usecols=['x'])['x']

    # only load relevant genes from reference dataset
    ref_data = []
    for chunk in track(
            pd.read_csv(ref_data_path, index_col=0, chunksize=step),
            description="loading..."):
        ref_data += [chunk.loc[chunk.index.isin(ref_variable_genes)]]

    # Revert ref_data
    ref_data = np.expm1(pd.concat(ref_data))
    print(ref_data.shape)
    new_data = pd.read_csv("", index_col=0)  # index is gene label

    # Load in reference variable genes, metadata, and sparse matrix
    ref_metadata = pd.read_csv(ref_meta_path, index_col=0)

    # Identify overlapping genes between reference dataset and query dataset
    (ref_data, new_data), g = intersect_index(ref_data, new_data)
    print(f'Using a common set of {len(g)} genes.')

    # normalize and convert to cell x gene matrices
    print('Normalizing sample matrix to sequencing depth per cell')
    ref_data = normalize_counts(
        ref_data, scalar=scale_factor, gene_x_cell=True, log=True).T
    print('Normalizing reference matrix to sequencing depth per cell')
    new_data = normalize_counts(
        new_data, scalar=scale_factor, gene_x_cell=True, log=True).T


if __name__ == '__main__':
    REF_DIR = Path(
        "/Users/ikogbonna/Documents/Lab/Cadwell Lab/Data/barcoded_tech_data/wang_ref_atlas")
    __main__()
