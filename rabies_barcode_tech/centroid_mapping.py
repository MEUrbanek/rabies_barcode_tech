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
        write_assignment_df: Path,
        write_corr_scores_df: Path,
        step: int = 1_000
):
    """
    Parameters
    ----------
    ref_counts :
        sparse reference matrix. Normalized and filtered to relevant genes.
        cell x gene matrix.
    new_counts :
        sparse sample matrix. Normalized and filtered to relevant genes.
        cell x gene matrix.
    ref_metadata :
    write_assignment_df : Path
        where to save cell type assignment dataframe to on local machine
    write_corr_scores_df : Path
        where to save correlation coefficient matrix on local machine
    step : int, optional

    Returns
    -------
    corr_Scores : pd.DataFrame
    type_assignment : pd.Series
    """
    start_time = time.time()
    print(start_time)

    if write_assignment_df.is_file() and write_corr_scores_df.is_file():
        corr_scores = pd.read_csv(write_corr_scores_df, index_col=0)
        type_assignment = pd.read_csv(
            write_assignment_df, index_col=0).squeeze()
        return corr_scores, type_assignment

    # filter to cells in count matrix and reference metadata
    ref_counts = ref_counts.loc[
        ref_counts.index.intersection(ref_metadata.index)]
    ref_counts['type_updated'] = ref_metadata.loc[
        ref_counts.index, 'type_updated']

    # build centroids table of mean for each variable gene within a cluster
    ref_counts = ref_counts.groupby("type_updated").mean()  # cluster x gene

    # For each cell, calculate correlations across all genes for each centroid
    corr_scores = []
    for i in track(range(0, new_counts.shape[0], step), description="corr..."):
        subset = new_counts.iloc[i:i+step]
        corr_scores += [pd.DataFrame(
            corr_2_matrix(subset.to_numpy(), ref_counts.to_numpy()),
            index=subset.index, columns=ref_counts.index)]

    corr_scores = pd.concat(corr_scores)
    corr_scores.to_csv(write_corr_scores_df)

    # Assign clusters based on highest value concordance (excluding any NaNs)
    type_assignment = corr_scores.idxmax(axis='columns')
    type_assignment.to_csv(write_assignment_df)

    print('Done assigning cell types! :)')
    print(f'This took... {time.time() - start_time} seconds')
    return corr_scores, type_assignment


def calculate_embeddings(
        ref_counts,
        reference_genes,
        new_counts,
        new_genes,
        reference_atlas,
        select_median: bool = True,
        knn: int = 10,
        write_assignment_dataframe: Path = None,
        step: int = 1_000
):
    """
    Parameters
    ----------
    ref_counts :
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
    assignment_positions : pd.DataFrame
    """
    start_time = time.time()
    print(start_time)
    if write_assignment_dataframe.is_file():
        assignment_positions = pd.read_csv(
            write_assignment_dataframe, index_col=0)
        return assignment_positions

    assignment_positions = {}
    # chunk over cells in new dataset
    for i in track(range(0, new_counts.shape[0], step), description="corr..."):
        sub_new = new_counts.iloc[i:i+step]
        corr_scores = []

        # chunk over cells in reference dataset
        for j in range(0, ref_counts.shape[0], step):
            sub_ref = ref_counts.iloc[j:j+step]

            # compute new_cells x ref_cells expression correlation matrix
            corr_scores += [pd.DataFrame(
                corr_2_matrix(sub_new.to_numpy(), sub_ref.to_numpy()),
                index=sub_new.index, columns=sub_ref.index)]

        # join ref_cell batches and iterate over new_cells
        for idx, r in pd.concat(corr_scores, axis=1).iterrows():
            # select knn reference cells based on largest correlations
            r = r.nlargest(knn)
            coords = ref_metadata.loc[r.index, ["umap0", "umap1"]].to_numpy()

            # take median/corr weighted average of reference UMAP coordinates
            coords = np.median(coords, axis=0) if select_median else np.average(
                coords, axis=0, weights=r.to_numpy())
            assignment_positions[idx] = coords

    assignment_positions = pd.DataFrame.from_dict(
        assignment_positions, orient="index", columns=["umap0", "umap1"])

    # Save output to local machine
    print(assignment_positions)
    assignment_positions.to_csv(write_assignment_dataframe)

    print('Done finding UMAP coords! :)')
    print(f'This took... {time.time() - start_time} seconds')
    return assignment_positions


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
