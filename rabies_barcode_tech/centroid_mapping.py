#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 1 09:00:00 2025
@author: ike
"""


import time
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

from pathlib import Path
from rich.table import Table
from rich.console import Console
from rich.progress import track
from rich.traceback import install


install(show_locals=True, width=120)


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


def format_str_index(
        index,
        replace: str = '_'
):
    """Format string index values to uppercase alphanumeric characters.

    Parameters
    ----------
    index : pd.Index | pd.Series
        Index of labels to format as strings.
    replace : str, optional
        Character used to replace non-alphanumeric characters in `index`.
        Defaults to '_'.

    Returns
    -------
    formatted_index : pd.Series
        `index` with non-alphanumeric characters replaced with `replace` and
        other characters converted to uppercase.
    """
    formatted_index = index.to_series().str.replace(
        r'[^0-9A-Za-z]', replace, regex=True).str.upper()
    return formatted_index


def normalize_counts(
        count_matrix: pd.DataFrame,
        counts: pd.Series = None,
        scalar: float = None,
        norm_seq_depth: bool = True,
        gene_x_cell: bool = True,
        log: bool = True
):
    """Normalize a count matrix to total UMI counts per cell.

    Parameters
    ----------
    count_matrix : pd.DataFrame
        One axis corresponds to identified genes. Other axis corresponds to
        sequenced cells.
    counts : pd.Series, optional
        Total UMI counts per cell_id (index). Useful if `count_matrix` is a
        filtered version of a larger DataFrame.
        Defaults to None, in which case counts are summed from `count_matrix`.
    scalar : float, optional
        Scale factor with which to adjust normalized counts.
        Defaults to None, in which case counts are not scaled.
    norm_seq_depth : bool, optional
    gene_x_cell : bool, optional
        If True, `count_matrix` is a [gene, cell] matrix. If False,
        `count_matrix` is a [cell, gene] matrix.
    log : bool, optional
        If True, take logarithm of counts + 1

    Returns
    -------
    new_counts : pd.DataFrame
        `count_matrix` normalized to total UMI counts per cell and scaled.

    Notes
    -----
    Genes counts for each cell are divided by the total UMIs for that cell,
    then multiplied by a scale factor.
    """
    new_counts = counts
    if norm_seq_depth:
        counts = counts if counts is not None else count_matrix.sum(
            numeric_only=True, axis=int(not gene_x_cell))
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
    type_assignment : pd.DataFrame
    """
    if write_assignment_df.is_file() and write_corr_scores_df.is_file():
        corr_scores = pd.read_csv(write_corr_scores_df, index_col=0)
        type_assignment = pd.read_csv(
            write_assignment_df, index_col=0).squeeze()
        return corr_scores, type_assignment

    # add cell type column
    ref_counts['type_updated'] = ref_metadata.loc[
        ref_counts.index, 'type_updated']

    # build centroids table of mean for each variable gene within a cluster
    ref_counts = ref_counts.groupby("type_updated").mean()  # cluster x gene

    # For each cell, calculate correlations across all genes for each centroid
    corr_scores = []
    for i in range(0, new_counts.shape[0], step):
        subset = new_counts.iloc[i:i+step]
        corr_scores += [pd.DataFrame(
            corr_2_matrix(subset.to_numpy(), ref_counts.to_numpy()),
            index=subset.index, columns=ref_counts.index)]

    corr_scores = pd.concat(corr_scores)
    corr_scores.to_csv(write_corr_scores_df)

    # Assign clusters based on highest value concordance (excluding any NaNs)
    type_assignment = corr_scores.idxmax(axis='columns').to_frame(
        name='celltype')
    type_assignment.to_csv(write_assignment_df)
    return corr_scores, type_assignment


def calculate_embeddings(
        ref_counts,
        new_counts,
        ref_umap,
        umap_cols,
        use_median: bool = True,
        knn: int = 10,
        write_assignment_dataframe: Path = None,
        step: int = 1_000
):
    """
    Parameters
    ----------
    ref_counts :
        sparse reference matrix
    new_counts :
        sparse sample matrix
    ref_umap :
        coordinates representing the embeddings for each cell from the
        reference atlas
    umap_cols :
        umap coord column names in `ref_umap`
    use_median : bool, optional
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
    if write_assignment_dataframe.is_file():
        return pd.read_csv(write_assignment_dataframe, index_col=0)

    # chunk over cells in new dataset
    assignment_positions = {}
    for i in range(0, new_counts.shape[0], step):
        sub_new = new_counts.iloc[i:i+step]
        corr_scores = []

        # chunk over cells in reference dataset
        for j in range(0, ref_counts.shape[0], step):
            sub_ref = ref_counts.iloc[j:j+step]

            # compute new_cells x ref_cells expression correlation matrix
            corr_scores += [pd.DataFrame(
                corr_2_matrix(sub_new.to_numpy(), sub_ref.to_numpy()),
                index=sub_new.index, columns=sub_ref.index)]

        # join batches by new cell id and iterate over new_cells
        for new_cell_id, row in pd.concat(corr_scores, axis=1).iterrows():
            # select knn reference cells based on largest correlations
            row = row.nlargest(knn)

            # extract umap coordinates for selected reference cells
            coords = ref_umap.loc[row.index, umap_cols].to_numpy()

            # take median/corr weighted average of reference UMAP coordinates
            coords = np.median(coords, axis=0) if use_median else np.average(
                coords, axis=0, weights=row.to_numpy())
            assignment_positions[new_cell_id] = coords

    assignment_positions = pd.DataFrame.from_dict(
        assignment_positions, orient="index", columns=list(umap_cols))

    # Save output to local machine
    assignment_positions.to_csv(write_assignment_dataframe)
    return assignment_positions


def __main__(
        ref_dir: Path,
        new_dir: Path,
        out_dir: Path,
        umap_cols: tuple = ("umap_1", "umap_2"),
        scale_factor: float = 10_000,
        norm_seq_depth: bool = True,
        use_median: bool = True,
        knn: int = 10,
        step: int = 1_000
):
    """"""
    """prepare reference dataset and output paths"""
    corr_dir = out_dir.joinpath("corr_scores")
    centroid_dir = out_dir.joinpath("centroid_assignments")
    coord_dir = out_dir.joinpath("umap_coordinates")
    plot_dir = out_dir.joinpath("umap_plots")
    for subdir in (corr_dir, centroid_dir, coord_dir, plot_dir):
        subdir.mkdir(exist_ok=True, parents=True)

    filtered_path = out_dir.joinpath("filtered_normed_wang_ref.csv")
    assignment_path = out_dir.joinpath("mapped_centroids.csv")

    # extract relevant genes and metadata from reference files
    ref_variable_genes = pd.read_csv(
        ref_dir.joinpath("ref_var_genes.csv"), usecols=['x'])['x']
    ref_metadata = pd.read_csv(
        ref_dir.joinpath("wang_metadata.csv"), index_col=0)

    # clean up cell ids to match those in raw data file
    ref_metadata.index = format_str_index(ref_metadata.index)
    ref_metadata['dataset_id'] = "wang et al"

    # only load relevant genes from reference dataset
    print(f'\nReference dataset file available? {filtered_path.is_file()}')
    if not filtered_path.is_file():
        ref_data = []
        ref_counts_path = ref_dir.joinpath("wang_ref.csv")
        counts = None
        start_time = time.time()
        for chunk in track(  # gene x cell matrix
                pd.read_csv(ref_counts_path, index_col=0, chunksize=step),
                description='    loading...', total=None):
            # clean cell ids, revert ref_data transformed with log1p
            chunk = np.expm1(chunk)
            chunk.columns = format_str_index(chunk.columns)

            # keep running total of UMI counts per cell, for all genes
            chunk_sum = chunk.sum(axis=0, numeric_only=True)
            counts = chunk_sum if counts is None else counts + chunk_sum

            # filter to only variable genes, cells with metadata
            chunk = chunk.loc[
                chunk.index.isin(ref_variable_genes),
                chunk.columns.intersection(ref_metadata.index)]
            ref_data += [chunk]

        # concat chunked raw data, normalize, save
        print(f'    load time: {time.time() - start_time:.2f}s')
        start_time = time.time()
        normalize_counts(
            pd.concat(ref_data), counts=counts, scalar=scale_factor,
            norm_seq_depth=norm_seq_depth, gene_x_cell=True, log=True).to_csv(
            filtered_path)
        print(f'    norm time: {time.time() - start_time:.2f}s')

    # load and convert to cell x gene matrix
    ref_data = pd.read_csv(filtered_path, index_col=0).T
    print(f'\nreference cell count: {ref_data.shape[0]}')
    print(f'variable genes count: {ref_data.shape[1]}')

    """run centroid_mapping, calculate_embeddings on each query dataset"""
    if not assignment_path.is_file():
        assignments = []
        logging = {}
        for query_path in track(
                list(new_dir.glob('[!.]*.csv')), description="mapping..."):
            """Prepare query dataset"""
            query_data = normalize_counts(
                pd.read_csv(query_path, index_col=0),  # index is gene
                scalar=scale_factor, norm_seq_depth=norm_seq_depth,
                gene_x_cell=True, log=True).T  # index is cell

            # Identify overlapping genes between reference and query datasets
            gg = ref_data.columns.intersection(query_data.columns)
            query_data = query_data.loc[:, gg]

            # map centroids
            centroid_time = time.time()
            corr, cell_type = centroid_mapping(
                ref_counts=ref_data.loc[:, gg],
                new_counts=query_data,
                ref_metadata=ref_metadata,
                write_assignment_df=centroid_dir.joinpath(query_path.name),
                write_corr_scores_df=corr_dir.joinpath(query_path.name))
            centroid_time = time.time() - centroid_time

            # get peak correlation for each cell
            corr = corr.max(axis=1).to_frame(name='high_score')
            corr[['dataset_id', 'cbc']] = corr.index.to_series().str.split(
                '_', n=1, expand=True)

            # add dataset id to cluster assignments
            name = query_path.stem
            corr['datasetid'] = name.split('_')[0] if '_' in name else name

            # extract umap coordinates
            embedding_time = time.time()
            coords = calculate_embeddings(
                ref_counts=ref_data.loc[:, gg],
                new_counts=query_data,
                ref_umap=ref_metadata,
                umap_cols=umap_cols,
                use_median=use_median,
                knn=knn,
                write_assignment_dataframe=coord_dir.joinpath(query_path.name))
            embedding_time = time.time() - embedding_time
            assignments += [pd.concat([corr, cell_type, coords], axis=1)]

            # look data from current loop
            logging[query_path.name] = (len(gg), centroid_time, embedding_time)

        # print logged metadata
        table = Table(title="Mapping Summary")
        table.add_column("dataset", style="bold cyan")
        table.add_column("n genes", style="magenta")
        table.add_column("t_centroid sec", style="green")
        table.add_column("t_embedding sec", style="green")
        for d, (g, t_c, t_e) in logging.items():
            table.add_row(d, str(g), f"{t_c:.2f}", f"{t_e:.2f}")

        Console().print(table)

        # index is cell id, cols high score, dataset_id, cbc, celltype, 2 umap
        assignments = pd.concat(assignments)
        mapping = {  # desired label: [dataset_ids]
            'SADB-19 cell': ['s1', 's2', 's3', 's4', 's5'],
            'CVS-N2c cell': ['c1', 'c2', 'c3', 'c4'],
            'CVS-N2c nuc': ['n1', 'n2', 'n3', 'n4'],
            'Uninfected cell': ['u1'], 'Uninfected nuc': ['u2']}
        for k, v in mapping.items():
            assignments.loc[assignments['dataset_id'].isin(v), 'cbc'] = k

        assignments.to_csv(out_dir.joinpath('filtered_assignments.csv'))

    # plot high score distribution
    assignments = pd.read_csv(assignment_path)
    plt.figure(figsize=(10, 5))
    plt.xticks(rotation=45)
    sns.violinplot(
        data=assignments, x='dataset_id', y='high_score', hue='dataset_id',
        legend=False, alpha=0.5)  #palette=dataset_palette)
    plt.show()

    # Drop all rows where high_score is <0.2
    total_cells = assignments.shape[0]
    assignments = assignments.query('high_score >= 0.2')
    print(f'Number of cells pre-thresholding: {total_cells}')
    print(f'Number of cells post-thresholding: {assignments.shape[0]}')
    print('Percentage of cells retained for analysis:')
    print(f'{100 * assignments.shape[0] / total_cells:.2f}%')

    # plot cell type umaps
    for cell_type in track(
            ref_metadata['type_updated'].unique(), description='plot...'):
        fig, ax = plt.subplots(figsize=(6, 6))
        sns.scatterplot(
            data=ref_metadata.query(f'type_updated == {cell_type}'),
            x=umap_cols[0], y=umap_cols[1], ax=ax, s=1, marker='.',
            legend=False, color="black")
        sns.scatterplot(
            data=assignments.query(f'type_updated == {cell_type}'),
            x=umap_cols[0], y=umap_cols[1], hue='dataset_id', ax=ax, s=1,
            marker='.', alpha=0.3, legend=True)
        plt.title(cell_type)
        plt.savefig(
            plot_dir.joinpath(f"{cell_type}.png"), dpi=300,
            bbox_inches="tight")
        plt.clf()


if __name__ == '__main__':
    REF_DIR = Path(
        "/Users/ikogbonna/Documents/Lab/Cadwell Lab/Data/barcoded_tech_data")
    if not REF_DIR.is_dir():
        REF_DIR = Path("/data/scratch/ike/barcoded_tech_data")

    __main__(
        ref_dir=REF_DIR.joinpath("wang_ref_atlas"),
        new_dir=REF_DIR.joinpath("sparse_matrices"),
        out_dir=REF_DIR.joinpath("outputs"))
