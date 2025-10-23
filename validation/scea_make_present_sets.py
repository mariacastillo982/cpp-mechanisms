#!/usr/bin/env python3
"""
Build per-cluster 'present gene' sets from Single Cell Expression Atlas (SCEA)
normalised counts + clustering file.

Inputs:
    - (A) MatrixMarket counts directory or archive containing:
            matrix.mtx[.gz], barcodes.tsv[.gz], genes.tsv[.gz]
        OR
    - (B) A wide TSV counts matrix (rows=genes (Ensembl IDs), cols=cells)

    - SCEA clustering TSV (maps cell_id -> cluster label/ID, possibly multiple resolutions).

    Outputs:
    - present_genes_by_cluster.csv  (cluster, n_cells, n_present, present_ensembl_ids)
    - cluster_gene_stats.csv        (cluster, gene_id, detection_fraction, mean_cpm, present)

    Presence rule (tunable):
    present == (detection_fraction >= --min-detect-frac) OR (mean_cpm >= --min-cpm)

    Example:
    python3 scea_make_present_sets.py \
        --counts-mtx /path/to/normalised_counts/ \
        --clustering  /path/to/clustering.tsv \
        --cluster-col cluster \
        --outdir out_eval/ \
        --min-detect-frac 0.1 \
        --min-cpm 1.0
"""
from __future__ import annotations
import argparse
import gzip
import os
from pathlib import Path
from typing import Optional, Tuple

import numpy as np
import pandas as pd
from scipy import sparse

# ---------------- I/O helpers ---------------- #

def _open_maybe_gzip(path: Path, mode: str = "rt"):
    if str(path).endswith(".gz"):
        return gzip.open(path, mode)
    return open(path, mode)

def _read_matrix_market(mtx_dir: Path) -> Tuple[sparse.csr_matrix, pd.Index, pd.Index]:
    """
    Read a MatrixMarket bundle (matrix.mtx[.gz], barcodes.tsv[.gz], genes.tsv[.gz]).
    Returns (X_csr (cells x genes) in CPM units if already normalised by SCEA,
                cell_ids (index),
                gene_ids (index of Ensembl IDs))
    """
    # allow either files or gz
    mtx = None
    for fname in ("matrix.mtx.gz", "matrix.mtx"):
        f = mtx_dir / fname
        if f.exists():
            mtx = f
            break
    if mtx is None:
        raise FileNotFoundError("Could not find matrix.mtx[.gz] in {}".format(mtx_dir))

    # genes file may be 'genes.tsv' or 'features.tsv'
    genes_file = mtx_dir / "genes.mtx_rows"
    # barcodes
    barcodes_file = mtx_dir / "barcodes.mtx_cols"
    
    # read MM
    from scipy.io import mmread
    X = mmread(str(mtx)).tocsr().astype(np.float32)  # cells x genes or genes x cells?
    # Heuristic: SCEA normalized counts often come as genes x cells (MatrixMarket standard).
    # We'll transpose to cells x genes if needed by comparing barcodes length.
    with _open_maybe_gzip(barcodes_file, "rt") as fh:
        barcodes = [line.strip().split("\t")[0] for line in fh if line.strip()]
    # genes file: first column is Ensembl ID
    with _open_maybe_gzip(genes_file, "rt") as fh:
        genes = [line.strip().split("\t")[0] for line in fh if line.strip()]

    # orient
    if X.shape[0] == len(genes) and X.shape[1] == len(barcodes):
        # genes x cells -> transpose
        X = X.T.tocsr()
    elif X.shape[0] == len(barcodes) and X.shape[1] == len(genes):
        # already cells x genes
        pass
    else:
        raise ValueError(f"Matrix shape {X.shape} doesn't match barcodes {len(barcodes)} and genes {len(genes)}")

    return X, pd.Index(barcodes, name="cell_id"), pd.Index(genes, name="gene_id")

def _read_counts_wide(tsv_path: Path) -> Tuple[pd.DataFrame, pd.Index, pd.Index]:
    """
    Read a wide counts TSV (rows=genes, cols=cells), returns (df_T (cells x genes), cell_ids, gene_ids).
    """
    df = pd.read_csv(tsv_path, sep=None, engine="python", index_col=0)
    # index are genes, columns are cells
    gene_ids = pd.Index(df.index.astype(str), name="gene_id")
    cell_ids = pd.Index(df.columns.astype(str), name="cell_id")
    # transpose to cells x genes
    return df.T, cell_ids, gene_ids

def _load_counts(args) -> Tuple[sparse.csr_matrix, pd.Index, pd.Index]:
    if args.counts_mtx:
        X, cells, genes = _read_matrix_market(Path(args.counts_mtx))
        # Normalised counts from SCEA are CPM (not log), so we can use as-is. :contentReference[oaicite:2]{index=2}
        return X, cells, genes
    elif args.counts_tsv:
        df_T, cells, genes = _read_counts_wide(Path(args.counts_tsv))
        # convert to sparse for efficiency
        X = sparse.csr_matrix(df_T.values.astype(np.float32))
        return X, cells, genes
    else:
        raise ValueError("Provide either --counts-mtx (MM directory) or --counts-tsv (wide TSV).")

# ---------------- Clustering ---------------- #

def _load_clustering(clustering_path: Path, cluster_col: str) -> pd.Series:
    """
    Read SCEA clustering TSV and return pd.Series mapping cell_id -> cluster_label.
    """
    cl = pd.read_csv(clustering_path, sep=None, engine="python", dtype=str)
    # Columns vary by experiment; typical include: 'cell', 'cluster', 'louvain', etc. SCEA lists this under Downloads. 
    # Try common id column names
    for cid in ("cell", "cell_id", "barcode", "run"):
        if cid in cl.columns:
            cell_col = cid
            break
    else:
        # take first column
        cell_col = cl.columns[0]
    if cluster_col not in cl.columns:
        raise ValueError(f"Cluster column '{cluster_col}' not found in clustering file. Available: {list(cl.columns)}")
    s = cl.set_index(cell_col)[cluster_col].astype(str)
    return s

# ---------------- Core logic ---------------- #

def main():
    ap = argparse.ArgumentParser(description="Make per-cluster present-gene sets from SCEA normalised counts + clustering")
    g = ap.add_mutually_exclusive_group(required=True)
    g.add_argument("--counts-mtx", help="Path to MatrixMarket counts directory (contains matrix.mtx[.gz], barcodes.tsv[.gz], genes.tsv[.gz])")
    g.add_argument("--counts-tsv", help="Path to WIDE TSV counts (rows=Ensembl gene IDs, cols=cell IDs)")
    ap.add_argument("--clustering", required=True, help="SCEA clustering TSV")
    ap.add_argument("--cluster-col", default="cluster", help="Column in clustering TSV with the cluster label to use (default: 'cluster')")
    ap.add_argument("--outdir", required=True, help="Output directory")
    ap.add_argument("--min-detect-frac", type=float, default=0.10, help="Min detection fraction to call present (default: 0.10)")
    ap.add_argument("--min-cpm", type=float, default=1.0, help="Min mean CPM to call present (default: 1.0)")
    ap.add_argument("--accept-prefix", action="append", default=["ENSG"], help="Accept Ensembl ID prefixes (repeatable). Default: ENSG")
    args = ap.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    # Load
    X, cell_ids, gene_ids = _load_counts(args)
    
    clustering = _load_clustering(Path(args.clustering), args.cluster_col)

    # Align cell order across counts and clustering
    # Keep only intersecting cells
    common = cell_ids.intersection(clustering.index)
    if common.empty:
        raise ValueError("No overlapping cell IDs between counts and clustering.")
    # reindex
    keep_mask = pd.Index(cell_ids).isin(common)
    X = X[keep_mask, :]
    cell_ids = cell_ids[keep_mask]
    clustering = clustering.loc[cell_ids]  # ordered to match cells

    # Filter genes by Ensembl prefix if requested
    if args.accept_prefix:
        pref_tup = tuple(p.upper() for p in args.accept_prefix)
        keep_genes = [i for i, g in enumerate(gene_ids) if str(g).upper().startswith(pref_tup)]
        if not keep_genes:
            raise ValueError("No genes matched the requested Ensembl prefixes.")
        X = X[:, keep_genes]
        gene_ids = gene_ids[keep_genes]

    clusters = pd.Index(sorted(clustering.unique()))
    n_cells_total = len(cell_ids)

    # Precompute boolean expression (X>0) per cell for detection
    X_bool = X.copy()
    X_bool.data = (X_bool.data > 0).astype(np.uint8)

    present_rows = []
    stats_rows = []

    for clab in clusters:
        mask = (clustering.values == clab)
        n_cells = int(mask.sum())
        if n_cells == 0:
            continue

        X_sub = X[mask, :]        # cells x genes (CPM)
        Xb_sub = X_bool[mask, :]  # cells x genes (bool)

        # detection fraction per gene
        det_counts = np.asarray(Xb_sub.sum(axis=0)).ravel()
        det_frac = det_counts / float(n_cells)

        # mean CPM per gene
        mean_cpm = np.asarray(X_sub.mean(axis=0)).ravel()

        # call present
        present = (det_frac >= args.min_detect_frac) | (mean_cpm >= args.min_cpm)

        # collect stats (long form)
        for j, gid in enumerate(gene_ids):
            stats_rows.append({
                "cluster": clab,
                "gene_id": gid,
                "detection_fraction": float(det_frac[j]),
                "mean_cpm": float(mean_cpm[j]),
                "present": int(present[j])
            })

    pd.DataFrame(stats_rows).to_csv(outdir / "cluster_gene_stats.csv", index=False)

if __name__ == "__main__":
    main()

