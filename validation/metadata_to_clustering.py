#!/usr/bin/env python3
# metadata_to_clustering.py
import argparse, pandas as pd
from pathlib import Path

def main():
    ap = argparse.ArgumentParser(description="Make cell→cluster mapping from SCEA metadata TSV.")
    ap.add_argument("--metadata", required=True, help="Path to SCEA metadata TSV/CSV")
    ap.add_argument("--id_col", required=True,
                    help="Column to use as cell identifier (e.g., 'Cell ID' or 'barcode')")
    ap.add_argument("--cluster_col", required=True,
                    help="Column to use as the cluster label (e.g., 'inferred cell type - ontology labels')")
    ap.add_argument("--out", required=True, help="Output TSV (columns: cell_id, cluster)")
    args = ap.parse_args()

    df = pd.read_csv(args.metadata, sep=None, engine="python", dtype=str)
    # Be robust to column whitespace
    df.columns = [c.strip() for c in df.columns]

    if args.id_col not in df.columns or args.cluster_col not in df.columns:
        raise SystemExit(f"Columns not found. Available: {list(df.columns)}")

    out = (
        df[[args.id_col, args.cluster_col]]
        .rename(columns={args.id_col: "cell_id", args.cluster_col: "cluster"})
        .dropna()
    )
    out["cell_id"] = out["cell_id"].astype(str).str.strip()
    out["cluster"] = out["cluster"].astype(str).str.strip()
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(args.out, index=False)
    print(f"Wrote {len(out)} rows to {args.out}")

if __name__ == "__main__":
    main()


"""
python metadata_to_clustering.py \
    --metadata "cpp-mechanisms/validation/GTEx/E-ANND-2.cells.txt" \
    --id_col "Cell ID" \
    --cluster_col "inferred cell type - ontology labels" \
    --out "cpp-mechanisms/validation/GTEx/clustering.tsv"
"""
