# Validation

### Data (Single Cell Expression Atlas)

We validated the endocytosis knowledge base using the public GTEx: snRNAseq atlas from the Single Cell Expression Atlas (SCEA). The data used for this validation were:

First, run the file metadata_to_clustering.py to obtain the results of unsupervised Louvain clustering at a range of resolution values:

```
python metadata_to_clustering.py \
    --metadata "cpp-mechanisms/validation/GTEx/E-ANND-2.cells.txt" \
    --id_col "Cell ID" \
    --cluster_col "inferred cell type - ontology labels" \
    --out "cpp-mechanisms/validation/clustering.tsv"
```

Secondly, run the file create_gene_sets.py to transform the normalised counts files (matrix.mtx[.gz], genes.tsv[.gz], barcodes.tsv[.gz]) to cluster_gene_stats.csv that is a table (cluster, gene_id, detection_fraction, mean_cpm, present).

```
python create_gene_sets.py \
    --counts-mtx "cpp-mechanisms/validation/GTEx/E-ANND-2-normalised-files" \
    --clustering "cpp-mechanisms/validation/GTEx/clustering.tsv" \
    --cluster-col "cluster" \
    --outdir "cpp-mechanisms/validation/GTEx/" \
    --min-detect-frac 0.10 # detection_fraction ≥ 0.10 (≥10% cells)\
    --min-cpm 1.0 # genes with at least 1 count per million reads \
    --accept-prefix ENSG
```
These per-cluster gene sets are then scored against our KB to compute gate_core + regulators coverage for each endocytic pathway.

### Validation Overview

To reproduce the results obtained for validation of the knowledge base, run the file test_kb.py. 
In summary, the code works as follows:

- `PathwayKBValidator` loads expression, ground-truth labels, and pathway KB, then hands the KB rules to `PathwayScorer`.
- For each cluster–pathway pair, the scorer blends gate_core coverage with regulators coverage using the optimized weights and stores the scores in a table.
- `PULearningEvaluator` reuses those scores to report PU metrics (precision lower bound, recall on labeled positives, rank percentiles) and a simple CV F1 with reliable negatives that are being generated.
- `WeightTuner` can optionally grid-search the core/regulator weights; it reruns scoring and evaluation for every pair, keeps the best objective (average of the existing metrics), and leaves the validator ready for visualization/exports.

Each pathway was evaluated for precision (lower bound), recall on labeled positives, average rank percentile, and cross-validation F1-score using Positive–Unlabeled (PU) learning.

|Pathway |	Precision (LB)	| Recall	| Avg. Rank (%) |	CV F1 ± SD|
|---|---|---|---|---|
|Phagocytosis	| 0.226	| 1.000	| 74.0	| 0.771 ± 0.390|
|Macropinocytosis	| 0.322	| 1.000	| 48.3	| 0.541 ± 0.304|
|Clathrin-mediated endocytosis	| 0.254	| 1.000	| 72.7	| 1.000 ± 0.000|
|Caveolae-mediated endocytosis	| 0.397	| 1.000	| 67.9	| 0.938 ± 0.081|

Overall, all pathways achieved perfect recall on labeled positives, with average rank percentiles between 48–74%, confirming consistent detection of known active mechanisms across cell clusters, with low precision since there are no negative labels. It was prioritized to maximize recall.

### 🔍 Visualization of Results

1. Pathway Activation Heatmap
Activation scores across all cell clusters showing distinct pathway engagement patterns.
📊 activation_heatmap.png

2. Precision, Recall, and Rank Metrics
Aggregate PU-learning metrics per pathway.
📈 performance_metrics.png

3. Score Distributions
Histogram of activation score distributions per pathway (threshold = 0.5).
📉 score_distributions.png


