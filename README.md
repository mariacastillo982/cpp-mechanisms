## 🧬 Cell Penetrating Peptide Uptake Mechanisms (CPP-Mechanisms)

Version: 1.0

Created: 2025-10-19
Maintainer: Maria Gomez
Contact: maria.castillo@kaust.edu.sa

Base namespace: https://w3id.org/cpp/schema#
Dataset IRI: https://w3id.org/cpp/dataset/mechanisms

## 📘 Overview

The CPP-Mechanisms Knowledge Base is a structured RDF dataset describing the endocytic molecular mechanisms through which cell-penetrating peptides (CPPs) and related cargos enter cells, such as macropinocytosis, phagocytosis, clathrin-mediated, and caveolae-mediated endocytosis.

# Each mechanism aggregates:

Core and regulatory genes (Ensembl identifiers)

Relevant pathways (GO, Reactome, KEGG, MSigDB)

Known chemical inhibitors (mapped to ChEBI and MeSH)

Literature-based evidence (DOIs)

## 🔍 Data model

Main classes and relationships used:

Entity	Description	Example
```
ex:Mechanism	Cellular uptake mechanism (e.g. macropinocytosis)	mech:Macropinocytosis
ex:GateCoreEntry	Core gene enabling internalization	ensembl:ENSG00000121879 (PIK3CA)
ex:RegulatorEntry	Regulator gene that modulates pathway activity	ensembl:ENSG00000124181 (PLCG1)
ex:ContextNegativeEntry	Drug or chemical that inhibits the mechanism	chebi:CHEBI_47499 (Imipramine)
ex:hasExternalRef	Links to GO, KEGG, Reactome, MSigDB	obo:GO_0044351
ex:evidence	Literature reference (DOI or PubMed link)	<https://doi.org/10.1111/bph.14439>
```
All identifiers are dereferenceable (e.g., Ensembl, ChEBI, Reactome) following the Identifiers.org
 scheme.

## Validation

#### Validation data (Single Cell Expression Atlas)

We validated the endocytosis knowledge base using public GTEx: snRNAseq atlas from Single Cell Expression Atlas (SCEA), the data used for this validation was:

- Normalised counts files: Untransformed expression values, normalised to counts per million. (MatrixMarket: matrix.mtx[.gz], genes.tsv[.gz], barcodes.tsv[.gz]).
- Clustering file: Results of unsupervised louvain clustering at a range of resolution values.


1. Load CPM counts & clustering.
Read the Normalised counts (CPM) and the Clustering file; align barcodes/Cell IDs to get a cell-ID → cluster label for the resolution/label we use. 

2. Aggregate by cluster (cell type).
For each cluster label (e.g., “B cell”, “fat cell”, etc.), subset cells and compute per-gene:

Detection fraction = fraction of cells with CPM > 0
Mean CPM across the cluster

3. Call genes “present.”
A gene is considered expressed/present in a cluster if it meets either threshold (defaults):
- detection_fraction ≥ 0.10 (≥10% cells), OR
- mean CPM ≥ 1.0

4. Write outputs.
```
present_genes_by_cluster.csv: cluster, n_cells, present_ensembl_ids (semicolon-joined)
cluster_gene_stats.csv: long table (cluster, gene_id, detection_fraction, mean_cpm, present)
```

These per-cluster gene sets are then scored against our KB to compute core∪reg and important coverage for each endocytic pathway.

# Validation Overview
- `PathwayKBValidator` loads expression, ground-truth labels, and pathway KB, then hands the KB rules to `PathwayScorer`.
- For each cluster–pathway pair the scorer blends core-module activity with regulator support using the current weights and stores the scores in a tidy DataFrame.
- `PULearningEvaluator` reuses those scores to report PU metrics (precision lower bound, recall on labeled positives, rank percentiles) and a simple CV F1 when reliable negatives are available.
- `WeightTuner` can optionally grid-search the core/regulator weights; it reruns scoring and evaluation for every pair, keeps the best objective (average of the existing metrics), and leaves the validator ready for visualization/exports.

Each pathway was evaluated for precision (lower bound), recall on labeled positives, average rank percentile, and cross-validation F1-score using Positive–Unlabeled (PU) learning.

Pathway	Precision | (LB)	| Recall	| Avg. Rank (%) |	CV F1 ± SD
Phagocytosis	| 0.226	| 1.000	| 74.0	| 0.771 ± 0.390
Macropinocytosis	| 0.322	| 1.000	| 48.3	| 0.541 ± 0.304
Clathrin-mediated endocytosis	| 0.254	| 1.000	| 72.7	| 1.000 ± 0.000
Caveolae-mediated endocytosis	| 0.397	| 1.000	| 67.9	| 0.938 ± 0.081

Overall, all pathways achieved perfect recall on labeled positives, with average rank percentiles between 48–74% and high F1-scores, confirming consistent detection of known active mechanisms across cell clusters.

🔍 Visualization of Results

1. Pathway Activation Heatmap
Activation scores across all cell clusters showing distinct pathway engagement patterns.
📊 activation_heatmap.png

2. Precision, Recall, and Rank Metrics
Aggregate PU-learning metrics per pathway.
📈 performance_metrics.png

3. Score Distributions
Histogram of activation score distributions per pathway (threshold = 0.5).
📉 score_distributions.png


## 🧠 Example SPARQL queries

### 1. List all mechanisms and their GO IDs
```
PREFIX ex: <https://w3id.org/cpp/schema#>
SELECT ?mechanism ?go WHERE {
  ?mechanism a ex:Mechanism ;
             ex:hasExternalRef ?go .
  FILTER(CONTAINS(STR(?go), "obo/GO_"))
}
```

### 2. Retrieve inhibitors (context_negative) and their CHEBI IDs
```
PREFIX ex: <https://w3id.org/cpp/schema#>
SELECT ?mechanism ?drug ?chebi WHERE {
  ?mechanism ex:hasContextNegative [
      rdfs:label ?drug ;
      ex:chemical ?chebi
  ] .
}
```
