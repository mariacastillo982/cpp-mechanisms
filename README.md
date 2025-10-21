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
