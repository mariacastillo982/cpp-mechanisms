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

Genes (Ensembl identifiers)

Cellular uptake mechanism (GO)

Chemical inhibitors (ChEBI)

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
