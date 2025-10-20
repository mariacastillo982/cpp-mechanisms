🧬 Cell Penetrating Peptide Uptake Mechanisms (CPP-Mechanisms)

Version: 1.0
License: CC BY 4.0

Created: 2025-10-19
Maintainer: Maria Gomez
Contact: maria.castillo@kaust.edu.sa

Base namespace: https://w3id.org/cpp/schema#
Dataset IRI: https://w3id.org/cpp/dataset/mechanisms

📘 Overview

The CPP-Mechanisms Knowledge Base is a structured RDF dataset describing the endocytic molecular mechanisms through which cell-penetrating peptides (CPPs) and related cargos enter cells such as macropinocytosis, phagocytosis, clathrin mediated and caveolae mediated endocytosis.

Each mechanism aggregates:

Core and regulatory genes (Ensembl identifiers)

Relevant pathways (GO, Reactome, KEGG, MSigDB)

Known chemical inhibitors (mapped to ChEBI and MeSH)

Literature-based evidence (DOIs)

The dataset is published as FAIR, machine-readable RDF and can be queried with SPARQL or reused in ontology-driven analyses.

🗂️ Directory structure
cpp-mechanisms/
├── data/
│   ├── mechanisms.json          # Original JSON input
│   ├── mechanisms.ttl           # RDF Turtle serialization
│   ├── mechanisms.jsonld        # RDF JSON-LD serialization
├── metadata/
│   └── dataset.ttl              # DCAT/PROV metadata file
├── shapes/
│   └── mechanism.shacl.ttl      # SHACL validation shapes
├── src/
│   └── json_to_rdf.py           # JSON → RDF conversion script
├── docs/
│   └── index.html               # Optional landing page for GitHub Pages
└── README.md                    # This file

🔍 Data model

Main classes and relationships used:

Entity	Description	Example
ex:Mechanism	Cellular uptake mechanism (e.g. macropinocytosis)	mech:Macropinocytosis
ex:GateCoreEntry	Core gene enabling internalization	ensembl:ENSG00000121879 (PIK3CA)
ex:RegulatorEntry	Regulator gene that modulates pathway activity	ensembl:ENSG00000124181 (PLCG1)
ex:ContextNegativeEntry	Drug or chemical that inhibits the mechanism	chebi:CHEBI_47499 (Imipramine)
ex:hasExternalRef	Links to GO, KEGG, Reactome, MSigDB	obo:GO_0044351
ex:evidence	Literature reference (DOI or PubMed link)	<https://doi.org/10.1111/bph.14439>

All identifiers are dereferenceable (e.g., Ensembl, ChEBI, Reactome) following the Identifiers.org
 scheme.

⚙️ How it was generated

Source JSON manually curated from literature on endocytic uptake pathways.

Converted to RDF using src/json_to_rdf.py

python src/json_to_rdf.py --in data/mechanisms.json --out data/mechanisms.ttl --format turtle


Validation

riot --validate data/mechanisms.ttl
pyshACL -s shapes/mechanism.shacl.ttl -d data/mechanisms.ttl -m -f human


Metadata annotated using DCAT, Dublin Core, and PROV-O.

Published to GitHub Pages and linked via a persistent w3id.org
 URI.

🔗 FAIR implementation summary
FAIR Principle	Implementation
F1. Findable	Persistent URIs via https://w3id.org/cpp/; metadata in metadata/dataset.ttl
A1. Accessible	Publicly available on GitHub/Zenodo; HTTP(S) access; CC-BY 4.0 license
I1. Interoperable	RDF (Turtle, JSON-LD); use of standard vocabularies (GO, KEGG, Reactome, ChEBI, Ensembl, DCAT, PROV)
R1. Reusable	Open license, provenance metadata, clear SHACL validation shapes
🧠 Example SPARQL queries

1. List all mechanisms and their GO IDs

PREFIX ex: <https://w3id.org/cpp/schema#>
SELECT ?mechanism ?go WHERE {
  ?mechanism a ex:Mechanism ;
             ex:hasExternalRef ?go .
  FILTER(CONTAINS(STR(?go), "obo/GO_"))
}


2. Retrieve inhibitors (context_negative) and their CHEBI IDs

PREFIX ex: <https://w3id.org/cpp/schema#>
SELECT ?mechanism ?drug ?chebi WHERE {
  ?mechanism ex:hasContextNegative [
      rdfs:label ?drug ;
      ex:chemical ?chebi
  ] .
}

📚 References (examples)

Commisso et al., Cell, 2013 — Macropinocytosis and nutrient uptake.

Kerr & Teasdale, Traffic, 2009 — Regulation of macropinocytosis.

Palamidessi et al., Curr Biol, 2021 — Cdc42 control in endocytosis.

DOIs: 10.1111/bph.14439 · 10.1016/j.ceb.2025.102563 · 10.1073/pnas.1201573109

🧩 Validation & Quality Control

Syntax check: Apache Jena RIOT

Shape validation: pySHACL with shapes/mechanism.shacl.ttl

Ontology reuse check: Verified GO, KEGG, Reactome, ChEBI, and Ensembl IDs resolve through Identifiers.org.

🧾 License

This dataset is released under the Creative Commons Attribution 4.0 International License (CC BY 4.0).
You are free to share, adapt, and use it commercially with proper attribution.

© 2025 Xavier López — CC BY 4.0

🧮 Citation

If you use this dataset, please cite:

Gomez, M. (2025). CPP Endocytic Mechanisms Knowledge Base (v1.0).
RDF dataset published at https://w3id.org/cpp/dataset/mechanisms.
DOI: [Add Zenodo DOI here]

🧱 Contributing

Pull requests are welcome!
If you add new mechanisms or corrections:

Update data/mechanisms.json

Re-run src/json_to_rdf.py

Validate using SHACL

Update metadata version and changelog

