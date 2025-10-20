import json
import argparse
from urllib.parse import quote
from rdflib import Graph, Namespace, URIRef, BNode, Literal
from rdflib.namespace import RDF, RDFS, XSD

EX = Namespace("http://example.org/schema#")
MECH = Namespace("http://example.org/mechanism/")

# ---------- Helpers to normalize IDs -> resolvable IRIs ----------

def iri_for_go(go_token: str) -> URIRef:
    # Accepts "GO:0044351" -> OBO PURL
    go_id = go_token.strip().replace(" ", "")
    if ":" in go_id:
        go_id = go_id.split(":")[1]
    return URIRef(f"http://purl.obolibrary.org/obo/GO_{go_id}")

def iri_for_reactome(rid: str) -> URIRef:
    # e.g., "R-HSA-8854222"
    return URIRef(f"https://identifiers.org/reactome/{rid.strip()}")

def iri_for_kegg(token: str) -> URIRef:
    # Accepts:
    #   "hsa04810"          -> pathway
    #   "KO: K17917"        -> orthology (strip leading "KO:")
    t = token.strip()
    if t.upper().startswith("KO"):
        # handle "KO: K17917" or "KO:K17917"
        ko = t.split(":")[-1].strip()
        return URIRef(f"https://identifiers.org/kegg.orthology/{ko}")
    else:
        # route generic KEGG codes to identifiers.org as well
        return URIRef(f"https://identifiers.org/kegg:{t}")

def iri_for_msigdb(token: str) -> URIRef:
    t = token.strip()
    return URIRef(f"https://identifiers.org/msigdb/{quote(t)}")

def iri_for_ensembl(ens: str) -> URIRef:
    return URIRef(f"https://identifiers.org/ensembl/{ens.strip()}")

def iri_for_chebi(chebi_token: str) -> URIRef:
    # Accepts "CHEBI:47499" or "47499"
    t = chebi_token.strip().upper()
    if t.startswith("CHEBI:"):
        t = t.split(":")[1]
    return URIRef(f"http://purl.obolibrary.org/obo/CHEBI_{t}")

def iri_for_mesh(mesh_token: str) -> URIRef:
    # Accepts "mesh:C000719991" or "C000719991"
    t = mesh_token.strip()
    if ":" in t:
        t = t.split(":")[1]
    return URIRef(f"https://id.nlm.nih.gov/mesh/{t}")

def iri_for_doi_or_url(s: str) -> URIRef | Literal:
    if not s:
        return Literal("")
    t = s.strip()
    # very simple normalization for DOIs
    if t.lower().startswith("doi:"):
        t = t[4:].strip()
    if t.startswith("10."):
        return URIRef(f"https://doi.org/{t}")
    if t.startswith("http://") or t.startswith("https://"):
        return URIRef(t)
    # fall back to literal if not obviously an IRI
    return Literal(t)

def mech_uri(key: str) -> URIRef:
    # mint a stable mechanism URI from the top-level key
    return URIRef(MECH[str(quote(key.strip()))])

# ---------- Core mapping functions ----------

def add_external_refs(g: Graph, mech: URIRef, ext: dict):
    # ext = { "GO": [...], "Reactome": [...], "KEGG": [...], "MSigDB": [...] }
    if not isinstance(ext, dict):
        return
    for k, vals in ext.items():
        if not isinstance(vals, (list, tuple)):
            continue
        for v in vals:
            if not v: 
                continue
            try:
                if k.upper() == "GO":
                    g.add((mech, EX.hasExternalRef, iri_for_go(v)))
                elif k.lower() == "reactome":
                    g.add((mech, EX.hasExternalRef, iri_for_reactome(v)))
                elif k.upper() == "KEGG":
                    g.add((mech, EX.hasExternalRef, iri_for_kegg(v)))
                elif k.upper() == "MSIGDB":
                    g.add((mech, EX.hasExternalRef, iri_for_msigdb(v)))
                else:
                    # Unknown external type -> store as literal with type tag
                    bn = BNode()
                    g.add((mech, EX.hasExternalRef, bn))
                    g.add((bn, EX.source, Literal(k)))
                    g.add((bn, EX.identifier, Literal(str(v))))
            except Exception:
                # If parsing fails, keep as literal node
                bn = BNode()
                g.add((mech, EX.hasExternalRef, bn))
                g.add((bn, EX.source, Literal(k)))
                g.add((bn, EX.identifier, Literal(str(v))))

def add_gene_entry(g: Graph, mech: URIRef, slot_predicate: URIRef, entry: dict, entry_type: URIRef):
    bn = BNode()
    g.add((mech, slot_predicate, bn))
    g.add((bn, RDF.type, entry_type))

    # "gene" + optional "ensembl_id"
    gene = entry.get("gene")
    if gene:
        g.add((bn, EX.geneSymbol, Literal(gene)))

    ensembl = entry.get("ensembl_id") or entry.get("ensemble_id") or entry.get("ensemblId")
    if ensembl:
        g.add((bn, EX.ensembl, iri_for_ensembl(ensembl)))

    # common text fields
    if "why" in entry:
        g.add((bn, EX.why, Literal(entry["why"])))
    if "type" in entry:
        g.add((bn, EX.type, Literal(entry["type"])))
    if "dose" in entry:
        g.add((bn, EX.dose, Literal(entry["dose"])))

    # evidence: allow DOI/URL or literal
    if "evidence" in entry:
        ev = iri_for_doi_or_url(entry["evidence"])
        if isinstance(ev, URIRef):
            g.add((bn, EX.evidence, ev))
        else:
            g.add((bn, EX.evidence, ev))

def add_small_molecule_entry(g: Graph, mech: URIRef, slot_predicate: URIRef, entry: dict):
    bn = BNode()
    g.add((mech, slot_predicate, bn))
    g.add((bn, RDF.type, EX.ContextNegativeEntry))  # reuse this type for inhibitors/conditions

    # label/name
    name = entry.get("name")
    if name:
        g.add((bn, RDFS.label, Literal(name)))

    # identifiers (CHEBI, MESH), may be nested or direct string
    chebi = entry.get("CHEBI")
    if isinstance(chebi, dict):
        # e.g., {"Amiloride": "CHEBI:2639", "EIPA": "CHEBI:136538"}
        for k, v in chebi.items():
            if v:
                g.add((bn, EX.chemical, iri_for_chebi(v)))
    elif isinstance(chebi, str) and chebi.strip():
        g.add((bn, EX.chemical, iri_for_chebi(chebi)))

    mesh = entry.get("MESH") or entry.get("Mesh") or entry.get("mesh")
    if isinstance(mesh, str) and mesh.strip():
        g.add((bn, EX.mesh, iri_for_mesh(mesh)))

    # common text fields
    for field, pred in [
        ("type", EX.type),
        ("dose", EX.dose),
        ("why", EX.why),
    ]:
        if field in entry and entry[field]:
            g.add((bn, pred, Literal(entry[field])))

    # evidence (IRI or literal)
    if "evidence" in entry and entry["evidence"]:
        ev = iri_for_doi_or_url(entry["evidence"])
        if isinstance(ev, URIRef):
            g.add((bn, EX.evidence, ev))
        else:
            g.add((bn, EX.evidence, ev))

def handle_array(g: Graph, mech: URIRef, key: str, arr):
    """
    Routes arrays by role:
        - gate_core      -> EX.hasGateCore + EX.GateCoreEntry
        - regulators     -> EX.hasRegulator + EX.RegulatorEntry
        - context_negative -> EX.hasContextNegative + EX.ContextNegativeEntry
    Supports entries that are gene-based (have 'gene') or small molecules (have 'name').
    """
    if not isinstance(arr, (list, tuple)):
        return
    role_map = {
        "gate_core": (EX.hasGateCore, EX.GateCoreEntry),
        "regulators": (EX.hasRegulator, EX.RegulatorEntry),
        "context_negative": (EX.hasContextNegative, EX.ContextNegativeEntry),
    }
    pred, etype = role_map.get(key, (EX.hasEntry, EX.Entry))
    for item in arr:
        if not isinstance(item, dict):
            continue
        if "gene" in item:
            add_gene_entry(g, mech, pred, item, etype)
        else:
            # treat as small-molecule or generic entry (has 'name' or just descriptors)
            add_small_molecule_entry(g, mech, pred, item)

def json_to_graph(data: dict, base_mech_ns: Namespace = MECH) -> Graph:
    g = Graph()
    g.bind("ex", EX)
    g.bind("mech", base_mech_ns)
    g.bind("rdfs", RDFS)
    g.bind("xsd", XSD)

    # Each top-level key is a mechanism
    for mech_name, mech_obj in data.items():
        if not isinstance(mech_obj, dict):
            continue
        mech = mech_uri(mech_name)
        g.add((mech, RDF.type, EX.Mechanism))
        g.add((mech, RDFS.label, Literal(mech_name)))

        # simple fields
        if "id" in mech_obj:
            g.add((mech, EX.id, Literal(mech_obj["id"])))
        if "label" in mech_obj and mech_obj["label"] != mech_name:
            g.add((mech, RDFS.label, Literal(mech_obj["label"])))
        if "class" in mech_obj:
            g.add((mech, EX.class_, Literal(mech_obj["class"])))

        # external refs
        if "external_refs" in mech_obj:
            add_external_refs(g, mech, mech_obj["external_refs"])

        # array slots
        for role_key in ("gate_core", "regulators", "context_negative"):
            if role_key in mech_obj:
                handle_array(g, mech, role_key, mech_obj[role_key])

    return g

# ---------- CLI ----------

def main():
    ap = argparse.ArgumentParser(description="Convert mechanisms JSON to RDF.")
    ap.add_argument("--in", dest="infile", required=True, help="Input JSON file")
    ap.add_argument("--out", dest="outfile", required=True, help="Output RDF file")
    ap.add_argument("--format", dest="fmt", default="turtle",
                    help="RDF format (turtle|xml|ntriples|nt|json-ld|trig|nquads)")
    ap.add_argument("--graph", dest="graph_uri", default=None,
                    help="Optional named graph URI (wraps output in one graph for Trig/N-Quads)")
    args = ap.parse_args()

    with open(args.infile, "r", encoding="utf-8") as f:
        data = json.load(f)

    g = json_to_graph(data)

    # Optional: write into a named graph (requires TriG/N-Quads serialization)
    if args.graph_uri and args.fmt.lower() in ("trig", "nquads"):
        from rdflib import ConjunctiveGraph
        cg = ConjunctiveGraph()
        cg.bind("ex", EX); cg.bind("mech", MECH); cg.bind("rdfs", RDFS); cg.bind("xsd", XSD)
        ng = cg.get_context(URIRef(args.graph_uri))
        for t in g:
            ng.add(t)
        cg.serialize(destination=args.outfile, format=args.fmt)
    else:
        g.serialize(destination=args.outfile, format=args.fmt)

if __name__ == "__main__":
    main()
