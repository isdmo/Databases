import streamlit as st
import requests
import pandas as pd
import matplotlib.pyplot as plt

# --- Hardcoded mappings ---
gene_to_ensembl = {
    "TP53": "ENSG00000141510",
    "BRCA1": "ENSG00000012048",
    "EGFR": "ENSG00000146648"
}
ensembl_to_gene = {v: k for k, v in gene_to_ensembl.items()}

# --- Database links ---
database_categories = {
    "Gene / Genome": {
        "Ensembl": "https://www.ensembl.org/Homo_sapiens/Gene/Summary?g={ensembl_id}",
        "GeneCards": "https://www.genecards.org/cgi-bin/carddisp.pl?gene={gene}",
        "NCBI Gene": "https://www.ncbi.nlm.nih.gov/gene/?term={gene}"
    },
    "Protein & Localization": {
        "UniProt": "https://www.uniprot.org/uniprotkb?query={gene}",
        "Protein Atlas": "https://www.proteinatlas.org/{ensembl_id}",
        "GTEx": "https://gtexportal.org/home/gene/{gene}"
    },
    "Population, Variation & Disease": {
        "gnomAD": "https://gnomad.broadinstitute.org/gene/{gene}",
        "ClinVar": "https://www.ncbi.nlm.nih.gov/clinvar/?term={gene}",
        "dbSNP": "https://www.ncbi.nlm.nih.gov/snp/?term={gene}",
        "OMIM": "https://omim.org/search?search={gene}",
        "cBioPortal": "https://www.cbioportal.org/results/cancerTypesSummary?geneId={gene}"
    },
    "Expression & Functional": {
        "OpenCell": "https://opencell.czbiohub.org/gene/{ensembl_id}",
        "DepMap": "https://depmap.org/portal/gene/{gene}"
    },
    "Pathways & Interactions": {
        "STRING": "https://string-db.org/cgi/network?identifiers={gene}&species=9606",
        "Reactome": "https://reactome.org/content/query?q={gene}",
        "UCSC Genome": "https://genome.ucsc.edu/cgi-bin/hgGene?db=hg38&hgg_gene={gene}",
    },
    "Structure": {
        "AlphaFold DB": "https://alphafold.ebi.ac.uk/search/text/{gene}",
        "PDB": "https://www.rcsb.org/search?request={%22query%22:{%22type%22:%22terminal%22,%22service%22:%22text%22,%22parameters%22:{%22value%22:%22{gene}%22}},%22return_type%22:%22entry%22}"
    }
}

# --- Fetch helpers ---
def fetch_from_ensembl_symbol(gene_symbol, species="human"):
    url = f"https://rest.ensembl.org/xrefs/symbol/{species}/{gene_symbol}?"
    r = requests.get(url, headers={"Content-Type": "application/json"})
    if r.ok:
        for entry in r.json():
            if entry.get("type") == "gene":
                return entry["id"]
    return None

def fetch_symbol_from_ensembl(ensembl_id):
    url = f"https://rest.ensembl.org/lookup/id/{ensembl_id}?content-type=application/json"
    r = requests.get(url)
    if r.ok:
        return r.json().get("display_name")
    return None

def fetch_ensembl_metadata(ensembl_id):
    url = f"https://rest.ensembl.org/lookup/id/{ensembl_id}?content-type=application/json"
    r = requests.get(url)
    if r.ok:
        return r.json()
    return None

def fetch_uniprot_metadata(gene):
    url = f"https://rest.uniprot.org/uniprotkb/search?query=gene:{gene}+AND+organism_id:9606&format=json"
    r = requests.get(url)
    if r.ok and r.json().get("results"):
        data = r.json()["results"][0]
        return {
            "UniProt ID": data["primaryAccession"],
            "Protein Name": data["proteinDescription"]["recommendedName"]["fullName"]["value"],
            "Organism": data["organism"]["scientificName"]
        }
    return None

# --- Streamlit app ---
st.title("🔎 Gene Database Lookup")

user_input = st.text_input("Enter a Gene Symbol or Ensembl ID (e.g., TP53, ENSG00000141510)")

if user_input:
    gene = None
    ensembl_id = None

    # Case 1: Hardcoded dict
    if user_input in gene_to_ensembl:
        gene = user_input.upper()
        ensembl_id = gene_to_ensembl[user_input]
    elif user_input in ensembl_to_gene:
        ensembl_id = user_input.upper()
        gene = ensembl_to_gene[user_input]
    else:
        # Case 2: Try Ensembl API
        if user_input.startswith("ENSG"):
            ensembl_id = user_input.upper()
            gene = fetch_symbol_from_ensembl(ensembl_id)
        else:
            gene = user_input.upper()
            ensembl_id = fetch_from_ensembl_symbol(gene)

    if gene and ensembl_id:
        st.subheader(f"Results for {gene} ({ensembl_id})")

        # Database Links as buttons grouped by category
        st.markdown("### 🔗 Database Links")
        for category, dbs in database_categories.items():
            st.markdown(f"**{category}**")
            cols = st.columns(4)
            idx = 0
            for db_name, url_pattern in dbs.items():
                url = url_pattern.format(gene=gene, ensembl_id=ensembl_id)
                # Use an anchor wrapping a simple button so the link opens in a new tab
                btn_html = (
                    f'<a href="{url}" target="_blank" style="text-decoration:none;">'
                    f'<button style="padding:1px 10px;border-radius:6px;cursor:pointer;">{db_name}</button>'
                    '</a>'
                )
                cols[idx % 4].markdown(btn_html, unsafe_allow_html=True)
                idx += 1

        # Ensembl metadata
        with st.expander("📊 Ensembl Metadata"):
            ensembl_meta = fetch_ensembl_metadata(ensembl_id)
            if ensembl_meta:
                meta_table = pd.DataFrame({
                    "Property": ["Gene Symbol", "Description", "Biotype", "Chromosome", "Start", "End"],
                    "Value": [
                        ensembl_meta.get("display_name"),
                        ensembl_meta.get("description"),
                        ensembl_meta.get("biotype"),
                        ensembl_meta.get("seq_region_name"),
                        ensembl_meta.get("start"),
                        ensembl_meta.get("end")
                    ]
                })
                # Ensure all values are strings so Streamlit/pyarrow can serialize the table
                meta_table["Value"] = meta_table["Value"].apply(lambda v: "" if v is None else str(v))
                st.table(meta_table)

        # UniProt metadata
        with st.expander("🧬 UniProt Metadata"):
            uniprot_meta = fetch_uniprot_metadata(gene)
            if uniprot_meta:
                rows = [(k, "" if v is None else str(v)) for k, v in uniprot_meta.items()]
                st.table(pd.DataFrame(rows, columns=["Property", "Value"]))

    else:
        st.error("❌ Could not resolve gene or Ensembl ID.")
