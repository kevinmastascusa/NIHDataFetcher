"""
This script demonstrates how to access and utilize data from various NIH databases.
"""

import GEOparse
from Bio import Entrez
import pandas as pd
import gzip


def fetch_geo_data(geo_id: str, destdir: str = './'):
    """
    Fetch a GEO dataset with a given ID.

    :param geo_id: The GEO dataset ID.
    :param destdir: The directory where the fetched dataset should be stored.
    """
    gse = GEOparse.get_GEO(geo=geo_id, destdir=destdir)
    print(f"Successfully fetched GEO dataset {geo_id} and stored it in {destdir}.")


def search_pubmed(term: str, email: str = "kzm24@drexel.edu"):
    """
    Search for a term in the PubMed database.

    :param term: The search term.
    :param email: The email to use with the Entrez API.
    Replace with your actual email when using the program.
    """
    Entrez.email = email
    handle = Entrez.esearch(db="pubmed", term=term)
    record = Entrez.read(handle)
    handle.close()
    idlist = record["IdList"]
    print(f"Found {len(idlist)} PubMed articles matching the term: {term}.")


def load_clinvar_data(filepath: str):
    """
    Load ClinVar data from a .vcf file.

    :param filepath: The path to the .vcf file.
    """
    data = pd.read_csv(filepath, comment='#', delimiter='\t')
    print(f"Successfully loaded ClinVar data from {filepath}.")


def load_vcf_file(file_path):
    """
    Load and decompress a VCF file.

    :param file_path: The path to the compressed VCF file (.vcf.gz).
    :return: The decompressed VCF file contents.
    """
    with gzip.open(file_path, "rt") as vcf_file:
        vcf_contents = vcf_file.read()
    return vcf_contents


if __name__ == "__main__":
    # Example test calls for each function
    fetch_geo_data("GSE2553")
    search_pubmed("asthma")
    vcf_file_path = "clinvar.vcf.gz"  # Replace with the path to your specific VCF file
    vcf_data = load_vcf_file(vcf_file_path)
    print(vcf_data)
    load_clinvar_data(vcf_file_path)
