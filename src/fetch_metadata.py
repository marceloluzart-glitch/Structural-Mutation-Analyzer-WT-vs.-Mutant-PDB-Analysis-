"""
GenBank Metadata Harvester
Description: Queries NCBI Entrez using BioPython to extract epidemiological metadata 
from FASTA segment files.
"""

from Bio import SeqIO, Entrez
import pandas as pd
import time
import os

# =====================================
# CONFIGURATION
# =====================================
Entrez.email = "your_email_here@gmail.com"  # Replace with your email

# Local FASTA input files dictionary
arquivos = {
    "L": "data/raw/sequence_SEGMENTO_L.fasta",
    "M": "data/raw/sequence_SEGMENTO_M.fasta",
    "S": "data/raw/sequence_SEGMENTO_S.fasta"
}

# =====================================
# EXTRACT ACCESSIONS
# =====================================
dados = []

for segmento, arquivo in arquivos.items():
    print(f"\n[INFO] Processing segment {segmento}...")
    
    if not os.path.exists(arquivo):
        print(f"[WARNING] File not found: {arquivo}. Skipping...")
        continue

    registros = SeqIO.parse(arquivo, "fasta")

    for registro in registros:
        accession = registro.id.split()[0]
        print(f"Fetching accession: {accession}")

        try:
            handle = Entrez.efetch(
                db="nucleotide",
                id=accession,
                rettype="gb",
                retmode="text"
            )

            gb_record = SeqIO.read(handle, "genbank")

            metadata = {
                "accession": accession,
                "segmento": segmento,
                "pais": "",
                "estado": "",
                "data_coleta": "",
                "hospedeiro": "",
                "strain": ""
            }

            for feature in gb_record.features:
                if feature.type == "source":
                    qualifiers = feature.qualifiers

                    # Country and location
                    country = qualifiers.get("country", [""])[0]
                    metadata["pais"] = country

                    # State/locality parsing
                    if ":" in country:
                        metadata["estado"] = country.split(":")[1].strip()

                    # Collection date
                    metadata["data_coleta"] = qualifiers.get(
                        "collection_date",
                        [""]
                    )[0]

                    # Host
                    metadata["hospedeiro"] = qualifiers.get(
                        "host",
                        [""]
                    )[0]

                    # Strain
                    metadata["strain"] = qualifiers.get(
                        "strain",
                        [""]
                    )[0]

            dados.append(metadata)
            time.sleep(0.4)  # Respect NCBI API rate limits

        except Exception as e:
            print(f"[ERROR] Failed on {accession}: {e}")

# =====================================
# DATAFRAME & EXPORTATION
# =====================================
if dados:
    df = pd.DataFrame(dados)
    
    # Ensure processed data directory exists
    os.makedirs("data/processed", exist_ok=True)
    output_csv = "data/processed/viral_metadata_summary.csv"
    df.to_csv(output_csv, index=False)

    print("\nSummary:")
    print(df.head())
    print(f"\n[SUCCESS] Metadata successfully saved to {output_csv}")
else:
    print("\n[WARNING] No data retrieved.")
