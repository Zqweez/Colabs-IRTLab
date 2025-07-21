"""
BLAST isolates against a database of annotated ASVs.
Database: ASVs PHE Annotated Database.
Runs blast through subprocess.
"""
import subprocess
from Bio import SeqIO
from pathlib import Path
import pandas as pd
from tqdm import tqdm

blast_db = Path("databases/ASVs/all_asvs_phe_annotated_db")
input_fasta = Path("data/sanger/formatted_consensus.fasta") #Path("results/ngs/all_asvs_phe.fasta")
output_dir = Path("results/ngs/isolate_ASVs_blast_results")

def run_blast():
    """
    Run BLAST on the input FASTA file against the ASVs PHE Annotated Database.
    """
    
    output_dir.mkdir(exist_ok=True)
    # Run BLAST
    all_hits = []
    sequences = list(SeqIO.parse(input_fasta, "fasta"))

    for record in tqdm(sequences, desc="Running BLAST"):
        query_id = record.id
        temp_query = output_dir / f"{query_id}.fasta"

        # Write the individual query FASTA
        SeqIO.write(record, temp_query, "fasta")
        # Prepare output file names
        temp_out_blast = output_dir / f"{query_id}_blast.tsv"

        # Prepare the BLAST command
        blast_command = [
            "blastn",
            "-query", str(temp_query),
            "-db", blast_db,
            "-out", str(temp_out_blast),
            "-outfmt", "6 qseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle",
            "-num_threads", "8"  # Adjust number of threads as needed
        ]

        try:
            subprocess.run(blast_command, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
            # Read result
            if temp_out_blast.exists() and temp_out_blast.stat().st_size > 0:
                df = pd.read_csv(temp_out_blast, sep="\t", header=None)
                df[0] = query_id  # Add query name in case qseqid is truncated
                #print(df.head())  # Debug: print first few rows
                all_hits.append(df)

        except subprocess.CalledProcessError:
            print(f"BLAST failed for {query_id}")

    # --- Combine results ---
    if all_hits:
        results_df = pd.concat(all_hits)
        results_df.columns = [
            "query", "pident", "length", "mismatch", "gapopen",
            "qstart", "qend", "sstart", "send", "evalue", "bitscore", "taxa"
        ]
        # Sort the results by query and bitscore
        results_df = (results_df
              .sort_values(by=["query", "bitscore"], ascending=[True, False])
              .reset_index(drop=True))
        results_df.to_csv(output_dir / "blast_results.tsv", sep="\t", index=False)
        print(f"BLAST results saved to {output_dir / 'blast_results.tsv'}")

        # Write top hit
        top_hits = results_df.groupby("query", as_index=False).head(1)
        top_hits.to_csv(output_dir / "blast_top_hit.tsv", sep="\t", index=False)
        print(f"Top hit saved to {output_dir / 'blast_top_hit.tsv'}")

        # Write top hits
        top_hits = results_df.groupby("query", as_index=False).head(5)
        top_hits.to_csv(output_dir / "blast_top_5_hit.tsv", sep="\t", index=False)
        print(f"Top 5 hits saved to {output_dir / 'blast_top_5_hit.tsv'}")
    else:
        print("No BLAST hits found.")

if __name__ == "__main__":
    run_blast()