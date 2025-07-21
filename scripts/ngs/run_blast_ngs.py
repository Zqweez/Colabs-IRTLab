"""
BLAST 16S rRNA sequences against a reference database.
Database: SILVA 16S rRNA database and ezBioCloud 16S rRNA database.
Runs blast though subprocess.
Requires:  blastn  (conda install -c bioconda blast)
"""
import subprocess
from Bio import SeqIO
from pathlib import Path
from tqdm import tqdm
import pandas as pd
import re

# --- Config ---
process_sample = "none" # "phe", "nap", "none"
input_fasta = Path(f"results/ngs/{process_sample}/all_asvs_{process_sample}.fasta")
blast_db = Path("databases/ezbiocloud/ezbio_db") # ezbiocloud/ezbio_db or SILVA/silva_nr99_db  # base name only, no .nsq etc.
vsearch_db = Path("databases/ezbiocloud/ezbiocloud.fa") # VSEARCH database, if using VSEARCH

output_dir = Path(f"results/ngs/{process_sample}/taxa_files")
output_dir.mkdir(exist_ok=True)
# Final combined output
combined_output = output_dir / "blast_results.tsv"
top_hits_output = output_dir / "blast_top_hit.tsv"
top_hits_compact_output = output_dir / "blast_top_hit_compact.tsv"
annotated_fasta_output = f"results/ngs/{process_sample}/all_asvs_{process_sample}_annotated.fasta"
sintax_combined_output = output_dir / "sintax_results.tsv"

def create_annotated_fasta():
    """
    Create an annotated FASTA file with top hit species added to headers.
    """
    # Read the top hits TSV file with pandas
    if top_hits_compact_output.exists():
        top_hits_df = pd.read_csv(top_hits_compact_output, sep="\t")
        # Create a dictionary for quick lookup of taxa by query ID
        taxa_dict = dict(zip(top_hits_df['query'], top_hits_df['taxa']))
    else:
        print(f"Warning: {top_hits_compact_output} not found. Creating annotated FASTA with 'No_hit_found' for all sequences.")
        taxa_dict = {}
    
    # Read original sequences and write annotated versions
    with open(annotated_fasta_output, 'w') as out_file:
        for record in SeqIO.parse(input_fasta, "fasta"):
            # Get the top hit taxa for this sequence
            top_hit_taxa = taxa_dict.get(record.id, "No_hit_found")
            
            # Extract species name from taxa string
            # Taxa format is typically like "Genus species strain_info" or similar
            species_name = extract_species_name(top_hit_taxa)
            
            # Create new header with original ID and species annotation
            new_header = f"{record.id}-{species_name}"
            
            # Write the annotated sequence
            out_file.write(f">{new_header}\n")
            out_file.write(f"{record.seq}\n")

def extract_species_name(taxa_string):
    """
    Extract species name from the taxa string.
    taxa_string: Full taxonomic description in format like:
        
    Returns:
        Species name in format "Genus_species"
    """
    # Extract genus: remove everything up to "g:" and everything after the next comma
    genus_match = re.search(r'g:([^,]+)', taxa_string)
    if genus_match:
        genus = genus_match.group(1)
    else:
        return "Unknown_species"
    
    # Extract species: remove everything up to "s:" and then take the species epithet (word after genus)
    species_match = re.search(r's:([^,]+)', taxa_string)
    if species_match:
        species_full = species_match.group(1)
        # Extract just the species epithet (second word after genus)
        species_words = species_full.split()
        if len(species_words) >= 2:
            species_epithet = species_words[1]
            return f"{genus} {species_epithet}"
        else:
            return f"{genus} sp."
    else:
        return f"{genus} sp."

# --- Process ---
method = "blast"  # blast or sintax
all_hits = []

if method == "blast":
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
            #"-max_target_seqs", "100", # If you want to limit the number of hits, might miss some good hits so use with caution
            "-num_threads", "8"  # Adjust number of threads as needed
        ]

        try:
            # Run BLAST
            if method == "blast":
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
              .sort_values(by=["query", "bitscore"],
                           ascending=[True, False])
              .reset_index(drop=True))
        results_df.to_csv(combined_output, sep="\t", index=False)
        print(f"\nCombined results written to {combined_output}")

        # Save only the top hit per query
        top_hits_df = results_df.groupby("query", as_index=False).head(3)
        top_hits_df.to_csv(top_hits_output, sep="\t", index=False)
        print(f"Top hits written to {top_hits_output}")

        # Only keep the query and taxa and save as tsv
        top_hits_df = results_df.groupby("query", as_index=False).first()[["query", "taxa"]]
        top_hits_df.to_csv(top_hits_compact_output, sep="\t", index=False)
        print(f"Top hits (taxa only) written to {top_hits_compact_output}")

        # Create annotated FASTA file with top hit species in headers
        create_annotated_fasta()
        print(f"Annotated FASTA written to {annotated_fasta_output}")
    else:
        print("No results collected.")
elif method == "sintax":
    # Run sintax classification
    vsearch_cmd = [
            "vsearch",
            "--sintax", str(input_fasta),
            "--db", str(vsearch_db),
            "--tabbedout", str(sintax_combined_output),
            "--threads", "4"
        ]
    subprocess.run(vsearch_cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
else:
    print("Invalid method specified. Use 'blast' or 'sintax'.")
    create_annotated_fasta()