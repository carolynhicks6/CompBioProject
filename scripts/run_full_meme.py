import sys
import pandas as pd
from Bio import SeqIO
import os

fasta_file = sys.argv[1]
motif_file = sys.argv[2]
msa_file = sys.argv[3]
output_file = sys.argv[4]

os.makedirs("results/meme", exist_ok=True)

#load data

motifs = [line.strip() for line in open(motif_file)]

sequences = {
    rec.id: str(rec.seq)
    for rec in SeqIO.parse(fasta_file, "fasta")
}

msa_df = pd.read_csv(msa_file, sep="\t")

#create strict motifs

def write_strict(motifs, out):
    alphabet = "ACDEFGHIKLMNPQRSTVWY"
    with open(out, "w") as f:
        f.write("MEME version 5\n")
        f.write(f"ALPHABET= {alphabet}\n\n")

        for i, motif in enumerate(motifs):
            f.write(f"MOTIF {motif} motif_{i+1}\n")
            f.write(f"letter-probability matrix: alength=20 w={len(motif)} nsites=1\n")

            for aa in motif:
                row = ["1.0" if aa == x else "0.0" for x in alphabet]
                f.write(" ".join(row) + "\n")

write_strict(motifs, "results/meme/strict.meme")

#build flexible

def build_flexible(msa_df):
    """Position-aware flexible dict matching old make_filtered_dict logic."""
    motif_dict = {}

    for _, row in msa_df.iterrows():
        motif = row["Motif"]
        aligned_motif = row["Aligned_Motif"]
        aligned_hit = row["Aligned_Hit"]

        if len(aligned_motif) < 4:
            continue

        motif_dict.setdefault(motif, {})[aligned_motif] = aligned_hit

    return motif_dict


def make_prob_matrix(motif, filtered_dict):
    alphabet = "ACDEFGHIKLMNPQRSTVWY"
    count_matrix = pd.DataFrame(0, index=range(len(motif)), columns=list(alphabet))

    # seed with the motif itself 5x more than substitutions
    for i, aa in enumerate(motif):
        if aa in count_matrix.columns:
            count_matrix.loc[i, aa] += 5

    # update positions based on aligned fragments
    for fragment, hit in filtered_dict.items():
        start = motif.find(fragment)
        if start != -1:
            for i, char in enumerate(hit):
                position = start + i
                if position < len(motif) and char in alphabet:
                    count_matrix.loc[position, char] += 1

    row_sums = count_matrix.sum(axis=1)
    return count_matrix.divide(row_sums, axis=0).round(6)


def write_flexible(motif_dict, out):
    alphabet = "ACDEFGHIKLMNPQRSTVWY"

    with open(out, "w") as f:
        f.write("MEME version 5\n")
        f.write(f"ALPHABET= {alphabet}\n\n")

        for i, (motif, filtered_dict) in enumerate(motif_dict.items()):
            probs = make_prob_matrix(motif, filtered_dict)

            f.write(f"MOTIF {motif} motif_{i+100}\n")
            f.write(f"letter-probability matrix: alength=20 w={len(motif)} nsites=1\n")
            for _, r in probs.iterrows():
                f.write(" ".join(f"{v:.6f}" for v in r) + "\n")
            f.write("\n")

flex_dict = build_flexible(msa_df)
write_flexible(flex_dict, "results/meme/flexible.meme")

import subprocess

#run CLI version of fimo

def run_fimo(meme_file, label):
    out_dir = f"results/meme/fimo_{label}"
    
    os.makedirs(out_dir, exist_ok=True)

    #running fimo on command line
    subprocess.run([
        "fimo",
        "--oc", out_dir,
        meme_file,
        fasta_file
    ], check=True)

    #load results
    fimo_file = os.path.join(out_dir, "fimo.tsv")

    df = pd.read_csv(fimo_file, sep="\t", comment="#")

    #rename columns to match pipeline
    df = df.rename(columns={
        "sequence_name": "Protein",
        "start": "Start",
        "stop": "End",
        "score": "Score",
        "p-value": "P_value",
        "matched_sequence": "Match_Sequence",
        "motif_id": "Motif"
    })

    return df

#run both motif sets (strict and flexible)
df_strict = run_fimo("results/meme/strict.meme", "strict")
df_flex = run_fimo("results/meme/flexible.meme", "flex")

#c
df = pd.concat([df_strict, df_flex])
df["Method"] = "MEME"

df = df[df["P_value"] < 1e-4]

df = df[df["Protein"] != "Nav1.5"]

df = df.rename(columns={"Match_Sequence": "Match"})
df = df[["Protein", "Motif", "Match", "P_value"]]   # keep P_value
df.to_csv(output_file, sep="\t", index=False)

df.to_csv(output_file, sep="\t", index=False)

print(f"MEME hits: {len(df)}")