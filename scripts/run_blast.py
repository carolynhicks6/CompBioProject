import pandas as pd
from Bio import SeqIO, Align
from Bio.Align import substitution_matrices
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import os
import sys

fasta_file = sys.argv[1]
motif_file = sys.argv[2]
outdir = sys.argv[3]

os.makedirs(outdir, exist_ok=True)

sequences = {rec.id: str(rec.seq) for rec in SeqIO.parse(fasta_file, "fasta")}

nav_seq = sequences["Nav1.5"]
muscle = {k: v for k, v in sequences.items() if k != "Nav1.5"}

SeqIO.write(SeqRecord(Seq(nav_seq), id="Nav1.5"), f"{outdir}/nav.fasta", "fasta")
SeqIO.write(
    [SeqRecord(Seq(v), id=k) for k, v in muscle.items()],
    f"{outdir}/muscle.fasta",
    "fasta"
)

os.system(f"makeblastdb -in {outdir}/muscle.fasta -dbtype prot -out {outdir}/db")
os.system(f"""
blastp -query {outdir}/nav.fasta -db {outdir}/db \
-out {outdir}/blast.tsv \
-outfmt "6 qseqid sseqid pident length qstart qend sstart send evalue bitscore"
""")

blast_df = pd.read_csv(
    f"{outdir}/blast.tsv",
    sep="\t",
    names=[
        "Query","Subject","Percent Identity","Alignment Length",
        "Query Start","Query End","Subject Start","Subject End",
        "E-value","Bit Score"
    ]
)

motifs = [line.strip() for line in open(motif_file)]

motif_positions = []
for m in motifs:
    idx = nav_seq.find(m)
    if idx != -1:
        motif_positions.append((m, idx, idx+len(m)))

def assign_motif(row):
    for m, start, end in motif_positions:
        if row["Query Start"] <= end and row["Query End"] >= start:
            return m
    return "None"

blast_df["Motif"] = blast_df.apply(assign_motif, axis=1)

def get_seq(row):
    seq = sequences[row["Subject"]]
    return seq[row["Subject Start"]-1:row["Subject End"]]

blast_df = blast_df[blast_df["Subject"] != "Nav1.5"]

blast_df["Match"] = blast_df.apply(get_seq, axis=1)
blast_df.to_csv(f"{outdir}/blast_results.tsv", sep="\t", index=False)

#alignment with gap penaltys + identity/similarity
aligner = Align.PairwiseAligner()
aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
aligner.open_gap_score = -11
aligner.extend_gap_score = -1
aligner.mode = "local"

blosum62 = substitution_matrices.load("BLOSUM62")
rows = []

for _, row in blast_df.iterrows():
    subject_name = row["Subject"]
    blast_hit_seq = row["Match"]

    for motif in motifs:
        aln = next(iter(aligner.align(motif, blast_hit_seq)))

        aligned_motif = aln[0]
        aligned_hit = aln[1]

        comparable = sum(1 for a, b in zip(aligned_motif, aligned_hit) if a != "-" and b != "-")
        identities = sum(a == b for a, b in zip(aligned_motif, aligned_hit) if a != "-" and b != "-")
        similarities = sum(
            1 for a, b in zip(aligned_motif, aligned_hit)
            if a != "-" and b != "-" and blosum62.get((a, b), -99) > 0
        )
        pct_identity = round(100.0 * identities / comparable, 1) if comparable else 0.0
        pct_similarity = round(100.0 * similarities / comparable, 1) if comparable else 0.0

        rows.append({
            "Protein": subject_name,
            "BLAST_Hit_Region": blast_hit_seq,
            "Motif": motif,
            "Aligned_Motif": aligned_motif,      # restored
            "Aligned_Hit": aligned_hit,           # restored
            "Identity(%)": pct_identity,          # restored
            "Similarity(%)": pct_similarity,      # restored
            "Alignment_Score": round(aln.score, 1)
        })

msa_df = pd.DataFrame(rows)

msa_df = msa_df[msa_df["Protein"] != "Nav1.5"]

msa_df.to_csv(f"{outdir}/msa.tsv", sep="\t", index=False)