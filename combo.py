##################
# IMPORT PACKAGES
##################
import os
from io import StringIO
import requests
import numpy as np
import pandas as pd
from Bio import Align, SeqIO
from Bio.Align import substitution_matrices
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from pymemesuite.common import MotifFile, Sequence
from pymemesuite.fimo import FIMO

##############################
# SHARED BLAST + MEME SETUP
##############################

#define list of proteins and uniprot accessions
PROTEINS = {
    "Nav1.5": "P15389",
    "ABLIM1": "O14639",
    "MYBPC": "Q14896",
    "MYL2":  "P10916",
}

#uniProt accessions for just muscle proteins (for blast and meme)
MUSCLE_PROTEINS = {k: v for k, v in PROTEINS.items() if k != "Nav1.5"}

#Nav1.5 motifs we will search for
MOTIFS = ["GTVLSDIIQKYFF", "LSDLIQ", "LSSLIQ", "LSDIIQ", "GLALSDLIQKYFF", "FKVGHGLAC", "VNILAKI"] #added in phase 2 motifs


#shared- get sequences from uniprot using provided accession
def fetch_fasta(uniprot_id: str) -> str:
    url = f"https://rest.uniprot.org/uniprotkb/{uniprot_id}.fasta"
    return requests.get(url).text

#shared- download sequences
def load_sequences(protein_dict: dict) -> dict:
    #download sequences for every protein in dict
    seqs = {}
    for name, uid in protein_dict.items():
        fasta_text = fetch_fasta(uid)
        record = SeqIO.read(StringIO(fasta_text), "fasta")
        seqs[name] = str(record.seq)
    return seqs


#download all sequences once (shared between methods)
sequences = load_sequences(PROTEINS)

#shared function
def find_motif_positions(query_seq: str, motifs: list) -> list:
    #locate each motif within nav1.5
    #return list of dictionaries with motif, start, end
    positions = []
    for motif in motifs:
        idx = query_seq.find(motif)
        if idx != -1:
            positions.append({
                "Motif": motif,
                "Start": idx + 1,           # convert to 1-based
                "End":   idx + len(motif),
            })
    return positions


#pre-compute motif positions in Nav1.5 (used by the blast assign_motif step)
nav_seq = sequences["Nav1.5"]
motif_positions = find_motif_positions(nav_seq, MOTIFS)

#########################
# BLAST PIPELINE
#########################

#writing sequence files for blast

#write the Nav1.5 sequence to its own fasta file
SeqIO.write(
    SeqRecord(Seq(sequences["Nav1.5"]), id="Nav1.5"),
    "Nav1.5.fasta",
    "fasta",
)

#write the muscle-protein sequences to a single multi-FASTA for blast db
muscle_records = [
    SeqRecord(Seq(sequences[name]), id=name)
    for name in MUSCLE_PROTEINS
]
SeqIO.write(muscle_records, "muscle_proteins.fasta", "fasta")

#build blast database and run blastp
os.system("makeblastdb -in muscle_proteins.fasta -dbtype prot -out muscle_db")

os.system(
    "blastp -query Nav1.5.fasta -db muscle_db "
    "-out nav_vs_muscle.tsv "
    "-outfmt '6 qseqid sseqid pident length qstart qend sstart send evalue bitscore'"
)

#load blast results
blast_df = pd.read_csv(
    "nav_vs_muscle.tsv",
    sep="\t",
    names=[
        "Query", "Subject",
        "Percent Identity", "Alignment Length",
        "Query Start", "Query End",
        "Subject Start", "Subject End",
        "E-value", "Bit Score",
    ],
)

#annotate BLAST hits with overlapping Nav1.5 motifs

def assign_motif(row):
    #return the first motif whose position in Nav1.5 overlaps the query alignment span in this BLAST hit
    for mp in motif_positions:
        if (row["Query Start"] <= mp["End"]) and (row["Query End"] >= mp["Start"]):
            return mp["Motif"]
    return "No motif"


blast_df["Motif"] = blast_df.apply(assign_motif, axis=1)

#extract the subject (muscle-protein) subsequence for each hit

def matched_seq(row):
    #slice the subject protein sequence to match the blast alignment region
    subj_seq = sequences[row["Subject"]]
    start = row["Subject Start"] - 1   #convert to 0-based
    end   = row["Subject End"]
    return subj_seq[start:end]


blast_df["Motif Match"] = blast_df.apply(matched_seq, axis=1)
blast_df.to_csv("nav_vs_muscle_results.tsv", sep="\t", index=False)

#pairwise alignment of each motif vs blast hit region

#configure local Smith-Waterman aligner with BLOSUM62
aligner = Align.PairwiseAligner()
aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
aligner.open_gap_score   = -11
aligner.extend_gap_score = -1
aligner.mode             = "local"

blosum62 = substitution_matrices.load("BLOSUM62")
msa_rows = []

for _, row in blast_df.iterrows():
    subject_name  = row["Subject"]
    blast_hit_seq = row["Motif Match"]

    for motif in MOTIFS:
        #align short motif against the BLAST-hit region
        alignment = next(iter(aligner.align(motif, blast_hit_seq)))

        aligned_motif = alignment[0]
        aligned_hit   = alignment[1]

        #count aligned (non-gap) positions
        comparable   = sum(1 for a, b in zip(aligned_motif, aligned_hit) if a != "-" and b != "-")
        identities   = sum(a == b for a, b in zip(aligned_motif, aligned_hit) if a != "-" and b != "-")
        similarities = sum(
            1 for a, b in zip(aligned_motif, aligned_hit)
            if a != "-" and b != "-" and blosum62.get((a, b), -99) > 0
        )

        pct_identity   = round(100.0 * identities   / comparable, 1) if comparable else 0.0
        pct_similarity = round(100.0 * similarities / comparable, 1) if comparable else 0.0

        msa_rows.append({
            "Subject":         subject_name,
            "BLAST_Hit_Region": blast_hit_seq,
            "Motif":           motif,
            "Aligned_Motif":   aligned_motif,
            "Aligned_Hit":     aligned_hit,
            "Identity(%)":     pct_identity,
            "Similarity(%)":   pct_similarity,
            "Alignment_Score": round(alignment.score, 1),
        })

msa_df = pd.DataFrame(msa_rows)
msa_df.to_csv("msa_results.tsv", sep="\t", index=False)

##########################
# MEME SUITE PIPELINE!
#########################

#write muscle-protien sequences to fasta for FIMO (using sequences dict we already have)
fasta_blocks = []
for name, uid in MUSCLE_PROTEINS.items():
    fasta_blocks.append(fetch_fasta(uid))   # raw FASTA text (header + sequence)

#join blocks with a single newline so FIMO sees one valid multi-FASTA
with open("muscle_seq.fasta", "w") as f:
    f.write("\n".join(fasta_blocks))


#############################
# BUILDING MEME MOTIF FILE
#############################
#each position in the motif gets probability 1.0 for the actual residue, 0.0 for everything else - overridden later by substitution-aware version, but left here as the original

def write_motif_file(motifs, out_file):
    """Write a MEME-format motif file with strict identity probability matrices."""
    alphabet = "ACDEFGHIKLMNPQRSTVWY"

    with open(out_file, "w") as f:
        f.write("MEME version 5\n")
        f.write(f"ALPHABET= {alphabet}\n\n")

        for i, residue in enumerate(motifs):
            f.write(f"MOTIF {residue} motif_{i+1}\n")
            f.write(
                f"letter-probability matrix: alength= {len(alphabet)} "
                f"w= {len(residue)} nsites= 1\n"
            )
            #one row per position in the motif
            for aa_pos in residue:
                row = ["1.0000" if aa_pos == aa else "0.0000" for aa in alphabet]
                f.write("  " + "  ".join(row) + "\n")


write_motif_file(MOTIFS, "motifs_input.meme")

########################
# MEME SUITE FIMO
########################

#parse muscle protein sequences into PyMEMEsuite Sequence objects
fimo_sequences = [
    Sequence(str(record.seq), name=record.id.encode())
    for record in SeqIO.parse("muscle_seq.fasta", "fasta")
]

#keep a plain dict for slicing subsequences later
seq_id_dict = {
    record.id: str(record.seq)
    for record in SeqIO.parse("muscle_seq.fasta", "fasta")
}

fimo = FIMO(both_strands=False)

#store FIMO hits in a list of dicts so we can build a DataFrame afterwards
meme_results_raw = []

#dict mapping protein id -> (start, stop) of its last FIMO hit (same as original script)
sequence_pts = {}
pdb_names_ids = []

with open("meme_suite_output.txt", "w") as out:

    with MotifFile("motifs_input.meme") as motif_file:

        for motif in motif_file:
            out.write(f"{motif.name.decode()} {motif.consensus}\n")
            header = "Accession\tstart\tstop\tstrand\tscore\tpvalue\tmatch_seq\n"
            out.write(header)

            #score every sequence against this motif
            pattern = fimo.score_motif(motif, fimo_sequences, motif_file.background)

            for m in pattern.matched_elements:
                match_name    = m.source.accession.decode()
                full_seq      = seq_id_dict[match_name]
                sub_seq       = full_seq[m.start - 1 : m.stop]

                line = (
                    f"{match_name}\t{m.start}\t{m.stop}\t"
                    f"{m.strand}\t{m.score}\t{m.pvalue}\t{sub_seq}\n"
                )
                out.write(line + "\n")

                #collect for combined analysis DataFrame (NEW)
                meme_results_raw.append({
                    "Method":        "MEME",
                    "Protein":       match_name,
                    "Motif":         motif.name.decode(),
                    "Start":         m.start,
                    "End":           m.stop,
                    "Score":         m.score,
                    "P_value":       m.pvalue,
                    "Match_Sequence": sub_seq,
                })

                #update position dict and generate PyMOL scripts (unchanged)
                sequence_pts[match_name] = (m.start, m.stop)
                sequence_keys = list(sequence_pts.keys())

                for protein_id in sequence_keys:
                    #ABLIM1 — no public PDB; use AlphaFold model downloaded manually - put this file in the github?
                    if protein_id == "sp|O14639|ABLM1_HUMAN":
                        pymol_script_1 = [
                            "load /Users/siennad0304/Documents/Senior Year/Comp Biology /AF-O14639-F1-model_v6.pdb",
                            "as cartoon",
                            f"select i. {m.start}-{m.stop}",
                            "set_name sele, motif",
                            "hide everything",
                            "show sticks, motif",
                            "fetch 6zu0",
                            "align 6zu0, motif",
                            "orient motif",
                        ]
                        with open(f"{protein_id}.pml", "w") as l:
                            l.write("\n".join(pymol_script_1))

                    # MYBPC — PDB 1PD6
                    if protein_id == "sp|Q14896|MYPC3_HUMAN":
                        pymol_script_2 = [
                            "fetch 1PD6",
                            "as cartoon",
                            f"select i. {m.start}-{m.stop}",
                            "set_name sele, motif",
                            "hide everything",
                            "show sticks, motif",
                            "fetch 6zu0",
                            "align 6zu0, motif",
                            "orient motif",
                        ]
                        with open(f"{protein_id}.pml", "w") as x:
                            x.write("\n".join(pymol_script_2))

                    # MYL2 — PDB 5TBY
                    if protein_id == "sp|P10916|MLRV_HUMAN":
                        pymol_script_3 = [
                            "fetch 5TBY",
                            "as cartoon",
                            f"select i. {m.start}-{m.stop}",
                            "set_name sele, motif",
                            "hide everything",
                            "show sticks, motif",
                            "fetch 6zu0",
                            "align 6zu0, motif",
                            "orient motif",
                        ]
                        with open(f"{protein_id}.pml", "w") as z:
                            z.write("\n".join(pymol_script_3))


#########################
# COMBINED RESULTS ANALYSIS
##########################

#merge blast and meme results into single dataframe, rank protein-motif combinations with composite score

#format BLAST+MSA results into the shared schema
#each row represents one (protein, motif) alignment from alignment step
blast_msa_results = []
for _, row in msa_df.iterrows():
    blast_msa_results.append({
        "Method":         "BLAST+MSA",
        "Protein":        row["Subject"],
        "Motif":          row["Motif"],
        "Start":          None,   # BLAST doesn't give exact motif-level positions
        "End":            None,
        "MEME_Score":     None,
        "MEME_Pvalue":    None,
        "MSA_Score":      row["Alignment_Score"],
        "MSA_Identity":   row["Identity(%)"],
        "MSA_Similarity": row["Similarity(%)"],
        "Match_Sequence": row["Aligned_Hit"].replace("-", ""),   # strip alignment gaps
    })

#format MEME/FIMO results into the shared structure
#meme_results_raw populated during the FIMO loop above

meme_results_formatted = []
for result in meme_results_raw:
    meme_results_formatted.append({
        "Method":         result["Method"],
        "Protein":        result["Protein"],
        "Motif":          result["Motif"],
        "Start":          result["Start"],
        "End":            result["End"],
        "MEME_Score":     result["Score"],
        "MEME_Pvalue":    result["P_value"],
        "MSA_Score":      None,
        "MSA_Identity":   None,
        "MSA_Similarity": None,
        "Match_Sequence": result["Match_Sequence"],
    })

#put both method outputs into one dataframe
all_results = blast_msa_results + meme_results_formatted
combined_df = pd.DataFrame(all_results)

#get full combined table for downstream use
combined_df.to_csv("all_methods_combined.tsv", sep="\t", index=False)


#compute a composite score and rank protein-motif pairs

def identify_best_matches(df: pd.DataFrame) -> pd.DataFrame:
    #checking every unique protein + motif pair, combine score from whichever methods detected a hit a produce a single composite score
    #weights for composite score: 
        #MEME  : -log10(p-value) × 0.40   — statistical significance
        #alignment score / 10 × 0.40  — structural fit
        #alignment identity (%) / 100   × 0.20  — sequence identity
    
    best_matches = []

    #iterate over every unique protein-motif combination found by any method
    for (protein, motif), group in df.groupby(["Protein", "Motif"]):

        composite_score  = 0.0
        score_components = []

        #MEME FIMO
        #use lowest p-value if multiple FIMO hits exist for pair
        meme_rows = group[group["Method"] == "MEME"]
        if not meme_rows.empty:
            best_pval = meme_rows["MEME_Pvalue"].dropna().min()
            if pd.notna(best_pval) and best_pval > 0:
                meme_contrib = -np.log10(best_pval) * 0.40
                composite_score += meme_contrib
                score_components.append(f"MEME: {meme_contrib:.2f}")
            elif pd.notna(best_pval) and best_pval == 0:
                # p-value of exactly 0 treated as extremely significant
                composite_score += 10 * 0.40
                score_components.append("MEME: 4.00 (p=0)")

        #alignment
        #use the best alignment score if multiple blast hits exist for pair
        msa_rows = group[group["Method"] == "BLAST+MSA"]
        if not msa_rows.empty:
            best_msa_score    = msa_rows["MSA_Score"].dropna().max()
            best_msa_identity = msa_rows["MSA_Identity"].dropna().max()

            if pd.notna(best_msa_score):
                msa_contrib = (best_msa_score / 10) * 0.40
                composite_score += msa_contrib
                score_components.append(f"MSA_score: {msa_contrib:.2f}")

            if pd.notna(best_msa_identity):
                id_contrib = (best_msa_identity / 100) * 0.20
                composite_score += id_contrib
                score_components.append(f"ID: {id_contrib:.2f}")

        #collect unique non-null match sequences across all methods
        unique_seqs = "; ".join(
            seq for seq in group["Match_Sequence"].dropna().unique() if seq
        )

        best_matches.append({
            "Protein":          protein,
            "Motif":            motif,
            "Methods_Used":     ", ".join(group["Method"].unique()),
            "Num_Methods":      group["Method"].nunique(),
            "Composite_Score":  round(composite_score, 3),
            "Score_Components": "; ".join(score_components),
            "Match_Sequences":  unique_seqs,
        })

    result_df = pd.DataFrame(best_matches).sort_values(
        "Composite_Score", ascending=False
    ).reset_index(drop=True)

    return result_df


best_matches_df = identify_best_matches(combined_df)

#print summary and save outputs

print(f"\nFound {len(best_matches_df)} unique protein-motif combinations")
print("\nTop 10 best overall matches:")
print(best_matches_df.head(10).to_string(index=False))

best_matches_df.to_csv("best_overall_matches.tsv", sep="\t", index=False)

print("\nANALYSIS COMPLETE!")
print("Results saved to:")
print("nav_vs_muscle_results.tsv (raw BLAST hits)")
print("msa_results.tsv (pairwise alignments)")
print("meme_suite_output.txt (FIMO motif scan)")
print("all_methods_combined.tsv (every hit from all methods)")
print("best_overall_matches.tsv (ranked composite scores)")

#to run
# cd CompBioProject
# python combined_motif_analysis.py