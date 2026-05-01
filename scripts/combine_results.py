import pandas as pd #dataframe
import numpy as np #scoring
import sys #command line
import os #files/paths

#input files
blast_file = sys.argv[1]
meme_file = sys.argv[2]
outdir = sys.argv[3]

os.makedirs(outdir, exist_ok=True)

#load data
blast = pd.read_csv(blast_file, sep="\t")
meme = pd.read_csv(meme_file, sep="\t")

#make sure Nav1.5 isn't there
blast = blast[blast["Protein"] != "Nav1.5"]
meme = meme[meme["Protein"] != "Nav1.5"]

#build blast rows
blast_rows = []
for _, row in blast.iterrows():
    blast_rows.append({
        "Method": "BLAST+MSA",
        "Protein": row["Protein"],
        "Motif": row["Motif"],
        "MSA_Score": row["Alignment_Score"],
        "MSA_Identity": row["Identity(%)"],
        "MSA_Similarity": row["Similarity(%)"],
        "MEME_Pvalue": None,
        "Match": row.get("Aligned_Hit", row.get("Match", ""))
    })

#build meme rows
meme_rows = []
for _, row in meme.iterrows():
    meme_rows.append({
        "Method": "MEME",
        "Protein": row["Protein"],
        "Motif": row["Motif"],
        "MSA_Score": None,
        "MSA_Identity": None,
        "MSA_Similarity": None,
        "MEME_Pvalue": row.get("P_value", None),
        "Match": str(row["Match"])
    })

#combine both methods
combined = pd.DataFrame(blast_rows + meme_rows)

#remove sequences shorter than 4 AAS
combined = combined[combined["Match"].str.len() >= 4]

#print output for check
print("Total rows:", len(combined))
print(combined["Method"].value_counts())

#unique matches per protein
simple = (
    combined.groupby("Protein")["Match"]
    .apply(lambda x: list(set(x.dropna())))
    .reset_index()
)
simple.to_csv(f"{outdir}/simple_unique_hits.tsv", sep="\t", index=False)

#scoring per unique found sequence
def score_sequence(group):
    composite = 0.0
    components = []

    #meme scoring
    meme_rows = group[group["Method"] == "MEME"]
    if len(meme_rows) > 0:
        pval = meme_rows.iloc[0]["MEME_Pvalue"]
        if pd.notna(pval) and pval > 0:
            s = -np.log10(pval) * 0.4
            composite += s
            components.append(f"MEME:{s:.2f}")

    #BLAST scoring
    msa_rows = group[group["Method"] == "BLAST+MSA"]
    if len(msa_rows) > 0:
        msa_score = msa_rows.iloc[0]["MSA_Score"]
        msa_id = msa_rows.iloc[0]["MSA_Identity"]
        if pd.notna(msa_score):
            s = (msa_score / 10) * 0.3
            composite += s
            components.append(f"MSA:{msa_score}")
        if pd.notna(msa_id):
            s = (msa_id / 100) * 0.2
            composite += s
            components.append(f"ID:{msa_id}%")

    return pd.Series({
        "Composite_Score": round(composite, 3),
        "Num_Methods": group["Method"].nunique(),
        "Score_Components": "; ".join(components),
        "Motifs_Matched": "; ".join(group["Motif"].dropna().unique())
    })

#group by protein + found sequence instead of protein + motif
scores = (
    combined.groupby(["Protein", "Match"])
    .apply(score_sequence)
    .reset_index()
    .sort_values("Composite_Score", ascending=False)
)

# keep top 5 per protein
scores = scores.groupby("Protein").head(5)
#save output
scores.to_csv(f"{outdir}/ranked_results.tsv", sep="\t", index=False)