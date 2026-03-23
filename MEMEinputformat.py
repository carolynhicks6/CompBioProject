from Bio import SeqIO
from Bio.Align import substitution_matrices
import numpy as np

# given fasta files and list of motifs
# generate correct format for inputs

# first, need to update headers of muscle protein fasta file
with open("MEME_proteins_input.fasta", "w") as output:
    for record in SeqIO.parse("muscle_proteins.fasta", "fasta"):
        record.description = ""
        record.id = record.id.replace("old", "new")
        SeqIO.write(record, output, "fasta")

# function to format input for motifs
# edit to include background frequency lines
def write_motif_file(motifs, out_file):
    alphabet = "ACDEFGHIKLMNPQRSTVWY"
    
    with open(out_file, "w") as f:
        # write header
        f.write("MEME version 5\n")
        f.write(f"ALPHABET= {alphabet}\n\n")

        for i, seq in enumerate(motifs):
            # motif identifier and alternate name
            f.write(f"MOTIF {seq} motif_{i+1}\n")

            # letter-probability matrix header
            f.write(f"letter-probability matrix: alength= {len(alphabet)} w= {len(seq)} nsites= 1\n")
            
            # design probability matrix
            # shows the probabilities of each amino acid at each position in motif
            # how to include scoring / chemical group similarity? this would modify code to account for weak motifs
            for i in seq:
                row = []
                for aa in alphabet:
                    if i == aa:
                        row.append("1.0000")
                    else:
                        row.append("0.0000")
                f.write("  " + "  ".join(row) + "\n")

        #blosum62 = substitution_matrices.load("BLOSUM62") 
        #trimmed_matrix = blosum62.select(alphabet)
        #f.write(str(trimmed_matrix))
    
motifs = ["GLALSDLIQKYFF", "LSDLIQ", "LSSLIQ"]
write_motif_file(motifs, "motifs_input.meme")

