import Bio.SeqIO
from Bio import SeqIO
from Bio.Seq import Seq
from io import StringIO


sequence_list = []
ablim_list = []
myl2_list = []
mybpc_list = []

for record in Bio.SeqIO.parse("muscle_seq.fna", "fasta"):
     #sequence order: ABLIM1, MYBPC, MYL2
    sequence_list.append(str(record.seq))

#saving muscle sequences as variables for later 
ablim1_seq = sequence_list[0]
mybpc_seq = sequence_list[1]
myl2_seq = sequence_list[2]


#calling the generated best motifs file: best_overall_matches 
motif_dic = {}
remove_list = "[],'"

with open ("simple_unique_hits.tsv") as f: 
    for i in f:
        var = i.strip()
        x = str.maketrans("", "", remove_list)
        var = var.translate(x)
        motif_var = var.split()
        muscle_name = motif_var[0]
        motif = (motif_var[1:])
        motif_dic[muscle_name] = motif 



ablim = motif_dic.get("ABLIM1")
ablim_list = ablim

mybpc = motif_dic.get("MYBPC")
mybpc_list = mybpc

myl2 = motif_dic.get("MYL2")
myl2_list = myl2

ablim_list = [m for m in motif_dic.get("ABLIM1", []) if m]
mybpc_list = [m for m in motif_dic.get("MYBPC", []) if m]
myl2_list = [m for m in motif_dic.get("MYL2", []) if m]

count = 0 
count2 = 0 
count3 = 0 
ablim_chai_motifs = []
mybpc_chai_motifs = []
myl2_chai_motifs = []

for i in range (len(ablim1_seq)): 
            for x in ablim_list:

                if ablim1_seq[i:i+len(x)] == x:
                    ablim_chai_motifs.append(ablim1_seq[i-20:i+len(x)+20])
                     #checking to see if motif is in sequence 
                    start_site = i+1 # start location for motif  
                    end_site = i+len(x) # end location for motif 
                    count += 1
                
                #generating a Script for PyMOl ** the script made can be added to pyMOL 
                    pymol_script_1 = ["load /Users/siennad0304/Documents/Senior Year/Comp Biology /ABLIM1.pdb", 
                                  #your own directory must be used due to visualize strucutres 
                                  # protein PDB files in GitHUb 
                                          "as cartoon", 
                                          (f"select i. {start_site}-{end_site}"), 
                                          "set_name sele, motif",
                                          "color lightblue, motif",
                                          "fetch 6uz0", # sodium channel 
                                          "select i. 1606-1619", #location of known motif in sodium channel 
                                          "set_name sele, navmotif", 
                                          "color green, navmotif",
                                          "hide everything",
                                          "enable navmotif",
                                          "show cartoon, navmotif",
                                          "enable motif",
                                          "show cartoon, motif",
                                          "align navmotif, motif",
                                          "orient navmotif",
                                          "disable motif"]
                
                    with open(f"ABLIM1_{count}.pml", "w") as l:
                        l.write("\n".join(pymol_script_1))
                
#MYBPC3 
for i in range (len(mybpc_seq)): 
        for x in mybpc_list: 

            if mybpc_seq[i:i+len(x)] == x: 
                mybpc_chai_motifs.append(mybpc_seq[i-20:i+len(x)+20])

                start_site = i+1
                end_site = i+len(x)
                count2 += 1

                #generating a Script for PyMOl ** the script made can be added to pyMOL 
                pymol_script_2 = ["load /Users/siennad0304/Documents/Senior Year/Comp Biology /MYBC3.pdb",
                                  #your own directory must be used due to visualize strucutres 
                                  # protein PDB files in GitHUb 
                                          "as cartoon", 
                                          (f"select i. {start_site}-{end_site}"), 
                                          "set_name sele, motif",
                                          "color purple, motif",
                                          "fetch 6uz0", # sodium channel 
                                          "select i. 1606-1619", #location of known motif in sodium channel 
                                          "set_name sele, navmotif", 
                                          "color green, navmotif",
                                          "hide everything",
                                          "enable navmotif",
                                          "show cartoon, navmotif",
                                          "enable motif",
                                          "show cartoon, motif",
                                          "align navmotif, motif",
                                          "orient navmotif",
                                          "disable motif"]
                
                with open(f"MYBPC_{count2}.pml", "w") as x:
                    x.write("\n".join(pymol_script_2))

#MYL2 
for i in range (len(myl2_seq)): 
        for x in myl2_list: 
            if myl2_seq[i:i+len(x)] == x: 
                myl2_chai_motifs.append(myl2_seq[i-20:i+len(x)+20])

                start_site = i+1
                end_site = i+len(x)
                count3 += 1

                # you don't need to download anything PDB number can be used instead 
                pymol_script_3 = ["fetch 5TBY",
                                    "as cartoon", 
                                    (f"select i. {start_site}-{end_site} and chain A+C+E"), 
                                    "set_name sele, motif",
                                    "color orange, motif",
                                    "fetch 6uz0", # sodium channel 
                                    "as cartoon",
                                    "select i. 1606-1619", #location of known motif in sodium channel 
                                    "set_name sele, navmotif", 
                                    "color green, navmotif",
                                    "hide everything",
                                    "enable navmotif",
                                    "show cartoon, navmotif",
                                    "enable motif",
                                    "show cartoon, motif",
                                    "align motif, navmotif",
                                    "orient navmotif",
                                    "disable motif"]

                        
                with open(f"MYL2_{count3}.pml", "w") as z:
                    z.write("\n".join(pymol_script_3))


print(ablim_chai_motifs)
print(mybpc_chai_motifs)
print(myl2_chai_motifs)

with open(f"Chai_Lab_Inputs.txt", "w") as s:
        s.write("\nABLIM1: \n")
        s.write("\n".join(ablim_chai_motifs))
        s.write("\n")

        s.write("\nMYBPC3: \n")
        s.write("\n".join(mybpc_chai_motifs))
        s.write("\n")

        s.write("\nMYL2: \n")
        s.write("\n".join(myl2_chai_motifs))


        