## Code Overview
Our project uses a snakefile to run a combination of shell commands and python scripts to fetch protein sequences and then run BLAST, multiple sequence alignment, and the MEME Suite FIMO tool to generate candidate amino acid motif matches from an input list of UniProt IDs and motifs of interest. Using a log-based scoring system, we used identified motifs to generate PyMOL scripts and Chai-1 input sequences to investigate putative structural interaction with the scorpion toxin LqhIII. 

## Project Setup
First, clone this repository into your environment using the following command in the terminal:
```
git clone https://github.com/annalundeen/CompBioProject
```
Then navigate into the project directory: 
```
cd CompBioProject/
```

## Conda Environment and Running Code
This code requires a conda environment for package management. This assumes that your machine has miniconda or anaconda installed previously. Documentation for installing miniconda can be found here (it can be installed globally):
https://www.anaconda.com/docs/getting-started/miniconda/main.

Once miniconda is installed, navigate into your project directory. Then, simply use the following command in the terminal to create the conda environment: 
```
conda env create --environment-spec environment.yml
```
Once the conda environment has been established in your project directory, activate it using the command:
```
conda activate motif_pipeline
```
If you want to deactivat the conda environment, simply use the command:
```
conda deactivate
```
When the virtual environment is active, the terminal prompt should show: 
```
(.motif_pipeline) user/path_to_project_directory
```
To run the snakemake file, run the command:
```
snakemake --cores 1
```
To run everything again, run the clean rule:
```
snakemake clean --cores 1
```

## Input
The initial input file to run the Snakefile is shown in config/input.txt, and it contains a list of proteins, UniProt IDs, and motifs of interest for this project. The inputs and outputs for each individual section of analysis are indicated in the the rules of the Snakefile. 

## Expected Outputs
The Snakefile provided will output two final motif charts in the results folder under final. The first ranked_results.tsv is a table ranking our identified motifs based off P-values and motif sequence length. The second chart, simple_unique_hits.tsv provides a simplified format of the motifs for future analysis. 

## Structural Prediction and Alignment 
### Chai-1
The PKH lab provided a Google Colab script to run Chai-1, a protein prediction and modeling software. The Chai script can be easily modified for other toxins and motifs. 
[Chai-1 Google Colab Script](https://colab.research.google.com/drive/1aFoc2KdYG9In8NDpSrkjNgxc08ujqWI1?usp=sharing) 
Our previously described pipeline outputs Chai-1 input sequences (results/structure/Chai_Lab_Inputs.txt) for visualization.

Protein modeling was done to compare structural similarity between the identified muscle protein motif and the input NAV1.5, as well as to determine binding potential between the identified motifs and LqhIII toxin. 20 amino acids upstream and downstream of the motif in the respective protein were used for visualization purposes. Aggregate confidence scores in Chai-1 were not significantly changed by increasing or decreasing sequence length.

After opening the link to the Google Colab script, the first step is to connect to the specified run-time. In the upper right corner, click on the icon displaying "RAM Disk" and select "Connect to a hosted runtime: T4". Example shown below:
<img width="2872" height="1372" alt="ConnectToRunTime" src="https://github.com/user-attachments/assets/6d9d3549-8d4a-4a62-8f7c-133c889e6ca2" />

When running the following code chunk:
```
!pip install chai-lab
```
An error message may pop up that requires the user to restart the runtime session. Follow this error message and run the code chunk again. It should finish running quickly without error the second time. 

In the Chai script, there is a fasta-format section where the sequences for the LqhIII toxin and the candidate protein motifs can be inserted. To run the script, use the following input data: 

LqhIII toxin: VRDGYIAQPENCVYHCFPGSSGCDTLCKEKGGTSGHCGFKVGHGLACWCNALPDNVGIIVEGEKCHS
ABLIM1 motif: CKGEVLRVQTKHFHIKCFTCKVCGCDLAQGGFFIKNGEYLCTLDYQRMYGTRC

The motif can be substituted for other proteins and other motifs. Each line in the file results/structure/Chai_Lab_Inputs.txt represents a unique muscle protein motif that can be plugged into the "peptide" section.

Example shown below: 
<img width="1212" height="484" alt="ExampleFastaRevised" src="https://github.com/user-attachments/assets/16e82e9f-ad4c-41d8-8842-6a94477e9b65" />


Then, run each cell by clicking the play buttons in each cell. If you encounter errors, restart the session, reconnect to the runtime, and run the cells again. 

After running code, the resulting structural prediction files (.cif format) can be downloaded from the "Files" tab, indicated by a folder icon on the left side of the google colab interface for the chai-1 script. Example shown below:
<img width="2876" height="1378" alt="CIF_Files" src="https://github.com/user-attachments/assets/db5110fe-2bbf-4b74-a12d-f7ea0dc5af41" />

### PyMOL
As an additional tool for structural visualization, the script pymol_file_generator.py is present in scripts/. If you wish to use this additional functionality for visualization, the first step is to install PyMOL locally on your device. Please refer to the PyMOL documentation for installation: https://github.com/schrodinger/pymol-open-source 
Note that you do not need to purchase a license to download pymol, the functions available within the free version is enough for the provided files. 

After successfully installing PyMOL, navigate to the original project directory and activate the Conda environment.
```
cd CompBioProject/
conda activate motif_pipeline
```
Then, download the .pdb files from data/pdb/ to your local machine. Navigate into the file pymol_file_generator.py and change the file paths to the paths to the newly downloaded .pdb files on your local machine. PyMOL requires local file paths to generate structures when PDB codes are not available. Due to incomplete PDB sequence for ABLIM1 and MYBPC3, each of these files will need to be downloaded onto your local machine.

Before running pymol_file_generator.py, two manual changes must be made to ensure you output accurate pymol scripts. Using the previously downloaded muscle .pdb files, copy your local directory under the respective muscle protein section within the script. The location and example path for both ABLIM1 and MYBPC3 are shown below.

#Example ABLIM1 Pathway and Location
<img width="1022" height="474" alt="Screenshot 2026-04-30 at 7 37 50 PM" src="https://github.com/user-attachments/assets/9b000273-8a12-4d56-8a91-bc850f625884" />

#Example MYBPC3 Pathway and Location 
<img width="986" height="344" alt="Screenshot 2026-04-30 at 8 51 39 PM" src="https://github.com/user-attachments/assets/922ec760-856c-432e-8191-9839f4096f2e" />
