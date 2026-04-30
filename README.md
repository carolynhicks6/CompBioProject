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
<img width="2868" height="1094" alt="ChaiConnectRuntime" src="https://github.com/user-attachments/assets/8420ecb2-819e-4d6f-aeae-a92ea6789ffa" />

In the Chai script, there is a fasta-format section where the sequences for the LqhIII toxin and the candidate protein motifs can be inserted. To run the script, use the following input data: 

LqhIII toxin: VRDGYIAQPENCVYHCFPGSSGCDTLCKEKGGTSGHCGFKVGHGLACWCNALPDNVGIIVEGEKCHS
ABLIM1 motif: CKGEVLRVQTKHFHIKCFTCKVCGCDLAQGGFFIKNGEYLCTLDYQRMYGTRC

The motif can be substituted for other proteins and other motifs. Each line in the file results/structure/Chai_Lab_Inputs.txt represents a unique muscle protein motif that can be plugged into the "peptide" section.

Example shown below: 
<img width="1212" height="484" alt="ExampleFastaRevised" src="https://github.com/user-attachments/assets/16e82e9f-ad4c-41d8-8842-6a94477e9b65" />


Then, run each cell by clicking the play buttons in each cell. If you encounter errors, restart the session, reconnect to the runtime, and run the cells again. 

After running code, the resulting structural prediction files (.cif format) can be downloaded from the "Files" tab, indicated by a folder icon on the left side of the google colab interface for the chai-1 script. Example shown below:
<img width="2866" height="1346" alt="ChaiCIF" src="https://github.com/user-attachments/assets/ad3189b6-ea9c-4e66-8380-ab247eea077d" />

### PyMOL
As an additional tool for structural visualization, the script pymol_file_generator.py is present in scripts/. If you wish to use this additional functionality for visualization, the first step is to install PyMOL locally on your device. Please refer to the PyMOL documentation for installation: https://github.com/schrodinger/pymol-open-source

After successfully installing PyMOL, navigate to the original project directory and activate the Conda environment.
```
cd CompBioProject/
conda activate motif_pipeline
```
Then, download the .pdb files from data/pdb/ to your local machine. Navigate into the file pymol_file_generator.py and change the file paths to the paths to the newly downloaded .pdb files on your local machine. PyMOL requires local file paths. 
