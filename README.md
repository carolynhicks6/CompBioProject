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
When the virtual environment is active, the terminal prompt should show: 
```
(.motif_pipeline) user/path_to_project_directory
```
To run the snakemake file, run the command:
```
snakemake Snakefile --cores 1
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
To determine how well the identified motifs interact with the LqhIII toxin, we used protein modeling and predicting software. Our pipeline creates PyMOL scripts (results/structure/) as well as Chai-1 input sequences (results/structure/Chai_Lab_Inputs.txt) for visualization. An example run of Chai-1 can be shown here:
https://colab.research.google.com/drive/1aFoc2KdYG9In8NDpSrkjNgxc08ujqWI1?usp=sharing.

This can be easily modified for other following comments on the colab script. Protein modeling was done to compare structural similarity between the identified muscle protein motif and the input NAV1.5, as well as to determine binding potential between the identified motifs and LqhIII toxin. 20 amino acids upstream and downstream of the motif in the respective protein were used for visualization purposes. Aggregate confidence scores in Chai-1 were not significantly changed by increasing or decreasing sequence length. 
