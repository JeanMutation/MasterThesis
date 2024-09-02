# MasterThesis
This repository conains all important files and scripts mentioned and used in my master thesis.

The python scripts were run with the python version 3.11. The R-Scripts were run with the R-version 4.3.1.

## Diversity analysis
Diversity analysis is performed on the constructed OTU tables. 
There are scripts availabe for:
- alpha-diversity
- top 10 OTU calculation
- beta-diversity

## *Albugo* and *Dioszegia* interactions
Here, the process can be performed in the following order:
1. Create input files (with the wanted filtering parameter and fitted to the format of Fastspar)
2. Run Fastspar (with the provided bash scripts)
3. Create filtered networks from Fastspar data using the after processing script
4. evtl. create gephi visualization (dataframes for gephi are provided)

In addition, hub analysis was performed on the networks (Script: hub_analysis.R)

## Keystone analysis

The scripts for the keystone analysis contain the following workflow:
![image](https://github.com/user-attachments/assets/cc6c5672-0d4e-461a-b3c3-8e1a7267cba1)


## Genome assembly of *A. laibachii* Nc14
Most of the scripts in this process are from the *"Beginners tutorial for genome assembly"* by Kim and Kim (2022)

Link: https://star-protocols.cell.com/protocols/1799#fig2

The scripts for the genome assembly contain the following workflow:
![image](https://github.com/user-attachments/assets/f06df118-d818-4f56-89d9-9648bea8c92f)

The data from assembly, the repeat masking, the genome annotation and the functional annotation are provided.
Due to the big size, the raw reads are provided [here](https://drive.google.com/file/d/1oAYSxDM1tzf0sKEnzZ8eniy8DEMaPR7L/view?usp=sharing)

