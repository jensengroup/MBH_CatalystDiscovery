# MBH_CatalystDiscovery

This repository includes code for the workflow and data analysis described in 

* Maria H. Rasmussen, Julius Seumer, and Jan H. Jensen. "De Novo Catalyst Discovery: Fast Identification of New Catalyst
Candidates", _ChemRxiv_ (June 2023)


## Meta-MD worflow
Code for running and extracting reactions based on meta-MD simulations can be found in ```metaMD_worfklow```.
The worfklow is set up to run on a cluster using the slurm workload manager by running
```
python control_metadyn_runs.py input_smiles.csv
```
where ```input_smiles.csv``` should contain a column labelled ```smiles``` with the reactant systems represented by mapped SMILES.
Hyperparameters for the meta-MD runs can be changed in ```control_metadyn_runs.py```
This will generate a directory for each reactant system (named by the index value in ```input_smiles.csv```
After the calculations have finished, running 
```
python combine_runs.py input_smiles.csv
```
will generated .csv files containing the found mapped reaction SMILES and count for each reactant system.


## Gibbs free energies for intermediates
The Gibbs free energies for intermediates are obtained using Gaussian 16 - the code for the workflow is found in ```GibbsFreeEnergies_worklow```
The calculations are done on individual molecules/fragments of the intermediates by running
```
python control_conf_search.py fragments.csv
```
```fragments.csv``` should contain a column named ```can_smiles``` containing canonical SMILES for all molecules/fragments 
needed to obtain the intermediate energies. Subsequently the Gibbs free energies for the individual molecules are obtained running 
```
python extract_free_energies.py fragments.csv
```
which also checks for changes in the adjacency matrix and errors in embedding/optimization. This information is collected in a generated .csv file. 


## Data availability

Data generated as part of the work described in the above paper is available [here](https://sid.erda.dk/sharelink/C4RVLJdhC5)

The data repository includes

* .csv files containing the mapped reaction SMILES and count for the three meta-MD iterations (metaMD_reactions)
* .xyz files for fragments optimized at the B3LYP/6-31+G(d,p)/SMD(methanol) level of theory (fragment_energies_dft)
* .xyz file for the transition state structures found (TS_structures) at the same level of theory
* output files for the five genetic algorithm runs (GA_output)

## Data Analysis
The directory ```DataAnalysis``` contains example code for generating networks like those presented in the paper. Here we only look at reactions found 
in the first iteration of meta-MD reactions (information about reaction collected in ```step1```).
* ```nodes_dictionary.csv``` contains the SMILES to integer mapping used to generate the networks. 
* ```node_energies_step1+2+3.csv``` contains energies for the intermediates calculated relative to the original reactant system. The energies are based on the Gibbs free energies 
calculated for individual molecules using the above described code.


## Dependencies

Information about the conda environment used to run the metaMD workflow is contained in ```MBH_environment.yml```
The reaction discovery workflow relies on the meta-dynamics method in the xtb program. Specifically the results in the paper are obtained with xtb version 6.1.
The DFT Gibbs free energies (for both intermediates and transition states) are obtained using Gaussian 16.
For the genetic algorithm runs we use code published as part of Seumer et al., _Angew. Chem.Int. Ed._ **2023**, _62_: [GA GitHub repository](https://github.com/jensengroup/mbh_catalyst_ga) 
