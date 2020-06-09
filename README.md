# Primer_TaxoCheck
In sillico validation of primers for diverses taxonomic groups

## General information


To use this repositery, clone it!

## PrimerTree

PrimerTree is a package developped by Jim Hester (https://github.com/jimhester/primerTree). This script is a wrapper around diverses functions in this package, inspired by the one of Bylemans et al. (2018) Ecology and Evolution 
 https://onlinelibrary.wiley.com/doi/full/10.1002/ece3.4387

### How to run PrimerTree

1. Download clustal-omega : http://www.clustal.org/omega/#Download, and note the path to the executable file, you will need this information.

2. If you want it to run faster, an NCBI API key can be added, see https://ncbiinsights.ncbi.nlm.nih.gov/2017/11/02/new-api-keys-for-the-e-utilities/
 to get one.

3. Open the PrimersCheck_PrimerTree.R script update the *Path* and *Data* section.

You will need the "00_Data/Primers.xlsx" excel sheet, with Primer, Group, Sequence.F, Sequence.R columns filled. 


### Note

The first part of PrimerTree is to use NCBI Primer Blast (https://www.ncbi.nlm.nih.gov/tools/primer-blast/). The database used is nr, without environmental data, and Organism is " ". Note that the maximum primer length is 36.

Default parameters:

num_permutations = 25 to limit the number of tested degenerated primer paires



## PrimerMiner

Primer miner is a package develop by Vasco Elbrecht (https://github.com/VascoElbrecht/PrimerMiner) intend to evaluate primer in sillico. This script is a wrapper around diverse function in We write a script to ease the use of this package.

### How to run PrimerMiner

1. Create an aligment (.fasta) of sequences for each locus taxonomic group of interest. This can be performed with the script **01_Code/Subset_TaxoGroup.R**.

### Reference :
https://besjournals.onlinelibrary.wiley.com/doi/full/10.1111/2041-210X.12687

