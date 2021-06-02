# IAVreassortment
## About
IAVreassortment is a comprehensive analysis to detect reassortment from IAV genomes. The order of execution is DataProcess, iqtree, ReassortmentDetection.
## Usage
### Data
The sequence data used in the programs was downloaded from https://ftp.ncbi.nih.gov/genomes/INFLUENZA/ (as of October 13, 2020), including the following files:
- genomeset.dat
- influenza.faa
- influenza.cds
- influenza.dat
### Mutiple sequence alighment
The aligned sequences used in the programs were done mafft-7.037-win64。
### Phylogentic tree
The phylogentic trees used in the programs were constructed by IQ-TREE.
<<<<<<< HEAD
## Host classification
We divided the hosts of IAVs into human, swine and avian. The avian host was subdivided into shorebirds, waterfowl, land birds and domestic birds, which was shown in /utils/avianClassification. Specifically, we gave the host classification of each avian influenza virus in /utils/genomeAvianHostsClassification.
## Country classification
We also divided the locations into 22 areas according to the administrative division, including North America, Western Europe and East Asia et al., which was shown in /utils/countryContinetContinentDetail.
