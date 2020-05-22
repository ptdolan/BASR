# BASR
Pipeline for Ancestral State Reconstruction


Analysis Tools for Initial Bayesian Ancestral State Reconstruction

## Workflow:

### 1. Alignment

BASR/Scripts/BasicAlignment.R - NT align sequence with DECIPHER and generate tree in initial ML tree (guide tree?) in FastTree. 

BASR/Scripts/CodonAlignment.R - Codon align CDS with DECIPHER and generate tree in initial ML tree (guide tree?) in FastTree. 

### 2. BEAST analysis - in BEAST software

### 3. Analysis of trees file from Bayesian phylogenetic analysis
ASR_TRiC2020.ipynb - iPython notebook that parses tree and collects ancestral state sequences, does basic analysis on the sequence strings. 

### 4. Visualization of Data. 
ASRAnalysis.R - Visualization and stats on ASR data file (.csv from ipynb) 
 
