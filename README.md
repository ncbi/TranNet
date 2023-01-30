# TranNet
TranNet: A computational method infers the gene-gene inteaction across normal and cancer tissues Tumors

#### Pre-required installations before running TranNet

Python libraries math, pandas, numpy and sklearn are prerequired to be installed before running TranNet

#### Input Data to TranNet

Sample matched gene expression matrices in Normal and Tumor respectively as cvs file (Column names are genes whereas row names are samples) 

```
nEXP=pd.read_csv('NormalData.csv', index_col=0) #Load Normal expression data
tEXP=pd.read_csv('TumorData.csv', index_col=0)  #Load Tumor expression data
```
#### Run TranNet
```
M=ProjectedGradient(tEXP, nEXP)  #Find optimal transition map M (sparse matrix presenting the interaction between genes across 2 conditions)

```
#### Output of TranNet
```
(M.T).to_csv('Transition_Matrix_M.csv') # Save transition matrix M as cvs file (Direction of interactions from row to columns in saved matrix) 
RP=(abs(M)).sum(axis=0).sort_values(ascending=False)  #Compute scores of regulatory potentials for genes (Total outgoing effects from each genes)
```
#### Results of the analysis on five types of cancer data sets.
The following files proves the results for the analysis on BRCA, LUAD, LUSC, PRAD and LIHC cancer data 

* The TranNet method is implemented in python and the codes are available as [Python Code](TranNet Code/TranNet.py) and [Jupyter Notebook module](TranNet Code/TranNet.ipynb).
