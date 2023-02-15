# TranNet
TranNet: A computational method infers the gene-gene interaction across normal and cancer tissues Tumors
#### Pre-required installations before running TranNet
Python libraries math, pandas, numpy and sklearn are prerequired to be installed before running TranNet
#### Input Data to TranNet
Two matrices represnting sample matched gene expression in Normal and Tumor respectively as cvs file (Column names are genes whereas row names are samples) 
We recommend users to have their data standardized for each genes before running the code.
```
nEXP=pd.read_csv('NormalData.csv', index_col=0) #Load Normal expression data
tEXP=pd.read_csv('TumorData.csv', index_col=0)  #Load Tumor expression data
```
#### To Run TranNet
ProjectedGradient procedure takes two matrices described above and return the transition weight matrix as output
```
M=ProjectedGradient(tEXP, nEXP)  #Find optimal transition map M (sparse matrix presenting the interaction between genes across 2 conditions)
```
#### Output of TranNet
Save result as the transitition weight matrix and genes ranked based on their regulatory potential scores
```
(M.T).to_csv('Transition_Matrix_M.csv') # Save transition matrix M as cvs file (Direction of interactions from row to columns in saved matrix) 
RP=(abs(M)).sum(axis=0).sort_values(ascending=False)  #Compute scores of regulatory potentials for genes (Total outgoing effects from each genes)
```
#### Python [Package](code) 
* The TranNet method is implemented in python and the codes are available as [Python Code](code/TranNet.py) and [Jupyter Notebook module](code/TranNet.ipynb).

#### Results of the analysis on five cancer data sets [Result](result).
The following files proves the results for the analysis on BRCA, LUAD, LUSC, PRAD and LIHC cancer data. 
* Lists of the gene and principal components sorted based on their score of regulatory potentials [Regulatory Potentials](result/Genes_Regulatory_Potentials.xlsx). The genes ranked in higher order have more contribution effect on the transition network.
* GO terms enriched for the sorted list of predictor genes [GO terms for predictors](result/GO_terms_Predictor_Genes.xlsx). The enriched GO terms are assumed as the representation of biological processes, molecular functions and cellular components involved in or contributing to the transition.
* GO terms enriched for the sorted list of predicted genes [GO terms for targets](result/GO_terms_Well_Predicted_Genes.xlsx). The enriched GO terms are assumed as the representation of biological processes, molecular functions and cellular components influenced by the transition.
