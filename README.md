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
As result, save the transitition weight matrix and list of genes ranked based on their regulatory potential scores
```
(M.T).to_csv('Transition_Matrix_M.csv') # Save transition matrix M as cvs file (Direction of interactions from row to columns in saved matrix) 
RP=(abs(M)).sum(axis=0).sort_values(ascending=False)  #Compute scores of regulatory potentials for genes (Total outgoing effects from each genes)
```
#### Python [Package](code) 
* The TranNet method is implemented in python and the codes are available as [Python Code](code/TranNet.py) and [Jupyter Notebook](code/TranNet.ipynb) modules.

#### Data sets for five cancers [Data](data)
For each of the five cancer data sets (BRCA, LUAD, LUSC, PRAD and LIHC), the gene expression matrix is standardized with z-score for each gene (columns) in control and tumor samples respectively. For the gene expression for each cancer, first half of samples (first half of rows in the gene expression matrix are control samples) represents control samples while the last half represents tumor samples (last half of rows in the gene expression matrix are tumor samples).

#### Results of the analysis on five cancer data sets [Result](result)
The following files proves the results for the analysis on BRCA, LUAD, LUSC, PRAD and LIHC cancer data. 
* Lists of the gene and principal components sorted based on their score of regulatory potentials [Regulatory Potentials](result/Genes_Regulatory_Potentials.xlsx). The genes ranked in higher order have more contribution effect on the transition network.
* GO terms enriched for the sorted list of predictor genes [GO terms for predictors](result/GO_terms_Predictor_Genes.xlsx). The enriched GO terms are assumed as the representation of biological processes, molecular functions and cellular components involved in or contributing to the transition.
* GO terms enriched for the list of the genes sorted by their prediction errors [GO terms for well predicted genes](result/GO_terms_Well_Predicted_Genes.xlsx). GO term biological processes, molecular functions and cellular components are provided in the file.
* KEGG pathways enriched for the target genes of the top 10 predictors [KEGG pathways for targets of top 10 predictors](result/KEGG_Pathways_Top-10-Genes_Targets.xlsx)
* Lists of genes whose expression in normal samples is identified as the predictor its corresponding expression in tumor [Self-Predictors](result/Self-Predictors.xlsx). GO term biological processes, molecular function and cellular components enriched for the self-predictor genes are provided[GO terms for Self-Predictors](result/GO_terms_Self-Predictor_Genes.xlsx)
