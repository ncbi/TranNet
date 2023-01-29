# TranNet
TranNet: A computational method infers the gene-gene inteaction across normal and cancer tissues Tumors

#### Pre-required installations before running TranNet

Python libraries math, pandas, numpy and sklearn are prerequired to be installed before running TranNet

#### Input to TranNet

Gene expression data matrices in Normal and Tumor respectively, 

```
nEXP=pd.read_csv('NormalData.csv', index_col=0) #Load Normal expression data
tEXP=pd.read_csv('TumorData.csv', index_col=0)  #Load Tumor expression data

```
