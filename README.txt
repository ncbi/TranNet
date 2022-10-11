--------------------------------------
Introduction
--------------------------------------

This README file includes the description of TranNet package,
which take gene expression matrices for normal and tumor tissue samples for cancer patients as an input.
It first constructs a weighted-directed bipartite network representing transition map from normal to tumor expressions and 
compute regulatory potential scores for genes   

References:

TranNet: TranNet Reveals Associations of Tumor Progression with the State  of Control Tissue  in Solid Tumors
Bayarbaatar Amgalan and Teresa M. Przytycka


Questions? bayarbaatar.amgalan@nih.gov

--------------------------------------
How to run the codes
--------------------------------------

* Construct a weighted-directed bipartite network representing transition map from normal to tumor expressions, 
  then compute scores of regulatory potentials of the genes

Input:  Two matrices denoting normal and tumor tissue samples of cancer patients. Their rows denote samples while the columns denote genes.
        "NormalData.csv": gene expression collected from normal tissues of the cancer patients 
        "TumorData.csv": gene expression collected from tumor tissues of the cancer patients 


Output: A weighted-directed bipartite network representing transition map from normal to tumor expressions
        "Transition_Matrix_M.csv": Non-zero entries m_ij represent the transition effect from regulator i-th regulator to j-th its target 
        A list representing genes with their contribution scores to the transition 
        "Regulatory_Potential.csv": Genes are sorted with respect to their scores for regulatory potential


--------------------------------------
Source Files
--------------------------------------

* TranNet.ipynb: a weighted-directed bipartite network representing transition map from normal to tumor expressions, and save the transition matrix 
  as "Transition_Matrix_M.csv", then compute the scores of regulatory potentials for the genes, and save the genes with their scores in "Regulatory_Potential.csv".
  The genes are sorted in the list, with respect to scores of their regulatory potential

