## This file contains detailed description of tidy data containt

1. **expressionset_matrix.csv**
    - Contains actual log normalized probeset intensity for all the samples (80 chips) in wide format
    - it will be handy if one wants to show actual gene expression instead of fold change
2. **foldchange_df.csv**
    - this csv holds differential gene expression analysis for all 10 comparison as stated in Readme of the project
    - the ID column represent drosophila gene symbol, other columns are self explainatory
3. **diff_fc.RData**
    - it contains output generated from the script *03_descriptive_plots.R*
    - **updown**
        - it contains contingency table for number of genes which are significantly up-regulated, down-regulated for each comparison 
    - **sign_fc**
        - sign_fc object holds list of genes which shows significant differential gene expression. 
        - this list has equal number of genes as represented in *updown* table
    - **fc_2793**
        - this holds list of 2793 genes for each comparison (27930).
        - there are 2793 genes which shows differential gene expression in at least one group comparison
    - **fc_1620**
        - it has list of 2793 genes for each comparsion 
        - this list is obtained after excluding *f0_female* group, since heatmap of *fc_2793* shows that all the expression signature is overshadowed by *f0_female* group (this single group has 1845 diff expressed genes)
        