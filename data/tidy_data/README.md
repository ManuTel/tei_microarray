### This file contains detailed description of tidy data folder

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
4. **f0_diet_sex_interaction_fc.csv**
    - it's output of the script *04_f0_sex_diet_interaction.R*
    - it has differential gene expression of F0 flies calculated with *model.matrix(~sex\*diet)*
    - three groups in group column denotes
        - male (HS vs LS diet)
        - female (HS vs LS diet)
        - HSmale (interaction of HS diet and Male sex)
    - conclusion: Only male data has relevant gene expression which reflects effect of high sugar
    
5. **F0_interaction_GO.xlsx**
    - This file is the output of the scripts *05_GO_analysis_part_1.R* and *05_GO_analysis_part_2.R*
    - It processes the data *f0_diet_sex_interaction_fc.csv* for Gene-ontology over-representation analysis for biological pathway only
    - for each group (male, female, HSmale), we formed three subsets of genelist (mixed, up, down) by applying filter (fold-change >= 1.3 and p-value < 0.05) and looked for the GO over-representation analysis using ```enrichGO``` from package ```clusterProfiler``` and annotation used from the package ```drosophila2.db``` - Affymetrix annnotation package at probe level information
    - for time-being p-value cutoff is not defined (tentative is unadjusted p-val < 0.05)
    - conclusion: Only male data has relevant gene expression which reflects effect of high sugar
6. **sign_f0_DietSex_intr_web_gsea.csv**
    - intermittant file (delelte later)