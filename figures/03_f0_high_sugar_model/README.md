# 03_f0_high_sugar_model: list of figures

### It should contain following figures
1. High sugar diet schematic
2. Results
    - body wt
    - hemolymph trehalose (normalized or original)
    - hemolymph glucose (normalized or original)
    - triglycerides level
    - F0 microarray related figures (heatmap, raw cel images, PCA)
    
### F0 microarray analysis design
1. basic design is finding the pairwise differences between high sugar diet (HSD) vs control diet(CD) in each sex, separately
    - HSD males vs CD males
    - HSD females vs CD females
2. nested design (Refer limma userguide section 9.5.3 A Nested Interaction Formula)
    - 
3. classic interaction design (Refer limma userguide section 9.5.4 Classic Interaction Models)
    - in case of F0 flies, there is no sex:diet interaction effect seen in triglycerides levels (we see increase in triglycerides irrespective of flies sex)
    - so it's better to use above two design methods to find differential gene expression in  HSD vs CD