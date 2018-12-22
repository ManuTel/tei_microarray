# basic images will go here

1. gsea_ls_only image is the heatmap of all samples reared on LS diet (except F0)
    - The command for plot is available in *07_enrichGo_drosophila2_db_single_grp_purrr_it.R*
    - The gsea results (minimum geneset size = 50) were filtered for following criteria
        - p value < 0.05
        - only LS samples (6 samples: 2 F0, 2F1, 2F2)
        - at least has significant Normalized Enrichment Score (NES) in 5 samples out of six
        - red is downregulated, yellow is upregulated (you can look for the actual scale in the filtered dataframe)
2. gsea_hs_only is same as above but for HS diet sample
    - Nothing interesting here, svg and tiff are saved in 00_rough folder

3. gsea_ls_only_3 is re-paletted version of gsea_ls_only
    - the gsea_ls_only truncated from downside where there is NA for HS male group for better viewing
  