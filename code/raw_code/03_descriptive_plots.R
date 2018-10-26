## several plots for microarray analysis ####

## load necessary R packages ####
library(affy)
library(limma)
library(drosophila2.db)
library(drosophila2cdf)
library(dplyr)
library(tidyr)
library(purrr)
library(broom)
library(ggplot2)
library(ggsci)
library(ggpubr)
library(gridExtra)
library(stringr)

## load expressionset matrix and foldchange dataframe ####
## expressionset matrix contains the expression value for each gene while
## foldchange dataframe contains 10 groups (case vs control) differential gene expression profiles
eset <- read.csv("./data/tidy_data/expressionset_matrix.csv")
fc   <- read.csv("./data/tidy_data/foldchange_df.csv")

updown <- fc %>% 
              group_by(group) %>%
              filter(P.Value <  0.05 &
                       2^abs(logFC) >= 1.3) %>%
              summarise(n = n(), up = sum(logFC > 0),
                               down = sum(logFC < 0))

## find the total number of genes that are significantly (p<0.05) differentially expressed (FC>=1.3) in 
## at least one comparison 
fc %>% 
  group_by(group) %>%
  filter(P.Value <  0.05  &
           2^abs(logFC) >= 1.3) %>% # & !str_detect(group, "HS") if we want only LS diet comparison #group != "f0_female"
  ungroup() %>%
  summarise(n = n_distinct(probe_id)) 
# 2793 probesets are differentially expressed in at least one comparison ##
# 2683 probesets if we consider only LS diet in F1 F2
# 2306 probesets if we consider only HS diet in F1 F2
sign_fc <- fc %>%
              group_by(group) %>%
              filter(P.Value <  0.05  &
                       2^abs(logFC) >= 1.3) %>%
              ungroup()
length(unique(sign_fc$probe_id))
length(unique(sign_fc$ID))
fc_2793 <- fc %>%
              filter(probe_id %in% unique(sign_fc$probe_id)) %>%
              select(probe_id, ID, logFC, group) %>%
              spread(group, logFC)

## Save outputs in RData
save(sign_fc, fc_2793, fc_1620, updown, file = "./data/tidy_data/diff_fc.RData")
heatmap(as.matrix(fc_1620[3:11]), Colv = NA)
  
