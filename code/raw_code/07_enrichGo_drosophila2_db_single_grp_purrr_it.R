## Aim of this script is to make better functions for gene-enrichment analysis for each groups
## since we have to perform same functions over time and time for 10 groups, 20 lists (up, down)
## use of functions from purrr package seems relevant 

## Load necessary packages ####
source('code/raw_code/00_packages.R')

## Load the foldchange data from csv file of interest ####
fcdata <- read.csv("./data/tidy_data/foldchange_df.csv")

## ClusterProfiler: Over-representation analysis using Affymetrix drosophila2.db annotations ####
library(clusterProfiler)
data(geneList) ## NOTE: the genelist is ordered foldchange vector with entrezid as names
gene <- names(geneList)[abs(geneList) > 2]
## 00: try with single data then use it as a template to do for all the data

## 01: prepare the data for HS vs LS male group
fcmale <- fcdata %>% filter(group == "f0_male") %>% droplevels()
head(fcmale); summary(fcmale)

## 02: filter the significant gene list and arrange by desc(logFC)
sign_fc <- fcmale %>% 
              filter(P.Value <  0.05 & 2^abs(logFC) >= 1.3) %>%
                arrange(desc(logFC))
## 03: make genelist for enrichgo function
genelist <- sign_fc$probe_id

## 04: perform enrichment analysis for up-down genes
## genelist contains all the genes which are upregulated as well as downregulated

egobp <- enrichGO(gene          = genelist,
                  universe      = fcmale$probe_id,
                  OrgDb         = drosophila2.db,
                  keyType       = "PROBEID", 
                  ont           = "BP",
                  pAdjustMethod = "BH",
                  pvalueCutoff  =  1,
                  qvalueCutoff  =  1,
                  readable      = TRUE)
##05a: select probe_id for upregulated genes
upgene   <- sign_fc %>% filter(logFC > 0) %>% pull(probe_id)

##05b: perform enrichment for up regulated genes
egobp_up <- enrichGO(gene         = upgene,
                     universe      = fcmale$probe_id,
                     OrgDb         = drosophila2.db,
                     keyType       = "PROBEID", 
                     ont           = "BP",
                     pAdjustMethod = "BH",
                     pvalueCutoff  =  1,
                     qvalueCutoff  =  1,
                     readable      = TRUE)

##06a: select probe_id for downregulated genes
dngene   <- sign_fc %>% filter(logFC < 0) %>% pull(probe_id)

##06b: perform enrichment for down regulated genes
egobp_dn <- enrichGO(gene          = dngene,
                     universe      = fcmale$probe_id,
                     OrgDb         = drosophila2.db,
                     keyType       = "PROBEID", 
                     ont           = "BP",
                     pAdjustMethod = "BH",
                     pvalueCutoff  =  1,
                     qvalueCutoff  =  1,
                     readable      = TRUE)

## extract gene symbol and entrezid from drosophila2.db for respective probes
keytypes(drosophila2.db)
