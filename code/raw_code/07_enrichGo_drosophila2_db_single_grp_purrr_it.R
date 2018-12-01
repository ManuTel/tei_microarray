## Aim of this script is to make better functions for gene-enrichment analysis for each groups
## since we have to perform same functions over time and time for 10 groups, 20 lists (up, down)
## use of functions from purrr package seems relevant 

## Load necessary packages ####
source('code/raw_code/00_packages.R')
rm(list =ls())
## Load the foldchange data from csv file of interest ####
fcdata <- read.csv("./data/tidy_data/foldchange_df.csv")

## ClusterProfiler: Over-representation analysis using Affymetrix drosophila2.db annotations ####
library(clusterProfiler)
data(geneList) ## NOTE: the genelist is ordered foldchange vector with entrezid as names
gene <- names(geneList)[abs(geneList) > 2]
## 00: try with single data then use it as a template to do for all the data ####

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
## Done ##
## extract gene symbol and entrezid from drosophila2.db for respective probes
keytypes(drosophila2.db)

## purrr it ####
fc_nest <- fcdata %>% nest(-group)

## 3rd column as signfc: significant foldchange subset ####
sign_filter <- function(df) { df %>% 
                              filter(P.Value < 0.05 & 2^abs(logFC) >=1.3) %>%
                              arrange(desc(logFC))
                                }
fc_nest <- fc_nest %>%
              mutate(signfc = map(data, sign_filter))
glimpse(fc_nest$signfc)
## 4th column as genelist
fc_nest <- fc_nest %>%
              mutate(genelist = map(signfc, function(df) pull(df, probe_id)))
object.size(fc_nest)
## 5th column: perform gene enrichment analysis for updown genes
df_enrichgo <- function(x, df) {
                           enrichGO(gene    = x,
                              universe      = df$probe_id,
                              OrgDb         = drosophila2.db,
                              keyType       = "PROBEID", 
                              ont           = "BP",
                              pAdjustMethod = "BH",
                              pvalueCutoff  =  1,
                              qvalueCutoff  =  1,
                              readable      = TRUE)}
fc_nest <- fc_nest %>% mutate(gobp = map2(genelist, data, df_enrichgo))

## 6th column: up_signfc: significant foldchange up regulated only
fc_nest <- fc_nest %>%
              mutate(up_signfc = map(signfc, function(df) filter(df, logFC > 0)))
map(fc_nest$up_signfc, dim)
## 7th column as up genelist
fc_nest <- fc_nest %>%
              mutate(uplist = map(up_signfc, function(df) pull(df, probe_id)))


## 8th column: dn_signfc: significant foldchange down regulated only
fc_nest <- fc_nest %>%
              mutate(dn_signfc = map(signfc, function(df) filter(df, logFC < 0)))
map(fc_nest$dn_signfc, dim)
## 9th column as down genelist
fc_nest <- fc_nest %>%
  mutate(dnlist = map(dn_signfc, function(df) pull(df, probe_id)))

## 10th column: go enrichment for up regulated genes
fc_nest <- fc_nest %>% mutate(up_gobp = map2(uplist, data, df_enrichgo))

## 11th column: go enrichment for down regulated genes
fc_nest <- fc_nest %>% mutate(dn_gobp = map2(dnlist, data, df_enrichgo))
save(fc_nest,file = "./data/tidy_data/enrichgobp_all.RData")
load("./data/tidy_data/enrichgobp_all.RData")

## use unnest to extract 3 columns for each group ####
## gobp, up_gobp, dn_gopb
names(fc_nest)

updn_gobp <- fc_nest %>% 
              mutate(gobpdf = map(gobp, function(x) as.data.frame(x))) %>%
                select(group, gobpdf) %>%
                  unnest(gobpdf)
up_gobp   <- fc_nest %>% 
              mutate(up_gobpdf = map(up_gobp, function(x) as.data.frame(x))) %>%
                select(group, up_gobpdf) %>%
                  unnest(up_gobpdf)
dn_gobp   <- fc_nest %>% 
              mutate(dn_gobpdf = map(dn_gobp, function(x) as.data.frame(x))) %>%
                select(group, dn_gobpdf) %>%
                  unnest(dn_gobpdf)
