## Load necessary packages ####
source('code/raw_code/00_packages.R')

## Load the foldchange data from csv file of interest ####
fcdata <- read.csv("./data/tidy_data/f0_diet_sex_interaction_fc.csv")

## ClusterProfiler: Over-representation analysis using Affymetrix drosophila2.db annotations ####
library(clusterProfiler)
data(geneList) ## NOTE: the genelist is ordered foldchange vector with entrezid as names

## prepare the data for HS vs LS male group
fcmale <- fcdata %>% filter(group == "male") %>% droplevels()
head(fcmale); summary(fcmale)

## extract gene symbol and entrezid from drosophila2.db for respective probes
keytypes(drosophila2.db)
fclist <- bitr(fcmale$probe_id, fromType="PROBEID", toType=c("ENTREZID", "SYMBOL"), OrgDb="drosophila2.db")
head(fclist)
fclist <- fclist %>% 
              left_join(select(fcmale, probe_id, logFC, P.Value), 
                        by = c("PROBEID" = "probe_id")) %>% 
              arrange(desc(logFC)) ## IMPORTANT

# NOTE: One probe_id to Many Entrezid mapping occurred during bitr command
fclist %>% group_by(SYMBOL) %>% summarise(n = n()) %>% filter(n > 1) %>% arrange(desc(n))
fclist %>% filter(SYMBOL == "mod(mdg4)")

gene <- fclist %>% 
            filter(P.Value <  0.05 & 2^abs(logFC) >= 1.3) %>% 
            select(PROBEID)
## Note the gene vector is already arranged by decreasing foldchange #refer to last fclist operation 
## gene contains all the genes which are upregulated as well as downregulated

egobp <- enrichGO(gene          = gene$PROBEID,
                  universe      = fclist$PROBEID,
                  OrgDb         = drosophila2.db,
                  keyType       = "PROBEID", 
                  ont           = "BP",
                  pAdjustMethod = "BH",
                  pvalueCutoff  =  1,
                  qvalueCutoff  =  1,
                  readable      = TRUE)

# save start ####
library(XLConnect)
wb <- loadWorkbook("./data/tidy_data/F0_interaction_GO.xlsx", create = TRUE)
createSheet(wb, name = "details")

#--- SAVE-1 ---
createSheet(wb, name = "BP_enrich_male_updown")
writeWorksheet(wb, as.data.frame(egobp), sheet = "BP_enrich_male_updown", startRow = 2)

#ego2 <- dropGO(egobp, level = 6)
#enrichMap(egobp)
#dotplot(egobp,  showCategory = 10)
#barplot(egobp, showCategory = 10)

## UPREGULATED GENES ONLY ####
up <- fclist %>% 
        filter(P.Value <  0.05 & 2^abs(logFC) >= 1.3) %>% 
        filter(logFC > 0) %>%
        select(PROBEID)


upegobp <- enrichGO(gene          = up$PROBEID,
                    universe      = fclist$PROBEID,
                    OrgDb         = drosophila2.db,
                    keyType       = "PROBEID", 
                    ont           = "BP",
                    pAdjustMethod = "BH",
                    pvalueCutoff  =  1,
                    qvalueCutoff  =  1,
                    readable      = TRUE)
dim(upegobp)
enrichMap(upegobp)

#--- SAVE-2 ---
createSheet(wb, name = "BP_enrich_male_up")
writeWorksheet(wb, as.data.frame(upegobp), sheet = "BP_enrich_male_up", startRow = 2)
#ego3 <- dropGO(upegobp, level =8)
#enrichMap(ego3)

## DOWNREGULATED GENES ONLY ####
down <- fclist %>% 
  filter(P.Value <  0.05 & 2^abs(logFC) >= 1.3) %>% 
  filter(logFC < 0) %>%
  select(PROBEID)


dnegobp <- enrichGO(gene          = down$PROBEID,
                    universe      = fclist$PROBEID,
                    OrgDb         = drosophila2.db,
                    keyType       = "PROBEID", 
                    ont           = "BP",
                    pAdjustMethod = "BH",
                    pvalueCutoff  =  1,
                    qvalueCutoff  =  1,
                    readable      = TRUE)
dim(dnegobp)
enrichMap(dnegobp)

#--- SAVE-3 ---
createSheet(wb, name = "BP_enrich_male_dn")
writeWorksheet(wb, as.data.frame(dnegobp), sheet = "BP_enrich_male_dn", startRow = 2)
#ego4 <- dropGO(dnegobp, level =8)
#enrichMap(ego4)
saveWorkbook(wb)
#wb <- loadWorkbook("./data/tidy_data/F0_interaction_GO.xlsx")

## prepare the data for HS vs LS female group
fcfemale <- fcdata %>% filter(group == "female") %>% droplevels()
head(fcfemale); summary(fcfemale)

## extract gene symbol and entrezid from drosophila2.db for respective probes
keytypes(drosophila2.db)
fclist <- bitr(fcfemale$probe_id, fromType="PROBEID", toType=c("ENTREZID", "SYMBOL"), OrgDb="drosophila2.db")
head(fclist)
fclist <- fclist %>% 
  left_join(select(fcfemale, probe_id, logFC, P.Value), 
            by = c("PROBEID" = "probe_id")) %>% 
  arrange(desc(logFC)) ## IMPORTANT

# NOTE: One probe_id to Many Entrezid mapping occurred during bitr command
fclist %>% group_by(SYMBOL) %>% summarise(n = n()) %>% filter(n > 1) %>% arrange(desc(n))
fclist %>% filter(SYMBOL == "mod(mdg4)")

gene <- fclist %>% 
  filter(P.Value <  0.05 & 2^abs(logFC) >= 1.3) %>% 
  select(PROBEID)
## Note the gene vector is already arranged by decreasing foldchange #refer to last fclist operation 
## gene contains all the genes which are upregulated as well as downregulated

egobp <- enrichGO(gene          = gene$PROBEID,
                  universe      = fclist$PROBEID,
                  OrgDb         = drosophila2.db,
                  keyType       = "PROBEID", 
                  ont           = "BP",
                  pAdjustMethod = "BH",
                  pvalueCutoff  =  1,
                  qvalueCutoff  =  1,
                  readable      = TRUE)

dim(egobp)
#--- SAVE-4 ---
createSheet(wb, name = "BP_enrich_female_updown")
writeWorksheet(wb, as.data.frame(egobp), sheet = "BP_enrich_female_updown", startRow = 2)

#ego2 <- dropGO(egobp, level = 6)
#enrichMap(egobp)
#dotplot(egobp,  showCategory = 10)
#barplot(egobp, showCategory = 10)

## UPREGULATED GENES ONLY ####
up <- fclist %>% 
  filter(P.Value <  0.05 & 2^abs(logFC) >= 1.3) %>% 
  filter(logFC > 0) %>%
  select(PROBEID)


upegobp <- enrichGO(gene          = up$PROBEID,
                    universe      = fclist$PROBEID,
                    OrgDb         = drosophila2.db,
                    keyType       = "PROBEID", 
                    ont           = "BP",
                    pAdjustMethod = "BH",
                    pvalueCutoff  =  1,
                    qvalueCutoff  =  1,
                    readable      = TRUE)
dim(upegobp)
enrichMap(upegobp)

#--- SAVE-5 ---
createSheet(wb, name = "BP_enrich_female_up")
writeWorksheet(wb, as.data.frame(upegobp), sheet = "BP_enrich_female_up", startRow = 2)
#ego3 <- dropGO(upegobp, level =8)
#enrichMap(ego3)

## DOWNREGULATED GENES ONLY ####
down <- fclist %>% 
  filter(P.Value <  0.05 & 2^abs(logFC) >= 1.3) %>% 
  filter(logFC < 0) %>%
  select(PROBEID)


dnegobp <- enrichGO(gene          = down$PROBEID,
                    universe      = fclist$PROBEID,
                    OrgDb         = drosophila2.db,
                    keyType       = "PROBEID", 
                    ont           = "BP",
                    pAdjustMethod = "BH",
                    pvalueCutoff  =  1,
                    qvalueCutoff  =  1,
                    readable      = TRUE)
dim(dnegobp)
enrichMap(dnegobp)

#--- SAVE-6 ---
createSheet(wb, name = "BP_enrich_female_dn")
writeWorksheet(wb, as.data.frame(dnegobp), sheet = "BP_enrich_female_dn", startRow = 2)
#ego4 <- dropGO(dnegobp, level =8)
#enrichMap(ego4)
saveWorkbook(wb)

#############################
## prepare the data for HS:male interaction group group

fchsmale <- fcdata %>% filter(group == "HSmale") %>% droplevels()
head(fchsmale); summary(fchsmale)

## extract gene symbol and entrezid from drosophila2.db for respective probes
keytypes(drosophila2.db)
fclist <- bitr(fchsmale$probe_id, fromType="PROBEID", toType=c("ENTREZID", "SYMBOL"), OrgDb="drosophila2.db")
head(fclist)
fclist <- fclist %>% 
  left_join(select(fchsmale, probe_id, logFC, P.Value), 
            by = c("PROBEID" = "probe_id")) %>% 
  arrange(desc(logFC)) ## IMPORTANT

# NOTE: One probe_id to Many Entrezid mapping occurred during bitr command
fclist %>% group_by(SYMBOL) %>% summarise(n = n()) %>% filter(n > 1) %>% arrange(desc(n))
fclist %>% filter(SYMBOL == "mod(mdg4)")

gene <- fclist %>% 
  filter(P.Value <  0.05 & 2^abs(logFC) >= 1.3) %>% 
  select(PROBEID)
## Note the gene vector is already arranged by decreasing foldchange #refer to last fclist operation 
## gene contains all the genes which are upregulated as well as downregulated

egobp <- enrichGO(gene          = gene$PROBEID,
                  universe      = fclist$PROBEID,
                  OrgDb         = drosophila2.db,
                  keyType       = "PROBEID", 
                  ont           = "BP",
                  pAdjustMethod = "BH",
                  pvalueCutoff  =  1,
                  qvalueCutoff  =  1,
                  readable      = TRUE)

dim(egobp)
#--- SAVE-7 ---
createSheet(wb, name = "BP_enrich_HSxmale_updown")
writeWorksheet(wb, as.data.frame(egobp), sheet = "BP_enrich_HSxmale_updown", startRow = 2)

#ego2 <- dropGO(egobp, level = 6)
#enrichMap(egobp)
#dotplot(egobp,  showCategory = 10)
#barplot(egobp, showCategory = 10)

## UPREGULATED GENES ONLY ####
up <- fclist %>% 
  filter(P.Value <  0.05 & 2^abs(logFC) >= 1.3) %>% 
  filter(logFC > 0) %>%
  select(PROBEID)


upegobp <- enrichGO(gene          = up$PROBEID,
                    universe      = fclist$PROBEID,
                    OrgDb         = drosophila2.db,
                    keyType       = "PROBEID", 
                    ont           = "BP",
                    pAdjustMethod = "BH",
                    pvalueCutoff  =  1,
                    qvalueCutoff  =  1,
                    readable      = TRUE)
dim(upegobp)
enrichMap(upegobp)

#--- SAVE-8 ---
createSheet(wb, name = "BP_enrich_HSxmale_up")
writeWorksheet(wb, as.data.frame(upegobp), sheet = "BP_enrich_HSxmale_up", startRow = 2)
#ego3 <- dropGO(upegobp, level =8)
#enrichMap(ego3)

## DOWNREGULATED GENES ONLY ####
down <- fclist %>% 
  filter(P.Value <  0.05 & 2^abs(logFC) >= 1.3) %>% 
  filter(logFC < 0) %>%
  select(PROBEID)


dnegobp <- enrichGO(gene          = down$PROBEID,
                    universe      = fclist$PROBEID,
                    OrgDb         = drosophila2.db,
                    keyType       = "PROBEID", 
                    ont           = "BP",
                    pAdjustMethod = "BH",
                    pvalueCutoff  =  1,
                    qvalueCutoff  =  1,
                    readable      = TRUE)
dim(dnegobp)
enrichMap(dnegobp)

#--- SAVE-9 ---
createSheet(wb, name = "BP_enrich_HSxmale_dn")
writeWorksheet(wb, as.data.frame(dnegobp), sheet = "BP_enrich_HSxmale_dn", startRow = 2)
#ego4 <- dropGO(dnegobp, level =8)
#enrichMap(ego4)
saveWorkbook(wb)

