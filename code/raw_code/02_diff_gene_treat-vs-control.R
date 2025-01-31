## Perform differential gene expression analysis for each control-treatment group in each generation (10 numbers) ####
## load necessary R packages ####
library(affy)
library(limma)
library(drosophila2.db)
library(drosophila2cdf)
library(dplyr)
library(tidyr)
library(purrr)
library(broom)
library(tibble)
library(ggplot2)
library(ggsci)
library(ggpubr)
library(gridExtra)
## make function for getting differential gene expression ####
fold_change <- function(exps, sex = c("Female", "Male"), gen = c("F0","F1", "F2"), diet = c("CD", "HSD") ){
                  exps_subset   <- exps[, exps$sex %in% sex & exps$generation %in% gen & exps$diet %in% diet]
                  print(exps_subset$group) #check levels
                  print(sampleNames(exps_subset))
                  exps_rma      <- rma(exps_subset)
                  design        <- model.matrix(~group, data = pData(exps_subset))
                  print(design)
                  efitlm        <- exps_rma %>% lmFit(design) %>% eBayes()
                  print(topTable(efitlm, coef = 2, sort.by = "p", genelist = link))
                  return(topTable(efitlm, coef = 2, sort.by = "none", genelist = link, number = Inf))
                  }

## Read .cel files with there phenodata ####
cel_files <- list.files(path="data/raw_data/", pattern = '.*CEL')
phen      <- read.AnnotatedDataFrame(filename = "data/raw_data/phenodata.txt")
all       <- ReadAffy(celfile.path = "data/raw_data/", phenoData = phen)
all$sex          ## double check the labels
all$name 
table(all$sex, all$diet, all$generation)
all$group <- factor(all$group, levels = c('control', 'case'))
all$group

link      <- unlist(mget(featureNames(all), envir =  drosophila2SYMBOL))

## calculate fold change for each group with simple design of CASE vs CONTROL ####

f0_male     <- fold_change(all, "Male", "F0");          f0_male["group"]     <- "f0_male"
f0_female   <- fold_change(all, "Female", "F0");        f0_female["group"]   <- "f0_female"

f1_maleLS   <- fold_change(all, "Male", "F1", "CD");    f1_maleLS["group"]   <- "f1_maleLS"
f1_maleHS   <- fold_change(all, "Male", "F1", "HSD");   f1_maleHS["group"]   <- "f1_maleHS"
f1_femaleLS <- fold_change(all, "Female", "F1", "CD");  f1_femaleLS["group"] <- "f1_femaleLS"
f1_femaleHS <- fold_change(all, "Female", "F1", "HSD"); f1_femaleHS["group"] <- "f1_femaleHS"

f2_maleLS   <- fold_change(all, "Male", "F2", "CD");    f2_maleLS["group"]   <- "f2_maleLS"
f2_maleHS   <- fold_change(all, "Male", "F2", "HSD");   f2_maleHS["group"]   <- "f2_maleHS"
f2_femaleLS <- fold_change(all, "Female", "F2", "CD");  f2_femaleLS["group"] <- "f2_femaleLS"
f2_femaleHS <- fold_change(all, "Female", "F2", "HSD"); f2_femaleHS["group"] <- "f2_femaleHS"

## combine all the foldchange ####
all_df   <- lapply(ls(pattern = "^f[0-2].*"), get)
all_df   <- map2(all_df, "probe_id", rownames_to_column)
all_fc   <- bind_rows(all_df)
write.csv(all_fc, file = "./data/tidy_data/foldchange_df.csv", row.names = FALSE)

