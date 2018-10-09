## install necessary packages
source('code/raw_code/00_packages.R')

## load necessary R packages ####
library(affy)
library(limma)
library(drosophila2.db)
library(drosophila2cdf)
library(dplyr)
library(purrr)
library(broom)
library(ggplot2)
library(ggsci)
library(ggpubr)
library(gridExtra)

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

## make function to extract expression set for each comparison ####
exprs_mtrx <- function(exps, sex = c("Female", "Male"), gen = c("F0","F1", "F2"), diet = c("CD", "HSD") ){
                        exps_subset   <- exps[, exps$sex %in% sex & exps$generation %in% gen & exps$diet %in% diet]
                        print(exps_subset$group)
                        print(sampleNames(exps_subset))
                        mtrx    <- rma(exps_subset) %>% exprs()
                        colnames(mtrx) <- exps_subset$name
                        print(head(mtrx))
                        return(mtrx)
                        }

f0_male   <- exprs_mtrx(all, "Male", "F0")
f0_female <- exprs_mtrx(all, "Female", "F0")

f1_male     <- exprs_mtrx(all, "Male", "F1", "CD")
f1_maleHS   <- exprs_mtrx(all, "Male", "F1", "HSD")
f1_femaleLS <- exprs_mtrx(all, "Female", "F1", "CD")
f1_femaleHS <- exprs_mtrx(all, "Female", "F1", "HSD")

f2_maleLS   <- exprs_mtrx(all, "Male", "F2", "CD")
f2_maleHS   <- exprs_mtrx(all, "Male", "F2", "HSD")
f2_femaleLS <- exprs_mtrx(all, "Female", "F2", "CD")
f2_femaleHS <- exprs_mtrx(all, "Female", "F2", "HSD")

## combine all the expression set ####
all_df   <- lapply(ls(pattern = "^f.*"), get)
all_mtrx <- do.call("cbind.data.frame", all_df)
write.csv(all_mtrx, file = "./data/tidy_data/expressionset_matrix.csv")
