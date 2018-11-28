## F0 Diet + Sex + DietxSex interaction ####

## this script is related to all the microarray analysis done on F0 samples to check 
## independent effect of diet, sex and their interaction #

## first use model 
design <- model.matrix(~sex*diet, data = pData(f0)) # baseline sex doesn't affect the diet effect on DE genes
## and then
design2 <- model.matrix(~sex + diet, data = pData(f0)) #check baseline sex whether affect the diet-effect on DEgenes
## extract diet coefficient from both analysis and check for 
source('code/raw_code/00_packages.R')
library(tibble)

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
## take out the only F0 data and make male baseline
f0      <- all[, all$generation == "F0"]
f0$sex  <- factor(f0$sex, levels = c("Male", "Female")) 
f0$diet <- factor(f0$diet, levels = c("CD", "HSD"))

f0_rma <- rma(f0)

design1 <- model.matrix(~ sex + diet, data = pData(f0))
design1

fit <- lmFit(f0_rma, design1)
fit1 <- eBayes(fit)
print(topTable(fit1, coef = 3, sort.by = "p", genelist = link, number = 10))
f0_hs     <- topTable(fit1, coef = 3, sort.by = "none", genelist = link, number = Inf)
summary(f0_hs)
library(tibble)
f0_hs <- rownames_to_column(f0_hs, "probe_id")
