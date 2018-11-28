## F0 Diet + Sex + DietxSex interaction ####

## this script is related to all the microarray analysis done on F0 samples to check 
## independent effect of diet, sex and their interaction #

## first use model 
design <- model.matrix(~sex*diet, data = pData(f0)) # baseline sex affect the diet effect on DE genes 
# so its decided that male will be baseline since we are more interested in effect of diet on male 
# (not the interaction)
## and then
design2 <- model.matrix(~sex + diet, data = pData(f0)) # baseline sex doesn't affect the diet-effect on DEgenes
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
## Model 0: subset of male on hs diet####
male <- f0_rma[, f0_rma$sex == "Male"]

design0 <- model.matrix(~ diet, data = pData(male))
design0

fit0 <- lmFit(male, design0)
fit0 <- eBayes(fit0)
print(topTable(fit0, coef = 2, sort.by = "p", genelist = link, number = 10))

## Model 1: diet*sex interaction design to study effect of diet on male ####
design1 <- model.matrix(~ diet*sex, data = pData(f0_rma))
design1

fit1 <- lmFit(f0_rma, design1)
fit1 <- eBayes(fit1)
print(topTable(fit1, coef = 2, sort.by = "p", genelist = link, number = 10))

## Model 2: diet +  sex design to study effect of diet regressed out for sex
design2 <- model.matrix(~ diet + sex, data = pData(f0_rma))
design2

fit2 <- lmFit(f0_rma, design2)
fit2 <- eBayes(fit2)
print(topTable(fit2, coef = 2, sort.by = "p", genelist = link, number = 10))

## Conclusion: Model 0 (baseline) and Model 1(diet*sex interaction) yields 
## the same fold change but different significance
## Model 2 differs in the fold change too, since we regress out the sex
mod_base    <- topTable(fit0, coef = 2, sort.by = "none", genelist = link, number = Inf)
mod_base$group <- "base"
mod_int     <- topTable(fit1, coef = 2, sort.by = "none", genelist = link, number = Inf)
mod_int$group <- "int"
mod_plus    <- topTable(fit2, coef = 2, sort.by = "none", genelist = link, number = Inf)
mod_plus$group <- "plus"

## combine all the foldchange ####
all_df   <- lapply(ls(pattern = "mod_*"), get)
all_df   <- map2(all_df, "probe_id", rownames_to_column)
all_mod  <- bind_rows(all_df)

updown <- all_mod %>% 
  group_by(group) %>%
  filter(P.Value <  0.05 &
           2^abs(logFC) >= 1.3)
updown
mat_updn <- updown %>% 
              select(probe_id, ID, logFC, group) %>% 
                spread(group, logFC)
