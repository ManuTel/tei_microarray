### The goal of this analysis is to formulate simple pipeline to get differentially regulated genes for
# given subset of samples (for eg. HS male vs LS male,  LS f1male from HS parent vs LS f1male from LS parent)
### the initial problem arised during normalizing all the samples with RMA
# so now it is required that first read all the cel files using ReadAffy and then
# do RMA only for the selected two groups which we want to compare 

library(affy)
library(limma)
library(dplyr)
library(drosophila2.db)
setwd('~/Desktop/trans all data microarray/all_cel/')

cel_files          <- list.files(path="./", pattern = '.*CEL')
phen               <- read.AnnotatedDataFrame(filename = "phenodata.txt")

all <- ReadAffy(phenoData = phen)
phenoData(all)$sex          ## double check the labels
phenoData(all)$diet 
phenoData(all)$generation
phenoData(all)$group <- factor(phenoData(all)$group, levels = c('control', 'case'))
phenoData(all)$group

link            <- unlist(mget(featureNames(all), 
                               envir =  drosophila2SYMBOL))
#### select F0 --- Male --- HS vs LS --- Do RMA --- Diff genes ####
one <- all[ ,phenoData(all)$sex == 'Male' & phenoData(all)$generation == "F0"]
phenoData(one)$diet
phenoData(one)$group

abc <- rma(one)
design <- model.matrix(~phenoData(abc)$group)
f0_rma <- abc %>% lmFit(design) %>% eBayes()
topTable(f0_rma, coef = 2, sort.by = "p", genelist = link)

#### now compare with direct justRMA of f0 generation ####
# cel_files2 <- list.celfiles(path = "../F0_male/") 
# setwd("../F0_male/")
# phen2      <- read.AnnotatedDataFrame(filename = "../F0_male/onlyf0.txt")
# two        <- justRMA(phenoData = phen2)
# two$group  <- factor(two$group, levels = c('control', 'case'))
# design     <- model.matrix(~phenoData(two)$group)
# link       <- unlist(mget(featureNames(two), 
#                         envir =  drosophila2SYMBOL))
# f0_justrma <- two %>% lmFit(design) %>% eBayes()
# setequal(f0_rma, f0_justrma) ## it was TRUE
#### conclusion : the process Readaffy() --> rma() --> lmFit() --> eBayes() is same as  ####
####                                    justRMA()  --> lmFit() --> eBayes()             #

make_rma_function <- function(exps, sex, gen){
 exps_subset   <- exps[, phenoData(exps)$sex == sex & phenoData(exps)$generation == gen]
 print(phenoData(exps_subset)$group)
 print(sampleNames(exps_subset))
 exps_rma      <- rma(exps_subset)
 design        <- model.matrix(~phenoData(exps_rma)$group)
 efitlm        <- exps_rma %>% lmFit(design) %>% eBayes()
 print(topTable(efitlm, coef = 2, sort.by = "p", genelist = link))
 return(efitlm)
}
a <- make_rma_function(all, 'Male', 'F0')
# setequal(exprs(a), exprs(abc))
topTable(a, coef = 2, sort.by = "p", genelist = link)
