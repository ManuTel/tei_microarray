## Perform differential gene expression analysis for each control-treatment group in each generation (10 numbers) ####

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
                  return(efitlm)
                        
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
## Later we will focus on sex*interaction design at least in case of F0. since all the samples of F0 were processed in one batch
f0_male   <- fold_change(all, "Male", "F0")
f0_female <- exprs_mtrx(all, "Female", "F0")

f1_male     <- exprs_mtrx(all, "Male", "F1", "CD")
f1_maleHS   <- exprs_mtrx(all, "Male", "F1", "HSD")
f1_femaleLS <- exprs_mtrx(all, "Female", "F1", "CD")
f1_femaleHS <- exprs_mtrx(all, "Female", "F1", "HSD")

f2_maleLS   <- exprs_mtrx(all, "Male", "F2", "CD")
f2_maleHS   <- exprs_mtrx(all, "Male", "F2", "HSD")
f2_femaleLS <- exprs_mtrx(all, "Female", "F2", "CD")
f2_femaleHS <- exprs_mtrx(all, "Female", "F2", "HSD")
