## Perform differential gene expression analysis for each control-treatment group in each generation (10 numbers) ####

## make function for getting differential gene expression ####
fold_change <- function(exps, sex = c("Female", "Male"), gen = c("F0","F1", "F2"), diet = c("CD", "HSD") ){
                  exps_subset   <- exps[, exps$sex %in% sex & exps$generation %in% gen & exps$diet %in% diet]
                  print(exps_subset$group) #check levels
                  print(sampleNames(exps_subset))
                  exps_rma      <- rma(exps_subset)
                  design        <- model.matrix(exps_subset$group)
                  print(design)
                  efitlm        <- exps_rma %>% lmFit(design) %>% eBayes()
                  print(topTable(efitlm, coef = 2, sort.by = "p", genelist = link))
                  return(efitlm)
                        
        }
#decide if we want efitlm object or dataframe of topTable
?treat()
