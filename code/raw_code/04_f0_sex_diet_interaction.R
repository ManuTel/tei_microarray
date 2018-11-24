## F0 Diet x Sex interaction ####

## this script is related to all the microarray analysis done on F0 samples to check for diet and sex interaction #
source('code/raw_code/00_packages.R')
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
## take out the only F0 data and make FEMALE baseline
f0      <- all[, all$generation == "F0"]
f0$sex  <- factor(f0$sex, levels = c("Female", "Male")) 
f0$diet <- factor(f0$diet, levels = c("CD", "HSD"))

f0_rma <- rma(f0)
design <- model.matrix(~sex*diet, data = pData(f0))
design
## Reference: limma userguide
## 9.5.4 Classic Interaction Models: page number 45 

# Coefficient          Comparison                                    Interpretation
# Intercept           Female.LS                                     Baseline level of female on LS
# sexMale             Male.LS - Female.LS                           Difference between sexes ##not of interest
# dietHSD             Female.HS - Female.LS                         Diet effect for female
# sexMale:dietHSD    (Male.HS-Male.LS)-(Female.HS-Famel.LS)         Interaction


fit <- lmFit(f0_rma, design)
fit_1bayes <- eBayes(fit)
print(topTable(fit_1bayes, coef = 2, sort.by = "p", genelist = link, number = 10))

## Notice that one of our comparisons of interest, Male.HS-Male.LS, is not represented 
## and instead the comparison Male.LS-Female.LS, which might not be of direct interest, is included
cont.matrix <- cbind(HSvsLSinFemale = c(0,0,1,0), ## same as 3rd coefficient in fit_1bayes model 
                     HSvsLSinMale   = c(0,0,1,1), ## this was absent in fit_1bayes model
                     Diff           = c(0,0,0,1)) ## same as 4th interaction coef in fit_1bayes model 
cont.matrix
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)
print(topTable(fit2, coef = 3, sort.by = "p", genelist = link, number = 10))
## all three coef from fit2 model are of important for us to study since they represent HSvsLs in female, male and interaction of diet with sex
