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
## all three coef from fit2 model are of important for us to study 
## since they represent HSvsLs in female, male and interaction of diet with sex

## Effect of HS diet in Female
f0_hs_female   <- topTable(fit2, coef = 1, sort.by = "none", genelist = link, number = Inf)
f0_hs_female$group <- "female"

## Effect of HS diet in Male
f0_hs_male     <- topTable(fit2, coef = 2, sort.by = "none", genelist = link, number = Inf)
f0_hs_male$group <- "male"

## Interaction of HS diet and male sex
f0_hs_male_int <- topTable(fit2, coef = 3, sort.by = "none", genelist = link, number = Inf) 
f0_hs_male_int$group <- "HSmale"
all_df   <- lapply(ls(pattern = "^f0_hs*"), get)
library(tibble)
all_df   <- map2(all_df, "probe_id", rownames_to_column)
all_fc   <- bind_rows(all_df)

#write.csv(all_fc, file = "./data/tidy_data/f0_diet_sex_interaction_fc.csv", row.names = FALSE)

## Exploratory data analysis ####
summary(all_fc)
updown <- all_fc %>% 
            group_by(group) %>%
              filter(P.Value <  0.05 &  ## NOTE: adjusted p value used
                      2^abs(logFC) >= 1.3) %>% ## NOTE: fold change 1.3
              summarise(n = n(), up = sum(logFC > 0),
                               down = sum(logFC < 0))
## stringent criteria of 2-fold change and adjusted p-value cut-off 0.05

sign_fc <- all_fc %>%
            group_by(group) %>%
            filter(adj.P.Val < 0.05 &
                     2^abs(logFC) > 2)
        

Reduce(intersect, split(abc$ID, abc$group))
## compare with normal HS vs Ls comparison per sex
fc   <- read.csv("./data/tidy_data/foldchange_df.csv")
f0_in <- read.csv("./data/tidy_data/f0_diet_sex_interaction_fc.csv")
# the fold change is similar across both the model (interaction and normal HS vs LS)
# but interaction model has more significant results as compare to normal model
## do check for F1 and F2 generation
