## several plots for microarray analysis ####

## load necessary R packages ####
library(affy)
library(limma)
library(drosophila2.db)
library(drosophila2cdf)
library(dplyr)
library(tidyr)
library(purrr)
library(broom)
library(ggplot2)
library(ggsci)
library(ggpubr)
library(gridExtra)

## load expressionset matrix and foldchange dataframe ####
## expressionset matrix contains the expression value for each gene while
## foldchange dataframe contains 10 groups (case vs control) differential gene expression profiles
eset <- read.csv("./data/tidy_data/expressionset_matrix.csv")
fc   <- read.csv("./data/tidy_data/foldchange_df.csv")

## F0 - PCA, Heatmap, Volcano
eset_f0 <- eset %>%
            select(X, contains("_f0_")) %>%
            rename(probe_id = X)
              
str(fc)
## Plot MA for all the expressionset
eset <- read.csv("./data/tidy_data/expressionset_matrix.csv")
## till here

ematrix <- as.matrix(eset[-1])
rownames(ematrix) <- eset$X

x <- t(ematrix)
pc <- prcomp(x)
names(pc)
plot(pc$x[, 1], pc$x[, 2], main = "PCA", xlab = "PC1", ylab = "PC2")

library(rafalib)
gender <- as.fumeric(as.character(a:b))
plot(pc$x[, 1], pc$x[, 2], main = "PCA", xlab = "PC1", ylab = "PC2", col = gender)