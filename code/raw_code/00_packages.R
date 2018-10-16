##install necessary packages if not installed

list.of.packages <- c("ggplot2", "dplyr", "purrr", "ggpubr", "ggsci", "gridExtra",
                      "affy", "limma", "drosophila2.db", "drosophila2cdf", "broom", "tidyr")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages,verbose=FALSE,quiet=TRUE)

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