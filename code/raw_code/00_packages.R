##install necessary packages if not installed

list_packages <- c("ggplot2", "dplyr", "purrr", "ggpubr", "ggsci", "gridExtra",
                       "broom", "tidyr", "tibble")
new_packages <- list_packages[!(list_packages %in% installed.packages()[,"Package"])]
if(length(new_packages)) install.packages(new_packages,verbose=FALSE,quiet=TRUE)

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
library(tibble)