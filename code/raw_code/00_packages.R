##install necessary packages if not installed

list.of.packages <- c("ggplot2", "dplyr", "purrr", "ggpubr", "ggsci", "gridExtra",
                      "affy", "limma", "drosophila2.db", "drosophila2cdf", "broom")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages,verbose=FALSE,quiet=TRUE)
