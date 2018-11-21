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
library(stringr)

## load expressionset matrix and foldchange dataframe ####
## expressionset matrix contains the expression value for each gene while
## foldchange dataframe contains 10 groups (case vs control) differential gene expression profiles
eset <- read.csv("./data/tidy_data/expressionset_matrix.csv")
fc   <- read.csv("./data/tidy_data/foldchange_df.csv")

updown <- fc %>% 
              group_by(group) %>%
              filter(P.Value <  0.05 &
                       2^abs(logFC) >= 1.3) %>%
              summarise(n = n(), up = sum(logFC > 0),
                               down = sum(logFC < 0))

## find the total number of genes that are significantly (p<0.05) differentially expressed (FC>=1.3) in 
## at least one comparison 
fc %>% 
  group_by(group) %>%
  filter(P.Value <  0.05  &
           2^abs(logFC) >= 1.3) %>% # & !str_detect(group, "HS") if we want only LS diet comparison #group != "f0_female"
  ungroup() %>%
  summarise(n = n_distinct(probe_id)) 
# 2793 probesets are differentially expressed in at least one comparison ##
# 2683 probesets if we consider only LS diet in F1 F2
# 2306 probesets if we consider only HS diet in F1 F2
sign_fc <- fc %>%
              group_by(group) %>%
              filter(P.Value <  0.05  &
                       2^abs(logFC) >= 1.3) %>%
              ungroup()
length(unique(sign_fc$probe_id))
length(unique(sign_fc$ID))
fc_2793 <- fc %>%
              filter(probe_id %in% unique(sign_fc$probe_id)) %>%
              select(probe_id, ID, logFC, group) %>%
              spread(group, logFC)

## Save outputs in RData
# save(sign_fc, fc_2793, fc_1620, updown, file = "./data/tidy_data/diff_fc.RData")

## load RData
load("./data/tidy_data/diff_fc.RData")
abc <- sign_fc %>%
          select(probe_id, ID, logFC, group) %>%
              spread(group, logFC) %>%
              mutate_at(vars(starts_with("f")), funs(ifelse( . > 0, 1, -1))) %>%
              mutate_at(vars(starts_with("f")), funs(replace(., is.na(.), 0)))

m <- as.matrix(abc[3:12])
rownames(m) <- abc$probe_id
breaks    <- c(-2,-1,0,1)
library(gplots)
gradient2 <- c("green", "grey55", "red")

dev.new(width=4, height=12)
heatmap.2(m,scale="none",breaks=breaks,col=gradient2, Colv = FALSE,
          trace="none", dendrogram = 'none',cexRow=0.6,cexCol=0.8,lhei = c(0.05,5),
          margin=c(7,10),lwid=c(4,2.0))

abc2 <- abc %>% select(-contains("HS"))
m2 <- as.matrix(abc2[4:8])
m2 <- m2[rowSums(m2) !=0, ]

## heatmap conclusion ####
## it doesn't give any insides into data structure and useless ##


## volcano plots ##
# add new boolean column to colour significant FC differently
fc <- mutate(fc, signFC = as.factor(P.Value < 0.05 & 2^abs(logFC) >= 1.3))
table(fc$group, fc$signFC)

#f0 plots 
gg <- fc %>% 
        filter(group == "f0_male") %>%
          ggplot(aes(x = logFC, y = -log10(P.Value), col = signFC)) +
            geom_point(alpha = 0.5, size = 2, shape = 16) +
            xlim(c(-8,8)) + ylim(c(0,15)) + 
            xlab("log2 fold change") + ylab("-log10 p-value") + 
            scale_color_manual(values = c("grey50", "red3")) +
            theme_minimal() +
            theme(legend.position = "none",
                  axis.title     = element_text(size = 14, color = "black"),
                  axis.text     = element_text(size = 12, color = "black"))
gg
gg2 <- gg %+% filter(fc, group == "f0_female")
# f1-f2 plots        minFC         maxFC
# f1_femaleHS       -2.30         1.07 
# f1_femaleLS       -1.14         2.52 
# f1_maleHS         -2.04         0.795
# f1_maleLS         -1.49         2.57 
gg3 <- fc %>% 
  filter(group == "f2_maleHS") %>%
  ggplot(aes(x = logFC, y = -log10(P.Value), col = signFC)) +
  geom_point(alpha = 0.5, size = 2, shape = 16) +
  xlim(c(-4,4)) + ylim(c(0,8)) + 
  xlab("log2 fold change") + ylab("-log10 p-value") + 
  scale_color_manual(values = c("grey50", "red3")) +
  theme_minimal() +
  theme(legend.position = "none",
        axis.title     = element_text(size = 14, color = "black"),
        axis.text     = element_text(size = 12, color = "black"))
gg4 <- gg3 %+% filter(fc, group == "f2_femaleHS")
gg4
grid.arrange(gg3, gg4, ncol = 2)

