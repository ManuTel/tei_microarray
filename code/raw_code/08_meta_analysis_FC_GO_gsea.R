## do all the meta analysis figures for F0_F1_F2 ####
## Load necessary packages ####
source('code/raw_code/00_packages.R')
rm(list =ls())

##load foldchange data -- Make it nest, gsea data
fcdata  <- read.csv("./data/tidy_data/foldchange_df.csv")
fc_nest <- fcdata %>% nest(-group)
head(fcdata)

load("./data/tidy_data/gsea_short.RData")
fcgsea <- fc_gsea_short
rm(fc_gsea_short)
names(fcgsea); str(fcgsea$num_genelist); head(fcgsea$gsebp[[1]])[,1:9]

## volcano plots (increased font size) #### 
# refer the script 03_descriptive_plots.R
# add new boolean column to colour significant FC differently
fcdata <- mutate(fcdata, signFC = as.factor(P.Value < 0.05 & 2^abs(logFC) >= 1.3))
table(fcdata$group, fcdata$signFC)

## f0 plots 

gg <- fcdata %>% 
  filter(group == "f1_maleHS") %>%
  ggplot(aes(x = logFC, y = -log10(P.Value), col = signFC)) +
  geom_point(alpha = 0.5, size = 2, shape = 16) +
  xlim(c(-3,3)) + ylim(c(0,7)) + 
  xlab("log2 fold change") + ylab("-log10 p-value") + 
  scale_color_manual(values = c("grey50", "red3")) +
  theme_minimal() +
  theme(legend.position = "none",
        axis.title     = element_text(size = 18, color = "black"),
        axis.text     = element_text(size = 16, color = "black"))

gg
gg2 <- gg %+% filter(fcdata, group == "f1_femaleHS")
grid.arrange(gg, gg2, ncol = 2)

## F1_F2 plots

gg3 <- fcdata %>% 
  filter(group == "f2_maleLS") %>%
  ggplot(aes(x = logFC, y = -log10(P.Value), col = signFC)) +
  geom_point(alpha = 0.5, size = 2, shape = 16) +
  xlim(c(-2,2)) + ylim(c(0,8)) + 
  xlab("log2 fold change") + ylab("-log10 p-value") + 
  scale_color_manual(values = c("grey50", "red3")) +
  theme_minimal() +
  theme(legend.position = "none",
        axis.title     = element_text(size = 18, color = "black"),
        axis.text     = element_text(size = 16, color = "black"))
gg4 <- gg3 %+% filter(fcdata, group == "f2_femaleLS")
gg4
grid.arrange(gg3, gg4, ncol = 2)

## Done still here ####

