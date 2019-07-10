library(ggplot2)
library(reshape)
library(tidyverse)
source("scripts/qqplot.R")

gdf <- read.table("data/crispr_diff_expr_6genes.txt",h=T,stringsAsFactors = F)
p_summary <- c(0.8294005,	0.9162661	,0.8296861,	0.2811773,	0.9675479,	0.6502330)
pltDf <- melt(gdf)
label <- tibble(
  Phenotype = Inf,
  value = Inf,
  variable = unique(pltDf$variable),
  label = str_c("Anova P = ", round(p_summary,3))
)

p <- ggplot(pltDf,aes(x=Phenotype,y=value)) + geom_boxplot() + geom_point(position = position_dodge(width=0.8) , color="red")+ 
  geom_text(aes(label = label), data = label, vjust= 2, hjust = 1.8, size =4) 
p + facet_wrap(. ~ variable ,nrow = 2, scales = "free_y") +ylab("log2(TPM+1)") + xlab("") +
  scale_y_continuous(expand = expand_scale(mult = c(0, .3))) + theme_bw()

ggsave(filename='figures/SF17a_diff_exression_6genes.pdf',width=8,height=8)

#Volcano plot
result <- read.table('data/single_clones_average_expression_volcano_results.txt.gz',header=T,stringsAsFactors = F)

result$logP <- -log10(result$Pval)

# plot QQ
qqunif.plot(out$P)

ggplot() +
  geom_point(
    data = result, 
    mapping = aes_string(x = "log2FC", y = "logP"),
    size = 0.5, stroke = 0.1, shape = 21
  )  + xlim(-2,2) + 
  geom_hline(yintercept = -log10(0.05/nrow(result)),  color="darkred", linetype = 2,size=0.5) +
  geom_vline(aes(xintercept=1.5), color="darkred", linetype=2, size=0.5) +
  geom_vline(aes(xintercept=-1.5), color="darkred", linetype=2, size=0.5) +
  labs(x = "log2(FC)", y = "-log10(P-value)" ) +
  theme_classic(base_size = 12) +
  theme(panel.grid = element_blank())

ggsave(filename='figures/SF17b_volcano.pdf',width=8,height=8)
