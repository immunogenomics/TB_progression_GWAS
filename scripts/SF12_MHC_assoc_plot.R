library(ggplot2)
library(dplyr)
library(ggrepel)
library(RColorBrewer)

gwas<-read.table("data/LIMAA_MHC_ALL.assoc.txt.gz",h=T,stringsAsFactors = F)
names(gwas)<-c("CHR","SNP","BP","n_miss","A1","A0","AF","beta","se","l_remle","l_mle","P_wald","P","P_score")

#assigning to HLA genes
gwas$GENE<-"intergenic"

match_index<-function(gene='A',data=gwas,buffer=500){
  
  patterns <- c(paste("SNP",gene,sep="_"), paste("HLA",gene,sep="_"), paste("AA",gene,sep="_"))
  
  matches <- unique (grep(paste(patterns,collapse="|"), 
                          gwas$SNP, value=FALSE))
  
  idx<-gwas$BP>=(gwas[matches[1],]$BP-buffer) & 
    gwas$BP<=(gwas[matches[length(matches)],]$BP+buffer)
  return(idx)
}

genes<-c("A","B","C","DQA1","DQB1","DRB1")

for (gene in genes){
  gwas[match_index(gene=gene,buffer=10000),]$GENE<-gene
}

### Plotting
snpsOfInterest<-gwas[gwas$P==min(gwas$P),]$SNP
for (gene in genes){
  df<-gwas[gwas$GENE==gene,]
  snpsOfInterest<-c(snpsOfInterest,df[df$P==min(df$P),]$SNP[1])
}

#manually defining this
snpsOfInterest <- c('AA_A_73_29910750_exon2_.I','AA_B_97_31324201_exon3_LST')

# Prepare the dataset
don <- gwas %>% 
  
  # Compute chromosome size
  group_by(CHR) %>% 
  summarise(chr_len=max(BP)) %>% 
  
  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(chr_len)-chr_len) %>%
  select(-chr_len) %>%

  # Add this info to the initial dataset
  left_join(gwas, ., by=c("CHR"="CHR")) %>%
  
  # Add a cumulative position of each SNP
  arrange(CHR, BP) %>%
  mutate( BPcum=BP+tot) %>%
  
  # Add highlight and annotation information
  mutate( is_highlight=ifelse(SNP %in% snpsOfInterest, "yes", "no")) %>%
  mutate( is_annotate=ifelse(SNP %in% snpsOfInterest, "yes", "no")) 

cols = c(brewer.pal(6, 'Set1'), '#777777')

# Make the plot
plt<-ggplot(don, aes(x=BPcum, y=-log10(P), color=GENE)) +
  # Show all points
  geom_point()+
  theme_bw() + xlab("Chromosome 6 Position (BP)")+ ylab('-log10(P)')+
  geom_hline(yintercept=5,color="red",linetype="dashed")+
  # Add highlighted points
  geom_point(data=subset(don, is_highlight=="yes"), size=2) +
  
  # Add label using ggrepel to avoid overlapping
  geom_label_repel( data=subset(don, is_annotate=="yes"), aes(label=c("HLA-A Position 73","HLA-B Position 97")), size=2)+
  scale_color_manual("HLA gene", values = cols )
plt 

ggsave("figures/SF12_MHC_association.png",width=10,height=6)
