library(ggplot2)
library(ggrepel)

hits<-read.table("data/top_snps.txt",h=T,stringsAsFactors = F)
header<-c("chr:pos_ref_alt","SNP","phenotypeID",
          "P","beta","Bonferroni.p.value","FDR","alt_allele_frequenc","std.error_of_beta")

### Monocytes
mono_k4me1<-read.table("data/mono_K4ME1_hit.txt.gz",h=F,stringsAsFactors = F)
names(mono_k4me1)<-header
dim(mono_k4me1)
dat<-mono_k4me1
dat$CHR<-as.numeric(lapply(dat[,1],function(x) unlist(strsplit(x,":"))[[1]]))

pos<-lapply(dat[,1],function(x) unlist(strsplit(x,":"))[[2]])

dat$BP<-as.numeric(lapply(pos, function(x) unlist(strsplit(x,"_"))[[1]]))

#check whether the p-val is the lowest among all k4me1 regiones
dat$is_min<-0

rsids<-unique(dat[dat$SNP %in% hits$rsID,]$SNP)

for (rsid in rsids){
  print(rsid)
  df<-dat[dat$SNP==rsid,]
  df2<-df[df$P==min(df$P),]
  dat[dat$SNP==rsid & dat$phenotypeID==df2$phenotypeID,]$is_min<-1
}
# Prepare the dataset for plotting
don_mono <- dat %>% 
  
  # Compute chromosome size
  group_by(CHR) %>% 
  summarise(chr_len=max(BP)) %>% 
  
  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(chr_len)-chr_len) %>%
  select(-chr_len) %>%
  
  # Add this info to the initial dataset
  left_join(dat, ., by=c("CHR"="CHR")) %>%
  
  # Add a cumulative position of each SNP
  arrange(CHR, BP) %>%
  mutate( BPcum=BP+tot) %>%
  
  # Add highlight and annotation information
  mutate( is_highlight=ifelse(SNP %in% hits$rsID, "yes", "no")) %>%
  mutate( is_annotate=ifelse(is_min==1, "yes", "no")) 

## Neutrophil

neut_k4me1<-read.table("data/neut_K4ME1_hit.txt.gz",h=F,stringsAsFactors = F)
names(neut_k4me1)<-header
dim(neut_k4me1)
dat<-neut_k4me1
dat$CHR<-as.numeric(lapply(dat[,1],function(x) unlist(strsplit(x,":"))[[1]]))

pos<-lapply(dat[,1],function(x) unlist(strsplit(x,":"))[[2]])

dat$BP<-as.numeric(lapply(pos, function(x) unlist(strsplit(x,"_"))[[1]]))

#check whether the p-val is the lowest among all k4me1 regiones
dat$is_min<-0

rsids<-unique(dat[dat$SNP %in% hits$rsID,]$SNP)

for (rsid in rsids){
  print(rsid)
  df<-dat[dat$SNP==rsid,]
  df2<-df[df$P==min(df$P),]
  dat[dat$SNP==rsid & dat$phenotypeID==df2$phenotypeID,]$is_min<-1
}

# Prepare the dataset for plotting
don_neut <- dat %>% 
  
  # Compute chromosome size
  group_by(CHR) %>% 
  summarise(chr_len=max(BP)) %>% 
  
  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(chr_len)-chr_len) %>%
  select(-chr_len) %>%
  
  # Add this info to the initial dataset
  left_join(dat, ., by=c("CHR"="CHR")) %>%
  
  # Add a cumulative position of each SNP
  arrange(CHR, BP) %>%
  mutate( BPcum=BP+tot) %>%
  
  # Add highlight and annotation information
  mutate( highlight=ifelse(SNP %in% hits$rsID, "yes", "no")) %>%
  mutate( is_annotate=ifelse(is_min==1, "yes", "no")) 



# CD4+ T cells
tcel_k4me1<-read.table("data/tcel_K4ME1_hit.txt.gz",h=F,stringsAsFactors = F)
names(tcel_k4me1)<-header
dat<-tcel_k4me1
dat$CHR<-as.numeric(lapply(dat[,1],function(x) unlist(strsplit(x,":"))[[1]]))

pos<-lapply(dat[,1],function(x) unlist(strsplit(x,":"))[[2]])

dat$BP<-as.numeric(lapply(pos, function(x) unlist(strsplit(x,"_"))[[1]]))

#check whether the p-val is the lowest among all k4me1 regiones
dat$is_min<-0

rsids<-unique(dat[dat$SNP %in% hits$rsID,]$SNP)

for (rsid in rsids){
  print(rsid)
  df<-dat[dat$SNP==rsid,]
  df2<-df[df$P==min(df$P),]
  dat[dat$SNP==rsid & dat$phenotypeID==df2$phenotypeID,]$is_min<-1
}

# Prepare the dataset for plotting
don_tcel <- dat %>% 
  
  # Compute chromosome size
  group_by(CHR) %>% 
  summarise(chr_len=max(BP)) %>% 
  
  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(chr_len)-chr_len) %>%
  select(-chr_len) %>%
  
  # Add this info to the initial dataset
  left_join(dat, ., by=c("CHR"="CHR")) %>%
  
  # Add a cumulative position of each SNP
  arrange(CHR, BP) %>%
  mutate( BPcum=BP+tot) %>%
  
  # Add highlight and annotation information
  mutate( is_highlight=ifelse(SNP %in% hits$rsID, "yes", "no")) %>%
  mutate( is_annotate=ifelse(is_min==1, "yes", "no")) 

# Make combined plot
g_combined<-ggplot(don_tcel, aes(x=BPcum, y=-log10(P))) + 
  
  # Show all points
  geom_point( color="skyblue", alpha=0.9, size=1) +
  
  #add another layer
  # neutrophile
  geom_point(data=don_neut, colour="darkgrey", alpha=0.9, size=1) +
  
  # monocyte
  geom_point(data=don_mono, colour="rosybrown2", alpha=0.9, size=1) +
  
  # custom X axis:
  #scale_x_continuous( label = axisdf$CHR, breaks= axisdf$center ) +
  scale_y_continuous(expand = c(0, 1) ) +     # remove space between plot area and x axis
  
  # Add highlighted points
  geom_point(data=subset(don_tcel, is_annotate=="yes"), color="blue", size=2) +
  
  geom_point(data=subset(don_neut, is_annotate=="yes"), color="black", size=2) +
  
  geom_point(data=subset(don_mono, is_annotate=="yes"), color="brown3", size=2) +
  
  # Add label using ggrepel to avoid overlapping
  geom_label_repel( data=subset(don_tcel, is_annotate=="yes"), aes(label=SNP), size=2) +
  
  # Custom the theme:
  theme_bw() + xlab("Position") + 
  theme( 
    legend.position="topleft",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  )

g_combined+ annotate(geom="text", x=141210000, y=22, col="rosybrown", label="Monocytes", parse=T)+
  annotate(geom="text", x=141210000, y=21, col="grey", label="Neutrophil", parse=T)+
  annotate(geom="text", x=141212000, y=20, col="skyblue", label="CD4+ T-cell", parse=T)

ggsave("figures/SF13_chromQTL.png",width=10,height=6)
