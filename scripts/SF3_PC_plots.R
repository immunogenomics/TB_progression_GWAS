#Author: yangluo
#Date:   2017-12-15T10:42:44-05:00

library(ggplot2)
theme_set(theme_bw(base_size = 16))

# To plot witin sample only PCs -- this is also used in association studies, where including the first 4 PCs as covariates

y<-read.table("data/LIMAA.eigenvec",h=F,stringsAsFactors = F)
eigenval<-scan("data/LIMAA.eigenval",as.numeric())
val<-eigenval/sum(eigenval)*100

df<-data.frame(y)
pheno<-read.table("data/TB.fam",h=F,stringsAsFactors = F)
df$Pheno<-ifelse(pheno[match(df$V2,pheno$V2),6]==1,"Controls","Cases")

ggplot(df,aes(x=V3,y=V4,colour=Pheno))+geom_point()+
  xlab(paste("PC1 (",round(val[1],2),"%)"))+ylab(paste("PC2 (",round(val[2],2),"%)"))+theme(legend.position ="top")+ scale_colour_discrete(name = "Phenotype")
ggsave("figures/SF3b_LIMAA-PCs.pdf",width=8,height=7)


# To plot global PCs merged with 1000 Genomes Phase 1 data
x<-read.table("data/g1k-limaa.eigenvec",h=F,stringsAsFactors = F)
eigenval<-scan("data/g1k-limaa.eigenval",as.numeric())
val<-eigenval/sum(eigenval)*100
df<-data.frame(x)

pop<-read.table("data/g1k.samples",h=T,stringsAsFactors = F)
df$POP<-pop[match(df$V1,pop$sample),]$pop
df$REGION<-pop[match(df$V1,pop$sample),]$super_pop
df[is.na(df$POP),]$REGION<-"LIMAA"
df[is.na(df$POP),]$POP<-"LIMAA"
ggplot(df,aes(x=V3,y=V4,colour=REGION))+geom_point()+
  xlab(paste("PC1 (",round(val[1],2),"%)"))+ylab(paste("PC2 (",round(val[2],2),"%)"))+
  scale_color_manual(values=c("dodgerblue4", "darkolivegreen4","darkorchid3", "goldenrod1","darkred","pink"),
    name="Populations",labels=c("African", "American", "East Asian","European","TB GWAS samples","South Asian"))+
  theme(legend.position="top")
ggsave(filename="figures/SF3a_Global-PCs.pdf",width=8,height=7)


#PCs with pre-caluclated PCs from SNPWT
pcs<-read.table("data/LIMAA-1000g.predpc",h=F,stringsAsFactors = F)
eigenval<-scan("data/LIMAA-1000g.predpc.eval",as.numeric())
val<-eigenval[1:20]/sum(eigenval[1:20])*100
x<-read.table("data/LIMAA-projected.predpc",h=F,stringsAsFactors = F)
g1k<-read.table("data/g1k.samples",h=T,stringsAsFactors = F)
pcs[pcs$V1 %in% g1k[g1k$pop=="PEL",1],6]<-"PEL"
names(pcs)[6]<-"Population"

ggplot(pcs,aes(x=V2,y=V3,colour=Population))+geom_point()+
  xlab(paste("PC1 (",round(val[1],2),"%)"))+ylab(paste("PC2 (",round(val[2],2),"%)"))+
  scale_color_manual(values=c("dodgerblue4", "darkolivegreen4",
                              "darkorchid3", "goldenrod1","darkred","pink"),
                     labels=c("African", "American", "East Asian","European","TB GWAS samples","South Asian"))+
  geom_point(data = x, aes(x=V4,y=V5), colour = alpha("darkred",1))+
  theme(legend.position = "top")
ggsave("figures/SF3c_LIMAA-1000g-snpwt-projected-pcs.pdf",width=8,height=7)

