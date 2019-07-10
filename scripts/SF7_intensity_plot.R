library(ggplot2)

df <- read.table("data/rs73226617_intensity.txt",header = T, stringsAsFactors = F)
#all

ggplot(df[df$dosage!=-1,],aes(x=A,y=B,colour=as.factor(dosage),shape=as.factor(dosage)))+geom_point(size=2.5)+xlim(0,3500)+ylim(0,3500)+
  theme_bw()+theme(legend.position=c(0.8,0.3),legend.text = element_text(size=14),legend.title=element_text(size=16),plot.title = element_text(hjust=0.5, face="bold",size=16))+
  scale_color_manual(name="Genotypes",values=c("red","orange","blue"),labels=c("AA","AG","GG"))+
  guides(colour = guide_legend(override.aes=list(shape=c(16,17,15),size=3))) + scale_shape(guide=FALSE)

ggsave(filename='figures/SF7_rs73226617_intensity_plot.pdf',width=10,height=7)

