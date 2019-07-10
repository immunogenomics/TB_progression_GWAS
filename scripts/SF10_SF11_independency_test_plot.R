library(ggplot2)

theme_set(theme_bw(base_size=20))

df <- read.table("data/beta_correlations.txt.gz",h=T, stringsAsFactors = F)


p<- ggplot()+geom_point(data=df,aes(x=beta.x,y=beta.y),colour="black")+xlab("case-control beta")+ylab("case-only beta")+theme_bw()+ 
  geom_point(data=df[df$rsid=="AX-14134648",],aes(x=beta.x,y=beta.y),colour="red",size=2)

p<-p+ annotate("text", x = -0.08, y = 0.25, label = paste("Pearson correlation = ",round(cor(df$beta.x,df$beta.y),3)),size=3.5)
p+annotate("text", x = 0.2, y = 0.09, label = "rs73226617",size=3.5,colour="red")
ggsave("figures/SF10_beta-corr.png",width=6,height=6)


perm<-read.table("data/case_perm_assoc.txt.gz",h=T,stringsAsFactors = F)

true<-data.frame(beta=8.438318e-02,se=3.487196e-02)
pval <- 1-sum(abs(true$beta)>abs(perm$beta))/nrow(perm)

p_beta<-ggplot(perm,aes(x=beta))+geom_histogram(binwidth=.005,color="black", fill="white")

p_beta<-p_beta+ geom_vline(aes(xintercept=true$beta),
                           color="red", size=1)+xlim(-0.2,0.2)
p_beta+theme_bw()+xlab("Effect size") + ylab("Counts")+
  annotate("text", x = 0.065, y = 600, label = paste("p =",round(pval,3)),size=4.5,colour="red")

ggsave("figures/SF11_permutation_pvalue.pdf",width=10,height=6.5)

