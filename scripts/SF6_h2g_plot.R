library(ggplot2)
library(ggsignif)

# original estimates
limaa_h2<-data.frame(h2=0.212,se=0.08,N=3212)
limaa_clustered_h2<-data.frame(h2=0.221,se=0.06,N=2535)
russian_ldsc<-data.frame(h2=0.155,se=0.04,N=11137)
russian_gcta<-data.frame(h2=0.178,se=0.02,N=11137)

df<-rbind(limaa_h2,limaa_clustered_h2,russian_ldsc,russian_gcta)

df<-cbind(df,Study=c("LIMAA TB progression (3,212)","LIMAA rapid TB progression (2,535)",
                     "Russian LDSC(11,137)","Russian GCTA(11,137)"))

df$pval<-2*pnorm(-abs(df$h2/df$se))

#plot 
myColors <- c("darkred","orange","mediumpurple","navyblue")
names(myColors) <- levels(df$Study)
colScale <- scale_fill_manual(name = "Study",values = myColors)

p<-ggplot(df,aes(x=Study,y=h2,fill=Study))+geom_bar(stat="identity",position=position_dodge())+
  geom_errorbar(aes(ymin=h2-se,ymax=h2+se),width=.2,position=position_dodge(.9))
p1<-p+theme_bw()+scale_fill_manual(values=myColors)+ylim(0,.5)+
  ylab(expression(h[g]^2))+theme(legend.position="none")+xlab("")

my_comparisons <- list( c("LIMAA TB progression (3,212)","LIMAA rapid TB progression (2,535)"), 
                        c("LIMAA TB progression (3,212)", "Russian GCTA(11,137)") )

#define p-value function for comparing two h2g estimates
compare<-function(df=stats){
  #standard deviation
  df$sd1i<-with(df,se1i*sqrt(n1i))
  df$sd2i<-with(df,se2i*sqrt(n2i))
  
  df$d<-with(df,m1i-m2i)
  df$s<-with(df,sqrt(((n1i-1)*sd1i^2 + (n2i-1)*sd2i^2)/(n1i + n2i - 2)))
  
  df$v<-with(df,s^2*(1/n1i+1/n2i))
  df$v2<-with(df,sd1i^2/n1i+sd2i^2/n2i)
  df$se<-with(df,sqrt(v)) #assuming same population deviance
  df$se2<-with(df,sqrt(v2)) #assuming different deviance
  df$tval<-with(df,d/se)
  df$tval2<-with(df,d/se2)
  df$pval <- 2 * with(df, pt(abs(d/se), df = n1i + n2i - 2, lower.tail=FALSE))
  df$pval2 <- 2 * with(df, pt(abs(d/se2), df = n1i + n2i - 2, lower.tail=FALSE))
  df$sig <- ifelse(df$pval<0.05/3,1,0)
  df$sig2<-ifelse(df$pval2<0.05/3,1,0)
  
  return(df)
}

# progression vs rapid progression (clustered cases only)
stats<-data.frame(Comparison="progression vs rapid progression",
                  m1i=df[1,]$h2,se1i=df[1,]$se,
                  m2i=df[2,]$h2,se2i=df[2,]$se,
                  n1i=df[1,]$N,n2i=df[2,]$N)

s1<-compare(stats)

# TB progression vs Russian GWAS GCTA
stats<-data.frame(Comparison="LIMAA vs Russian",
                     m1i=df[df$Study=="LIMAA TB progression (3,212)",]$h2,
                  se1i=df[df$Study=="LIMAA TB progression (3,212)",]$se,
                           m2i=df[4,]$h2,se2i=df[4,]$se, 
                  n1i=df[1,]$N,n2i=df[4,]$N)

s2<-compare(stats)

p1+ 
  geom_signif(y_position=c(0.35, 0.45), xmin=c(1, 1), xmax=c(2, 3),
              annotation=c(paste("t-test P-value =",round(s1$pval2,2)), paste("t-test P-value = ",round(s2$pval2,2))), tip_length=.1)

ggsave(filename='figures/SF6_heritability_estimates.pdf',width=10,height=7)
