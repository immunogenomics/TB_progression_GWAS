# some pre-settings for plotting using pophelper package (http://royfrancis.github.io/pophelper)
#install.packages(c("Cairo","devtools","ggplot2","gridExtra","gtable","tidyr"),dependencies=T)
#devtools::install_github('royfrancis/pophelper')
# load library
library(pophelper)

#install.packages("seriation")
library(seriation)
library(gridExtra)


#plotting parameters
cb_paired=c("darkred","#1F78B4","#B2DF8A","#33A02C","#FB9A99","#E31A1C","#FDBF6F","#FF7F00","#CAB2D6","#6A3D9A","#FFFF99","#B15928","#F2F3F4","#222222","#F3C300")

##reordering clusters (mannually for input files)
# for (k in c(3:15)){
#   f<-read.table(paste("data/clean-TB.newids.fixed-Reich.newids.fixed-1kg-pruned.",k,".Q",sep=""),h=F,stringsAsFactors = F)
#   tb_index<-which(grepl("M00",fam$V2)==TRUE)
#   write.table(f[,order(colMeans(f[tb_index,]),decreasing = T)],file=paste("data/admixture.",k,".Q",sep=""),row.names = F,col.names = F,quote=F)
# }  


# read in cluster files 
qfiles <- list.files(path = "data/",pattern = "admixture",full.names = T)
slist<-readQ(file=qfiles[c(7:13,1:6)],filetype="basic")                     # automatically detects input filetype

# create labels in the same order as the cluster file
fam<-read.table("data/SF1_admxiture_samples.txt",h=T,sep="\t",stringsAsFactors = F)

#rename
fam$super_pop<-ifelse(fam$super_pop=="Equatorial-Tucanoan" | fam$super_pop=="Ge-Pano-Carib","Southern\nAmerindian",fam$super_pop)
fam$super_pop<-ifelse(fam$super_pop=="Central-Amerind","Central\nAmerindian",fam$super_pop)
fam$super_pop<-ifelse(fam$super_pop=="Northern Amerind","Northern\nAmerindian",fam$super_pop)

pop_label<-data.frame(Cohort=fam$cohort,Region=fam$super_pop,pop=fam$pop,status=fam$V6,stringsAsFactors = F)
sapply(pop_label,is.character)


# assign individual ids
if(length(unique(sapply(slist,nrow)))==1) slist <- lapply(slist,"rownames<-",fam$V2)

  
#selected samples only
selected_pop<-which(attributes(slist[[1]])$row.names %in% fam[fam$super_pop!="others" & fam$super_pop!="TB",2])
selected_TB<-sample(which(fam$super_pop=="TB"),1013)
selected_index<-c(selected_pop,selected_TB)
slist1 <- sapply(slist,function(x) x[c(selected_index),])
pop_label1<-pop_label[selected_index,]
dim(pop_label1)

#grouping pops to super groups 
pop_label1[pop_label1$pop%in% c("Cree","Ojibwa","Algonquin"),]$pop<-"Other\nNorthern\nAmerindian"
pop_label1[pop_label1$Region=="Northern\nAmerindian",]$pop<-"Northern\nAmerindian"
pop_label1[pop_label1$Region=="Central\nAmerindian",]$pop<-"Central\nAmerindian"
pop_label1[pop_label1$Region=="Southern\nAmerindian",]$pop<-"Southern\nAmerindian"
pop_label1[pop_label1$Region=="Andean" & pop_label1$pop!="Quechua"  & pop_label1$pop!="Aymara",]$pop<-"Andean\n"
pop_label1[pop_label1$pop=="Quechua",]$pop<-"Andean"
pop_label1[pop_label1$pop=="Aymara",]$pop<-"Andean"
pop_label1[pop_label1$Cohort=="Reich et al",]$Cohort<-"NAT"
pop_label1[pop_label1$pop=="TB" & pop_label1$status==2,]$pop<-"Cases"
pop_label1[pop_label1$pop=="TB" & pop_label1$status==1,]$pop<-"Controls"

#single plot
selected_pops<-c("CHB","YRI","CEU","PUR","CLM","MXL","PEL","Northern\nAmerindian","Central\nAmerindian","Southern\nAmerindian",
                 "Andean","Cases","Controls")

nspl=round(sum(pop_label1$pop%in%selected_pops)/1)

#multiple plots
for (i in c(1:13)){
  assign(paste("p",i,sep=""),plotQMultiline(slist1[i],spl=nspl,useindlab=F,showindlab=F,basesize=10,width=100,height=20,
                                            barsize=1,grplab=pop_label1[,c(1,3)],ordergrp=TRUE,selgrp="pop",subsetgrp=selected_pops,
                                            grplabsize = 20,sort="all",clustercol = cb_paired,returnplot=F,exportplot=T,
                                            showtitle=T,showsubtitle=F,titlelab=paste0("ADMIXTURE K=",sapply(slist1[i],ncol)),titlesize=36,titlehjust = 0.5) )
  
}
#g<-grid.arrange(p1$plot[[1]][[1]],p2$plot[[1]][[1]],p3$plot[[1]][[1]],p4$plot[[1]][[1]],p5$plot[[1]][[1]],p6$plot[[1]][[1]],
#                p7$plot[[1]][[1]],p8$plot[[1]][[1]],p9$plot[[1]][[1]],p10$plot[[1]][[1]],p11$plot[[1]][[1]],p12$plot[[1]][[1]],p13$plot[[1]][[1]],nrow=13,
#                top="ADMIXTURE plots")
#ggsave("figures/SF1_LIMAA_ADMIXTURE.png",g,height=40,width=50,limitsize=FALSE)

#for displaying in the main figure, pick K=6
# same length as the reference in previous step
#selected_TB<-sample(which(fam$super_pop=="TB"),1013)

i=4
k6<-plotQMultiline(slist1[i],spl=1013,useindlab=F,showindlab=F,basesize=10,width=80,height=20,
                   barsize=1,grplab=pop_label1[,c(1,3)],ordergrp=TRUE,selgrp="pop",subsetgrp=selected_pops,
                   grplabsize = 24,sort="all",clustercol = cb_paired,returnplot=F,exportplot = T,
                   showtitle=T,showsubtitle=F,titlelab=paste0("ADMIXTURE K=",sapply(slist1[i],ncol)),titlesize=36,titlehjust = 0.5,imgtype="pdf",titlespacer =1.5) 
#ggsave("figures/F2a_LIMAA_ADMIXTURE.png",k6,height=40,width=50,limitsize=FALSE)

