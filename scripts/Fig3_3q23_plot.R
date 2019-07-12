library('Sushi')
library('biomaRt')
source("scripts/plotManhattan2.R")

# define zoom region
chrom1="chr3" 
chromstart1=141383625-60000
chromend1=141406933+60000

# A) zoom in of the 3q23 regionls
data<-read.table("data/3q23.assoc",h=T,stringsAsFactors = F)
df<-data[,c(1,3,2,13)]
names(df)<-c("CHR","BP","SNP","P")
gwas<-data.frame(chr.hg19=paste("chr",df$CHR,sep=""),pos.hg19=as.numeric(df$BP),pos.hg19.1=as.numeric(df$BP),rsid=df$SNP,pval=df$P)

# ld info
ld<-read.table("data/rs73226617.ld",h=T,stringsAsFactors = F)
gwas$ld<-0
gwas$ld<-abs(ld[match(gwas$rsid,ld$rsid),2])
#gwas[is.na(gwas$ld),]$ld<-0
#add gene names
chrom_biomart = gsub("chr","",chrom1)
mart=useMart(host='grch37.ensembl.org', biomart='ENSEMBL_MART_ENSEMBL',
             dataset='hsapiens_gene_ensembl')
a<-getBM(attributes=c("chromosome_name","external_gene_name","ensembl_transcript_id",
                      "transcript_length","transcript_gencode_basic","gene_biotype"),filters= c("chromosome_name","start","end"),
         values=list(chrom_biomart,chromstart1,chromend1), mart)
b<-na.omit(a[a$transcript_gencode_basic==1 & (a$gene_biotype=="lincRNA" |a$gene_biotype=="protein_coding"),])
selected_transcripts<-aggregate(transcript_length~external_gene_name+ensembl_transcript_id,data=b,max)$ensembl_transcript_id
b <- getBM(c("ensembl_transcript_id", "exon_chrom_start","exon_chrom_end","strand"),
           "ensembl_transcript_id",selected_transcripts, mart)
d<-merge(a,b)

geneinfobed<-data.frame(chrom=d$chromosome_name,start=d$exon_chrom_start,stop=d$exon_chrom_end,gene=d$external_gene_name,score='.',strand=d$strand)
geneinfobed[,1] = paste("chr",geneinfobed[,1],sep="")

#zoom in functional region
chromstart2=141383625-1000
chromend2=141406933+1000

#Finemap
finemap<-read.table("data/finemap.bed",h=F,stringsAsFactors = F)

#all 400 features 
features <- read.table("data/EnhancerRegionTiled_overlap.txt", header = T, stringsAsFactors = F)
features.bedgraph<-cbind(region[,1:3],features[match(region$V2,features$start),2])

names(features.bedgraph)<-c("chrom","start","end","value")


#IMPACT
region<-read.table("data/impact_region.txt",h=F,stringsAsFactors = F)

impact<-read.table("data/Enhancer_Region_IMPACTAnnotations_CEPBP_values.txt",h=T,stringsAsFactors = F)
region$cebpb<-impact[match(region$V2,impact$start),2]

impact_bed<-data.frame(chrom=region$V1,start=region$V2,end=region$V3,cebpb=region$cebpb)

#define 11 risk variants
rs73226617=c(141400653,141400654)
indel<-c(141401146,141401151)

rs73239724<-c(141383625,141383626)
rs73226608<-c(141388147,141388147)
rs58538713<-c(141389447,141389448)
rs11710569<-c(141390125,141390126)
rs11714221<-c(141390219,141390220)
rs189348793<-c(141391318,141391319)
rs73226619<-c(141403952,141403953)
rs112304167<-c(141404369,141404370)
rs146526750<-c(141406933,141406934)

par(family = '')

pdf("figures/Fig3a.pdf",height=14,width=16)
layout(matrix(c(1,
                1,
                1,
                1,
                2,
                3,
                3,
                4,
                4,
                4,
                5,
                5,
                5,6
),nrow=14, byrow = TRUE))

#regional manhattan
par(mar=c(1,6,5, 2 ) )
cex.val <- 1.2 + (gwas$pval < 1E-5)*.5 - (gwas$pval > 1E-5)*.5
plotManhattan2(bedfile=gwas,chrom=chrom1,chromstart = chromstart1,chromend = chromend1,
               pvalues=(gwas$pval),ld=gwas$ld,col="blue",cex=cex.val,ylim=c(0,8))

axis(side=2,las=2,tcl=.5,cex.axis=1.5,at=c(0,2,4,6,8),labels = c("1",expression(10^{-2}),expression(10^{-4}) ,expression(10^{-6}),expression(10^{-8})))
mtext("P-value",line=3.5,side=2,cex=1.5,font=2)
#zoom in signal region
labelplot("(a). Regional Manhattan\n",letteradj=-.06,lettercex=1.5)
zoombox(passthrough =TRUE,topextend = -1,lty = 1)
zoomsregion(region=c(chromstart2,chromend2),chrom=chrom1,extend=c(0.05,1),col="#BEBEBE4C")
#zoom in signal region
par(mar=c(4,6,0, 2 ) )
pg = plotGenes(geneinfobed,chrom1,chromstart1,chromend1 ,labeloffset=0.3,wigglefactor=0.01,
               labeltext=TRUE,maxrows=1,bheight=0.05,plotgenetype="box",bentline=FALSE,
               fontsize=1.8,labelat="end")
#labelgenome( chrom1, chromstart1,chromend1,n=3,scale="Mb",scalecex=1.5,chromcex=1.5)

#zoom in signal region
#par(mar=c(0,4,0, 2 ) )
zoombox(zoomregion=c(chromstart2,chromend2),lty = 1)
zoomsregion(region=c(chromstart2,chromend2),chrom=chrom1,
            extend=c(1,-0.1),wideextend = 1.5,offsets=c(0.02,0.02),col="#BEBEBE4C")



#Finemap
par(mar=c(4,6,3, 4) )
plotBedgraph(signal=finemap,chrom1,chromstart2,chromend2,color = "black",lwd=3,range=c(0,0.6),cex=1.5,cex.axis=2,cex.lab=2)
axis(side=2,las=2,tcl=.5,cex.axis=2,at=c(0,0.3,0.6),labels=c("0","0.3","0.6"))
par(mar=c(2,5,3, 4) )
mtext("Probability",side=2,line=3.5,cex=1.5,font=2)
zoombox(passthrough = TRUE,topextend = -1,lty = 1)
zoomsregion(region=rs73226617,chrom=chrom1,genome=NULL, extend=c(0.07,2),offsets=c(0,0),col="#FA80727F",highlight = TRUE,wideextend=0) #red
#zoomsregion2(region=indel,chrom=chrom1,genome=NULL, extend=c(0.07,2),offsets=c(0,0.0),highlight = FALSE,wideextend=0,lwd=1,zoomborder = NULL)#,col="#ADD8E67F") #blue
zoomsregion3(region=indel,top=finemap[finemap$V2<indel[2] & finemap$V3>indel[1],4],chrom=chrom1,highlight=TRUE, extend=c(0.07,2),offsets=c(0,0.5))#,col="#ADD8E67F") #blue

zoomsregion3(region=rs58538713,top=finemap[finemap$V2<rs58538713[2] & finemap$V3>rs58538713[1],4],chrom=chrom1,genome=NULL, extend=c(0.07,2),offsets=c(0,0),highlight = TRUE)#,col="#FFDAB97F")
zoomsregion3(region=rs73239724,top=finemap[finemap$V2<rs73239724[2] & finemap$V3>rs73239724[1],4],chrom=chrom1,genome=NULL, extend=c(0.07,2),offsets=c(0,0),highlight = TRUE)
zoomsregion3(region=rs73226608,top=finemap[finemap$V2<=rs73226608[2] & finemap$V3>=rs73226608[1],4],chrom=chrom1,genome=NULL, extend=c(0.07,2),offsets=c(0,0),highlight = TRUE)
zoomsregion3(region=rs11710569,top=finemap[finemap$V2<rs11710569[2] & finemap$V3>rs11710569[1],4],chrom=chrom1,genome=NULL, extend=c(0.07,2),offsets=c(0,0),highlight = TRUE)
zoomsregion3(region=rs11714221,top=finemap[finemap$V2<rs11714221[2] & finemap$V3>rs11714221[1],4],chrom=chrom1,genome=NULL, extend=c(0.07,2),offsets=c(0,0),highlight = TRUE)
zoomsregion3(region=rs189348793,top=finemap[finemap$V2<rs189348793[2] & finemap$V3>rs189348793[1],4],chrom=chrom1,genome=NULL, extend=c(0.07,2),offsets=c(0,0),highlight = TRUE)
zoomsregion3(region=rs73226619,top=finemap[finemap$V2<rs73226619[2] & finemap$V3>rs73226619[1],4],chrom=chrom1,genome=NULL, extend=c(0.07,2),offsets=c(0,0),highlight = TRUE)
zoomsregion3(region=rs112304167,top=finemap[finemap$V2<rs112304167[2] & finemap$V3>rs112304167[1],4],chrom=chrom1,genome=NULL, extend=c(0.07,2),offsets=c(0,0),highlight = TRUE)
zoomsregion3(region=rs146526750,top=finemap[finemap$V2<rs146526750[2] & finemap$V3>rs146526750[1],4],chrom=chrom1,genome=NULL, extend=c(0.07,2),offsets=c(0,0),highlight = TRUE)

labelplot("(b). Finemap\n",letteradj=-.06,lettercex=1.5)


par(mar=c(4,6,3, 4) )
#ALL epigenetic features
irf.bedgraph<-features.bedgraph[features.bedgraph$value>0,]
irf.bedgraph$value<-log10(features.bedgraph[features.bedgraph$value>0,]$value)
plotBedgraph(signal=irf.bedgraph,chrom1,chromstart2,chromend2,transparency=0.5,color = "grey",lwd=1,range=c(0,log10(800)))
zoombox(zoomregion = c(chromstart2,chromend2),passthrough = TRUE,lty = 1)
zoomsregion(region=rs73226617,chrom=chrom1,highlight= TRUE, extend=c(2,2),offsets=c(0,0),col="#FA80727F") #red
zoomsregion2(region=indel,chrom=chrom1,highlight=TRUE, extend=c(2,2),offsets=c(0,0.5))#,col="#ADD8E67F") #blue
labelplot("(c). Epigenomic features overlap",letteradj=-.075,lettercex=1.5)
zoomsregion2(region=rs58538713,chrom=chrom1,genome=NULL, extend=c(2,2),offsets=c(0,0),highlight = TRUE)#,col="#FFDAB97F")
zoomsregion2(region=rs73239724,chrom=chrom1,genome=NULL, extend=c(2,2),offsets=c(0,0),highlight = TRUE)
zoomsregion2(region=rs73226608,chrom=chrom1,genome=NULL, extend=c(2,2),offsets=c(0,0),highlight = TRUE)
zoomsregion2(region=rs11710569,chrom=chrom1,genome=NULL, extend=c(2,2),offsets=c(0,0),highlight = TRUE)
zoomsregion2(region=rs11714221,chrom=chrom1,genome=NULL, extend=c(2,2),offsets=c(0,0),highlight = TRUE)
zoomsregion2(region=rs189348793,chrom=chrom1,genome=NULL, extend=c(2,2),offsets=c(0,0),highlight = TRUE)
zoomsregion2(region=rs73226619,chrom=chrom1,genome=NULL, extend=c(2,2),offsets=c(0,0),highlight = TRUE)
zoomsregion2(region=rs112304167,chrom=chrom1,genome=NULL, extend=c(2,2),offsets=c(0,0),highlight = TRUE)
zoomsregion2(region=rs146526750,chrom=chrom1,genome=NULL, extend=c(2,2),offsets=c(0,0),highlight = TRUE)

axis(side=2,las=2,tcl=.5,cex.axis=2,at=c(0,1,2,log10(515)),labels = c("0","10","100","515"))
par(mar=c(2,5,3, 4) )
mtext("No. of features",side=2,line=3.5,cex=1.5,font=2)

#IMPACT
par(mar=c(2,6,3, 4) )
plotBedgraph(signal=impact_bed,chrom1,chromstart2,chromend2,transparency=.50,color = "grey",linecolor = NA)
zoombox(zoomregion = c(chromstart2,chromend2),passthrough = FALSE,lty = 1)
zoomsregion(region=rs73226617,chrom=chrom1,highlight= TRUE, extend=c(2,0),offsets=c(0,0),col="#FA80727F") #red
zoomsregion2(region=indel,chrom=chrom1,highlight=TRUE, extend=c(2,0),offsets=c(0,0.5))#,col="#ADD8E67F") #blue
zoomsregion2(region=rs58538713,chrom=chrom1,genome=NULL, extend=c(2,0),offsets=c(0,0),highlight = TRUE)#,col="#FFDAB97F")
zoomsregion2(region=rs73239724,chrom=chrom1,genome=NULL, extend=c(2,0),offsets=c(0,0),highlight = TRUE)
zoomsregion2(region=rs73226608,chrom=chrom1,genome=NULL, extend=c(2,0),offsets=c(0,0),highlight = TRUE)
zoomsregion2(region=rs11710569,chrom=chrom1,genome=NULL, extend=c(2,0),offsets=c(0,0),highlight = TRUE)
zoomsregion2(region=rs11714221,chrom=chrom1,genome=NULL, extend=c(2,0),offsets=c(0,0),highlight = TRUE)
zoomsregion2(region=rs189348793,chrom=chrom1,genome=NULL, extend=c(2,0),offsets=c(0,0),highlight = TRUE)
zoomsregion2(region=rs73226619,chrom=chrom1,genome=NULL, extend=c(2,0),offsets=c(0,0),highlight = TRUE)
zoomsregion2(region=rs112304167,chrom=chrom1,genome=NULL, extend=c(2,0),offsets=c(0,0),highlight = TRUE)
zoomsregion2(region=rs146526750,chrom=chrom1,genome=NULL, extend=c(2,0),offsets=c(0,0),highlight = TRUE)

labelplot("(d). Cell type (monocyte) specific activity",letteradj=-.085,lettercex=1.5)

axis(side=2,las=2,tcl=.5,at=c(0,0.5,1),labels = c("0","0.5","1"),cex.axis=2)
par(mar=c(2,5,3, 4) )
mtext("IMPACT score",side=2,line=3,cex=1.5,font=2)


par(mar=c(1,5,4, 4) )
labelgenome( chrom1, chromstart2,chromend2,n=5,scale="Kb",scalecex =1.5,chromcex = 1.5,cex.axis=2)

dev.off()

