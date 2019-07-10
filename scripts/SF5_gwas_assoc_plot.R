## Manhattan and qqplot ###
source("scripts/qqplot.R")
source("scripts/manhattan.R")

#### rare association study
load("data/rare-plot.RData")
png("figures/SF5d_rare-qqplot.png",width=800,height=800)
qqunif.plot(out$P)
dev.off()

png("figures/SF5c_rare-manhattan.png",width=1000,height=650)
manhattan(out)
dev.off()

#### common association study stats have been deposit to GWAS Catalog, 
#### plotting script is similar as above

x <- read.table("data/all.assoc.txt.gz",h=F,stringsAsFactors = F)
names(x)<-c("CHR","BP","P")

png("figures/SF18a_withancestry-qqplot.png",width=800,height=800)
qqunif.plot(x$P)
dev.off()

png("figures/SF18b_withancestry-manhattan.png",width=1000,height=650)
manhattan(x)
dev.off()