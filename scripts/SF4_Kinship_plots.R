######REAP Estimate ######

reap<-scan("data/reap_kinship.txt.gz",as.numeric())


#######PLINK Estimate#####

plink<-scan("data/plink_kinship.txt.gz",as.numeric())

#####PLINK vs REAP ###
png("figures/SF4_plink-vs-reap-kinship.png",width=1000,height=800)
plot(plink,reap,pch=20,xlab="PLINK kinship estimates",ylab="REAP kinship estimates",xlim=c(0,0.5),ylim=c(0,0.5))

lines(c(0.5,0.5),c(0,0.5),col="red",lty=3,lwd=3)
lines(c(0.5,0),c(0.5,0.5),col="red",lty=3,lwd=3)

lines(c(0.25,0.25),c(0,0.25),col="orange",lty=3,lwd=3)
lines(c(0.25,0),c(0.25,0.25),col="orange",lty=3,lwd=3)

lines(c(0.125,0.125),c(0,0.125),col="blue",lty=3,lwd=3)
lines(c(0.125,0),c(0.125,0.125),col="blue",lty=3,lwd=3)

lines(c(0.0625,0.0625),c(0,0.0625),col="green",lty=3,lwd=3)
lines(c(0.0625,0),c(0.0625,0.0625),col="green",lty=3,lwd=3)

lines(c(0.03125,0.03125),c(0,0.03125),col="yellow",lty=3,lwd=3)
lines(c(0.03125,0),c(0.03125,0.03125),col="yellow",lty=3,lwd=3)

abline(h = 0, col = "black", lty = 1)
abline(v = 0, col = "black", lty = 1)
dev.off()
