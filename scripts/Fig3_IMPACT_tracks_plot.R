library(data.table)
library(readr)
library(scales)

TFs <- c("Tbet","Gata3","Stat3","Foxp3","Stat5","IRF5","IRF1_Mac","IRF1","CEBPB","PAX5","HNF4A","PolII","TCF7L2","RXRA","Rest","CUX1","adaptive","innate","Smad2")
fancyTFs <- c("Th1 (T-BET)","Th2 (GATA3)","Th17 (STAT3)", "Treg (FOXP3)","Treg (STAT5)","Macrophage (IRF5)","Macrophage (IRF1)","Monocyte (IRF1)","Monocyte (CEBPB)",
              "B cell (PAX5)","Liver (HNF4A)","Blood (PolII)","Pancreas (TCF7L2)","Brain (RXRA)","Fetal Brain (REST)", "K562 (CUX1)","adaptive","innate","Cardio (SMAD2)")


TFs <- c("Tbet","Gata3","Stat3","Foxp3","Stat5","IRF5","IRF1_Mac","IRF1","CEBPB","PAX5","HNF4A","PolII","TCF7L2","RXRA","Rest")
fancyTFs <- c("Th1 (T-BET)","Th2 (GATA3)","Th17 (STAT3)", "Treg (FOXP3)","Treg (STAT5)","Macrophage (IRF5)","Macrophage (IRF1)","Monocyte (IRF1)","Monocyte (CEBPB)",
              "B cell (PAX5)","Liver (HNF4A)","Blood (PolII)","Pancreas (TCF7L2)","Brain (RXRA)","Fetal Brain (REST)")

start <- 141383625


load("data/IMPACT_enh_regions.Rdata")

variant_positions <- c(141383625, 141388147, 141389447, 141390125, 141390219, 141391318, 141400653, 141401146, 141403952, 141404369, 141406933)
variant_names <- c("rs73239724", "rs73226608", "rs58538713", "rs11710569", "rs11714221", "rs189348793", "rs73226617", "rs148722713", "rs73226619", "rs112304167", "rs146526750")

label_positions <- (variant_positions - start)/3
label_positions <- c(0.000 ,1507.333, 1940.667, 2166.667, 
                     2198.00 + 100,2564.333 ,
                     5676.000, 5840.333 ,6775.667, 6914.667+100, 7769.333)

pdf("figures/Fig3b.pdf")
par(mfrow = c(length(TFs),1))
par(cex = 1)
par(mar = c(0, 0, 0, 0), oma = c(4, 4, 0.5, 0.5))
par(tcl = -0.25)
par(mgp = c(2, 0.6, 0))

TFcolors <- hue_pal()(length(TFs))
for (i in 1:length(TFs)){
  dat <- get(paste0(TFs[i],"_dat"))
  plot(c(0,dat$V1,0), axes = F, type = "n", ylim = c(0,1.1), ylab = "")
  polygon(c(0,dat$V1,0), col = TFcolors[i], border = NA)
  #mtext(fancyTFs[i], side = 3, line = -1.5, adj = 1, cex = 1)
  legend("topleft",fancyTFs[i], bty = "n", cex = 0.8)
  axis(2, col = "black", col.axis = "black", at = c(0.1,1), cex.axis = 0.5, labels = c(0,1))
  abline(v = (variant_positions - start)/3, lty = 4, lwd = 1) #segments where variants are
  box(col = "black")
}
axis(1,at=label_positions[1:6],labels=variant_names[1:6],tck=.01,cex.axis=0.7, srt=60,las=2, col.ticks = "grey")
axis(1,at=label_positions[7],labels=variant_names[7],tck=.01,font=2,cex.axis=0.7, srt=60,las=2, col.ticks = "grey")
axis(1,at=label_positions[8:11],labels=variant_names[8:11],tck=.01,font=1,cex.axis=0.7, srt=60,las=2, col.ticks = "grey")
mtext("IMPACT score (p(cell-state regulatory element))",side=2,outer=TRUE,padj = -2)
dev.off()


