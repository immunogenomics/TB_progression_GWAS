# Some of the data is too large to load, script only for reproducing SF18 

#source("https://bioconductor.org/biocLite.R")
#biocLite("Gviz","EnsDb.Hsapiens.v75")

library(Gviz)
library(GenomicRanges)
library(EnsDb.Hsapiens.v75)

#Load data : class = GRanges
mono <- read.table("~/Downloads//GSE74912_ATACseq_All_Counts.txt.gz", sep = "\t", header = T, stringsAsFactors = FALSE)
monocyte_names <- c("X6792.7A","X6792.7B","Donor7256.7A.141106","Donor7256.7B.141106")
dataset_names <- names(mono)
mono_dat <- mono[,c(1:3,match(monocyte_names,dataset_names))]

locus_start <- 141383625
locus_end <- 141406933

window <- 5000

#Annotation track, title ="Monocytes"
monoDat <- mono_dat[mono_dat$Chr == "chr3" & mono_dat$Start >= locus_start - window ,]
mono <- monoDat[monoDat$End <= locus_end+window,]
atrack <- AnnotationTrack(start=mono$Start,end=mono$End, data = as.vector(mono$Donor7256.6A.141106),name = "Primary monocytes",genome="hg19", chromosome = 'chr3')
#plotTracks(atrack)

## genomic coordinates
gtrack <- GenomeAxisTrack()
plotTracks(list(gtrack,atrack))

#Ideogram track
itrack <- IdeogramTrack(genome = 'hg19', chromosome = 'chr3')
plotTracks(list(itrack, gtrack, atrack))

#add gene names
## Retrieving a Gviz compatible GRanges object with all genes
## encoded on chromosome 3.
edb <- EnsDb.Hsapiens.v75
seqlevelsStyle(edb) <- "UCSC"

gr <- getGeneRegionTrackForGviz(edb, chromosome = "chr3",
                                start = locus_start-window, end = locus_end+window)
genes <- unique(gr$symbol)

transcripts.plots <- c("ENSE00002023189","ENST00000286364","ENST00000480908","ENST00000474979","ENST00000286371")
grtrack <- GeneRegionTrack(gr[gr$transcript %in% transcripts.plots,], name = 'Genes',
                           background.title = "white")


#Data Tracks
# THP1 H3K4me (data too big to load)
thp1_h3k4me1 <- DataTrack(range = "~/Downloads/GSM3514950_NM-L5_S5_L001_R1_001_filt_sorted.bw",
                          genome = 'hg19',
                          type="h",
                          chromosome = 'chr3',
                          name = "THP-1 ChIP-Seq \n (H3K4me1)",
                          col = "orange",ylim=c(0,50),
                          background.title = "orange")


#peaks from ATAC-Seq data
THP1_ATAC_peak <- AnnotationTrack(range = bed_to_granges('~/Downloads/GSM3731848_R1_StATAC_THP1_PMA_Rep2A_24hr_peaks.broadPeak.gz'),
                                  genome = 'hg19',
                                  name = '',
                                  chromosome = 'chr3',
                                  fill='red',ylim=c(0,50))

dat <- read.table("~/Downloads/GSM3731848_R1_StATAC_THP1_PMA_Rep2A_24hr_peaks.broadPeak.gz",h=F,stringsAsFactors = F)
thp1Dat <- dat[dat$V1 == "chr3" & dat$V2 >= locus_start - window ,]
thp1 <- thp1Dat[thp1Dat$V3 <= locus_end + window,]
thp1ATACtrack <- DataTrack(start=thp1$V2,end=thp1$V3, data = (thp1$V9),
                           chromosome = 'chr3', genome = 'hg19',
                           name = "THP1 ATAC-seq \n (peaks)",
                           type='histogram',fill="orange",
                           ylim=c(0,80),
                           background.title = "orange")
#plotTracks(thp1ATACtrack)

#Pirmary Monocytes
# ENCODE ChIP-Seq CD14+ monocyte
encodetrack <- DataTrack(range = "~/Downloads/GSM1003535_hg19_wgEncodeBroadHistoneMonocd14ro1746H3k04me1Sig.bigWig",
                         genome = 'hg19',
                         type="h",
                         chromosome = 'chr3',
                         name = "Monocyte ChIP-Seq \n (H3K4me1)",
                         col = "darkblue",
                         #background.panel = "#FFFEDB",
                         background.title = "darkblue")


# ATAC-seq
monoATAC <- DataTrack(range = "~/Downloads/GSM2325687_ATAC_Mono_8444.bw",
                      genome = 'hg19',
                      type = 'h',
                      chromosome = 'chr3',
                      name = "Monocyte ATAC-seq",
                      col = "steelblue")


monoATACpeak <- DataTrack(start=mono_dat$Start,end=mono_dat$End, 
                          data = (mono_dat$X6792.5B),
                          chromosome = 'chr3', genome = 'hg19',
                          name = "Monocyte ATAC-seq \n (peaks)",
                          type='histogram',fill="darkblue",
                          #background.panel = "#FFFEDB",
                          background.title = "darkblue")
#Plot
# variant position
variant_positions <- c(141383625, 141388147, 141389447, 141390125, 141390219, 141391318, 141400653, 141401146, 141403952, 141404369, 141406933)
end_positions <- c(141383625, 141388147, 141389447, 141390125, 141390219, 141391318, 141400653 + 100, 141401146, 141403952, 141404369, 141406933)

#highlight 
ht<- HighlightTrack(trackList = list(grtrack,encodetrack,thp1_h3k4me1,monoATACpeak,thp1ATACtrack), 
                    start = variant_positions, end = end_positions ,chromosome = 'chr3')
plotTracks(list(itrack, gtrack,ht),from=locus_start-window,to=locus_end+window,
           transcriptAnnotation = "symbol",fontsize = 16)

