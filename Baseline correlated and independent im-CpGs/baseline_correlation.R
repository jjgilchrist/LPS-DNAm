rm(list = ls(all=TRUE))

library(ggplot2)
library(ggrepel)
library(RColorBrewer)
library(XGR)
library(GenomicRanges)

RData.location <- "http://galahad.well.ox.ac.uk/bigdata.backup"

cols <- brewer.pal(9, "Set1")


#read in delta methylation values for 7350 CpGs differentially methylation by LPS in monocytes
delta <- read.table("~/Documents/fairfax_lab/cd14_rnaseq/mQTL/final/code_for_github/delta_methylation.txt",header = T, row.names = 1)

#read in baseline methylation values for 7350 CpGs differentially methylation by LPS in monocytes
baseline<- read.table("~/Documents/fairfax_lab/cd14_rnaseq/mQTL/final/code_for_github/baseline_methylation.txt",header = T, row.names = 1)

stats <- read.table("~/Documents/fairfax_lab/cd14_rnaseq/mQTL/final/code_for_github/LPS_diff_CpG_stats.txt", header = T, row.names = 1)


#correlate baseline methylation and LPS-induced change in methylation at each CpG

cor.p <- c()
cor.r <- c()
for (i in c(1:7350)){
cor.p[i] <- cor.test(delta[,i], baseline[,i], method = 'pearson')$p.value
cor.r[i] <- cor.test(delta[,i], baseline[,i], method = 'pearson')$estimate
}


#add in logFC estimates for LPS stimulation at each CpG
for_plot <- cbind(stats, cor.p, cor.r)

#label baseline correlated CpGs
for_plot$is.cor <- 0
for_plot$is.cor[which(for_plot$cor.p<0.05)] <- 1

#label CpGs with LPS-induced gain in methylation
for_plot$is.meth <- 0
for_plot$is.meth[which(for_plot$logFC>0)] <- 1

table(for_plot$is.meth, for_plot$is.cor)

#scatter plot highlighting baselin correlated and uncorrelated CpGs

plot <- ggplot(for_plot, aes(x=cor.r, y=logFC)) +
  geom_point(alpha = 0.72, size = 5) +
  aes(fill = factor(is.cor), col=factor(is.cor)) + scale_fill_manual(values = cols[c(1,2)]) + scale_colour_manual(values = cols[c(1,2)]) +
  ylab("LPS-induced change in methylation: log(fold change)") +
  xlab("Pearson's r: baseline methylation vs LPS-induced methylation change") +
  #scale_x_continuous(breaks=c(0,0.5,1.0), limits = c(0,1)) +
  theme_bw() +
  theme(legend.position = "none", axis.text=element_text(size=20),
        axis.title=element_text(size=30),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),axis.line = element_line(colour = "grey"))

plot

ggplot2::ggsave("scatter_plot.baseline_correlation.jpg",
width = 17,
height = 14,
dpi = 300
)

#explore differences in baseline methylation and LPS effect by baseline correlation
for_plot2 <- cbind(stats, for_plot)
for_plot3 <- for_plot2[,c("AveExpr", "t", "is.cor", "is.meth")]

for_plot3$corr <- "correlated"
for_plot3$corr[for_plot3$is.cor==0] <- "independent"

for_plot3$corr <- factor(for_plot3$corr)

for_plot3$response <- "demethylated"
for_plot3$response[for_plot3$is.meth==1] <- "methylated"

for_plot3$response <- factor(for_plot3$response)

# Comparing distribution of baseline methylation of demethylated and methylated correlated and non correlated CpGs
#demeth
ks.test(x = subset(for_plot3, is.cor==0 & is.meth==0)$AveExpr, y = subset(for_plot3, is.cor==1 & is.meth==0)$AveExpr)
#meth
ks.test(x = subset(for_plot3, is.cor==0 & is.meth==1)$AveExpr, y = subset(for_plot3, is.cor==1 & is.meth==1)$AveExpr)

# Comparing distribution of effect size of demethylated and methylated correlated and non correlated CpGs
#demeth
ks.test(x = subset(for_plot3, is.cor==0 & is.meth==0)$t, y = subset(for_plot3, is.cor==1 & is.meth==0)$t)
#meth
ks.test(x = subset(for_plot3, is.cor==0 & is.meth==1)$t, y = subset(for_plot3, is.cor==1 & is.meth==1)$t)

#direction of effect for response
#demeth
wilcox.test(x = subset(for_plot3, is.cor==0 & is.meth==0)$t, y = subset(for_plot3, is.cor==1 & is.meth==0)$t)
#meth
wilcox.test(x = subset(for_plot3, is.cor==0 & is.meth==1)$t, y = subset(for_plot3, is.cor==1 & is.meth==1)$t)

#for demeth - BI CpGs have greater effect size of LPS change (p<2.2e-16)
#for meth - BI CpGs have lesser effect size of LPS change (p=0.002)

# plot baseline methylation densities for baseline correlated and independent CpGs

baseline.dens <- ggplot(for_plot3,aes(AveExpr))+
  geom_density(aes(colour=factor(corr)), size=3)+
  scale_colour_manual(values=c("steelblue4","firebrick4"))+
  facet_wrap(~response)+
  theme_bw()+
  theme(axis.text.y = element_text(size = 20),
        axis.text.x = element_text(size = 20),
        axis.title = element_text(size = 25),
        strip.text.x = element_text(size = 30, face = "bold"),
        legend.title = element_blank(),
        legend.text = element_text(size = 25),
        legend.position="top")+
  scale_x_continuous("Average beta", breaks = c(0,0.25,0.5,0.75), labels = c("0","0.25","0.5","0.75"))

# plot LPS effect size for baseline coorrelated and independent CpGs

resp.dens <- ggplot(for_plot3,aes(t))+
  geom_density(aes(colour=factor(corr)), size=3)+
  scale_colour_manual(values=c("steelblue4","firebrick4"))+
  facet_wrap(~response, scales="free")+
  theme_bw()+
  theme(axis.text.y = element_text(size = 20),
        axis.text.x = element_text(size = 20),
        axis.title = element_text(size = 25),
        strip.text.x = element_text(size = 30, face = "bold"),
        legend.title = element_blank(),
        legend.text = element_text(size = 25),
        legend.position="top")+
  scale_x_continuous("t statistic")


dens_plot_out <- (baseline.dens|resp.dens) 

ggplot2::ggsave("dens_plots.jpg",
  width = 28,
  height = 10,
  dpi = 300
)

#Assess enrichment of baseline correlated and uncorrelated LPS-modulated CpGs among genomic features

#read in background CpGs
b.ground <- read.table("~/Documents/fairfax_lab/cd14_rnaseq/mQTL/final/code_for_github/background_probes.txt", header = T, row.names = 1)

#create Granges object for background CpGs tested
pGR <- GRanges(
  seqnames=Rle(b.ground$chr),
  ranges = IRanges(start=as.numeric(b.ground$start), end=as.numeric(b.ground$end), names=rownames(b.ground)),
  strand = Rle(rep('*',nrow(b.ground)))
)

#create Granges object for demethylated LPS-modulated CpGs
ind <- match(rownames(subset(for_plot3, is.meth==0)), names(pGR))
demethGR <- pGR[ind[!is.na(ind)]]

#create Granges object for hypermethylated LPS-modulated CpGs
ind <- match(rownames(subset(for_plot3, is.meth==1)), names(pGR))
methGR <- pGR[ind[!is.na(ind)]]

#create Granges object for all LPS-modulated CpGs
ind <- match(rownames(for_plot3), names(pGR))
sigGR <- pGR[ind[!is.na(ind)]]

#create Granges object for baseline correlated LPS-modulated CpGs
ind <- match(rownames(subset(for_plot3, is.cor==1)), names(pGR))
corrGR <- pGR[ind[!is.na(ind)]]

#create Granges object for baseline independent LPS-modulated CpGs
ind <- match(rownames(subset(for_plot3, is.cor==0)), names(pGR))
indGR <- pGR[ind[!is.na(ind)]]

#create Granges object for baseline correlated demethylated LPS-modulated CpGs
ind <- match(rownames(subset(for_plot3, is.cor==1 & is.meth==0)), names(pGR))
demeth.corr.GR <- pGR[ind[!is.na(ind)]]

#create Granges object for baseline independent demethylated LPS-modulated CpGs
ind <- match(rownames(subset(for_plot3, is.cor==0 & is.meth==0)), names(pGR))
demeth.ind.GR <- pGR[ind[!is.na(ind)]]

#create Granges object for baseline correlated hypermethylated LPS-modulated CpGs
ind <- match(rownames(subset(for_plot3, is.cor==1 & is.meth==1)), names(pGR))
meth.corr.GR <- pGR[ind[!is.na(ind)]]

#create Granges object for baseline independent hypermethylated LPS-modulated CpGs
ind <- match(rownames(subset(for_plot3, is.cor==0 & is.meth==1)), names(pGR))
meth.ind.GR <- pGR[ind[!is.na(ind)]]

#read in intervals for transcription factor binding sites in K562 cells
all_gr <- xRDataLoader(RData.customised="Uniform_TFBS", RData.location=RData.location)
ind <- grep('K562', names(all_gr))
anno_gr <- all_gr[ind]

# Remove Ifn treated samples
ind2 <- grep('Ifn', names(anno_gr))
anno_gr <- anno_gr[-ind2]

### use XGR to compute K562 TFBS enrichment in:
#baseline correlated, LPs-demethylated CpGs
eTerm_tfbs_CorDemeth <- xGRviaGenomicAnno(data.file=demeth.corr.GR, annotation.file=anno_gr, background.file=pGR, format.file="GRanges", resolution="hybrid", RData.location=RData.location, p.tail = "two-tails")
xEnrichViewer(eTerm_tfbs_CorDemeth, 50)

#baseline independent, LPs-demethylated CpGs
eTerm_tfbs_IndDemeth <- xGRviaGenomicAnno(data.file=demeth.ind.GR, annotation.file=anno_gr, background.file=pGR, format.file="GRanges", resolution="hybrid", RData.location=RData.location, p.tail = "two-tails")
xEnrichViewer(eTerm_tfbs_IndDemeth, 50)

#baseline correlated, LPs-hypermethylated CpGs
eTerm_tfbs_CorMeth <- xGRviaGenomicAnno(data.file=meth.corr.GR, annotation.file=anno_gr, background.file=pGR, format.file="GRanges", resolution="hybrid", RData.location=RData.location, p.tail = "two-tails")
xEnrichViewer(eTerm_tfbs_CorMeth, 50)

#baseline independent, LPs-hypermethylated CpGs
eTerm_tfbs_IndMeth <- xGRviaGenomicAnno(data.file=meth.ind.GR, annotation.file=anno_gr, background.file=pGR, format.file="GRanges", resolution="hybrid", RData.location=RData.location, p.tail = "two-tails")
xEnrichViewer(eTerm_tfbs_IndMeth, 50)

#plot out enrichment
tfbs.plot <- rbind(xEnrichViewer(eTerm_tfbs_CorMeth, 50), xEnrichViewer(eTerm_tfbs_CorDemeth, 50), xEnrichViewer(eTerm_tfbs_IndDemeth, 50), xEnrichViewer(eTerm_tfbs_IndMeth, 50))

tfbs.plot$facet <- c(rep("BC_meth",50),
                     rep("BC_demeth",50),
                     rep("BI_demeth",50),
                     rep("BI_meth",50))

tfbs.plot$facet <- factor(tfbs.plot$facet, levels = c("BC_meth", "BI_meth", "BC_demeth", "BI_demeth"))

tfbs.plot$is.annot <- 0
tfbs.plot$is.annot[c(1,2,3, 51,52,53,101,102,103)] <- 1

top_labelled<-tfbs.plot[which(tfbs.plot$is.annot==1),]

tfbs_enrich<-ggplot(tfbs.plot,aes(zscore,fc))+
  geom_point(data = subset(tfbs.plot, adjp>0.05),pch=21,colour="grey30", fill="grey30", size = 1)+
  geom_point(data = subset(tfbs.plot, adjp<0.05 & fc!=0), aes(fill=-log10(pvalue),size=nOverlap),pch=21,colour="grey30")+
  scale_fill_distiller(palette = "Spectral") +
  theme_bw()+ 
  theme(axis.title=element_text(size = 15),
        axis.text=element_text(size = 15),
        legend.title=element_text(size=10),
        legend.text=element_text(size = 10),
        strip.text = element_blank(),
        plot.title = element_text(size = 25, face = "bold")) +
  scale_size(range=c(2,15),"nOverlap")+
  scale_y_continuous("fold change")+
  facet_wrap(~facet, ncol=2) +
  #guides(size=guide_legend(),fill=guide_legend())+
  geom_text_repel(data = top_labelled, 
                  mapping = aes(label = name),
                  size = 7,
                  fontface = 'italic',
                  color = 'black')

tfbs_enrich

ggplot2::ggsave("tfbs_enrich_corr_indep.jpg",
                 height = 7,
                 width = 7,
                 dpi = 300
)

#Point towards intervals for ChromHMM chromatin states (15 state model) 
GR.annotation <- "EpigenomeAtlas_15Segments_E029"

### use XGR to compute ChromHMM chromatin state enrichment in:
#baseline correlated, LPS-demethylated CpGs
eTerm_chrom_CorDemeth <- xGRviaGenomicAnno(data.file=demeth.corr.GR, GR.annotation=GR.annotation, background.file=pGR, format.file="GRanges", resolution="regions", RData.location=RData.location, p.tail = "two-tails")
xEnrichViewer(eTerm_chrom_CorDemeth, 50)

#baseline independent, LPS-demethylated CpGs
eTerm_chrom_IndDemeth <- xGRviaGenomicAnno(data.file=demeth.ind.GR, GR.annotation=GR.annotation, background.file=pGR, format.file="GRanges", resolution="regions", RData.location=RData.location, p.tail = "two-tails")
xEnrichViewer(eTerm_chrom_IndDemeth, 50)

#baseline correlated, LPS-hypermethylated CpGs
eTerm_chrom_CorMeth <- xGRviaGenomicAnno(data.file=meth.corr.GR, GR.annotation=GR.annotation, background.file=pGR, format.file="GRanges", resolution="regions", RData.location=RData.location, p.tail = "two-tails")
xEnrichViewer(eTerm_chrom_CorMeth, 50)

#baseline independent, LPS-hypermethylated CpGs
eTerm_chrom_IndMeth <- xGRviaGenomicAnno(data.file=meth.ind.GR, GR.annotation=GR.annotation, background.file=pGR, format.file="GRanges", resolution="regions", RData.location=RData.location, p.tail = "two-tails")
xEnrichViewer(eTerm_chrom_IndMeth, 50)

#plot out chromatin enrichment
BC_demeth_chrom<-data.frame(xEnrichViewer(eTerm_chrom_CorDemeth, 16))
BI_demeth_chrom<-data.frame(xEnrichViewer(eTerm_chrom_IndDemeth, 16))
BC_meth_chrom<-data.frame(xEnrichViewer(eTerm_chrom_CorMeth, 16))
BI_meth_chrom<-data.frame(xEnrichViewer(eTerm_chrom_IndMeth, 16))

BC_demeth_chromDf<-data.frame(BC_demeth_chrom)
BI_demeth_chromDf<-data.frame(BI_demeth_chrom)
BC_meth_chromDf<-data.frame(BC_meth_chrom)
BI_meth_chromDf<-data.frame(BI_meth_chrom)

a<-unlist(strsplit(BC_demeth_chromDf[,1],"\\("))
names2<-a[grep("\\)",a)]
b<-unlist(strsplit(names2,"\\)"))
BC_demeth_chromDf$names<-b
BC_demeth_chromDf$names<-factor(BC_demeth_chromDf$names, levels=c("Active TSS","Flanking Active TSS","Transcr. at gene 5' and 3'","Strong transcription","Weak transcription","Genic enhancers","Enhancers","ZNF genes & repeats","Heterochromatin","Bivalent/Poised TSS","Flanking Bivalent TSS/Enh","Bivalent Enhancer","Repressed PolyComb","Weak Repressed PolyComb","Quiescent/Low"))

a<-unlist(strsplit(BI_demeth_chromDf[,1],"\\("))
names2<-a[grep("\\)",a)]
b<-unlist(strsplit(names2,"\\)"))
BI_demeth_chromDf$names<-b
BI_demeth_chromDf$names<-factor(BI_demeth_chromDf$names, levels=c("Active TSS","Flanking Active TSS","Transcr. at gene 5' and 3'","Strong transcription","Weak transcription","Genic enhancers","Enhancers","ZNF genes & repeats","Heterochromatin","Bivalent/Poised TSS","Flanking Bivalent TSS/Enh","Bivalent Enhancer","Repressed PolyComb","Weak Repressed PolyComb","Quiescent/Low"))

a<-unlist(strsplit(BC_meth_chromDf[,1],"\\("))
names2<-a[grep("\\)",a)]
b<-unlist(strsplit(names2,"\\)"))
BC_meth_chromDf$names<-b
BC_meth_chromDf$names<-factor(BC_meth_chromDf$names, levels=c("Active TSS","Flanking Active TSS","Transcr. at gene 5' and 3'","Strong transcription","Weak transcription","Genic enhancers","Enhancers","ZNF genes & repeats","Heterochromatin","Bivalent/Poised TSS","Flanking Bivalent TSS/Enh","Bivalent Enhancer","Repressed PolyComb","Weak Repressed PolyComb","Quiescent/Low"))

a<-unlist(strsplit(BI_meth_chromDf[,1],"\\("))
names2<-a[grep("\\)",a)]
b<-unlist(strsplit(names2,"\\)"))
BI_meth_chromDf$names<-b
BI_meth_chromDf$names<-factor(BI_meth_chromDf$names, levels=c("Active TSS","Flanking Active TSS","Transcr. at gene 5' and 3'","Strong transcription","Weak transcription","Genic enhancers","Enhancers","ZNF genes & repeats","Heterochromatin","Bivalent/Poised TSS","Flanking Bivalent TSS/Enh","Bivalent Enhancer","Repressed PolyComb","Weak Repressed PolyComb","Quiescent/Low"))

BC_demeth_chromDf$sig <- 0
BC_demeth_chromDf$sig[which(BC_demeth_chromDf$adjp<0.05)] <- 1
BI_demeth_chromDf$sig <- 0
BI_demeth_chromDf$sig[which(BI_demeth_chromDf$adjp<0.05)] <- 1
BC_meth_chromDf$sig <- 0
BC_meth_chromDf$sig[which(BC_meth_chromDf$adjp<0.05)] <- 1
BI_meth_chromDf$sig <- 0
BI_meth_chromDf$sig[which(BI_meth_chromDf$adjp<0.05)] <- 1


BC_demeth_chrom <- ggplot(BC_demeth_chromDf,aes(names,zscore, fill=factor(sig)))+
  geom_bar(stat="identity",position="dodge",colour="grey30")+
  scale_fill_manual(values=c("grey70", "steelblue2")) +
  ylim(-100, 100) +
  ylab("z score") +
  theme_bw()+
  theme(axis.title.y=element_text(size = 30), axis.text.y=element_text(size = 30), axis.text.x=element_text(angle = 90,vjust=c(0.5), size=30), axis.title.x = element_blank(), legend.position="none")

BI_demeth_chrom <- ggplot(BI_demeth_chromDf,aes(names,zscore, fill=factor(sig)))+
  geom_bar(stat="identity",position="dodge",colour="grey30")+
  scale_fill_manual(values=c("grey70", "firebrick2")) +
  ylim(-100, 100) +
  ylab("zscore") +
  theme_bw()+
  theme(axis.title.y=element_text(size = 30), axis.text.y=element_text(size = 30), axis.text.x=element_text(angle = 90,vjust=c(0.5), size=30), axis.title.x = element_blank(), legend.position="none")

BC_meth_chrom <- ggplot(BC_meth_chromDf,aes(names,zscore, fill=factor(sig)))+
  geom_bar(stat="identity",position="dodge",colour="grey30")+
  scale_fill_manual(values=c("grey70", "steelblue2")) +
  ylim(-20, 20) +
  ylab("z score") +
  theme_bw()+
  theme(axis.title.y=element_text(size = 30), axis.text.y=element_text(size = 30), axis.text.x=element_blank(), axis.title.x = element_blank(), legend.position="none")

BI_meth_chrom <- ggplot(BI_meth_chromDf,aes(names,zscore, fill=factor(sig)))+
  geom_bar(stat="identity",position="dodge",colour="grey30")+
  scale_fill_manual(values=c("grey70", "firebrick2")) +
  ylim(-20, 20) +
  ylab("zscore") +
  theme_bw()+
  theme(axis.title.y=element_text(size = 30), axis.text.y=element_text(size = 30), axis.text.x=element_blank(), axis.title.x = element_blank(), legend.position="none")


chrom.out <- (BC_meth_chrom|BI_meth_chrom)/(BC_demeth_chrom|BI_demeth_chrom)

  
ggplot2::ggsave("chrom_enrich_corr_indep.jpg",
                height = 14,
                width = 14,
                dpi = 300
)

