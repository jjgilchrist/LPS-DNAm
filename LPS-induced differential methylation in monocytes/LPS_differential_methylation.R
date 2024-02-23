rm(list = ls(all=TRUE))

library(ggplot2)
library(RColorBrewer)
library(dplyr)
library(ggrepel)
library(LDlinkR)
library(ggplotify)
library(patchwork)
library(ggbio)
library(data.table)
library(GenomicRanges)
library(LDlinkR)
library(biomaRt)
library(bacon)
library(limma)
library(qvalue)

library(XGR)
RData.location <- "http://galahad.well.ox.ac.uk/bigdata_dev"

cols <- brewer.pal(8,"Set2")
cols2 <- brewer.pal(8,"Dark2")
cols3 <- brewer.pal(8,"Paired")

#read in differential methylation output - as calculated with pairwise linear models on normalized data using the lmFit and eBayes function in limma
d1 <- read.table("diff_meth.out", header = T, row.names = 1)

#estimate bias and inflation with bacon from t statistics
bc <- bacon(d1$t)
inflation(bc)
bias(bc)

#calculate corrected p-values from adjusted t statistics (189 degrees of freedom)
corr.pvals <- pt(abs(tstat(bc)),189, lower.tail = FALSE)*2

#adjust new bacon-corrected p-values with q-value
corr.qvals <- qvalue(p = corr.pvals)$qvalues

#add corrected p/qvalues to dataframe
d.corr <- cbind(d1, corr.pvals, corr.qvals)

#plot manhattans

d2 <- d.corr %>% 
  
  group_by(chr) %>% 
  summarise(chr_len=max(bp)) %>% 
  
  mutate(tot=cumsum(as.numeric(chr_len))-chr_len) %>%
  dplyr::select(-chr_len) %>%
  
  left_join(d.corr, ., by=c("chr"="chr")) %>%
  arrange(chr, bp) %>%
  mutate( BPcum=bp+tot)


axisdf = d2 %>% group_by(chr) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )

total.diff.manh <- ggplot(d2, aes(x=BPcum, y=-log10(corr.pvals))) +
  
  geom_point(data = d2, aes(color=as.factor(chr)), size=3, alpha= 0.75) +
  scale_color_manual(values = rep(c("grey70","grey30"), 11 )) +
  
  scale_x_continuous( label = axisdf$chr[c(1:14,16,18,20,22)], breaks= axisdf$center[c(1:14,16,18,20,22)] ) +
  scale_y_continuous(expand = c(0, 0) ) + 
  ylim(NA,100) +
  xlab("chromosome") +
  
  theme_bw() +
  theme( 
    legend.position="none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.title=element_text(size=30),
    axis.text=element_text(size=25)
  )


ggplot2::ggsave(
  "total_manh.jpg",
  width = 14,
  height = 7,
  dpi = 100
)

meth.diff.manh <- ggplot(d2, aes(x=BPcum, y=-log10(corr.pvals))) +
  
  geom_point(data = subset(d2, logFC>0), aes(color=as.factor(chr)), size=3, alpha= 0.75) +
  scale_color_manual(values = rep(c("firebrick2","grey30"), 11 )) +
  
  scale_x_continuous( label = axisdf$chr[c(1:14,16,18,20,22)], breaks= axisdf$center[c(1:14,16,18,20,22)] ) +
  scale_y_continuous(expand = c(0, 0) ) + 
  ylim(NA,100) +
  xlab("chromosome") +
  
  theme_bw() +
  theme( 
    legend.position="none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.title=element_text(size=30),
    axis.text=element_text(size=25)
  )


ggplot2::ggsave(
  "meth_manh.jpg",
  width = 14,
  height = 7,
  dpi = 100
)

demeth.diff.manh <- ggplot(d2, aes(x=BPcum, y=-log10(corr.pvals))) +
  
  geom_point(data = subset(d2, logFC<0), aes(color=as.factor(chr)), size=3, alpha= 0.75) +
  scale_color_manual(values = rep(c("steelblue2","grey30"), 11 )) +
  
  scale_x_continuous( label = axisdf$chr[c(1:14,16,18,20,22)], breaks= axisdf$center[c(1:14,16,18,20,22)] ) +
  scale_y_continuous(expand = c(0, 0) ) + 
  ylim(NA,100) +
  xlab("chromosome") +
  
  theme_bw() +
  theme( 
    legend.position="none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.title=element_text(size=30),
    axis.text=element_text(size=25)
  )


ggplot2::ggsave(
  "demeth_manh.jpg",
  width = 14,
  height = 7,
  dpi = 100
)

#qq plot - comparing uncorrected and bacon corrected p-values

df <- data.frame(observed = -log10(sort(na.omit(d.corr$P.Value))),
                 corr = -log10(sort(na.omit(d.corr$corr.pvals))),
                   expected = -log10(ppoints(length(na.omit(d.corr$P.Value)))))

log10Pe <- expression(paste("Expected -log"[10], plain(P)))
log10Po <- expression(paste("Observed -log"[10], plain(P)))

meth.qq <- ggplot(df) +
  geom_abline(intercept = 0, slope = 1, size = 1, colour = cols3[6]) +
  geom_point(aes(expected, observed), size = 3, colour = cols3[2]) +
  geom_point(aes(expected, corr), size = 3, colour = cols3[8]) +
  xlab(log10Pe) +
  ylab(log10Po) +
  theme_bw() +
  theme( 
    legend.position="none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),axis.text=element_text(size=15),
    axis.title=element_text(size=20), plot.title=element_text(size=30, face = "bold", hjust = 0.5))

ggplot2::ggsave(
  "meth.qq.jpg",
  width = 7,
  height = 7,
  dpi = 100
)

#dataframe of high-confidence differentially methylated CpGs
sig.lps <- subset(d.corr, corr.qvals<0.05)

#find nearest genes to CpGs
UCSC_genes <- load("/Users/jamesgilchrist/Documents/fairfax_lab/cd14_rnaseq/mQTL/dis_ont/UCSC_genes.RData")
symbols <- mcols(UCSC_genes)$Symbol

d.corr$chr <- paste('chr',d.corr$chr, sep='')

pGR <- GRanges(
  seqnames=Rle(d.corr$chr),
  ranges = IRanges(start=as.numeric(d.corr$bp), end=(as.numeric(d.corr$bp)+50), names=rownames(d.corr)),
  strand = Rle(rep('*',nrow(d.corr)))
)

ind <- rownames(d.corr)
bGR <- pGR[ind[!is.na(ind)]]

## link UCSC_genes to probes
query <- pGR[ind[!is.na(ind)]]
subject <- UCSC_genes
maxgap <- -1L # 1MB
minoverlap <- 1L # 1b overlaps
q2r <- as.matrix(suppressWarnings(GenomicRanges::findOverlaps(query=query, subject=subject, maxgap=maxgap, minoverlap=minoverlap, type="any", select="all", ignore.strand=T)))
ind_query <- q2r[,1]
ind_subject <- q2r[,2]
probes_symbols <- data.frame(probes=names(query)[ind_query], symbols=symbols[ind_subject])


d.corr$gene <- NA
d.corr$gene <- probes_symbols$symbols[match(rownames(d.corr), probes_symbols$probes)]

head(d.corr)
str(probes_symbols)


#Analysis of chromatn usage
#First define GR objects for total, hyper and hypomethylated probes
total <- d.corr
meth <- subset(d.corr, logFC>0)
demeth <- subset(d.corr, logFC<0)

#used probes
ind <- match(rownames(total), names(pGR))
usedGR <- pGR[ind[!is.na(ind)]]
str(usedGR)

#hyperGR
meth.sig <- subset(meth, corr.qvals<0.05)
meth.nonsig <- subset(meth, corr.qvals>0.05)

ind <- match(rownames(meth.sig), names(pGR))
hyperGR <- pGR[ind[!is.na(ind)]]


#hypoGR
demeth.sig <- subset(demeth, corr.qvals<0.05)
demeth.nonsig <- subset(demeth, corr.qvals>0.05)

ind <- match(rownames(demeth.sig), names(pGR))
hypoGR <- pGR[ind[!is.na(ind)]]

#Assess enrichment with XGR against chromatin usage
GR.annotation <- "EpigenomeAtlas_15Segments_E029"

eTerm_segment_hyper <- xGRviaGenomicAnno(data.file=hyperGR, background.file=usedGR, format.file="GRanges", GR.annotation=GR.annotation, resolution="regions", RData.location=RData.location,p.tail="two-tails")
xEnrichViewer(eTerm_segment_hyper, 15)

eTerm_segment_hypo <- xGRviaGenomicAnno(data.file=hypoGR, background.file=usedGR, format.file="GRanges", GR.annotation=GR.annotation, resolution="regions", RData.location=RData.location,p.tail="two-tails")
xEnrichViewer(eTerm_segment_hypo, 15)


list_eTerm <- list(eTerm_segment_hyper, eTerm_segment_hypo)
names(list_eTerm) <- c('Hyper-methylated CpG', 'Hypo-methylated CpG')
bp_segment <- xEnrichCompare(list_eTerm, displayBy="fc", FDR.cutoff=1e-2, bar.label=T, signature=F)
bp_segment + scale_fill_manual(values=c("skyblue","pink1"))


hyperAll<-data.frame(xEnrichViewer(eTerm_segment_hyper, 16))
hypoAll<-data.frame(xEnrichViewer(eTerm_segment_hypo, 16))

hypoAllDf<-data.frame(hypoAll)
hyperAllDf<-data.frame(hyperAll)

a<-unlist(strsplit(hypoAllDf[,1],"\\("))
names2<-a[grep("\\)",a)]
b<-unlist(strsplit(names2,"\\)"))


hypoAllDf$names<-b

hypoAllDf$names<-factor(hypoAllDf$names, levels=c("Active TSS","Flanking Active TSS","Transcr. at gene 5' and 3'","Strong transcription","Weak transcription","Genic enhancers","Enhancers","ZNF genes & repeats","Heterochromatin","Bivalent/Poised TSS","Flanking Bivalent TSS/Enh","Bivalent Enhancer","Repressed PolyComb","Weak Repressed PolyComb","Quiescent/Low"))

a<-unlist(strsplit(hyperAllDf[,1],"\\("))
names2<-a[grep("\\)",a)]
b<-unlist(strsplit(names2,"\\)"))

hyperAllDf$names<-b

hyperAllDf$names<-factor(hyperAllDf$names, levels=c("Active TSS","Flanking Active TSS","Transcr. at gene 5' and 3'","Strong transcription","Weak transcription","Genic enhancers","Enhancers","ZNF genes & repeats","Heterochromatin","Bivalent/Poised TSS","Flanking Bivalent TSS/Enh","Bivalent Enhancer","Repressed PolyComb","Weak Repressed PolyComb","Quiescent/Low"))

hypoAllDf$sig <- 0
hypoAllDf$sig[which(hypoAllDf$adjp<0.05)] <- 1

hyperAllDf$sig <- 0
hyperAllDf$sig[which(hyperAllDf$adjp<0.05)] <- 1

#pdf("../figures/2.10.regulatedRegionsAll.pdf",width=10,height=8)
hypo_enrich <- ggplot(hypoAllDf,aes(names,zscore, fill=factor(sig)))+
  geom_bar(stat="identity",position="dodge",colour="grey30")+
  scale_fill_manual(values=c("grey70", "steelblue2")) +
  ylim(-100, 100) +
  ylab("z score") +
  #geom_hline(yintercept=c(-3,3), linetype="dashed", color = "grey40")+
  #geom_text(aes(label=ifelse(adjp<0.01,signif(-log(adjp,10),2),"n.s.")),nudge_y = -5,angle=45)+
  theme_bw()+
  theme(axis.title.y=element_text(size = 20), axis.text.y=element_text(size = 20), axis.text.x=element_blank(), axis.title.x = element_blank(), legend.position="none")

hyper_enrich <- ggplot(hyperAllDf,aes(names,zscore, fill=factor(sig)))+
  geom_bar(stat="identity",position="dodge",colour="grey30")+
  scale_fill_manual(values=c("grey70", "firebrick2")) +
  ylim(-100, 100) +
  ylab("zscore") +
  #geom_hline(yintercept=c(-3,3), linetype="dashed", color = "grey40")+
  #geom_text(aes(label=ifelse(adjp<0.01,signif(-log(adjp,10),2),"n.s.")),nudge_y = -5,angle=45)+
  theme_bw()+
  theme(axis.title.y=element_text(size = 20), axis.text.y=element_text(size = 20), axis.text.x=element_text(angle = 90,vjust=c(0.5), size=20), axis.title.x = element_blank(), legend.position="none")


chrom.out <- (hypo_enrich/hyper_enrich)

ggplot2::ggsave(
  "chrom_enrich.jpg",
  width = 7,
  height = 14,
  dpi = 100
)


#Next transcription factor usage - again in XGR

all_gr <- xRDataLoader(RData.customised="Uniform_TFBS", RData.location=RData.location)
ind <- grep('K562', names(all_gr))
anno_gr <- all_gr[ind]

ind2 <- grep('Ifn', names(anno_gr))
anno_gr <- anno_gr[-ind2]

eTerm_tfbs_hyper <- xGRviaGenomicAnno(data.file=hyperGR, annotation.file=anno_gr, background.file=usedGR, format.file="GRanges", resolution="hybrid", RData.location=RData.location,p.tail = "two-tails")
write.table(xEnrichViewer(eTerm_tfbs_hyper, 100), "~/Documents/fairfax_lab/cd14_rnaseq/mQTL/revisions/inflation/txn_f_usage.hyper.txt", quote = F)

eTerm_tfbs_hypo <- xGRviaGenomicAnno(data.file=hypoGR, annotation.file=anno_gr, background.file=usedGR, format.file="GRanges", resolution="hybrid", RData.location=RData.location,p.tail = "two-tails")
xEnrichViewer(eTerm_tfbs_hypo, 50)

colnames(eTerm_tfbs_hypo)[6] <- "P"
colnames(eTerm_tfbs_hyper)[6] <- "P"

top_labelled<-eTerm_tfbs_hypo[eTerm_tfbs_hypo$zscore>20,]

hypoMethTf<-ggplot(eTerm_tfbs_hypo[eTerm_tfbs_hypo$adjp<0.05,],aes(zscore,fc))+
  geom_point(aes(fill=-log10(P),size=-nOverlap),pch=21,colour="grey30")+
  scale_fill_distiller(palette = "Spectral") +
  theme_bw()+ 
  theme(axis.title=element_text(size = 25),
        axis.text=element_text(size = 25),
        legend.title=element_text(size=20),
        legend.text=element_text(size = 20),
        plot.title = element_text(size = 25, face = "bold")) +
  ggtitle("demethylated CpGs") +
  scale_size(range=c(2,15),"nOverlap")+
  scale_y_continuous("fold change")+
  #guides(size=guide_legend(),fill=guide_legend())+
  geom_text_repel(data = top_labelled, 
                  mapping = aes(label = gsub("Tfbs_Sydh_K562_","",gsub("Tfbs_Uchicago_K562_E","",gsub("Tfbs_Haib_K562_","",gsub("ab50322_Iggrab","",gsub("sc183_V0416101","",gsub("_Pcr1x","",gsub("_Iggrab","",gsub("ab50322","",gsub("sc150_V0422111","",name)))))))))),
                  size = 10,
                  fontface = 'italic',
                  color = 'grey40',
                  box.padding = unit(1.5, "lines"),
                  point.padding = unit(0.75, "lines"))

top_labelled<-eTerm_tfbs_hyper[eTerm_tfbs_hyper$zscore>3.5| eTerm_tfbs_hyper$zscore<(-3),]

hyperMethTf<-ggplot(eTerm_tfbs_hyper[eTerm_tfbs_hyper$adjp<0.05,],aes(zscore,fc))+
  geom_point(aes(fill=-log10(P),size=-nOverlap),pch=21,colour="grey30")+
  scale_fill_distiller(palette = "Spectral") +
  theme_bw()+ 
  theme(axis.title=element_text(size = 25),
        axis.text=element_text(size = 25),
        legend.title=element_text(size=20),
        legend.text=element_text(size = 20),
        plot.title = element_text(size = 25, face = "bold")) +
  scale_size(range=c(2,15),"nOverlap")+
  scale_y_continuous("fold change")+
  ggtitle("methylated CpGs") +
  #guides(size=guide_legend(),fill=guide_legend())+
  geom_text_repel(data = top_labelled, 
                  mapping = aes(label = gsub("301772a","",gsub("_V0416101","",gsub("Tfbs_Haib_K562_","",gsub("Tfbs_Uta_K562_","",gsub("Tfbs_Broad_K562_","",gsub("Tfbs_Sydh_K562_","",gsub("_Iggrab","",gsub("_Pcr1x","",name))))))))),
                  size = 10,
                  fontface = 'italic',
                  color = 'grey40',
                  box.padding = unit(1.5, "lines"),
                  point.padding = unit(0.75, "lines"))

txn_f.plot <- (hypoMethTf/hyperMethTf)

ggplot2::ggsave(
  "txn_factor_usage.jpg",
  width = 14,
  height = 14,
  dpi = 100
)

#enrichment in promoter capture HiC

ig <- xRDataLoader(RData.customised="ig.PCHiC.Monocytes", RData.location=RData.location)
df_edges <- igraph::get.data.frame(ig, what="edges")
gr_prey <- xGR(data=df_edges$to, format="chr:start-end", RData.location=RData.location)
annotation.file <- list(hyper=hyperGR, hypo=hypoGR)
eTerm_prey <- xGRviaGenomicAnno(data.file=gr_prey, annotation.file=annotation.file, background.file=usedGR, format.file="GRanges", resolution="region", RData.location=RData.location)
xEnrichViewer(eTerm_prey, 2)

gr_bait <- xGR(data=df_edges$from, format="chr:start-end", RData.location=RData.location)
eTerm_bait <- xGRviaGenomicAnno(data.file=gr_bait, annotation.file=annotation.file, background.file=usedGR, format.file="GRanges", resolution="region", RData.location=RData.location)
xEnrichViewer(eTerm_bait, 2)

bp_hic <- xEnrichCompare(list_eTerm, displayBy="zscore", FDR.cutoff=1, bar.label=F, signature=F)



query<-list("ig.PCHiC.Monocytes","ig.PCHiC.Macrophages_M0","ig.PCHiC.Macrophages_M1","ig.PCHiC.Macrophages_M2","ig.PCHiC.Neutrophils","ig.PCHiC.Megakaryocytes","ig.PCHiC.Endothelial_precursors","ig.PCHiC.Erythroblasts","ig.PCHiC.Fetal_thymus","ig.PCHiC.Naive_CD4_T_cells","ig.PCHiC.Total_CD4_T_cells","ig.PCHiC.Activated_total_CD4_T_cells","ig.PCHiC.Nonactivated_total_CD4_T_cells","ig.PCHiC.Naive_CD8_T_cells","ig.PCHiC.Total_CD8_T_cells","ig.PCHiC.Naive_B_cells","ig.PCHiC.Total_B_cells")

dataAll<-data.frame()
for(i in 1:length(query)){
  ig <- xRDataLoader(RData.customised=query[i], RData.location=RData.location)
  df_edges <- igraph::get.data.frame(ig, what="edges")
  gr_prey <- xGR(data=df_edges$to, format="chr:start-end", RData.location=RData.location)
  annotation.file <- list(hyper=hyperGR, hypo=hypoGR)
  eTerm_prey <- xGRviaGenomicAnno(data.file=gr_prey, annotation.file=annotation.file, background.file=bGR, format.file="GRanges", resolution="region",
                                  RData.location=RData.location,background.annotatable.only = TRUE)
  data<-data.frame(xEnrichViewer(eTerm_prey, 2))
  data$name<-print(query[[i]])
  dataAll<-rbind(dataAll,data)
}

dataAll$direction<-rep(c("demethylated","methylated"),17)
dataAll$name<-gsub("ig.PCHiC.","",dataAll$name)
dataAllPlot<-dataAll[dataAll$direction=="demethylated",]
dataAllPlot$name<-factor(dataAllPlot$name,levels=c(dataAllPlot[order(dataAllPlot$zscore),c("name")]))

colnames(dataAllPlot)[6] <- "P"

promoterEnhancerPlot<-
  ggplot(dataAllPlot,aes(name,zscore))+
  geom_point(aes(fill=-log10(P),size=nOverlap),pch=21,colour="grey30")+
  scale_fill_distiller(palette = "Spectral") +
  coord_flip()+
  theme_bw()+scale_size(range=c(2,20),"nOverlap")+
  theme(legend.text = element_text(size = 20),
        legend.title = element_text(size = 20),
        legend.position="topleft",
        axis.title.y = element_blank(),
        axis.title.x = element_text(size = 20),
        axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20)) +
  scale_x_discrete("Cell type")+
  #guides(size=guide_legend(),fill=guide_legend() )+
  scale_y_continuous("zscore",limits = c(0,9))
promoterEnhancerPlot

ggplot2::ggsave(
  "~/Documents/fairfax_lab/cd14_rnaseq/mQTL/revisions/inflation/HiC_enhancer_enrichment.jpg",
  width = 7,
  height = 12,
  dpi = 100
)

#enrichment of genes proximal to methylated/demthyalted probes for REACTOME & Disease Ontology terms
demeth.sig <- subset(total, corr.qvals<0.05 & t<0)
meth.sig <- subset(total, corr.qvals<0.05 & t>0)
# Demethylated probes
background <- unique(as.character(na.omit(total$gene)))
data <- unique(as.character(na.omit(demeth.sig$gene)))

eTerm<-xEnricherGenes(data=data, background=background, ontology="DO", ontology.algorithm = "lea")
eTermDemeth <- data.frame(xEnrichViewer(eTerm,250, details = T))

eTerm<-xEnricherGenes(data=data, background=background, ontology="REACTOME_ImmuneSystem")
eTermDeMeth.REACTOME_ImmuneSystem <- data.frame(xEnrichViewer(eTerm,100, details = T))

demeth.reactome <- subset(eTermDeMeth.REACTOME_ImmuneSystem, adjp<0.05)
demeth.reactome <- head(demeth.reactome[order(demeth.reactome$pvalue),],15)

demeth.reactome$name[11] <- "Regulation of actin for phagocytic cup formation"
demeth.reactome$name[7] <- "FCGR dependent phagocytosis"
demeth.reactome$name[13] <- "DDX58/IFIH1-mediated induction of IFN-alpha/beta"

demeth.reactome$name <- factor(demeth.reactome$name, levels = c(demeth.reactome$name[order(demeth.reactome$fc)]))

outReactPlot <- ggplot(demeth.reactome, aes(x = name, y = fc)) +
  geom_point(aes(fill = -log10(pvalue), size = nOverlap), pch = 21) +
  scale_fill_distiller(palette = "Spectral") +
  expand_limits(y=c(1,4)) +
  scale_size(range=c(5,25),"nOverlap")+
  scale_y_continuous(breaks=c(0,1,2,3,4,5)) +
  ylab("Fold change") +
  theme_bw() +
  theme(axis.text.y=element_text(size=30), 
        axis.text.x=element_text(size=30), 
        axis.title=element_text(size=30), 
        axis.title.y = element_blank(),
        legend.text = element_text(size = 25),
        legend.title = element_text(size = 30),
        legend.key.size = unit(2, 'cm')) +
  coord_flip() +
  scale_color_manual(values = c("dodgerblue", "firebrick2"))
outReactPlot


demeth.do <- subset(eTermDemeth, adjp<0.05)
demeth.do <- head(demeth.do[order(demeth.do$pvalue),],15)

demeth.do$name <- factor(demeth.do$name, levels = c(demeth.do$name[order(demeth.do$fc)]))

outDOPlot <- ggplot(demeth.do, aes(x = name, y = fc)) +
  geom_point(aes(fill = -log10(pvalue), size = nOverlap), pch = 21) +
  scale_fill_distiller(palette = "Spectral") +
  expand_limits(y=c(1,4)) +
  scale_size(range=c(5,25),"nOverlap")+
  scale_y_continuous(breaks=c(0,1,2,3,4,5)) +
  ylab("Fold change") +
  theme_bw() +
  theme(axis.text.y=element_text(size=30), 
        axis.text.x=element_text(size=30), 
        axis.title=element_text(size=30), 
        axis.title.y = element_blank(),
        legend.text = element_text(size = 25),
        legend.title = element_text(size = 30),
        legend.key.size = unit(2, 'cm')) +
  coord_flip() +
  scale_color_manual(values = c("dodgerblue", "firebrick2"))
outDOPlot


react_do_enrich <- (outReactPlot|outDOPlot)

ggplot2::ggsave(
  "REACTOME_DO_enrichment.jpg",
  width =42,
  height = 14,
  dpi = 100
) 

#Methylated probes
background <- unique(as.character(na.omit(total$gene)))
data <- unique(as.character(na.omit(meth.sig$gene)))

eTerm<-xEnricherGenes(data=data, background=background, ontology="DO", ontology.algorithm = "lea")
eTermMeth <- data.frame(xEnrichViewer(eTerm,100, details = T))

eTerm<-xEnricherGenes(data=data, background=background, ontology="REACTOME_ImmuneSystem")
eTermMeth.REACTOME_ImmuneSystem <- data.frame(xEnrichViewer(eTerm,100, details = T))

#Transcription factor binding site enrichment in LCLs

#LCL data as annotation
all_gr <- xRDataLoader(RData.customised="Uniform_TFBS", RData.location=RData.location)
ind<- grep('_Gm12878_', names(all_gr))
anno_gr <- all_gr[ind]


eTerm_tfbs_hyper <- xGRviaGenomicAnno(data.file=hyperGR, annotation.file=anno_gr, background.file=usedGR, format.file="GRanges", resolution="hybrid", RData.location=RData.location,p.tail = "two-tails")
xEnrichViewer(eTerm_tfbs_hyper, 50)

eTerm_tfbs_hypo <- xGRviaGenomicAnno(data.file=hypoGR, annotation.file=anno_gr, background.file=usedGR, format.file="GRanges", resolution="hybrid", RData.location=RData.location,p.tail = "two-tails")
xEnrichViewer(eTerm_tfbs_hypo, 50)




top_labelled<-eTerm_tfbs_hypo[eTerm_tfbs_hypo$zscore>20|eTerm_tfbs_hypo$zscore<(-20),]
top_labelled$name<-gsub("Tfbs_Haib_Gm12878_","",top_labelled$name)
top_labelled$name<-gsub("Tfbs_Sydh_Gm12878_","",top_labelled$name)

LCLhypoPlot<-ggplot(eTerm_tfbs_hypo[eTerm_tfbs_hypo$adjp<0.05,],aes(zscore,fc))+
  geom_point(aes(fill=-log(adjp,10),size=-log(adjp,10)),pch=21,colour="grey30")+
  scale_fill_gradient2(low="steelblue2",mid="white",high="firebrick2",midpoint = 100,"-log adj. P")+
  theme_bw()+ 
  scale_size(range=c(1,10),"-log adj. P")+
  scale_y_continuous("fold change")+
  guides(size=guide_legend(),fill=guide_legend())+
  geom_text_repel(data = top_labelled, 
                  mapping = aes(label = gsub("sc81335_V0422111","",gsub("iknucla","",gsub("sc137065","",gsub("ab50322_Iggrab","",gsub("_Pcr1x","",gsub("sc6059_Pcr1x","",gsub("_Iggrab","",gsub("_Iggmus","",gsub("sc101553_V0422111","",gsub("_Pcr2x","",name))))))))))),
                  size = 4,
                  fontface = 'italic',
                  color = 'grey40',
                  box.padding = unit(1.4, "lines"),
                  point.padding = unit(1.2, "lines"))

top_labelled<-eTerm_tfbs_hyper[eTerm_tfbs_hyper$zscore>3|eTerm_tfbs_hyper$zscore<(-6),]
top_labelled$name<-gsub("Tfbs_Haib_Gm12878_","",top_labelled$name)
top_labelled$name<-gsub("Tfbs_Sydh_Gm12878_","",top_labelled$name)
LCLhyperPlot<-ggplot(eTerm_tfbs_hyper[eTerm_tfbs_hyper$adjp<0.05,],aes(zscore,fc))+
  geom_point(aes(fill=-log(adjp,10),size=-log(adjp,10)),pch=21,colour="grey30")+
  scale_fill_gradient2(low="steelblue2",mid="white",high="firebrick2",midpoint = 6,"-log adj. P")+
  theme_bw()+ 
  scale_size(range=c(1,10),"-log adj. P")+
  scale_y_continuous("fold change")+
  guides(size=guide_legend(),fill=guide_legend())+
  geom_text_repel(data = top_labelled, 
                  mapping = aes(label = gsub("Tfbs_Uta_Gm12878_","",gsub("sc15914c20","",gsub("Tfbs_Haib_K562_","",gsub("sc137065","",gsub("sc137065_Pcr1x","",gsub("Tfbs_Sydh_K562_","",gsub("_Iggrab","",gsub("_Pcr2x","",name))))))))),
                  size = 4,
                  fontface = 'italic',
                  color = 'grey40',
                  box.padding = unit(1.5, "lines"),
                  point.padding = unit(0.75, "lines"))

p1 <- (LCLhypoPlot|LCLhyperPlot)

ggplot2::ggsave(
    "LCL_tfbs_enrich.jpg",
    width = 14,
    height = 7,
    dpi = 100
  )
