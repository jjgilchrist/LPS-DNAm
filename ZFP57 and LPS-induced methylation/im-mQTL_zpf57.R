rm(list = ls(all=TRUE))
library(ggplot2)
library(RColorBrewer)
library(dplyr)
library(ggrepel)
library(LDlinkR)
library(ggplotify)
library(patchwork)
library(biomaRt)
library(coloc)

#example of an LPS-modulated mQTL at cg16885113
#read in genotypes and untreated, LPS-treated and delta methylation (LPS-treated minus untreated) for cg16885113
total <- read.table("cg16885113.txt", header = T)

#plot out box plot of per genotype delta methylation at cg16885113
p_labels = data.frame(expt = c("base"), label = c("italic(P)==2.24%*%10^-7"))
cols <- brewer.pal(9,"Set1")
total$state <- "delta"
box.delta.cg16885113 = ggplot(total, aes(x=factor(round(rs3129058,0)), y=delta)) +
geom_dotplot(binaxis="y", binwidth=0.002, stackdir="center", alpha = 0.72) + geom_boxplot(alpha = 0.2) +
scale_x_discrete(breaks=c(0,1,2), labels = c("CC", "CT", "TT")) + aes(fill = factor(round(rs3129058,0)), col=factor(round(rs3129058,0))) + scale_fill_manual(values = cols[c(4,4,4)]) + scale_colour_manual(values = cols[c(4,4,4)]) +
ylab("cg16885113") + xlab("rs3129058") +
theme_bw() + theme(legend.position = "none", axis.text=element_text(size=12),
                             axis.title=element_text(size=15), strip.background = element_rect(color="black", fill=cols[9], size=1.5, linetype="solid"), strip.text.x = element_text(
                               size = 15, color = "white", face = "bold")) + ylim(-0.12, 0.12) +
facet_wrap( ~ state, ncol=1) +
geom_text(x=2, y=0.11, aes(label=label), data=p_labels, parse=TRUE, inherit.aes=F, size = 5)
box.delta.cg16885113

#plot out box plot of per genotype baseline and LPS-treated methylation at cg02724909
box.facet <- data.frame(rbind(cbind(total$rs3129058, total$ut, "UT"),
                   cbind(total$rs3129058, total$lps, "LPS")))
colnames(box.facet) <- c("rs3129058", "probe", "state")
box.facet$rs3129058 <- as.numeric(box.facet$rs3129058)
box.facet$probe <- as.numeric(box.facet$probe)
box.facet$state <- ordered(box.facet$state, levels = c("UT", "LPS"))
p_labels = data.frame(state = c("UT", "LPS"), label = c("italic(P)==1.75%*%10^-37", "italic(P)==6.46%*%10^-33"))
p_labels$state <- ordered(p_labels$state, levels = c("UT", "LPS"))

facet.cg16885113.box <- ggplot(box.facet, aes(x=factor(round(rs3129058,0)), y=probe)) +
geom_dotplot(binaxis="y", binwidth=0.003, stackdir="center", alpha = 0.72) + geom_boxplot(alpha = 0.2) +
scale_x_discrete(breaks=c(0,1,2), labels = c("CC", "CT", "TT")) + aes(fill = state, col=state) + scale_fill_manual(values = cols[c(3,1)]) + scale_colour_manual(values = cols[c(3,1)]) +
ylab("cg16885113") + xlab("rs3129058") +
theme_bw() + theme(legend.position = "none", axis.text=element_text(size=12),
                             axis.title=element_text(size=15)) + 
facet_wrap( ~ state, ncol=2) + ylim(NA, 1) +
  geom_text(x=2, y=0.97, aes(label=label), data=p_labels, parse=TRUE, inherit.aes=F, size = 5) +
theme(strip.background = element_rect(color="black", fill=cols[9], size=1.5, linetype="solid"), strip.text.x = element_text(
    size = 15, color = "white", face = "bold"))
facet.cg16885113.box


p1 <- (facet.cg16885113.box|box.delta.cg16885113) + plot_layout(widths = c(2, 1))

ggplot2::ggsave(
  "cg16885113_box_plots.jpg",
  width = 6,
  height = 7,
  dpi = 300
) 

#colocalisation plot for LPS-modulated mQTL at cg16885113 and Lung Cancer
#read in summary stats and r2 to rs3129058 for region (+/-250kb)
total.stats <- read.table("cg16885113_Lung_CA_colocalisation.txt", header = T)

#split r2 values into bins
total.stats$bin_r2 <- 1
total.stats$bin_r2[which(total.stats$r2>0.2 & total.stats$r2 <= 0.5)] <- 2
total.stats$bin_r2[which(total.stats$r2>0.5 & total.stats$r2 <= 0.8)] <- 3
total.stats$bin_r2[which(total.stats$r2>0.8)] <- 4

#plot colocalisation
cols2 <- brewer.pal(11,"Spectral")
for_facet <- data.frame(cbind(total.stats$cpg.p,total.stats$lca.p,total.stats$bin_r2,"Lung Ca"))
colnames(for_facet) <- c("cpg.p","gwas.p", "bin_r2","trait")
for_facet$cpg.p <- as.numeric(for_facet$cpg.p)
for_facet$gwas.p <- as.numeric(for_facet$gwas.p)
for_facet$bin_r2 <- as.numeric(for_facet$bin_r2)
pp4_labels = data.frame(trait = c("Lung Ca"),
                        label1 = c("PP4==0.94"))
for_facet$trait <- factor(for_facet$trait)

lung_ca.coplot <- ggplot(for_facet, aes(x=-log10(gwas.p), y=-log10(cpg.p))) + 
  geom_point(size=3) + 
  aes(fill = factor(bin_r2), col=factor(bin_r2)) + 
  scale_fill_manual(values = c(cols[9],cols2[c(5,3,1)])) + scale_colour_manual(values = c(cols[9],cols2[c(5,3,1)])) +
  ylab("cg16885113 imQTL (-logP)") + 
  xlab("GWAS (-logP)") + 
  facet_wrap( ~ trait, ncol=1, scales = "free_x") + geom_text(x=0, y=7, aes(label=label1), data=pp4_labels, parse=TRUE, inherit.aes=F, size=5, hjust=0) +
  theme_bw() + theme(legend.position = "none", axis.text=element_text(size=12),
                     axis.title=element_text(size=12)) + ylim(NA,8) +
  theme(strip.background = element_rect(color="black", fill=cols[9], size=1.5, linetype="solid"), strip.text.x = element_text(
    size = 12, color = "white", face = "bold"), plot.title = element_text(size = 20, face = "bold"))
lung_ca.coplot

ggsave(
  "cg16885113_lung_ca_coplots.jpg",
  width = 2,
  height = 5,
  dpi = 300
) 

#RA plot for LPS-modulated mQTL at cg16885113
#read in summary stats and r2 to rs3129058 for region (+/-250kb)
cpg.stats <- read.table("cg16885113_RA.txt", header = T)

#split r2 values into bins
cpg.stats$bin_r2 <- 1
cpg.stats$bin_r2[which(cpg.stats$r2>0.2 & cpg.stats$r2 <= 0.5)] <- 2
cpg.stats$bin_r2[which(cpg.stats$r2>0.5 & cpg.stats$r2 <= 0.8)] <- 3
cpg.stats$bin_r2[which(cpg.stats$r2>0.8)] <- 4

#RA plot
#set colours
cols3 <- brewer.pal(8,"Paired")
cols2 <- brewer.pal(11,"Spectral")
cols <- brewer.pal(9,"Set1")

#1st plot regional association
cpg.stats$annotate <- 0
cpg.stats$annotate[c(which(cpg.stats$rsid=="rs3129058"))] <- 1


ra.plot <- ggplot(cpg.stats, aes(x=var_from, y=-log10(nom_pval))) + 
  xlim(min(cpg.stats$var_from), max(cpg.stats$var_from)) +
  annotate("rect", xmin = 29646507, xmax = 29650507, ymin = 0, ymax = 9.5,
           alpha = 1,fill = cols[4]) +
  geom_point(data=subset(cpg.stats, bin_r2==1), color=cols[9], size=3) + 
  geom_point(data=subset(cpg.stats, bin_r2==2), color=cols2[5], size=3) + 
  geom_point(data=subset(cpg.stats, bin_r2==3), color=cols2[3], size=3) + 
  geom_point(data=subset(cpg.stats, bin_r2==4), color=cols2[1], size=3) +
  scale_y_continuous(name="-log P-value", breaks=c(0,2,4,6,8,10)) +
  xlab(NULL) + 
  theme_bw() + theme(legend.position = "none", axis.text=element_text(size=12),
                     axis.title=element_text(size=12)) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  geom_text_repel( data=subset(cpg.stats, annotate==1), aes(label=rsid), size=5, col = c("black"), nudge_x = 50000) +
  annotate("text", x = (min(cpg.stats$var_from)+20000), y = 10, label = quote(r^2), hjust = 0.5,size=3) +
  annotate("point", x = (min(cpg.stats$var_from)+60000), y = 9.5, size = 3, colour = cols2[1]) +
  annotate("text", x = (min(cpg.stats$var_from)+20000), y = 9.5, label = c("0.8-1.0"), hjust = 0.5,size=3) +
  annotate("point", x = (min(cpg.stats$var_from)+60000), y =9, size = 3, colour = cols2[2]) +
  annotate("text", x = (min(cpg.stats$var_from)+20000), y = 9, label = c("0.5-0.8"), hjust = 0.5,size=3) +
  annotate("point", x = (min(cpg.stats$var_from)+60000), y = 8.5, size = 3, colour = cols2[3]) +
  annotate("text", x = (min(cpg.stats$var_from)+20000), y = 8.5, label = c("0.2-0.5"), hjust = 0.5,size=3) +
  annotate("point", x = (min(cpg.stats$var_from)+60000), y = 8, size = 3, colour = cols2[5]) +
  annotate("text", x = (min(cpg.stats$var_from)+20000), y = 8, label = c("0.1-0.3"), hjust = 0.5,size=3) +
  annotate("point", x = (min(cpg.stats$var_from)+60000), y = 7.5, size = 3, colour = cols[9]) +
  annotate("text", x = (min(cpg.stats$var_from)+20000), y = 7.5, label = c("<0.2"), hjust = 0.5,size=3) +
  annotate("text", x = 29648507, y = 10, size=5, label = c("cg16885113"), hjust = 0.5)
ra.plot

#2nd plot local genes
gene.ensembl <- useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl", GRCh = 37)
sel.chr=6
sel.pos=min(cpg.stats$var_from)+(max(cpg.stats$var_from)-min(cpg.stats$var_from))/2
range=1000000

out.bm.genes.region <- getBM(
  attributes = c('start_position','end_position','ensembl_gene_id','external_gene_name', 'gene_biotype', 'strand'), 
  filters = c('chromosome_name','start','end'), 
  values = list(sel.chr, sel.pos - range, sel.pos + range), 
  mart = gene.ensembl)

out.bm.genes.region$mid <- out.bm.genes.region$start_position+(out.bm.genes.region$end_position-out.bm.genes.region$start_position)/2
genes <- subset(out.bm.genes.region, gene_biotype=="protein_coding")
genes$start <- genes$start_position
genes$end <- genes$end_position
genes$start[which(genes$strand==-1)] <- genes$end_position[which(genes$strand==-1)]
genes$end[which(genes$strand==-1)] <- genes$start_position[which(genes$strand==-1)]


plot.range <- c(min(cpg.stats$var_from), max(cpg.stats$var_from))
genes$order <- rep(seq(1:1),100)[c(1:length(genes$end_position))]
genes.plot <- ggplot(genes, aes(x=start, y=order+1)) + 
  geom_point(size=0) +
  xlim(min(cpg.stats$var_from), max(cpg.stats$var_from)) +
  ylim(c(1.9,2.1)) +
  geom_segment(data = genes[seq(1,45,2),],
               aes(x=start, xend=end, y=order+1, yend=order+1), size = 1, colour = cols3[2],
               arrow = arrow(length = unit(0.25, "cm"))) +
  geom_segment(data = genes[seq(2,45,2),],
               aes(x=start, xend=end, y=order+1.1, yend=order+1.1), size = 1, colour = cols3[2],
               arrow = arrow(length = unit(0.25, "cm"))) +
  geom_text_repel( data = genes[seq(1,45,2),], aes(x=mid,y=order+1, label=external_gene_name), size=2, col = c("black"),
                   nudge_y =-0.05, segment.color = NA) +
  geom_text_repel( data = genes[seq(2,45,2),], aes(x=mid,y=order+1.1, label=external_gene_name), size=2, col = c("black"),
                   nudge_y =-0.05, segment.color = NA) +
  theme_bw() +
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
genes.plot

#3rd plot out recombination rate
recomb <- read.table("genetic_map_chr6_combined_b37.txt", header = T)
recomb <- subset(recomb, position>min(cpg.stats$var_from) & position<max(cpg.stats$var_from))

recomb_rate <- ggplot(recomb, aes(x=position, y=COMBINED_rate.cM.Mb.)) + 
  geom_line() +
  theme_bw() +
  ylab("cM/Mb") +
  xlab("chromosome 1") +
  scale_x_continuous(breaks=c(29000000,29500000,30000000,30500000),
                     labels=c("29.0Mb", "29.5Mb", "30.0Mb", "30.5Mb")) +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"))
recomb_rate

#add association, genes and recombination plots together
p1 <- (ra.plot/genes.plot/recomb_rate) + plot_layout(heights = c(3, 0.5, 0.5))

ggsave(
  "cg16885113_im-mQTL_RA_plot.jpg",
  width = 7,
  height = 5,
  dpi = 300
) 

##LPS-induced mQTL at cg16885113 colocalises with an eQTL for zfp57 in untreated monocytes
#read in genotype and gene expression data for zfp57
eqtl <- read.table("zfp57_eqtl.txt", header = T)

#plot out box plot of per genotype (rs3129051) gene expression of ZFP57 in untreated and LPS-stimuated monocytes
for.plot <- data.frame(rbind(cbind(eqtl$geno, eqtl$zfp57.ut, "UT"),
                  cbind(eqtl$geno, eqtl$zfp57.lps24, "LPS - 24h")))
colnames(for.plot) <- c("geno", "exp", "trait")
for.plot$geno <- as.numeric(for.plot$geno)
for.plot$exp <- as.numeric(for.plot$exp)
for.plot$trait <- factor(for.plot$trait, levels = c("UT", "LPS - 24h"))
p_labels = data.frame(trait = c("UT", "LPS - 24h"), label = c("italic(P)==1.25%*%10^-19", "NS"))
p_labels$trait <- factor(p_labels$trait, levels = c("UT", "LPS - 24h"))
cols <- brewer.pal(9,"Set1")

zfp57.eqtl.plot <- ggplot(for.plot, aes(x=factor(round(geno,0)), y=exp)) + 
geom_dotplot(binaxis="y", binwidth=0.004, stackdir="center", alpha = 0.72) + geom_boxplot(alpha = 0.2) +
scale_x_discrete(breaks=c(0,1,2), labels = c("CC", "CT", "TT")) + aes(fill = factor(trait), col=factor(trait)) + scale_fill_manual(values = cols[c(3,1)]) + scale_colour_manual(values = cols[c(3,1)]) +
ylab(expression(italic(ZFP57)~ expression)) + xlab("rs3129058") +
theme_bw() + theme(legend.position = "none", axis.text=element_text(size=12),
                             axis.title=element_text(size=15)) +
  facet_wrap(~trait) +
  theme(strip.background = element_rect(color="black", fill=cols[9], size=1.5, linetype="solid"), strip.text.x = element_text(
    size = 15, color = "white", face = "bold")) +
geom_text(x=2, y=7.65, aes(label=label), data=p_labels, parse=TRUE, inherit.aes=F, size = 7)
zfp57.eqtl.plot

ggplot2::ggsave(
  "eqtl_zfp57.box.jpg",
  width = 5,
  height = 7,
  dpi = 300
) 

#RA plot for zfp57 eqtl (untreated)
#read in summary stats and r2 to rs417764 for region (+/-250kb)
for.ra <- read.table("zfp57_monoycte_eqtl_UT.txt", header = T)

#split r2 values into bins
for.ra$bin_r2 <- 1
for.ra$bin_r2[which(for.ra$r2>0.2 & for.ra$r2 <= 0.5)] <- 2
for.ra$bin_r2[which(for.ra$r2>0.5 & for.ra$r2 <= 0.8)] <- 3
for.ra$bin_r2[which(for.ra$r2>0.8)] <- 4

#RA plot
#set colours
cols3 <- brewer.pal(8,"Paired")
cols2 <- brewer.pal(11,"Spectral")
cols <- brewer.pal(9,"Set1")

#1st plot regional association
for.ra$annotate <- 0
for.ra$annotate[c(which(for.ra$id=="rs417764"))] <- 1
for.ra$annotate[c(which(for.ra$id=="rs375984"))] <- 1
for.ra$annotate[c(which(for.ra$id=="rs3129058"))] <- 1

ra.plot <- ggplot(for.ra, aes(x=pos, y=-log10(p))) + 
  xlim(min(for.ra$pos), max(for.ra$pos)) +
  geom_point(data=subset(for.ra, bin_r2==1), color=cols[9], size=3) + 
  geom_point(data=subset(for.ra, bin_r2==2), color=cols2[5], size=3) + 
  geom_point(data=subset(for.ra, bin_r2==3), color=cols2[3], size=3) + 
  geom_point(data=subset(for.ra, bin_r2==4), color=cols2[1], size=3) +
  scale_y_continuous(name="-log P-value", breaks=c(0,5,10,15,20,25)) +
  xlab(NULL) + 
  theme_bw() + theme(legend.position = "none", axis.text=element_text(size=12),
                     axis.title=element_text(size=12)) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  geom_text_repel( data=subset(for.ra, annotate==1), aes(label=id), size=4, col = c("black"), nudge_x = 75000) +
  annotate("text", x = (min(for.ra$pos)+20000), y = 24, label = quote(r^2), hjust = 0.5,size=3) +
  annotate("point", x = (min(for.ra$pos)+60000), y = 23, size = 3, colour = cols2[1]) +
  annotate("text", x = (min(for.ra$pos)+20000), y = 23, label = c("0.8-1.0"), hjust = 0.5,size=3) +
  annotate("point", x = (min(for.ra$pos)+60000), y =22, size = 3, colour = cols2[2]) +
  annotate("text", x = (min(for.ra$pos)+20000), y = 22, label = c("0.5-0.8"), hjust = 0.5,size=3) +
  annotate("point", x = (min(for.ra$pos)+60000), y = 21, size = 3, colour = cols2[3]) +
  annotate("text", x = (min(for.ra$pos)+20000), y = 21, label = c("0.2-0.5"), hjust = 0.5,size=3) +
  annotate("point", x = (min(for.ra$pos)+60000), y = 20, size = 3, colour = cols2[5]) +
  annotate("text", x = (min(for.ra$pos)+20000), y = 20, label = c("0.1-0.3"), hjust = 0.5,size=3) +
  annotate("point", x = (min(for.ra$pos)+60000), y = 19, size = 3, colour = cols[9]) +
  annotate("text", x = (min(for.ra$pos)+20000), y = 19, label = c("<0.2"), hjust = 0.5,size=3)
ra.plot


#2nd plot out location of genes in region
library(biomaRt)
gene.ensembl <- useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl", GRCh = 37)
sel.chr=6
sel.pos=min(for.ra$pos)+(max(for.ra$pos)-min(for.ra$pos))/2
range=1000000
out.bm.genes.region <- getBM(
  attributes = c('start_position','end_position','ensembl_gene_id','external_gene_name', 'gene_biotype', 'strand'), 
  filters = c('chromosome_name','start','end'), 
  values = list(sel.chr, sel.pos - range, sel.pos + range), 
  mart = gene.ensembl)
out.bm.genes.region$mid <- out.bm.genes.region$start_position+(out.bm.genes.region$end_position-out.bm.genes.region$start_position)/2
genes <- subset(out.bm.genes.region, gene_biotype=="protein_coding")
genes$start <- genes$start_position
genes$end <- genes$end_position
genes$start[which(genes$strand==-1)] <- genes$end_position[which(genes$strand==-1)]
genes$end[which(genes$strand==-1)] <- genes$start_position[which(genes$strand==-1)]

plot.range <- c(min(for.ra$pos), max(for.ra$pos))
genes$order <- rep(seq(1:1),100)[c(1:length(genes$end_position))]
genes.plot <- ggplot(genes, aes(x=start, y=order+1)) + 
  geom_point(size=0) +
  xlim(min(for.ra$pos), max(for.ra$pos)) +
  ylim(c(1.9,2.1)) +
  geom_segment(data = genes[seq(14,23,2),],
               aes(x=start, xend=end, y=order+1, yend=order+1), size = 1, colour = cols3[2],
               arrow = arrow(length = unit(0.25, "cm"))) +
  geom_segment(data = genes[seq(13,23,2),],
               aes(x=start, xend=end, y=order+1.1, yend=order+1.1), size = 1, colour = cols3[2],
               arrow = arrow(length = unit(0.25, "cm"))) +
  geom_text_repel( data = genes[seq(14,23,2),], aes(x=mid,y=order+1, label=external_gene_name), size=2, col = c("black"),
                   nudge_y =-0.05, segment.color = NA) +
  geom_text_repel( data = genes[seq(13,23,2),], aes(x=mid,y=order+1.1, label=external_gene_name), size=2, col = c("black"),
                   nudge_y =-0.05, segment.color = NA) +
  theme_bw() +
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
genes.plot

#3rd plot out local recombination rate
recomb <- read.table("genetic_map_chr6_combined_b37.txt", header = T)
recomb <- subset(recomb, position>min(for.ra$pos) & position<max(for.ra$pos))
recomb_rate <- ggplot(recomb, aes(x=position, y=COMBINED_rate.cM.Mb.)) + 
  geom_line() +
  theme_bw() +
  ylab("cM/Mb") +
  xlab("chromosome 6") +
  scale_x_continuous(breaks=c(29000000,29500000,30000000,30500000),
                     labels=c("29.0Mb", "29.5Mb", "30.0Mb", "30.5Mb")) +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"))
recomb_rate

#add association, genes and recombination plots together
p1 <- (ra.plot/genes.plot/recomb_rate) + plot_layout(heights = c(3, 0.5, 0.5))

ggplot2::ggsave(
  "zfp57.eqtl_ra_plot.jpg",
  width = 5,
  height = 5,
  dpi = 300
) 

#circos plot of effect of ZFP57 eQTL on genome-wide methylation in monocytes
library(GenomicRanges)
library(ggbio)
library(stringr)

#read in data
zfp57 <- read.table("baseline_ut_zfp57.txt", header = T)

#subset into unimprinted and imprinted genetic regions and make Granges objects
bed.non.imp <- data.frame(subset(zfp57, imprinted==0)[,c("chr", "start", "end")])
colnames(bed.non.imp) <- c("chr", "start", "end")
bed.non.imp$chr.b <- str_remove_all(bed.non.imp$chr, "chr")
bed.non.imp$chr.b <- as.numeric(bed.non.imp$chr.b)
bed.non.imp.sort <- bed.non.imp[order(bed.non.imp$chr.b, bed.non.imp$start),]
non_imp.gr <- with(bed.non.imp.sort,GRanges(chr.b, IRanges(start, end)))

bed.imp <- data.frame(subset(zfp57, imprinted==1)[,c("chr", "start", "end")])
colnames(bed.imp) <- c("chr", "start", "end")
bed.imp$chr.b <- str_remove_all(bed.imp$chr, "chr")
bed.imp$chr.b <- as.numeric(bed.imp$chr.b)
bed.imp.sort <- bed.imp[order(bed.imp$chr.b, bed.imp$start),]
imp.gr <- with(bed.imp.sort,GRanges(chr.b, IRanges(start, end)))

data("hg19Ideogram", package = "biovizBase")
seqlevelsStyle(hg19Ideogram) <- "NCBI"
seqlevels(hg19Ideogram, pruning.mode="coarse") <- c(seq(1:22), "X", "Y")

gr2 <- GRanges(c("6"), IRanges(29648137, width = 1))
seqlevels(gr2, pruning.mode="coarse") <- c(seq(1:22), "X", "Y")
seqlengths(gr2) <- seqlengths(hg19Ideogram)

values(imp.gr)$to.gr <- gr2
values(non_imp.gr)$to.gr <- gr2

seqlevels(imp.gr, pruning.mode="coarse") <- c(seq(1:22), "X", "Y")
seqlengths(imp.gr) <- seqlengths(hg19Ideogram)

seqlevels(non_imp.gr, pruning.mode="coarse") <- c(seq(1:22), "X", "Y")
seqlengths(non_imp.gr) <- seqlengths(hg19Ideogram)

#plot circos
zfp57.circos <- ggplot() + layout_circle(hg19Ideogram, geom = "ideo", fill = "gray70",
                                         radius = 30, trackWidth = 4)
zfp57.circos <- zfp57.circos + layout_circle(hg19Ideogram, geom = "scale", size = 2, radius = 35, trackWidth= 2)
zfp57.circos <- zfp57.circos + layout_circle(hg19Ideogram, geom = "text", aes(label = seqnames), vjust = 0, radius = 38, trackWidth = 7)
zfp57.circos <- zfp57.circos + layout_circle(non_imp.gr, geom = "link", linked.to = "to.gr", color = "dark grey", radius=27, trackwidth=6, size=2)
zfp57.circos <- zfp57.circos + layout_circle(imp.gr, geom = "link", linked.to = "to.gr", color = cols[5], radius=27, trackwidth=6, size=2)
zfp57.circos

ggplot2::ggsave(
  "zfp57.circos.jpg",
  width = 4,
  height = 4,
  dpi = 300
) 

#correlation between ZFP57 eQTL-associated delta methylation and baseline
zfp57.baseline <- read.table("zfp57_baseline_delta_methylation.txt", header = T)

#plot baseline:delta methylation scatter plot (n=151 CpGs)
zfp57_baseline.scatter <- ggplot(zfp57.baseline, aes(x=base.med, y=beta)) +
  geom_point(alpha = 0.72, col = cols[4]) +
  geom_smooth(method=lm, fill=cols[4], col="black") +
  ylab("rs417764:methylation beta") +
  xlab("baseline methylation") +
  scale_x_continuous(breaks=c(0,0.5,1.0), limits = c(0,1)) +
  theme_bw() +
  theme(legend.position = "none", axis.text=element_text(size=8),
        axis.title=element_text(size=12),#panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),axis.line = element_line(colour = "grey"))
zfp57_baseline.scatter

ggplot2::ggsave(
  "zpf57_beta_scatter.jpg",
  width = 2.5,
  height = 2.5,
  dpi = 300
) 

#baseline CpG methylation densities per ZFP57 eQTL genotype
#read in baseline methylation at ZFP57-associated CPGs
zfp57.baseline <- read.table("zpf57_associated_cpgs.baseline_cpg_methylation.txt", header = T)

#split out into estimated genotypes
zfp57.baseline.geno0 <- subset(zfp57.baseline, round(geno)==0)
zfp57.baseline.geno1 <- subset(zfp57.baseline, round(geno)==1)
zfp57.baseline.geno2 <- subset(zfp57.baseline, round(geno)==2)

base.med.geno0 <- NA
for (i in c(1:151)){
  base.med.geno0[i] <- median(zfp57.baseline.geno0[,i])
}
base.med.geno1 <- NA
for (i in c(1:151)){
  base.med.geno1[i] <- median(zfp57.baseline.geno1[,i])
}
base.med.geno2 <- NA
for (i in c(1:151)){
  base.med.geno2[i] <- median(zfp57.baseline.geno2[,i])
}

#density plot
library(ggridges)
for.dens.plot <- data.frame(rbind(cbind(base.med.geno0, "GG"),
                                  cbind(base.med.geno1, "GA"),
                                  cbind(base.med.geno2, "AA")))

colnames(for.dens.plot) <- c("beta", "rs417764")
for.dens.plot$beta <- as.numeric(for.dens.plot$beta)
for.dens.plot$rs417764 <- factor(for.dens.plot$rs417764, levels = c("GG", "GA", "AA"))

cols.bupu <- brewer.pal(9, "BuPu")

zfp57.dens.plot <- ggplot(for.dens.plot, aes(x = beta, y = rs417764, group = rs417764)) + 
  geom_density_ridges() +
  aes(fill = rs417764, col=rs417764) + scale_fill_manual(values = cols.bupu[c(5,6,7)]) + scale_colour_manual(values = cols.bupu[c(5,6,7)]) +
  scale_x_continuous(breaks=c(0,0.5,1.0)) +
  theme_bw() + theme(legend.position = "none", axis.text=element_text(size=8),
                     axis.title=element_text(size=12),panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(),
                     panel.border = element_blank(),
                     panel.background = element_blank(),axis.line = element_line(colour = "grey"))
zfp57.dens.plot
ggplot2::ggsave(
  "zpf57_baseline_densities.jpg",
  width = 2.5,
  height = 2.5,
  dpi = 300
) 

#effect of ZFP57 eQTL on LPS-induced methylation in trans
#read in per genotype delta methylation values at 7 significantly associated CpGs
for.plot <- read.table("delta_cpg_zfp57.txt", header = T)

#plot effect of ZFP57 eQTL genotype on LPS-induced methylation in trans at 7 CpGs
cols <- brewer.pal(9,"Set1")

cpg.plot <- ggplot(for.plot, aes(x=factor(round(geno,0)), y=d.meth)) + 
  geom_dotplot(binaxis="y", binwidth=0.002, stackdir="center", alpha = 0.72) + geom_boxplot(alpha = 0.2) +
  scale_x_discrete(breaks=c(0,1,2), labels = c("GG", "GA", "AA")) + aes(fill = factor(cpg), col=factor(cpg)) + scale_fill_manual(values = cols[rep(4,7)]) + scale_colour_manual(values = cols[rep(4,7)]) +
  ylab("delta methylation") + xlab("rs417764") +
  theme_bw() + theme(legend.position = "none", axis.text=element_text(size=7),
                     axis.title=element_text(size=15)) +
  facet_wrap(~cpg, ncol=7) +
  theme(strip.background = element_rect(color="black", fill=cols[9], size=1.5, linetype="solid"), strip.text.x = element_text(
    size = 6, color = "white", face = "bold"))
cpg.plot
ggplot2::ggsave(
  "zpf57_delta_methylation.jpg",
  width = 2.5,
  height = 2.5,
  dpi = 300
) 









