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

#baseline enhanced LPS-induced demethylation at cg19906672 
#read in genotypes and untreated, LPS-treated and delta methylation (LPS-treated minus untreated) for cg19906672
total <- read.table("cg19906672.txt", header = T)

#plot out box plot of per genotype delta methylation at cg19906672
p_labels = data.frame(expt = c("base"), label = c("italic(P)==1.31%*%10^-7"))
cols <- brewer.pal(9,"Set1")
total$state <- "delta"

box.delta.cg19906672 = ggplot(total, aes(x=factor(round(rs6446553,0)), y=delta)) +
geom_dotplot(binaxis="y", binwidth=0.007, stackdir="center", alpha = 0.72) + geom_boxplot(alpha = 0.2) +
scale_x_discrete(breaks=c(0,1,2), labels = c("TT", "TG", "GG")) + aes(fill = factor(round(rs6446553,0)), col=factor(round(rs6446553,0))) + scale_fill_manual(values = cols[c(4,4,4)]) + scale_colour_manual(values = cols[c(4,4,4)]) +
ylab("cg19906672") + xlab("rs6446553") +
theme_bw() + theme(legend.position = "none", axis.text=element_text(size=12),
                             axis.title=element_text(size=15), strip.background = element_rect(color="black", fill=cols[9], size=1.5, linetype="solid"), strip.text.x = element_text(
                               size = 15, color = "white", face = "bold")) + ylim(NA, 0.1) +
facet_wrap( ~ state, ncol=1) +
geom_text(x=2, y=0.075, aes(label=label), data=p_labels, parse=TRUE, inherit.aes=F, size = 5)

box.delta.cg19906672

#plot out box plot of per genotype baseline and LPS-treated methylation at cg19906672
box.facet <- data.frame(rbind(cbind(total$rs6446553, total$ut, "UT"),
                   cbind(total$rs6446553, total$lps, "LPS")))
colnames(box.facet) <- c("rs6446553", "probe", "state")
box.facet$rs6446553 <- as.numeric(box.facet$rs6446553)
box.facet$probe <- as.numeric(box.facet$probe)
box.facet$state <- ordered(box.facet$state, levels = c("UT", "LPS"))

facet.cg19906672.box <- ggplot(box.facet, aes(x=factor(round(rs6446553,0)), y=probe)) +
geom_dotplot(binaxis="y", binwidth=0.002, stackdir="center", alpha = 0.72) + geom_boxplot(alpha = 0.2) +
scale_x_discrete(breaks=c(0,1,2), labels = c("TT", "TG", "GG")) + aes(fill = state, col=state) + scale_fill_manual(values = cols[c(3,1)]) + scale_colour_manual(values = cols[c(3,1)]) +
ylab("cg19906672") + xlab("rs6446553") +
theme_bw() + theme(legend.position = "none", axis.text=element_text(size=12),
                             axis.title=element_text(size=15)) + 
facet_wrap( ~ state, ncol=2) +
theme(strip.background = element_rect(color="black", fill=cols[9], size=1.5, linetype="solid"), strip.text.x = element_text(
    size = 15, color = "white", face = "bold"))
facet.cg19906672.box


p1 <- (facet.cg19906672.box/box.delta.cg19906672)

ggplot2::ggsave(
  "cg19906672_boxes.jpg",
  width = 3,
  height = 5,
  dpi = 300
) 

#RA plot for LPS-modulated mQTL at cg19906672
#read in summary stats and r2 to rs6446553 for region (+/-250kb) & summary stats for regional associations with:
#White Cell count, Neutrophil count, Neutrophil%, Lymphocyte%, MCV, MSCV in UK biobank (Neale lab GWAS v2)
total.stats<- read.table("RA_coloc_cg19906672.txt", header = T)

#split r2 values into bins
total.stats$bin_r2 <- 1
total.stats$bin_r2[which(total.stats$r2>0.2 & total.stats$r2 <= 0.5)] <- 2
total.stats$bin_r2[which(total.stats$r2>0.5 & total.stats$r2 <= 0.8)] <- 3
total.stats$bin_r2[which(total.stats$r2>0.8)] <- 4

#RA plot
#1st plot regional association
total.stats$annotate <- 0
total.stats$annotate[which(total.stats$rsid=="rs6446553")] <- 0
cols2 <- brewer.pal(11,"Spectral")
cols3 <- brewer.pal(8,"Paired")
ra.plot <- ggplot(total.stats, aes(x=var_from, y=-log10(cpg.p))) + 
  xlim(min(total.stats$var_from), max(total.stats$var_from)) +
  annotate("rect", xmin = 6918868, xmax = 6923868, ymin = 0, ymax = 7.5,
           alpha = 1,fill = cols[4]) +
  geom_point(data=subset(total.stats, bin_r2==1), color=cols[9], size=3) + 
  geom_point(data=subset(total.stats, bin_r2==2), color=cols2[5], size=3) + 
  geom_point(data=subset(total.stats, bin_r2==3), color=cols2[3], size=3) + 
  geom_point(data=subset(total.stats, bin_r2==4), color=cols2[1], size=3) +
  scale_y_continuous(name="-log P-value", breaks=c(0,2,4,6,8)) +
  xlab(NULL) + 
  theme_bw() + theme(legend.position = "none", axis.text=element_text(size=12),
                     axis.title=element_text(size=12)) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  geom_text_repel( data=subset(total.stats, annotate==1), aes(label=rsid), size=5, col = c("black"), nudge_x = 50000) +
  annotate("text", x = (min(total.stats$var_from)+20000), y = 8, label = quote(r^2), hjust = 0.5,size=3) +
  annotate("point", x = (min(total.stats$var_from)+60000), y = 7.5, size = 3, colour = cols2[1]) +
  annotate("text", x = (min(total.stats$var_from)+20000), y = 7.5, label = c("0.8-1.0"), hjust = 0.5,size=3) +
  annotate("point", x = (min(total.stats$var_from)+60000), y =7, size = 3, colour = cols2[2]) +
  annotate("text", x = (min(total.stats$var_from)+20000), y = 7, label = c("0.5-0.8"), hjust = 0.5,size=3) +
  annotate("point", x = (min(total.stats$var_from)+60000), y = 6.5, size = 3, colour = cols2[3]) +
  annotate("text", x = (min(total.stats$var_from)+20000), y = 6.5, label = c("0.2-0.5"), hjust = 0.5,size=3) +
  annotate("point", x = (min(total.stats$var_from)+60000), y = 6, size = 3, colour = cols2[5]) +
  annotate("text", x = (min(total.stats$var_from)+20000), y = 6, label = c("0.1-0.3"), hjust = 0.5,size=3) +
  annotate("point", x = (min(total.stats$var_from)+60000), y = 5.5, size = 3, colour = cols[9]) +
  annotate("text", x = (min(total.stats$var_from)+20000), y = 5.5, label = c("<0.2"), hjust = 0.5,size=3) +
  annotate("text", x = 6918868, y = 8, size=5, label = c("cg19906672"), hjust = 0.5)
ra.plot

#2nd plot out genes
gene.ensembl <- useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl", GRCh = 37)
sel.chr=4
sel.pos=min(total.stats$var_from)+(max(total.stats$var_from)-min(total.stats$var_from))/2
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

plot.range <- c(min(total.stats$var_from), max(total.stats$var_from))
genes$order <- rep(seq(1:1),100)[c(1:length(genes$end_position))]
genes.plot <- ggplot(genes, aes(x=start, y=order+1)) + 
  geom_point(size=0) +
  xlim(min(cpg.stats$var_from), max(cpg.stats$var_from)) +
  ylim(c(1.9,2.2)) +
  geom_segment(data = genes[seq(1,18,3),],
               aes(x=start, xend=end, y=order+1, yend=order+1), size = 1, colour = cols3[2],
               arrow = arrow(length = unit(0.25, "cm"))) +
  geom_segment(data = genes[seq(2,18,3),],
               aes(x=start, xend=end, y=order+1.1, yend=order+1.1), size = 1, colour = cols3[2],
               arrow = arrow(length = unit(0.25, "cm"))) +
  geom_segment(data = genes[seq(3,18,3),],
               aes(x=start, xend=end, y=order+1.2, yend=order+1.2), size = 1, colour = cols3[2],
               arrow = arrow(length = unit(0.25, "cm"))) +
  geom_text_repel( data = genes[seq(1,18,3),], aes(x=mid,y=order+1, label=external_gene_name), size=2, col = c("black"),
                   nudge_y =-0.05, segment.color = NA) +
  geom_text_repel( data = genes[seq(2,18,3),], aes(x=mid,y=order+1.1, label=external_gene_name), size=2, col = c("black"),
                   nudge_y =-0.05, segment.color = NA) +
  geom_text_repel( data = genes[seq(3,18,3),], aes(x=mid,y=order+1.2, label=external_gene_name), size=2, col = c("black"),
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
recomb <- read.table("genetic_map_chr4_combined_b37.txt", header = T)
recomb <- subset(recomb, position>min(total.stats$var_from) & position<max(total.stats$var_from))

recomb_rate <- ggplot(recomb, aes(x=position, y=COMBINED_rate.cM.Mb.)) + 
  geom_line() +
  theme_bw() +
  ylab("cM/Mb") +
  xlab("chromosome 4") +
  scale_x_continuous(breaks=c(6700000,6800000,6900000,7000000,7100000),
                     labels=c("6.7Mb", "6.8Mb", "6.9Mb", "7.0Mb", "7.1Mb")) +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"))

#add association, genes and recombination plots together
p1 <- (ra.plot/genes.plot/recomb_rate) + plot_layout(heights = c(3, 1.5, 0.5))

ggplot2::ggsave(
  "cg19906672_ra_plot.jpg",
  width = 5,
  height = 5,
  dpi = 300
) 

#plot out colocalisation of LPS-induced mQTL and GWAS of haematological traits in UK Biobank
for_facet <- data.frame(rbind(cbind(total.stats$cpg.p,total.stats$wbc.p,total.stats$bin_r2,"WBCC"),
                              cbind(total.stats$cpg.p,total.stats$mcv.p,total.stats$bin_r2,"MCV"),
                              cbind(total.stats$cpg.p,total.stats$neut.p,total.stats$bin_r2,"Neut"),
                              cbind(total.stats$cpg.p,total.stats$l_per.p,total.stats$bin_r2,"Lymph%"),
                              cbind(total.stats$cpg.p,total.stats$n_per.p,total.stats$bin_r2,"Neut%"),
                              cbind(total.stats$cpg.p,total.stats$mscv.p,total.stats$bin_r2,"MSCV")))

colnames(for_facet) <- c("cpg.p","gwas.p", "bin_r2","trait")

for_facet$cpg.p <- as.numeric(for_facet$cpg.p)
for_facet$gwas.p <- as.numeric(for_facet$gwas.p)
for_facet$bin_r2 <- as.numeric(for_facet$bin_r2)

pp4_labels = data.frame(trait = c("WBCC", "MCV", "Neut", "Lymph%", "Neut%", "MSCV"),
                        label1 = c("PP4==0.97", "PP4==0.93", "PP4==0.97", "PP4==0.98", "PP4==0.97", "PP4==0.96"))

pp4_labels$trait <- factor(pp4_labels$trait, levels = c("WBCC", "Neut", "Neut%", "Lymph%", "MCV","MSCV"))

for_facet$trait <- factor(for_facet$trait, levels = c("WBCC", "Neut", "Neut%", "Lymph%", "MCV","MSCV"))

facet.coplot <- ggplot(for_facet, aes(x=-log10(gwas.p), y=-log10(cpg.p))) + 
  geom_point(size=2) + 
  aes(fill = factor(bin_r2), col=factor(bin_r2)) + 
  scale_fill_manual(values = c(cols[9],cols2[c(5,3,1)])) + scale_colour_manual(values = c(cols[9],cols2[c(5,3,1)])) +
  ylab("cg19906672 imQTL (-logP)") + 
  xlab("GWAS (-logP)") + 
  facet_wrap( ~ trait, ncol=3) + geom_text(x=0, y=7.5, aes(label=label1), data=pp4_labels, parse=TRUE, inherit.aes=F, size=3, hjust=0) +
  theme_bw() + theme(legend.position = "none", axis.text=element_text(size=10),
                     axis.title=element_text(size=12)) + ylim(NA, 8) +
  theme(strip.background = element_rect(color="black", fill=cols[9], size=1.5, linetype="solid"), strip.text.x = element_text(
    size = 12, color = "white", face = "bold"), plot.title = element_text(size = 20, face = "bold"))

facet.coplot

ggplot2::ggsave(
  "cg19906672_co_plots.jpg",
  width = 3,
  height = 5,
  dpi = 300
) 



#plot out forest plot of rs6446553:G effect on haematological traits to illustrate effect direction
label <- c("WBCC", "Neut", "Neut%", "Lymph%", "MCV","MSCV")
mean  <- c(total.stats$wbc.beta[which(total.stats$rsid=="rs6446553")],
           total.stats$neut.beta[which(total.stats$rsid=="rs6446553")],
           total.stats$n_per.beta[which(total.stats$rsid=="rs6446553")], 
           total.stats$l_per.beta[which(total.stats$rsid=="rs6446553")],
           total.stats$mcv.beta[which(total.stats$rsid=="rs6446553")],
           total.stats$mscv.beta[which(total.stats$rsid=="rs6446553")])
se  <- c(total.stats$wbc.se[which(total.stats$rsid=="rs6446553")],
         total.stats$neut.se[which(total.stats$rsid=="rs6446553")],
         total.stats$n_per.se[which(total.stats$rsid=="rs6446553")], 
         total.stats$l_per.se[which(total.stats$rsid=="rs6446553")],
         total.stats$mcv.se[which(total.stats$rsid=="rs6446553")],
         total.stats$mscv.se[which(total.stats$rsid=="rs6446553")])
lower <- mean-1.96*se
upper <- mean+1.96*se
facet <- "rs6446553:G"

df <- data.frame(label, mean, lower, upper, facet)

cg19906672.fp <- ggplot(data=df, aes(x=label, y=mean, ymin=lower, ymax=upper)) +
  geom_pointrange(fatten=2.5, size = 1) + aes(fill = label, col=label) + scale_fill_manual(values = c(cols[rep(4,6)])) + scale_colour_manual(values = c(cols[rep(4,6)])) +
  geom_hline(yintercept=0, lty=2) +  # add a dotted line at x=1 after flip
  coord_flip() +  # flip coordinates (puts labels on y axis)
  ylab("beta") + scale_x_discrete(expand = expand_scale(add = 1)) +
  scale_y_continuous(breaks=c(-0.02,0,0.02),
                     labels=c("-0.02", "0","0.02")) +
  theme_bw() + theme(legend.position = "none", axis.title.y = element_blank(), axis.text=element_text(size=15),
                     axis.title=element_text(size=15)) +
  facet_wrap(~facet, ncol = 1) +
  theme(strip.background = element_rect(color="black", fill=cols[9], size=1.5, linetype="solid"), strip.text.x = element_text(
    size = 12, color = "white", face = "bold"))
cg19906672.fp

ggplot2::ggsave(
  "cg19906672_haem_fp.jpg",
  width = 2.5,
  height = 5,
  dpi = 300
) 




#LPS-induced mQTL at cg19906672 colocalises with an LPS-induced eQTL (also in monocytes)
#read in genotype and gene expression data for TBC1D14
eqtl <- read.table("eqtl.txt", header = T)

#plot out box plot of per genotype (rs6446553) gene expression of TBC1D14 in untreated and LPS-stimuated monocytes
for.plot <- data.frame(rbind(cbind(eqtl$geno, eqtl$tbc1d14.ut, "UT"),
                  cbind(eqtl$geno, eqtl$tbc1d14.lps24, "LPS - 24h")))

colnames(for.plot) <- c("geno", "exp", "trait")
for.plot$geno <- as.numeric(for.plot$geno)
for.plot$exp <- as.numeric(for.plot$exp)
for.plot$trait <- factor(for.plot$trait, levels = c("UT", "LPS - 24h"))

p_labels = data.frame(trait = c("UT", "LPS - 24h"), label = c("NS","italic(P)==1.64%*%10^-13"))

p_labels$trait <- factor(p_labels$trait, levels = c("UT", "LPS - 24h"))

cols <- brewer.pal(9,"Set1")
exp.plot <- ggplot(for.plot, aes(x=factor(round(geno,0)), y=exp)) + 
geom_dotplot(binaxis="y", binwidth=0.02, stackdir="center", alpha = 0.72) + geom_boxplot(alpha = 0.2) +
scale_x_discrete(breaks=c(0,1,2), labels = c("TT", "TG", "GG")) + aes(fill = factor(trait), col=factor(trait)) + scale_fill_manual(values = cols[c(3,1,1)]) + scale_colour_manual(values = cols[c(3,1,1)]) +
ylab(expression(italic(TBC1D14)~ expression)) + xlab("rs6446553") +
theme_bw() + theme(legend.position = "none", axis.text=element_text(size=8),
                             axis.title=element_text(size=15)) +
  facet_wrap(~trait) + ylim(NA,12) +
  theme(strip.background = element_rect(color="black", fill=cols[9], size=1.5, linetype="solid"), strip.text.x = element_text(
    size = 9, color = "white", face = "bold")) +
geom_text(x=2, y=11.75, aes(label=label), data=p_labels, parse=TRUE, inherit.aes=F, size = 3)
exp.plot

ggplot2::ggsave(
  "eqtl_tbc1d14.box.jpg",
  width = 3.5,
  height = 5,
  dpi = 300
) 

