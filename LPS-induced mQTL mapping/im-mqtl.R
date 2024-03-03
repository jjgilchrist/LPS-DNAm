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

#read in summary statistics describing LPS-modulated mQTL in monocytes (n=7,369 CpGs)
#to limit file size, SNP:CpG associations are limited to p<0.01 and a random subset (n=100,000) of SNP:CpG associations where p>0.01
manh.data <- read.table("cis_im-mqtl.txt.gz", header = T)

#construct manhattan pot of cis mapping
manh.data$chr <- factor(manh.data$chr, levels = c(1:22))

don <- manh.data %>% 
  
  # Compute chromosome size
  group_by(chr) %>% 
  summarise(chr_len=max(bp)) %>% 
  
  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(as.numeric(chr_len))-chr_len) %>%
  dplyr::select(-chr_len) %>%
  
  # Add this info to the initial dataset
  left_join(manh.data, ., by=c("chr"="chr")) %>%
  
  # Add a cumulative position of each SNP
  arrange(chr, bp) %>%
  mutate( BPcum=bp+tot)

axisdf = don %>% group_by(chr) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )

mqtl.manh <- ggplot(don, aes(x=BPcum, y=-log10(pval))) +
  
  # Show all points
  
  # Add highlighted points
  geom_point(data=subset(don, cols=="demeth.sig"), color="steelblue2", size=3) +
  geom_point(data=subset(don, cols=="meth.sig"), color="firebrick2", size=3) +
  geom_point(data = subset(don, sig==0), aes(color=as.factor(chr)), size=3, alpha= 0.75) +
  scale_color_manual(values = rep(c("grey70","grey30"), 11 )) +
  
  # custom X axis:
  scale_x_continuous( label = axisdf$chr[c(1:14,16,18,20,22)], breaks= axisdf$center[c(1:14,16,18,20,22)] ) +
  scale_y_continuous(expand = c(0, 0) ) + 
  ylim(NA, 50) +     # remove space between plot area and x axis
  xlab("chromosome") +
  
  # Custom the theme:
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
  "im-mqtl_manhattan.jpg",
  width = 14,
  height = 7,
  dpi = 300
)

#example of a de novo LPS-modulated mQTL at cg02724909
#read in genotypes and untreated, LPS-treated and delta methylation (LPS-treated minus untreated) forcg02724909
total <- read.table("cg02724909.txt", header = T)

#plot out box plot of per genotype delta methylation at cg02724909
p_labels = data.frame(expt = c("base"), label = c("italic(P)==1.63%*%10^-29"))
cols <- brewer.pal(9,"Set1")
total$state <- "delta"

box.delta.cg02724909 = ggplot(total, aes(x=factor(round(rs869191,0)), y=delta)) +
  geom_dotplot(binaxis="y", binwidth=0.003, stackdir="center", alpha = 0.72) + geom_boxplot(alpha = 0.2) +
  scale_x_discrete(breaks=c(0,1,2), labels = c("CC", "CG", "GG")) + aes(fill = factor(round(rs869191,0)), col=factor(round(rs869191,0))) + scale_fill_manual(values = cols[c(4,4,4)]) + scale_colour_manual(values = cols[c(4,4,4)]) +
  ylab("cg02724909") + xlab("rs869191") +
  theme_bw() + theme(legend.position = "none", axis.text=element_text(size=12),
                     axis.title=element_text(size=15), strip.background = element_rect(color="black", fill=cols[9], size=1.5, linetype="solid"), strip.text.x = element_text(
                       size = 15, color = "white", face = "bold")) + ylim(NA, 0.1) +
  facet_wrap( ~ state, ncol=1) +
  geom_text(x=2, y=0.075, aes(label=label), data=p_labels, parse=TRUE, inherit.aes=F, size = 5)

box.delta.cg02724909


#plot out box plot of per genotype baseline and LPS-treated methylation at cg02724909
box.facet <- data.frame(rbind(cbind(total$rs869191, total$ut, "UT"),
                              cbind(total$rs869191, total$lps, "LPS")))
colnames(box.facet) <- c("rs869191", "probe", "state")
box.facet$rs869191 <- as.numeric(box.facet$rs869191)
box.facet$probe <- as.numeric(box.facet$probe)

box.facet$state <- ordered(box.facet$state, levels = c("UT", "LPS"))

p_labels = data.frame(state = c("UT", "LPS"), label = c("NS", "italic(P)==2.65%*%10^-31"))
p_labels$state <- ordered(p_labels$state, levels = c("UT", "LPS"))

facet.cg02724909.box <- ggplot(box.facet, aes(x=factor(round(rs869191,0)), y=probe)) +
  geom_dotplot(binaxis="y", binwidth=0.002, stackdir="center", alpha = 0.72) + geom_boxplot(alpha = 0.2) +
  scale_x_discrete(breaks=c(0,1,2), labels = c("CC", "CG", "GG")) + aes(fill = state, col=state) + scale_fill_manual(values = cols[c(3,1)]) + scale_colour_manual(values = cols[c(3,1)]) +
  ylab("cg02724909") + xlab("rs869191") +
  theme_bw() + theme(legend.position = "none", axis.text=element_text(size=12),
                     axis.title=element_text(size=15)) + 
  facet_wrap( ~ state, ncol=2) +
  ylim(NA, 0.9) +
  geom_text(x=2, y=0.85, aes(label=label), data=p_labels, parse=TRUE, inherit.aes=F, size = 5) +
  theme(strip.background = element_rect(color="black", fill=cols[9], size=1.5, linetype="solid"), strip.text.x = element_text(
    size = 15, color = "white", face = "bold"))


facet.cg02724909.box

p1 <- (facet.m.box|box.cg02724909) + plot_layout(widths = c(2, 1))

ggplot2::ggsave(
  "cg02724909_de_novo_box.jpg",
  width = 7,
  height = 7,
  dpi = 300
) 

#RA plot for LPS-modulated mQTL at cg02724909
#read in summary stats and r2 to rs869191 for region (+/-250kb)
total.stats <- read.table("cg02724909_region.txt", header = T)

#split r2 values into bins
total.stats$bin_r2 <- 1
total.stats$bin_r2[which(total.stats$r2>0.1 & total.stats$r2 <= 0.3)] <- 2
total.stats$bin_r2[which(total.stats$r2>0.3 & total.stats$r2 <= 0.5)] <- 3
total.stats$bin_r2[which(total.stats$r2>0.5 & total.stats$r2 <= 0.8)] <- 4
total.stats$bin_r2[which(total.stats$r2>0.8)] <- 5

#RA plot
#set colours
cols3 <- brewer.pal(8,"Paired")
cols2 <- brewer.pal(11,"Spectral")
cols <- brewer.pal(9,"Set1")

#1st plot regional association
total.stats$annotate <- 0
total.stats$annotate[which(total.stats$rsid=="rs869191")] <- 0

ra.plot <- ggplot(total.stats, aes(x=var_from, y=-log10(cpg.p))) + 
  xlim(min(total.stats$var_from), max(total.stats$var_from)) +
  annotate("rect", xmin = 149625104, xmax = 149630104, ymin = 0, ymax = 37,
           alpha = 1,fill = cols[4]) +
  geom_point(data=subset(total.stats, bin_r2==1), color=cols[9], size=3) + 
  geom_point(data=subset(total.stats, bin_r2==2), color=cols2[5], size=3) + 
  geom_point(data=subset(total.stats, bin_r2==3), color=cols2[3], size=3) + 
  geom_point(data=subset(total.stats, bin_r2==4), color=cols2[2], size=3) +
  geom_point(data=subset(total.stats, bin_r2==5), color=cols2[1], size=3) +
  scale_y_continuous(name="-log P-value", breaks=c(0,5,10,15,20,25,30,35,40)) +
  xlab(NULL) + 
  theme_bw() + theme(legend.position = "none", axis.text=element_text(size=12),
                     axis.title=element_text(size=12)) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  geom_text_repel( data=subset(total.stats, annotate==1), aes(label=rsid), size=5, col = c("black"), nudge_x = 50000) +
  annotate("text", x = (min(total.stats$var_from)+20000), y = 35, label = quote(r^2), hjust = 0.5,size=4) +
  annotate("point", x = (min(total.stats$var_from)+60000), y = 33, size = 4, colour = cols2[1]) +
  annotate("text", x = (min(total.stats$var_from)+20000), y = 33, label = c("0.8-1.0"), hjust = 0.5,size=4) +
  annotate("point", x = (min(total.stats$var_from)+60000), y =31, size = 4, colour = cols2[2]) +
  annotate("text", x = (min(total.stats$var_from)+20000), y = 31, label = c("0.5-0.8"), hjust = 0.5,size=4) +
  annotate("point", x = (min(total.stats$var_from)+60000), y = 29, size = 4, colour = cols2[3]) +
  annotate("text", x = (min(total.stats$var_from)+20000), y = 29, label = c("0.3-0.5"), hjust = 0.5,size=4) +
  annotate("point", x = (min(total.stats$var_from)+60000), y = 27, size = 4, colour = cols2[5]) +
  annotate("text", x = (min(total.stats$var_from)+20000), y = 27, label = c("0.1-0.3"), hjust = 0.5,size=4) +
  annotate("point", x = (min(total.stats$var_from)+60000), y = 25, size = 4, colour = cols[9]) +
  annotate("text", x = (min(total.stats$var_from)+20000), y = 25, label = c("<0.1"), hjust = 0.5,size=4) +
  annotate("text", x = 149625104, y = 38, size=5, label = c("cg02724909"), hjust = 0.5)
ra.plot

#2nd plot out genes
gene.ensembl <- useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl", GRCh = 37)
sel.chr=5
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
  geom_segment(data = genes[seq(1,63,3),],
               aes(x=start, xend=end, y=order+1, yend=order+1), size = 1, colour = cols3[2],
               arrow = arrow(length = unit(0.25, "cm"))) +
  geom_segment(data = genes[seq(2,63,3),],
               aes(x=start, xend=end, y=order+1.1, yend=order+1.1), size = 1, colour = cols3[2],
               arrow = arrow(length = unit(0.25, "cm"))) +
  geom_segment(data = genes[seq(3,63,3),],
               aes(x=start, xend=end, y=order+1.2, yend=order+1.2), size = 1, colour = cols3[2],
               arrow = arrow(length = unit(0.25, "cm"))) +
  geom_text_repel( data = genes[seq(1,63,3),], aes(x=mid,y=order+1, label=external_gene_name), size=2, col = c("black"),
                   nudge_y =-0.05, segment.color = NA) +
  geom_text_repel( data = genes[seq(2,63,3),], aes(x=mid,y=order+1.1, label=external_gene_name), size=2, col = c("black"),
                   nudge_y =-0.05, segment.color = NA) +
  geom_text_repel( data = genes[seq(3,63,3),], aes(x=mid,y=order+1.2, label=external_gene_name), size=2, col = c("black"),
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
recomb <- read.table("genetic_map_chr5_combined_b37.txt", header = T)
recomb <- subset(recomb, position>min(total.stats$var_from) & position<max(total.stats$var_from))

recomb_rate <- ggplot(recomb, aes(x=position, y=COMBINED_rate.cM.Mb.)) + 
  geom_line() +
  theme_bw() +
  ylab("cM/Mb") +
  xlab("chromosome 5") +
  scale_x_continuous(breaks=c(149400000,149500000,149600000,149700000,149800000),
                     labels=c("149.4Mb", "149.5Mb", "149.6Mb", "149.7Mb", "149.8Mb")) +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"))

#add association, genes and recombination plots together
p1 <- (ra.plot/genes.plot/recomb_rate) + plot_layout(heights = c(3, 1.5, 0.5))

ggplot2::ggsave(
  "cg02724909_de_novo_RA.jpg",
  width = 7,
  height = 7,
  dpi = 300
) 

#example of a baseline enhanced LPS-modulated mQTL at cg17462560
#read in genotypes and untreated, LPS-treated and delta methylation (LPS-treated minus untreated) for cg17462560
total <- read.table("cg17462560.txt", header = T)

p_labels = data.frame(expt = c("base"), label = c("italic(P)==1.06%*%10^-19"))

cols <- brewer.pal(9,"Set1")

#plot out box plot of per genotype delta methylation at cg17462560
total$state <- "delta"
box.delta.cg17462560 = ggplot(total, aes(x=factor(round(rs6875879,0)), y=delta)) +
  geom_dotplot(binaxis="y", binwidth=0.005, stackdir="center", alpha = 0.72) + geom_boxplot(alpha = 0.2) +
  scale_x_discrete(breaks=c(0,1,2), labels = c("GG", "GA", "AA")) + aes(fill = factor(round(rs6875879,0)), col=factor(round(rs6875879,0))) + scale_fill_manual(values = cols[c(4,4,4)]) + scale_colour_manual(values = cols[c(4,4,4)]) +
  ylab("cg17462560") + xlab("rs6875879") +
  theme_bw() + theme(legend.position = "none", axis.text=element_text(size=12),
                     axis.title=element_text(size=15), strip.background = element_rect(color="black", fill=cols[9], size=1.5, linetype="solid"), strip.text.x = element_text(
                       size = 15, color = "white", face = "bold")) + ylim(NA, 0.1) +
  facet_wrap( ~ state, ncol=1) +
  geom_text(x=2, y=0.075, aes(label=label), data=p_labels, parse=TRUE, inherit.aes=F, size = 5)
box.delta.cg17462560

#plot out box plot of per genotype baseline and LPS-treated methylation at cg17462560
box.facet <- data.frame(rbind(cbind(total$rs6875879, total$ut, "UT"),
                              cbind(total$rs6875879, total$lps, "LPS")))
colnames(box.facet) <- c("rs6875879", "probe", "state")
box.facet$rs6875879 <- as.numeric(box.facet$rs6875879)
box.facet$probe <- as.numeric(box.facet$probe)

box.facet$state <- ordered(box.facet$state, levels = c("UT", "LPS"))
p_labels = data.frame(state = c("UT", "LPS"), label = c("italic(P)==5.91%*%10^-5", "italic(P)==2.13%*%10^-24"))
p_labels$state <- ordered(p_labels$state, levels = c("UT", "LPS"))

facet.cg17462560.box <- ggplot(box.facet, aes(x=factor(round(rs6875879,0)), y=probe)) +
  geom_dotplot(binaxis="y", binwidth=0.003, stackdir="center", alpha = 0.72) + geom_boxplot(alpha = 0.2) +
  scale_x_discrete(breaks=c(0,1,2), labels = c("GG", "GA", "AA")) + aes(fill = state, col=state) + scale_fill_manual(values = cols[c(3,1)]) + scale_colour_manual(values = cols[c(3,1)]) +
  ylab("cg17462560") + xlab("rs6875879") +
  theme_bw() + theme(legend.position = "none", axis.text=element_text(size=12),
                     axis.title=element_text(size=15)) + 
  facet_wrap( ~ state, ncol=2) +
  ylim(NA, 1.0) +
  geom_text(x=2, y=0.95, aes(label=label), data=p_labels, parse=TRUE, inherit.aes=F, size = 5) +
  theme(strip.background = element_rect(color="black", fill=cols[9], size=1.5, linetype="solid"), strip.text.x = element_text(
    size = 15, color = "white", face = "bold"))
facet.cg17462560.box

p1 <- (facet.cg17462560.box|box.delta.cg17462560) + plot_layout(widths = c(2, 1))

ggplot2::ggsave(
  "cg17462560_baseline_enhanced_box.jpg",
  width = 7,
  height = 7,
  dpi = 300
) 

#RA plot for LPS-modulated mQTL at cg17462560
#read in summary stats and r2 to rs6875879 for region (+/-250kb)
total.stats <- read.table("cg17462560_region.txt", header = T)

#split r2 values into bins
total.stats$bin_r2 <- 1
total.stats$bin_r2[which(total.stats$r2>0.1 & total.stats$r2 <= 0.3)] <- 2
total.stats$bin_r2[which(total.stats$r2>0.3 & total.stats$r2 <= 0.5)] <- 3
total.stats$bin_r2[which(total.stats$r2>0.5 & total.stats$r2 <= 0.8)] <- 4
total.stats$bin_r2[which(total.stats$r2>0.8)] <- 5

#RA plot
#set colours
cols3 <- brewer.pal(8,"Paired")
cols2 <- brewer.pal(11,"Spectral")
cols <- brewer.pal(9,"Set1")

#1st plot regional association
total.stats$annotate <- 0
total.stats$annotate[which(total.stats$rsid=="rs6875879")] <- 0

ra.plot <- ggplot(total.stats, aes(x=var_from, y=-log10(cpg.p))) + 
  xlim(min(total.stats$var_from), max(total.stats$var_from)) +
  annotate("rect", xmin = 77814503, xmax = 77819503, ymin = 0, ymax = 26,
           alpha = 1,fill = cols[4]) +
  geom_point(data=subset(total.stats, bin_r2==1), color=cols[9], size=3) + 
  geom_point(data=subset(total.stats, bin_r2==2), color=cols2[5], size=3) + 
  geom_point(data=subset(total.stats, bin_r2==3), color=cols2[3], size=3) + 
  geom_point(data=subset(total.stats, bin_r2==4), color=cols2[2], size=3) +
  geom_point(data=subset(total.stats, bin_r2==5), color=cols2[1], size=3) +
  scale_y_continuous(name="-log P-value", breaks=c(0,5,10,15,20,25)) +
  xlab(NULL) + 
  theme_bw() + theme(legend.position = "none", axis.text=element_text(size=12),
                     axis.title=element_text(size=12)) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  geom_text_repel( data=subset(total.stats, annotate==1), aes(label=rsid), size=5, col = c("black"), nudge_x = 50000) +
  annotate("text", x = (min(total.stats$var_from)+20000), y = 25, label = quote(r^2), hjust = 0.5,size=4) +
  annotate("point", x = (min(total.stats$var_from)+60000), y = 23.5, size = 4, colour = cols2[1]) +
  annotate("text", x = (min(total.stats$var_from)+20000), y = 23.5, label = c("0.8-1.0"), hjust = 0.5,size=4) +
  annotate("point", x = (min(total.stats$var_from)+60000), y =22, size = 4, colour = cols2[2]) +
  annotate("text", x = (min(total.stats$var_from)+20000), y = 22, label = c("0.5-0.8"), hjust = 0.5,size=4) +
  annotate("point", x = (min(total.stats$var_from)+60000), y = 20.5, size = 4, colour = cols2[3]) +
  annotate("text", x = (min(total.stats$var_from)+20000), y = 20.5, label = c("0.3-0.5"), hjust = 0.5,size=4) +
  annotate("point", x = (min(total.stats$var_from)+60000), y = 19, size = 4, colour = cols2[5]) +
  annotate("text", x = (min(total.stats$var_from)+20000), y = 19, label = c("0.1-0.3"), hjust = 0.5,size=4) +
  annotate("point", x = (min(total.stats$var_from)+60000), y = 17.5, size = 4, colour = cols[9]) +
  annotate("text", x = (min(total.stats$var_from)+20000), y = 17.5, label = c("<0.1"), hjust = 0.5,size=4) +
  annotate("text", x = 77814503, y = 27, size=5, label = c("cg02724909"), hjust = 0.5)
ra.plot

#2nd plot out location of genes in region
#plot genes
gene.ensembl <- useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl", GRCh = 37)
sel.chr=5
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
  geom_segment(data = genes[seq(1,13,3),],
               aes(x=start, xend=end, y=order+1, yend=order+1), size = 1, colour = cols3[2],
               arrow = arrow(length = unit(0.25, "cm"))) +
  geom_segment(data = genes[seq(2,13,3),],
               aes(x=start, xend=end, y=order+1.1, yend=order+1.1), size = 1, colour = cols3[2],
               arrow = arrow(length = unit(0.25, "cm"))) +
  geom_segment(data = genes[seq(3,13,3),],
               aes(x=start, xend=end, y=order+1.2, yend=order+1.2), size = 1, colour = cols3[2],
               arrow = arrow(length = unit(0.25, "cm"))) +
  geom_text_repel( data = genes[seq(1,13,3),], aes(x=mid,y=order+1, label=external_gene_name), size=2, col = c("black"),
                   nudge_y =-0.05, segment.color = NA) +
  geom_text_repel( data = genes[seq(2,13,3),], aes(x=mid,y=order+1.1, label=external_gene_name), size=2, col = c("black"),
                   nudge_y =-0.05, segment.color = NA) +
  geom_text_repel( data = genes[seq(3,13,3),], aes(x=mid,y=order+1.2, label=external_gene_name), size=2, col = c("black"),
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
recomb <- read.table("genetic_map_chr5_combined_b37.txt", header = T)
recomb <- subset(recomb, position>min(total.stats$var_from) & position<max(total.stats$var_from))

recomb_rate <- ggplot(recomb, aes(x=position, y=COMBINED_rate.cM.Mb.)) + 
  geom_line() +
  theme_bw() +
  ylab("cM/Mb") +
  xlab("chromosome 5") +
  scale_x_continuous(breaks=c(77600000,77700000,77800000,77900000,78000000),
                     labels=c("77.6Mb", "77.7Mb", "77.8Mb", "77.9Mb", "78.0Mb")) +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"))

#add association, genes and recombination plots together
p1 <- (ra.plot/genes.plot/recomb_rate) + plot_layout(heights = c(3, 1.5, 0.5))
ggplot2::ggsave(
  "cg17462560_baseline_enhanced_RA.jpg",
  width = 7,
  height = 7,
  dpi = 300
) 

