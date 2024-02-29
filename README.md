# Characterisation of the genetic determinants of context specific DNA methylation in primary monocytes

**Abstract**
DNA methylation (DNAm) has pervasive effects on gene expression and associations with ageing-related traits. Here we describe monocyte DNAm responses to inflammatory stimuli across 190 individuals. We find that, unlike the similarly widespread changes in gene expression elicited by LPS and IFNγ, DNAm is markedly more sensitive to LPS. Exposure to LPS caused differential methylation at 20,858 immune modulated CpGs (imCpGs) which display distinct genomic localisation and transcription factor usage, dependent upon whether methylation is lost or gained. Demethylated imCpGs are profoundly enriched for enhancers, and are over-represented by genes implicated in human dis eases, most notably cancer. We find LPS-induced demethylation follows hydroxymethylation and for most sites the degree of demethylation correlates with baseline signal. Notably, we find LPS exposure triggers gain in epigenetic age by approximately 6 months, identifying a potential cause of accelerated epigentic aging which has diverse negative health associations. Finally, we explore the effect of genetic variation on LPS induced changes in DNAm in 188 individuals, identifying 234 imCpGs under genetic control. Exploring shared causal loci between LPS induced DNAm responses and human disease traits highlights examples of human disease associated loci that also modulate imCpG formation. In summary, our findings suggest innate immune activity continually remodels DNAm in a highly punctate, enhancer enriched fashion that is under tight genetic control and predominantly involves genes commonly mutated in cancer.

[Preprint](https://www.biorxiv.org/content/10.1101/2023.05.17.541041v1)

**Overview of repository**
* LPS-induced differential methylation: contains R script and source data to reproduce analysis and figures exploring changes in DNAm in monocytes in response to LPS stimulation. Includes; control of inflation and bias using bacon [Maarten van Itersen, *et al*. "Controlling bias and inflation in epigenome- and transcriptome-wide association studies using the empirical null distribution." *Genome Biology.* 2017.](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-016-1131-9), enrichment analyses exploring immune-modulated CpG genomic context, transcription factor binding site overlap, and local genic enrichment for biological and disease pathways.
* LPS-induced epigenetic age acceleration: contains R script and source data to reproduce analysis and figures describing changes in epigenetic age in monocytes following LPS stimulation. Epigenetic age estimates are calculated from normalised methylation data using ENmix [Zongli Xu, *et al*. "The ENmix DNA methylation analysis pipeline for Illumina BeadChip and comparisons with seven other preprocessing pipelines" *Clinical Epigenetics.* 2021.](https://clinicalepigeneticsjournal.biomedcentral.com/articles/10.1186/s13148-021-01207-1) and dnaMethyAge [Yuchen Wang, *et al*. "Insights into ageing rates comparison across tissues from recalibrating cerebellum DNA methylation clock" *Geroscience.* 2024.](https://pubmed.ncbi.nlm.nih.gov/37597113).
