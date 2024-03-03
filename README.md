# Characterisation of the genetic determinants of context specific DNA methylation in primary monocytes

**Abstract**
DNA methylation (DNAm) has pervasive effects on gene expression and associations with ageing-related traits. Here we describe monocyte DNAm responses to inflammatory stimuli across 190 individuals. We find that, unlike the similarly widespread changes in gene expression elicited by LPS and IFN$\gamma$, DNAm is markedly more sensitive to LPS. Exposure to LPS caused differential methylation at 7,359 high confidence immune-modulated CpGs (imCpGs) which differed in genomic localisation and transcription factor usage according to methylation loss or gain. Demethylated imCpGs are profoundly enriched for enhancers, and are over-represented by genes implicated in human diseases, most notably cancer. We find LPS-induced demethylation follows hydroxymethylation and for most sites the degree of demethylation correlates with baseline signal. Notably, we find 24 hours of LPS exposure triggers gain in epigenetic age by approximately 6 months, directly linking the adverse health associated phenomena of accelerated epigenetic ageing with innate immune activity. Finally, we explore the effect of genetic variation on LPS-induced changes in DNAm in 188 individuals, identifying 234 imCpGs under genetic control. Exploring shared causal loci between LPS-induced DNAm responses and human disease traits highlights examples of human disease associated loci that also modulate imCpG formation. In summary, our findings demonstrate innate immune activity remodels DNAm in a highly punctate, enhancer-enriched fashion that is influenced by genetic variation and predominantly involves genes frequently mutated in cancer.

[Preprint](https://www.biorxiv.org/content/10.1101/2023.05.17.541041v1)

**Overview of repository**
* LPS-induced differential methylation in monocytes: contains R script and source data to reproduce analysis and figures exploring changes in DNAm in monocytes in response to LPS stimulation. Includes; control of inflation and bias using bacon [Maarten van Itersen, *et al*. "Controlling bias and inflation in epigenome- and transcriptome-wide association studies using the empirical null distribution." *Genome Biology.* 2017.](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-016-1131-9), enrichment analyses exploring immune-modulated CpG genomic context, transcription factor binding site overlap, and local genic enrichment for biological and disease pathways.
* LPS-induced epigenetic age acceleration: contains R script and source data to reproduce analysis and figures describing changes in epigenetic age in monocytes following LPS stimulation. Epigenetic age estimates are calculated from normalised methylation data using ENmix [Zongli Xu, *et al*. "The ENmix DNA methylation analysis pipeline for Illumina BeadChip and comparisons with seven other preprocessing pipelines" *Clinical Epigenetics.* 2021.](https://clinicalepigeneticsjournal.biomedcentral.com/articles/10.1186/s13148-021-01207-1) and dnaMethyAge [Yuchen Wang, *et al*. "Insights into ageing rates comparison across tissues from recalibrating cerebellum DNA methylation clock" *Geroscience.* 2024.](https://pubmed.ncbi.nlm.nih.gov/37597113).
* Baseline correlated and independent im-CpGs: contains R script and source data to reproduce analysis and figures describing correlation of LPS-induced changes in CpG methylation with the baseline methylation state.
* LPS-induced mQTL mapping: contains R script and source data to reproduce analysis and figures describing eQTL mapping of LPS-induced methylation changes in monocytes. eQTL mapping performed in QTLtools: [Olivier Delaneau, *et al*. "A complete tool set for molecular QTL discovery and analysis." *Nature Communications.* 2017.](https://www.nature.com/articles/ncomms15452)
* LPS-induced mQTL colocalisation: contains R script and source data to reproduce analysis and figures describing an LPS-induced mQTL, which colocalises with GWAS signals for haematological traits and RNA expression eQTL in monocytes. Colocalisation analysis was performed in coloc: [Claudia Giambartolomei, *et al*. "Bayesian Test for Colocalisation between Pairs of Genetic Association Studies Using Summary Statistics." *PLOS Genetics.* 2014.](https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1004383). GWAS summary statistics for haematological traits in the UK Biobank are taken from [Neale lab GWAS round 2](http://www.nealelab.is/uk-biobank/). eQTL summary statistics of gene expression in untreated and LPS-stimulation monocytes are taken from: [Benjamin Fairfax, *et al*. "Innate immune activity conditions the effect of regulatory variants upon monocyte gene expression." *Science.* 2014.](https://pubmed.ncbi.nlm.nih.gov/24604202/).
* ZFP57 and LPS-induced methylation: contains R script and source data to reproduce analysis and figures describing the effect of an eQTL for ZPF57 on local and distal methylation patterns in monocytes, which result in differential sensivitity of CpGs to LPS stimulation.
