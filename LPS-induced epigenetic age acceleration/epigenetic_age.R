rm(list = ls(all=TRUE))

library(ggplot2)
library(RColorBrewer)
library(dplyr)
library(ggrepel)
library(patchwork)

cols <- brewer.pal(8,"Set2")
cols2 <- brewer.pal(8,"Dark2")
cols3 <- brewer.pal(8,"Paired")


#read in sample info with epigenetic age estimates
age_meth <- read.table("age_clocks.txt", header = T)

######epigenetic age estimates calculated with ENlist and dnaMethyage

#calculate differential epigenetic/chronological age
age_meth$delta_horvath<-age_meth$dnam_horvath-age_meth$Age
age_meth$delta_hannum<-age_meth$dnam_hannum-age_meth$Age
age_meth$delta_phenoage<-age_meth$dnam_phenoage-age_meth$Age
age_meth$delta_PCgrimAge<-age_meth$dnam_PCgrimAge-age_meth$Age

#dichotomise samples into high and low baseline age acceleration
increasedId_horvath<-age_meth[age_meth$delta_horvath>median(age_meth$delta_horvath)&age_meth$Treatment=="UT",c("Pool_ID")]
reducedId_horvath<-age_meth[age_meth$delta_horvath<median(age_meth$delta_horvath)&age_meth$Treatment=="UT",c("Pool_ID")]

increasedId_hannum<-age_meth[age_meth$delta_hannum>median(age_meth$delta_hannum)&age_meth$Treatment=="UT",c("Pool_ID")]
reducedId_hannum<-age_meth[age_meth$delta_hannum<median(age_meth$delta_hannum)&age_meth$Treatment=="UT",c("Pool_ID")]

increasedId_phenoage<-age_meth[age_meth$delta_phenoage>median(age_meth$delta_phenoage)&age_meth$Treatment=="UT",c("Pool_ID")]
reducedId_phenoage<-age_meth[age_meth$delta_phenoage<median(age_meth$delta_phenoage)&age_meth$Treatment=="UT",c("Pool_ID")]

increasedId_PCgrimAge<-age_meth[age_meth$delta_PCgrimAge>median(age_meth$delta_PCgrimAge)&age_meth$Treatment=="UT",c("Pool_ID")]
reducedId_PCgrimAge<-age_meth[age_meth$delta_PCgrimAge<median(age_meth$delta_PCgrimAge)&age_meth$Treatment=="UT",c("Pool_ID")]

#calculate LPS-induced age-acceleration for Horvath clock

panCohort_horvath<-wilcox.test(
  age_meth[age_meth$Treatment=="LPS24",c("dnam_horvath")],
  age_meth[age_meth$Treatment=="UT",c("dnam_horvath")],
  paired = TRUE)
panCohort_horvath

median(age_meth[age_meth$Treatment=="LPS24",c("dnam_horvath")]-age_meth[age_meth$Treatment=="UT",c("dnam_horvath")])

high_baseline_horvath<-wilcox.test(
  age_meth[age_meth$Pool_ID%in%increasedId_horvath&age_meth$Treatment=="LPS24",c("dnam_horvath")],
  age_meth[age_meth$Pool_ID%in%increasedId_horvath&age_meth$Treatment=="UT",c("dnam_horvath")],
  paired = TRUE)
high_baseline_horvath

median(age_meth[age_meth$Pool_ID%in%increasedId_horvath&age_meth$Treatment=="LPS24",c("dnam_horvath")]-age_meth[age_meth$Pool_ID%in%increasedId_horvath&age_meth$Treatment=="UT",c("dnam_horvath")])

low_baseline_horvath<-wilcox.test(
  age_meth[age_meth$Pool_ID%in%reducedId_horvath&age_meth$Treatment=="LPS24",c("dnam_horvath")],
  age_meth[age_meth$Pool_ID%in%reducedId_horvath&age_meth$Treatment=="UT",c("dnam_horvath")],
  paired = TRUE)
low_baseline_horvath

median(age_meth[age_meth$Pool_ID%in%reducedId_horvath&age_meth$Treatment=="LPS24",c("dnam_horvath")]-age_meth[age_meth$Pool_ID%in%reducedId_horvath&age_meth$Treatment=="UT",c("dnam_horvath")])

#calculate LPS-induced age-acceleration for Hannum clock

panCohort_hannum<-wilcox.test(
  age_meth[age_meth$Treatment=="LPS24",c("dnam_hannum")],
  age_meth[age_meth$Treatment=="UT",c("dnam_hannum")],
  paired = TRUE)
panCohort_hannum

median(age_meth[age_meth$Treatment=="LPS24",c("dnam_hannum")]-age_meth[age_meth$Treatment=="UT",c("dnam_hannum")])

high_baseline_hannum<-wilcox.test(
  age_meth[age_meth$Pool_ID%in%increasedId_hannum&age_meth$Treatment=="LPS24",c("dnam_hannum")],
  age_meth[age_meth$Pool_ID%in%increasedId_hannum&age_meth$Treatment=="UT",c("dnam_hannum")],
  paired = TRUE)
high_baseline_hannum

median(age_meth[age_meth$Pool_ID%in%increasedId_hannum&age_meth$Treatment=="LPS24",c("dnam_hannum")]-age_meth[age_meth$Pool_ID%in%increasedId_hannum&age_meth$Treatment=="UT",c("dnam_hannum")])

low_baseline_hannum<-wilcox.test(
  age_meth[age_meth$Pool_ID%in%reducedId_hannum&age_meth$Treatment=="LPS24",c("dnam_hannum")],
  age_meth[age_meth$Pool_ID%in%reducedId_hannum&age_meth$Treatment=="UT",c("dnam_hannum")],
  paired = TRUE)
low_baseline_hannum

median(age_meth[age_meth$Pool_ID%in%reducedId_hannum&age_meth$Treatment=="LPS24",c("dnam_hannum")]-age_meth[age_meth$Pool_ID%in%reducedId_hannum&age_meth$Treatment=="UT",c("dnam_hannum")])

#calculate LPS-induced age-acceleration for PhenoAge clock

panCohort_phenoage<-wilcox.test(
  age_meth[age_meth$Treatment=="LPS24",c("dnam_phenoage")],
  age_meth[age_meth$Treatment=="UT",c("dnam_phenoage")],
  paired = TRUE)
panCohort_phenoage

median(age_meth[age_meth$Treatment=="LPS24",c("dnam_phenoage")]-age_meth[age_meth$Treatment=="UT",c("dnam_phenoage")])

high_baseline_phenoage<-wilcox.test(
  age_meth[age_meth$Pool_ID%in%increasedId_phenoage&age_meth$Treatment=="LPS24",c("dnam_phenoage")],
  age_meth[age_meth$Pool_ID%in%increasedId_phenoage&age_meth$Treatment=="UT",c("dnam_phenoage")],
  paired = TRUE)
high_baseline_phenoage

median(age_meth[age_meth$Pool_ID%in%increasedId_phenoage&age_meth$Treatment=="LPS24",c("dnam_phenoage")]-age_meth[age_meth$Pool_ID%in%increasedId_phenoage&age_meth$Treatment=="UT",c("dnam_phenoage")])

low_baseline_phenoage<-wilcox.test(
  age_meth[age_meth$Pool_ID%in%reducedId_phenoage&age_meth$Treatment=="LPS24",c("dnam_phenoage")],
  age_meth[age_meth$Pool_ID%in%reducedId_phenoage&age_meth$Treatment=="UT",c("dnam_phenoage")],
  paired = TRUE)
low_baseline_phenoage

median(age_meth[age_meth$Pool_ID%in%reducedId_phenoage&age_meth$Treatment=="LPS24",c("dnam_phenoage")]-age_meth[age_meth$Pool_ID%in%reducedId_phenoage&age_meth$Treatment=="UT",c("dnam_phenoage")])

#calculate LPS-induced age-acceleration for PCGrimAge clock

panCohort_PCgrimAge<-wilcox.test(
  age_meth[age_meth$Treatment=="LPS24",c("dnam_PCgrimAge")],
  age_meth[age_meth$Treatment=="UT",c("dnam_PCgrimAge")],
  paired = TRUE)
panCohort_PCgrimAge

median(age_meth[age_meth$Treatment=="LPS24",c("dnam_PCgrimAge")]-age_meth[age_meth$Treatment=="UT",c("dnam_PCgrimAge")])

high_baseline_PCgrimAge<-wilcox.test(
  age_meth[age_meth$Pool_ID%in%increasedId_PCgrimAge&age_meth$Treatment=="LPS24",c("dnam_PCgrimAge")],
  age_meth[age_meth$Pool_ID%in%increasedId_PCgrimAge&age_meth$Treatment=="UT",c("dnam_PCgrimAge")],
  paired = TRUE)
high_baseline_PCgrimAge

median(age_meth[age_meth$Pool_ID%in%increasedId_PCgrimAge&age_meth$Treatment=="LPS24",c("dnam_PCgrimAge")]-age_meth[age_meth$Pool_ID%in%increasedId_PCgrimAge&age_meth$Treatment=="UT",c("dnam_PCgrimAge")])

low_baseline_PCgrimAge<-wilcox.test(
  age_meth[age_meth$Pool_ID%in%reducedId_PCgrimAge&age_meth$Treatment=="LPS24",c("dnam_PCgrimAge")],
  age_meth[age_meth$Pool_ID%in%reducedId_PCgrimAge&age_meth$Treatment=="UT",c("dnam_PCgrimAge")],
  paired = TRUE)
low_baseline_PCgrimAge

median(age_meth[age_meth$Pool_ID%in%reducedId_PCgrimAge&age_meth$Treatment=="LPS24",c("dnam_PCgrimAge")]-age_meth[age_meth$Pool_ID%in%reducedId_PCgrimAge&age_meth$Treatment=="UT",c("dnam_PCgrimAge")])

#plot LPS-induced age acceleration for Horvath clock - overall
age_meth$Treatment <- factor(age_meth$Treatment, levels = c("UT", "LPS24"))

ann_text.horvath<-data.frame(dnam_horvath=c(70),Treatment=1.5,p.value=c(panCohort_horvath[[3]]))

plot_Horvath<-ggplot(age_meth, aes(Treatment, dnam_horvath))+
  geom_violin(aes(fill=Treatment))+
  geom_boxplot(colour="grey45",alpha=0.25,position=position_dodge(1),width=0.64)+
  geom_text(aes(label=paste0("P=",signif(p.value,2))), data = ann_text.horvath) +
  scale_y_continuous("methylation age (yrs)")+
  scale_x_discrete("Treatment")+
  ggtitle("Horvath") +
  theme_bw() +
  theme(plot.title = element_text(size = 25),strip.background = element_rect(color="black", fill="dark grey", size=1.5, linetype="solid"), strip.text.x = element_text(
    size = 15, color = "white", face = "bold"))
plot_Horvath

#plot LPS-induced age acceleration for Horvath clock - facetted by baseline age acceleration
ann_text.horvath_facet<-data.frame(high_low_horvath=c("high", "low"), dnam_horvath=c(70),Treatment=1.5,p.value=c(high_baseline_horvath[[3]], low_baseline_horvath[[3]]))

plot_Horvath.facet<-ggplot(age_meth, aes(Treatment, dnam_horvath))+
  geom_violin(aes(fill=Treatment))+
  geom_boxplot(colour="grey45",alpha=0.25,position=position_dodge(1),width=0.64)+
  geom_text(aes(label=paste0("P=",signif(p.value,2))), data = ann_text.horvath_facet) +
  scale_y_continuous("methylation age (yrs)")+
  scale_x_discrete("Treatment")+
  facet_wrap(~high_low_horvath) +
  theme_bw() +
  theme(strip.background = element_rect(color="black", fill="dark grey", size=1.5, linetype="solid"), strip.text.x = element_text(
    size = 15, color = "white", face = "bold"))
plot_Horvath.facet



#plots - Hannum

panCohort<-wilcox.test(
  age_meth[age_meth$Sample_Group=="LPS24",c("dnam_hannum")],
  age_meth[age_meth$Sample_Group=="UT",c("dnam_hannum")],
  paired = TRUE)

age_meth$Sample_Group <- factor(age_meth$Sample_Group, levels = c("UT", "LPS24"))
ann_text.1<-data.frame(dnam_hannum=c(70),Sample_Group=1.5,p.value=c(panCohort[[3]]))

plot_Hannum<-ggplot(age_meth, aes(Sample_Group, dnam_hannum))+
  geom_violin(aes(fill=Sample_Group))+
  geom_boxplot(colour="grey45",alpha=0.25,position=position_dodge(1),width=0.64)+
  geom_text(aes(label=paste0("P=",signif(p.value,2))), data = ann_text.1) +
  scale_y_continuous("methylation age (yrs)")+
  scale_x_discrete("Treatment")+
  ggtitle("Hannum") +
  theme_bw() +
  theme(plot.title = element_text(size = 25),strip.background = element_rect(color="black", fill="dark grey", size=1.5, linetype="solid"), strip.text.x = element_text(
    size = 15, color = "white", face = "bold"))
plot_Hannum

age_meth$deltaAge<-age_meth$dnam_hannum-age_meth$Age
age_meth$baseline<-ifelse(age_meth$Age<age_meth$dnam_hannum,"reduced","increased")

high_acc<-age_meth[age_meth$deltaAge>median(age_meth$deltaAge)&age_meth$Sample_Group=="UT",c("Pool_ID")]
low_acc<-age_meth[age_meth$deltaAge<median(age_meth$deltaAge)&age_meth$Sample_Group=="UT",c("Pool_ID")]

high_baseline<-wilcox.test(
  age_meth[age_meth$Pool_ID%in%high_acc&age_meth$Sample_Group=="LPS24",c("dnam_hannum")],
  age_meth[age_meth$Pool_ID%in%high_acc&age_meth$Sample_Group=="UT",c("dnam_hannum")],
  paired = TRUE)

low_baseline<-wilcox.test(
  age_meth[age_meth$Pool_ID%in%low_acc&age_meth$Sample_Group=="LPS24",c("dnam_hannum")],
  age_meth[age_meth$Pool_ID%in%low_acc&age_meth$Sample_Group=="UT",c("dnam_hannum")],
  paired = TRUE)

age_meth$high_low_hannum <- "high"
age_meth$high_low_hannum[which(age_meth$Pool_ID %in% low_acc)] <- "low"

age_meth$high_low_hannum <- factor(age_meth$high_low_hannum, levels = c("high", "low"))

age_meth$Sample_Group <- factor(age_meth$Sample_Group, levels = c("UT", "LPS24"))
ann_text.1<-data.frame(high_low_hannum=c("high", "low"), dnam_hannum=c(70),Sample_Group=1.5,p.value=c(high_baseline[[3]], low_baseline[[3]]))

plot_Hannum.facet<-ggplot(age_meth, aes(Sample_Group, dnam_hannum))+
  geom_violin(aes(fill=Sample_Group))+
  geom_boxplot(colour="grey45",alpha=0.25,position=position_dodge(1),width=0.64)+
  geom_text(aes(label=paste0("P=",signif(p.value,2))), data = ann_text.1) +
  scale_y_continuous("methylation age (yrs)")+
  scale_x_discrete("Treatment")+
  facet_wrap(~high_low_hannum) +
  theme_bw() +
  theme(strip.background = element_rect(color="black", fill="dark grey", size=1.5, linetype="solid"), strip.text.x = element_text(
    size = 15, color = "white", face = "bold"))
plot_Hannum.facet



#plots - PhenoAge
str(age_meth)

panCohort<-wilcox.test(
  age_meth[age_meth$Sample_Group=="LPS24",c("dnam_phenoage")],
  age_meth[age_meth$Sample_Group=="UT",c("dnam_phenoage")],
  paired = TRUE)

age_meth$Sample_Group <- factor(age_meth$Sample_Group, levels = c("UT", "LPS24"))
ann_text.1<-data.frame(dnam_phenoage=c(75),Sample_Group=1.5,p.value=c(panCohort[[3]]))

plot_phenoage<-ggplot(age_meth, aes(Sample_Group, dnam_phenoage))+
  geom_violin(aes(fill=Sample_Group))+
  geom_boxplot(colour="grey45",alpha=0.25,position=position_dodge(1),width=0.64)+
  geom_text(aes(label=paste0("P=",signif(p.value,2))), data = ann_text.1) +
  scale_y_continuous("methylation age (yrs)")+
  scale_x_discrete("Treatment")+
  ggtitle("PhenoAge") +
  theme_bw() +
  theme(plot.title = element_text(size = 25), strip.background = element_rect(color="black", fill="dark grey", size=1.5, linetype="solid"), strip.text.x = element_text(
    size = 15, color = "white", face = "bold"))
plot_phenoage


age_meth$deltaAge<-age_meth$dnam_phenoage-age_meth$Age
age_meth$baseline<-ifelse(age_meth$Age<age_meth$dnam_phenoage,"reduced","increased")

high_acc<-age_meth[age_meth$deltaAge>median(age_meth$deltaAge)&age_meth$Sample_Group=="UT",c("Pool_ID")]
low_acc<-age_meth[age_meth$deltaAge<median(age_meth$deltaAge)&age_meth$Sample_Group=="UT",c("Pool_ID")]

high_baseline<-wilcox.test(
  age_meth[age_meth$Pool_ID%in%high_acc&age_meth$Sample_Group=="LPS24",c("dnam_phenoage")],
  age_meth[age_meth$Pool_ID%in%high_acc&age_meth$Sample_Group=="UT",c("dnam_phenoage")],
  paired = TRUE)

low_baseline<-wilcox.test(
  age_meth[age_meth$Pool_ID%in%low_acc&age_meth$Sample_Group=="LPS24",c("dnam_phenoage")],
  age_meth[age_meth$Pool_ID%in%low_acc&age_meth$Sample_Group=="UT",c("dnam_phenoage")],
  paired = TRUE)

age_meth$high_low_phenoage <- "high"
age_meth$high_low_phenoage[which(age_meth$Pool_ID %in% low_acc)] <- "low"

age_meth$high_low_phenoage <- factor(age_meth$high_low_phenoage, levels = c("high", "low"))

age_meth$Sample_Group <- factor(age_meth$Sample_Group, levels = c("UT", "LPS24"))
ann_text.1<-data.frame(high_low_phenoage=c("high", "low"), dnam_phenoage=c(80),Sample_Group=1.5,p.value=c(high_baseline[[3]], low_baseline[[3]]))

plot_phenoage.facet<-ggplot(age_meth, aes(Sample_Group, dnam_phenoage))+
  geom_violin(aes(fill=Sample_Group))+
  geom_boxplot(colour="grey45",alpha=0.25,position=position_dodge(1),width=0.64)+
  geom_text(aes(label=paste0("P=",signif(p.value,2))), data = ann_text.1) +
  scale_y_continuous("methylation age (yrs)")+
  scale_x_discrete("Treatment")+
  facet_wrap(~high_low_phenoage) +
  theme_bw() +
  theme(strip.background = element_rect(color="black", fill="dark grey", size=1.5, linetype="solid"), strip.text.x = element_text(
    size = 15, color = "white", face = "bold"))
plot_phenoage.facet



#plots - PCGrimAge
str(age_meth)

panCohort<-wilcox.test(
  age_meth[age_meth$Sample_Group=="LPS24",c("PCGrimAge")],
  age_meth[age_meth$Sample_Group=="UT",c("PCGrimAge")],
  paired = TRUE)

age_meth$Sample_Group <- factor(age_meth$Sample_Group, levels = c("UT", "LPS24"))
ann_text.1<-data.frame(PCGrimAge=c(80),Sample_Group=1.5,p.value=c(panCohort[[3]]))

plot_grimage<-ggplot(age_meth, aes(Sample_Group, PCGrimAge))+
  geom_violin(aes(fill=Sample_Group))+
  geom_boxplot(colour="grey45",alpha=0.25,position=position_dodge(1),width=0.64)+
  geom_text(aes(label=paste0("P=",signif(p.value,2))), data = ann_text.1) +
  scale_y_continuous("methylation age (yrs)")+
  scale_x_discrete("Treatment")+
  ggtitle("PCGrimAge") +
  theme_bw() +
  theme(plot.title = element_text(size = 25), strip.background = element_rect(color="black", fill="dark grey", size=1.5, linetype="solid"), strip.text.x = element_text(
    size = 15, color = "white", face = "bold"))
plot_grimage


age_meth$deltaAge<-age_meth$PCGrimAge-age_meth$Age
age_meth$baseline<-ifelse(age_meth$Age<age_meth$PCGrimAge,"reduced","increased")

high_acc<-age_meth[age_meth$deltaAge>median(age_meth$deltaAge)&age_meth$Sample_Group=="UT",c("Pool_ID")]
low_acc<-age_meth[age_meth$deltaAge<median(age_meth$deltaAge)&age_meth$Sample_Group=="UT",c("Pool_ID")]

high_baseline<-wilcox.test(
  age_meth[age_meth$Pool_ID%in%high_acc&age_meth$Sample_Group=="LPS24",c("PCGrimAge")],
  age_meth[age_meth$Pool_ID%in%high_acc&age_meth$Sample_Group=="UT",c("PCGrimAge")],
  paired = TRUE)

low_baseline<-wilcox.test(
  age_meth[age_meth$Pool_ID%in%low_acc&age_meth$Sample_Group=="LPS24",c("PCGrimAge")],
  age_meth[age_meth$Pool_ID%in%low_acc&age_meth$Sample_Group=="UT",c("PCGrimAge")],
  paired = TRUE)

age_meth$high_low_grimage <- "high"
age_meth$high_low_grimage[which(age_meth$Pool_ID %in% low_acc)] <- "low"

age_meth$high_low_grimage <- factor(age_meth$high_low_grimage, levels = c("high", "low"))

age_meth$Sample_Group <- factor(age_meth$Sample_Group, levels = c("UT", "LPS24"))
ann_text.1<-data.frame(high_low_grimage=c("high", "low"), PCGrimAge=c(80),Sample_Group=1.5,p.value=c(high_baseline[[3]], low_baseline[[3]]))

plot_grimage.facet<-ggplot(age_meth, aes(Sample_Group, PCGrimAge))+
  geom_violin(aes(fill=Sample_Group))+
  geom_boxplot(colour="grey45",alpha=0.25,position=position_dodge(1),width=0.64)+
  geom_text(aes(label=paste0("P=",signif(p.value,2))), data = ann_text.1) +
  scale_y_continuous("methylation age (yrs)")+
  scale_x_discrete("Treatment")+
  facet_wrap(~high_low_grimage) +
  theme_bw() +
  theme(strip.background = element_rect(color="black", fill="dark grey", size=1.5, linetype="solid"), strip.text.x = element_text(
    size = 15, color = "white", face = "bold"))
plot_grimage.facet





p.out <- ((plot_Horvath|plot_Horvath.facet)/(plot_Hannum|plot_Hannum.facet)/(plot_phenoage|plot_phenoage.facet)/(plot_grimage|plot_grimage.facet)) + plot_layout(widths = c(1,2))

ggplot2::ggsave("supp_age_ests.jpg",
                width = 10,
                height = 14,
                dpi = 300
)




