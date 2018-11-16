library(ggplot2)
library(ggsci)
library(dplyr)
library(cowplot)
library(plyr)
library(tidyr)
#setwd("/home/tavshalom/Desktop/temp")
setwd("/home/avshalom/projects/def-wyeth/avshalom/clinvar_anno")

figtheme <- theme_bw() + 
  theme(axis.text.x = element_text(size=18),
        axis.text.y = element_text(size=18),
        axis.title.x = element_text(size=20),
        axis.title.y = element_text(size=20),
        legend.text = element_text(size=18),
        legend.title = element_text(size=20)
  )
theme_set(figtheme)

clnsig_binner <- function(str){
  if (grepl("pathogenic", str, ignore.case = TRUE)) {
    return("pathogenic")
  }
  else if (grepl("benign", str, ignore.case = TRUE)) {
    return ("benign")
  } 
  else {
    return("other")
  }
}

full_df <- read.csv("parse.clinvar.anno.tsv", sep= "\t", header = TRUE)
#full_df <- read.csv("check.tsv", sep= "\t", header = TRUE)
full_df$X.CHROM <- as.factor(full_df$X.CHROM)
full_df$CLNSIG <- simplify2array(lapply(full_df$CLNSIG, FUN=clnsig_binner))

MC <- full_df %>% select(CLNSIG, MC)
MC <- na.omit(MC)
MC <- MC %>% mutate(D=gsub("\\,.*","",MC))
MC <- separate(data = MC, col = D, into = c("transcript", "MC_name"), sep = "\\|")
MC_freq <- count(MC, c('MC_name', 'CLNSIG'))
MC_freq <- MC_freq[MC_freq$CLNSIG != "other",]
MC_plot <- ggplot(MC_freq, aes(MC_name, freq, fill=CLNSIG))+
  geom_bar(stat="identity") +
  labs(y="Count", x="MC")  +
  theme(axis.text.x = element_text(angle = 40, hjust = 1))+
  scale_fill_manual(values=c("#0072B5FF", "#BC3C29FF"))
save_plot("MC_bar.png", MC_plot,base_height = 9,base_width = 16)

##INFO=<ID=ANN,Number=.,Type=String,Description="Functional annotations: 'Allele | Annotation | Annotation_Impact | Gene_Name | Gene_ID | Feature_Type | Feature_ID | Transcript_BioType | Rank | HGVS.c | HGVS.p | cDNA.pos / cDNA.length | CDS.pos / CDS.length | AA.pos / AA.length | Distance | ERRORS / WARNINGS / INFO' ">
ANN <- full_df %>% select(CLNSIG, ANN)
ANN <- ANN %>% mutate(ANN=gsub("\\,.*","",ANN))
ANN <- separate(data = ANN, col = ANN, into=c("Allele", "Annotation", "Annotation_Impact", "Gene_Name", 'Gene_ID', "Feature_Type", "Feature_ID", "Transcript_BioType", "Rank", "HGVS.c", "HGVS.p", "cDNA.pos", 'CDS.pos', "AA.pos" , "Distance", "ERRORS"), sep = "\\|")
gene <- count(ANN, c('Gene_Name', 'CLNSIG'))
gene <- gene[gene$CLNSIG == "pathogenic",]
gene <- gene[order(-gene$freq),]
gene_top20 <- gene[order(-gene$freq),][0:20,]

ANN_hist <- ggplot(gene, aes(x=freq))+
  geom_histogram(color="black", fill="#BC3C29FF", binwidth=0.05)+
  labs(y="Count", x="Variants/Gene")+ 
  geom_vline(aes(xintercept=median(freq)),
             color="#0072B5FF", linetype="dashed", size=1)+
  xlim(0,1000)+
  scale_x_continuous(trans='log10')
save_plot("ann_hist.png", ANN_hist,base_height = 9,base_width = 16)

ANN_plot <- ggplot(gene_top20, aes(Gene_Name, freq, fill=CLNSIG))+
  geom_bar(stat="identity") +
  labs(y="Count", x="Gene Name")  +
  theme(axis.text.x = element_text(angle = 40, hjust = 1))+
  scale_fill_manual(values=c("#BC3C29FF"))
save_plot("ann_bar.png", ANN_plot,base_height = 9,base_width = 16)
##CADD############
CADD <- full_df %>% select(X.CHROM,CLNSIG, CADD)
CADD <- na.omit(CADD)
CADD_boxplot_jitter <- ggplot(CADD, aes(CLNSIG, CADD)) + 
  geom_boxplot(aes(fill = CLNSIG), alpha=0.2) +
  geom_jitter(aes(color = CLNSIG), shape=16, position=position_jitter(0.4), size = 0.5, alpha=0.8)+ 
  labs(y="CADD Score", x="Pathogenicity") +
  theme(legend.position="none")+
  scale_fill_nejm()+
  scale_color_nejm()+
  geom_hline(yintercept=20, size=1)
save_plot("CADD_jitter.png", CADD_boxplot_jitter,base_height = 9,base_width = 16)

##gnomad_exome_af_global############
CADD <- full_df %>% select(X.CHROM,CLNSIG, gnomad_exome_af_global)
CADD <- na.omit(CADD)
length(unique(CADD$CLNSIG))
CADD_boxplot_jitter <- ggplot(CADD, aes(CLNSIG, gnomad_exome_af_global)) + 
  geom_boxplot(aes(fill = CLNSIG), alpha=0.2) +
  geom_jitter(aes(color = CLNSIG), shape=16, position=position_jitter(0.4), size = 0.5, alpha=0.8)+ 
  labs(y="gnomad_exome_af_global Score", x="Pathogenicity") +
  theme(legend.position="none")+
  scale_fill_nejm()+
  scale_color_nejm()+
  geom_hline(yintercept=0.01, size=1)+ 
  scale_y_continuous(trans='log10')
save_plot("gnomad_exome_af_global_jitter.png", CADD_boxplot_jitter,base_height = 9,base_width = 16)



##gnomad_exome_hom_global############
CADD <- full_df %>% select(X.CHROM,CLNSIG, gnomad_exome_hom_global)
CADD <- na.omit(CADD)
CADD_boxplot_jitter <- ggplot(CADD, aes(CLNSIG, gnomad_exome_hom_global)) + 
  geom_boxplot(aes(fill = CLNSIG), alpha=0.2) +
  geom_jitter(aes(color = CLNSIG), shape=16, position=position_jitter(0.4), size = 0.5, alpha=0.8)+ 
  labs(y="gnomad_exome_hom_global Score", x="Pathogenicity") +
  theme(legend.position="none")+
  scale_fill_nejm()+
  scale_color_nejm()+
  geom_hline(yintercept=10, size=1)+ 
  scale_y_continuous(trans='log10')
save_plot("gnomad_exome_hom_global_jitter.png", CADD_boxplot_jitter,base_height = 9,base_width = 16)



##gnomad_genome_af_global############
CADD <- full_df %>% select(X.CHROM,CLNSIG, gnomad_genome_af_global)
CADD <- na.omit(CADD)
CADD_boxplot_jitter <- ggplot(CADD, aes(CLNSIG, gnomad_genome_af_global)) + 
  geom_boxplot(aes(fill = CLNSIG), alpha=0.2) +
  geom_jitter(aes(color = CLNSIG), shape=16, position=position_jitter(0.4), size = 0.5, alpha=0.8)+ 
  labs(y="gnomad_genome_af_global Score", x="Pathogenicity") +
  theme(legend.position="none")+
  scale_fill_nejm()+
  scale_color_nejm()+
  geom_hline(yintercept=0.01, size=1)+ 
  scale_y_continuous(trans='log10')
save_plot("gnomad_genome_af_global_jitter.png", CADD_boxplot_jitter,base_height = 9,base_width = 16)



##gnomad_genome_hom_global############
CADD <- full_df %>% select(X.CHROM,CLNSIG, gnomad_genome_hom_global)
CADD <- na.omit(CADD)
CADD_boxplot_jitter <- ggplot(CADD, aes(CLNSIG, gnomad_genome_hom_global)) + 
  geom_boxplot(aes(fill = CLNSIG), alpha=0.2) +
  geom_jitter(aes(color = CLNSIG), shape=16, position=position_jitter(0.4), size = 0.5, alpha=0.8)+ 
  labs(y="gnomad_genome_hom_global Score", x="Pathogenicity") +
  theme(legend.position="none")+
  scale_fill_nejm()+
  scale_color_nejm()+
  geom_hline(yintercept=10, size=1)+ 
  scale_y_continuous(trans='log10')
save_plot("gnomad_genome_hom_global_jitter.png", CADD_boxplot_jitter,base_height = 9,base_width = 16)
