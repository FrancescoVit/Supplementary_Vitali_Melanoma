#################################################################################################################
### Script to perform data analysis and obtain figure of Vitali et. al 2020 manuscript                  #########
### Author: Francesco Vitali                                                                            #########
### Licence: This work is licensed under a Creative Commons Attribution 4.0 International License       #########
#################################################################################################################

#Load Package
library(tidyverse)
library(colorblindr)
library(reshape2)
library(ggsci)
library(vegan)
library(RColorBrewer)
library(NbClust)
library(FactoMineR)
library(metagenomeSeq)
library(DESeq2)
library(ggpubr)
library(microbiome) 
library(phyloseq)
library(scales)
library(doBy)
library(factoextra)
library(plotly)
library(xtable)
library(knitr)
library(calibrate)
library(data.table)
library(selbal)
library(mvabund)
library(cowplot)
library(awtools)
library(dutchmasters)
library(plyr)
library(picante)
library("randomForest")
library("plyr") # for the "arrange" function
library("rfUtilities") # to test model significance
library("caret") # to get leave-one-out cross-validation accuracies and also contains the nearZeroVar function 


#################
## FILE IMPORT ##
#################

otus_16S_fp <- "./otutable16S_melanoma.csv"
meta_16S_fp <- "./METADATA_MELANOMA_16S.csv"
tax_16S_fp <- "./taxtable16S_melanoma.csv"
#
tax_16S_table <- read.csv(tax_16S_fp, header = F, row.names = 1, sep = "\t")
tax_16S_tt <- tax_table(as.matrix(tax_16S_table))
#
otu_16S_table <- read.csv(otus_16S_fp, header = T, row.names = 1, sep = "\t")
otu_16S_ot <- otu_table(as.matrix(otu_16S_table), taxa_are_rows = T)
#
metadata_16S_table <- read.csv(meta_16S_fp, header = T, row.names = 1, sep = "\t")
metadata_16S_sd <- sample_data(metadata_16S_table)

# mount the phyloseq object
qiime_file_16S_initial <- merge_phyloseq(otu_16S_ot, tax_16S_tt, metadata_16S_sd)
sample_names(qiime_file_16S_initial)


# rename rank
colnames(tax_table(qiime_file_16S_initial)) <- c(Rank1 = "Kingdom", Rank2 = "Phylum", Rank3 = "Class", Rank4 = "Order", Rank5 = "Family", Rank6 = "Genus")

# check for mitochondria and chloroplast 
taxa <- as.data.frame(tax_table(qiime_file_16S_initial)@.Data)

na.omit(taxa[taxa$Family == "Mitochondria",])
na.omit(taxa[taxa$Family == "Chloroplast",])

#qiime_file_16S_initial <- subset_taxa(qiime_file_16S_initial, Family  != "Mitochondria")
#qiime_file_16S_initial <- subset_taxa(qiime_file_16S_initial, Family  != "Chloroplast")
#qiime_file_16S_initial

# Remove C17
qiime_file_16S_initial <- subset_samples(qiime_file_16S_initial, Label != "C17F")
summary(sample_sums(qiime_file_16S_initial)) # Variation interval 31736 - 172771
sum(sample_sums(qiime_file_16S_initial))

# Reomve ML sample (after initial MS revision)
qiime_file_16S_initial <- subset_samples(qiime_file_16S_initial, Pathology != "Leucoderma_Melanoma") %>% # remove otus = 0
  prune_taxa(taxa_sums(.) > 0, .)
summary(sample_sums(qiime_file_16S_initial)) # Variation interval 31736 - 172771
sum(sample_sums(qiime_file_16S_initial))


## otus different from 0
# converting otu table to 1/0
otu_16S_10 <- ifelse(otu_16S_ot >0,1,0)
summary(sample_sums(otu_16S_10))


## obtaining rarefaction 
qiime_file_16S_raref <- rarefy_even_depth(qiime_file_16S_initial, 
                                          rngseed=1234, 
                                          sample.size=0.99*min(sample_sums(qiime_file_16S_initial)), 
                                          replace=F)


# CSS scaling
# singletons remove
summary(taxa_sums(qiime_file_16S_initial))
boxplot(taxa_sums(qiime_file_16S_initial))
doubleton <- genefilter_sample(qiime_file_16S_initial, filterfun_sample(function(x) x > 1), A=1)
doubleton <- prune_taxa(doubleton, qiime_file_16S_initial) 
summary(taxa_sums(doubleton))
# transforming see also mixomics, I'll use CSS
data.metagenomeSeq = phyloseq_to_metagenomeSeq(doubleton)
p = cumNormStat(data.metagenomeSeq)
data.cumnorm = cumNorm(data.metagenomeSeq, p=p)
#data.cumnorm
data.CSS = MRcounts(data.cumnorm, norm=TRUE, log=TRUE)
dim(data.CSS)  # make sure the data are in a correct formal: number of samples in rows
qiime_file_16S_css <- qiime_file_16S_initial
otu_table(qiime_file_16S_css) <- otu_table(data.CSS, taxa_are_rows = T)

######### FUNGI ###########

otus_ITS_fp <- "./otutableITS_melanoma.csv"
meta_ITS_fp <- "./METADATA_MELANOMA_ITS.csv"
tax_ITS_fp <- "./taxtableITS_melanoma.csv"
#
tax_ITS_table <- read.csv(tax_ITS_fp, header = F, row.names = 1, sep = "\t")
tax_ITS_tt <- tax_table(as.matrix(tax_ITS_table))
#
otu_ITS_table <- read.csv(otus_ITS_fp, header = T, row.names = 1, sep = "\t")
otu_ITS_ot <- otu_table(as.matrix(otu_ITS_table), taxa_are_rows = T)
#
metadata_ITS_table <- read.csv(meta_ITS_fp, header = T, row.names = 1, sep = "\t")
metadata_ITS_sd <- sample_data(metadata_ITS_table)

# mount the phyloseq object
qiime_file_ITS_initial <- merge_phyloseq(otu_ITS_ot, tax_ITS_tt, metadata_ITS_sd)
colnames(tax_table(qiime_file_ITS_initial)) <- c(V2 = "Kingdom", V3= "Phylum", V4 = "Class", V5 = "Order", V6 = "Family", V7 = "Genus", V8 = "Species")

sample_names(qiime_file_16S_initial)

# Remove C17
qiime_file_ITS_initial <- subset_samples(qiime_file_ITS_initial, Label != "C17F")

summary(sample_sums(qiime_file_ITS_initial)) 
sum(sample_sums(qiime_file_ITS_initial))

# Reomve ML sample (after initial MS revision)
qiime_file_ITS_initial <- subset_samples(qiime_file_ITS_initial, Pathology != "Leucoderma_Melanoma")  %>% # remove otus = 0
  prune_taxa(taxa_sums(.) > 0, .)
summary(sample_sums(qiime_file_ITS_initial)) # Variation interval 31736 - 172771
sum(sample_sums(qiime_file_ITS_initial))

## otus different from 0
# converting otu table to 1/0
otu_ITS_10 <- ifelse(otu_ITS_ot >0,1,0)
summary(sample_sums(otu_ITS_10))


## obtaining rarefaction 
qiime_file_ITS_raref <- rarefy_even_depth(qiime_file_ITS_initial, 
                                          rngseed=1234, 
                                          sample.size=0.99*min(sample_sums(qiime_file_ITS_initial)), 
                                          replace=F)


# CSS scaling
doubleton <- genefilter_sample(qiime_file_ITS_initial, filterfun_sample(function(x) x > 1), A=1)
doubleton <- prune_taxa(doubleton, qiime_file_ITS_initial) 
# transforming 
data.metagenomeSeq = phyloseq_to_metagenomeSeq(doubleton)
p = cumNormStat(data.metagenomeSeq) #default is 0.5
data.cumnorm = cumNorm(data.metagenomeSeq, p=p)
#data.cumnorm
data.CSS = MRcounts(data.cumnorm, norm=TRUE, log=TRUE)
dim(data.CSS)  # make sure the data are in a correct formal: number of samples in rows
qiime_file_ITS_css <- qiime_file_ITS_initial
otu_table(qiime_file_ITS_css) <- otu_table(data.CSS, taxa_are_rows = T)


########### Joining ITS and 16S scaled tables

sample_data(qiime_file_16S_css)

#creo oggetto phyloseq con nomi campione = lable, in modo da avere nomi uguali tra 16S e ITS; phyloseq obj with names = label so they are same between 16S adn ITS
qiime_file_16S_css_join <- qiime_file_16S_css
sample_names(qiime_file_16S_css_join) <- sample_data(qiime_file_16S_css_join)$Label
qiime_file_ITS_css_join <- qiime_file_ITS_css
sample_names(qiime_file_ITS_css_join) <- sample_data(qiime_file_ITS_css_join)$Label
#check
sample_names(qiime_file_16S_css_join) == sample_names(qiime_file_ITS_css_join)

# extract tax and otu table
otu_16S_join <- otu_table(qiime_file_16S_css_join)
row.names(otu_16S_join) <- paste(row.names(otu_16S_join),"16S", sep = "|")
otu_ITS_join <- otu_table(qiime_file_ITS_css_join)
row.names(otu_ITS_join) <- paste(row.names(otu_ITS_join),"ITS", sep = "|")

tax_16S_join <- tax_table(qiime_file_16S_css_join)
row.names(tax_16S_join) <- paste(row.names(tax_16S_join),"16S", sep = "|")

#add variable species = genera to use tax_glom() at the species lavel
tax_16S_join <- as.data.frame(tax_16S_join)
tax_16S_join$Species <- tax_16S_join$Genus
tax_ITS_join <- tax_table(qiime_file_ITS_css_join)
row.names(tax_ITS_join) <- paste(row.names(tax_ITS_join),"ITS", sep = "|")
tax_ITS_join <- as.data.frame(tax_ITS_join)

# cbind
colnames(tax_16S_join)
colnames(tax_ITS_join)

taxa_16SITS <- rbind(tax_16S_join, tax_ITS_join)
otu_16SITS <- rbind(otu_16S_join, otu_ITS_join)
str(taxa_16SITS)
str(otu_16SITS)

# create joined phyloseq obj
otus <- otu_table(otu_16SITS,taxa_are_rows = T)
taxas <- tax_table(taxa_16SITS)
row.names(taxas) <- row.names(otus)
phyloseq_16SITS <- phyloseq(otu_table(otus), tax_table(taxas), sample_data(sample_data(qiime_file_16S_css_join)))

#####################################################
#####################################################

## Look at the reads!

### 16S
read.per.sample <- data.frame(sample = sample_names(qiime_file_16S_initial), reads = sample_sums(qiime_file_16S_initial))
read.per.sample.order <- read.per.sample[with(read.per.sample, order(reads)), ]
d <- density(read.per.sample.order$reads)
plot(d, type="n", main="Seq per sample Bacteria", xlim = c(0,300000))
polygon(d, col="lightgray", border="gray")
rug(read.per.sample$reads, col="red")
sdt = data.table(as(sample_data(qiime_file_16S_initial), "data.frame"),
                 TotalReads = sample_sums(qiime_file_16S_initial), keep.rownames = TRUE)
setnames(sdt, "rn", "SampleID")
#n tot reads
sums <- sample_sums(qiime_file_16S_initial)
min(sums)
max(sums)
sum(sums)
p_rich <- ggboxplot(data = sdt, x = "Pathology", y = "TotalReads", add = "jitter", fill = "Pathology", 
                    palette ="npg")
p_rich <- p_rich + stat_compare_means(label.y = 150000, method = "anova")
p_rich

### ITS
read.per.sample <- data.frame(sample = sample_names(qiime_file_ITS_initial), reads = sample_sums(qiime_file_ITS_initial))
read.per.sample.order <- read.per.sample[with(read.per.sample, order(reads)), ]
d <- density(read.per.sample.order$reads)
plot(d, type="n", main="Seq per sample Fungi",xlim = c(0,300000))
polygon(d, col="lightgray", border="gray")
rug(read.per.sample$reads, col="red")
sdt = data.table(as(sample_data(qiime_file_ITS_initial), "data.frame"),
                 TotalReads = sample_sums(qiime_file_ITS_initial), keep.rownames = TRUE)
setnames(sdt, "rn", "SampleID")
#n tot reads
sums <- sample_sums(qiime_file_ITS_initial)
min(sums)
max(sums)
sum(sums)
p_rich <- ggboxplot(data = sdt, x = "Pathology", y = "TotalReads", add = "jitter", fill = "Pathology", 
                    palette ="npg")
p_rich <- p_rich + stat_compare_means(label.y = 200000, method = "anova")
p_rich

############################################
#### DIVERSITY AND ORDINATION ####
############################################


mel_con <- subset_samples(phyloseq_16SITS, Pathology == "Melanoma" | Pathology == "Control")
leve <- c("Control", "Melanoma")
facto <- factor(sample_data(mel_con)$Pathology, levels = leve)
sample_data(mel_con)$Pathology <- facto
mel_con_noM6 <- subset_samples(mel_con, Label != "M6F")

mel_con_bact <- subset_samples(qiime_file_16S_css, Pathology == "Melanoma" | Pathology == "Control")
leve <- c("Control", "Melanoma")
facto <- factor(sample_data(mel_con_bact)$Pathology, levels = leve)
sample_data(mel_con_bact)$Pathology <- facto
mel_con_noM6_bact <- subset_samples(mel_con_bact, Label != "M6F")

mel_con_fung <- subset_samples(qiime_file_ITS_css, Pathology == "Melanoma" | Pathology == "Control")
leve <- c("Control", "Melanoma")
facto <- factor(sample_data(mel_con_fung)$Pathology, levels = leve)
sample_data(mel_con_fung)$Pathology <- facto
mel_con_noM6_fung <- subset_samples(mel_con_fung, Label != "M6F")

mel_nom6 <- subset_samples(mel_con_noM6, Pathology != "Control")


##################################
######## BUILDING FIGURE 1 #######
##################################

# BACT ORDINATION

ordination_bact <- ordinate(subset_samples(mel_con_noM6_bact, Pathology != "Leucoderma_Melanoma"), "PCoA", distance = "bray")
p_ordination_bact <- plot_ordination(physeq = subset_samples(mel_con_noM6_bact, Pathology != "Leucoderma_Melanoma"), ordination_bact, axes = c(1,2), color="Pathology")
p12 <- p_ordination_bact + geom_point(size = 7, shape = 20) 
p12 <- p12 + theme_bw() + 
  theme( panel.grid.major =  element_blank(),
         panel.grid.minor =  element_blank())
p12 <- p12 + scale_color_aaas() + xlab("PCoA 1 (14.0%)") +
  ylab("PCoA 2 (7.7%)")
p12 + stat_ellipse(type = "t", linetype = 2) + theme(legend.position="none")

set.seed(1233)
datas <- as(sample_data(subset_samples(mel_con_noM6_bact, Pathology != "Leucoderma_Melanoma")), "data.frame")
adonis2(distance(subset_samples(mel_con_noM6_bact, Pathology != "Leucoderma_Melanoma"), method = "bray") ~ Pathology, data = datas, permutations = 9999)


## BACT BARPLOT FAMILY

library(pals)
phyloseq_16S_noM6 <- subset_samples(qiime_file_16S_css, Label != "M6F" & Pathology != "Vitiligo")
paletta <- colorRampPalette(as.character(alphabet(26)))(30)
otu_table_collapsed <- merge_samples(qiime_file_16S_css, group = "Pathology")
otu_table_collapsed <- tax_glom(otu_table_collapsed, taxrank="Family")
qiime_file_proportional <- transform_sample_counts(otu_table_collapsed, function(x) 100 * x/sum(x))
qiime_file_proportional_oneperc <-  filter_taxa(qiime_file_proportional, function(x) max(x) > 1, TRUE)
p_phylum <- plot_bar(qiime_file_proportional_oneperc, fill= "Family")
p_phylum <- p_phylum + theme_linedraw() + scale_fill_manual(values= paletta) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position="right") + ylab("Abundance (%)") + xlab("")
p_phylum + guides(fill = guide_legend(nrow = 15)) + scale_x_discrete(limits=c("Control","Melanoma"))


# FUNGI ORDINATION

ordination_fung <- ordinate(subset_samples(mel_con_noM6_fung, Pathology != "Leucoderma_Melanoma"), "PCoA", distance = "bray")
p_ordination_fung <- plot_ordination(physeq = subset_samples(mel_con_noM6_fung, Pathology != "Leucoderma_Melanoma"), ordination_fung, axes = c(1,2), color="Pathology")
p12 <- p_ordination_fung + geom_point(size = 7, shape = 20) 
p12 <- p12 + theme_bw() + 
  theme( panel.grid.major =  element_blank(),
         panel.grid.minor =  element_blank())
p12 <- p12 + scale_color_aaas() + xlab("PCoA 1 (12.9%)") +
  ylab("PCoA 2 (11%)")
p12 + stat_ellipse(type = "t", linetype = 2) + theme(legend.position="none")

set.seed(1233)
datas <- as(sample_data(subset_samples(mel_con_noM6_fung, Pathology != "Leucoderma_Melanoma")), "data.frame")
adonis2(distance(subset_samples(mel_con_noM6_fung, Pathology != "Leucoderma_Melanoma"), method = "bray") ~ Pathology, data = datas, permutations = 9999)



## FUNGI BARPLOT 

library(pals)
phyloseq_ITS_noM6 <- subset_samples(qiime_file_ITS_css, Label != "M6F" & Pathology != "Vitiligo")
sost <- revalue(sample_data(phyloseq_ITS_noM6)$Pathology, c("Leucoderma_Melanoma" = "Leucoderma_Melanoma"))
sample_data(phyloseq_ITS_noM6)$Pathology <- sost
paletta <- colorRampPalette(as.character(alphabet2(26)))(40)
otu_table_collapsed <- merge_samples(phyloseq_ITS_noM6, group = "Pathology")
otu_table_collapsed <- tax_glom(otu_table_collapsed, taxrank="Order")
qiime_file_proportional <- transform_sample_counts(otu_table_collapsed, function(x) 100 * x/sum(x))
qiime_file_proportional_oneperc <-  filter_taxa(qiime_file_proportional, function(x) max(x) > 1, TRUE)
p_phylum <- plot_bar(qiime_file_proportional_oneperc, fill= "Order")
p_phylum <- p_phylum + theme_linedraw() + scale_fill_manual(values= paletta) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position="right") + ylab("Abundance (%)") + xlab("")
p_phylum + guides(fill = guide_legend(nrow = 15)) + scale_x_discrete(limits=c("Control","Melanoma"))



### BACTERIA AND FUNGI

ordination_fungBact <- ordinate(subset_samples(mel_con_noM6, Pathology != "Leucoderma_Melanoma"), "PCoA", distance = "bray")
p_ordination_fungBact <- plot_ordination(physeq = subset_samples(mel_con_noM6, Pathology != "Leucoderma_Melanoma"), ordination_fungBact, axes = c(1,2), color="Pathology")
p12 <- p_ordination_fungBact + geom_point(size = 7, shape = 20) 
p12 <- p12 + theme_bw() + 
  theme( panel.grid.major =  element_blank(),
         panel.grid.minor =  element_blank())
p12 <- p12 + scale_color_aaas() + xlab("PCoA 1 (13%)") +
  ylab("PCoA 2 (7.3%)")
p12 + stat_ellipse(type = "t", linetype = 2) + theme(legend.position="none")

set.seed(1233)
datas <- as(sample_data(subset_samples(mel_con_noM6, Pathology != "Leucoderma_Melanoma")), "data.frame")
adonis2(distance(subset_samples(mel_con_noM6, Pathology != "Leucoderma_Melanoma"), method = "bray") ~ Pathology, data = datas, permutations = 9999)


## ALL SAMPLES ORDINATION

mel_ML_con <- subset_samples(phyloseq_16SITS, Pathology != "Vitiligo")
sample_data(mel_ML_con)$Pathology
leve <- c("Control", "Melanoma", "Leucoderma_Melanoma")
facto <- factor(sample_data(mel_ML_con)$Pathology, levels = leve)
sample_data(mel_ML_con)$Pathology <- facto
mel_ML_con_noM6 <- subset_samples(mel_ML_con, Label != "M6F")
sample_data(mel_ML_con_noM6)$Pathology

mel_ML_con_bact <- subset_samples(qiime_file_16S_css, Pathology != "Vitiligo")
leve <- c("Control", "Melanoma", "MAL-melanoma")
facto <- factor(sample_data(mel_ML_con_bact)$Pathology, levels = leve)
sample_data(mel_ML_con_bact)$Pathology <- facto
mel_ML_con_noM6_bact <- subset_samples(mel_ML_con_bact, Label != "M6F")

mel_ML_con_fung <- subset_samples(qiime_file_ITS_css, Pathology != "Vitiligo")
leve <- c("Control", "Melanoma", "Leucoderma_Melanoma")
facto <- factor(sample_data(mel_ML_con_fung)$Pathology, levels = leve)
sample_data(mel_ML_con_fung)$Pathology <- facto
mel_ML_con_noM6_fung <- subset_samples(mel_ML_con_fung, Label != "M6F")


ordination_fungBact <- ordinate(mel_ML_con, "PCoA", distance = "bray")
p_ordination_fungBact <- plot_ordination(physeq = mel_ML_con, ordination_fungBact, axes = c(1,2), color="Pathology")
p12 <- p_ordination_fungBact + geom_point(size = 7, shape = 20) 
p12 <- p12 + theme_bw() + 
  theme( panel.grid.major =  element_blank(),
         panel.grid.minor =  element_blank())
p12 <- p12 + scale_color_aaas() + xlab("PCoA 1 (12.4%)") +
  ylab("PCoA 2 (7.6%)")
p12 + theme(legend.position="none")

#### ALPHA DIVERSITY #
## USING microbiome package


alphadiv <- microbiome::alpha(subset_samples(qiime_file_16S_raref, Pathology != "Vitiligo"))

df <- as(sample_data(subset_samples(qiime_file_16S_raref, Pathology != "Vitiligo")), "data.frame")

df$rich <- alphadiv$observed
df$shan <- alphadiv$diversity_shannon
df$eve <- alphadiv$evenness_pielou
df$sim <- alphadiv$diversity_inverse_simpson

df$Pathology <- factor(df$Pathology, levels = c("Control", "Melanoma"))

# Plot on Pathology
variabile_plot = "Pathology"
my_comparisons <- list(c("Control","Melanoma"))

p_rich_b <- ggboxplot(data = df, x = variabile_plot, y = "rich", add = "jitter", size = 0.5, fill = variabile_plot, outlier.shape = NA)
p_rich_b <- p_rich_b + stat_compare_means(method = "t.test", comparisons = my_comparisons)                   # Add global p-value
p_rich_b <- p_rich_b + ylab("Richness (n° of OTUs)") + scale_fill_aaas() + rremove(object = "legend")
p_shan <- ggdotplot(data = df, x = variabile_plot, y = "shan", add = "mean_sd", size = 0.5, fill = variabile_plot)
p_shan <- p_shan + stat_compare_means(method = "wilcox.test", comparisons = my_comparisons)                   # Add global p-value
p_shan <- p_shan + ylab("Shannon's Index") + scale_fill_aaas()
p_eve_b <- ggboxplot(data = df, x = variabile_plot, y = "eve", add = "jitter", size = 0.5, fill = variabile_plot, outlier.shape = NA)
p_eve_b <- p_eve_b + stat_compare_means(method = "t.test", comparisons = my_comparisons)                   # Add global p-value
p_eve_b <- p_eve_b + ylab("Evenness Index") + scale_fill_aaas() + rremove(object = "legend")
p_sim <- ggdotplot(data = df, x = variabile_plot, y = "sim", add = "mean_sd", size = 0.5, fill = variabile_plot)
p_sim <- p_sim + stat_compare_means(method = "wilcox.test", comparisons = my_comparisons)                   # Add global p-value
p_sim <- p_sim + ylab("Simpson Index") + scale_fill_aaas()
# 
 plot_grid(p_rich, p_eve,p_shan,p_sim, ncol = 4)

plot_grid(p_rich_b, p_eve_b, ncol = 2)

#### Fungi


alphadiv <- microbiome::alpha(qiime_file_ITS_raref, index = c("observed", "evenness_pielou", "diversity_shannon", "diversity_inverse_simpson"))

df <- as(sample_data(subset_samples(qiime_file_ITS_raref, Pathology != "Vitiligo")), "data.frame")

df$rich <- alphadiv$observed
df$shan <- alphadiv$diversity_shannon
df$eve <- alphadiv$evenness_pielou
df$sim <- alphadiv$diversity_inverse_simpson

df$Pathology <- factor(df$Pathology, levels = c("Control", "Melanoma", "Leucoderma_Melanoma"))


# Plot on Pathology
variabile_plot = "Pathology"
my_comparisons <- list(c("Control","Melanoma"))

p_rich_f <- ggboxplot(data = df, x = variabile_plot, y = "rich", add = "jitter", size = 0.5, fill = variabile_plot, outlier.shape = NA)
p_rich_f <- p_rich_f + stat_compare_means(method = "t.test", comparisons = my_comparisons)                   # Add global p-value
p_rich_f <- p_rich_f + ylab("Richness (n° of OTUs)") + scale_fill_aaas() + rremove("legend")
p_shan <- ggdotplot(data = df, x = variabile_plot, y = "shan", add = "mean_sd", size = 0.5, fill = variabile_plot)
p_shan <- p_shan + stat_compare_means(method = "wilcox.test", comparisons = my_comparisons)                   # Add global p-value
p_shan <- p_shan + ylab("Shannon's Index") + scale_fill_aaas()
p_eve_f <- ggboxplot(data = df, x = variabile_plot, y = "eve", add = "jitter", size = 0.5, fill = variabile_plot, outlier.shape = NA)
p_eve_f <- p_eve_f + stat_compare_means(method = "t.test", comparisons = my_comparisons)                   # Add global p-value
p_eve_f <- p_eve_f + ylab("Evenness Index") + scale_fill_aaas() + rremove("legend")
p_sim <- ggdotplot(data = df, x = variabile_plot, y = "sim", add = "mean_sd", size = 0.5, fill = variabile_plot)
p_sim <- p_sim + stat_compare_means(method = "wilcox.test", comparisons = my_comparisons)                   # Add global p-value
p_sim <- p_sim + ylab("Simpson Index")+ scale_fill_aaas()

plot_grid(p_rich_f, p_shan, p_eve_f, p_sim, ncol = 4)

plot_grid(p_rich_b, p_eve_b, p_rich_f, p_eve_f, ncol = 2)

###################################################
###### BUILDING SUPPLEMENTARY FIGURE 1 and 2 ######
###################################################

#BACTERIA
#basic dendrogram
xyz <- qiime_file_16S_css
sample_data(xyz)
sample_names(xyz) <- sample_data(xyz)$Label
d.bact = phyloseq::distance(xyz, "bray")
upgma <- hclust(d.bact, "complete")
plot(upgma)
library(dendextend)
library(RColorBrewer)
dend <- as.dendrogram(upgma)
# We hang the dendrogram a bit:
dend <- hang.dendrogram(dend)
dend <- set(dend, "labels_cex", 0.5)
# Create a vector of colors
data <- as.data.frame(sample_data(qiime_file_16S_css))
my_colors2 <- as.factor(data$Pathology)
pal_futu <- pal_aaas("default")(3)
my_colors2 <- revalue(my_colors2, c("Control" = pal_futu[1], "Melanoma" = pal_futu[2], "Leucoderma_Melanoma" = pal_futu[3]))
my_colors2 <- as.character(my_colors2)
# Make the dendrogram
par(mar = c(4,4,1,12))
plot(dend, horiz = T)
colored_bars(colors = my_colors2, dend = dend,rowLabels = "Pathology", horiz = T)

# barplot at family with same order as dendrorgam

x <- as.factor(labels(dend)) # bottom to top order


library(pals)
paletta <- colorRampPalette(as.character(glasbey(30)))(44)
otu_table_collapsed <- tax_glom(qiime_file_16S_css, taxrank="Family")
df <- as.data.frame(sample_data(qiime_file_16S_css))
df_ord <- df[match(x, df$Label),]
df_ord$Label == x
df_ord$Label <- factor(x = df_ord$Label, levels = df_ord$Label, ordered = T)
sample_data(mel_con_bact) <- df_ord
p_phylum <- plot_bar(otu_table_collapsed, fill= "Family", x = "Label")
p_phylum <- p_phylum + theme_linedraw() + scale_fill_manual(values= paletta) + coord_flip() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + ylab("Abundance CSS/log scaled") + xlab("")
p_phylum <- p_phylum + guides(fill = guide_legend(nrow = 22))



##############
#### FUNGI ###

#basic dendrogram
xyz <- qiime_file_ITS_css
sample_data(xyz)
sample_names(xyz) <- sample_data(xyz)$Label
d.fung = phyloseq::distance(xyz, "bray")
upgma <- hclust(d.fung, "complete")
plot(upgma)
library(dendextend)
library(RColorBrewer)
dend <- as.dendrogram(upgma)
# We hang the dendrogram a bit:
dend <- hang.dendrogram(dend)
dend <- set(dend, "labels_cex", 0.5)
# Create a vector of colors
data <- as.data.frame(sample_data(qiime_file_ITS_css))
my_colors2 <- as.factor(data$Pathology)
pal_futu <- pal_aaas("default")(3)
my_colors2 <- revalue(my_colors2, c("Control" = pal_futu[1], "Melanoma" = pal_futu[2], "Leucoderma_Melanoma" = pal_futu[3]))
my_colors2 <- as.character(my_colors2)
# Make the dendrogram
par(mar = c(4,4,1,12))
plot(dend, horiz = T)
colored_bars(colors = my_colors2, dend = dend,rowLabels = "Pathology", horiz = T)

# bbarplot at family with same order as dendrorgam

x <- as.factor(labels(dend)) # bottom to top order



library(pals)
paletta <- colorRampPalette(as.character(glasbey(32)))(38)
otu_table_collapsed <- tax_glom(mel_con_fung, taxrank="Order")
df <- as.data.frame(sample_data(mel_con_fung))
df_ord <- df[match(x, df$Label),]
df_ord$Label == x
df_ord$Label <- factor(x = df_ord$Label, levels = df_ord$Label, ordered = T)
sample_data(mel_con_fung) <- df_ord
p_phylum <- plot_bar(otu_table_collapsed, fill= "Order", x = "Label")
p_phylum <- p_phylum + theme_linedraw() + scale_fill_manual(values= paletta) + coord_flip() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + ylab("Abundance CSS/log scaled") + xlab("")
p_phylum <- p_phylum + guides(fill = guide_legend(nrow = 19))



### BARPLOT for sample category

library(pals)
paletta <- as.character(glasbey(30))
otu_table_collapsed <- merge_samples(qiime_file_16S_css, group = "Pathology")
otu_table_collapsed <- tax_glom(otu_table_collapsed, taxrank="Family")
qiime_file_proportional <- transform_sample_counts(otu_table_collapsed, function(x) 100 * x/sum(x))
qiime_file_proportional_oneperc <-  filter_taxa(qiime_file_proportional, function(x) max(x) > 1, TRUE)
p_phylum <- plot_bar(qiime_file_proportional_oneperc, fill= "Family")
p_phylum <- p_phylum + theme_linedraw() + scale_fill_manual(values= paletta) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + ylab("Abundance (%)") + xlab("")
p_phylum + guides(fill = guide_legend(nrow = 15))



library(pals)

taxxxxa <- as.data.frame(tax_table(qiime_file_ITS_css))
taxxxxas <- taxxxxa$Species
taxxxxas <- sub("\\|.*", "", taxxxxas)
taxxxxa$Species <- taxxxxas

tax_table(qiime_file_ITS_css) <- tax_table(as.matrix(taxxxxa))

paletta <- colorRampPalette(as.character(glasbey(30)))(33)
otu_table_collapsed <- merge_samples(qiime_file_ITS_css, group = "Pathology")
otu_table_collapsed <- tax_glom(otu_table_collapsed, taxrank="Order")
qiime_file_proportional <- transform_sample_counts(otu_table_collapsed, function(x) 100 * x/sum(x))
qiime_file_proportional_oneperc <-  filter_taxa(qiime_file_proportional, function(x) max(x) > 1, TRUE)
p_phylum <- plot_bar(qiime_file_proportional_oneperc, fill= "Order")
p_phylum <- p_phylum + theme_linedraw() + scale_fill_manual(values= paletta) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + ylab("Abundance (%)") + xlab("")
p_phylum + guides(fill = guide_legend(nrow = 17))



##################################
######## BUILDING FIGURE 3    ####  
#######  BUILDING SUPP. FIG 4 ####
##################################

# same as above but on the untrasformed data for deseq analysis
# Creating joined bacteria and fungi dataset

qiime_file_16S_NOcss_join <- qiime_file_16S_initial
sample_names(qiime_file_16S_NOcss_join) <- sample_data(qiime_file_16S_NOcss_join)$Label
qiime_file_ITS_NOcss_join <- qiime_file_ITS_initial
sample_names(qiime_file_ITS_NOcss_join) <- sample_data(qiime_file_ITS_NOcss_join)$Label

sample_names(qiime_file_ITS_NOcss_join) == sample_names(qiime_file_ITS_NOcss_join)

otu_16S_join_NOcss <- otu_table(qiime_file_16S_NOcss_join)
row.names(otu_16S_join_NOcss) <- paste(row.names(otu_16S_join_NOcss),"16S", sep = "|")
otu_ITS_join_NOcss <- otu_table(qiime_file_ITS_NOcss_join)
row.names(otu_ITS_join_NOcss) <- paste(row.names(otu_ITS_join_NOcss),"ITS", sep = "|")


tax_16S_join_NOcss <- tax_table(qiime_file_16S_NOcss_join)
row.names(tax_16S_join_NOcss) <- paste(row.names(tax_16S_join_NOcss),"16S", sep = "|")


tax_16S_join_NOcss <- as.data.frame(tax_16S_join_NOcss)
tax_16S_join_NOcss$Species <- tax_16S_join_NOcss$Genus
tax_ITS_join_NOcss <- tax_table(qiime_file_ITS_NOcss_join)
row.names(tax_ITS_join_NOcss) <- paste(row.names(tax_ITS_join_NOcss),"ITS", sep = "|")
tax_ITS_join_NOcss <- as.data.frame(tax_ITS_join_NOcss)


taxa_16SITS_NOcss <- rbind(tax_16S_join_NOcss, tax_ITS_join_NOcss)
otu_16SITS_NOcss <- rbind(otu_16S_join_NOcss, otu_ITS_join_NOcss)
str(taxa_16SITS_NOcss)
str(otu_16SITS_NOcss)


otus <- otu_table(otu_16SITS_NOcss,taxa_are_rows = T)
taxas <- tax_table(taxa_16SITS_NOcss)
row.names(taxas) <- row.names(otus)
phyloseq_16SITS_NOcss <- phyloseq(otu_table(otus), tax_table(taxas), sample_data(sample_data(qiime_file_16S_NOcss_join)))


mel_con_NOcss <- subset_samples(phyloseq_16SITS_NOcss, Pathology == "Melanoma" | Pathology == "Control")
leve <- c("Control", "Melanoma")
facto <- factor(sample_data(mel_con_NOcss)$Pathology, levels = leve)
sample_data(mel_con_NOcss)$Pathology <- facto
mel_con_noM6_NOcss <- subset_samples(mel_con_NOcss, Label != "M6F")


  ## TEST ON PATHOLOGY

# following for reference https://bioconductor.org/packages/devel/bioc/vignettes/phyloseq/inst/doc/phyloseq-mixture-models.html
kosticB <- mel_con_noM6_NOcss
tax_table(kosticB)
sample_names(kosticB)

diagdds <- phyloseq_to_deseq2(mel_con_noM6_NOcss, ~ Pathology)
# calculate geometric means prior to estimate size factors
gm_mean <- function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
geoMeans <- apply(counts(diagdds), 1, gm_mean)
diagdds <- estimateSizeFactors(diagdds, geoMeans = geoMeans)
diagdds <- DESeq(diagdds, fitType="local",test = "Wald")

res <- results(diagdds,cooksCutoff = F,pAdjustMethod = "BH")
res <- res[order(res$padj, na.last=NA), ]
alpha <- 0.05
sigtab <- res[which(res$padj < alpha), ]
sigtab <- as(sigtab, "data.frame")
sigtax <- as.data.frame(tax_table(kosticB)@.Data)

rownames(sigtab) %in% sigtax[rownames(sigtab),]
sigtab = res[which(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(kosticB)[rownames(sigtab), ], "matrix"))
sigtab

list_CvsM <- row.names(sigtab)

# SUPPLEMENTARY TABLE 5
#write.table(x = sigtab, "./tabella_DESEQ_melanomaVScontrol.csv", sep = ",") # 


# SUPPLEMENTARY FIGURE 4
library(EnhancedVolcano)
passing <- cbind(as(res, "data.frame"), as(tax_table(kosticB)[rownames(res), ], "matrix"))
relabel <- paste(rownames(passing),passing$ta7,sep = "|")
res2 <- res
rownames(res2) <- relabel


shaping <- grepl("ITS",x = rownames(res2))

res2$shape <- shaping

p1 <- EnhancedVolcano(res2,
                title = "",subtitle = "", 
                axisLabSize = 10,titleLabSize = 8, legendLabSize = 7, 
                legendPosition = "right", 
                lab = "",
                transcriptPointSize = 5,
                x = 'log2FoldChange',
                y = 'padj',
                pCutoff = 0.05,
                FCcutoff = 2,
                drawConnectors = TRUE,
                xlab = bquote(~Log[2]~ 'fold change'),
                ylab = bquote(~-Log[10]~adjusted~italic(P)))

p1 + geom_point(size = 3, aes(shape = shape)) + scale_shape_manual(values=c(3,4)) + 
  annotate("text", x = -15, y = 25, label = "Enriched in Controls") + annotate("text", x = 15, y = 25, label = "Enriched in Melanoma")


## TEST ON INVASIVITY

qiime_file_16S_NOcss_join <- qiime_file_16S_initial
sample_names(qiime_file_16S_NOcss_join) <- sample_data(qiime_file_16S_NOcss_join)$Label
qiime_file_ITS_NOcss_join <- qiime_file_ITS_initial
sample_names(qiime_file_ITS_NOcss_join) <- sample_data(qiime_file_ITS_NOcss_join)$Label

sample_names(qiime_file_ITS_NOcss_join) == sample_names(qiime_file_ITS_NOcss_join)


otu_16S_join_NOcss <- otu_table(qiime_file_16S_NOcss_join)
row.names(otu_16S_join_NOcss) <- paste(row.names(otu_16S_join_NOcss),"16S", sep = "|")
otu_ITS_join_NOcss <- otu_table(qiime_file_ITS_NOcss_join)
row.names(otu_ITS_join_NOcss) <- paste(row.names(otu_ITS_join_NOcss),"ITS", sep = "|")

 tax table
tax_16S_join_NOcss <- tax_table(qiime_file_16S_NOcss_join)
row.names(tax_16S_join_NOcss) <- paste(row.names(tax_16S_join_NOcss),"16S", sep = "|")

tax_16S_join_NOcss <- as.data.frame(tax_16S_join_NOcss)
tax_16S_join_NOcss$Species <- tax_16S_join_NOcss$Genus
tax_ITS_join_NOcss <- tax_table(qiime_file_ITS_NOcss_join)
row.names(tax_ITS_join_NOcss) <- paste(row.names(tax_ITS_join_NOcss),"ITS", sep = "|")
tax_ITS_join_NOcss <- as.data.frame(tax_ITS_join_NOcss)

#
taxa_16SITS_NOcss <- rbind(tax_16S_join_NOcss, tax_ITS_join_NOcss)
otu_16SITS_NOcss <- rbind(otu_16S_join_NOcss, otu_ITS_join_NOcss)
str(taxa_16SITS_NOcss)
str(otu_16SITS_NOcss)

#
otus <- otu_table(otu_16SITS_NOcss,taxa_are_rows = T)
taxas <- tax_table(taxa_16SITS_NOcss)
row.names(taxas) <- row.names(otus)
phyloseq_16SITS_NOcss <- phyloseq(otu_table(otus), tax_table(taxas), sample_data(sample_data(qiime_file_16S_NOcss_join)))


mel_con_NOcss <- subset_samples(phyloseq_16SITS_NOcss, Pathology == "Melanoma" | Pathology == "Control")
leve <- c("Control", "Melanoma")
facto <- factor(sample_data(mel_con_NOcss)$Pathology, levels = leve)
sample_data(mel_con_NOcss)$Pathology <- facto
mel_con_noM6_NOcss <- subset_samples(mel_con_NOcss, Label != "M6F")

###

phylo_lista2 <- subset_samples(mel_con_noM6_NOcss, Pathology == "Melanoma")

sample_data(phylo_lista2)

# following https://bioconductor.org/packages/devel/bioc/vignettes/phyloseq/inst/doc/phyloseq-mixture-models.html

kosticB <- phylo_lista2
tax_table(kosticB)
sample_names(kosticB)

diagdds <- phyloseq_to_deseq2(kosticB, ~ Invasivity)
# calculate geometric means prior to estimate size factors
gm_mean <- function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
geoMeans <- apply(counts(diagdds), 1, gm_mean)
diagdds <- estimateSizeFactors(diagdds, geoMeans = geoMeans)
diagdds <- DESeq(diagdds, fitType="local",test = "Wald")

res <- results(diagdds,cooksCutoff = F,pAdjustMethod = "BH")
res <- res[order(res$padj, na.last=NA), ]
alpha <- 0.05
sigtab <- res[which(res$padj < alpha), ]
sigtab <- as(sigtab, "data.frame")
sigtax <- as.data.frame(tax_table(kosticB)@.Data)

rownames(sigtab) %in% sigtax[rownames(sigtab),]
sigtab = res[which(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(kosticB)[rownames(sigtab), ], "matrix"))
sigtab

list_InvasivevsInsitu <- row.names(sigtab)

# SUPPLEMENTARY TABLE 7
#write.table(x = sigtab, "./tabella_DESEQ_invasiveVSinsitu.csv", sep = ",")

# SUPPLEMENTARY FIGURE 4
library(EnhancedVolcano)
passing <- cbind(as(res, "data.frame"), as(tax_table(kosticB)[rownames(res), ], "matrix"))
relabel <- paste(rownames(passing),passing$ta7,sep = "|")
res2 <- res
rownames(res2) <- relabel

shaping <- grepl("ITS",x = rownames(res2))

res2$shape <- shaping

p1 <- EnhancedVolcano(res2,
                title = "",subtitle = "", 
                axisLabSize = 10,titleLabSize = 8, legendLabSize = 7, 
                legendPosition = "right", 
                lab = "",
                transcriptPointSize = 5,
                x = 'log2FoldChange',
                y = 'padj',
                pCutoff = 0.05,
                FCcutoff = 2,
                drawConnectors = TRUE,
                xlab = bquote(~Log[2]~ 'fold change'),
                ylab = bquote(~-Log[10]~adjusted~italic(P)))
p1 + geom_point(size = 3, aes(shape = shape)) + scale_shape_manual(values=c(3,4)) + 
  annotate("text", x = -15, y = 20, label = "Enriched in In situ") + annotate("text", x = 15, y = 20, label = "Enriched in Invasive")


## TEST ON REGRESSION

phylo_lista2 <- subset_samples(mel_con_noM6_NOcss, Pathology == "Melanoma")

sample_data(phylo_lista2)

# following https://bioconductor.org/packages/devel/bioc/vignettes/phyloseq/inst/doc/phyloseq-mixture-models.html

kosticB <- phylo_lista2
tax_table(kosticB)
sample_names(kosticB)
sample_data(kosticB)

diagdds <- phyloseq_to_deseq2(kosticB, ~ Regression_generic)
# calculate geometric means prior to estimate size factors
gm_mean <- function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
geoMeans <- apply(counts(diagdds), 1, gm_mean)
diagdds <- estimateSizeFactors(diagdds, geoMeans = geoMeans)
diagdds <- DESeq(diagdds, fitType="local",test = "Wald")

res <- results(diagdds,cooksCutoff = F,pAdjustMethod = "BH")
res <- res[order(res$padj, na.last=NA), ]
alpha <- 0.05
sigtab <- res[which(res$padj < alpha), ]
sigtab <- as(sigtab, "data.frame")
sigtax <- as.data.frame(tax_table(kosticB)@.Data)

rownames(sigtab) %in% sigtax[rownames(sigtab),]
sigtab = res[which(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(kosticB)[rownames(sigtab), ], "matrix"))
sigtab

list_RegressedvsNonregressed <- row.names(sigtab)

# SUPPLEMENTARY TABLE 6
#write.table(x = sigtab, "./tabella_DESEQ_regressionVSnoregression.csv", sep = ",")


# SUPPLEMENTARY FIGURE 4
library(EnhancedVolcano)
passing <- cbind(as(res, "data.frame"), as(tax_table(kosticB)[rownames(res), ], "matrix"))
relabel <- paste(rownames(passing),passing$ta7,sep = "|")
res2 <- res
rownames(res2) <- relabel


shaping <- grepl("ITS",x = rownames(res2))

res2$shape <- shaping

p1 <- EnhancedVolcano(res2,
                      title = "",subtitle = "", 
                      axisLabSize = 10,titleLabSize = 8, legendLabSize = 7, 
                      legendPosition = "right", 
                      lab = "",
                      transcriptPointSize = 5,
                      x = 'log2FoldChange',
                      y = 'padj',
                      pCutoff = 0.05,
                      FCcutoff = 2,
                      drawConnectors = TRUE,
                      xlab = bquote(~Log[2]~ 'fold change'),
                      ylab = bquote(~-Log[10]~adjusted~italic(P)))
p1 + geom_point(size = 3, aes(shape = shape)) + scale_shape_manual(values=c(3,4)) + 
  annotate("text", x = -15, y = 22.5, label = "Not regressed") + annotate("text", x = 15, y = 22.5, label = "Regressed")



            
                     
#### BUILDING G PANEL

set1 <- list_InvasivevsInsitu
set2 <- list_CvsM
set3 <- list_RegressedvsNonregressed

# library

library(eulerr)
vd <- euler(c(Invasivity = 171, C_vs_M = 14, Regression = 12,
              "Invasivity&C_vs_M" = 1, "Invasivity&Regression" = 7, "C_vs_M&Regression" = 3,
              "Invasivity&C_vs_M&Regression" = 1), shape = "ellipse")
plot(vd, legend = T, quantities = T, key = TRUE, counts = TRUE,
     edges = list(lty = 1),
     labels = list(font = 2, cex = 1),
     fills = list(fill = c("#f7f7f7", "#bdbdbd", "#636363")))




#######################
## BUILDING TABLE 1 ##

mel_nom6 <- subset_samples(mel_con_noM6, Pathology != "Control")
mel_nom6_bact <- subset_samples(mel_con_noM6_bact, Pathology != "Control")
mel_nom6_fung <- subset_samples(mel_con_noM6_fung, Pathology != "Control")

#####################################
## test Bacteria and Fungi dataset ##
datas <- as(sample_data(mel_nom6), "data.frame")
datas <- datas[,12:length(colnames(datas))]
# remove inflammation as NA breaks the loop
datas <- datas[,-4]
variables <- colnames(datas)
result_bact_fung <- data.frame(matrix(nrow = length(variables), ncol = 4))
colnames(result_bact_fung) <- c("R_adonis","p_adonis","R_anosim","p_anosim")
rownames(result_bact_fung) <- colnames(datas)

options(digits=3)
for (i in 1:length(variables)) {
  set.seed(12387)
  ado <- adonis2(distance(mel_nom6, method = "bray") ~ datas[,i], permutations = 9999)
  set.seed(12387)
  test_group  <- get_variable(mel_nom6, variables[i])
  test_anosim <-  anosim(phyloseq::distance(mel_nom6, "bray"), test_group, permutations = 9999)
  result_bact_fung[i, 1] <- ado$R2[1]
  result_bact_fung[i, 2] <- ado$`Pr(>F)`[1]
  result_bact_fung[i, 3] <- test_anosim$statistic
  result_bact_fung[i, 4] <- test_anosim$signif
}  

result_bact_fung

#####################################
## test Bacteria  dataset ##

datas <- as(sample_data(mel_nom6_bact), "data.frame")
datas <- datas[,12:length(colnames(datas))]
# remove inflammation as NA breaks the loop
datas <- datas[,-4]
variables <- colnames(datas)
result_bact <- data.frame(matrix(nrow = length(variables), ncol = 4))
colnames(result_bact) <- c("R_adonis","p_adonis","R_anosim","p_anosim")
rownames(result_bact) <- colnames(datas)

options(digits=4)
for (i in 1:length(variables)) {
  set.seed(12387)
  ado <- adonis2(distance(mel_nom6_bact, method = "bray") ~ datas[,i], permutations = 9999)
  set.seed(12387)
  test_group  <- get_variable(mel_nom6_bact, variables[i])
  test_anosim <-  anosim(phyloseq::distance(mel_nom6_bact, "bray"), test_group, permutations = 9999)
  result_bact[i, 1] <- ado$R2[1]
  result_bact[i, 2] <- ado$`Pr(>F)`[1]
  result_bact[i, 3] <- test_anosim$statistic
  result_bact[i, 4] <- test_anosim$signif
}  

result_bact

#####################################
## test fungi  dataset ##

datas <- as(sample_data(mel_nom6_fung), "data.frame")
datas <- datas[,12:length(colnames(datas))]
# remove inflammation as NA breaks the loop
datas <- datas[,-4]
variables <- colnames(datas)
result_fung <- data.frame(matrix(nrow = length(variables), ncol = 4))
colnames(result_fung) <- c("R_adonis","p_adonis","R_anosim","p_anosim")
rownames(result_fung) <- colnames(datas)

options(digits=4)
for (i in 1:length(variables)) {
  set.seed(12387)
  ado <- adonis2(distance(mel_nom6_fung, method = "bray") ~ datas[,i], permutations = 9999)
  set.seed(12387)
  test_group  <- get_variable(mel_nom6_fung, variables[i])
  test_anosim <-  anosim(phyloseq::distance(mel_nom6, "bray"), test_group, permutations = 9999)
  result_fung[i, 1] <- ado$R2[1]
  result_fung[i, 2] <- ado$`Pr(>F)`[1]
  result_fung[i, 3] <- test_anosim$statistic
  result_fung[i, 4] <- test_anosim$signif
}  

result_fung

##### TEST ON INVASIVITY : FIGURE 3 PANEL A, B, C

# INVASIVITY on total community

sample_data(mel_con_noM6)
mel_nom6 <- subset_samples(mel_con_noM6, Pathology != "Control")

summary(sample_data(mel_nom6)$Invasivity)
summary(sample_data(mel_nom6)$Pathology)

ordination_fungBact_melanoma <- ordinate(mel_nom6, "PCoA", distance = "bray")
p_ordination_fungBact <- plot_ordination(physeq = mel_nom6, ordination_fungBact_melanoma, axes = c(1,2), color="Invasivity")
p12 <- p_ordination_fungBact + geom_point(size = 5, shape = 20) + geom_text(aes(label = Label), size = 3, vjust = 2)
p12 <- p12 + theme_bw() + 
  theme( panel.grid.major =  element_blank(),
         panel.grid.minor =  element_blank())
p12 <- p12 + scale_color_aaas() + xlab("PCoA 1 (13.2%)") +
  ylab("PCoA 2 (9.6%)")
ordination_fungBact_melanoma <- ordinate(mel_nom6, "PCoA", distance = "bray")
p_ordination_fungBact <- plot_ordination(physeq = mel_nom6, ordination_fungBact_melanoma, axes = c(1,3), color="Invasivity")
p13 <- p_ordination_fungBact + geom_point(size = 5, shape = 20) + geom_text(aes(label = Label), size = 3, vjust = 2)
p13 <- p13 + theme_bw() + 
  theme( panel.grid.major =  element_blank(),
         panel.grid.minor =  element_blank())
p13 <- p13 + scale_color_aaas() + xlab("PCoA 1 (13.2%)") +
  ylab("PCoA 3 (8.2%)")

plot_grid(p12, p13, labels = c('A', 'B'), label_size = 12)

#PERMANOVA
set.seed(12387)
datas <- as(sample_data(mel_nom6), "data.frame")
adonis2(distance(mel_nom6, method = "bray") ~ Invasivity, data = datas, permutations = 9999)
#ANOSIM
set.seed(12387)
df <- as(sample_data(mel_nom6), "data.frame")
test_group  <- get_variable(mel_nom6, "Invasivity")
test_anosim <-  anosim(phyloseq::distance(mel_nom6, "bray"), test_group, permutations = 9999)
summary(test_anosim)
plot(test_anosim)

############## #######################
## PLOT ALPHADIV BETWEEN INVASIVITY ##

#RICHNESS of bacteria
alphadiv <- microbiome::alpha(subset_samples(qiime_file_16S_css, Pathology == "Melanoma" & Label != "M6F"), index = "observed")
df <- as(sample_data(subset_samples(qiime_file_16S_css, Pathology == "Melanoma" & Label != "M6F")), "data.frame")
df$rich <- alphadiv$observed
df$Invasivity
# Plot on Invasivity
variabile_plot = "Invasivity"
my_comparisons <- list(c("in_situ","invasive"))

p_rich <- ggdotplot(data = df, x = variabile_plot, y = "rich", add = "boxplot", size = 0.8, fill = variabile_plot)
p_rich <- p_rich + stat_compare_means(method = "wilcox.test", comparisons = my_comparisons)                   # Add global p-value
p_rich <- p_rich + ylab("Richness (n° of OTUs)") + scale_fill_aaas()+ xlab("")+ theme(legend.position="none")
p_rich_bact <- p_rich

#RICHNESS of fungi
alphadiv <- microbiome::alpha(subset_samples(qiime_file_ITS_css, Pathology == "Melanoma" & Label != "M6F"), index = "observed")
df <- as(sample_data(subset_samples(qiime_file_ITS_css, Pathology == "Melanoma" & Label != "M6F")), "data.frame")
df$rich <- alphadiv$observed
df$Invasivity
# Plot on Invasivity
variabile_plot = "Invasivity"
my_comparisons <- list(c("in_situ","invasive"))

p_rich <- ggdotplot(data = df, x = variabile_plot, y = "rich", add = "boxplot", size = 0.8, fill = variabile_plot)
p_rich <- p_rich + stat_compare_means(method = "wilcox.test", comparisons = my_comparisons)                   # Add global p-value
p_rich <- p_rich + ylab("Richness (n° of OTUs)") + scale_fill_aaas() + xlab("")+theme(legend.position="none")
p_rich_fung <- p_rich

plot_grid(p_rich_bact, p_rich_fung, labels = c('A', 'B'), label_size = 12) + 
  geom_point(size = 6)

#########

#SHANNON of bacteria
alphadiv <- microbiome::alpha(subset_samples(qiime_file_16S_css, Pathology == "Melanoma" & Label != "M6F"), index = "shannon")
df <- as(sample_data(subset_samples(qiime_file_16S_css, Pathology == "Melanoma" & Label != "M6F")), "data.frame")
df$shan<- alphadiv$diversity_shannon
df$Invasivity

df[,c("Invasivity", "shan")]
# Plot on Invasivity
variabile_plot = "Invasivity"
my_comparisons <- list(c("in_situ","invasive"))

p_shan <- ggdotplot(data = df, x = variabile_plot, y = "shan", add = "boxplot", size = 0.8, fill = variabile_plot)
p_shan <- p_shan + stat_compare_means(method = "wilcox.test", comparisons = my_comparisons)                   # Add global p-value
p_shan <- p_shan + ylab("Shannon's Index") + scale_fill_aaas()+ xlab("")+theme(legend.position="none")
p_shan_bact <- p_shan

#SHANNON of fungi
alphadiv <- microbiome::alpha(subset_samples(qiime_file_ITS_css, Pathology == "Melanoma" & Label != "M6F"), index = "shannon")
df <- as(sample_data(subset_samples(qiime_file_ITS_css, Pathology == "Melanoma" & Label != "M6F")), "data.frame")
df$shan <- alphadiv$diversity_shannon
df$Invasivity
# Plot on Invasivity
variabile_plot = "Invasivity"
my_comparisons <- list(c("in_situ","invasive"))

p_shan <- ggdotplot(data = df, x = variabile_plot, y = "shan", add = "boxplot", size = 0.8, fill = variabile_plot)
p_shan <- p_shan + stat_compare_means(method = "wilcox.test", comparisons = my_comparisons)                   # Add global p-value
p_shan <- p_shan + ylab("Richness (n° of OTUs)") + scale_fill_aaas() + xlab("")+theme(legend.position="none")
p_shan_fung <- p_shan

plot_grid(p_rich_bact, p_rich_fung, p_shan_bact,p_shan_fung, nrow = 2, label_size = 12) + 
  geom_point(size = 6)


## REMOVE OUTLIER

#RICHNESS of bacteria
alphadiv <- microbiome::alpha(subset_samples(qiime_file_16S_css, Pathology == "Melanoma" & Label != "M6F" & Label != "M7F"), index = "observed")
df <- as(sample_data(subset_samples(qiime_file_16S_css, Pathology == "Melanoma" & Label != "M6F" & Label != "M7F")), "data.frame")
df$rich <- alphadiv$observed
df$Invasivity

# Plot on Invasivity
variabile_plot = "Invasivity"
my_comparisons <- list(c("in_situ","invasive"))

p_rich <- ggdotplot(data = df, x = variabile_plot, y = "rich", add = "boxplot", size = 0.8, fill = variabile_plot)
p_rich <- p_rich + stat_compare_means(method = "wilcox.test", comparisons = my_comparisons)                   # Add global p-value
p_rich <- p_rich + ylab("Richness (n° of OTUs)") + scale_fill_aaas()+ xlab("")+theme(legend.position="none")
p_rich_bact <- p_rich

#RICHNESS of fungi
alphadiv <- microbiome::alpha(subset_samples(qiime_file_ITS_css, Pathology == "Melanoma" & Label != "M6F" & Label != "M7F"), index = "observed")
df <- as(sample_data(subset_samples(qiime_file_ITS_css, Pathology == "Melanoma" & Label != "M6F" & Label != "M7F")), "data.frame")
df$rich <- alphadiv$observed
df$Invasivity
# Plot on Invasivity
variabile_plot = "Invasivity"
my_comparisons <- list(c("in_situ","invasive"))

p_rich <- ggdotplot(data = df, x = variabile_plot, y = "rich", add = "boxplot", size = 0.8, fill = variabile_plot)
p_rich <- p_rich + stat_compare_means(method = "wilcox.test", comparisons = my_comparisons)                   # Add global p-value
p_rich <- p_rich + ylab("Richness (n° of OTUs)") + scale_fill_aaas() + xlab("")+theme(legend.position="none")
p_rich_fung <- p_rich

#plot_grid(p_rich_bact, p_rich_fung, labels = c('A', 'B'), label_size = 12) + 
  geom_point(size = 6)

#########

#SHANNON of bacteria
alphadiv <- microbiome::alpha(subset_samples(qiime_file_16S_css, Pathology == "Melanoma" & Label != "M6F" & Label != "M7F"), index = "shannon")
df <- as(sample_data(subset_samples(qiime_file_16S_css, Pathology == "Melanoma" & Label != "M6F" & Label != "M7F")), "data.frame")
df$shan<- alphadiv$diversity_shannon
df$Invasivity
# Plot on Invasivity
variabile_plot = "Invasivity"
my_comparisons <- list(c("in_situ","invasive"))

p_shan <- ggdotplot(data = df, x = variabile_plot, y = "shan", add = "boxplot", size = 0.8, fill = variabile_plot)
p_shan <- p_shan + stat_compare_means(method = "wilcox.test", comparisons = my_comparisons)                   # Add global p-value
p_shan <- p_shan + ylab("Shannon's Index") + scale_fill_aaas()+ xlab("")+theme(legend.position="none")
p_shan_bact <- p_shan

#SHANNON of fungi
alphadiv <- microbiome::alpha(subset_samples(qiime_file_ITS_css, Pathology == "Melanoma" & Label != "M6F" & Label != "M7F"), index = "shannon")
df <- as(sample_data(subset_samples(qiime_file_ITS_css, Pathology == "Melanoma" & Label != "M6F" & Label != "M7F")), "data.frame")
df$shan <- alphadiv$diversity_shannon
df$Invasivity
# Plot on Invasivity
variabile_plot = "Invasivity"
my_comparisons <- list(c("in_situ","invasive"))

p_shan <- ggdotplot(data = df, x = variabile_plot, y = "shan", add = "boxplot", size = 0.8, fill = variabile_plot)
p_shan <- p_shan + stat_compare_means(method = "wilcox.test", comparisons = my_comparisons)                   # Add global p-value
p_shan <- p_shan + ylab("Richness (n° of OTUs)") + scale_fill_aaas() + xlab("")+theme(legend.position="none")
p_shan_fung <- p_shan

plot_grid(p_rich_bact, p_rich_fung, p_shan_bact,p_shan_fung, nrow = 2, label_size = 12) + 
  geom_point(size = 6)

## TEST ON REGRESSION : FIGURE 3 PANEL E, F

# REGRESSION on fungi community
summary(sample_data(mel_con_noM6_fung)$Pathology)

ordination_fungBact_melanoma <- ordinate(subset_samples(mel_con_noM6_fung, Pathology != "Control"), "PCoA", distance = "bray")
p_ordination_fungBact <- plot_ordination(physeq = subset_samples(mel_con_noM6_fung, Pathology != "Control"), ordination_fungBact_melanoma, axes = c(1,2), color="Regression_generic")
p12 <- p_ordination_fungBact + geom_point(size = 5, shape = 20) + geom_text(aes(label = Label), size = 3, vjust = 2)
p12 <- p12 + theme_bw() + 
  theme( panel.grid.major =  element_blank(),
         panel.grid.minor =  element_blank())
p12 <- p12 + scale_color_aaas() + xlab("PCoA 1 (16.3%)") +
  ylab("PCoA 2 (13.2%)")
ordination_fungBact_melanoma <- ordinate(subset_samples(mel_con_noM6_fung, Pathology != "Control"), "PCoA", distance = "bray")
p_ordination_fungBact <- plot_ordination(physeq = subset_samples(mel_con_noM6_fung, Pathology != "Control"), ordination_fungBact_melanoma, axes = c(1,3), color="Regression_generic")
p13 <- p_ordination_fungBact + geom_point(size = 5, shape = 20) + geom_text(aes(label = Label), size = 3, vjust = 2)
p13 <- p13 + theme_bw() + 
  theme( panel.grid.major =  element_blank(),
         panel.grid.minor =  element_blank())
p13 <- p13 + scale_color_aaas() + xlab("PCoA 1 (16.3%)") +
  ylab("PCoA 3 (10.9%)")

plot_grid(p12, p13, labels = c('A', 'B'), label_size = 12)

#ANOSIM
set.seed(12387)
df <- as(sample_data(subset_samples(mel_con_noM6_fung, Pathology != "Control")), "data.frame")
test_group  <- get_variable(subset_samples(mel_con_noM6_fung, Pathology != "Control"), "Invasivity")
test_anosim <-  anosim(phyloseq::distance(subset_samples(mel_con_noM6_fung, Pathology != "Control"), "bray"), test_group, permutations = 9999)
summary(test_anosim)
plot(test_anosim)

summary(df$Regression_generic)

#PERMANOVA
set.seed(12387)
datas <- as(sample_data(subset_samples(mel_con_noM6_fung, Pathology != "Control")), "data.frame")
adonis2(distance(subset_samples(mel_con_noM6_fung, Pathology != "Control"), method = "bray") ~ Invasivity, data = datas, permutations = 9999)


#RICHNESS of bacteria
alphadiv <- microbiome::alpha(subset_samples(qiime_file_16S_css, Pathology == "Melanoma" & Label != "M6F"), index = "observed")
df <- as(sample_data(subset_samples(qiime_file_16S_css, Pathology == "Melanoma" & Label != "M6F")), "data.frame")
df$rich <- alphadiv$observed
df$Regression_generic
# Plot on Regression
variabile_plot = "Regression_generic"
my_comparisons <- list(c("yes","no"))

p_rich <- ggdotplot(data = df, x = variabile_plot, y = "rich", add = "boxplot", size = 0.8, fill = variabile_plot)
p_rich <- p_rich + stat_compare_means(method = "wilcox.test", comparisons = my_comparisons)                   # Add global p-value
p_rich <- p_rich + ylab("Richness (n° of OTUs)") + scale_fill_aaas()
p_rich

#RICHNESS of fungi
alphadiv <- microbiome::alpha(subset_samples(qiime_file_ITS_css, Pathology == "Melanoma" & Label != "M6F"), index = "observed")
df <- as(sample_data(subset_samples(qiime_file_ITS_css, Pathology == "Melanoma" & Label != "M6F")), "data.frame")
df$rich <- alphadiv$observed
df$Invasivity
# Plot on Regression
variabile_plot = "Regression_generic"
my_comparisons <- list(c("yes","no"))

p_rich <- ggdotplot(data = df, x = variabile_plot, y = "rich", add = "boxplot", size = 0.8, fill = variabile_plot)
p_rich <- p_rich + stat_compare_means(method = "wilcox.test", comparisons = my_comparisons)                   # Add global p-value
p_rich <- p_rich + ylab("Richness (n° of OTUs)") + scale_fill_aaas()
p_rich


# REGRESSION on total community

ordination_fungBact_melanoma <- ordinate(mel_nom6, "PCoA", distance = "bray")
p_ordination_fungBact <- plot_ordination(physeq = mel_nom6, ordination_fungBact_melanoma, axes = c(1,2), color="Regression_generic")
p12 <- p_ordination_fungBact + geom_point(size = 5, shape = 20) + geom_text(aes(label = Label), size = 3, vjust = 2)
p12 <- p12 + theme_bw() + 
  theme( panel.grid.major =  element_blank(),
         panel.grid.minor =  element_blank())
p12 <- p12 + scale_color_aaas() + xlab("PCoA 1 (13.2%)") +
  ylab("PCoA 2 (9.6%)")
ordination_fungBact_melanoma <- ordinate(mel_nom6, "PCoA", distance = "bray")
p_ordination_fungBact <- plot_ordination(physeq = mel_nom6, ordination_fungBact_melanoma, axes = c(1,3), color="Regression_generic")
p13 <- p_ordination_fungBact + geom_point(size = 5, shape = 20) + geom_text(aes(label = Label), size = 3, vjust = 2)
p13 <- p13 + theme_bw() + 
  theme( panel.grid.major =  element_blank(),
         panel.grid.minor =  element_blank())
p13 <- p13 + scale_color_aaas() + xlab("PCoA 1 (13.2%)") +
  ylab("PCoA 2 (8.2%)")

plot_grid(p12, p13, labels = c('A', 'B'), label_size = 12)

#PERMANOVA
datas <- as(sample_data(mel_nom6), "data.frame")
adonis2(distance(mel_nom6, method = "bray") ~ Regression_generic, data = datas, permutations = 9999)
#ANOSIM
set.seed(12387)
df <- as(sample_data(mel_nom6), "data.frame")
test_group  <- get_variable(mel_nom6, "Regression_generic")
test_anosim <-  anosim(phyloseq::distance(mel_nom6, "bray"), test_group, permutations = 9999)
summary(test_anosim)
plot(test_anosim)

##################################
######## BUILDING FIGURE 4 #######
##################################

#### RANDOM FOREST
## following for reference https://github.com/LangilleLab/microbiome_helper/wiki/Random-Forest-Tutorial

remove_rare <- function( table , cutoff_pro ) {
  row2keep <- c()
  cutoff <- ceiling( cutoff_pro * ncol(table) )  
  for ( i in 1:nrow(table) ) {
    row_nonzero <- length( which( table[ i , ]  > 0 ) ) 
    if ( row_nonzero > cutoff ) {
      row2keep <- c( row2keep , i)
    }
  }
  return( table [ row2keep , , drop=F ])
}

## following for reference https://github.com/LangilleLab/microbiome_helper/wiki/Random-Forest-Tutorial

# set data; same obj as ordination analysis

taxxx <- tax_table(mel_nom6)@.Data

head(tax_table(mel_con_noM6))
tail(tax_table(mel_con_noM6))


phylo_RF_BF <- mel_nom6
sample_names(phylo_RF_BF)
otu_RF_BF <- as.data.frame(otu_table(phylo_RF_BF))
metadata_RF_BF <- as(sample_data(phylo_RF_BF), "data.frame")
summary(metadata_RF_BF$Pathology)
dim(otu_RF_BF)
dim(metadata_RF_BF)

## PRE-PROCESSING

#viz
otu_nonzero_counts <- apply(otu_RF_BF, 1, function(y) sum(length(which(y > 0))))
hist(otu_nonzero_counts, breaks=100, col="grey", main="", ylab="Number of OTUs", xlab="Number of Non-Zero Values", xlim= c(0,20), xaxt="n")
axis(1, at = seq(0, 20, by = 1), las=1)

#Typically researchers discard OTUs that are zero in greater 
#than 75-90% of samples although these cut-offs are somewhat arbitrary. 

# remove rare
otu_table_rare_removed_BF <- remove_rare(table=otu_RF_BF, cutoff_pro=0.15) 

#viz
otu_nonzero_counts <- apply(otu_table_rare_removed_BF, 1, function(y) sum(length(which(y > 0))))
hist(otu_nonzero_counts, breaks=100, col="grey", main="", ylab="Number of OTUs", xlab="Number of Non-Zero Values", xlim= c(0,20), xaxt="n")
axis(1, at = seq(0, 20, by = 1), las=1)

otu_table_rare_removed_BF <- as.data.frame(t(otu_table_rare_removed_BF))
dim(otu_table_rare_removed_BF)
tail(otu_table_rare_removed_BF)
str(otu_table_rare_removed_BF)

#### RF MODEL ON INVASIVITY
# get the Df in the right format
otu_table_rare_removed_BF$Invasivity <- metadata_RF_BF$Invasivity
set.seed(151)  
RF_invasivity_classify_BF <- randomForest(x=otu_table_rare_removed_BF[,1:(ncol(otu_table_rare_removed_BF)-1)] , 
                                          y=otu_table_rare_removed_BF[ , ncol(otu_table_rare_removed_BF)] , 
                                          ntree=1001, importance=TRUE, proximities=TRUE )
RF_invasivity_classify_BF
plot(RF_invasivity_classify_BF)

library(rfUtilities)
RF_invasivity_classify_BF.sig <- rf.significance( x=RF_invasivity_classify_BF ,  
                                                  xdata=otu_table_rare_removed_BF[,1:(ncol(otu_table_rare_removed_BF)-1)] ,
                                                  nperm=1000 , ntree=1001 )  
RF_invasivity_classify_BF.sig
plot(RF_invasivity_classify_BF.sig)

# following for reference https://stackoverflow.com/questions/30366143/how-to-compute-roc-and-auc-under-roc-after-training-using-caret-in-r
library(ROCR)
predictions=as.vector(RF_invasivity_classify_BF$votes[,2])
pred=prediction(predictions, otu_table_rare_removed_BF[ , ncol(otu_table_rare_removed_BF)])
perf_AUC=performance(pred,"auc") #Calculate the AUC value
performance(pred,"err")
AUC=perf_AUC@y.values[[1]]
perf_ROC=performance(pred,"tpr","fpr") #plot the actual ROC curve
plot(perf_ROC,col=2,lwd=2)
abline(a=0,b=1,lwd=2,lty=2,col="gray")
text(0.2,1,paste("AUC = ",format(AUC, digits=3, scientific=FALSE)))
text(0.2,0.8,paste("OOB = ",format(RF_invasivity_classify_BF.sig$test.OOB, digits=3, scientific=FALSE)))
text(0.2,0.9,paste("p-val = ",format(RF_invasivity_classify_BF.sig$pValue, digits=3, scientific=FALSE)))


#explore RF
library(randomForestExplainer)

explain_forest(RF_invasivity_classify_BF, interactions = FALSE, data = NULL,vars = NULL, no_of_pred_plots = 3, pred_grid = 100,measures = NULL)

RF_importance <- measure_importance(RF_invasivity_classify_BF)
important_variables(RF_importance,measures = c("accuracy_decrease", "gini_decrease"),
                    k = 10, ties_action = "all")

importance <- plot_multi_way_importance(RF_importance, 
                                        x_measure = "accuracy_decrease", 
                                        y = "gini_decrease")
importance + xlab("Accuracy Decrease") + 
  ylab("Gini Decrease") + theme_linedraw()

plot_importance_rankings(RF_importance, main = "Relations between rankings according to different measures")


dieffe_plot_RF <-data.frame(OTU1559 = otu_table_rare_removed_BF[,"DENOVO1559|16S", drop = F],
                            OTU231 = otu_table_rare_removed_BF[,"DENOVO231|16S", drop = F],
                            OTU530 = otu_table_rare_removed_BF[,"DENOVO530|16S", drop = F],
                            OTU392 = otu_table_rare_removed_BF[,"DENOVO392|16S", drop = F],
                            OTU2266 = otu_table_rare_removed_BF[,"DENOVO2266|16S", drop = F],
                            OTU232 = otu_table_rare_removed_BF[,"DENOVO232|16S", drop = F],
                            OTU1804 = otu_table_rare_removed_BF[,"DENOVO1804|16S", drop = F],
                            OTU240 = otu_table_rare_removed_BF[,"DENOVO240|16S", drop = F],
                            OTU517 = otu_table_rare_removed_BF[,"DENOVO517|16S", drop = F],
                            OTU197 = otu_table_rare_removed_BF[,"DENOVO197|16S", drop = F],
                            Invasivity = otu_table_rare_removed_BF$Invasivity,
                            ID = row.names(otu_table_rare_removed_BF))

dieffe_plot_RF_melt <- melt(dieffe_plot_RF, c("Invasivity","ID"))
colnames(dieffe_plot_RF_melt) <- c("Invasivity","ID","OTU","Abundance")

ggplot(data =dieffe_plot_RF_melt, aes(x = OTU, y = Abundance, fill = Invasivity, group = Invasivity) ) + 
  geom_bar(stat = "identity", position=position_dodge()) + theme_linedraw() + coord_flip() + scale_x_discrete(limits = rev(levels(dieffe_plot_RF_melt$OTU)))



##################################
######## BUILDING FIGURE 6 #######
##################################


qiime_file_16S_NOcss_join <- qiime_file_16S_initial
sample_names(qiime_file_16S_NOcss_join) <- sample_data(qiime_file_16S_NOcss_join)$Label
qiime_file_ITS_NOcss_join <- qiime_file_ITS_initial
sample_names(qiime_file_ITS_NOcss_join) <- sample_data(qiime_file_ITS_NOcss_join)$Label

sample_names(qiime_file_ITS_NOcss_join) == sample_names(qiime_file_ITS_NOcss_join)

otu_16S_join_NOcss <- otu_table(qiime_file_16S_NOcss_join)
row.names(otu_16S_join_NOcss) <- paste(row.names(otu_16S_join_NOcss),"16S", sep = "|")
otu_ITS_join_NOcss <- otu_table(qiime_file_ITS_NOcss_join)
row.names(otu_ITS_join_NOcss) <- paste(row.names(otu_ITS_join_NOcss),"ITS", sep = "|")

tax_16S_join_NOcss <- tax_table(qiime_file_16S_NOcss_join)
row.names(tax_16S_join_NOcss) <- paste(row.names(tax_16S_join_NOcss),"16S", sep = "|")
tax_16S_join_NOcss <- as.data.frame(tax_16S_join_NOcss)
tax_16S_join_NOcss$Species <- tax_16S_join_NOcss$Genus
tax_ITS_join_NOcss <- tax_table(qiime_file_ITS_NOcss_join)
row.names(tax_ITS_join_NOcss) <- paste(row.names(tax_ITS_join_NOcss),"ITS", sep = "|")
tax_ITS_join_NOcss <- as.data.frame(tax_ITS_join_NOcss)

taxa_16SITS_NOcss <- rbind(tax_16S_join_NOcss, tax_ITS_join_NOcss)
otu_16SITS_NOcss <- rbind(otu_16S_join_NOcss, otu_ITS_join_NOcss)
str(taxa_16SITS_NOcss)
str(otu_16SITS_NOcss)

otus <- otu_table(otu_16SITS_NOcss,taxa_are_rows = T)
taxas <- tax_table(taxa_16SITS_NOcss)
row.names(taxas) <- row.names(otus)
phyloseq_16SITS_NOcss <- phyloseq(otu_table(otus), tax_table(taxas), sample_data(sample_data(qiime_file_16S_NOcss_join)))

sample_data(phyloseq_16SITS_NOcss)$Pathology
mel_ML_NOcss <- subset_samples(phyloseq_16SITS_NOcss, Pathology == "Melanoma" | Pathology == "Leucoderma_Melanoma")
leve <- c("Leucoderma_Melanoma", "Melanoma")
facto <- factor(sample_data(mel_ML_NOcss)$Pathology, levels = leve)
sample_data(mel_ML_NOcss)$Pathology <- facto
mel_ML_NOm6_NOcss <- subset_samples(mel_ML_NOcss, Label != "M6F")


# ORDINATION

#SCALE CSS

summary(taxa_sums(mel_ML_NOm6_NOcss))
boxplot(taxa_sums(mel_ML_NOm6_NOcss))
doubleton <- genefilter_sample(mel_ML_NOm6_NOcss, filterfun_sample(function(x) x > 1), A=1)
doubleton <- prune_taxa(doubleton, mel_ML_NOm6_NOcss) 
summary(taxa_sums(doubleton))

data.metagenomeSeq = phyloseq_to_metagenomeSeq(doubleton)
p = cumNormStat(data.metagenomeSeq) #default is 0.5
data.cumnorm = cumNorm(data.metagenomeSeq, p=p)
data.CSS = MRcounts(data.cumnorm, norm=TRUE, log=TRUE)
dim(data.CSS)  # make sure the data are in a correct formal: number of samples in rows
MAL_M_css <- mel_ML_NOm6_NOcss
otu_table(MAL_M_css) <- otu_table(data.CSS, taxa_are_rows = T)


sample_data(MAL_M_css)

ordination_bact <- ordinate(MAL_M_css, "PCoA", distance = "bray")
p_ordination_bact <- plot_ordination(physeq = MAL_M_css, ordination_bact, axes = c(1,2), color="Pathology", shape = "Invasivity")
p12 <- p_ordination_bact + geom_point(size = 5) 
p12 <- p12 + theme_bw() + 
  theme( panel.grid.major =  element_blank(),
         panel.grid.minor =  element_blank())
p12 <- p12 + scale_color_aaas() + xlab("PCoA 1 (12.1%)") +
  ylab("PCoA 2 (9%)")
p12 + theme(legend.position="none") #+ stat_ellipse(type = "t", linetype = 2) 

set.seed(1233)
datas <- as(sample_data(MAL_M_css), "data.frame")
adonis2(distance(MAL_M_css, "bray") ~ Pathology * Invasivity, data = datas, permutations = 9999)


# following for reference https://bioconductor.org/packages/devel/bioc/vignettes/phyloseq/inst/doc/phyloseq-mixture-models.html

kosticB <- tax_glom(mel_ML_NOm6_NOcss, taxrank = "ta7")
tax_table(kosticB)
sample_names(kosticB)

diagdds <- phyloseq_to_deseq2(kosticB, ~ Pathology)
# calculate geometric means prior to estimate size factors
gm_mean <- function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
geoMeans <- apply(counts(diagdds), 1, gm_mean)
diagdds <- estimateSizeFactors(diagdds, geoMeans = geoMeans)
diagdds <- DESeq(diagdds, fitType="local",test = "Wald")

res <- results(diagdds,cooksCutoff = F,pAdjustMethod = "BH")
res <- res[order(res$padj, na.last=NA), ]
alpha <- 0.05
sigtab <- res[which(res$padj < alpha), ]
sigtab <- as(sigtab, "data.frame")
sigtax <- as.data.frame(tax_table(kosticB)@.Data)

rownames(sigtab) %in% sigtax[rownames(sigtab),]
sigtab = res[which(res$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(kosticB)[rownames(sigtab), ], "matrix"))
sigtab

df <-sigtab[order(sigtab$log2FoldChange),]
df$ta7 <- as.character(df$ta7)

df$ta7

##
ggdotchart(df, x = "ta7", y = "log2FoldChange", add = "segment", 
           color = "ta1", rotate = T, dot.size = 3) + 
  geom_hline(yintercept = 0,linetype="dashed") + 
  theme(legend.position="none") +
  annotate("text", y = 15, x = 15, label = "Enriched in M")+
  annotate("text", y = -20, x = 1, label = "Enriched in MAL")

# Define the taxa you don't want like this:
badTaxa = row.names(df)
allTaxa = taxa_names(kosticB)
allTaxa <- allTaxa[allTaxa %in% badTaxa]
ex1 = prune_taxa(allTaxa, transform(kosticB, "identity"))
# new phyloseq object with just the taxa you kept.
confronto_M_MAL <- ex1

sample_names(confronto_M_MAL)

taxa_names(confronto_M_MAL) <- as.data.frame(tax_table(confronto_M_MAL))$ta7
options(scipen=10000, digits=5)
# plot_heatmap(confronto_M_MAL, 
#              taxa.order = rev(df$ta7),
#              sample.order = sample_names(confronto_M_MAL),
#              low="#0072B2", high="#D55E00", na.value = "white") + 
#   theme_bw()+
#   theme(axis.title.y=element_blank(),
#         axis.text.y=element_blank(),
#         axis.ticks.y=element_blank(),
#         axis.text.x = element_text(angle = 90, hjust = 1),
#         axis.title.x=element_blank())
##

library(pheatmap)
inp <- otu_table(confronto_M_MAL)
taxa_names(confronto_M_MAL)

inp <- inp[match(df$ta7, row.names(inp)),]

inp[inp == 0] <- NA

pheatmap(inp,
         color = colorRampPalette(c("#0072B2", "#D55E00"))(20),
         cluster_cols=F, 
         cluster_rows = F, 
         show_rownames = F,
         border_color = "grey27",
         gaps_col=c(19),
         gaps_row = c(7),
         na_col = "white",
         scale = "row")
## NB Pheatmap with "row" scaling is removing entry n° 3 and 4 which are present only in one sample



##################################
######## BUILDING FIGURE 7 #######
##################################

#
metadata <- read.csv("./human_16S.sampleinfo.csv", sep = ",")
input_metaAnalisi_science_data <- sample_data(metadata)
sample_names(input_metaAnalisi_science_data) <- as.character(sample_data(input_metaAnalisi_science_data)$Sample)

in_metadata <- sample_data(input_metaAnalisi_science_data)

# construct otu table
input_metaAnalisi_science_otu <- read.csv("./otu_table.csv", sep = ";", row.names = 1, header = T)
head(input_metaAnalisi_science_otu)

in_otu<- otu_table(input_metaAnalisi_science_otu, taxa_are_rows = T)

# construct tax table
# .biom table is not in this github repository, but can be found at github repository of original manuscript https://github.com/cribioinfo/sci2017_analysis
#input_metaAnalisi_science <- biom_otu_tax <- import_biom("./human_16S.even13190.abs.full.biom")
out <- tax_table(input_metaAnalisi_science)
#write.csv(out,"./taxonomy_table.csv")

input_metaAnalisi_science_tax <- read.csv("./taxonomy_table_refined.csv", sep = ";", row.names = 1, header = T)
head(input_metaAnalisi_science_tax)

in_tax<- tax_table(as.matrix(input_metaAnalisi_science_tax))

## JOIN tables

metaAnalisi_obj <- merge_phyloseq(in_otu, in_tax, in_metadata)


## UNIFORMING objects

sample_sums(qiime_file_16S_initial)
# rarefaction
depth <- 13190
qiime_file_16S_raref <- rarefy_even_depth(physeq = qiime_file_16S_initial, sample.size = depth, rngseed = T)
sample_sums(qiime_file_16S_raref)


sample_sums(qiime_file_16S_raref)
sample_sums(metaAnalisi_obj)

## summarize to genus

mio_obj_genus <- tax_glom(qiime_file_16S_raref, taxrank = "Genus")
tax_table(mio_obj_genus)

science_obj_genus <- tax_glom(metaAnalisi_obj, taxrank = "Genus")
tax_table(science_obj_genus)
science_obj_genus <- subset_taxa(science_obj_genus, Genus  != " g__")

## join

taxs <- as.data.frame(tax_table(mio_obj_genus))
taxa_names(mio_obj_genus) <- taxs$Genus

# remove thing like k__ and the "species" column
taxs <- as.matrix(tax_table(science_obj_genus)@.Data)
# remove column 7
taxs <- taxs[,1:6]
# where only "g__" or similar change in NA
taxs <- gsub(taxs, pattern = ".*__$", replacement = "NA")

taxs <- gsub(taxs, pattern = "^k__", replacement = "")
taxs <- gsub(taxs, pattern = "^ p__", replacement = "")
taxs <- gsub(taxs, pattern = "^ c__", replacement = "")
taxs <- gsub(taxs, pattern = "^ o__", replacement = "")
taxs <- gsub(taxs, pattern = "^ f__", replacement = "")
taxs <- gsub(taxs, pattern = "^ g__", replacement = "")

taxs <- as.data.frame(taxs)
taxs$Genus[duplicated(taxs$Genus)]

# remove duplicated genera

taxs$Genus[taxs$Genus == "Ruminococcus"]
otu_table(science_obj_genus)["366352",]
otu_table(science_obj_genus)["583398",]

taxs$Genus[taxs$Genus == "Prevotella"]
otu_table(science_obj_genus)["840914",]
otu_table(science_obj_genus)["196947",] 


taxs$Genus[taxs$Genus == "Clostridium"]
otu_table(science_obj_genus)["558420",]
otu_table(science_obj_genus)["New.CleanUp.ReferenceOTU10523",] 


####
# Define the taxa you don't want like this:
badTaxa = c("New.CleanUp.ReferenceOTU10523")
allTaxa = taxa_names(science_obj_genus)
allTaxa <- allTaxa[!(allTaxa %in% badTaxa)]
ex1 = prune_taxa(allTaxa, science_obj_genus)
# new phyloseq object with just the taxa you kept.
science_obj_genus <- ex1

#### 
otu_table(science_obj_genus)["840914",] <- otu_table(science_obj_genus)["840914",] + otu_table(science_obj_genus)["196947",] 

# Define the taxa you don't want like this:
badTaxa = c("196947")
allTaxa = taxa_names(science_obj_genus)
allTaxa <- allTaxa[!(allTaxa %in% badTaxa)]
ex1 = prune_taxa(allTaxa, science_obj_genus)
# new phyloseq object with just the taxa you kept.
science_obj_genus <- ex1

#### 
otu_table(science_obj_genus)["366352",] <- otu_table(science_obj_genus)["366352",] + otu_table(science_obj_genus)["583398",] 

# Define the taxa you don't want like this:
badTaxa = c("583398")
allTaxa = taxa_names(science_obj_genus)
allTaxa <- allTaxa[!(allTaxa %in% badTaxa)]
ex1 = prune_taxa(allTaxa, science_obj_genus)
# new phyloseq object with just the taxa you kept.
science_obj_genus <- ex1


## ok
science_obj_genus
taxs2 <- as.matrix(tax_table(science_obj_genus)@.Data)
# remove column 7
taxs2 <- taxs2[,1:6]
# where only "g__" or similar change in NA
taxs2 <- gsub(taxs2, pattern = ".*__$", replacement = "NA")

taxs2 <- gsub(taxs2, pattern = "^k__", replacement = "")
taxs2 <- gsub(taxs2, pattern = "^ p__", replacement = "")
taxs2 <- gsub(taxs2, pattern = "^ c__", replacement = "")
taxs2 <- gsub(taxs2, pattern = "^ o__", replacement = "")
taxs2 <- gsub(taxs2, pattern = "^ f__", replacement = "")
taxs2 <- gsub(taxs2, pattern = "^ g__", replacement = "")

taxs2 <- as.data.frame(taxs2)
tax_table(science_obj_genus) <- as.matrix(taxs2)
taxa_names(science_obj_genus) <- taxs2$Genus

####  join

tax_table(science_obj_genus)
tax_table(mio_obj_genus)

metaAnalisi_obj <- merge_phyloseq(science_obj_genus, mio_obj_genus)
sample_sums(metaAnalisi_obj)
tax_table(metaAnalisi_obj)

plot_heatmap(metaAnalisi_obj, na.value = "white",low="#0072B2", high="#D55E00")

metadata_metaAnalisi <- as.data.frame(sample_data(metaAnalisi_obj))
c(as.character(metadata_metaAnalisi$Response[1:42]), as.character(metadata_metaAnalisi$Pathology[43:80]))
sample_data(metaAnalisi_obj)$MetaAnalisi <- c(as.character(metadata_metaAnalisi$Response[1:42]), as.character(metadata_metaAnalisi$Pathology[43:80]))

sample_data(metaAnalisi_obj)$LabelMeta <- sample_names(metaAnalisi_obj)

### ORDINATION (Figure 6)

metaAnalisi_obj_noV_noM6 <- subset_samples(metaAnalisi_obj, LabelMeta != "M6F.V3.V4")

relabel <- c(as.character(sample_data(metaAnalisi_obj_noV_noM6)$Sample[1:42]),
             as.character(sample_data(metaAnalisi_obj_noV_noM6)$Label_bis[43:length(sample_data(metaAnalisi_obj_noV_noM6)$Label_bis)]))

sample_names(metaAnalisi_obj_noV_noM6) <- relabel

sample_data(metaAnalisi_obj_noV_noM6)$Study <- c(rep("Matson_et_al.", 42),rep("This_study", 37))

xtx <- transform(metaAnalisi_obj_noV_noM6,"compositional", target = "OTU")

#library(metagMisc)
#meaAnalisi_filtered2 <- phyloseq_filter_prevalence(xtx, prev.trh = 1, abund.trh = 0.000001) 
# function broke if no abund.trh is given, last version do not filter

table(sample_data(xtx)$Study) # 42 samples from matson and 37 from this study
# can use 37 as treshold

# extract otu table
prevotu <- otu_table(xtx)
colSums(prevotu)
rowSums(prevotu)
# convert to 1/0
prevotu_01 <- as.data.frame(prevotu)
prevotu_01[prevotu_01!=0] <- 1
# select only otu in more than 37 samples
prevotu_01 <- prevotu_01[rowSums(prevotu_01) > 37,]
otu_to_keep <- row.names(prevotu_01)

# select those otus from initial object

meaAnalisi_filtered2 <- prune_taxa(taxa = otu_to_keep, xtx)

# remove Oscillospira
# Define the taxa you don't want like this:
badTaxa = c("Oscillospira")
allTaxa = taxa_names(meaAnalisi_filtered2)
allTaxa <- allTaxa[!(allTaxa %in% badTaxa)]
ex1 = prune_taxa(allTaxa, meaAnalisi_filtered2)
# new phyloseq object with just the taxa you kept.
meaAnalisi_filtered2 <- ex1

library(pheatmap)
inp <- otu_table(meaAnalisi_filtered2)
inp[inp == 0] <- NA

pheatmap(inp,
         clustering_distance_rows = "euclidean",
         fontsize_col = 8,
         fontsize_row = 8,
         color = colorRampPalette(c("white","#0072B2", "#D55E00"))(50),
         cluster_cols=F, 
         cluster_rows = F, 
         show_rownames = T,
         border_color = "grey27",
         gaps_col = c(42,57),
         na_col = "white",
         scale = "row")

colori <- factor(sample_data(meaAnalisi_filtered2)$MetaAnalisi, levels = c("Leucoderma_Melanoma", "Melanoma", "NonResponder", "Responder", "Control"))
sample_data(meaAnalisi_filtered2)$Colors <- colori

ordination_metanalisi <- ordinate(meaAnalisi_filtered2, "PCoA", distance = "bray")
p_ordination_metanalisi <- plot_ordination(physeq = meaAnalisi_filtered2, ordination_metanalisi, axes = c(1,2), color = "Colors", shape = "Study")
p12 <- p_ordination_metanalisi + geom_point(size = 4) 
p12 <- p12 + theme_bw() + scale_color_jama() + 
  theme( panel.grid.major =  element_blank(),
         panel.grid.minor =  element_blank())
p12 <- p12 + xlab("PCoA 1 (12.1%)") +
  ylab("PCoA 2 (9%)") + scale_shape_manual(values = c(15,16))
  #stat_ellipse(type = "t", linetype = 2) 



set.seed(12387)
df <- as(sample_data(meaAnalisi_filtered2), "data.frame")
d.bacteria = phyloseq::distance(meaAnalisi_filtered2, "bray")
adonis(d.bacteria ~ Colors, data = df, permutations = 9999)
library(pairwiseAdonis)
pairwise.adonis(d.bacteria, factors = df$Colors, perm = 9999, p.adjust.m = "bonferroni")


rowSums(otu_table(meaAnalisi_filtered2))
colSums(otu_table(meaAnalisi_filtered2))

colnames(otu_table(meaAnalisi_filtered2)[,1:42])
colnames(otu_table(meaAnalisi_filtered2)[,58:76])

media_matsen <- rowMeans2(otu_table(meaAnalisi_filtered2)[,1:42])
media_nostro <- rowMeans2(otu_table(meaAnalisi_filtered2)[,58:76])
sd_matsen <- rowSds(otu_table(meaAnalisi_filtered2)[,1:42])
sd_nostro <- rowSds(otu_table(meaAnalisi_filtered2)[,58:76])
df <- cbind(media_matsen,media_nostro, sd_matsen, sd_nostro)
row.names(df) <- row.names(otu_table(meaAnalisi_filtered2))
df <- as.data.frame(df)
df <- df[order(-df$media_nostro),]

df$diff <- abs(df$media_nostro-df$media_matsen)
df$label <- row.names(df)

df$abbundHiger <- ifelse(df$media_matsen > df$media_nostro, "Matsen","Vitali")

library(ggrepel)

top_n(df,n = 10, wt = diff)

df$label[which(df$diff < 0.01)] <- "" # remove label of otus with less than 5% variation in mean abb

ggplot(data = df, aes(x = media_nostro, y = media_matsen, label = label)) +  
  geom_pointrange(aes(ymin = media_matsen - sd_matsen, ymax = media_matsen + sd_matsen)) +
  geom_errorbarh(aes(xmax = media_nostro + sd_nostro, xmin = media_nostro - sd_nostro, height = 0)) +
  geom_point(aes(size = diff)) +
  geom_segment(aes(x = 0, y = 0, xend = 1, yend = 1), colour = "black", linetype  ="dotted") + 
  xlab("Mean relative abundance of this work") + ylab("Mean relative abundance of Matsen at al.") + 
  theme_cowplot() + 
  geom_label_repel(aes(color = abbundHiger),
  arrow = arrow(length = unit(0.01, "npc"), type = "closed", ends = "first"),
  force = 10)





