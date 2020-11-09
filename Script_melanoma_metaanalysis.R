#################################################################################################################
### Script to perform data analysis and obtain figure of Vitali et. al 2020 manuscript                  #########
### Author: Francesco Vitali                                                                            #########
### Licence: This work is licensed under a Creative Commons Attribution 4.0 International License       #########
#################################################################################################################


#Load Package
library(tidyverse)
library(ggsci)
library(vegan)
library(NbClust)
library(FactoMineR)
library(metagenomeSeq)
library(DESeq2)
library(ggpubr)
library(microbiome)
library(phyloseq)
library(data.table)

#################
## FILE IMPORT ##
#################

otus_16S_fp <- "./denovo_greedy_metaanalysis/otutable.csv"
meta_16S_fp <- "./denovo_greedy_metaanalysis/metadata_metaanalysis.csv"
tax_16S_fp <- "./denovo_greedy_metaanalysis/taxa.csv"
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

# rename phylogenetic rank
colnames(tax_table(qiime_file_16S_initial)) <- c(Rank1 = "Kingdom", Rank2 = "Phylum", Rank3 = "Class", Rank4 = "Order", Rank5 = "Family", Rank6 = "Genus")

# remove mito and chloro
taxa <- as.data.frame(tax_table(qiime_file_16S_initial)@.Data)

na.omit(taxa[taxa$Family == "Mitochondria",])
na.omit(taxa[taxa$Family == "Chloroplast",])

#qiime_file_16S_initial <- subset_taxa(qiime_file_16S_initial, Family  != "Mitochondria") # un ci sono... 
#qiime_file_16S_initial <- subset_taxa(qiime_file_16S_initial, Family  != "Chloroplast")
#qiime_file_16S_initial

summary(sample_sums(qiime_file_16S_initial)) # variation interval 31736 - 172771
sum(sample_sums(qiime_file_16S_initial))

## REMOVE ML
qiime_file_16S_initial <- subset_samples(qiime_file_16S_initial, Pathology != "ML")

# Plot reads to check for bad samples
sdt = data.table(as(sample_data(qiime_file_16S_initial), "data.frame"),
                 TotalReads = sample_sums(qiime_file_16S_initial), keep.rownames = TRUE)
setnames(sdt, "rn", "SampleID")
means <- aggregate(TotalReads ~  Pathology*Study, sdt, mean)
sum(sdt$TotalReads) # a total of 4034643 reads were obtained

min(sdt$TotalReads)

p_reads <- ggboxplot(data = sdt, x = "Pathology", y = "TotalReads", add = "segments", 
                     color = "Pathology", palette ="npg", rotate = T, dot.size = 4) + theme_cleveland() +
  stat_compare_means(method = "anova", label.y = 3000, label.x = 3.5)
p_reads

# rarefaction scaling

depth <- min(sdt$TotalReads)
  qiime_file_16S_raref <- rarefy_even_depth(physeq = qiime_file_16S_initial, sample.size = depth, rngseed = T)
sample_sums(qiime_file_16S_raref)

qiime_file_16S_raref_a <- qiime_file_16S_raref

qiime_file_16S_raref <- transform(qiime_file_16S_raref, transform = "compositional")


### Ordination analysis
sample_names(qiime_file_16S_raref)

mypal = pal_aaas("default", alpha = 1)(9)
mypal
library("scales")
show_col(mypal)

ordination_metanalisi <- ordinate(qiime_file_16S_raref, "PCoA", distance = "bray")
p_ordination_metanalisi <- plot_ordination(physeq = qiime_file_16S_raref, ordination_metanalisi, 
                                           axes = c(1,2), shape = "Study") + scale_shape_manual(values = c(21,24))
p12 <- p_ordination_metanalisi + stat_ellipse(linetype = 2) + 
  geom_point(size = 4, aes(fill = Pathology)) 
p12 <- p12 + theme_bw() + scale_fill_manual(values = c("#3B4992FF","#EE0000FF","#A20056FF")) + 
  theme( panel.grid.major =  element_blank(),
         panel.grid.minor =  element_blank())
p12 <- p12 + xlab("PCoA 1 (15.6%)") +
  ylab("PCoA 2 (7.6%)")
p12




ordination_metanalisi <- ordinate(qiime_file_16S_raref, "PCoA", distance = "bray")
p_ordination_metanalisi <- plot_ordination(physeq = qiime_file_16S_raref, ordination_metanalisi, 
                                           axes = c(1,2), shape = "Study") #+ scale_shape_manual(values = c(21,22))
p12 <- p_ordination_metanalisi + stat_ellipse(linetype = 2) + 
  geom_point(size = 4, aes(color = Plot)) 
p12 <- p12 + theme_bw() + scale_color_jama() + 
  theme( panel.grid.major =  element_blank(),
         panel.grid.minor =  element_blank())
p12 <- p12 + xlab("PCoA 1 (15.6%)") +
  ylab("PCoA 2 (7.6%)")
p12


### test 

set.seed(12387)
df <- as(sample_data(qiime_file_16S_raref), "data.frame")
d.bacteria = phyloseq::distance(qiime_file_16S_raref, "bray")
adonis(d.bacteria ~ Pathology, data = df, permutations = 9999)
library(pairwiseAdonis)
pairwise.adonis(d.bacteria, factors = df$Pathology, perm = 9999, p.adjust.m = "bonferroni")


## no control

ordination_metanalisi <- ordinate(subset_samples(qiime_file_16S_raref, Pathology != "C"), "PCoA", distance = "bray")
p_ordination_metanalisi <- plot_ordination(physeq = qiime_file_16S_raref, ordination_metanalisi, axes = c(1,2), shape = "Study")
p12 <- p_ordination_metanalisi + stat_ellipse(linetype = 2) + 
  geom_point(size = 4, aes(color = Varible)) 
p12 <- p12 + theme_bw() + scale_color_jama() + 
  theme( panel.grid.major =  element_blank(),
         panel.grid.minor =  element_blank())
p12 <- p12 + xlab("PCoA 1 (15.1%)") +
  ylab("PCoA 2 (8.3%)") + scale_shape_manual(values = c(15,16,17,18)) 
p12

library(ggrepel)

sample_data(qiime_file_16S_raref)$ID <- substr(sample_names(qiime_file_16S_raref),1,4)

ordination_metanalisi <- ordinate(subset_samples(qiime_file_16S_raref, Pathology != "C"), "PCoA", distance = "bray")
p_ordination_metanalisi <- plot_ordination(physeq = qiime_file_16S_raref, ordination_metanalisi, axes = c(1,2), shape = "Study")
p12 <- p_ordination_metanalisi + stat_ellipse(linetype = 2) + 
  geom_point(size = 4, aes(color = Varible)) 
p12 <- p12 + theme_bw() + scale_color_jama() + 
  theme( panel.grid.major =  element_blank(),
         panel.grid.minor =  element_blank())
p12 <- p12 + xlab("PCoA 1 (15.1%)") +
  ylab("PCoA 2 (8.3%)") + scale_shape_manual(values = c(15,16,17,18)) 
p12 + geom_text_repel(aes(label = ID))

##

set.seed(12387)
df <- as(sample_data(subset_samples(qiime_file_16S_raref, Pathology != "C")), "data.frame")
d.bacteria = phyloseq::distance(subset_samples(qiime_file_16S_raref, Pathology != "C"), "bray")
adonis(d.bacteria ~ Varible, data = df, permutations = 9999)
library(pairwiseAdonis)
pairwise.adonis(d.bacteria, factors = df$Varible, perm = 9999, p.adjust.m = "bonferroni")

### Alpha diversity analysis

alphadiv <- microbiome::alpha(qiime_file_16S_raref_a) # rarefied but not relative dataset

df <- as(sample_data(qiime_file_16S_raref_a), "data.frame")

df$rich <- alphadiv$observed
df$shan <- alphadiv$diversity_shannon
df$eve <- alphadiv$evenness_pielou
df$sim <- alphadiv$diversity_inverse_simpson
df$dom <- alphadiv$dominance_gini

#df_melanomi <- df[df$Pathology!="C",]

##' Plot and test on Invasivity

variabile_plot = "Pathology"

# get levels of comparison
lev <- levels(as.factor(df$Pathology))
my_comparisons <- list(c(lev[1], lev[2]),c(lev[1], lev[3]),
                       c(lev[2], lev[3]))

p_rich <- ggboxplot(data = df, x = variabile_plot, y = "rich", add = "jitter", size = 0.5, fill = variabile_plot, outlier.shape = NA)
p_rich <- p_rich + stat_compare_means(method = "t.test", comparisons = my_comparisons,  hide.ns = F, label = "p.signif")                   # Add global p-value
p_rich <- p_rich + ylab("Richness (n° of OTUs)") + scale_fill_manual(values = c("#3B4992FF","#EE0000FF","#A20056FF")) + rremove("legend") + rremove("xlab") + rremove("x.text") + rremove("x.ticks") #+ theme(axis.text.x = element_text(angle = 45, hjust = 1))
p_shan <- ggboxplot(data = df, x = variabile_plot, y = "shan", add = "jitter", size = 0.5, fill = variabile_plot, outlier.shape = NA)
p_shan <- p_shan +stat_compare_means(method = "t.test", comparisons = my_comparisons,  hide.ns = F, label = "p.signif")                  # Add global p-value
p_shan <- p_shan + ylab("Shannon's Index") +  scale_fill_manual(values = c("#3B4992FF","#EE0000FF","#A20056FF"))+ rremove("legend") + rremove("xlab") + rremove("x.text") + rremove("x.ticks") #+ theme(axis.text.x = element_text(angle = 45, hjust = 1))
p_eve <- ggboxplot(data = df, x = variabile_plot, y = "eve", add = "jitter", size = 0.5, fill = variabile_plot, outlier.shape = NA)
p_eve <- p_eve + stat_compare_means(method = "t.test", comparisons = my_comparisons,  hide.ns = F, label = "p.signif")                   # Add global p-value
p_eve <- p_eve + ylab("Evenness Index") +  scale_fill_manual(values = c("#3B4992FF","#EE0000FF","#A20056FF"))+ rremove("legend") + rremove("xlab") + rremove("x.text") + rremove("x.ticks") #+ theme(axis.text.x = element_text(angle = 45, hjust = 1))
p_sim <- ggboxplot(data = df, x = variabile_plot, y = "sim", add = "jitter", size = 0.5, fill = variabile_plot, outlier.shape = NA)
p_sim <- p_sim + stat_compare_means(method = "t.test", comparisons = my_comparisons, hide.ns = F, label = "p.signif")                 # Add global p-value
p_sim <- p_sim + ylab("Simpson Index") +  scale_fill_manual(values = c("#3B4992FF","#EE0000FF","#A20056FF"))+ rremove("legend") + rremove("xlab") + rremove("x.text") + rremove("x.ticks") #+ theme(axis.text.x = element_text(angle = 45, hjust = 1))

library(cowplot)
plot_alpha<- plot_grid(p_rich, p_eve, ncol = 1) + theme_void()

plot_alpha

##' Plot and test on Invasivity

variabile_plot = "Study"

# get levels of comparison

p_rich <- ggboxplot(data = df_melanomi, x = variabile_plot, y = "rich", add = "jitter", size = 0.5, fill = variabile_plot)
p_rich <- p_rich + stat_compare_means(method = "t.test",  hide.ns = T)                   # Add global p-value
p_rich <- p_rich + ylab("Richness (n° of OTUs)") + scale_fill_grey() + rremove("legend") + rremove("xlab") + theme(axis.text.x = element_text(angle = 45, hjust = 1))
p_shan <- ggboxplot(data = df_melanomi, x = variabile_plot, y = "shan", add = "jitter", size = 0.5, fill = variabile_plot)
p_shan <- p_shan +stat_compare_means(method = "t.test",  hide.ns = T)                  # Add global p-value
p_shan <- p_shan + ylab("Shannon's Index") + scale_fill_grey()+ rremove("legend") + rremove("xlab") + theme(axis.text.x = element_text(angle = 45, hjust = 1))
p_eve <- ggboxplot(data = df_melanomi, x = variabile_plot, y = "eve", add = "jitter", size = 0.5, fill = variabile_plot)
p_eve <- p_eve + stat_compare_means(method = "t.test", hide.ns = T)                   # Add global p-value
p_eve <- p_eve + ylab("Evenness Index") + scale_fill_grey()+ rremove("legend") + rremove("xlab") + theme(axis.text.x = element_text(angle = 45, hjust = 1))
p_sim <- ggboxplot(data = df_melanomi, x = variabile_plot, y = "sim", add = "jitter", size = 0.5, fill = variabile_plot)
p_sim <- p_sim + stat_compare_means(method = "t.test", hide.ns = T)                 # Add global p-value
p_sim <- p_sim + ylab("Simpson Index") + scale_fill_grey()+ rremove("legend") + rremove("xlab") + theme(axis.text.x = element_text(angle = 45, hjust = 1))

library(cowplot)
plot_alpha<- plot_grid(p_rich, p_eve,p_shan, p_sim, ncol = 4, labels = c("A","B","C","D"))

plot_alpha

# Heatmap

plot_heatmap(qiime_file_16S_raref)

otus_table <- as.data.frame(otu_table(tax_glom(qiime_file_16S_raref, taxrank = "Genus")))

pheatmap::pheatmap(as.matrix(otus_table), scale = "column",
                   clustering_distance_rows = "euclidean", 
                   clustering_distance_cols = "euclidean",
                   clustering_method = "ward.D2")
  
#### SIMPER analysis 

# aim to see main contributors in the difference between pre and post metastatic melanoma 

set.seed(12387)
df <- as(sample_data(subset_samples(qiime_file_16S_raref, Pathology != "C")), "data.frame")
d.bacteria = phyloseq::distance(subset_samples(qiime_file_16S_raref, Pathology != "C"), "bray")
adonis(d.bacteria ~ Varible, data = df, permutations = 9999)

set.seed(12387)
df <- as(sample_data(subset_samples(qiime_file_16S_raref, Pathology != "C")), "data.frame")
d.bacteria = phyloseq::distance(subset_samples(qiime_file_16S_raref, Pathology != "C"), "bray")
adonis(d.bacteria ~ Study, data = df, permutations = 9999)

# in the last adonis, we see that the "study" variable (i.e. early vs metastatic melanoma) is significant

ordination_metanalisi <- ordinate(subset_samples(qiime_file_16S_raref, Pathology != "C"), "PCoA", distance = "bray")
p_ordination_metanalisi <- plot_ordination(physeq = qiime_file_16S_raref, ordination_metanalisi, axes = c(1,2), shape = "Study")
p12 <- p_ordination_metanalisi + stat_ellipse(linetype = 2) + 
  geom_point(size = 4, aes(color = Study)) 
p12 <- p12 + theme_bw() + scale_color_jama() + 
  theme( panel.grid.major =  element_blank(),
         panel.grid.minor =  element_blank())
p12 <- p12 + xlab("PCoA 1 (15.1%)") +
  ylab("PCoA 2 (8.3%)") + scale_shape_manual(values = c(15,16,17,18)) 
p12

# and here we visualize it. As ordination and adonis are based on BC distance, we know that BC distance is able to differentiate
# between the two studies

simper_data <- tax_glom(subset_samples(qiime_file_16S_raref_a, Pathology != "C"), taxrank = "Genus")
any(taxa_sums(simper_data) == 0)
simper_data <- prune_taxa(taxa_sums(simper_data) > 1, simper_data)

simper_table <- as.data.frame(otu_table(simper_data))

simper_metadata <- as(sample_data(simper_data), "data.frame")

sim <- simper(comm = t(simper_table), group = simper_metadata$Study, permutations = 100)

sim

summary(sim)

## take vector of taxa from simper

simper_otu_res <- c("DENOVO4", "DENOVO1", "DENOVO5", "DENOVO2", "DENOVO14", "DENOVO17", "DENOVO26", "DENOVO42", "DENOVO25", "DENOVO7", "DENOVO9", "DENOVO13", "DENOVO19", "DENOVO36")

simper_phyloseq <- prune_taxa(simper_otu_res, simper_data)


ordination_metanalisi <- ordinate(simper_phyloseq, "PCoA", distance = "bray")
p_ordination_metanalisi <- plot_ordination(physeq = simper_phyloseq, ordination_metanalisi, axes = c(1,2), shape = "Study")
p12 <- p_ordination_metanalisi + stat_ellipse(linetype = 2) + 
  geom_point(size = 4, aes(color = Study)) 
p12 <- p12 + theme_bw() + scale_color_jama() + 
  theme( panel.grid.major =  element_blank(),
         panel.grid.minor =  element_blank())
p12 <- p12 + xlab("PCoA 1 (15.1%)") +
  ylab("PCoA 2 (8.3%)") + scale_shape_manual(values = c(15,16,17,18)) 
p12

ss <- transform_sample_counts(merge_samples(simper_phyloseq, group = "Invasivity"), function(x) 100 * x/sum(x))
plot_bar(ss, fill = "Genus") + scale_fill_igv()




  