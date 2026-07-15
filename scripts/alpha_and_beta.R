# Daniel Castaneda Mogollon, PhD
# April 19th, 2026
# This script evaluates alpha and beta diversity of the samples in the Parkinson's project. 

library(vegan)
library(ggplot2)
library(phyloseq)
library(devtools)

install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")

library(pairwiseAdonis)


metadata_df = read.csv("/Users/danielcm/Desktop/MIGUT/metadata/MIGUT_metadata_in_tab_format.txt", sep = "\t")
rownames(metadata_df) = metadata_df$ID
species_table_abs = read.csv("/Users/danielcm/Desktop/MIGUT/analysis_may2026/filtered_metaphlan4_abundances/merged_absolute_abundances_table_GTDB_filtered.txt", sep = "\t", row.names = 1)

#species_table_rel_f = read.csv("/Users/danielcm/Desktop/MIGUT/novaseq/metaphlan_merged_gtdb_kevin/merged_abundances_kevin_script/species_filtered_gtdb_prev_10_mean_00001.txt", sep = "\t", row.names = 1)

# TSS normalization: divide each sample (column) by its total sum
species_table_abs = species_table_abs[, colSums(species_table_abs) > 0]  # Remove zero-sum samples before normalizing
species_table_rel_f = sweep(species_table_abs, 2, colSums(species_table_abs), FUN = "/")

species_table_rel_f = species_table_rel_f[rowSums(species_table_rel_f) > 0, ]

ps_object_rel = phyloseq(otu_table(as.matrix(species_table_rel_f), taxa_are_rows = TRUE), sample_data(metadata_df))
ps_object_abs = phyloseq(otu_table(as.matrix(species_table_abs), taxa_are_rows = TRUE), sample_data(metadata_df))
sample_data(ps_object_rel)$Group_and_sex

#Group and sex
meta_rel = data.frame(sample_data(ps_object_rel))
distance = vegdist(t(otu_table(ps_object_rel)), method = "bray")
permanova = adonis2(distance ~ Group_and_sex, data = meta_rel)
permanova


View(meta_rel)

#Group


#PD only
ps_pd <- subset_samples(ps_object_rel, Group == "PD")
ps_pd <- prune_samples(sample_sums(ps_pd) > 0, ps_pd)

dist_pd <- phyloseq::distance(ps_pd, method = "bray")
meta_pd <- data.frame(sample_data(ps_pd))

adonis2(dist_pd ~ Sex_Male1_Female2, data = meta_pd)

#HC only
ps_hc <- subset_samples(ps_object_rel, Group == "HC")
ps_hc <- prune_samples(sample_sums(ps_hc) > 0, ps_hc)

dist_hc <- phyloseq::distance(ps_hc, method = "bray")
meta_hc <- data.frame(sample_data(ps_hc))
    
adonis2(dist_hc ~ Sex_Male1_Female2, data = meta_hc)

# Males only
ps_males = subset_samples(ps_object_rel, Sex_Male1_Female2 == 1)
ps_males = prune_samples(sample_sums(ps_males) > 0, ps_males)

dist_males = phyloseq::distance(ps_males, method = "bray")
meta_male = data.frame(sample_data(ps_males))

adonis2(dist_males ~ Group, data = meta_male)

# Females only
ps_females = subset_samples(ps_object_rel, Sex_Male1_Female2 == 2)
ps_females = prune_samples(sample_sums(ps_females) > 0, ps_females)

dist_females = phyloseq::distance(ps_females, method = "bray")
meta_female = data.frame(sample_data(ps_females))

adonis2(dist_females ~ Group, data = meta_female)

# Ordinating
ordination = ordinate(ps_object_rel, method = "PCoA", distance = "bray")


png("/Users/danielcm/Desktop/MIGUT/analysis_may2026/alpha_beta/pcoa_bray_prev10_by_diagnosis.png", width = 10, height = 10, units = "in", res = 600)

# Plotting all individuals by diagnosis
p1 = plot_ordination(ps_object_rel, ordination = ordination, color = "Group") +
    geom_point(size = 4) +
    stat_ellipse(aes(group = Group), type = "t", linetype = 1, linewidth = 1.5) +
    scale_color_manual(values = c("HC" = "#030202", "PD" = "#1F78B4")) +
    theme(
    axis.title.x = element_text(size = 22),
    axis.title.y = element_text(size = 22),
    axis.text.x = element_text(size = 22),
    axis.text.y = element_text(size = 22),
    legend.title = element_text(size = 22),
    legend.text = element_text(size = 22),
    plot.title = element_text(size = 22),
    aspect.ratio = 1  # Make PCA plot square
  )
p1

permanova = adonis2(distance ~ Group, data = meta_rel)
permanova

dev.off()



# Plotting all individuals divided by sex and diagnosis
x = as.data.frame(sample_data(ps_pd))
x$Sex_Male1_Female2


png("/Users/danielcm/Desktop/MIGUT/analysis_may2026/alpha_beta/pcoa_bray_prev10_by_parkinsons_only.png", width = 10, height = 10, units = "in", res = 600)

p2 = plot_ordination(ps_pd, ordination = ordinate(ps_pd,"PCoA","bray"), color = "Sex_Male1_Female2") +
    geom_point(size = 4) +
    stat_ellipse(aes(group = Sex_Male1_Female2), type = "t", linetype = 1, linewidth = 1.5) +
    scale_color_manual(values = c(1 = "red", 2 = "pink")) +
    theme(
    axis.title.x = element_text(size = 22),
    axis.title.y = element_text(size = 22),
    axis.text.x = element_text(size = 22),
    axis.text.y = element_text(size = 22),
    legend.title = element_text(size = 22),
    legend.text = element_text(size = 22),
    plot.title = element_text(size = 22),
    aspect.ratio = 1  # Make PCA plot square
  )
p2

png("/Users/danielcm/Desktop/MIGUT/analysis_may2026/alpha_beta/pcoa_bray_prev10_sex_and_diagnosis.png", width = 10, height = 10, units = "in", res = 600)


p3 = plot_ordination(ps_object_rel, ordination = ordination, color = "Group_and_sex") +
    geom_point(size = 4) +
    stat_ellipse(aes(group = Group_and_sex), type = "t", linetype = 1, linewidth = 1.5) +
    scale_color_manual(values = c("HC_Male" = "#030202", "HC_Female" = "#1F78B4", "PD_Male" = "#33A02C", "PD_Female" = "#E31A1C")) +
    theme(
    axis.title.x = element_text(size = 22),
    axis.title.y = element_text(size = 22),
    axis.text.x = element_text(size = 22),
    axis.text.y = element_text(size = 22),
    legend.title = element_text(size = 22),
    legend.text = element_text(size = 22),
    plot.title = element_text(size = 22),
    aspect.ratio = 1  # Make PCA plot square
  )
p3

dev.off()
permanova = adonis2(distance ~ Group_and_sex, data = meta_rel)
pairwise_results = pairwise.adonis2(x = distance ~ Group_and_sex, data = meta_rel, p.adjust.method = "BH")
pairwise_results


png("/Users/danielcm/Desktop/MIGUT/analysis_may2026/alpha_beta/pcoa_male_prev10_bray.png", width = 10, height = 10, units = "in", res = 600)


p4 = plot_ordination(ps_males, ordination = ordinate(ps_males,"PCoA","bray"), color = "Group_and_sex") +
    geom_point(size = 4) +
    stat_ellipse(aes(group = Group_and_sex), type = "t", linetype = 1, linewidth = 1.5) +
    scale_color_manual(values = c("HC_Male" = "purple", "PD_Male" = "brown")) +
    theme(
    axis.title.x = element_text(size = 22),
    axis.title.y = element_text(size = 22),
    axis.text.x = element_text(size = 22),
    axis.text.y = element_text(size = 22),
    legend.title = element_text(size = 22),
    legend.text = element_text(size = 22),
    plot.title = element_text(size = 22),
    aspect.ratio = 1  # Make PCA plot square
  )
p4

dev.off()

png("/Users/danielcm/Desktop/MIGUT/analysis_may2026/alpha_beta/pcoa_female_bray.png", width = 10, height = 10, units = "in", res = 600)


p5 = plot_ordination(ps_females, ordination = ordinate(ps_females,"PCoA","bray"), color = "Group_and_sex") +
    geom_point(size = 4) +
    stat_ellipse(aes(group = Group_and_sex), type = "t", linetype = 1, linewidth = 1.5) +
    scale_color_manual(values = c("HC_Female" = "pink", "PD_Female" = "darkgreen")) +
    theme(
    axis.title.x = element_text(size = 22),
    axis.title.y = element_text(size = 22),
    axis.text.x = element_text(size = 22),
    axis.text.y = element_text(size = 22),
    legend.title = element_text(size = 22),
    legend.text = element_text(size = 22),
    plot.title = element_text(size = 22),
    aspect.ratio = 1  # Make PCA plot square
  )
p5


dev.off()

sample_data(ps_pd)$H_and_Y_stage = factor(sample_data(ps_pd)$H_and_Y_stage)
p6 = plot_ordination(ps_pd, ordination = ordinate(ps_pd,"PCoA","bray"), color = "H_and_Y_stage") +
    geom_point(size = 4) +
    stat_ellipse(aes(group = factor(H_and_Y_stage)), type = "t", linetype = 1, linewidth = 1.5) +
    scale_color_manual(values = c("1" = "#030202", "2" = "#1F78B4", "3" = "#33A02C", "4" = "#E31A1C"),
                       labels = c("1","2","3","4"),
                       breaks = c("1","2","3","4")) +
    theme(
    axis.title.x = element_text(size = 22),
    axis.title.y = element_text(size = 22),
    axis.text.x = element_text(size = 22),
    axis.text.y = element_text(size = 22),
    legend.title = element_text(size = 22),
    legend.text = element_text(size = 22),
    plot.title = element_text(size = 22),
    aspect.ratio = 1  # Make PCA plot square
  )
p6

sample_data(ps_object_abs)

shannon_diversity = estimate_richness(ps_object_abs, measures = "Shannon")
simpson_diversity = estimate_richness(ps_object_abs, measures = "Simpson")
observed_diversity = estimate_richness(ps_object_abs, measures = "Observed")
meta <- data.frame(sample_data(ps_object_abs))
shannon_diversity$SampleID <- rownames(shannon_diversity)
shannon_diversity$Group_and_Sex = meta$Group_and_sex
simpson_diversity$SampleID <- rownames(simpson_diversity)
simpson_diversity$Group_and_Sex = meta$Group_and_sex
observed_diversity$SampleID <- rownames(observed_diversity)
observed_diversity$Group_and_Sex = meta$Group_and_sex


df_alpha <- merge(shannon_diversity,simpson_diversity, meta, by.x = "SampleID", by.y = "row.names")
df_alpha <- merge(df_alpha, observed_diversity, meta, by.x = "SampleID.y", by.y = "row.names")


colnames(df_alpha)
write.table(df_alpha, "/Users/danielcm/Desktop/MIGUT/analysis_may2026/alpha_diversity.txt", sep = "\t", row.names = FALSE)




