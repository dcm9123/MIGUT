# Daniel Castaneda Mogollon, PhD
# July 15th, 2026
# Purpose: This script will run ANCOM-BC2 on microbiome data sets

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("ANCOMBC")

library(ANCOMBC)
library(phyloseq)
library(dplyr)


# Reading the main files to work with, that is, ASV count tables and metadata files.    
dataset_calgary = as.data.frame(read.csv("/Users/danielcm/Desktop/MIGUT/analysis_may2026/results1/merged_tables_abundances/merged_absolute_abundances_table_GTDB.txt", sep = "\t", row.names = 1))
metadata_calgary = as.data.frame(read.csv("/Users/danielcm/Desktop/MIGUT/metadata/MIGUT_metadata_in_tab_format.txt", sep = "\t"))
rownames(metadata_calgary) = metadata_calgary$ID
dataset_vancouver = as.data.frame(read.csv("/Users/danielcm/Desktop/MIGUT/vancouver/merged_tables/gtdb_absolute_abundance_merged.txt", sep = "\t", row.names = 1))
metadata_vancouver = as.data.frame(read.csv("/Users/danielcm/Desktop/MIGUT/vancouver/vancouver_participant_metadata.csv", sep = ",", row.names = 1))
colnames(dataset_vancouver) = sapply(strsplit(colnames(dataset_vancouver), "_"), `[`, 1)  # Removes the "_S123_gtdb_profile" part so it matches the metadata's sample names


metadata_vancouver$Location = "Vancouver"
metadata_calgary$Location = "Calgary"

# Originally, there were 184 samples in the Calgary cohort metadata's file, but 10 of these were not sequenced (17, 18, 19, 20, 21, 22, 23, 24, 168 and 187)
# Two sample IDs that were sequenced do not have a match in the metadata file from Davide: 40, 70
# 18 samples sequenced had too few reads and were removed from the analysis: 199, 212, 50, 107, 216, 28, 81, 191, 218, 211, 36, 83, 192, 217, 26, 182, 112, 161

# From 184 to 

print(dim(dataset_calgary)) # 176 samples in total
print(dim(metadata_calgary)) # 184 samples in total. 


# Removing samples that were originally excluded for any reason. This assumes that the metadata has them all but not the dataset count table, so the latter is a subset of the metadata.   
excluding_samples = function(dataset_calgary, metadata_calgary, remove_flagged_samples, flagged_samples){
    samples_to_exclude = c()
    for(item in metadata_calgary$ID){
        if(!item %in% colnames(dataset_calgary)){
            samples_to_exclude = c(samples_to_exclude, item)
        }
    }
    
    metadata_calgary_f = metadata_calgary[!metadata_calgary$ID %in% samples_to_exclude, ]
    samples_to_exclude = c()

    for(col in colnames(dataset_calgary)){
        if(!col %in% metadata_calgary_f$ID){
            samples_to_exclude = c(samples_to_exclude, col)
        }
    }
    dataset_calgary_f = dataset_calgary[, !colnames(dataset_calgary) %in% samples_to_exclude]

    if(remove_flagged_samples == TRUE){
        metadata_calgary_f = metadata_calgary_f[!metadata_calgary_f$ID %in% flagged_samples, ]
        dataset_calgary_f = dataset_calgary_f[, !colnames(dataset_calgary_f) %in% flagged_samples]
    }

    return(list(metadata_calgary_f = metadata_calgary_f, dataset_calgary_f = dataset_calgary_f))
}

# Filters by taxonomic rank, most likely will be always species
filtering_dataset_by_taxa_rank = function(dataset, taxa_rank){
    taxa_to_remove = c()
    for(row in rownames(dataset)){
        if(!grepl(taxa_rank, row)){
            taxa_to_remove = c(taxa_to_remove, row)
        }
    }
    dataset_filtered = dataset[!rownames(dataset) %in% taxa_to_remove, ]
    print(paste0("The original dataset had ", nrow(dataset), " taxa abn_data and ", ncol(dataset), " samples."))
    print(paste0("The dataset has been filtered to keep only the taxa at the rank of ", taxa_rank, ". The new dataset has ", nrow(dataset_filtered), " rows and ", ncol(dataset_filtered), " columns."))
    print(paste0("Reduction of ", round(ceiling((nrow(dataset) - nrow(dataset_filtered))/nrow(dataset)*100)), "% of the original dataset."))
    return(dataset_filtered)
}

flagged_samples = c("PD_199P_profile","PD_107P_profile","PD_81P_profile","PD_218S_profile","PD_36P_profile","PD_192S_profile","PD_26P_profile","PD_112S_profile",
                    "PD_212S_profile","PD_216S_profile","PD_191P_profile","PD_211P_profile","PD_83P_profile","PD_217P_profile","PD_182S_profile","PD_161S_profile",
                    "PD_50P_profile","PD_28P_profile")

filtered_list = excluding_samples(dataset_calgary, metadata_calgary, TRUE, flagged_samples)  #This will remove the samples that were originally excluded from the Calgary cohort
metadata_calgary_f = filtered_list$metadata_calgary_f
dataset_calgary_f = filtered_list$dataset_calgary_f

metadata_calgary_f$Sex_Male1_Female2 = sapply(metadata_calgary_f$Sex, function(x) ifelse(x == 1, "Male", "Female"))  #This will create a new column in the metadata that has 1 for Male and 2 for

print(dim(dataset_calgary_f)) # 156 samples in total
print(dim(metadata_calgary_f)) # 156 samples in total

df_species_yyc = filtering_dataset_by_taxa_rank(dataset_calgary_f, "s__")  #This will filter the dataset to keep only the species-level taxa

# Main function to run ANCOM-BC2
ancombc2_yyc_results1 = ancombc2(data = df_species_yyc,
                            taxa_are_rows = TRUE,
                            fix_formula = "Group + Age + Sex_Male1_Female2",
                            p_adj_method = "BH",
                            meta_data = metadata_calgary_f,
                            alpha = 0.05,
                            prv_cut = 0.10,
                            verbose = TRUE)

ancombc2_yyc_results2 = ancombc2(data = df_species_yyc,
                            taxa_are_rows = TRUE,
                            fix_formula = "Group + Sex_Male1_Female2",
                            p_adj_method = "BH",
                            meta_data = metadata_calgary_f,
                            alpha = 0.05,
                            prv_cut = 0.10,
                            verbose = TRUE)


adding_metrics = function(ancombc2_results){
    results_res = ancombc2_results$res
    results_res$log2FC_PD = sapply(results_res$lfc_GroupPD, function(x) x = log2(exp(x)))
    #results_res$log2FC_Age = sapply(results_res$lfc_Age, function(x) x = log2(exp(x)))
    results_res$log2FC_Sex = sapply(results_res$lfc_Sex_Male1_Female2Male, function(x) x = log2(exp(x)))
    results_res$NegLog10qvalPD = sapply(results_res$q_GroupPD, function(x) x = -log10(x))
    #results_res$NegLog10qvalAge = sapply(results_res$q_Age, function(x) x = -log10(x))
    results_res$NegLog10qvalSex = sapply(results_res$q_Sex_Male1_Female2, function(x) x = -log10(x))
    results_res$Specific_taxa = sapply(results_res$taxon, function(x) strsplit(x, ";s__")[[1]][2])  #This will extract the species name from the taxon column, which has the format "k__Bacteria;p__Firmicutes;c__Clostridia;o__Clostridiales;f__Lachnospiraceae;g__Roseburia;s__intestinalis"
    return(results_res)
}
yyc_results1_res = adding_metrics(ancombc2_yyc_results1)
yyc_results2_res = adding_metrics(ancombc2_yyc_results2)



write.csv(yyc_results1_res, "/Users/danielcm/Desktop/MIGUT/vancouver/ancombc/calgary/ancombc2_results_calgary_species_level_prev10_sex_pd_age.csv", row.names = TRUE)
write.csv(yyc_results2_res, "/Users/danielcm/Desktop/MIGUT/vancouver/ancombc/calgary/ancombc2_results_calgary_species_level_prev10_sex_pd.csv", row.names = TRUE)

