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
reading_input_files = function(){
    dataset_calgary = as.data.frame(read.csv("/Users/danielcm/Desktop/MIGUT/analysis_may2026/results1/merged_tables_abundances/merged_absolute_abundances_table_GTDB.txt", sep = "\t", row.names = 1))
    metadata_calgary = as.data.frame(read.csv("/Users/danielcm/Desktop/MIGUT/metadata/MIGUT_metadata_in_tab_format.txt", sep = "\t"))
    rownames(metadata_calgary) = metadata_calgary$ID
    dataset_vancouver = as.data.frame(read.csv("/Users/danielcm/Desktop/MIGUT/vancouver/merged_tables/gtdb_absolute_abundance_merged.txt", sep = "\t", row.names = 1))
    metadata_vancouver = as.data.frame(read.csv("/Users/danielcm/Desktop/MIGUT/vancouver/vancouver_participant_metadata.csv", sep = ",", row.names = 1))
    metadata_vancouver$Location = "Vancouver"
    metadata_calgary$Location = "Calgary"
}

# Removing samples that were originally excluded for any reason. This assumes that the metadata has them all but not the dataset count table, so the latter is a subset of the metadata.   
excluding_samples = function(dataset_calgary, metadata_calgary){
    samples_to_exclude = c()                                            #There are some samples that I need to exclude from the Calgary metadata as they were very low on reads 
    for(item in colnames(dataset_calgary)){                             #Go over all the name of the samples in the ASV count table
        if(!item %in% rownames(metadata_calgary)){                      #If one of them is not present in the metadata, then save it in the vector and carry on with the rest
            samples_to_exclude = c(samples_to_exclude, item)
        }
    }

    dataset_calgary_f = dataset_calgary[, !colnames(dataset_calgary) %in% samples_to_exclude]       #Remove the samples that are not good for the analysis (just a sanity check here)
    metadata_calgary_f = metadata_calgary[!rownames(metadata_calgary) %in% samples_to_exclude, ]    #Remove the samples that are not good fot the analysis, from all metadata to just the ones of interest

    return(list(dataset_calgary_f = dataset_calgary_f,
                metadata_calgary_f = metadata_calgary_f))
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


df_samples_excluded = excluding_samples(dataset_calgary, metadata_calgary)[[1]]
df_metadata_calgary_f = excluding_samples(dataset_calgary, metadata_calgary)[[2]]
df_species = filtering_dataset_by_taxa_rank(df_samples_excluded, "s__")  #This will filter the dataset to keep only the species-level taxa


# Main function to run ANCOM-BC2
ancombc2_results = ancombc2(data = df_species,
                            taxa_are_rows = TRUE,
                            fix_formula = "Group + Age + Sex_Male1_Female2",
                            p_adj_method = "BH",
                            #struc_zero = TRUE,
                            meta_data = df_metadata_calgary_f,
                            group = "Group",
                            alpha = 0.05,
                            global = TRUE,
                            prv_cut = 0.10,
                            verbose = TRUE)

ancombc2_results_2 = ancombc2(data = df_species,
                            taxa_are_rows = TRUE,
                            fix_formula = "Group + Age + Sex_Male1_Female2",
                            p_adj_method = "BH",
                            #struc_zero = TRUE,
                            meta_data = df_metadata_calgary_f,
                            group = "Group",
                            alpha = 0.05,
                            global = TRUE,
                            prv_cut = 0.05,
                            verbose = TRUE)

ancombc2_results_3 = ancombc2(data = df_species,
                            taxa_are_rows = TRUE,
                            fix_formula = "Group",
                            p_adj_method = "BH",
                            #struc_zero = TRUE,
                            meta_data = df_metadata_calgary_f,
                            group = "Group",
                            alpha = 0.05,
                            global = TRUE,
                            prv_cut = 0.10,
                            verbose = TRUE)


write.csv(ancombc2_results$res, "/Users/danielcm/Desktop/MIGUT/vancouver/ancombc/calgary/ancombc2_results_calgary_species_level_prev10.csv", row.names = TRUE)
write.csv(ancombc2_results_2$res, "/Users/danielcm/Desktop/MIGUT/vancouver/ancombc/calgary/ancombc2_results_calgary_species_level_prev5.csv", row.names = TRUE)
write.csv(ancombc2_results_3$res, "/Users/danielcm/Desktop/MIGUT/vancouver/ancombc/calgary/ancombc2_results_calgary_species_level_prev10_only_diagnosis.csv", row.names = TRUE)