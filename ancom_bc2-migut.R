# Daniel Castaneda Mogollon, PhD
# July 15th, 2026
# Purpose: This script will run ANCOM-BC2 on microbiome data sets

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("ANCOMBC")

library(ANCOMBC)
#library(phyloseq)
library(dplyr)

library(tidyr); library(dplyr); library(tibble)


path = "/Users/danielcm/Desktop/MIGUT/vancouver/ancombc"
setwd(path)
# Reading the main files to work with, that is, ASV count tables and metadata files.    
dataset_calgary = as.data.frame(read.csv("/Users/danielcm/Desktop/MIGUT/analysis_may2026/results1/merged_tables_abundances/merged_absolute_abundances_table_GTDB.txt", sep = "\t", row.names = 1))
metadata_calgary = as.data.frame(read.csv("/Users/danielcm/Desktop/MIGUT/metadata/MIGUT_metadata_in_tab_format.txt", sep = "\t"))
rownames(metadata_calgary) = metadata_calgary$ID
dataset_vancouver = as.data.frame(read.csv("/Users/danielcm/Desktop/MIGUT/vancouver/merged_tables/gtdb_absolute_abundance_merged.txt", sep = "\t", row.names = 1))
metadata_vancouver = as.data.frame(read.csv("/Users/danielcm/Desktop/MIGUT/vancouver/vancouver_participant_metadata.csv", sep = ",", row.names = 1))
colnames(dataset_vancouver) = sapply(strsplit(colnames(dataset_vancouver), "_"), `[`, 1)  # Removes the "_S123_gtdb_profile" part so it matches the metadata's sample names

metadata_vancouver = metadata_vancouver %>% rename(Group = Disease)  #This will rename the "Disease" column to "Group" in the metadata_vancouver dataframe
metadata_vancouver = metadata_vancouver %>% rename(Constipation_No0_Yes1 = Const)
metadata_vancouver = metadata_vancouver %>% rename(UPDRS.TOTAL = updrstotal)
metadata_vancouver = metadata_vancouver %>% rename(Sex_Male1_Female2 = Female_1)
metadata_vancouver$Sex_Male1_Female2 = sapply(metadata_vancouver$Sex_Male1_Female2, function(x) ifelse(x == 1, 2, 1))  #This will change the coding
metadata_vancouver$ID = rownames(metadata_vancouver)
metadata_vancouver$Group = sapply(metadata_vancouver$Group, function(x) ifelse(x == "Control", "HC", "PD"))  #This will change the coding of the Group column to match the Calgary cohort

metadata_calgary = metadata_calgary %>% rename(updrs1total = `Part.1.Total`)  
metadata_calgary = metadata_calgary %>% rename(updrs2total = `Part.2.Total`) 
metadata_calgary = metadata_calgary %>% rename(updrs3total = `Part.3.Total`)  
metadata_calgary = metadata_calgary %>% rename(updrs4total = `Part.4.Total`)  

metadata_vancouver$Sex_Male1_Female2 = sapply(metadata_vancouver$Sex_Male1_Female2, function(x) ifelse(x == 1, "Male", "Female"))  #This will create a new column
#metadata_calgary$Sex_Male1_Female2 = as.integer(metadata_calgary$Sex_Male1_Female2)  #This will convert the Sex column to integer
metadata_calgary$Sex_Male1_Female2 = sapply(metadata_calgary$Sex_Male1_Female2, function(x) ifelse(x == 1, "Male", "Female"))  #This will create a new column in the metadata that has 1 for Male and 2 for

metadata_vancouver$location = "Vancouver"
metadata_calgary$location = "Calgary"

dataset_ab_bc = merge(dataset_calgary, dataset_vancouver, by = 0, all = TRUE)
dataset_ab_bc[is.na(dataset_ab_bc)] = 0  #This will replace all the NA values with 0, since these are taxa that are not present in the sample
rownames(dataset_ab_bc) = dataset_ab_bc$Row.names

dim(dataset_ab_bc)  # 455 samples in total from both locations 

# Originally, there were 184 samples in the Calgary cohort metadata's file, but 10 of these were not sequenced (17, 18, 19, 20, 21, 22, 23, 24, 168 and 187)
# Two sample IDs that were sequenced do not have a match in the metadata file from Davide: 40, 70
# 18 samples sequenced had too few reads and were removed from the analysis: 199, 212, 50, 107, 216, 28, 81, 191, 218, 211, 36, 83, 192, 217, 26, 182, 112, 161

# From 184 to 

#print(dim(dataset_calgary)) # 176 samples in total
#print(dim(metadata_calgary)) # 184 samples in total. 


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
        print(row)
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

dim(dataset_calgary_f) # 156 samples in total

##
dataset_ab_bc = merge(dataset_calgary_f, dataset_vancouver, by = 0, all = TRUE)
rownames(dataset_ab_bc) = dataset_ab_bc$Row.names
dataset_ab_bc$Row.names = NULL                 # (1) drop the leftover column
dataset_ab_bc[is.na(dataset_ab_bc)] = 0

metadata_ab_bc = dplyr::full_join(
    metadata_calgary_f, metadata_vancouver,
    by = c("ID","Age","Sex_Male1_Female2","Group",
           "Constipation_No0_Yes1","UPDRS.TOTAL","location", "updrs1total", "updrs2total", "updrs3total", "updrs4total"))
rownames(metadata_ab_bc) = metadata_ab_bc$ID   # (2) restore names for matching


pd_subset_merged_metadata = metadata_ab_bc[metadata_ab_bc$Group == "PD", ]
pd_subset_merged_healthy = metadata_ab_bc[metadata_ab_bc$Group == "HC", ]
pd_subset_data_merged = dataset_ab_bc[, colnames(dataset_ab_bc) %in% pd_subset_merged_metadata$ID]
pd_subset_data_healthy = dataset_ab_bc[, colnames(dataset_ab_bc) %in% pd_subset_merged_healthy$ID]


df_species_bc_ab = filtering_dataset_by_taxa_rank(dataset_ab_bc, "s__")
df_species_bc_ab_pd_only = filtering_dataset_by_taxa_rank(pd_subset_data_merged, "s__")
df_species_bc_ab_healthy_only = filtering_dataset_by_taxa_rank(pd_subset_data_healthy, "s__")

##


#View(dataset_ab_bc_f)

print(dim(dataset_calgary_f)) # 156 samples in total
print(dim(metadata_calgary_f)) # 156 samples in total

df_species_yyc = filtering_dataset_by_taxa_rank(dataset_calgary_f, "s__")  #This will filter the dataset to keep only the species-level taxa
df_species_yvr = filtering_dataset_by_taxa_rank(dataset_vancouver, "s__")  #This will filter the dataset to keep only the species-level taxa
df_species_bc_ab = filtering_dataset_by_taxa_rank(dataset_ab_bc, "s__")  #This will filter the dataset to keep only the species-level taxa

# Main function to run ANCOM-BC2
ancombc2_yyc_results1 = ancombc2(data = df_species_yyc,
                            taxa_are_rows = TRUE,
                            fix_formula = "Group + Age + Sex_Male1_Female2",
                            p_adj_method = "BH",
                            meta_data = metadata_calgary_f,
                            alpha = 0.05,
                            prv_cut = 0.10,
                            struc_zero = TRUE,
                            group = "Group",
                            verbose = TRUE)

ancombc2_yyc_results2 = ancombc2(data = df_species_yyc,
                            taxa_are_rows = TRUE,
                            fix_formula = "Group + Sex_Male1_Female2",
                            p_adj_method = "BH",
                            meta_data = metadata_calgary_f,
                            alpha = 0.05,
                            prv_cut = 0.10,
                            verbose = TRUE)

ancombc_yvr_results1 = ancombc2(data = df_species_yvr,
                            taxa_are_rows = TRUE,
                            fix_formula = "Group + Age + Sex_Male1_Female2",
                            p_adj_method = "BH",
                            meta_data = metadata_vancouver,
                            alpha = 0.05,
                            prv_cut = 0.10,
                            struc_zero = FALSE,
                            verbose = TRUE)

ancombc_yvr_results2 = ancombc2(data = df_species_yvr,
                            taxa_are_rows = TRUE,
                            fix_formula = "Group + Sex_Male1_Female2",
                            p_adj_method = "BH",
                            meta_data = metadata_vancouver,
                            alpha = 0.05,
                            prv_cut = 0.10,
                            struc_zero = FALSE,
                            verbose = TRUE)

ancombc_yvr_results3 = ancombc2(data = df_species_yvr,
                            taxa_are_rows = TRUE,
                            fix_formula = "Group + Sex_Male1_Female2 + Age",
                            p_adj_method = "BH",
                            meta_data = metadata_vancouver,
                            alpha = 0.05,
                            prv_cut = 0.10,
                            struc_zero = FALSE,
                            verbose = TRUE)

ancombc_yvr_results4 = ancombc2(data = df_species_yvr,
                            taxa_are_rows = TRUE,
                            fix_formula = "Group + Sex_Male1_Female2",
                            p_adj_method = "BH",
                            meta_data = metadata_vancouver,
                            alpha = 0.05,
                            prv_cut = 0.10,
                            struc_zero = FALSE,
                            verbose = TRUE)



# This one runs the stratified PD across AB and BC cohorts, with location as a covariate, age, and sex
ancombc_bc_ab_results_by_city_in_pd = ancombc2(data = df_species_bc_ab_pd_only,
                            taxa_are_rows = TRUE,
                            fix_formula = "Age +Sex_Male1_Female2 + location",
                            p_adj_method = "BH",
                            meta_data = pd_subset_merged_metadata,
                            alpha = 0.05,
                            prv_cut = 0.10,
                            struc_zero = FALSE,
                            verbose = TRUE) 

ancombc_bc_ab_results_by_city_in_healthy = ancombc2(data = df_species_bc_ab_healthy_only,
                            taxa_are_rows = TRUE,
                            fix_formula = "Age + Sex_Male1_Female2 + location",
                            p_adj_method = "BH",
                            meta_data = pd_subset_merged_healthy,
                            alpha = 0.05,
                            prv_cut = 0.10,
                            struc_zero = FALSE,
                            verbose = TRUE) 

ancombc_bc_ab_results_diagnosis_and_interaction = ancombc2(data = df_species_bc_ab,
                            taxa_are_rows = TRUE,
                            fix_formula = "Group * location + Age + Sex_Male1_Female2",
                            p_adj_method = "BH",
                            meta_data = metadata_ab_bc,
                            alpha = 0.05,
                            prv_cut = 0.10,
                            struc_zero = FALSE,
                            verbose = TRUE)

ancombc_bc_ab_results_diagnosis_and_interaction2 = ancombc2(data = df_species_bc_ab,
                            taxa_are_rows = TRUE,
                            fix_formula = "Group * location + Sex_Male1_Female2",
                            p_adj_method = "BH",
                            meta_data = metadata_ab_bc,
                            alpha = 0.05,
                            prv_cut = 0.10,
                            struc_zero = FALSE,
                            verbose = TRUE)

dim(ancombc_bc_ab_results_diagnosis_and_interaction2$bias_correct_log_table)


adding_metrics = function(ancombc2_results, merged){
    results_res = ancombc2_results$res

    if(merged == TRUE){
        results_res$log2FC_location = sapply(results_res$lfc_locationVancouver, function(x) x = log2(exp(x)))
        results_res$NegLog10qval_location = sapply(results_res$q_locationVancouver, function(x) x = -log10(x))
        results_res$Significant_location = ifelse(results_res$q_locationVancouver > 0.05, "No evidence of enrichment",
        ifelse(abs(results_res$log2FC_location) > 1 & results_res$q_locationVancouver < 0.05, "Enriched",
        ifelse(abs(results_res$log2FC_location) > 0.585 & abs(results_res$log2FC_location) <= 1 & results_res$q_locationVancouver < 0.05, "Potentially enriched",
        "No evidence of enrichment")))
    }

    results_res$log2FC_PD = sapply(results_res$lfc_GroupPD, function(x) x = log2(exp(x)))
    #results_res$log2FC_Age = sapply(results_res$lfc_Age, function(x) x = log2(exp(x)))
    results_res$log2FC_Sex = sapply(results_res$lfc_Sex_Male1_Female2Male, function(x) x = log2(exp(x)))
    results_res$NegLog10qvalPD = sapply(results_res$q_GroupPD, function(x) x = -log10(x))
    #results_res$NegLog10qvalAge = sapply(results_res$q_Age, function(x) x = -log10(x))
    results_res$NegLog10qvalSex = sapply(results_res$q_Sex_Male1_Female2, function(x) x = -log10(x))
    results_res$Specific_taxa = sapply(results_res$taxon, function(x) strsplit(x, ";s__")[[1]][2])  #This will extract the species name from the taxon column, which has the format "k__Bacteria;p__Firmicutes;c__Clostridia;o__Clostridiales;f__Lachnospiraceae;g__Roseburia;s__intestinalis"
    results_res$SignificantPD = ifelse(results_res$q_GroupPD > 0.05, "No evidence of enrichment", 
    ifelse(abs(results_res$log2FC_PD) > 1 & results_res$q_GroupPD < 0.05, "Enriched",
    ifelse(abs(results_res$log2FC_PD) > 0.585 & abs(results_res$log2FC_PD) <= 1 & results_res$q_GroupPD < 0.05, "Potentially enriched",
    "No evidence of enrichment")))
    return(results_res)
}

yyc_results1_res = adding_metrics(ancombc2_yyc_results1, merged = FALSE)
yyc_results2_res = adding_metrics(ancombc2_yyc_results2, merged = FALSE)
yvr_results1_res = adding_metrics(ancombc_yvr_results1, merged = FALSE)
yvr_results2_res = adding_metrics(ancombc_yvr_results2, merged = FALSE)
yvr_results3_res = adding_metrics(ancombc_yvr_results3, merged = FALSE)
yvr_results4_res = adding_metrics(ancombc_yvr_results4, merged = FALSE)
    
ab_bc_results_by_city_in_pd_res = adding_metrics(ancombc_bc_ab_results_by_city_in_pd, merged = TRUE)
ab_bc_results_by_city_in_healthy_res = adding_metrics(ancombc_bc_ab_results_by_city_in_healthy, merged = TRUE)
ab_bc_results_diagnosis_and_interaction_res = adding_metrics(ancombc_bc_ab_results_diagnosis_and_interaction, merged = TRUE)
ab_bc_results_diagnosis_and_interaction2_res = adding_metrics(ancombc_bc_ab_results_diagnosis_and_interaction2, merged = TRUE)


write.csv(yyc_results1_res, "/Users/danielcm/Desktop/MIGUT/vancouver/ancombc/calgary/ancombc2_results_calgary_species_level_prev10_sex_pd_age.csv", row.names = TRUE)
write.csv(yyc_results2_res, "/Users/danielcm/Desktop/MIGUT/vancouver/ancombc/calgary/ancombc2_results_calgary_species_level_prev10_sex_pd.csv", row.names = TRUE)
write.csv(yvr_results1_res, "/Users/danielcm/Desktop/MIGUT/vancouver/ancombc/ancombc2_results_vancouver_species_level_prev10_sex_pd_age.csv", row.names = TRUE)
write.csv(yvr_results2_res, "/Users/danielcm/Desktop/MIGUT/vancouver/ancombc/ancombc2_results_vancouver_species_level_prev10_sex_pd.csv", row.names = TRUE)
write.csv(yvr_results3_res, "/Users/danielcm/Desktop/MIGUT/vancouver/ancombc/ancombc2_results_vancouver_species_level_prev10_sex_age.csv", row.names = TRUE)    
write.csv(yvr_results4_res, "/Users/danielcm/Desktop/MIGUT/vancouver/ancombc/ancombc2_results_vancouver_species_level_prev10_sex_pd.csv", row.names = TRUE)
write.csv(ab_bc_results_by_city_in_pd_res, "/Users/danielcm/Desktop/MIGUT/vancouver/ancombc/ancombc2_results_merged_species_level_prev10_sex_location_age_pd_stratified.csv", row.names = TRUE)
write.csv(ab_bc_results_by_city_in_healthy_res, "/Users/danielcm/Desktop/MIGUT/vancouver/ancombc/ancombc2_results_merged_species_level_prev10_sex_location_hc_stratified.csv", row.names = TRUE)
write.csv(ab_bc_results_diagnosis_and_interaction_res, "/Users/danielcm/Desktop/MIGUT/vancouver/ancombc/ancombc2_results_merged_species_level_prev10_sex_location_age_pd_interaction.csv", row.names = TRUE)
write.csv(ab_bc_results_diagnosis_and_interaction2_res, "/Users/danielcm/Desktop/MIGUT/vancouver/ancombc/ancombc2_results_merged_species_level_prev10_sex_location_pd_interaction.csv", row.names = TRUE)

### Getting Spearman Correlation coefficients

clr_merged_data = ancombc_bc_ab_results_diagnosis_and_interaction2$bias_correct_log_table #CLR corrected table

#View(bc)
metadata_pd_ab_bc = metadata_ab_bc[metadata_ab_bc$Group == "PD", ]
metadata_pd_ab_bc = metadata_pd_ab_bc %>% arrange(ID) 
yvr_metadata = metadata_pd_ab_bc %>% filter(location == "Vancouver" | location == "Calgary")
yvr_dataset = clr_merged_data %>% select(all_of(yvr_metadata$ID))
updrs_metadata_vancouver = yvr_metadata %>% select(ID,updrs1total,updrs2total,updrs3total,updrs4total)
rownames(updrs_metadata_vancouver) = NULL

dim(yvr_dataset)  # 105 samples in total (All PD for YYC) # 197 for Vancouver

x = filtering_dataset_by_taxa_rank(yvr_dataset, "s__")  #Sanity check, there should be a reduction of 0%, but just in case ...

### This block allows me to remove the rest of the taxonomic rank and only keep the species. It also allows me to format it in a way that I can 
### Generate a pearson/spearman correlation between a variable of interest (UPDRS) and CLR of species across samples in PD

yvr_dataset$species = sapply(rownames(yvr_dataset), function(x) strsplit(x, ";s__")[[1]][2])  #This will remove the "k__Bacteria;p__Firmicutes;c__Clostridia;o__Clostridiales;f__Lachnospiraceae;g__Roseburia;s__" part of the taxa names, so only the species name remains
yvr_dataset = yvr_dataset %>% relocate(species)  #This will move the species column to the first column of the dataset
yvr_dataset = t(yvr_dataset)
yvr_dataset = as.data.frame(yvr_dataset) #This will check the class of the dataset, which should be a matrix
colnames(yvr_dataset) = yvr_dataset[1,]  #This will set the first row of the dataset as the column names
yvr_dataset = yvr_dataset[-1,]  #This will remove the first row of the dataset, which is now redundant since it has been set as the column names
yvr_dataset$ID = rownames(yvr_dataset)
rownames(yvr_dataset) = NULL  #This will remove the row names of the dataset, which are now redundant since they have been set as the ID column
yvr_dataset = yvr_dataset %>% relocate(ID)  #This will move the ID column to the first column of the dataset
yvr_dataset = yvr_dataset[, -42] # Removing an NA column with no species assigned
yvr_dataset = yvr_dataset[, -119] # Removing an NA column with no species assigned

species_and_updrs_yyc = full_join(updrs_metadata_vancouver, yvr_dataset, by = "ID")  #This will join the two data frames by the ID column, which is the sample ID
na_counts = colSums(is.na(species_and_updrs_yyc))  #This will count the number of NA values in each column of the data frame
dim(species_and_updrs_yyc)  #This will check the dimensions of the data frame, which should be 105 samples and 4 UPDRS columns + 1 ID column + 119 species columns = 124 columns

keep_metadata = c("ID","updrs1total","updrs2total","updrs3total","updrs4total")
species_cols = setdiff(colnames(species_and_updrs_yyc), keep_metadata)  #This will create a vector of the column names that are not in the metadata columns
good_species = species_cols[colMeans(is.na(species_and_updrs_yyc[species_cols])) <= 1]  #This will create a vector of the column names that have less than 70% NA values
species_and_updrs_yyc_clean = species_and_updrs_yyc %>%
                            select(all_of(c(keep_metadata, good_species)))

rownames(species_and_updrs_yyc_clean) = NULL 

calculating_pearson = function(dataset, good_species, updrs_col){
  updrs_values = as.numeric(dataset[[updrs_col]])
  res = do.call(rbind, lapply(good_species, function(col){
    taxa_values = as.numeric(dataset[[col]])
    ok = is.finite(taxa_values) & is.finite(updrs_values)
    if(sum(ok) < 10) return(NULL)
    ct = cor.test(taxa_values[ok], updrs_values[ok], method = "pearson")
    data.frame(Species = col,
               Pearson_Correlation = unname(ct$estimate),
               P_value = ct$p.value)
  }))
  res$Q_value = p.adjust(res$P_value, method = "BH")   # correct across all taxa
  rownames(res) = NULL
  return(res)
}


### CHECK AGAIN THIS PART TO SEE IF IT'S WORTH IT

updrs1_pearson_results = calculating_pearson(species_and_updrs_yyc_clean, good_species, updrs_col = "updrs1total") #no significant q-vals
updrs2_pearson_results = calculating_pearson(species_and_updrs_yyc_clean, good_species, updrs_col = "updrs2total") #no significant q-vals
updrs3_pearson_results = calculating_pearson(species_and_updrs_yyc_clean, good_species, updrs_col = "updrs3total") #
updrs4_pearson_results = calculating_pearson(species_and_updrs_yyc_clean, good_species, updrs_col = "updrs4total")  


View(updrs2_pearson_results)

View(species_and_updrs_yyc_clean$`Hominicoprocola sp900317525`)
