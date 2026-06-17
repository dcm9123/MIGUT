BiocManager::install("devtools")
BiocManager::install("remotes")
BiocManager::install("biobakery/maaslin3")

library(maaslin3)
library(ggplot2)

path = "/Users/danielcm/Desktop/MIGUT"
setwd(path)

maaslin_input = read.table(paste0(path,"/analysis_may2026/filtered_metaphlan4_abundances/merged_absolute_abundances_table_GTDB_filtered.txt"), header = TRUE, sep = "\t", row.names = 1)
metadata_df = read.csv(paste0(path,"/metadata","/MIGUT_metadata_in_tab_format.txt"), header = TRUE, sep = "\t")
metadata_df = as.data.frame(metadata_df)

#humann_data = read.csv("final_results_raw/merged_pathabundance_unstratified_samples_only.csv", header = TRUE, row.names = 1)
#humann_data_f = read.csv("final_results_raw/merged_pathabundance_unstratified_samples_filtered.csv", header = TRUE, row.names = 1)
#humann_data_semi_f = read.csv("final_results_raw/merged_pathabundance_unstratified_samples_only_all_classified.csv", header = TRUE, row.names = 1)
rownames(metadata_df) <- metadata_df$ID
transformation_method = "LOG"
normalization_to_use = "TSS"


#Make sure the IDs match. Use the IDs as rownames of the metadata file so it matches the colnames of the input file
ids_input_file = colnames(maaslin_input)

ids_input_file
ids_encoded = grepl("S_profile", ids_input_file)


#Adjusting sex as a categorical variable by changing the 0 and 1 coding
metadata_df$Sex_Male1_Female2 <- factor(metadata_df$Sex_Male1_Female2,
                          levels = c(1, 2),
                          labels = c("Male", "Female"))

metadata_df$Constipation_No0_Yes1 <- factor(metadata_df$Constipation_No0_Yes1,
                          levels = c(0, 1),
                          labels = c("No", "Yes"))

calling_maaslin3 = function(maaslin_input, metadata_df, output_dir, median_comparison_abundance, median_comparison_prevalence, formula){
    x = maaslin3(input_data = maaslin_input, input_metadata = metadata_df,
        output = output_dir, formula = formula, normalization = normalization_to_use,
        transform = transformation_method, augment = TRUE, standardize = TRUE, max_significance = 0.05,
        median_comparison_abundance = TRUE, median_comparison_prevalence = FALSE, max_pngs = 600, cores = 8)
}

# Subsetting metadata


metadata_male_only = metadata_df[metadata_df$Sex_Male1_Female2 == "Male", ]
metadata_female_only = metadata_df[metadata_df$Sex_Male1_Female2 == "Female", ]
metadata_parkinsons_only = metadata_df[metadata_df$Group == "PD", ]
metadata_healthy_only = metadata_df[metadata_df$Group == "HC", ]
metadata_parkinsons_only_hfiltered = metadata_parkinsons_only[metadata_parkinsons_only$H_and_Y_stage != "", ]
metadata_parkinsons_only_hfiltered = metadata_parkinsons_only_hfiltered[!is.na(metadata_parkinsons_only_hfiltered$H_and_Y_stage), ]
metadata_parkinsons_only_hfiltered$H_and_Y_stage
metadata_parkinsons_levodopa = metadata_parkinsons_only[metadata_parkinsons_only$Levodopa_no_0_yes_1 == "Yes", ]

reference_value = NULL

rownames(metadata_df) = metadata_df$ID

# PD vs HC - Abundance Yes, Prevalence No
pd_results = calling_maaslin3(maaslin_input = maaslin_input, metadata_df = metadata_df, output_dir = paste0(path,"/analysis_may2026/maaslin3/PD_vs_HC"), 
                                formula = "~ Group ", median_comparison_abundance = TRUE, 
                                median_comparison_prevalence = FALSE)

#formula = "~ Group + Sex_Male1_Female2 + Age + BMI + Constipation_No0_Yes1",
# Parkinson's vs Healthy only Males
pd_male_only_results = calling_maaslin2(maaslin_input, metadata_male_only, output_file = paste0(path,"/analysis_may2026/maaslin2/PD_vs_HC_in_males"), fixed_effects_to_use = "Males_only",
                    transformation_method = transformation_method, normalization_to_use = normalization_to_use, random_to_use = random_to_use, 
                    reference_value = "HC_Male")

# Parkinson's vs Healthy only Females
pd_female_only_results = calling_maaslin2(maaslin_input, metadata_female_only, output_file = paste0(path,"/analysis_may2026/maaslin2/PD_vs_HC_in_females"), fixed_effects_to_use = "Females_only",
                    transformation_method = transformation_method, normalization_to_use = normalization_to_use, random_to_use = random_to_use, 
                    reference_value = "HC_Female")

# Healthy ctrls males vs females
pd_sex_results = calling_maaslin2(maaslin_input, metadata_healthy_only, output_file = paste0(path,"/analysis_may2026/maaslin2/HC_males_vs_females"), fixed_effects_to_use = "Sex_Male1_Female2",
                    transformation_method = transformation_method, normalization_to_use = normalization_to_use, random_to_use = random_to_use, 
                    reference_value = "Male")

# Parkinson's males vs females
hc_sex_results = calling_maaslin2(maaslin_input, metadata_parkinsons_only, output_file = paste0(path,"/analysis_may2026/maaslin2/PD_males_vs_females"), fixed_effects_to_use = "Sex_Male1_Female2",
                    transformation_method = transformation_method, normalization_to_use = normalization_to_use, random_to_use = random_to_use, 
                    reference_value = "Male")









pd_parkinsons_only_results = calling_maaslin2(maaslin_input, metadata_parkinsons_only, output_file = paste0(path,"/analysis_may2026/maaslin2/pd_only"), fixed_effects_to_use = "Sex_male1_female2",
                    transformation_method = transformation_method, normalization_to_use = normalization_to_use, random_to_use = random_to_use, 
                    reference_value = "Female")


hc_parkinsons_only_results = calling_maaslin2(maaslin_input, metadata_healthy_only, output_file = "maaslin2/hc_only", fixed_effects_to_use = "Sex_male1_female2",
                    transformation_method = transformation_method, normalization_to_use = normalization_to_use, random_to_use = random_to_use, 
                    reference_value = "Female")

parkinsons_constipation = calling_maaslin2(maaslin_input, metadata_parkinsons_only, output_file = "maaslin2/pd_constipation", fixed_effects_to_use = "Constipation_no_0_yes_1",
                    transformation_method = transformation_method, normalization_to_use = normalization_to_use, random_to_use = random_to_use, 
                    reference_value = "Yes")

h_and_y_results = calling_maaslin2(maaslin_input, metadata_parkinsons_only_hfiltered, output_file = "maaslin2/pd_h_and_y", fixed_effects_to_use = "H_and_Y_stage",
                    transformation_method = transformation_method, normalization_to_use = normalization_to_use, random_to_use = random_to_use, 
                    reference_value = "H_and_Y_stage,4")

updrs_total_results = calling_maaslin2(maaslin_input, metadata_parkinsons_only, output_file = "maaslin2/pd_updrs_total", fixed_effects_to_use = "UPDRS_Total",
                    transformation_method = transformation_method, normalization_to_use = normalization_to_use, random_to_use = random_to_use, 
                    reference_value = NULL)

humann_pd_results = calling_maaslin2(humann_data, metadata_df, output_file = "maaslin2/humann_pd_vs_no_pd", fixed_effects_to_use = "Group",     
                    transformation_method = transformation_method, normalization_to_use = normalization_to_use, random_to_use = random_to_use, 
                    reference_value = reference_value)

human_pd_f_results = calling_maaslin2(humann_data, metadata_df, output_file = "maaslin2/humann_pd_and_no_pd_f", fixed_effects_to_use = "Group",
                    transformation_method = transformation_method, normalization_to_use = normalization_to_use, random_to_use = random_to_use, 
                    reference_value = reference_value)

human_pd_semi_f_results = calling_maaslin2(humann_data_semi_f, metadata_df, output_file = "maaslin2/humann_pd_and_no_pd_semi_f", fixed_effects_to_use = "Group",
                    transformation_method = transformation_method, normalization_to_use = normalization_to_use, random_to_use = random_to_use, 
                    reference_value = reference_value)
