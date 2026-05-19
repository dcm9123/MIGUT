library(Maaslin2)
library(ggplot2)

setwd("/Users/danielcm/Desktop/MIGUT/novaseq/")
maaslin_input = read.table("metaphlan_merged_gtdb_kevin/merged_abundances_kevin_script/species_filtered_gtdb_prev_10_mean_00001.txt", header = TRUE, sep = "\t", row.names = 1)
metadata_df = read.csv("encoded_metadata.txt", header = TRUE, sep = "\t")
metadata_df = as.data.frame(metadata_df)

humann_data = read.csv("final_results_raw/merged_pathabundance_unstratified_samples_only.csv", header = TRUE, row.names = 1)
humann_data_f = read.csv("final_results_raw/merged_pathabundance_unstratified_samples_filtered.csv", header = TRUE, row.names = 1)
humann_data_semi_f = read.csv("final_results_raw/merged_pathabundance_unstratified_samples_only_all_classified.csv", header = TRUE, row.names = 1)
rownames(metadata_df) <- metadata_df$ID
transformation_method = "LOG"
normalization_to_use = "TSS"
random_to_use = NULL


#Make sure the IDs match. Use the IDs as rownames of the metadata file so it matches the colnames of the input file
ids_input_file = colnames(maaslin_input)
ids_humann_file = colnames(humann_data)
#ids_humann_file
#ids_input_file

#Adjusting sex as a categorical variable by changing the 0 and 1 coding
metadata_df$Sex_male1_female2 <- factor(metadata_df$Sex_male1_female2,
                          levels = c(1, 2),
                          labels = c("Male", "Female"))

metadata_df$Constipation_no_0_yes_1 <- factor(metadata_df$Constipation_no_0_yes_1,
                          levels = c(0, 1),
                          labels = c("No", "Yes"))

calling_maaslin2 = function(maaslin_input, metadata_df, output_file, fixed_effects_to_use, transformation_method, normalization_to_use, random_to_use, reference_value){
    x = Maaslin2(input_data = maaslin_input, input_metadata = metadata_df, output=output_file, # nolint
        fixed_effects = fixed_effects_to_use, transform = transformation_method, 
        normalization = normalization_to_use, random_effects = random_to_use, reference = reference_value,
        cores = 6, plot_scatter = FALSE)
}

# Subsetting metadata


metadata_male_only = metadata_df[metadata_df$Sex_male1_female2 == "Male", ]

metadata_female_only = metadata_df[metadata_df$Sex_male1_female2 == "Female", ]

metadata_parkinsons_only = metadata_df[metadata_df$Group == "PD", ]

metadata_healthy_only = metadata_df[metadata_df$Group == "HC", ]

metadata_parkinsons_only_hfiltered = metadata_parkinsons_only[metadata_parkinsons_only$H_and_Y_stage != "", ]
metadata_parkinsons_only_hfiltered = metadata_parkinsons_only_hfiltered[!is.na(metadata_parkinsons_only_hfiltered$H_and_Y_stage), ]
metadata_parkinsons_only_hfiltered$H_and_Y_stage

metadata_parkinsons_levodopa = metadata_parkinsons_only[metadata_parkinsons_only$Levodopa_no_0_yes_1 == "Yes", ]

reference_value = NULL

pd_results = calling_maaslin2(maaslin_input, metadata_df, output_file = "maaslin2/pd_vs_no_pd", fixed_effects_to_use = "Group",     
                    transformation_method = transformation_method, normalization_to_use = normalization_to_use, random_to_use = random_to_use, 
                    reference_value = reference_value)

ps_sex_results = calling_maaslin2(maaslin_input, metadata_df, output_file = "maaslin2/pd_and_sex", fixed_effects_to_use = "Group_and_sex",
                    transformation_method = transformation_method, normalization_to_use = normalization_to_use, random_to_use = random_to_use, 
                    reference_value = "Group_and_sex,HC male")


pd_male_only_results = calling_maaslin2(maaslin_input, metadata_male_only, output_file = "maaslin2/pd_and_sex_male", fixed_effects_to_use = "Males_only",
                    transformation_method = transformation_method, normalization_to_use = normalization_to_use, random_to_use = random_to_use, 
                    reference_value = "HC male")


pd_female_only_results = calling_maaslin2(maaslin_input, metadata_female_only, output_file = "maaslin2/pd_and_sex_female", fixed_effects_to_use = "Females_only",
                    transformation_method = transformation_method, normalization_to_use = normalization_to_use, random_to_use = random_to_use, 
                    reference_value = "HC female")


pd_parkinsons_only_results = calling_maaslin2(maaslin_input, metadata_parkinsons_only, output_file = "maaslin2/pd_only", fixed_effects_to_use = "Sex_male1_female2",
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
