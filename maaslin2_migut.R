library(Maaslin2)
library(ggplot2)

setwd("/Users/danielcm/Desktop/MIGUT/novaseq/")
maaslin_input = read.table("metaphlan_merged_gtdb_kevin/merged_abundances_kevin_script/species_filtered_gtdb_prev_10_mean_00001.txt", header = TRUE, sep = "\t", row.names = 1)
metadata_df = read.csv("encoded_metadata.txt", header = TRUE, sep = "\t")
metadata_df = as.data.frame(metadata_df)
rownames(metadata_df) <- metadata_df$ID
transformation_method = "LOG"
normalization_to_use = "TSS"
random_to_use = NULL


#Make sure the IDs match. Use the IDs as rownames of the metadata file so it matches the colnames of the input file
ids_input_file = colnames(maaslin_input)
#ids_input_file

#Adjusting sex as a categorical variable by changing the 0 and 1 coding
metadata_df$Sex_male1_female2 <- factor(metadata_df$Sex_male1_female2,
                          levels = c(1, 2),
                          labels = c("Male", "Female"))

calling_maaslin2 = function(maaslin_input, metadata_df, output_file, fixed_effects_to_use, transformation_method, normalization_to_use, random_to_use, reference_value){
    x = Maaslin2(input_data = maaslin_input, input_metadata = metadata_df, output=output_file, # nolint
        fixed_effects = fixed_effects_to_use, transform = transformation_method, 
        normalization = normalization_to_use, random_effects = random_to_use, reference = reference_value,
        cores = 6, plot_scatter = FALSE)
}

reference_value = NULL

pd_results = calling_maaslin2(maaslin_input, metadata_df, output_file = "maaslin2/pd_vs_no_pd", fixed_effects_to_use = "Group",     
                    transformation_method = transformation_method, normalization_to_use = normalization_to_use, random_to_use = random_to_use, 
                    reference_value = reference_value)

ps_sex_results = calling_maaslin2(maaslin_input, metadata_df, output_file = "maaslin2/pd_and_sex", fixed_effects_to_use = "Group_and_sex",
                    transformation_method = transformation_method, normalization_to_use = normalization_to_use, random_to_use = random_to_use, 
                    reference_value = "Group_and_sex,HC male")

metadata_male_only = metadata_df[metadata_df$Sex_male1_female2 == 1, ]

rownames(metadata_male_only)

pd_male_only_results = calling_maaslin2(maaslin_input, metadata_male_only, output_file = "maaslin2/pd_and_sex_male", fixed_effects_to_use = "Males_only",
                    transformation_method = transformation_method, normalization_to_use = normalization_to_use, random_to_use = random_to_use, 
                    reference_value = "HC male")
