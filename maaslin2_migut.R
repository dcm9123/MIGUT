BiocManager::install("devtools")
BiocManager::install("remotes")
BiocManager::install("biobakery/maaslin3")
BiocManager::install("biobakery/Maaslin2")


library(Maaslin2)
library(ggplot2)

path = "/Users/danielcm/Desktop/MIGUT"
setwd(path)

maaslin_input = read.table(paste0(path,"/analysis_may2026/results1/merged_tables_abundances/merged_absolute_abundances_table_GTDB.txt"), header = TRUE, sep = "\t", row.names = 1)
metadata_df = read.csv(paste0(path,"/metadata","/MIGUT_metadata_in_tab_format.txt"), header = TRUE, sep = "\t")
metadata_df = as.data.frame(metadata_df)


#humann_data = read.csv("final_results_raw/merged_pathabundance_unstratified_samples_only.csv", header = TRUE, row.names = 1)
#humann_data_f = read.csv("final_results_raw/merged_pathabundance_unstratified_samples_filtered.csv", header = TRUE, row.names = 1)
#humann_data_semi_f = read.csv("final_results_raw/merged_pathabundance_unstratified_samples_only_all_classified.csv", header = TRUE, row.names = 1)
rownames(metadata_df) <- metadata_df$ID
transformation_method = "LOG"
normalization_to_use = "TSS"

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

flagged_samples = c("PD_199P_profile","PD_107P_profile","PD_81P_profile","PD_218S_profile","PD_36P_profile","PD_192S_profile","PD_26P_profile","PD_112S_profile",
                    "PD_212S_profile","PD_216S_profile","PD_191P_profile","PD_211P_profile","PD_83P_profile","PD_217P_profile","PD_182S_profile","PD_161S_profile",
                    "PD_50P_profile","PD_28P_profile")


metadata_filtered = excluding_samples(maaslin_input, metadata_df, TRUE, flagged_samples)$metadata_calgary_f
maaslin_input_filtered = excluding_samples(maaslin_input, metadata_df, TRUE, flagged_samples)$dataset_calgary_f

ncol(maaslin_input_filtered) # 156, matches the same as ANCOM-BC2
#Make sure the IDs match. Use the IDs as rownames of the metadata file so it matches the colnames of the input file
ids_input_file = colnames(maaslin_input_filtered)

ids_input_file
ids_encoded = grepl("S_profile", ids_input_file)

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


#Adjusting sex as a categorical variable by changing the 0 and 1 coding
metadata_filtered$Sex_Male1_Female2 <- factor(metadata_filtered$Sex_Male1_Female2,
                          levels = c("1", "2"),
                          labels = c("Male", "Female"))
View(metadata_filtered)
metadata_filtered$Constipation_No0_Yes1 <- factor(metadata_filtered$Constipation_No0_Yes1,
                          levels = c(0, 1),
                          labels = c("No", "Yes"))

calling_maaslin3 = function(maaslin_input, metadata_df, output_dir, median_comparison_abundance, median_comparison_prevalence, formula=NULL, fixed_effects = NULL,
group_effects = NULL){
    x = maaslin3(input_data = maaslin_input, input_metadata = metadata_df,
        output = output_dir, formula = formula, fixed_effects = fixed_effects, group_effects = group_effects, normalization = normalization_to_use,
        transform = transformation_method, augment = TRUE, standardize = TRUE, max_significance = 0.05, evaluate_only = "abundance",
        median_comparison_abundance = TRUE, median_comparison_prevalence = FALSE, max_pngs = 600, cores = 8, warn_prevalence = FALSE)
}

calling_maaslin2 = function(maaslin_input, metadata_df, output_dir, fixed_effects_to_use,
                            transformation_method, normalization_to_use,random_to_use = NULL, reference_value){

    x = Maaslin2(input_data = maaslin_input, input_metadata = metadata_df, output = output_dir,
        fixed_effects = fixed_effects_to_use, transform = transformation_method,
        normalization = normalization_to_use, random_effects = random_to_use,
        reference = reference_value, cores = 8, plot_scatter = FALSE, min_prevalence = 0.05,
        standardize = FALSE)

    df_in = read.csv(file.path(output_dir, "all_results.tsv"), header = TRUE, sep = "\t")

    df_in$NegLog10qval  = -log10(df_in$qval)
    df_in$Specific_taxa = sub(".*s__", "", df_in$feature)

    df_in$qval_within_level = ave(df_in$pval, df_in$metadata, df_in$value,
                                  FUN = function(p) p.adjust(p, method = "BH"))
    
    df_in$NegLog10qval_within_level = -log10(df_in$qval_within_level)
    df_in$Specific_taxa_within_level = sub(".*s__", "", df_in$feature)
    df_in$Significant <- ifelse(
        df_in$qval_within_level >= 0.05, "No evidence of enrichment",
        ifelse(abs(df_in$coef) >= 1 & df_in$qval_within_level < 0.05, "Enriched",
        ifelse(abs(df_in$coef) > 0.585 & df_in$qval_within_level < 0.05, "Potentially enriched",
        "Lolz")))

    write.csv(df_in, file.path(output_dir, "all_results_with_neglog10qval.csv"), row.names = FALSE)

}

# Subsetting metadata


#metadata_male_only = metadata_df[metadata_df$Sex_Male1_Female2 == "Male", ]
#metadata_female_only = metadata_df[metadata_df$Sex_Male1_Female2 == "Female", ]
#metadata_parkinsons_only = metadata_df[metadata_df$Group == "PD", ]
#metadata_healthy_only = metadata_df[metadata_df$Group == "HC", ]
#metadata_parkinsons_only_hfiltered = metadata_parkinsons_only[metadata_parkinsons_only$H_and_Y_stage != "", ]
#metadata_parkinsons_only_hfiltered = metadata_parkinsons_only_hfiltered[!is.na(metadata_parkinsons_only_hfiltered$H_and_Y_stage), ]
#metadata_parkinsons_only_hfiltered$H_and_Y_stage
#metadata_parkinsons_levodopa = metadata_parkinsons_only[metadata_parkinsons_only$Levodopa_no_0_yes_1 == "Yes", ]

maaslin2_input_f_species = filtering_dataset_by_taxa_rank(maaslin_input_filtered, ";s__")

View(maaslin2_input_f_species)
View(metadata_filtered)

### MAASLN 2 ###
# PD vs HC
pd_results_2_with_sex = calling_maaslin2(maaslin2_input_f_species, metadata_filtered, output_dir = paste0(path,"/analysis_may2026/results1/maaslin2/PD_vs_HC_accounting_sex_prev05"), fixed_effects_to_use = c("Group", "Sex_Male1_Female2"),
                    transformation_method = transformation_method, normalization_to_use = normalization_to_use,
                    reference_value = "PD")

pd_results_3_with_sex_and_age = calling_maaslin2(maaslin2_input_f_species, metadata_filtered, output_dir = paste0(path,"/analysis_may2026/results1/maaslin2/PD_vs_HC_accounting_sex_and_age"), fixed_effects_to_use = c("Group", "Sex_Male1_Female2", "Age"),
                    transformation_method = transformation_method, normalization_to_use = normalization_to_use,
                    reference_value = "PD")     

pd_results_4_with_sex_bmi_age = calling_maaslin2(maaslin_input_f_species, metadata_filtered, output_dir = paste0(path,"/analysis_may2026/maaslin2/PD_vs_HC_accounting_sex_bmi_age"), fixed_effects_to_use = c("Group", "Sex_Male1_Female2", "BMI", "Age"),
                    transformation_method = transformation_method, normalization_to_use = normalization_to_use,
                    reference_value = "PD")

pd_only_by_sex = calling_maaslin2(maaslin_input_f_species, metadata_filtered[metadata_filtered$Group == "PD", ], output_dir = paste0(path,"/analysis_may2026/maaslin2/PD_only"), fixed_effects_to_use = "Group_and_sex",
                    transformation_method = transformation_method, normalization_to_use = normalization_to_use,
                    reference_value = "PD_Male") #PD male vs PD female

sex_differences_pd_only = calling_maaslin2(maaslin_input_f_species, metadata_filtered[metadata_filtered$Sex_Male1_Female2 == "Male", ], output_dir = paste0(path,"/analysis_may2026/maaslin2/sex_differences_only"), fixed_effects_to_use = "Group_and_sex",
                    transformation_method = transformation_method, normalization_to_use = normalization_to_use,
                    reference_value = "PD_Male") #PD male vs HC male

pd_results_5_with_sex_and_laxative = calling_maaslin2(maaslin_input_f_species, metadata_filtered, output_dir = paste0(path,"/analysis_may2026/maaslin2/PD_vs_HC_accounting_for_laxative_and_sex"), fixed_effects_to_use = c("Group", "Sex_Male1_Female2", "Laxatives"),
                    transformation_method = transformation_method, normalization_to_use = normalization_to_use,
                    reference_value = "PD")


### MAASLIN2 STRATIFIED ###

pd_results_stratified = calling_maaslin2(maaslin_input_f_species, metadata_filtered[metadata_filtered$Sex_Male1_Female2 == "Male", ], output_dir = paste0(path,"/analysis_may2026/maaslin2/PD_vs_HC_in_males_prev10"), fixed_effects_to_use = "Group",
                    transformation_method = transformation_method, normalization_to_use = normalization_to_use,
                    reference_value = "PD_Male")

female_results_stratified = calling_maaslin2(maaslin_input_f_species, metadata_filtered[metadata_filtered$Sex_Male1_Female2 == "Female", ], output_dir = paste0(path,"/analysis_may2026/maaslin2/PD_vs_HC_in_females_prev10"), fixed_effects_to_use = "Group",
                    transformation_method = transformation_method, normalization_to_use = normalization_to_use,
                    reference_value = "PD_Female")


### TESTING UNIVARIATE MODELS FIRST ###

pd_results_2 = calling_maaslin2(maaslin_input_filtered, metadata_filtered, output_dir = paste0(path,"/analysis_may2026/maaslin2/PD_vs_HC_prev10"), fixed_effects_to_use = "Group",
                    transformation_method = transformation_method, normalization_to_use = normalization_to_use,
                    reference_value = "PD")

age_differences = calling_maaslin2(maaslin_input_filtered, metadata_filtered, output_dir = paste0(path,"/analysis_may2026/maaslin2/age_differences"), fixed_effects_to_use = "Age",
                    transformation_method = transformation_method, normalization_to_use = normalization_to_use,
                    reference_value = NULL) #Age as a continuous variable

sex_differences = calling_maaslin2(maaslin_input_f_species, metadata_df, output_dir = paste0(path,"/analysis_may2026/maaslin2/sex_differences"), fixed_effects_to_use = "Sex_Male1_Female2",
                    transformation_method = transformation_method, normalization_to_use = normalization_to_use,
                    reference_value = "Male") #Male vs Female

constipation_differences = calling_maaslin2(maaslin_input_f_species, metadata_df, output_dir = paste0(path,"/analysis_may2026/maaslin2/constipation_differences"), fixed_effects_to_use = "Constipation_No0_Yes1",
                    transformation_method = transformation_method, normalization_to_use = normalization_to_use,
                    reference_value = "Yes") #No constipation vs Yes constipation

bmi_differences = calling_maaslin2(maaslin_input_f_species, metadata_df, output_dir = paste0(path,"/analysis_may2026/maaslin2/bmi_differences"), fixed_effects_to_use = "BMI",
                    transformation_method = transformation_method, normalization_to_use = normalization_to_use,
                    reference_value = NULL) #BMI as a continuous variable

levodopa_differences = calling_maaslin2(maaslin_input_f_species, metadata_parkinsons_only, output_dir = paste0(path,"/analysis_may2026/maaslin2/levodopa_differences"), fixed_effects_to_use = "Levodopa_dosage",
                    transformation_method = transformation_method, normalization_to_use = normalization_to_use,
                    reference_value = NULL) #Levodopa yes vs no in Parkinson's patients only

ledd_differences = calling_maaslin2(maaslin_input_f_species, metadata_parkinsons_only, output_dir = paste0(path,"/analysis_may2026/maaslin2/ledd_differences"), fixed_effects_to_use = "LEDD",
                    transformation_method = transformation_method, normalization_to_use = normalization_to_use,
                    reference_value = NULL) #LEDD as a continuous variable in Parkinson's patients only

updrs_total_differences = calling_maaslin2(maaslin_input_f_species, metadata_parkinsons_only, output_dir = paste0(path,"/analysis_may2026/maaslin2/updrs_total_differences"), fixed_effects_to_use = "UPDRS.TOTAL",
                    transformation_method = transformation_method, normalization_to_use = normalization_to_use,
                    reference_value = NULL) #UPDRS total as a continuous variable in Parkinson's patients only

antidepressants_differences = calling_maaslin2(maaslin_input_f_species, metadata_df, output_dir = paste0(path,"/analysis_may2026/maaslin2/antidepressants_differences"), fixed_effects_to_use = "Antidepressants_coded",
                    transformation_method = transformation_method, normalization_to_use = normalization_to_use,
                    reference_value = NULL) #Antidepressants yes vs no in Parkinson's patients only

laxatives_differences = calling_maaslin2(maaslin_input_f_species, metadata_df, output_dir = paste0(path,"/analysis_may2026/maaslin2/laxatives_differences"), fixed_effects_to_use = "Laxatives",
                    transformation_method = transformation_method, normalization_to_use = normalization_to_use,
                    reference_value = NULL) #Laxatives yes vs no in Parkinson's patients only

benzos_differences = calling_maaslin2(maaslin_input_f_species, metadata_df, output_dir = paste0(path,"/analysis_may2026/maaslin2/benzos_differences"), fixed_effects_to_use = "Benzos_coded",
                    transformation_method = transformation_method, normalization_to_use = normalization_to_use,
                    reference_value = NULL) #Benzos yes vs no in Parkinson's patients only

### MAASLIN3 ###
















# PD vs HC - Abundance Yes, Prevalence No
pd_results = calling_maaslin3(maaslin_input = maaslin_input, metadata_df = metadata_df, output_dir = paste0(path,"/analysis_may2026/maaslin3/PD_vs_HC"), 
                                formula = "~ Group ", median_comparison_abundance = TRUE, 
                                median_comparison_prevalence = FALSE)

dim(maaslin_input)

pd_results = calling_maaslin3(maaslin_input = maaslin_input, metadata_df = metadata_df, output_dir = paste0(path,"/analysis_may2026/maaslin3/PD_vs_HC_Group_and_Sex"),
                                fixed_effects = c("Group"), median_comparison_abundance = TRUE,
                                median_comparison_prevalence = FALSE)
                                

#formula = "~ Group + Sex_Male1_Female2 + Age + BMI + Constipation_No0_Yes1",
# Parkinson's vs Healthy only Males
pd_male_only_results = calling_maaslin2(maaslin_input_f_species, metadata_male_only, output_file = paste0(path,"/analysis_may2026/maaslin2/PD_vs_HC_in_males"), fixed_effects_to_use = "Males_only",
                    transformation_method = transformation_method, normalization_to_use = normalization_to_use, random_to_use = random_to_use, 
                    reference_value = "HC_Male")

# Parkinson's vs Healthy only Females
pd_female_only_results = calling_maaslin2(maaslin_input_f_species, metadata_female_only, output_file = paste0(path,"/analysis_may2026/maaslin2/PD_vs_HC_in_females"), fixed_effects_to_use = "Females_only",
                    transformation_method = transformation_method, normalization_to_use = normalization_to_use, random_to_use = random_to_use, 
                    reference_value = "HC_Female")

# Healthy ctrls males vs females
pd_sex_results = calling_maaslin2(maaslin_input_f_species, metadata_healthy_only, output_file = paste0(path,"/analysis_may2026/maaslin2/HC_males_vs_females"), fixed_effects_to_use = "Sex_Male1_Female2",
                    transformation_method = transformation_method, normalization_to_use = normalization_to_use, random_to_use = random_to_use, 
                    reference_value = "Male")

# Parkinson's males vs females
hc_sex_results = calling_maaslin2(maaslin_input_f_species, metadata_parkinsons_only, output_file = paste0(path,"/analysis_may2026/maaslin2/PD_males_vs_females"), fixed_effects_to_use = "Sex_Male1_Female2",
                    transformation_method = transformation_method, normalization_to_use = normalization_to_use, random_to_use = random_to_use, 
                    reference_value = "Male")









pd_parkinsons_only_results = calling_maaslin2(maaslin_input_f_species, metadata_parkinsons_only, output_file = paste0(path,"/analysis_may2026/maaslin2/pd_only"), fixed_effects_to_use = "Sex_male1_female2",
                    transformation_method = transformation_method, normalization_to_use = normalization_to_use, random_to_use = random_to_use, 
                    reference_value = "Female")


hc_parkinsons_only_results = calling_maaslin2(maaslin_input_f_species, metadata_healthy_only, output_file = "maaslin2/hc_only", fixed_effects_to_use = "Sex_male1_female2",
                    transformation_method = transformation_method, normalization_to_use = normalization_to_use, random_to_use = random_to_use, 
                    reference_value = "Female")

parkinsons_constipation = calling_maaslin2(maaslin_input_f_species, metadata_parkinsons_only, output_file = "maaslin2/pd_constipation", fixed_effects_to_use = "Constipation_no_0_yes_1",
                    transformation_method = transformation_method, normalization_to_use = normalization_to_use, random_to_use = random_to_use, 
                    reference_value = "Yes")

h_and_y_results = calling_maaslin2(maaslin_input_f_species, metadata_parkinsons_only_hfiltered, output_file = "maaslin2/pd_h_and_y", fixed_effects_to_use = "H_and_Y_stage",
                    transformation_method = transformation_method, normalization_to_use = normalization_to_use, random_to_use = random_to_use, 
                    reference_value = "H_and_Y_stage,4")

updrs_total_results = calling_maaslin2(maaslin_input_f_species, metadata_parkinsons_only, output_file = "maaslin2/pd_updrs_total", fixed_effects_to_use = "UPDRS_Total",
                    transformation_method = transformation_method, normalization_to_use = normalization_to_use, random_to_use = random_to_use, 
                    reference_value = NULL)

humann_pd_results = calling_maaslin2(humann_data_f_species, metadata_df, output_file = "maaslin2/humann_pd_vs_no_pd", fixed_effects_to_use = "Group",     
                    transformation_method = transformation_method, normalization_to_use = normalization_to_use, random_to_use = random_to_use, 
                    reference_value = reference_value)

human_pd_f_results = calling_maaslin2(humann_data_f_species, metadata_df, output_file = "maaslin2/humann_pd_and_no_pd_f", fixed_effects_to_use = "Group",
                    transformation_method = transformation_method, normalization_to_use = normalization_to_use, random_to_use = random_to_use, 
                    reference_value = reference_value)

human_pd_semi_f_results = calling_maaslin2(humann_data_semi_f_species, metadata_df, output_file = "maaslin2/humann_pd_and_no_pd_semi_f", fixed_effects_to_use = "Group",
                    transformation_method = transformation_method, normalization_to_use = normalization_to_use, random_to_use = random_to_use, 
                    reference_value = reference_value)
