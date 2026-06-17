# Daniel Castaneda Mogollon, PhD
# May 26th, 2026
# This script filters the Metaphlan4 abundances to keep only whatever samples/features the user wants based
# on a range of options

library(dplyr)

input_dir = "/Users/danielcm/Desktop/MIGUT/analysis_may2026/merged_tables_abundances"
output_dir = "/Users/danielcm/Desktop/MIGUT/analysis_may2026/filtered_metaphlan4_abundances"

samples_to_remove = c("PD_199P_profile","PD_107P_profile","PD_81P_profile","PD_218S_profile","PD_36P_profile","PD_192S_profile","PD_26P_profile","PD_112S_profile",
                        "PD_212S_profile","PD_216S_profile","PD_191P_profile","PD_211P_profile","PD_83P_profile","PD_217P_profile","PD_182S_profile","PD_161S_profile",
                        "PD_50P_profile","PD_28P_profile")

setwd(input_dir)

print(paste0("A total of ",length(samples_to_remove)," samples will be removed from the data frame."))
df_in = read.csv("merged_absolute_abundances_table_GTDB.txt", header = TRUE, sep = "\t", row.names = 1)

print(paste0("The original data frame has ",nrow(df_in)," rows and ",ncol(df_in)," columns."))

df_filt = df_in %>%
    slice(2:nrow(df_in)) %>% #This removes the 'UNCLASSIFIED' row from the data frame
    select(-all_of(samples_to_remove))

print(paste0("After filtering by a cutoff of 900,000 reads, the data frame has ",nrow(df_filt)," rows and ",ncol(df_filt)," columns."))

df_filt_species = df_filt %>%
    filter(grepl("s__", rownames(df_filt)))

print(paste0("After filtering to keep only species-level features, the data frame has ",nrow(df_filt_species)," rows and ",ncol(df_filt_species)," columns."))
print(paste0("A total of ",nrow(df_filt)-nrow(df_filt_species)," features were removed by filtering to keep only species-level features."))

proportion_of_zeros_to_remove = 1.0




df_filt2 = as.data.frame(df_filt_species)

df_filt2$zero_count = rowSums(df_filt2 == 0)
df_filt2$zero_percentage = df_filt2$zero_count/ncol(df_filt2)
df_filt2$zero_percentage

df_filt2 = df_filt2 %>%
    filter(zero_percentage <= proportion_of_zeros_to_remove)


print(paste0("After filtering to keep only features with less than ", proportion_of_zeros_to_remove*100, "% zeros, the data frame has ",nrow(df_filt2)," rows and ",ncol(df_filt2)," columns."))
print(paste0("A total of ",nrow(df_filt_species)-nrow(df_filt2)," features were removed by filtering to keep only features with less than ", proportion_of_zeros_to_remove*100, "% zeros."))

write.table(df_filt2[,1:(ncol(df_filt2)-2)], file.path(output_dir, "merged_absolute_abundances_table_GTDB_filtered.txt"), sep = "\t", row.names = TRUE, quote = FALSE)
