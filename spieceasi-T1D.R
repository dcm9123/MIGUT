# Daniel Castaneda Mogollon, PhD
# May 27th, 2026
# This script will run the network co-association between taxa by using spiceasi

library("SpiecEasi")
library("phyloseq")
library("Matrix")
library("igraph")
library("qgraph")
library("stringr")
#library("ComplexHeatmap")
library("pheatmap")


set.seed(23)

path = "/Users/danielcm/Desktop/SickKids/"
setwd(path)

#f_input = read.table(paste0(path,"/analysis_may2026/filtered_metaphlan4_abundances/merged_absolute_abundances_table_GTDB_filtered.txt"), header = TRUE, sep = "\t", row.names = 1)
# Generating the ps object

consortium = "ns1"

raw_input = read.table(paste0(path,"PICRUSt2.6/",consortium,"_input/",consortium,"_filtered_asvs_count600_len400_prev20_f_sra.tsv"), header = TRUE, sep = "\t", row.names = 1)
metadata_df = read.csv(paste0(path,"/Metadata/Danska_diabetes_metadata364_20260409.txt"), header = TRUE, sep = "\t")

metadata_df = as.data.frame(metadata_df)
colnames(metadata_df)

otu_table = as.matrix(raw_input)
colnames(otu_table) = lapply(colnames(otu_table), function(name) gsub("X","",name))


taxa_table = as.matrix(read.csv(paste0(path,"/Phyloseq2/taxa_ps_mice_samples.tsv"), header = TRUE, sep = "\t"))
rownames(taxa_table) = taxa_table[,1]
taxa_table = taxa_table[,-1]
taxa_table = as.matrix(taxa_table)
taxa_df = as.data.frame(taxa_table)
taxa_df$Genus_and_species = paste(taxa_df$genus_final, taxa_df$species_final, sep = "_")
taxa_matrix = as.matrix(taxa_df)


sample_data = metadata_df[,c("Id","Sex","Consortium","Control","Inoculum","Timepoint","Merged.Weeks","Sex.and.Consortium","Sex.and.Timepoint","Week.and.Consortium")]
rownames(sample_data) = sample_data$Id
sample_data$Id = lapply(sample_data$Id, function(id) gsub("Plate\\d_","",id))
sample_data$Id = lapply(sample_data$Id, function(id) gsub("_S\\d+_L001","",id))
sample_data = subset(sample_data, Consortium == toupper(consortium))
sample_data = sample_data[!(sample_data$Id %in% c("Positive_Control","Negative_Control", paste0(toupper(consortium), "_Inoculum"))),]
rownames(sample_data) = sample_data$Id
rownames(sample_data)

ps = phyloseq(otu_table(otu_table, taxa_are_rows = TRUE), tax_table(taxa_matrix), sample_data(sample_data))
#ps = merge_taxa(ps, "species_final")


#SPIEC-EASI
# This software works by inferring microbial association networks from composituional microbiome data
# It does not simply generate a matrix of Pearson correlation coefficients, but instead it tries to determine or 
# generate a newtork of conditional dependencies between taxa. Standard correlations using relative abundances may provide 
# false associations simply by having different absolute abundances.

# It asks if taxa X and Y are associated after taking into account the rest of the taxa in the system/dataset.
# Rows = samples, columns = taxa. Use raw or filtered counts.




# Filters applied: Keep only species level taxa, remove those present less than 5% of samples, remove samples with no abundances, 
# removes taxa with a mean relative abundance of less than 0.001%.

applying_phyloseq_filters = function(ps_object, sample_prevalence_filter, taxa_mean_prevalence, category_to_analyze, category_column) {
    if(!(category_column %in% colnames(sample_data(ps_object)))) {
        stop(paste("Metadata's column specified:", category_column, "is not present in the metadata file"))
    }

    metadata = as.data.frame(sample_data(ps_object))
    samples_to_keep = !is.na(metadata[[category_column]]) & metadata[[category_column]] == category_to_analyze
    names(samples_to_keep) = rownames(metadata)

    ps = prune_samples(samples_to_keep, ps_object) # Subset on samples of interest from the metadata (i.e. PD only, HC_male)
    print(ps)
    ps_f = subset_taxa(ps, !is.na(Species) & Species != "") # Removes everything that is not speciated
    ps_f = prune_samples(sample_sums(ps_f) > 0, ps_f) # Removes samples with no abundances (0)
    print(ps_f)    
    ps_f = filter_taxa(ps_f, function(x) sum(x >0) > (sample_prevalence_filter*nsamples(ps_f)), TRUE) # Removes taxa present in less than 5% of samples or whatever sample_pevalence_filter was set to.
    
    ps_f_relative_a = transform_sample_counts(ps_f, function(x) x/sum(x)) # Transforms to relative abundances
    ps_f_relative_a = filter_taxa(ps_f_relative_a, function(x) mean(x) > taxa_mean_prevalence, TRUE) # Removes taxa with a mean relative abundance of less than 0.001% or whatever the user set it to 
    
    ps_f_final = prune_taxa(taxa_names(ps_f_relative_a), ps_f) # Prune the original ps object to keep the same taxa as the filtered relative abundance object, but with the original counts. This is the one we will use for SPIEC-EASI
    ps_f_final = prune_taxa(taxa_sums(ps_f_final) > 0, ps_f_final) # Remove taxa with no counts after all the filtering. This is just in case some taxa had a mean relative abundance of less than 0.001% but still had some counts in some samples, and thus were not removed by the previous filter.
    ps_f_final = prune_samples(sample_sums(ps_f_final) > 0, ps_f_final) # Remove samples with no counts after all the filtering. This is just in case some samples had some taxa with a mean relative abundance of less than 0.001% but still had some counts in some taxa, and thus were not removed by the previous filter. 

    print(paste("Analyzing the ps object with samples from the category:", category_to_analyze))
    print(paste("The number of samples before filtering was:", nsamples(ps_object)))
    print(paste("The number of samples after filtering was:", nsamples(ps_f_final)))
    print(paste("The number of taxa before filtering was:", ntaxa(ps_object)))
    print(paste("The number of taxa after filtering was:", ntaxa(ps_f_final)))
    return(ps_f_final)
}


# Log-ratio transforms the data, selects the best model with 'pulsar using stars', and fits the final estimate with the method of choice:
# Both 'mb' and 'glasso' are methods to infer conditional dependency between taxa. The first is neighborhood selection, and it 
# uses the lasso regression ().

# Step 1: Data input as columns taxa, and rows samples. 
# Step 2: Filter data if needed (i.e. sample prevalence, or rare taxa).
# Step 3: CLR transformation of each taxon. This gives us a CLR matrix. 
# Step 4: Use LASSO regression to minimize the RSS + lambda(sum(betas)), where each beta is a coefficient for each taxa. We use one taxon at a time
# as the response variable (the CLR value) and the rest of the taxa*beta as the predictors. We iterate each beta to minimize the RSS + lambda(sum(betas)) as much
# as we can while keeping the simplest model possible. All samples are used at once to iterate over all betas for all the response variables. Once the best model is found,
# then repeat with a new Lambda value. A coefficient is kept if the reduction of the RSS is worth the increase in penalty.
# SPIEC-EASI chooses the sparsest network that remains stable across many subsampling with many lambda values. If it sees an edge is kept, then we know it holds.
# The tradeoff between RSS reduction and model simplicity is determined by the lambda parameter. The higher the lambda, the more we penalize the model for having more edges, and thus the sparser the network.

# Neighborhood selection (mb) first, glasso second, and sparcc third

network_generation = function(ps_object, method_to_use, lambda_min_ratio, n_lambda, pulsar_rep_num, pulsar_ncores, sparcc_method){
    asv_table = t(otu_table(ps_object))
    spic_network = spiec.easi(asv_table, method = method_to_use, lambda.min.ratio = lambda_min_ratio, nlambda = n_lambda, pulsar.params = list(rep.num = pulsar_rep_num, ncores = pulsar_ncores))
    
    sparcc_network = sparcc(asv_table)
    sparcc_graph = abs(sparcc_network$Cor) >= 0.5
    diag(sparcc_graph) = 0
    sparcc_graph = Matrix(sparcc_graph, sparse = TRUE)
    if(sparcc_method == FALSE){
        return(spic_network)
    }
    else{
        return(sparcc_graph)
    }
}

  # sparcc = a correlation method that does not try to remove indirect associations, and it estimates it from 
  # compositional data, not dependency-based correlations from sparcity data.
  # assumptions: the true underlying abundances are not compositional, and most taxa pairs are not strongly correlated.

make_tax_rank_colors = function(ps_object, tax_rank){
    taxa_df = as.data.frame(tax_table(ps_object))
    tax_rank_prefix = paste0(tolower(substr(tax_rank, 1, 1)), "__")
    tax_rank_values = sub(paste0("^", tax_rank_prefix), "", taxa_df[, tax_rank])
    tax_rank_values = sort(unique(tax_rank_values[!is.na(tax_rank_values) & tax_rank_values != ""]))
    setNames(
        grDevices::hcl.colors(length(tax_rank_values), palette = "Dark 3"),
        tax_rank_values)
}


plotting_networks = function(network_to_plot, spic_network, category, ps_object, tax_rank, plot_file_name, 
                        tax_rank_colors = NULL, multicommunity, labels, color_for_vertex, shared_list, hei=8, wid=14){
    taxa_matrix = tax_table(ps_object)
    taxa_df = as.data.frame(taxa_matrix)
    species_names = taxa_df$Genus_and_species
    
    matrix_to_work = as.matrix(otu_table(ps_object))
    matrix_to_work = as.data.frame(matrix_to_work)
    matrix_to_work$non_zero_count = rowSums(matrix_to_work>0)/ncol(matrix_to_work)
    non_zero_count = matrix_to_work$non_zero_count
    names(non_zero_count) = rownames(matrix_to_work)

    if(spic_network == "mb" | spic_network == "glasso"){
        adj_matrix = as.matrix(getRefit(network_to_plot))
        igraph_object = adj2igraph(adj_matrix, rmEmptyNodes = FALSE)
    }
    else{
        igraph_object = adj2igraph(network_to_plot)
    }
    
    taxa_df = as.data.frame(tax_table(ps_object))
    edge_list = ends(igraph_object, E(igraph_object), names = TRUE)

    if(spic_network == "mb"){
        beta_coefficients = getOptBeta(network_to_plot)
        beta_symmetric = symBeta(beta_coefficients, mode = "maxabs")
        E(igraph_object)$weight = beta_symmetric[cbind(edge_list[,1], edge_list[,2])]
    }
    else{
        cov_matrix = cov2cor(getOptCov(network_to_plot))
        #cov_matrix = getOptCov(network_to_plot)
        weights = adj_matrix * cov_matrix
        E(igraph_object)$weight = weights[cbind(edge_list[,1], edge_list[,2])]
    }

    E(igraph_object)$abs_weight <- abs(E(igraph_object)$weight)
    E(igraph_object)$sign <- ifelse(E(igraph_object)$weight > 0, "positive", "negative")

    tax_rank_prefix = paste0(tolower(substr(tax_rank, 1, 1)), "__")
    V(igraph_object)$alpha = non_zero_count[V(igraph_object)$name]
    V(igraph_object)$tax_rank = taxa_df[V(igraph_object)$name, tax_rank]
    V(igraph_object)$tax_rank = sub(paste0("^", tax_rank_prefix), "", V(igraph_object)$tax_rank)
    V(igraph_object)$label = species_names
    V(igraph_object)$shape = "circle"
    V(igraph_object)$color = "black"

     # Store the original OTU ID key in an attribute so we can find it after vertex deletions
    V(igraph_object)$otu_id = V(igraph_object)$name
    V(igraph_object)$name = species_names

    cutoff <- quantile(E(igraph_object)$abs_weight, 0, na.rm = TRUE)

    igraph_filt <- delete_edges(
        igraph_object,
        E(igraph_object)[abs_weight < cutoff]
    )

    igraph_filt <- delete_vertices(
        igraph_filt,
        V(igraph_filt)[degree(igraph_filt) < 0]
    )

    E(igraph_filt)$layout_weight <- E(igraph_filt)$abs_weight
    if(is.null(tax_rank_colors)){
        tax_rank_colors = make_tax_rank_colors(ps_object, tax_rank)
    }
    
    if(color_for_vertex == TRUE){
        classes = sort(unique(V(igraph_filt)$tax_rank))
        class_colors = tax_rank_colors[classes]
        V(igraph_filt)$color = class_colors[V(igraph_filt)$tax_rank]
    }
    else{
        if(is.null(shared_list) == FALSE){
            for(item in shared_list){
                V(igraph_filt)[label == item]$color = "blue"
                V(igraph_filt)[label == item]$shape = "square"
            }
        }
    }
    set.seed(123)

    layout_fr <- layout_with_fr(
        igraph_filt,
        weights = E(igraph_filt)$layout_weight
    )

    if(labels == TRUE){
        v_label = V(igraph_filt)$label
    } else {
        v_label = NA
    }

    vertex_size = log10(apply(otu_table(ps_object), 1, median, na.rm = TRUE) + 1)
    vertex_size = scales::rescale(vertex_size, to = c(1, 30))

    print(min(vertex_size))
    print(max(vertex_size))

    # 1. Compute your 0.2 to 0.8 transparency limits
    raw_alpha_vector <- as.numeric(V(igraph_filt)$alpha)
    raw_alpha_vector[is.na(raw_alpha_vector)] <- 0  # Replace NAs with zero
    final_alpha <- as.numeric(scales::rescale(raw_alpha_vector, to = c(0, 1)))

    # 2. Extract original colors and guarantee zero NA values
    base_colors <- as.character(V(igraph_filt)$color)
    base_colors[is.na(base_colors) | base_colors == "" | base_colors == "NA"] <- "black"

    # 3. Apply alpha mapping using the reliable scales::alpha function
    # This completely replaces adjustcolor() and avoids the d == c(4L, 4L) check entirely
    final_vertex_colors <- scales::alpha(base_colors, alpha = final_alpha)
    final_frame_colors  <- scales::alpha("black", alpha = final_alpha)

    # 4. Render the PNG network graph cleanly
    png(plot_file_name, height = hei, width = wid, units = "in", res = 600)

    plot(igraph_filt,
         vertex.size = vertex_size,
         vertex.label = v_label,
         main = category,
         edge.width = scales::rescale(E(igraph_filt)$abs_weight, to = c(1, 3)),
         edge.color = ifelse(E(igraph_filt)$sign == "positive", "black", "red"),

         # Pass the safe hex strings generated by scales::alpha
         vertex.color = final_vertex_colors,
         vertex.frame.color = final_frame_colors,
         layout = layout_fr)

    if(color_for_vertex == TRUE){
        legend(
            "topleft",
            legend = names(class_colors),
            col = class_colors,
            pch = 19,
            pt.cex = 1.5,
            bty = "n",
            cex = 0.8
        )
    }
    dev.off()
    return(igraph_object)
}



#Generating best model and plotting_networks, SPIC method
ps_w9w10 = prune_samples(sample_data(ps)$Timepoint == "week_9" | sample_data(ps)$Timepoint == "week_10", ps)
#40 taxa, 52 samples
ps_w9w10 = prune_taxa(taxa_sums(otu_table(ps_w9w10)) > 0, ps_w9w10)
#40 taxa, 52 samples for ns1
#39 taxa, 47 samples for s2

ps_w5 = prune_samples(sample_data(ps)$Timepoint == "week_5", ps)    
ps_w5
ps_w9w10
#40 taxa, 41 samples for ps
#39 taxa, 50 samples for s2


###         MANUAL ADJUSTMENT                 ###
# Run only for manual curation of known species #
#################################################
# S2 NEEDS THIS!
print("Running manual adjustment of species, week 5")
tax_df = as.data.frame(tax_table(ps_w5))
tax_df$Genus_and_species[tax_df$genus_final == "Parabacteroides"] <- "Parabacteroides_distasonis"
tax_df$species_final[tax_df$Genus_and_species == "Parabacteroides_distasonis"] <- "distasonis"
tax_df$Genus_and_species[tax_df$genus_final == "Hungatella"] <- "Hungatella_effluvi"
tax_df$species_final[tax_df$Genus_and_species == "Hungatella_effluvi"] <- "effluvi"
tax_table(ps_w5) = as.matrix(tax_df)

tax_df = as.data.frame(tax_table(ps_w9w10))
tax_df$Genus_and_species[tax_df$genus_final == "Parabacteroides"] <- "Parabacteroides_distasonis"
tax_df$species_final[tax_df$Genus_and_species == "Parabacteroides_distasonis"] <- "distasonis"
tax_df$Genus_and_species[tax_df$genus_final == "Hungatella"] <- "Hungatella_effluvi"
tax_df$species_final[tax_df$Genus_and_species == "Hungatella_effluvi"] <- "effluvi"
tax_table(ps_w9w10) = as.matrix(tax_df)

################################################

ps_w9w10 = tax_glom(ps_w9w10, taxrank = "Genus_and_species")
ps_w9w10 #20 taxa, 52 samples (ns1), 19 taxa, 47 samples, s2
ps_w5 = tax_glom(ps_w5, taxrank = "Genus_and_species")
ps_w5 #20 taxa, 41 samples (ns1), 19 taxa, 50 samples, s2

#View(tax_table(ps_w9w10))
otu_table_ps = as.data.frame(otu_table(ps_w5))
otu_table_ps$rowmedian = (apply(otu_table_ps, 1, median, na.rm = TRUE))
#View(otu_table_ps)



ps_w9w10_network_mb = network_generation(ps_w9w10, method_to_use = "mb", lambda_min_ratio = 1e-2, n_lambda = 30, pulsar_rep_num = 50, pulsar_ncores = 8, sparcc_method = FALSE)
ps_w5_network_mb = network_generation(ps_w5, method_to_use = "mb", lambda_min_ratio = 1e-2, n_lambda = 30, pulsar_rep_num = 50, pulsar_ncores = 8, sparcc_method = FALSE)
ps_w9w10_network_lasso = network_generation(ps_w9w10, method_to_use = "glasso", lambda_min_ratio = 1e-2, n_lambda = 30, pulsar_rep_num = 50, pulsar_ncores = 8, sparcc_method = FALSE)
ps_w5_network_lasso = network_generation(ps_w5, method_to_use = "glasso", lambda_min_ratio = 1e-2, n_lambda = 30, pulsar_rep_num = 50, pulsar_ncores = 8, sparcc_method = FALSE)


# Plotting plotting_networks
species_colors = make_tax_rank_colors(ps_w9w10, "Genus_and_species")
shared_microbes = c("Akkermansia_muciniphila", "Alistipes_finegoldii", "Alistipes_onderdonkii",
                    "Escherichia_coli", "Bifidobacterium_longum", "Eggerthella_lenta", "Eisenbergiella_tayi",
                    "Enterocloster_bolteae","Flavonifractor_plautii","Hungatella_effluvi","Parabacteroides_distasonis",
                    "Phocaeicola_vulgatus","Staphylococcus_hominis")
length(shared_microbes)



### WEEK 9 AND WEEK 10 FIRST ###
adj_matrix_w9w10 = as_adjacency_matrix(plotting_networks(ps_w9w10_network_mb, spic_network = "mb", category = paste0(toupper(consortium), " w9w10"), ps_object = ps_w9w10, 
                  tax_rank = "Genus_and_species", plot_file_name = paste0(consortium, "_w9w10_mb_network_no_labels.png"),
                  tax_rank_colors = NULL, labels = FALSE, color_for_vertex = FALSE, shared_list = shared_microbes))

plotting_networks(ps_w9w10_network_mb, spic_network = "mb", category = paste0(toupper(consortium), " w9w10"), ps_object = ps_w9w10, 
                  tax_rank = "Genus_and_species", plot_file_name = paste0(consortium, "_w9w10_mb_network_labels.png"),
                  tax_rank_colors = NULL, labels = TRUE, color_for_vertex = FALSE, shared_list = shared_microbes, hei = 4, wid = 7)

plotting_networks(ps_w9w10_network_mb, spic_network = "mb", category = paste0(toupper(consortium), " w9w10"), ps_object = ps_w9w10, 
                  tax_rank = "Genus_and_species", plot_file_name = paste0(consortium, "mini_w9w10_mb_network_no_labels.png"),
                  tax_rank_colors = NULL, labels = FALSE, color_for_vertex = FALSE, shared_list = shared_microbes, hei = 4, wid = 7)

glasso_adj_matrix_w9w10 = as_adjacency_matrix(plotting_networks(ps_w9w10_network_lasso, spic_network = "glasso", category = paste0(toupper(consortium), " w9w10"), ps_object = ps_w9w10, 
                  tax_rank = "Genus_and_species", plot_file_name = paste0(consortium, "_w9w10_glasso_network_no_labels.png"),
                  tax_rank_colors = NULL, labels = FALSE, color_for_vertex = FALSE, shared_list = shared_microbes))

plotting_networks(ps_w9w10_network_lasso, spic_network = "glasso", category = paste0(toupper(consortium), " w9w10"), ps_object = ps_w9w10, 
                  tax_rank = "Genus_and_species", plot_file_name = paste0(consortium, "_w9w10_glasso_network_labels.png"),
                  tax_rank_colors = NULL, labels = TRUE, color_for_vertex = FALSE, shared_list = shared_microbes)

plotting_networks(ps_w9w10_network_lasso, spic_network = "glasso", category = paste0(toupper(consortium), " w9w10"), ps_object = ps_w9w10, 
                  tax_rank = "Genus_and_species", plot_file_name = paste0(consortium, "mini_w9w10_glasso_network_no_labels.png"),
                  tax_rank_colors = NULL, labels = FALSE, color_for_vertex = FALSE, shared_list = shared_microbes, hei = 4, wid = 7)

x = adj2igraph(glasso_adj_matrix_w9w10)
V(x)$edges



View(otu_table(ps_w5))

### WEEK 5 NOW ###
adj_matrix_w5 = as_adjacency_matrix(plotting_networks(ps_w5_network_mb, spic_network = "mb", category = paste0(toupper(consortium), " w5"), ps_object = ps_w5, 
                  tax_rank = "Genus_and_species", plot_file_name = paste0(consortium, "_w5_mb_network_no_labels.png"),
                  tax_rank_colors = NULL, labels = FALSE, color_for_vertex = FALSE, shared_list = shared_microbes))

plotting_networks(ps_w5_network_mb, spic_network = "mb", category = paste0(toupper(consortium), " w5"), ps_object = ps_w5, 
                  tax_rank = "Genus_and_species", plot_file_name = paste0(consortium, "_w5_mb_network_labels.png"),
                  tax_rank_colors = NULL, labels = TRUE, color_for_vertex = FALSE, shared_list = shared_microbes, hei = 4, wid = 7)

plotting_networks(ps_w5_network_mb, spic_network = "mb", category = paste0(toupper(consortium), " w5"), ps_object = ps_w5, 
                  tax_rank = "Genus_and_species", plot_file_name = paste0(consortium, "mini_w5_mb_network_no_labels.png"),
                  tax_rank_colors = NULL, labels = FALSE, color_for_vertex = FALSE, shared_list = shared_microbes, hei = 4, wid = 7)

glasso_adj_matrix_w5 = as_adjacency_matrix(plotting_networks(ps_w5_network_lasso, spic_network = "glasso", category = paste0(toupper(consortium), " w5"), ps_object = ps_w5, 
                  tax_rank = "Genus_and_species", plot_file_name = paste0(consortium, "_w5_glasso_network_no_labels.png"),
                  tax_rank_colors = NULL, labels = FALSE, color_for_vertex = FALSE, shared_list = shared_microbes))

plotting_networks(ps_w5_network_lasso, spic_network = "glasso", category = paste0(toupper(consortium), " w5"), ps_object = ps_w5, 
                  tax_rank = "Genus_and_species", plot_file_name = paste0(consortium, "_w5_glasso_network_labels.png"),
                  tax_rank_colors = NULL, labels = TRUE, color_for_vertex = FALSE, shared_list = shared_microbes)

plotting_networks(ps_w5_network_lasso, spic_network = "glasso", category = paste0(toupper(consortium), " w5"), ps_object = ps_w5, 
                  tax_rank = "Genus_and_species", plot_file_name = paste0(consortium, "mini_w5_glasso_network_no_labels.png"),
                  tax_rank_colors = NULL, labels = FALSE, color_for_vertex = FALSE, shared_list = shared_microbes, hei = 4, wid = 7)

log10(apply(otu_table(ps_w9w10), 1, median, na.rm = TRUE))
cbind(otu_table(ps_w9w10), tax_table(ps_w9w10))




# CONNECTIONS
psw9w10_network_matrix_lasso = as.matrix(getRefit(ps_w9w10_network_lasso))
total_connections_psw9w10_lasso = sum(psw9w10_network_matrix_lasso)/2
total_connections_psw9w10_lasso

psw5_network_matrix_lasso = as.matrix(getRefit(ps_w5_network_lasso))
total_connections_psw5_lasso = sum(psw5_network_matrix_lasso)/2
total_connections_psw5_lasso


# HARMONIC CENTRALITY / CLOSENNES ()
psw9w10_network_adjmatrix_lasso = adj2igraph(psw9w10_network_matrix_lasso, rmEmptyNodes = TRUE,vertex.attr = list(name = colnames(ps_w9w10)))
psw9w10_closeness_lasso = harmonic_centrality(psw9w10_network_adjmatrix_lasso, mode = "all", normalized = FALSE, cutoff = -1)
psw9w10_closeness_lasso

psw5_network_adjmatrix_lasso = adj2igraph(psw5_network_matrix_lasso, rmEmptyNodes = TRUE,vertex.attr = list(name = colnames(ps_w5)))
psw5_closeness_lasso = harmonic_centrality(psw5_network_adjmatrix_lasso, mode = "all", normalized = FALSE, cutoff = -1)
psw5_closeness_lasso

for(item in psw9w10_closeness_lasso){
  cat((log10(item)),"\n")
}

for(item in psw5_closeness_lasso){
  cat((log10(item)),"\n")
}



# DENSITY
density_psw9w10_lasso = edge_density(psw9w10_network_adjmatrix_lasso, loops = FALSE)
density_psw9w10_lasso

density_psw5_lasso = edge_density(psw5_network_adjmatrix_lasso, loops = FALSE)
density_psw5_lasso


# BETWEENNESS
ps_w9w10_betweenness_lasso = betweenness(psw9w10_network_adjmatrix_lasso,  directed = FALSE, normalized = TRUE, cutoff = -1)
ps_w5_betweenness_lasso = betweenness(psw5_network_adjmatrix_lasso,  directed = FALSE, normalized = TRUE, cutoff = -1)
ps_w9w10_betweenness_lasso
ps_w5_betweenness_lasso

ps_w9w10_betweenness_lasso_raw = betweenness(psw9w10_network_adjmatrix_lasso,  directed = FALSE, normalized = FALSE, cutoff = -1)
ps_w5_betweenness_lasso_raw = betweenness(psw5_network_adjmatrix_lasso,  directed = FALSE, normalized = FALSE, cutoff = -1)
ps_w9w10_betweenness_lasso_raw
ps_w5_betweenness_lasso_raw


ps_w9w10_mean_distance = mean_distance(psw9w10_network_adjmatrix_lasso, directed = FALSE)
ps_w9w10_mean_distance
ps_w5_mean_distance = mean_distance(psw5_network_adjmatrix_lasso, directed = FALSE)
ps_w5_mean_distance

ps_w9w10_distances = distances(psw9w10_network_adjmatrix_lasso)[upper.tri(distances(psw9w10_network_adjmatrix_lasso))]
ps_w9w10_distances = ps_w9w10_distances[is.finite(ps_w9w10_distances)]
mean(ps_w9w10_distances)
sd(ps_w9w10_distances)
length(ps_w9w10_distances)

ps_w5_distances = distances(psw5_network_adjmatrix_lasso)[upper.tri(distances(psw5_network_adjmatrix_lasso))]
ps_w5_distances = ps_w5_distances[is.finite(ps_w5_distances)]
mean(ps_w5_distances)
sd(ps_w5_distances)
length(ps_w5_distances)


taxa_degrees_ps_w9w10 = rowSums(psw9w10_network_matrix_lasso)
zero_connections_ps_w9w10 = sum(taxa_degrees_ps_w9w10 == 0)
zero_connections_ps_w9w10

taxa_degrees_ps_w5 = rowSums(psw5_network_matrix_lasso)
zero_connections_ps_w5 = sum(taxa_degrees_ps_w5 == 0)
zero_connections_ps_w5
