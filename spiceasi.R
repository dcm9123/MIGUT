# Daniel Castaneda Mogollon, PhD
# May 27th, 2026
# This script will run the network co-association between taxa by using spiceasi

library("SpiecEasi")
library("phyloseq")
library("Matrix")
library("igraph")
library("qgraph")

set.seed(23)

path = "/Users/danielcm/Desktop/MIGUT"
setwd(path)

#f_input = read.table(paste0(path,"/analysis_may2026/filtered_metaphlan4_abundances/merged_absolute_abundances_table_GTDB_filtered.txt"), header = TRUE, sep = "\t", row.names = 1)
# Generating the ps object

raw_input = read.table(paste0(path,"/analysis_may2026/merged_tables_abundances/merged_absolute_abundances_table_GTDB.txt"), header = TRUE, sep = "\t", row.names = 1)
metadata_df = read.csv(paste0(path,"/metadata","/MIGUT_metadata_in_tab_format.txt"), header = TRUE, sep = "\t")

metadata_df = as.data.frame(metadata_df)

otu_table = as.matrix(raw_input)
taxa_table = as.matrix(read.csv(paste0(path,"/analysis_may2026/merged_tables_abundances/GTDB_taxa_table.csv"), header = TRUE))
rownames(taxa_table) = taxa_table[,1]

sample_data = metadata_df[,c("Study_Number","ID","Age","Sex_Male1_Female2","Group","BMI","Constipation_No0_Yes1","Group","Group_and_sex","Males_only","Females_only")]
rownames(sample_data) = sample_data$ID
ps = phyloseq(otu_table(otu_table, taxa_are_rows = TRUE), tax_table(taxa_table), sample_data(sample_data))


#SPIEC-EASI
# This software works by inferring microbial association networks from composituional microbiome data
# It does not simply generate a matrix of Pearson correlation coefficients, but instead it tries to determine or 
# generate a newtork of conditional dependencies between taxa. Standard correlations using relative abundances may provide 
# false associations simply by having different absolute abundances.

# It asks if taxa X and Y are associated after taking into account the rest of the taxa in the system/dataset.
# Rows = samples, columns = taxa. Use raw or filtered counts.




# Filters applied: Keep only species level taxa, remove those present less than 5% of samples, remove samples with no abundances, 
# removes taxa with a mean relative abundance of less than 0.001%.

applying_phyloseq_filters = function(ps_object, sample_prevalence_filter, taxa_mean_prevalence, category_to_analyze, category_column){
  
  if(!(category_column %in% colnames(sample_data(ps_object)))){
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
    tax_rank_values
  )
}

plotting_networks = function(network_to_plot, spic_network, category, ps_object, tax_rank, plot_file_name, tax_rank_colors = NULL){
  if(spic_network == TRUE){
    igraph_object = adj2igraph(getRefit(network_to_plot), rmEmptyNodes = TRUE)
  }
  else{
    igraph_object = adj2igraph(network_to_plot)
  }
  
  taxa_df = as.data.frame(tax_table(ps_object))

  beta_coefficients = getOptBeta(network_to_plot)
  beta_symmetric = symBeta(beta_coefficients, mode = "maxabs")
  edge_list = ends(igraph_object, E(igraph_object), names = TRUE)
  
  E(igraph_object)$weight = beta_symmetric[cbind(edge_list[,1], edge_list[,2])]
  E(igraph_object)$abs_weight <- abs(E(igraph_object)$weight)
  E(igraph_object)$sign <- ifelse(E(igraph_object)$weight > 0, "positive", "negative")

  tax_rank_prefix = paste0(tolower(substr(tax_rank, 1, 1)), "__")
  V(igraph_object)$tax_rank = taxa_df[V(igraph_object)$name, tax_rank]
  V(igraph_object)$tax_rank = sub(paste0("^", tax_rank_prefix), "", V(igraph_object)$tax_rank)

  cutoff <- quantile(E(igraph_object)$abs_weight, 0.50, na.rm = TRUE)

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
  classes = sort(unique(V(igraph_filt)$tax_rank))
  class_colors = tax_rank_colors[classes]
  V(igraph_filt)$color = class_colors[V(igraph_filt)$tax_rank]

  print(head(V(igraph_filt)$tax_rank))


  set.seed(123)

  layout_fr <- layout_with_fr(
    igraph_filt,
    weights = E(igraph_filt)$layout_weight
  )
  
  vertex_size = clr(rowMeans(otu_table(ps_object),1)) + 4
  png(plot_file_name, height = 8, width = 14, units = "in", res = 600)
  plot(igraph_filt,  vertex.size = vertex_size, vertex.label = NA, main = category,
        edge.width = scales::rescale(E(igraph_filt)$abs_weight, to = c(1, 3)), 
        edge.color = ifelse(E(igraph_filt)$sign == "positive", "black", "red"),
        layout = layout_fr, vertex.color = V(igraph_filt)$color)

  legend(
      "topleft",
      legend = names(class_colors),
      col = class_colors,
      pch = 19,
      pt.cex = 1.5,
      bty = "n",
      cex = 0.8
  )
  dev.off()
  return(igraph_object)
  }


# Getting ps objects
ps_pd_only = applying_phyloseq_filters(ps, sample_prevalence_filter = 0.10, taxa_mean_prevalence = 0.00000 ,category_to_analyze = "PD", category_column = "Group")
ps_hc_only = applying_phyloseq_filters(ps, sample_prevalence_filter = 0.10, taxa_mean_prevalence = 0.00000 ,category_to_analyze = "HC", category_column = "Group")
#ps_male_only = applying_phyloseq_filters(ps, sample_prevalence_filter = 0.05, taxa_mean_prevalence = 0.00001 ,category_to_analyze = "1", category_column = "Sex_Male1_Female2")
#ps_female_only = applying_phyloseq_filters(ps, sample_prevalence_filter = 0.05, taxa_mean_prevalence = 0.00001 ,category_to_analyze = "2", category_column = "Sex_Male1_Female2")
#ps_pd_male_only = applying_phyloseq_filters(ps, sample_prevalence_filter = 0.05, taxa_mean_prevalence = 0.00001 ,category_to_analyze = "PD_Male", category_column = "Group_and_sex")
#ps_pd_female_only = applying_phyloseq_filters(ps, sample_prevalence_filter = 0.05, taxa_mean_prevalence = 0.00001 ,category_to_analyze = "PD_Female", category_column = "Group_and_sex")

# Generating best model and plotting_networks, SPIC method
pd_network = network_generation(ps_pd_only, method_to_use = "mb", lambda_min_ratio = 1e-2, n_lambda = 30, pulsar_rep_num = 50, pulsar_ncores = 8, sparcc_method = FALSE)
hc_network = network_generation(ps_hc_only, method_to_use = "mb", lambda_min_ratio = 1e-2, n_lambda = 30, pulsar_rep_num = 50, pulsar_ncores = 8, sparcc_method = FALSE)
#male_network = network_generation(ps_male_only, method_to_use = "mb", lambda_min_ratio = 1e-2, n_lambda = 30, pulsar_rep_num = 50, pulsar_ncores = 8, sparcc_method = FALSE)
#female_network = network_generation(ps_female_only, method_to_use = "mb", lambda_min_ratio = 1e-2, n_lambda = 30, pulsar_rep_num = 50, pulsar_ncores = 8, sparcc_method = FALSE)
#pd_male_network = network_generation(ps_pd_male_only, method_to_use = "mb", lambda_min_ratio = 1e-2, n_lambda = 30, pulsar_rep_num = 50, pulsar_ncores = 8, sparcc_method = FALSE)
#pd_female_network = network_generation(ps_pd_female_only, method_to_use = "mb", lambda_min_ratio = 1e-2, n_lambda = 30, pulsar_rep_num = 50, pulsar_ncores = 8, sparcc_method = FALSE)

ps_pd_only
ps_hc_only
#pd_network$refit
#pd_network$est
#pd_network$fun
#pd_network$select
#pd_network$lambda

#class(pd_network)
# Generating best model and plotting networks, SPARCC method

# CONNECTIONS
pd_network_matrix = as.matrix(getRefit(pd_network))
total_connections_pd = sum(pd_network_matrix/2)
total_connections_pd

hc_network_matrix = as.matrix(getRefit(hc_network))
total_connections_hc = sum(hc_network_matrix/2)
total_connections_hc

# CLOSENNES ()
pd_network_graph = adj2igraph(getRefit(pd_network), rmEmptyNodes = TRUE,vertex.attr = list(name = colnames(ps_pd_only)))
pd_closeness = harmonic_centrality(pd_network_graph, mode = "all", normalized = FALSE, cutoff = -1)
pd_closeness

hc_network_graph = adj2igraph(getRefit(hc_network), rmEmptyNodes = TRUE,vertex.attr = list(name = colnames(ps_hc_only)))
hc_closeness = harmonic_centrality(hc_network_graph, mode = "all", normalized = FALSE, cutoff = -1)
hc_closeness

for(item in pd_closeness){
  cat((log10(item)),"\n")
}

for(item in hc_closeness){
  cat((log10(item)),"\n")
}

# DENSITY
density_pd = edge_density(pd_network_graph, loops = FALSE)
density_pd

density_hc = edge_density(hc_network_graph, loops = FALSE)
density_hc


# BETWEENNESS
pd_betweenness = betweenness()


pd_betweenness_norm = betweenness(pd_network_graph,  directed = FALSE, normalized = TRUE, cutoff = -1)
hc_betweenness_norm = betweenness(hc_network_graph,  directed = FALSE, normalized = TRUE, cutoff = -1)

pd_betweenness_raw = betweenness(pd_network_graph,  directed = FALSE, normalized = FALSE, cutoff = -1)
hc_betweenness_raw = betweenness(hc_network_graph,  directed = FALSE, normalized = FALSE, cutoff = -1)


for(item in pd_betweenness_raw){
  cat((item),"\n")
}

for(item in hc_betweenness_raw){
  cat((item),"\n")
}

pd_mean_distance = mean_distance(pd_network_graph, directed = FALSE)
pd_mean_distance
hc_mean_distance = mean_distance(hc_network_graph, directed = FALSE)
hc_mean_distance

pd_distances = distances(pd_network_graph)[upper.tri(distances(pd_network_graph))]
pd_distances = pd_distances[is.finite(pd_distances)]
mean(pd_distances)
sd(pd_distances)
length(pd_distances)
as.matrix(pd_distances)

hc_distances = distances(hc_network_graph)[upper.tri(distances(hc_network_graph))]
hc_distances = hc_distances[is.finite(hc_distances)]
mean(hc_distances)
sd(hc_distances)
length(hc_distances)

components(pd_network_graph)$no
components(hc_network_graph)$no
components(pd_network_graph)$csize
components(hc_network_graph)$csize
components(pd_network_graph)


possible_pd_connections_max = (ntaxa(ps_pd_only) * (ntaxa(ps_pd_only)-1))/2
possible_hc_connections_max = (ntaxa(ps_hc_only) * (ntaxa(ps_hc_only)-1))/2

taxa_degrees_pd = rowSums(pd_network_matrix)
zero_connections_pd = sum(taxa_degrees_pd == 0)
zero_connections_pd

taxa_degrees_hc = rowSums(hc_network_matrix)
zero_connections_hc = sum(taxa_degrees_hc == 0)
zero_connections_hc



# Plotting plotting_networks
class_colors = make_tax_rank_colors(ps, "Class")
order_colors = make_tax_rank_colors(ps, "Order")
phylum_colors = make_tax_rank_colors(ps, "Phylum")

adj_matrix_pd = as_adjacency_matrix(plotting_networks(pd_network, spic_network = TRUE, category = "PD", ps_object = ps_pd_only, 
                  tax_rank = "Class", plot_file_name = "pd_top_0.75_network.png",
                  tax_rank_colors = class_colors))

rownames(adj_matrix_pd) = colnames(ps_pd_only)
colnames(adj_matrix_pd) = colnames(ps_pd_only)

adj_matrix_hc = as_adjacency_matrix(plotting_networks(hc_network, spic_network = TRUE, category = "HC", ps_object = ps_hc_only, 
                  tax_rank = "Class", plot_file_name = "hc_top_0.75_network.png",
                  tax_rank_colors = class_colors))

rownames(adj_matrix_hc) = colnames(ps_hc_only)
colnames(adj_matrix_hc) = colnames(ps_hc_only)

taxa_degrees_pd = sort(rowSums(adj_matrix_pd), decreasing = TRUE)
taxa_degrees_hc = sort(rowSums(adj_matrix_hc), decreasing = TRUE)

for(item in taxa_degrees_pd){
  cat(item,"\n")
}

for(item in taxa_degrees_hc){
  cat(item,"\n")
}


plotting_networks(pd_network, spic_network = TRUE, category = "PD", ps_object = ps_pd_only, 
                  tax_rank = "Phylum", plot_file_name = "pd_top_0.50_phylum_network.png",
                  tax_rank_colors = phylum_colors)

plotting_networks(hc_network, spic_network = TRUE, category = "HC", ps_object = ps_hc_only, 
                  tax_rank = "Phylum", plot_file_name = "hc_top_0.50_phylum_network.png",
                  tax_rank_colors = phylum_colors)


tax_table(ps_pd_only) == Species

#plotting_networks(pd_male_network, vertex_extra = 2, spic_network = TRUE, category = "PD_Male", ps_object = ps_pd_male_only)
