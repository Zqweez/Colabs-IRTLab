# Percent Identity (PID) Based Clustering -------------------------------------------------------
library(ape)
library(seqinr)
library(Biostrings)
library(ggplot2)
library(stringr)
library(dplyr)


# Load MSA file
msa_file <- "data/sanger/03_ab1_assembled_sequences/ALL_SAMPLES_consensus.fasta"
dna_aln <- read.alignment(msa_file, format = "fasta")

# Process sequences
clean_seqs <- sapply(dna_aln$seq, function(seq) gsub("\\s+", "", seq))
sample_names <- sub("_consensus", "", dna_aln$nam)
df <- data.frame(
  sample = sub("_$", "", sample_names),
  sequence = clean_seqs,
  stringsAsFactors = FALSE
)
df$sample <- sub("^(Takasu|Kasper)_", "", df$sample)
df$group <- sub("-[12]$", "", df$sample)
df <- df[order(df$sample), ]
rownames(df) <- NULL
df$ambiguous_percent <- str_count(df$sequence, "n") / nchar(df$sequence)
df$sequence_length <- nchar(df$sequence)

dna <- DNAStringSet(df$sequence) 
names(dna) <- sprintf("%s", df$sample)
writeXStringSet(dna, filepath = "data/sanger/formatted_consensus.fasta", format   = "fasta")

# Prepare the taxa
taxa_strains <- read.csv("data/sanger/blast_top_hit_ezbio.tsv", sep = "\t", header = TRUE)
taxa_strains <- taxa_strains %>%
  mutate(
    genus = sub(",s:.*", "", sub(".*g:", "", taxa)),
    species = paste(genus, sub("^.* ", "", sub(".*s:", "", taxa))),
    query = sub("^(Takasu|Kasper)_", "", query),
    pident =  round(pident, 1)
  ) %>%
  select(query, genus, species, pident)

write.csv(taxa_strains, file = "data/sanger/taxa-genus-species-pid.csv")

# Setup Alignment
alignment_type <- "local" # or "global"
gap_opening <- -5
gap_extension <- -2

sub_df <- df
dna_set <- DNAStringSet(sub_df$sequence)
names(dna_set) <- sub_df$sample

# Create empty matrix
alignment_results <- data.frame()
n <- length(dna_set)

sequence_lengths <- setNames(nchar(df$sequence), df$sample) # To calculate weighed PID score
# Fill with pairwise percent identity use either the df for both local and global or matrix for only one
pb <- progress_bar$new(total = n * n)
for (i in 1:n) {
  for (j in 1:n) {
    pb$tick() # Progress bar
    aln  <- pairwiseAlignment(dna_set[[i]], dna_set[[j]], type = alignment_type, substitutionMatrix = NULL, 
                                    gapOpening = gap_opening, gapExtension = gap_extension)
    
    # To weigh the alignment by the coverage, only for local, returns 1 for global
    min_length <- min(sequence_lengths[[names(dna_set)[i]]], sequence_lengths[[names(dna_set)[j]]])
    coverage <- nchar(gsub("-", "", as.character(alignedPattern(aln)))) / min_length
    
    # Print for debugging
    # print(sprintf("Sample %s - %s coverage: %s, pid: %s, w_pid: %s", names(dna_set)[i], names(dna_set)[j], coverage, pid(aln_local), pid(aln)*coverage))
    
    alignment_results <- bind_rows(alignment_results, data.frame(
      Sample1 = names(dna_set)[i],
      Sample2 = names(dna_set)[j],
      AlignmentType = alignment_type,
      Score = score(aln),
      PID = pid(aln),
      pid_weighted = pid(aln)*coverage,
      AlignedWidth =width(alignedPattern(aln))
    )) 
  }
}

# Convert df to a matrix for distance matrix downsteams
identity_matrix <- alignment_results %>%
  select(Sample1, Sample2, pid_weighted) %>%
  pivot_wider(names_from = Sample2, values_from = pid_weighted) %>%
  tibble::column_to_rownames("Sample1")  # make Sample1 the rownames

# Heatmap
library(viridis)
heat_colors <- c("#000000", "#120054", "#71047a", "#cf2649", "#fe9564", "#ffdda1", "#ffffda") # c("white", "white", "lightblue", "steelblue", "grey10")
ggplot(alignment_results, aes(Sample1, Sample2, fill = pid_weighted)) +
  geom_tile(color = "white") +
  scale_fill_gradientn(colors = heat_colors, values = scales::rescale(c(0, 30, 40, 60, 80, 97, 100)), limits = c(0, 100)) +
  #scale_fill_viridis(option = "G", name = "Abundance", direction = 1, na.value = "white") +
  xlab("") + ylab("") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6), axis.text.y = element_text(size = 6)) +
  labs(title = "Pairwise Aligned (local) Sequence Identity Heatmap", subtitle = "Weighted by alignment length / shortest sequence", fill = "% Identity")

ggsave("figures/Heatmap-local-weighted-align-fire.pdf", width = 8, height = 8)

## Hierarcheal clustering of the pid data
distance_matrix <- 100 - identity_matrix  # Low pid => high distance
distance_obj <- as.dist(distance_matrix)
hc <- hclust(distance_obj, method = "average") # "average", "complete", "ward.D2", "single", "centroid", "median"

pdf("figures/Hierarchical-clustering-avg-weighted-local-all-isolates.pdf", width = 8, height = 5)
plot(hc, main = "Hierarchical Clustering of Sequences", #hang = -1,
     xlab = "Samples", cex = 0.8)
mtext("Local weighted alignment", side = 3, line = 0.5)
dev.off()

## Combined heatmap
library(pheatmap)       # clustered heatmaps
library(dichromat)
library(ComplexHeatmap)
library(dendextend)
library(circlize)
library(grid)
library(RColorBrewer)
library(colorspace) 
library(Polychrome)

pid_mat <- as.matrix(identity_matrix)
dist_obj <- as.dist(1 - pid_mat / 100)   # values 0-1; 0 = identical
combined_cols <- colorRamp2(
  c(0, 90, 100),
  c("#152525", "#c1c6c6", "#fefefe")             # light grey → almost-black
)
#colorRampPalette(c("#000000", "#fffddd"))(100)

pheatmap(
  pid_mat,
  color                = combined_cols,
  cluster_rows         = TRUE,
  cluster_cols         = TRUE,
  clustering_distance_rows = dist_obj,   # use custom distances
  clustering_distance_cols = dist_obj,
  clustering_method    = "average",      # UPGMA; try "complete", "ward.D2",
  display_numbers      = TRUE,          # or TRUE if you want the PID printed
  number_format = "%.0f",
  treeheight_row  = 70,   # default = 50  increase for tall row trees
  treeheight_col  = 70,   # same for column tree
  cellwidth       = NA,    # let pheatmap shrink cells automatically
  cellheight      = NA,
  border_color    = NA,
  fontsize = 6,
  legend_breaks        = c(0, 50, 100),
  legend_labels        = c("0", "50", "100 % PID"),
  main                 = "Pairwise PID clustered heat-map"
)

# Complex heatmap
ht <- Heatmap(
  pid_mat,
  name        = "PID",
  col         = combined_cols,
  clustering_distance_rows    = dist_obj,
  clustering_distance_columns = dist_obj,
  clustering_method_rows      = "average",
  clustering_method_columns   = "average",
  # ── move both dendrogram and labels to the LEFT ──
  row_dend_side  = "left",
  row_names_side = "left",
  row_dend_width = unit(15, "mm"),
  
  # keep column dendrogram at top (default) but you could also move:
  #column_dend_side  = "bottom",
  #column_names_side = "bottom",
  #column_dend_height = unit(40, "mm"),
  show_column_dend = F,
  
  cell_fun = function(j, i, x, y, width, height, fill) {
    identity_i <- pid_mat[i, j]
    color_ide <- ifelse(identity_i < 52, "#bbb", "#000")
    grid.text(sprintf("%.0f", identity_i),        # one decimal place
              x, y,
              gp = gpar(fontsize = 3, col = color_ide))# change to white if fill is dark
  },
  
  # cosmetic
  row_names_gp = gpar(fontsize = 7),      # shrink row text if many samples
  column_names_gp = gpar(fontsize = 7),
  heatmap_legend_param = list(
    at = c(0, 50, 100),
    labels = c("0 %", "50 %", "100 %")
  ),
  column_title      = "Combined Heatmap of Percent Identity (%) \n ",
  column_title_gp   = gpar(fontsize = 12, fontface = "bold")
)

ht_draw <- draw(ht)
dend <- row_dend(ht_draw)
k <- 15
row_clu <- cutree(dend, k = k)   # named vector: sample → cluster_id

clu_genus <- tibble(sample = names(row_clu),
                    cluster = row_clu) %>%
  left_join(taxa_strains, by = c("sample" = "query")) %>%
  count(cluster, genus, name = "n") %>%
  group_by(cluster) %>%
  mutate(prop = n / sum(n)) %>%
  slice_max(order_by = prop, n = 1, with_ties = FALSE) %>%   # top genus
  mutate(genus_lab = if_else(prop >= 0.7, genus, "Mixed"))   # threshold 60 %

# define colors for genera
cluster_levels <- as.character(unique(clu_genus$genus_lab))
n_clusters     <- length(cluster_levels)
# Colors using polychrome
set.seed(0)
gen_cols <- createPalette(n_clusters, c("#1b9e77", "#d95f02", "#7570b3"))
pastelise <- function(cols, desat = 0.6, lighten = 0.3) {
  cols %>%
    desaturate(amount = desat) %>%   # pull chroma toward grey
    lighten  (amount = lighten)      # push luminance up toward white
}
gen_cols <- pastelise(gen_cols)
names(gen_cols) <- c(cluster_levels)

# Order the annotation in a nice order
cluster_to_genus <- setNames(clu_genus$genus_lab, clu_genus$cluster)
genus_per_row <- cluster_to_genus[as.character(row_clu)]
row_ord <- row_order(ht_draw)             # numeric indices (top→bottom)
samp_ord <- rownames(pid_mat)[row_ord]    # sample IDs in display order
clust_ord <- row_clu[samp_ord]            # cluster ID per displayed sample
genus_ord <- cluster_to_genus[as.character(clust_ord)]
legend_levels <- unique(genus_ord)        # first time each genus appears

row_genus_vec <- factor(genus_per_row, levels = legend_levels)
gen_cols      <- gen_cols[legend_levels]
row_anno <- rowAnnotation(
  Genus = row_genus_vec,
  col   = list(Genus = gen_cols),
  width = unit(4, "mm"),
  show_legend = TRUE
)

ht_full <- ht + row_anno

pdf("figures/Heatmap-genus-complete.pdf", width = 9, height = 8, useDingbats = FALSE)  # safer for non-standard fonts
draw(ht_full)
grid.text(
  "Local weighted alignment, average linkage",
  x = unit(0.5, "npc"),                    # centred
  y = unit(1, "npc") - unit(24, "pt"),                   
  gp = gpar(fontsize = 10)
)
dev.off()


## Do PCA with the same distance matrix as above
pca <- prcomp(distance_matrix, scale. = TRUE)
summary(pca)

pca_var <- pca$sdev^2
pca_var <- pca_var / sum(pca_var)

pca_df <- as.data.frame(pca$x)
pca_df$Sample <- rownames(pca_df)
pca_df$group = df$group


ggplot(pca_df, aes(x = PC1, y = PC2, color = group, label = Sample)) +
  geom_point(size = 2) +
  geom_text(vjust = -1, size = 1.5, , show.legend = FALSE) +
  theme_minimal() +
  xlab(sprintf("PC1 (%s)%%", round(100 * pca_var[1], 1))) +
  ylab(sprintf("PC2 (%s)%%", round(100 * pca_var[2], 1))) +
  labs(title = "PCA of Pairwise Identity local weighted")

ggsave("figures/PCA-PID-weighted-PC1-PC2.pdf", width = 8, height = 5)

## K-means from the identity matrix and plotted with PCA
set.seed(42)
k = 8
kmeans_result <- kmeans(identity_matrix, centers = k)  # adjust k
kmeans_result <- kmeans(pca_df[,1:6], centers = k)  # k-means built with the top x PCs

pca_df$Cluster <- as.factor(kmeans_result$cluster)

# Do an elbow plot if wanted to choose K
wss <- sapply(1:10, function(k) {
  kmeans(pca_df[,1:6], centers = k, nstart = 10)$tot.withinss # update to either matrix or pca_df
})
plot(1:10, wss, type = "b", main = "Elbow Method for K", xlab = "Number of Clusters")

ggplot(pca_df, aes(x = PC1, y = PC2, color = Cluster, label = Sample)) +
  geom_point(size = 2) +
  geom_text(vjust = -1, size = 1.5, show.legend = FALSE) +
  theme_minimal() +
  labs(title = "PCA of weighted PID + K-means Clustering", subtitle = sprintf("Based on PC1-6, K=%s", k))
ggsave("figures/K-means-PC1-6-local-wieghted.pdf", width = 8, height = 5)
