# install.packages(c("ape", "seqinr"))
library(ape)
library(seqinr)
library(Biostrings)
library(ggplot2)
library(stringr)
library(tidyr)
library(dplyr)
#install.packages("progress")
library(progress)
library(dendextend)
library(doSNOW)
library(foreach)


fasta_A <- "data/sanger/formatted_consensus.fasta"
fasta_B <- "results/ngs/top_ASVs_per_taxa.fasta"

# Read the alignment
seqs_A <- read.alignment(fasta_A, format = "fasta")
seqs_B <- read.alignment(fasta_B, format = "fasta")

# Extract sample names
clean_names <- function(x) {
  x %>%
    str_remove("_consensus$") %>%   # drop trailing "_consensus"
    str_remove("_$")               # drop trailing underscore, if any
}
names_A <- clean_names(seqs_A$nam)
names_B <- clean_names(seqs_B$nam)

seqs_A <- sapply(seqs_A$seq, function(seq) gsub("\\s+", "", seq))
seqs_B <- sapply(seqs_B$seq, function(seq) gsub("\\s+", "", seq))

# Create dataframe
df <- data.frame(
  sample = c(names_A, names_B),
  sequence = c(seqs_A, seqs_B),
  stringsAsFactors = FALSE
)

df_clean <- df %>%
  filter(!str_ends(sample, "KS-02-2")) %>%
  mutate(base_id = str_remove(sample, "-[12]$")) %>% 
  group_by(base_id) %>%
  dplyr::slice(1) %>%                                           
  ungroup() %>%
  select(-sample) %>%                           
  dplyr::rename(sample = base_id) 

# Extract common prefix by removing the last "-1" or "-2"
df_clean <- df_clean[order(df_clean$sample),]
rownames(df_clean) <- NULL # Reset the indexing, (rownames)

# Count the number of ambiguous bases
df_clean$sequence_length <- nchar(df_clean$sequence)
df_clean$ambiguous_percent <- str_count(df_clean$sequence, "n")/nchar(df_clean$sequence)

# ---- Compute the Pairwise alignment of each of the consensus generating a large table -----

# Convert to DNAStringSet use only a subset of the df if desired
sub_df <- df_clean
dna_set <- DNAStringSet(sub_df$sequence)
names(dna_set) <- sub_df$sample

# Create empty matrix
alignment_results <- data.frame()
n <- length(dna_set)

alignment_type <- "local" # or "global"
gap_opening <- -5
gap_extension <- -2

sequence_lengths <- setNames(nchar(sub_df$sequence), sub_df$sample) # To calculate weighed PID score
# Fill with pairwise percent identity use either the df for both local and global or matrix for only one
pb <- progress_bar$new(total = n * n)
for (i in 1:n) {
  for (j in 1:n) {
    pb$tick() # Progress bar
    aln  <- pwalign::pairwiseAlignment(dna_set[[i]], dna_set[[j]], type = alignment_type, substitutionMatrix = NULL, 
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

## ---- ##
num_cores <- detectCores() - 4
if (is.na(num_cores) || num_cores < 1) {
  num_cores <- 1
} # Failsafe for single-core machines
cl <- makeCluster(num_cores)
registerDoSNOW(cl)
cat(paste("Registered", num_cores, "cores for parallel processing.\n"))
# Generate all unique pairs of indices to iterate over
pairs <- combn(1:n, 2, simplify = FALSE)

pb_parallel <- txtProgressBar(max = length(pairs), style = 3)
progress <- function(n) setTxtProgressBar(pb_parallel, n)
opts <- list(progress = progress)

# Use foreach to process in parallel
alignment_results_parallel <- foreach(
  p = pairs,
  .combine = 'bind_rows',
  .packages = c('pwalign', 'dplyr'),
  .export = c("dna_set", "alignment_type", "gap_opening", "gap_extension", "sequence_lengths"),
  .options.snow = opts
) %dopar% {
  # Get the indices for the current pair
  i <- p[1]
  j <- p[2]
  
  # Perform pairwise alignment (same logic as before)
  aln <- pwalign::pairwiseAlignment(dna_set[[i]], dna_set[[j]],
                                    type = alignment_type,
                                    substitutionMatrix = NULL,
                                    gapOpening = gap_opening,
                                    gapExtension = gap_extension
  )
  
  # Calculate weighted PID
  min_length <- min(sequence_lengths[[names(dna_set)[i]]], sequence_lengths[[names(dna_set)[j]]])
  coverage <- nchar(gsub("-", "", as.character(alignedPattern(aln)))) / min_length
  
  # Return a data frame for this single pair. `foreach` will collect and combine them.
  data.frame(
    Sample1 = names(dna_set)[i],
    Sample2 = names(dna_set)[j],
    AlignmentType = alignment_type,
    Score = score(aln),
    PID = pid(aln),
    pid_weighted = pid(aln) * coverage,
    AlignedWidth = width(alignedPattern(aln))
  )
}

# Stop the cluster when its done
close(pb_parallel)
stopCluster(cl)
cat("Parallel Processing Finished!")

mirrored_results <- alignment_results_parallel %>%
  rename(
    Sample1_temp = Sample1,
    Sample2_temp = Sample2
  ) %>%
  rename(
    Sample1 = Sample2_temp,
    Sample2 = Sample1_temp
  )
diagonal_results <- data.frame(
  Sample1 = names(dna_set),
  Sample2 = names(dna_set),
  pid_weighted = 100
)

alignment_results <- bind_rows(alignment_results_parallel, mirrored_results, diagonal_results)

# Get KS and ST taxa
taxa_strains <- read.csv("data/sanger/taxa-genus-species-pid.csv", sep = ",", header = TRUE)

taxa_strains <- taxa_strains %>%
  filter(!str_ends(query, "KS-02-2")) %>%
  mutate(query = sub("-[12]$", "", query)) %>%
  group_by(query) %>%
  dplyr::slice(1) %>%                                           
  ungroup() %>%
  rename(ASV = query)

taxa_asvs <- read.csv("results/ngs/taxa_files/blast_top_hit_compact.tsv", sep = "\t", header = TRUE)
taxa_asvs <- taxa_asvs %>% 
  mutate(
    genus = sub(",s:.*", "", sub(".*g:", "", taxa)),
    species = paste(genus, sub("^.* ", "", sub(".*s:", "", taxa)))
  ) %>%
  rename(ASV = query)

combinded_taxa <- bind_rows(taxa_asvs, taxa_strains) %>%
  select(ASV, species)
# Add taxa to the alignment results
alignment_taxa <- alignment_results %>%
  left_join(combinded_taxa, by = c("Sample1" = "ASV")) %>%
  rename(Taxa1 = species) %>%
  left_join(combinded_taxa, by = c("Sample2" = "ASV")) %>%
  rename(Taxa2 = species) %>%
  mutate(
    Sample1 = paste0(Sample1, "-", Taxa1),
    Sample2 = paste0(Sample2, "-", Taxa2)
  )

# Convert df to a matrix for distance matrix downsteams
identity_matrix <- alignment_taxa %>%
  select(Sample1, Sample2, pid_weighted) %>%
  pivot_wider(names_from = Sample2, values_from = pid_weighted) %>%
  tibble::column_to_rownames("Sample1")  # make Sample1 the rownames

# Make sure row- and col-names are identical AND in the same order
ids <- intersect(rownames(identity_matrix), colnames(identity_matrix))
identity_matrix <- identity_matrix[ids, ids] 

## Hierarcheal clustering of the pid data
identity_matrix_sub <- identity_matrix

distance_matrix <- 100 - identity_matrix_sub  # Low pid => high distance
distance_obj <- as.dist(distance_matrix)
hc <- hclust(distance_obj, method = "complete") # "average", "complete", "ward.D2", "single", "centroid", "median"

# Save current plotting parameters
old_par <- par(no.readonly = TRUE)

# Method 1: Vertical dendrogram with colors
dend <- as.dendrogram(hc)
labels_dend <- labels(dend)
label_colors <- ifelse(startsWith(labels_dend, "ASV"), "black", "red2")

dend_colored <- dend %>%
  set("labels_col", value = label_colors) %>%
  set("labels_cex", 0.3)

# Horizontal dendrogram with colors
pdf("results/ngs/figures/Hierarchical-clustering-avg-weighted-local-ASVs-horizontal.pdf", width = 8, height = 12)
par(mar = c(5.1, 4.1, 4.1, 8))
plot(dend_colored, main = "Hierarchical Clustering of Sequences", 
     horiz = TRUE, xlab = "Distance")
mtext("Local weighted alignment", side = 3, line = 0.2, cex = 0.7)
dev.off()

# Reset plotting parameters
par(old_par)


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

ggsave("figures/PCA-PID-weighted-PC1-PC2-KS-only.pdf", width = 8, height = 5)

## K-mers from the identity matrix and plotted with PCA
set.seed(42)
k = 6
kmeans_result <- kmeans(identity_matrix, centers = k)  # adjust k
kmeans_result <- kmeans(pca_df[,1:6], centers = k)  # k-means built with the top x PCs

pca_df$Cluster <- as.factor(kmeans_result$cluster)

# Do an elbow plot if wanted to choose K
# wss <- sapply(1:10, function(k) {
#   kmeans(pca_df[,1:6], centers = k, nstart = 10)$tot.withinss # update to either matrix or pca_df
# })
# plot(1:10, wss, type = "b", main = "Elbow Method for K", xlab = "Number of Clusters")

ggplot(pca_df, aes(x = PC1, y = PC2, color = Cluster, label = Sample)) +
  geom_point(size = 2) +
  geom_text(vjust = -1, size = 1.5, show.legend = FALSE) +
  theme_minimal() +
  labs(title = "PCA of weighted PID + K-means Clustering", subtitle = sprintf("Based on PC1-6, K=%s", k))
ggsave("figures/K-means-PC1-6-local-wieghted.pdf", width = 8, height = 5)
