# K-mer Based Clustering -------------------------------------------------------
library(ape)
library(seqinr)
library(Biostrings)
library(ggplot2)
library(stringr)
library(dplyr)

# Load MSA file
msa_file <- "02_ab1_assembled_sequences/ALL_SAMPLES_consensus.fasta"
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

# K-mer frequency matrix
k-mers <- 4 # Try different K-mer values
cleaned_sequences <- gsub("N", "", df$sequence)
dna_set <- DNAStringSet(cleaned_sequences)
kmer_matrix <- oligonucleotideFrequency(dna_set, width = k-mers)
rownames(kmer_matrix) <- df$sample

# PCA
pca <- prcomp(kmer_matrix, scale. = TRUE)
summary(pca)
pca_var <- pca$sdev^2
pca_var <- pca_var / sum(pca_var)

pca_df <- as.data.frame(pca$x)
pca_df$sample <- df$sample
pca_df$group = df$group


ggplot(pca_df, aes(x = PC1, y = PC2, color = group, label = sample)) +
  geom_point(size = 2) +
  geom_text(vjust = -1, size = 1.5, , show.legend = FALSE) +
  theme_minimal() +
  xlab(sprintf("PC1 (%s)%%", round(100 * pca_var[1], 1))) +
  ylab(sprintf("PC2 (%s)%%", round(100 * pca_var[2], 1))) +
  labs(title = sprintf("PCA of DNA %s-mer Frequency", nchar(colnames(kmer_matrix)[1])))

ggsave(sprintf("figures/PCA-PC1-PC2-%s-mer.pdf", nchar(colnames(kmer_matrix)[1])), width = 8, height = 5)

# Hierarchical clustering
d <- dist(kmer_matrix, method = "euclidean")
clust_method <- "average"
hc <- hclust(d, method = clust_method) # try "average", "complete", "ward.D2"

# Plot and save dendrogram
pdf(sprintf("figures/Hierarchical-clustering-%smer-%s.pdf", nchar(colnames(kmer_matrix)[1]), clust_method), width = 8, height = 5)
plot(hc, labels = rownames(kmer_matrix), main = sprintf("Hierarchical Clustering of %s-mers", nchar(colnames(kmer_matrix)[1])), xlab = "", cex = 0.7)
dev.off()

# K-means clustering
set.seed(42)
k <- 6  # adjust as needed
km <- kmeans(kmer_matrix, centers = k)
pca_df$cluster <- as.factor(km$cluster)

ggplot(pca_df, aes(PC1, PC2, color = Cluster, label = sample)) +
  geom_point(size = 2) +
  geom_text(vjust = 1.5, size = 2, show.legend = FALSE) +
  theme_minimal() +
  ggtitle(sprintf("PCA of %s-mer Frequencies with K-means Clusters, K=%s", nchar(colnames(kmer_matrix)[1]), k))

ggsave(sprintf("figures/K-means-clusters-PCA-of-%s-mers-K=%s.pdf", nchar(colnames(kmer_matrix)[1]), k), width = 8, height = 5)

