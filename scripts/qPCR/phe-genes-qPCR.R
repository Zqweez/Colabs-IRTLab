## This script is to compare the nidA and phnAc genes in all samples
# Since two runs were made we first recalculate the Cq values to make them comparable
# The approach: use nidA run as reference to calculate thresholds, then apply to phnAc run
library(tidyverse)
library(readxl)
library(ggplot2)
library(rstatix)

# Load nidA data (reference run)
nidA_amp_raw <- read_excel("data/phe_growth/qPCR/nidA/Kasper_2025-07-30 14-15-02_BR005661 -  Quantification Amplification Results.xlsx", sheet = "SYBR")
nidA_cq_raw <- read_excel("data/phe_growth/qPCR/nidA/Kasper_2025-07-30 14-15-02_BR005661 -  Quantification Cq Results.xlsx", sheet = "cq")

# For phnAc data - use same files for now (will be updated when phnAc run is available)
# TODO: Update these paths when phnAc run is completed
phnAc_amp_raw <- read_excel("data/phe_growth/qPCR/phnAc/Kasper_2025-07-30 16-08-15_BR005661 -  Quantification Amplification Results.xlsx", sheet = "SYBR")
phnAc_cq_raw <- read_excel("data/phe_growth/qPCR/phnAc/Kasper_2025-07-30 16-08-15_BR005661 -  Quantification Cq Results.xlsx", sheet = "cq")

# Function to calculate threshold from amplification data and Cq
calculate_threshold <- function(amp_data, sample_col, cq_value) {
  # Get the amplification curve for this sample
  cycle_col <- amp_data[[1]]  # "Sample" column contains cycle numbers (starting from row 2)
  amp_values <- as.numeric(amp_data[[sample_col]])
  
  # Remove the first row which contains well identifiers
  cycles <- as.numeric(cycle_col[-1])
  fluorescence <- amp_values[-1]
  
  # Remove any NA values
  valid_idx <- !is.na(cycles) & !is.na(fluorescence)
  cycles <- cycles[valid_idx]
  fluorescence <- fluorescence[valid_idx]
  
  if(length(cycles) == 0 || length(fluorescence) == 0) {
    return(NA)
  }
  
  # Find the fluorescence value at the Cq
  # Linear interpolation between the cycle before and after Cq
  if(cq_value < min(cycles) || cq_value > max(cycles)) {
    return(NA)
  }
  
  # Find the cycles around the Cq value
  before_idx <- which(cycles <= cq_value)
  after_idx <- which(cycles > cq_value)
  
  if(length(before_idx) == 0 || length(after_idx) == 0) {
    return(NA)
  }
  
  before_cycle <- max(cycles[before_idx])
  after_cycle <- min(cycles[after_idx])
  before_fluor <- fluorescence[cycles == before_cycle][1]
  after_fluor <- fluorescence[cycles == after_cycle][1]
  
  # Linear interpolation
  threshold <- before_fluor + (after_fluor - before_fluor) * 
    (cq_value - before_cycle) / (after_cycle - before_cycle)
  
  return(threshold)
}
# Function to calculate Cq from amplification data given a threshold
calculate_cq_from_threshold <- function(amp_data, sample_col, threshold) {
  # Get the amplification curve for this sample
  cycle_col <- amp_data[[1]]  # "Sample" column contains cycle numbers (starting from row 2)
  amp_values <- as.numeric(amp_data[[sample_col]])
  
  # Remove the first row which contains well identifiers
  cycles <- as.numeric(cycle_col[-1])
  fluorescence <- amp_values[-1]
  
  # Remove any NA values
  valid_idx <- !is.na(cycles) & !is.na(fluorescence)
  cycles <- cycles[valid_idx]
  fluorescence <- fluorescence[valid_idx]
  
  if(length(cycles) == 0 || length(fluorescence) == 0) {
    return(NA)
  }
  
  # Find where fluorescence crosses the threshold
  # Look for the first point where fluorescence exceeds threshold
  cross_idx <- which(fluorescence >= threshold)
  
  if(length(cross_idx) == 0) {
    return(NA)  # Never crosses threshold
  }
  
  first_cross <- min(cross_idx)
  
  if(first_cross == 1) {
    return(cycles[1])  # Crosses at first cycle
  }
  
  # Linear interpolation between the cycle before and at crossing
  before_idx <- first_cross - 1
  after_idx <- first_cross
  
  before_cycle <- cycles[before_idx]
  after_cycle <- cycles[after_idx]
  before_fluor <- fluorescence[before_idx]
  after_fluor <- fluorescence[after_idx]
  
  # Linear interpolation to find exact crossing point
  cq <- before_cycle + (after_cycle - before_cycle) * 
    (threshold - before_fluor) / (after_fluor - before_fluor)
  
  return(cq)
}

# Function to create sample name to column mapping
create_sample_mapping <- function(cq_data, amp_data) {
  # Get sample names from Cq data and column names from amplification data
  sample_names <- cq_data$Sample
  amp_columns <- colnames(amp_data)[-1]  # Exclude the "Sample" column
  
  # Create mapping
  mapping <- data.frame(
    Sample = sample_names,
    Well = cq_data$Well,
    Column = NA,
    stringsAsFactors = FALSE
  )
  
  # Match samples to columns
  for(i in seq_len(nrow(mapping))) {
    sample_name <- mapping$Sample[i]
    
    # Clean sample name for matching (remove hyphens, standardize format)
    clean_sample <- gsub("-", "", sample_name)
    
    # Find matching column
    for(col in amp_columns) {
      clean_col <- gsub("-", "", col)
      if(clean_sample == clean_col) {
        mapping$Column[i] <- col
        break
      }
    }
  }
  
  return(mapping)
}

cat("Processing nidA data for threshold calculation...\n")

# Create sample mapping for nidA
nidA_mapping <- create_sample_mapping(nidA_cq_raw, nidA_amp_raw)

# Calculate thresholds from nidA reference samples
cat("Calculating thresholds from nidA samples...\n")

threshold_results <- data.frame(
  Sample = character(),
  Well = character(),
  Original_Cq = numeric(),
  Calculated_Threshold = numeric(),
  Column_Used = character(),
  stringsAsFactors = FALSE
)

for(i in seq_len(nrow(nidA_mapping))) {
  sample_name <- nidA_mapping$Sample[i]
  well <- nidA_mapping$Well[i]
  column_name <- nidA_mapping$Column[i]
  
  if(!is.na(column_name)) {
    # Get Cq value for this sample
    cq_value <- nidA_cq_raw$Cq[nidA_cq_raw$Sample == sample_name][1]
    
    if(!is.na(cq_value)) {
      threshold <- calculate_threshold(nidA_amp_raw, column_name, cq_value)
      
      threshold_results <- rbind(threshold_results, data.frame(
        Sample = sample_name,
        Well = well,
        Original_Cq = cq_value,
        Calculated_Threshold = threshold,
        Column_Used = column_name,
        stringsAsFactors = FALSE
      ))
    }
  } else {
    cat(sprintf("Warning: Could not find matching column for sample %s\n", sample_name))
  }
}

# Calculate final threshold statistics
valid_thresholds <- threshold_results$Calculated_Threshold[!is.na(threshold_results$Calculated_Threshold)]

if(length(valid_thresholds) > 0) {
  mean_threshold <- mean(valid_thresholds, na.rm = TRUE)
  median_threshold <- median(valid_thresholds, na.rm = TRUE)
  sd_threshold <- sd(valid_thresholds, na.rm = TRUE)
  
  cat(sprintf("\nThreshold statistics from %d valid samples:\n", length(valid_thresholds)))
  cat(sprintf("Mean threshold: %.3f\n", mean_threshold))
  cat(sprintf("Median threshold: %.3f\n", median_threshold))
  cat(sprintf("Standard deviation: %.3f\n", sd_threshold))
  cat(sprintf("Range: %.3f to %.3f\n", min(valid_thresholds), max(valid_thresholds)))
  
  # Use median threshold as final threshold (more robust to outliers)
  final_threshold <- median_threshold
  cat(sprintf("\nUsing median threshold: %.3f\n", final_threshold))
} else {
  cat("Error: No valid thresholds calculated!\n")
  stop("Cannot proceed without valid threshold calculations")
}

# Now apply the threshold to both nidA and phnAc samples to get corrected Cq values
cat("\nApplying threshold to all samples...\n")

# Create mappings for both datasets
nidA_mapping <- create_sample_mapping(nidA_cq_raw, nidA_amp_raw)
phnAc_mapping <- create_sample_mapping(phnAc_cq_raw, phnAc_amp_raw)

# Calculate corrected Cq values for nidA samples
corrected_nidA <- data.frame(
  Sample = character(),
  Well = character(),
  Original_Cq_nidA = numeric(),
  Corrected_Cq_nidA = numeric(),
  stringsAsFactors = FALSE
)

for(i in seq_len(nrow(nidA_mapping))) {
  sample_name <- nidA_mapping$Sample[i]
  well <- nidA_mapping$Well[i]
  column_name <- nidA_mapping$Column[i]
  
  if(!is.na(column_name)) {
    original_cq <- nidA_cq_raw$Cq[nidA_cq_raw$Sample == sample_name][1]
    corrected_cq <- calculate_cq_from_threshold(nidA_amp_raw, column_name, final_threshold)
    
    corrected_nidA <- rbind(corrected_nidA, data.frame(
      Sample = sample_name,
      Well = well,
      Original_Cq_nidA = original_cq,
      Corrected_Cq_nidA = corrected_cq,
      stringsAsFactors = FALSE
    ))
  }
}

# Calculate corrected Cq values for phnAc samples
corrected_phnAc <- data.frame(
  Sample = character(),
  Well = character(),
  Original_Cq_phnAc = numeric(),
  Corrected_Cq_phnAc = numeric(),
  stringsAsFactors = FALSE
)

for(i in seq_len(nrow(phnAc_mapping))) {
  sample_name <- phnAc_mapping$Sample[i]
  well <- phnAc_mapping$Well[i]
  column_name <- phnAc_mapping$Column[i]
  
  if(!is.na(column_name)) {
    original_cq <- phnAc_cq_raw$Cq[phnAc_cq_raw$Sample == sample_name][1]
    corrected_cq <- calculate_cq_from_threshold(phnAc_amp_raw, column_name, final_threshold)
    
    corrected_phnAc <- rbind(corrected_phnAc, data.frame(
      Sample = sample_name,
      Well = well,
      Original_Cq_phnAc = original_cq,
      Corrected_Cq_phnAc = corrected_cq,
      stringsAsFactors = FALSE
    ))
  }
}

# Merge the results to create the final CSV with all 5 columns
final_results <- merge(corrected_nidA, corrected_phnAc, by = "Sample", all = TRUE, suffixes = c("_nidA", "_phnAc"))

# Clean up the merged data and select only the required columns
final_output <- final_results %>%
  select(Sample, Original_Cq_nidA, Original_Cq_phnAc, Corrected_Cq_nidA, Corrected_Cq_phnAc) %>%
  arrange(Sample)

# Display summary of results
cat("\nSummary of final results:\n")
cat(sprintf("Total samples processed: %d\n", nrow(final_output)))
cat(sprintf("Samples with nidA data: %d\n", sum(!is.na(final_output$Original_Cq_nidA))))
cat(sprintf("Samples with phnAc data: %d\n", sum(!is.na(final_output$Original_Cq_phnAc))))

# Show first few rows
cat("\nFirst 10 rows of final results:\n")
print(head(final_output, 10))

# Calculate and display correction statistics
nidA_corrections <- final_output$Corrected_Cq_nidA - final_output$Original_Cq_nidA
phnAc_corrections <- final_output$Corrected_Cq_phnAc - final_output$Original_Cq_phnAc

cat("\nnidA Cq correction statistics:\n")
valid_nidA_corrections <- nidA_corrections[!is.na(nidA_corrections)]
if(length(valid_nidA_corrections) > 0) {
  cat(sprintf("Mean correction: %.3f\n", mean(valid_nidA_corrections)))
  cat(sprintf("Range: %.3f to %.3f\n", min(valid_nidA_corrections), max(valid_nidA_corrections)))
}

cat("\nphnAc Cq correction statistics:\n")
valid_phnAc_corrections <- phnAc_corrections[!is.na(phnAc_corrections)]
if(length(valid_phnAc_corrections) > 0) {
  cat(sprintf("Mean correction: %.3f\n", mean(valid_phnAc_corrections)))
  cat(sprintf("Range: %.3f to %.3f\n", min(valid_phnAc_corrections), max(valid_phnAc_corrections)))
}

# Save all results
write.csv(threshold_results, "results/growth_phe/qPCR/genes/qPCR_nidA_thresholds.csv", row.names = FALSE)
write.csv(final_output, "results/growth_phe/qPCR/genes/qPCR_corrected_Cq_values.csv", row.names = FALSE)

cat(sprintf("\nResults saved to:\n"))
cat("- results/growth_phe/qPCR/genes/qPCR_nidA_thresholds.csv (threshold calculations)\n")
cat("- results/growth_phe/qPCR/genes/qPCR_corrected_Cq_values.csv (final 5-column output)\n")

cat("\nScript completed successfully!\n")
cat("Note: Update phnAc file paths when the phnAc run data becomes available.\n")


## Now make plots to show the Cq values in a chart but first filter out the outliers.

# Prepare data for plotting - create proper structure for both genes
genes_plot_data <- final_output %>%
  filter(!(Sample %in% c("NC-2-2"))) %>%
  mutate(
    bio_rep = sub("-[123]$", "", Sample),
    Combination = case_when(
      grepl("^NC-[12]", Sample) ~ "NC",  # Only NC-1 and NC-2, not NC-samp
      TRUE ~ sub("-[123]$", "", bio_rep)
    )
  )

# Create long format data for both genes
genes_long <- genes_plot_data %>%
  pivot_longer(
    cols = c(Corrected_Cq_nidA, Corrected_Cq_phnAc),
    names_to = "Gene",
    values_to = "Cq"
  ) %>%
  mutate(
    Gene = case_when(
      Gene == "Corrected_Cq_nidA" ~ "nidA",
      Gene == "Corrected_Cq_phnAc" ~ "phnAc",
      TRUE ~ Gene
    )
  )

# Calculate NC reference values for normalization (only NC-1 and NC-2, not NC-samp)
nc_reference <- genes_long %>%
  filter(grepl("^NC-[12]", Sample)) %>%  # Only NC-1 and NC-2
  group_by(Gene) %>%
  summarise(
    nc_mean_cq = mean(Cq, na.rm = TRUE),
    .groups = "drop"
  )

cat("\nNC reference values for normalization:\n")
print(nc_reference)

# Normalize Cq values and calculate relative abundance
genes_long_normalized <- genes_long %>%
  left_join(nc_reference, by = "Gene") %>%
  mutate(
    # Normalize by dividing by NC, then take reciprocal for relative abundance
    # Higher values = higher abundance
    relative_abundance = 1 / (Cq / nc_mean_cq)
  )

# Calculate summary statistics for plotting (normalized values)
genes_summary <- genes_long_normalized %>%
  group_by(bio_rep, Combination, Gene) %>%
  summarise(
    mean_abundance = mean(relative_abundance, na.rm = TRUE),
    sd_abundance = sd(relative_abundance, na.rm = TRUE),
    n = n(),
    .groups = "drop"
  ) %>%
  # Handle cases where sd is NA (only one replicate)
  mutate(sd_abundance = ifelse(is.na(sd_abundance), 0, sd_abundance))

# Define colors for genes
gene_colors <- c('nidA' = '#e56562', 'phnAc' = '#51c558')

# Jitter position for data points
pos_jd <- position_jitterdodge(
  dodge.width = 0.8,
  jitter.width = 0.10
)

# First plot: Bio_rep level with both genes (normalized)
ggplot(genes_summary,
       aes(x = bio_rep, y = mean_abundance, fill = Gene)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  geom_point(data = genes_long_normalized,
             aes(x = bio_rep, y = relative_abundance, fill = Gene),
             position = pos_jd,
             shape = 21,
             size = 0.7,
             stroke = 0.6,
             alpha = 0.7) +
  geom_errorbar(aes(ymin = mean_abundance - sd_abundance,
                    ymax = mean_abundance + sd_abundance),
                position = position_dodge(width = 0.8),
                width = 0.25) +
  facet_wrap(~ Combination, scales = "free_x", nrow = 1) +
  scale_fill_manual(values = gene_colors, name = "Gene") +
  labs(title = "Relative Gene Abundance for Biological Replicates",
       subtitle = "Normalized to NC controls (higher bars = higher abundance)",
       x = "Biological Replicate",
       y = "Relative Abundance (normalized to NC)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
        strip.text = element_text(size = 8, angle = 0))

ggsave("results/growth_phe/qPCR/genes/Relative-abundance-per-biorep-nidA-phnAc.pdf", width = 12, height = 6)

## Do significance testing per gene to assess if there is a significant difference between the samples
# Test each gene independently for differences between Combinations using normalized abundance values

cat("\nPerforming statistical analysis on normalized abundance values...\n")

# Perform Wilcoxon tests for each gene separately using relative abundance
stat_test_nidA <- genes_long_normalized %>%
  filter(Gene == "nidA") %>%
  wilcox_test(relative_abundance ~ Combination, p.adjust.method = "BH") %>%
  add_significance("p.adj") %>%
  mutate(Gene = "nidA")

stat_test_phnAc <- genes_long_normalized %>%
  filter(Gene == "phnAc") %>%
  wilcox_test(relative_abundance ~ Combination, p.adjust.method = "BH") %>%
  add_significance("p.adj") %>%
  mutate(Gene = "phnAc")

# Combine results
stat_test_combined <- bind_rows(stat_test_nidA, stat_test_phnAc)

# Display results
cat("\nStatistical test results (relative abundance):\n")
print(stat_test_combined %>% select(Gene, group1, group2, p, p.adj, p.adj.signif))

# Custom function to generate significance letters
generate_significance_letters <- function(stat_results, alpha = 0.05) {
  # Get unique groups
  groups <- unique(c(stat_results$group1, stat_results$group2))
  groups <- sort(groups)
  
  # Create adjacency matrix for significant comparisons
  n_groups <- length(groups)
  sig_matrix <- matrix(FALSE, nrow = n_groups, ncol = n_groups)
  rownames(sig_matrix) <- groups
  colnames(sig_matrix) <- groups
  
  # Fill matrix with significant comparisons (p.adj < alpha)
  for(i in seq_len(nrow(stat_results))) {
    if(stat_results$p.adj[i] < alpha) {
      g1 <- stat_results$group1[i]
      g2 <- stat_results$group2[i]
      sig_matrix[g1, g2] <- TRUE
      sig_matrix[g2, g1] <- TRUE
    }
  }
  
  # Generate letters using a simple algorithm
  letters_assigned <- character(n_groups)
  names(letters_assigned) <- groups
  current_letter <- 1
  
  for(i in seq_len(n_groups)) {
    if(letters_assigned[i] == "") {
      # Find all groups that are NOT significantly different from this one
      same_group <- c(i)
      for(j in (i+1):n_groups) {
        if(j <= n_groups && !sig_matrix[i, j]) {
          same_group <- c(same_group, j)
        }
      }
      
      # Assign the same letter to all groups in this set
      letter <- letters[current_letter]
      for(idx in same_group) {
        if(letters_assigned[idx] == "") {
          letters_assigned[idx] <- letter
        } else {
          # If already has a letter, combine them
          if (!grepl(letter, letters_assigned[idx])) {
            letters_assigned[idx] <- paste0(letters_assigned[idx], letter)
          }
        }
      }
      current_letter <- current_letter + 1
    }
  }
  
  return(data.frame(
    Combination = names(letters_assigned),
    letters = letters_assigned,
    stringsAsFactors = FALSE
  ))
}

# Generate significance letters for each gene
cld_nidA <- generate_significance_letters(stat_test_nidA, alpha = 0.05) %>%
  mutate(Gene = "nidA")
cld_phnAc <- generate_significance_letters(stat_test_phnAc, alpha = 0.05) %>%
  mutate(Gene = "phnAc")

# Combine letters
cld_combined <- bind_rows(cld_nidA, cld_phnAc)

cat("\nSignificance letters for nidA (relative abundance):\n")
print(cld_nidA)
cat("\nSignificance letters for phnAc (relative abundance):\n")
print(cld_phnAc)

# Calculate combination-level summary for plotting (normalized values)
combination_summary <- genes_long_normalized %>%
  group_by(Combination, Gene) %>%
  summarise(
    mean_abundance = mean(relative_abundance, na.rm = TRUE),
    sd_abundance = sd(relative_abundance, na.rm = TRUE),
    n = n(),
    .groups = "drop"
  ) %>%
  mutate(sd_abundance = ifelse(is.na(sd_abundance), 0, sd_abundance))

# Add significance letters to summary
combination_summary_labeled <- combination_summary %>%
  left_join(cld_combined, by = c("Combination", "Gene")) %>%
  mutate(letter_y = mean_abundance + sd_abundance + 0.1)  # Position letters above error bars

# Define colors for combinations (you can adjust these as needed)
combination_colors <- c('KS8-od' = '#e56562', 'EHC-od' = '#51c558', 'EHC' = '#67b4ef', 
                       'KS8-10' = "#e6ab02", 'KS8-100' = "#fed996", 'KS3-10' = "#1b9e77", 
                       'KS3-100' = "#7570b3", 'KS3-KS8' = "#fbd3e1", 'NC' = "#a6cee3",
                       'NC-samp' = "#ff7f00")

# Plot by combination with significance letters, faceted by gene (normalized values)
ggplot(combination_summary_labeled,
       aes(x = Combination, y = mean_abundance, fill = Combination)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  geom_point(data = genes_long_normalized,
             aes(x = Combination, y = relative_abundance, fill = Combination),
             position = position_jitter(width = 0.2),
             shape = 21,
             size = 0.8,
             stroke = 0.6,
             alpha = 0.7) +
  geom_errorbar(aes(ymin = mean_abundance - sd_abundance,
                    ymax = mean_abundance + sd_abundance),
                position = position_dodge(width = 0.8),
                width = 0.25) +
  geom_text(aes(x = Combination, y = letter_y, label = letters),
            size = 3, fontface = "bold") +
  facet_wrap(~ Gene, scales = "free_y", nrow = 1) +
  scale_fill_manual(values = combination_colors, name = "Combination") +
  labs(title = "Relative Gene Abundance per Combination",
       subtitle = "Normalized to NC controls - Letters indicate statistical significance (p<0.05)",
       x = "Combination",
       y = "Relative Abundance (normalized to NC)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
        strip.text = element_text(size = 10, face = "bold"),
        legend.position = "bottom")

ggsave("results/growth_phe/qPCR/genes/Relative-abundance-per-combination-letters-genes.pdf", width = 12, height = 6)

cat("\nPlotting completed successfully!\n")
cat("Generated plots:\n")
cat("1. Relative-abundance-per-biorep-nidA-phnAc.pdf - Biological replicate level with normalized abundance\n")
cat("2. Relative-abundance-per-combination-letters-genes.pdf - Combination level with significance testing\n")
cat("\nNormalization approach:\n")
cat("- Cq values normalized by dividing by NC control mean (NC-1 and NC-2 only)\n")
cat("- Relative abundance calculated as 1/(normalized Cq) - higher bars = higher abundance\n")
cat("- NC controls excluded from normalization reference (only NC-1 and NC-2 used)\n")
cat("\nStatistical analysis summary:\n")
cat("- Performed pairwise Wilcoxon tests on normalized abundance values for each gene\n")
cat("- Applied Benjamini-Hochberg correction for multiple comparisons\n")
cat("- Generated significance letters for each gene independently\n")
cat("- Higher relative abundance values indicate higher gene expression\n")
