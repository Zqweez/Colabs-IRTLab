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
  ) %>%
  filter(!is.na(Cq))

# Calculate EHC reference values for normalization (EHC samples, excluding EHC-od)
ehc_reference <- genes_long %>%
  filter(grepl("^EHC-[12]", Sample)) %>%  # Only EHC-1 and EHC-2, not EHC-od
  group_by(Gene) %>%
  summarise(
    ehc_mean_cq = mean(Cq, na.rm = TRUE),
    .groups = "drop"
  )

cat("\nEHC reference values for normalization:\n")
print(ehc_reference)

# Calculate fold change in copy number using PCR efficiency
# Efficiency = 70% = 0.7, so each cycle multiplies by (1 + 0.7) = 1.7
pcr_efficiency <- 0.7
amplification_factor <- 1 + pcr_efficiency  # 1.7

genes_long_normalized <- genes_long %>%
  left_join(ehc_reference, by = "Gene") %>%
  mutate(
    # Calculate delta Cq (difference from EHC reference)
    delta_cq = Cq - ehc_mean_cq,
    # Calculate fold change: amplification_factor^(-delta_cq)
    # Negative delta_cq means lower Cq (more copies), so fold change > 1
    # Positive delta_cq means higher Cq (fewer copies), so fold change < 1
    fold_change = amplification_factor^(-delta_cq)
  )

# Calculate summary statistics for plotting (fold change values)
genes_summary <- genes_long_normalized %>%
  group_by(bio_rep, Combination, Gene) %>%
  summarise(
    mean_Cq_change = mean(Cq, na.rm = TRUE),
    sd_Cq_change = sd(Cq, na.rm = TRUE),
    n = n(),
    .groups = "drop"
  ) %>%
  # Handle cases where sd is NA (only one replicate)
  mutate(sd_Cq_change = ifelse(is.na(sd_Cq_change), 0, sd_Cq_change))

# Define colors for genes
gene_colors <- c('nidA' = '#e56562', 'phnAc' = '#51c558')

# Jitter position for data points
pos_jd <- position_jitterdodge(
  dodge.width = 0.8,
  jitter.width = 0.10
)

# First plot: Bio_rep level with both genes (fold change on log scale)
ggplot(genes_summary,
       aes(x = bio_rep, y = mean_Cq_change, fill = Gene)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  geom_point(data = genes_long_normalized,
             aes(x = bio_rep, y = Cq, fill = Gene),
             position = pos_jd,
             shape = 21,
             size = 0.7,
             stroke = 0.6,
             alpha = 0.7) +
  geom_errorbar(aes(ymin = pmax(mean_Cq_change - sd_Cq_change, 0.001),  # Prevent negative values on log scale
                    ymax = mean_Cq_change + sd_Cq_change),
                position = position_dodge(width = 0.8),
                width = 0.25) +
  facet_wrap(~ Combination, scales = "free_x", nrow = 1) +
  scale_fill_manual(values = gene_colors, name = "Gene") +
  labs(title = "Quantification Cycle for Biological Replicates",
       subtitle = "Technical replicates shown as points",
       x = "Biological Replicate",
       y = "Quantification Cycle (Cq)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
        strip.text = element_text(size = 8, angle = 0))

ggsave("results/growth_phe/qPCR/genes/Fold-change-per-biorep-nidA-phnAc.pdf", width = 12, height = 6)

## Do significance testing per gene to assess if there is a significant difference between the samples
# Test each gene independently for differences between Combinations using fold change values

cat("\nPerforming statistical analysis on fold change values...\n")

# Perform Wilcoxon tests for each gene separately using fold change
stat_test_nidA <- genes_long_normalized %>%
  filter(Gene == "nidA") %>%
  wilcox_test(fold_change ~ Combination, p.adjust.method = "BH") %>%
  add_significance("p.adj") %>%
  mutate(Gene = "nidA")

stat_test_phnAc <- genes_long_normalized %>%
  filter(Gene == "phnAc") %>%
  wilcox_test(fold_change ~ Combination, p.adjust.method = "BH") %>%
  add_significance("p.adj") %>%
  mutate(Gene = "phnAc")

# Combine results
stat_test_combined <- bind_rows(stat_test_nidA, stat_test_phnAc)

# Display results
cat("\nStatistical test results (fold change):\n")
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

cat("\nSignificance letters for nidA (fold change):\n")
print(cld_nidA)
cat("\nSignificance letters for phnAc (fold change):\n")
print(cld_phnAc)

# Calculate combination-level summary for plotting (fold change values)
combination_summary <- genes_long_normalized %>%
  group_by(Combination, Gene) %>%
  summarise(
    mean_fold_change = mean(fold_change, na.rm = TRUE),
    sd_fold_change = sd(fold_change, na.rm = TRUE),
    n = n(),
    .groups = "drop"
  ) %>%
  mutate(sd_fold_change = ifelse(is.na(sd_fold_change), 0, sd_fold_change))

# Add significance letters to summary
combination_summary_labeled <- combination_summary %>%
  left_join(cld_combined, by = c("Combination", "Gene")) %>%
  mutate(
    # Position letters above error bars, accounting for log scale
    letter_y = mean_fold_change + sd_fold_change + 0.2
  )

# Define colors for combinations (you can adjust these as needed)
combination_colors <- c('#e56562', '#51c558', '#67b4ef',   "#e6ab02", "#fed996", "#1b9e77", "#7570b3", "#fbd3e1",
  "#a6cee3", "#b2df8a"
)
# combination_colors <- c('KS8-od' = '#e56562', 'EHC-od' = '#51c558', 'EHC' = '#67b4ef', 
#                        'KS8-10' = "#e6ab02", 'KS8-100' = "#fed996", 'KS3-10' = "#1b9e77", 
#                        'KS3-100' = "#7570b3", 'KS3-KS8' = "#fbd3e1", 'NC' = "#a6cee3",
#                        'NC-samp' = "#ff7f00")

# Plot by combination with significance letters, faceted by gene (fold change on log scale)
y_max <- ceiling(max(combination_summary_labeled$mean_fold_change +
                       combination_summary_labeled$sd_fold_change)*1.2)
ggplot(combination_summary_labeled,
       aes(x = Combination, y = mean_fold_change, fill = Combination)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  geom_point(data = genes_long_normalized,
             aes(x = Combination, y = fold_change, fill = Combination),
             position = position_jitter(width = 0.2),
             shape = 21,
             size = 0.8,
             stroke = 0.6,
             alpha = 0.7) +
  geom_errorbar(aes(ymin = pmax(mean_fold_change - sd_fold_change, 0.001),  # Prevent negative values on log scale
                    ymax = mean_fold_change + sd_fold_change),
                position = position_dodge(width = 0.8),
                width = 0.25) +
  geom_text(aes(x = Combination, y = letter_y, label = letters),
            size = 3, fontface = "bold") +
  facet_wrap(~ Gene, nrow = 1) +
  scale_fill_manual(values = combination_colors, name = "Combination") +
  scale_y_continuous(
    breaks        = seq(0, y_max, by = 1),      # major grid every 1
    minor_breaks  = seq(0, y_max, by = 0.5)    # minor grid every 0.5
  ) +
  labs(title = "Fold Change in Gene Copy Number per Combination",
       subtitle = "Normalized to EHC controls, 70% PCR efficiency - Letters indicate significance (p<0.05)",
       x = "Combination",
       y = "Fold Change in Copy Number") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
        strip.text = element_text(size = 10, face = "bold"),
        legend.position = "none",
        panel.grid.major.x = element_blank(),                 # keep x grid clean
        panel.grid.major.y = element_line(colour = "grey40",  # darker major lines
                                          linewidth = 0.4),
        panel.grid.minor.y = element_line(colour = "grey70",  # extra horizontal lines
                                          linewidth = 0.25),
        panel.grid.minor.x = element_blank()) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "gray50", alpha = 0.7)

ggsave("results/growth_phe/qPCR/genes/nidA-Fold-change-per-combination-letters-genes.pdf", width = 12, height = 8)

cat("\nPlotting completed successfully!\n")
cat("Generated plots:\n")
cat("1. Fold-change-per-biorep-nidA-phnAc.pdf - Biological replicate level with fold change\n")
cat("2. Fold-change-per-combination-letters-genes.pdf - Combination level with significance testing\n")
cat("\nNormalization approach:\n")
cat("- Cq values normalized to EHC control mean (EHC-1 and EHC-2, excluding EHC-od)\n")
cat("- Fold change calculated using 70% PCR efficiency: (1.7)^(-ΔCq)\n")
cat("- ΔCq = Sample_Cq - EHC_mean_Cq\n")
cat("- Values > 1 indicate more copies than EHC, values < 1 indicate fewer copies\n")
cat("- Log scale used for visualization to accommodate wide range of fold changes\n")
cat("\nStatistical analysis summary:\n")
cat("- Performed pairwise Wilcoxon tests on fold change values for each gene\n")
cat("- Applied Benjamini-Hochberg correction for multiple comparisons\n")
cat("- Generated significance letters for each gene independently\n")
cat("- Dashed line at y=1 represents no change compared to EHC reference\n")
