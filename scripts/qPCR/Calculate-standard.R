library(tidyverse)
library(readxl)
library(ggplot2)
library(gridExtra)  # For combining plots
library(rstatix)    # For statistical analysis

## Load the sample data both Cq and Amplification
sample_data_raw <- read_excel("data/phe_growth/qPCR/Kasper_2025-07-28 14-45-16_BR005661 -  Quantification Amplification Results.xlsx", sheet = "SYBR")
sample_cq_raw <- read_excel("data/phe_growth/qPCR/Kasper_2025-07-28 14-45-16_BR005661 -  Quantification Cq Results.xlsx", sheet = "0")

## Load the Standard data both Cq and Amplification
std_data_raw <- read_excel("data/phe_growth/qPCR/Standard-curve.xlsx", sheet = "amplification")
std_cq_raw <- read_excel("data/phe_growth/qPCR/Standard-curve.xlsx", sheet = "cq")

# Function to calculate threshold from amplification data and Cq
calculate_threshold <- function(amp_data, sample_col, cq_value) {
  # Get the amplification curve for this sample
  if(sample_col == "copies/ml") {
    # For standard data, the first column is "copies/ml" containing cycle info
    cycle_col <- amp_data[[1]]
    amp_values <- as.numeric(amp_data[[sample_col]])
  } else {
    # For sample data, the first column is "Sample" containing cycle info
    cycle_col <- amp_data[[1]]
    amp_values <- as.numeric(amp_data[[sample_col]])
  }
  
  # Remove the first row which contains headers/sample names
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
  if(sample_col == "copies/ml") {
    # For standard data, the first column is "copies/ml" containing cycle info
    cycle_col <- amp_data[[1]]
    amp_values <- as.numeric(amp_data[[sample_col]])
  } else {
    # For sample data, the first column is "Sample" containing cycle info
    cycle_col <- amp_data[[1]]
    amp_values <- as.numeric(amp_data[[sample_col]])
  }
  
  # Remove the first row which contains headers/sample names
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

# Calculate thresholds for all samples
cat("Calculating thresholds for sample data...\n")

# Prepare sample data for analysis
sample_cq_clean <- sample_cq_raw %>%
  filter(!is.na(Cq)) %>%
  select(Sample, Cq)

# Get all sample column names from amplification data (excluding "Sample" column)
sample_columns <- colnames(sample_data_raw)[-1]

# Calculate threshold for each sample
threshold_results <- data.frame(
  Sample = character(),
  Original_Cq = numeric(),
  Calculated_Threshold = numeric(),
  Column_Matched = character(),
  stringsAsFactors = FALSE
)

for(i in seq_len(nrow(sample_cq_clean))) {
  sample_name <- sample_cq_clean$Sample[i]
  cq_value <- sample_cq_clean$Cq[i]
  
  # Find corresponding column in amplification data
  # Clean sample name for matching (remove hyphens and standardize)
  clean_sample_name <- gsub("-", "", sample_name)
  
  # Find matching column by checking if sample name appears in column name
  matching_col <- NULL
  for(col_name in sample_columns) {
    clean_col_name <- gsub("-", "", col_name)
    if(clean_sample_name == clean_col_name) {
      matching_col <- col_name
      break
    }
  }
  
  if(!is.null(matching_col)) {
    threshold <- calculate_threshold(sample_data_raw, matching_col, cq_value)
    
    threshold_results <- rbind(threshold_results, data.frame(
      Sample = sample_name,
      Original_Cq = cq_value,
      Calculated_Threshold = threshold,
      Column_Matched = matching_col,
      stringsAsFactors = FALSE
    ))
  } else {
    cat(sprintf("Warning: Could not find matching column for sample %s\n", sample_name))
  }
}

# Display threshold results with more detail
cat("\nCalculated thresholds for samples:\n")
print(head(threshold_results, 10))
cat(sprintf("... and %d more samples\n", nrow(threshold_results) - 10))

# Show threshold distribution
cat("\nThreshold distribution:\n")
valid_thresholds <- threshold_results$Calculated_Threshold[!is.na(threshold_results$Calculated_Threshold)]
cat(sprintf("Min: %.3f, Max: %.3f, Range: %.3f\n", 
           min(valid_thresholds), max(valid_thresholds), 
           max(valid_thresholds) - min(valid_thresholds)))

# Calculate mean threshold from samples (excluding outliers)
mean_threshold <- mean(valid_thresholds, na.rm = TRUE)
median_threshold <- median(valid_thresholds, na.rm = TRUE)
sd_threshold <- sd(valid_thresholds, na.rm = TRUE)

cat(sprintf("\nThreshold statistics:\n"))
cat(sprintf("Mean threshold: %.3f\n", mean_threshold))
cat(sprintf("Median threshold: %.3f\n", median_threshold))
cat(sprintf("Standard deviation: %.3f\n", sd_threshold))

# Remove outliers (beyond 2 standard deviations)
outlier_threshold <- 2
clean_thresholds <- valid_thresholds[abs(valid_thresholds - mean_threshold) <= outlier_threshold * sd_threshold]
final_threshold <- mean(clean_thresholds)

cat(sprintf("Final threshold (outliers removed): %.3f\n", final_threshold))
cat(sprintf("Number of samples used: %d out of %d\n", length(clean_thresholds), length(valid_thresholds)))

# Now apply this threshold to the standard curve data
cat("\nApplying threshold to standard curve data...\n")

# Create a mapping between wells and columns in standard amplification data
# The standard data has columns named by concentration, but we need to map wells to concentrations
std_columns <- colnames(std_data_raw)[-1]  # Remove "copies/ml" column

cat("Available standard columns:", paste(std_columns, collapse = ", "), "\n")
cat("Standard curve wells and concentrations:\n")
print(std_cq_raw[, c("Well", "copies/ml")])

# Create a manual mapping based on the data structure we observed
# From the output, we can see the mapping of wells to columns
well_to_column_map <- list(
  "H11" = "100000000...2",    # 100,000,000 copies/ml
  "A12" = "10000000",         # 10,000,000 copies/ml  
  "B12" = "1000000",          # 1,000,000 copies/ml
  "C12" = "100000...5",       # 100,000 copies/ml
  "D12" = "10000",            # 10,000 copies/ml
  "E12" = "1000...7",         # 1,000 copies/ml
  "F12" = "100000000...8",    # 100,000,000 copies/ml (duplicate)
  "G12" = "100000...9",       # 100,000 copies/ml (duplicate)
  "H12" = "1000...10"         # 1,000 copies/ml (duplicate)
)

# Calculate corrected Cq values for standards
corrected_std_results <- data.frame(
  Sample = character(),
  Well = character(),
  Copies_per_ml = numeric(),
  Log10_copies = numeric(),
  Original_Cq = numeric(),
  Corrected_Cq = numeric(),
  Threshold_Used = numeric(),
  Column_Used = character(),
  stringsAsFactors = FALSE
)

for(i in seq_len(nrow(std_cq_raw))) {
  well <- std_cq_raw$Well[i]
  copies_ml <- std_cq_raw$`copies/ml`[i]
  log10_copies <- std_cq_raw$log10[i]
  original_cq <- std_cq_raw$cq[i]
  
  # Find corresponding column using our mapping
  matching_col <- well_to_column_map[[well]]
  
  if(!is.null(matching_col) && matching_col %in% std_columns) {
    cat(sprintf("Processing well %s (%g copies/ml) using column %s\n", well, copies_ml, matching_col))
    
    corrected_cq <- calculate_cq_from_threshold(std_data_raw, matching_col, final_threshold)
    
    corrected_std_results <- rbind(corrected_std_results, data.frame(
      Sample = paste0(well, " (", format(copies_ml, scientific = FALSE), " copies/ml)"),
      Well = well,
      Copies_per_ml = copies_ml,
      Log10_copies = log10_copies,
      Original_Cq = original_cq,
      Corrected_Cq = corrected_cq,
      Threshold_Used = final_threshold,
      Column_Used = matching_col,
      stringsAsFactors = FALSE
    ))
  } else {
    cat(sprintf("Warning: Could not find matching column for well %s\n", well))
  }
}

# Display corrected standard results
cat("\nCorrected Cq values for standard curve:\n")
print(corrected_std_results)

# Calculate the difference between original and corrected Cq values
corrected_std_results$Cq_Difference <- corrected_std_results$Corrected_Cq - corrected_std_results$Original_Cq

# Filter out NA values for summary statistics
valid_corrections <- corrected_std_results[!is.na(corrected_std_results$Cq_Difference), ]

cat("\nSummary of Cq corrections:\n")
if(nrow(valid_corrections) > 0) {
  cat(sprintf("Number of successful corrections: %d out of %d\n", nrow(valid_corrections), nrow(corrected_std_results)))
  cat(sprintf("Mean Cq difference: %.3f\n", mean(valid_corrections$Cq_Difference)))
  cat(sprintf("Range of differences: %.3f to %.3f\n", 
             min(valid_corrections$Cq_Difference),
             max(valid_corrections$Cq_Difference)))
} else {
  cat("No successful corrections were made.\n")
}

# Display any failed calculations
failed_corrections <- corrected_std_results[is.na(corrected_std_results$Cq_Difference), ]
if(nrow(failed_corrections) > 0) {
  cat("\nFailed corrections:\n")
  print(failed_corrections[, c("Well", "Copies_per_ml", "Original_Cq", "Column_Used")])
}

# Save results to CSV files for review
write.csv(threshold_results, "results/growth_phe/qPCR/qPCR_sample_thresholds.csv", row.names = FALSE)
write.csv(corrected_std_results, "results/growth_phe/qPCR/qPCR_corrected_standards.csv", row.names = FALSE)

cat("\nResults saved to:\n")
cat("- results/growth_phe/qPCR/qPCR_sample_thresholds.csv\n")
cat("- results/growth_phe/qPCR/qPCR_corrected_standards.csv\n")

# Final validation: Show the relationship between log concentration and Cq values
cat("\nStandard curve validation:\n")
cat("Log10(copies/ml) vs Original Cq vs Corrected Cq:\n")
comparison_table <- corrected_std_results[, c("Log10_copies", "Original_Cq", "Corrected_Cq", "Cq_Difference")]
comparison_table <- comparison_table[order(comparison_table$Log10_copies, decreasing = TRUE), ]
print(comparison_table)

# Calculate R-squared for both original and corrected standard curves
if(nrow(valid_corrections) >= 3) {
  original_lm <- lm(Original_Cq ~ Log10_copies, data = valid_corrections)
  corrected_lm <- lm(Corrected_Cq ~ Log10_copies, data = valid_corrections)
  
  cat(sprintf("\nStandard curve quality:\n"))
  cat(sprintf("Original curve R² = %.4f\n", summary(original_lm)$r.squared))
  cat(sprintf("Corrected curve R² = %.4f\n", summary(corrected_lm)$r.squared))
  cat(sprintf("Original slope = %.3f\n", coef(original_lm)[2]))
  cat(sprintf("Corrected slope = %.3f\n", coef(corrected_lm)[2]))
  
  # Ideal qPCR slope should be around -3.32 (100% efficiency)
  original_efficiency <- (10^(-1/coef(original_lm)[2]) - 1) * 100
  corrected_efficiency <- (10^(-1/coef(corrected_lm)[2]) - 1) * 100
  
  cat(sprintf("Original PCR efficiency = %.1f%%\n", original_efficiency))
  cat(sprintf("Corrected PCR efficiency = %.1f%%\n", corrected_efficiency))
}




# Create plots for both original and corrected standard curves
cat("\nGenerating standard curve plots...\n")

# Prepare data for plotting (remove duplicates by averaging)
plot_data <- valid_corrections %>%
  group_by(Log10_copies) %>%
  summarise(
    Original_Cq_mean = mean(Original_Cq, na.rm = TRUE),
    Corrected_Cq_mean = mean(Corrected_Cq, na.rm = TRUE),
    Original_Cq_sd = sd(Original_Cq, na.rm = TRUE),
    Corrected_Cq_sd = sd(Corrected_Cq, na.rm = TRUE),
    n_replicates = n(),
    .groups = "drop"
  ) %>%
  # Replace NAs in sd with 0 for single measurements
  mutate(
    Original_Cq_sd = ifelse(is.na(Original_Cq_sd), 0, Original_Cq_sd),
    Corrected_Cq_sd = ifelse(is.na(Corrected_Cq_sd), 0, Corrected_Cq_sd)
  )

# Create linear models for plotting
original_lm_plot <- lm(Original_Cq_mean ~ Log10_copies, data = plot_data)
corrected_lm_plot <- lm(Corrected_Cq_mean ~ Log10_copies, data = plot_data)

# Calculate statistics for annotation
original_stats <- list(
  r2 = summary(original_lm_plot)$r.squared,
  slope = coef(original_lm_plot)[2],
  efficiency = (10^(-1/coef(original_lm_plot)[2]) - 1) * 100
)

corrected_stats <- list(
  r2 = summary(corrected_lm_plot)$r.squared,
  slope = coef(corrected_lm_plot)[2],
  efficiency = (10^(-1/coef(corrected_lm_plot)[2]) - 1) * 100
)

# Prepare combined data for plotting
combined_data <- plot_data %>%
  pivot_longer(cols = c(Original_Cq_mean, Corrected_Cq_mean), 
               names_to = "Curve_Type", values_to = "Cq_Value") %>%
  mutate(Curve_Type = ifelse(Curve_Type == "Original_Cq_mean", "Original", "Corrected"))

# Add error bars separately for each curve type
combined_data <- combined_data %>%
  left_join(
    plot_data %>% 
      select(Log10_copies, Original_Cq_sd, Corrected_Cq_sd) %>%
      pivot_longer(cols = c(Original_Cq_sd, Corrected_Cq_sd),
                   names_to = "sd_type", values_to = "Cq_sd") %>%
      mutate(Curve_Type = ifelse(sd_type == "Original_Cq_sd", "Original", "Corrected")) %>%
      select(-sd_type),
    by = c("Log10_copies", "Curve_Type")
  )

# Create legend labels with statistics
legend_labels <- c(
  sprintf("Original (R² = %.4f, Slope = %.3f, Eff = %.1f%%)", 
          original_stats$r2, original_stats$slope, original_stats$efficiency),
  sprintf("Corrected (R² = %.4f, Slope = %.3f, Eff = %.1f%%)", 
          corrected_stats$r2, corrected_stats$slope, corrected_stats$efficiency)
)

# Create the combined plot with statistics in legend
p_combined <- ggplot(combined_data, aes(x = Log10_copies, y = Cq_Value, color = Curve_Type)) +
  geom_point(size = 3, alpha = 0.7) +
  geom_errorbar(aes(ymin = Cq_Value - Cq_sd, ymax = Cq_Value + Cq_sd), 
                width = 0.1, alpha = 1) +
  geom_line(stat="smooth", method = "lm", alpha = 0.4, linewidth = 1, linetype ="dashed") +
  scale_color_manual(values = c("Original" = "blue", "Corrected" = "red"), 
                     labels = legend_labels) +
  guides(color = guide_legend(
    override.aes = list(size = 1.2),   # icon / key size
    nrow = 2 
  )) +
  labs(
    title = "qPCR Standard Curve Comparison",
    subtitle = "Original vs Corrected Cq Values",
    x = expression(Log[10](copies/mL)),
    y = "Cq Value",
    color = "Curve Type"
  ) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    legend.text = element_text(size = 10)
  )

# Display the plot
print(p_combined)

# Save the combined plot
ggsave("results/growth_phe/qPCR/standard_curve_comparison.pdf", p_combined, width = 12, height = 8)

# Calculate sample concentrations using both curves
cat("\nCalculating sample concentrations...\n")

# Prepare sample data for quantification
sample_quantification <- sample_cq_clean %>%
  mutate(
    # Calculate copies/mL using original standard curve
    Log10_copies_original = (Cq - coef(original_lm)[1]) / coef(original_lm)[2],
    Copies_per_mL_original = 10^Log10_copies_original,
    
    # Calculate copies/mL using corrected standard curve  
    Log10_copies_corrected = (Cq - coef(corrected_lm)[1]) / coef(corrected_lm)[2],
    Copies_per_mL_corrected = 10^Log10_copies_corrected,
    
    # Calculate the difference
    Fold_difference = Copies_per_mL_corrected / Copies_per_mL_original,
    Log2_fold_difference = log2(Fold_difference)
  ) %>%
  # Add threshold information
  left_join(threshold_results[, c("Sample", "Calculated_Threshold", "Column_Matched")], 
            by = "Sample") %>%
  # Rename for clarity
  rename(
    Sample_Name = Sample,
    Original_Cq = Cq,
    Threshold_Used = Calculated_Threshold
  ) %>%
  # Select and order columns
  select(
    Sample_Name,
    Original_Cq,
    Threshold_Used,
    Column_Matched,
    Log10_copies_original,
    Copies_per_mL_original,
    Log10_copies_corrected,
    Copies_per_mL_corrected,
    Fold_difference,
    Log2_fold_difference
  )

sample_quant <- sample_quantification %>% 
  select(Sample_Name, Original_Cq, Log10_copies_original, Log10_copies_corrected) %>%
  filter(!(Sample_Name %in% c("EHC-1-1", "KS8-10-2-3", "KS3-10-2-3", "EHC-od-2-3"))) %>%
  mutate(
    bio_rep = sub("-[123]$", "", Sample_Name),
    Combination = sub("-[123]$", "", bio_rep),
  ) 

sample_quant_summary <- sample_quant %>%
  group_by(bio_rep, Combination) %>% # bio_rep,
  summarise(
    mean_log10_ori = mean(Log10_copies_original, na.rm = TRUE),
    sd_log10_ori = sd(Log10_copies_original, na.rm = TRUE),
    mean_log10_corr = mean(Log10_copies_corrected, na.rm = TRUE),
    sd_log10_corr = sd(Log10_copies_corrected, na.rm = TRUE),
    .groups = "drop"
  )

# Jitter for datapoints
pos_jd <- position_jitterdodge(
  dodge.width  = 0.8,  # same width as for the bars
  jitter.width = 0.10  # sideways spread
)
Cq_colors = c('#e56562', '#51c558', '#67b4ef',   "#e6ab02", "#fed996", "#1b9e77", "#7570b3", "#fbd3e1",
              "#a6cee3", "#b2df8a", '#d95f02', '#f781be'
)
# Barplot
ggplot(sample_quant_summary,
       aes(x = bio_rep, y = mean_log10_corr,
           fill = factor(Combination),    
           group = factor(Combination))) + # makes sure bars & points share groups
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  geom_point(data = sample_quant,
             aes(x = bio_rep, y = Log10_copies_corrected,
                 fill  = factor(Combination),   # same mapping as the bars
                 group = factor(Combination)),
             position = pos_jd,
             shape = 21,                 # hollow circle whose fill we control
             size  = 0.7,
             stroke = 0.6,
             alpha = 0.7) +
  geom_errorbar(aes(ymin = mean_log10_corr - sd_log10_corr,
                    ymax = mean_log10_corr + sd_log10_corr),
                position = position_dodge(width = 0.8),
                width = 0.25) +
  facet_wrap(~ Combination, scales = "free_x", nrow = 1) +
  scale_fill_manual(values = Cq_colors, name = "Sample") +
  labs(title = "Cells/mL for Biological Duplicates",
       subtitle = "Mean of technical replicates with SD as error bars",
       x = "Biological Duplicate",
       y = expression(Log[10](cells/mL))) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
        strip.text = element_text(size = 8, angle = 0))

ggsave("results/growth_phe/qPCR/cells-ml-samples.pdf", width = 8, height = 5)

# Add statistical analysis to compare combinations using Wilcoxon tests
cat("\nPerforming statistical analysis on log10 corrected values...\n")

# Load rstatix for easier statistical testing (if not already loaded)
if (!require(rstatix, quietly = TRUE)) {
  install.packages("rstatix")
  library(rstatix)
}

# Prepare data for statistical analysis - using corrected values
sample_quant_stats <- sample_quant %>%
  select(Sample_Name, bio_rep, Combination, Log10_copies_corrected) %>%
  filter(!is.na(Log10_copies_corrected))

# Perform Kruskal-Wallis test first
# H₀: All groups come from the same distribution
# H₁: At least one group differs in median
kw_test <- sample_quant_stats %>% 
  kruskal_test(Log10_copies_corrected ~ Combination)
cat("Kruskal-Wallis test results:\n")
print(kw_test)

# Post-hoc pairwise comparisons using Wilcoxon tests
stat_test <- sample_quant_stats %>%
  wilcox_test(Log10_copies_corrected ~ Combination, p.adjust.method = "BH") %>%
  add_significance("p.adj")

cat("\nPairwise Wilcoxon test results:\n")
print(stat_test)

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

# Generate significance letters
cld_df <- generate_significance_letters(stat_test, alpha = 0.05)
cat("\nSignificance letters:\n")
print(cld_df)

# Create summary with letters for the updated plot
sample_quant_summary_labeled <- sample_quant %>%
  group_by(Combination) %>% # bio_rep,
  summarise(
    mean_log10_corr = mean(Log10_copies_corrected, na.rm = TRUE),
    sd_log10_corr = sd(Log10_copies_corrected, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  left_join(cld_df, by = "Combination") %>%
  mutate(letter_y = mean_log10_corr + sd_log10_corr + 0.3)  # Position letters above error bars

# Also create a simpler plot showing combinations without faceting for easier comparison
y_max <- ceiling(max(combination_summary_labeled$mean_fold_change +
                       combination_summary_labeled$sd_fold_change)*1.2)
ggplot(sample_quant_summary_labeled,
       aes(x = Combination, y = mean_log10_corr,
           fill = factor(Combination))) +
  geom_col(width = 0.7) +
  geom_point(data = sample_quant_stats,
             aes(x = Combination, y = Log10_copies_corrected,
                 fill = factor(Combination)),
             position = position_jitter(width = 0.2),
             shape = 21,
             size = 1.5,
             stroke = 0.6,
             alpha = 0.7) +
  geom_errorbar(aes(ymin = mean_log10_corr - sd_log10_corr,
                    ymax = mean_log10_corr + sd_log10_corr),
                width = 0.25) +
  geom_text(aes(x = Combination, y = letter_y, label = letters),
            size = 4, fontface = "bold") +
  scale_fill_manual(values = Cq_colors, name = "Combination") +
  scale_y_continuous(
    breaks        = seq(0, y_max, by = 1),      # major grid every 1
    minor_breaks  = seq(0, y_max, by = 0.5)    # minor grid every 0.5
  ) +
  labs(title = "Cells/mL Comparison Across Combinations",
       subtitle = "Letters indicate statistical significance (p<0.05, Wilcoxon test with BH correction)",
       x = "Consortium Combination",
       y = expression(Log[10](cells/mL))) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        legend.position = "none",
        panel.grid.major.x = element_blank(),                 # keep x grid clean
        panel.grid.major.y = element_line(colour = "grey40",  # darker major lines
                                          linewidth = 0.4),
        panel.grid.minor.y = element_line(colour = "grey70",  # extra horizontal lines
                                          linewidth = 0.25),
        panel.grid.minor.x = element_blank())  # Remove legend since it's redundant with x-axis

ggsave("results/growth_phe/qPCR/cells-ml-combinations-comparison-stats.pdf", width = 10, height = 6)

# Save quantification results
write.csv(sample_quantification, "results/growth_phe/qPCR/sample_quantification_results.csv", row.names = FALSE)
