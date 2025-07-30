library(tidyverse)
library(ggpubr)
library(ggplot2)
library(readxl)
library(rstatix)  # For easier statistical analysis
# Use the Cq values to do analysis and statistical tests.

Cq_data_raw <- read_excel("data/phe_growth/qPCR/Kasper_2025-07-28 14-45-16_BR005661 -  Quantification Cq Results.xlsx", sheet = "0")

Cq_data <- Cq_data_raw %>% 
  select(Well, Fluor, Sample, Cq) %>%
  filter(!(Sample %in% c("EHC-1-1", "KS8-10-2-3", "KS3-10-2-3", "EHC-od-2-3"))) %>%
  mutate(
    bio_rep = sub("-[123]$", "", Sample),
    Combination = sub("-[123]$", "", bio_rep),
    Cq = as.numeric(Cq)
  )

Cq_summary <- Cq_data %>%
  group_by(Combination) %>% # bio_rep,
  summarise(
    mean_cq = mean(Cq, na.rm = TRUE),
    sd_cq = sd(Cq, na.rm = TRUE),
    .groups = "drop"
  )

## Very small sample size so difficult to determine parametric or non-parametric data, Assume non-parametric

# Perform Kruskal-Wallis test
# H₀: All groups come from the same distribution
# H₁: At least one group differs in median.
kw_test <- Cq_data %>% kruskal_test(Cq ~ Combination)
cat("Kruskal-Wallis test results:\n")
print(kw_test)

# Post-hoc pairwise comparisons using rstatix
stat_test <- Cq_data %>%
  wilcox_test(Cq ~ Combination, p.adjust.method = "BH") %>%
  add_significance("p.adj")

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

# Create summary with letters
Cq_summary_labeled <- Cq_summary %>%
  left_join(cld_df, by = "Combination") %>%
  mutate(letter_y = mean_cq + sd_cq + 2)  # adjust position above error bar

# Define variables for plotting
pos_jd <- position_jitterdodge(
  dodge.width  = 0.8,  # same width as for the bars
  jitter.width = 0.10  # sideways spread
)

Cq_colors = c('#e56562', '#51c558', '#67b4ef', "#e6ab02", "#fed996", "#1b9e77", 
              "#7570b3", "#fbd3e1", "#a6cee3", "#b2df8a", '#d95f02', '#f781be')



## Plot the barchart with the significance levels
ggplot(Cq_summary_labeled,
       aes(x = Combination, y = mean_cq,
           fill = factor(Combination),    
           group = factor(Combination))) + # makes sure bars & points share groups
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  geom_point(data = Cq_data,
             aes(x = Combination, y = Cq,
                 fill  = factor(Combination),   # same mapping as the bars
                 group = factor(Combination)),
             position = pos_jd,
             shape = 21,                 # hollow circle whose fill we control
             size  = 0.7,
             stroke = 0.6,
             alpha = 0.7) +
  geom_errorbar(aes(ymin = mean_cq - sd_cq,
                    ymax = mean_cq + sd_cq),
                position = position_dodge(width = 0.8),
                width = 0.25) +
  geom_text(aes(x = Combination, y = letter_y, label = letters),
            size = 3, fontface = "bold") +
  scale_fill_manual(values = Cq_colors, name = "Sample") +
  labs(title = "Quantification Cycle per Consortium",
       subtitle = "Letters indicate statistical significance (p<0.05)",
       x = "Consortium Combination",
       y = "Quantification Cycle (Cq)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
        strip.text = element_text(size = 8, angle = 0))

ggsave("results/growth_phe/qPCR/Cq-mean-SD-per-combination-letters.pdf", width = 8, height = 5)



