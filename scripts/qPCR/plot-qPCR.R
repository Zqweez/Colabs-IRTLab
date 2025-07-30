library(ggplot2)
library(ggpattern)
library(janitor)
library(tidyverse)
library(colorspace)  # for desaturate() and lighten()
library(Polychrome)

qpcr_data_raw <- read_excel("data/phe_growth/qPCR/Kasper_2025-07-28 14-45-16_BR005661 -  Quantification Amplification Results.xlsx", sheet = "SYBR")

# Filter out the outliers that are obviously resulting from problems with the machine
qpcr_data_raw <- qpcr_data_raw %>% 
  select(-`EHC-1-1`, -`KS8-10-2-3`, -`KS3-10-2-3`, -`EHC-od2-3`)

# Convert to long form
qpcr_long <- qpcr_data_raw %>%
  slice(-1) %>%
  pivot_longer(
    cols = -Sample,
    names_to = "replicate",
    values_to = "fluorescence"
  ) %>%
  rename(Cycle = Sample) %>%
  mutate(
    bio_rep = sub("-[123]$", "", replicate),
    fluorescence = as.numeric(fluorescence),
    Cycle = as.numeric(Cycle)
  )

rep_levels <- as.character(unique(qpcr_long$bio_rep))
n_clusters     <- length(rep_levels)

set.seed(0)
anchor_pastel <- c("#FF0000", "#00FF00", "#0000FF")
main_colors <- createPalette(n_clusters, anchor_pastel)

pastelise <- function(cols, desat = 0.2, lighten = 0.1) {
  cols %>%
    desaturate(amount = desat) %>%
    lighten(amount = lighten)
}
main_colors <- pastelise(main_colors)
final_palette <- main_colors
names(final_palette) <- rep_levels

# Plot all data no grouping or means
ggplot(qpcr_long, aes(x = Cycle, y = fluorescence, group = replicate, color = bio_rep)) +
  geom_line(alpha = 0.8) +
  geom_point(size  = 0.5,
             alpha = 0.7) +
  scale_color_manual(name = "Sample", values = final_palette) +
  guides(color = guide_legend(
    override.aes = list(size = 1.2),   # icon / key size
    ncol = 2 
  )) +
  labs(title = "qPCR Amplification Curves",
       subtitle = "All replicates",
       y = "Fluorescence",
       x = "Cycle") +
  theme_minimal()
ggsave("results/growth_phe/qPCR/amp-curve-all-replicates.pdf", width = 8, height = 5)

qpcr_summary <- qpcr_long %>%
  group_by(Cycle, bio_rep) %>%
  summarise(
    mean_fluo = mean(fluorescence, na.rm = TRUE),
    sd_fluo = sd(fluorescence, na.rm = TRUE),
    .groups = "drop"
  )

ggplot(qpcr_summary, aes(x = Cycle, y = mean_fluo, color = bio_rep)) +
  geom_line(size = 0.7) +
  geom_errorbar(aes(ymin = mean_fluo - sd_fluo, ymax = mean_fluo + sd_fluo), width = 0.3, alpha = 0.5) +
  geom_point(size = 1.2) +
  scale_color_manual(name = "Sample", values = final_palette) +
  guides(color = guide_legend(
    override.aes = list(size = 1.2),   # icon / key size
    ncol = 2 
  )) +
  labs(title = "Mean Amplification Curves with Standard Deviation", 
       y = "Mean Fluorescence",
       x = "Cycle") +
  theme_minimal()
ggsave("results/growth_phe/qPCR/amp-curve-mean-SD-per-biorep.pdf", width = 12, height = 5)


## ---- Moving on to using the Cq values instead of the full amplification curve ---
# Load the other xlsx file and show the Cq values in bar charts grouped by samples
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
  group_by(bio_rep,Combination) %>%
  summarise(
    mean_cq = mean(Cq, na.rm = TRUE),
    sd_cq = sd(Cq, na.rm = TRUE),
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
ggplot(Cq_summary,
       aes(x = bio_rep, y = mean_cq,
           fill = factor(Combination),    
           group = factor(Combination))) + # makes sure bars & points share groups
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  geom_point(data = Cq_data,
             aes(x = bio_rep, y = Cq,
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
  facet_wrap(~ Combination, scales = "free_x", nrow = 1) +
  scale_fill_manual(values = Cq_colors, name = "Sample") +
  labs(title = "Quantification Cycle for Biological Duplicates",
       subtitle = "Mean of technical replicates with SD as error bars",
       x = "Biological Duplicate",
       y = "Quantification Cycle (Cq)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
        strip.text = element_text(size = 8, angle = 0))

ggsave("results/growth_phe/qPCR/Cq-mean-SD-per-biorep.pdf", width = 8, height = 5)

