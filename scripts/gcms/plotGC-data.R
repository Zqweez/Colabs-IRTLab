library(dplyr)
library(ggplot2)
library(ggpattern)
library(gtools) 
library(readxl)
library(lubridate)     # nicer date handling

## Load the CSV formatted in excel
# file_path <- "250528-GCMS-v1.csv"

# Read the CSV file with comma as decimal separator
# data <- read.csv(file_path, dec = ",", comment.char="", quote="", sep=";", header=TRUE)
# Complete-GCMS_data.xlsx
data <- read_excel("data/gcms_raw/L-tubes.xlsx", sheet = "Samples with phe")

# Convert necessary columns to numeric# Convert necessary columns to numeric# Convert necessary columns to numeric
data$Ret.Time <- as.numeric(data$`Ret Time`)

data <- data[order(data$Sample),]

# data$Sample <- as.numeric(data$Sample)
data <- data %>%
  mutate(
    Date = ymd(sub("-.*", "", data$`File Name`)),
    Sample = toupper(as.factor(Sample)),
    Combination = sub("-[12]$", "", Sample),
    Duplicate = sub(".*-([12])$", "\\1", Sample),
    SampleDate = interaction(Sample, Date, sep = "_")
  )

# DF with analytic metrics
summary_df <- data %>%
  filter(Analyte == 'Phe') %>%
  group_by(Sample, Combination, Duplicate, Date) %>%
  summarise(
    mean = mean(`Phe/decane`, na.rm = TRUE),
    se = sd(`Phe/decane`, na.rm = TRUE) / sqrt(n()),
    n = n(),
    sd = sd(`Phe/decane`, na.rm = TRUE),
    .groups   = "drop"
  )
summary_df <- summary_df %>%
  mutate(
    Sample = factor(Sample, levels = mixedsort(unique(Sample))), # 1,2,3,…
    Date   = factor(Date,   levels = sort(unique(Date)))         # oldest → newest
  )
# Df for naphtol
summary_nap_df <- data %>%
  filter(Analyte == 'Naphtol') %>%
  group_by(Sample, Combination, Duplicate, Date) %>%
  summarise(
    mean = mean(`Nap/decane`, na.rm = TRUE),
    se = sd(`Nap/decane`, na.rm = TRUE) / sqrt(n()),
    n = n(),
    sd = sd(`Nap/decane`, na.rm = TRUE),
    .groups   = "drop"
  )
summary_nap_df <- summary_nap_df %>%
  mutate(
    Sample = factor(Sample, levels = mixedsort(unique(Sample))), # 1,2,3,…
    Date   = factor(Date,   levels = sort(unique(Date)))         # oldest → newest
  )

# Columnplot for each sample
date_colors = c('#e56562', '#51c558', '#67b4ef',   "#e6ab02", "#fed996", "#1b9e77", "#7570b3", "#f781be",
                "#a6cee3", "#b2df8a", '#d95f02', '#fbd3e1'
)

# Columnplot for each combination color by combination
ggplot(summary_df, aes(x = Duplicate, y = mean, fill = factor(Combination))) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  geom_errorbar(
    aes(ymin = mean - sd, ymax = mean + sd),
    position = position_dodge(width = 0.8),
    width = 0.25
  ) +
  geom_jitter(data = data,
              aes(x = Duplicate, y = `Phe/decane`),
              width = 0.1,
              size = 0.5,
              alpha = 0.7,
              inherit.aes = FALSE) +
  facet_wrap(~ Combination, scales = "free_x", nrow = 1) +  # group by Isolate
  scale_fill_manual(values = date_colors, name = "Taxa") +
  labs(
    title = "Phenanthrene Ratios Across Samples",
    subtitle = "Grouped by Sample",
    x = "Sample",
    y = "Phenanthrene/NC2 (mean ± SD)",
    fill = "Date"
  ) +
  theme_minimal() +
  theme(legend.position = "bottom", axis.text.x = element_text(angle = 0, size = 7))
ggsave("results/gcms_plots/l-tubes.pdf", width = 8, height = 5)

# Color by date for Phenanthrene
ggplot(summary_df, aes(x = Duplicate, y = mean, fill = factor(Date))) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  geom_jitter(data = data,
              aes(x = Duplicate, y = `Phe/decane`),
              width = 0.1,
              size = 0.5,
              alpha = 0.7,
              inherit.aes = FALSE) +
  geom_errorbar(
    aes(ymin = mean - sd, ymax = mean + sd),
    position = position_dodge(width = 0.8),
    width = 0.25
  ) +
  facet_wrap(~ Combination, scales = "free_x", nrow = 1) +  # group by Isolate
  scale_fill_manual(values = date_colors, name = "Date") +
  labs(
    title = "Phenanthrene/Decane Ratios Across Samples",
    x = "Duplicate",
    y = "Phenanthrene/Decane (mean ± SD)",
    fill = "Date"
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")
ggsave("results/gcms_plots/l-tubes-phe-dates.pdf", width = 8, height = 5)

# Color by date for Naphtol
ggplot(summary_nap_df, aes(x = Duplicate, y = mean, fill = factor(Date))) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  geom_jitter(data = data,
              aes(x = Duplicate, y = `Nap/decane`),
              width = 0.1,
              size = 0.5,
              alpha = 0.7,
              inherit.aes = FALSE) +
  geom_errorbar(
    aes(ymin = mean - sd, ymax = mean + sd),
    position = position_dodge(width = 0.8),
    width = 0.25
  ) +
  facet_wrap(~ Combination, scales = "free_x", nrow = 1) +  # group by Isolate
  scale_fill_manual(values = date_colors, name = "Date") +
  labs(
    title = "Naphtol/Decane Ratios Across Samples",
    x = "Duplicate",
    y = "Naphtol/Decane (mean ± SD)",
    fill = "Date"
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")
ggsave("results/gcms_plots/l-tubes-nap-dates.pdf", width = 8, height = 5)
