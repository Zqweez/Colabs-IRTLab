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
data <- read_excel("data/gcms_raw/naphtol-calibration.xlsx", sheet = "Samples with phe")

# Convert necessary columns to numeric# Convert necessary columns to numeric# Convert necessary columns to numeric
data$Ret.Time <- as.numeric(data$`Ret Time`)

data <- data[order(data$Sample),]

# data$Sample <- as.numeric(data$Sample)
data <- data %>%
  mutate(
    Date = ymd(sub("-.*", "", data$`File Name`)),
    Sample = as.numeric(Sample),
    SampleDate = interaction(Sample, Date, sep = "_")
  )
# Df for naphtol
summary_nap_df <- data %>%
  filter(Analyte == 'Naphtol') %>%
  group_by(Sample, Date) %>%
  summarise(
    mean = mean(`Nap/decane`, na.rm = TRUE),
    se = sd(`Nap/decane`, na.rm = TRUE) / sqrt(n()),
    n = n(),
    sd = sd(`Nap/decane`, na.rm = TRUE),
    .groups   = "drop"
  )

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
    title = "Naphthol/Decane Ratios Across Samples",
    x = "Duplicate",
    y = "Naphtol/Decane (mean ± SD)",
    fill = "Date"
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")

ggplot(summary_nap_df, aes(x = Sample, y = mean, group = Date)) +
    geom_line(colour = '#e56562') +
    geom_point(data = data,
               aes(x = Sample, y = `Nap/decane`),
               size = 0.7,
               alpha = 1,
               stroke = 1,
               shape = 4) +
    geom_errorbar(
      aes(ymin = mean - sd, ymax = mean + sd),
      position = position_dodge(width = 0.8),
      width = 0.005,
      linewidth = 0.3
    ) +
    coord_cartesian(xlim = c(0, 0.2), ylim = c(0, 1)) +
    labs(
      title = "Naphthol/Decane Ratios Across Samples",
      x = "Naphthol Concentartion (g/L)",
      y = "Naphthol/Decane (mean ± SD)",
      fill = "Date"
    ) +
    theme_minimal() +
    theme(legend.position = "bottom")

ggsave("results/gcms_plots/naphthol.pdf", width = 8, height = 5)
