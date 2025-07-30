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

# To correct for volume change, initial volume 5 mL, after first sample ~4 mL, before final ~1 mL
# Thus divide all concentrations for the last date by 2 to be a bit conservative
data <- data %>%
  mutate(
    `Phe/decane` = if_else(
      Date == as.Date("2025-07-22") & !Sample %in% c("NC-1", "NC-2"),   # rows to alter
      if_else(Sample %in% c("EHC-OD-1", "KS3-100-1"), `Phe/decane` / 1.25, `Phe/decane` / 2),
      `Phe/decane`),
    `Nap/decane`= if_else(
      Date == as.Date("2025-07-22"),   # rows to alter
      `Nap/decane` / 2,
      `Nap/decane`)
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

# Color by date for Phenanthrene
# Helper for jitter
pos_jd <- position_jitterdodge(
  dodge.width  = 0.8,  # same width as for the bars
  jitter.width = 0.10  # sideways spread
)
## For phenanthrene
ggplot(summary_df,
       aes(x = Duplicate,
           y = mean,
           fill = factor(Date),     # colour code by date
           group = factor(Date))) + # makes sure bars & points share groups
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  geom_point(data = data,
             aes(x = Duplicate,
                 y = `Phe/decane`,
                 fill  = factor(Date),   # same mapping as the bars
                 group = factor(Date)),
             position = pos_jd,
             shape = 21,                 # hollow circle whose fill we control
             size  = 0.7,
             stroke = 0.6,
             alpha = 0.7) +
  geom_errorbar(aes(ymin = mean - sd,
                    ymax = mean + sd),
                position = position_dodge(width = 0.8),
                width = 0.25) +
  
  facet_wrap(~ Combination, scales = "free_x", nrow = 1) +
  scale_fill_manual(values = date_colors, name = "Date") +
  labs(title = "Phenanthrene/Decane Ratios Across Samples",
       subtitle = "Corrected",
       x = "Duplicate",
       y = "Phenanthrene/Decane (mean ± SD)") +
  theme_minimal() +
  theme(legend.position = "bottom")
ggsave("results/gcms_plots/l-tubes-phe-corrected-final.pdf", width = 8, height = 5)

# Color by date for Naphthol
ggplot(summary_nap_df,
       aes(x = Duplicate,
           y = mean,
           fill = factor(Date),     # colour code by date
           group = factor(Date))) + # makes sure bars & points share groups
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  geom_point(data = data,
             aes(x = Duplicate,
                 y = `Nap/decane`,
                 fill  = factor(Date),   # same mapping as the bars
                 group = factor(Date)),
             position = pos_jd,
             shape = 21,                 # hollow circle whose fill we control
             size  = 0.7,
             stroke = 0.6,
             alpha = 0.7) +
  geom_errorbar(aes(ymin = mean - sd,
                    ymax = mean + sd),
                position = position_dodge(width = 0.8),
                width = 0.25) +
  
  facet_wrap(~ Combination, scales = "free_x", nrow = 1) +
  scale_fill_manual(values = date_colors, name = "Date") +
  labs(title = "Phenanthrene/Decane Ratios Across Samples",
       x = "Duplicate",
       y = "Phenanthrene/Decane (mean ± SD)") +
  theme_minimal() +
  theme(legend.position = "bottom")
ggsave("results/gcms_plots/l-tubes-nap-corrected-final.pdf", width = 8, height = 5)
