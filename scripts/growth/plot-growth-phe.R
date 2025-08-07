library(gplots)
library(tidyverse)
library(lubridate)
library(slider)
# ---- Samples are
# This part is to modify the raw data in dt files and save it as CSV files with a correct header
samples <- c("20250702-m1.dt1", "20250702-m2.dt2", "20250702-m4.dt4")
# samples <- c("20250714-m1.dt1", "20250714-m2.dt2", "20250714-m4.dt4")
sample <- samples[3]

growth <- read.table(paste0("./data/phe_growth/",sample), sep = ",", header = FALSE)
sample <- sub("\\.dt[0-9]+$", "", sample)  # Remove the .dt1, .dt2, etc. from the sample name
# Load the header from the txt file
header <- read.csv(paste0("./data/phe_growth/",sub(".*-","",sample), ".txt"), sep=",", header = F)
colnames(growth) <- header[1,]

# The Biophotorecorder is a bit unreliable when it comes to when it starts to measure either it starts at 0 or below 0
# We correct for this by either removing or keeping the first datapoint 
growth <- growth %>% 
  slice(-1) %>%
  mutate(point = point - 2) %>%             
  mutate(
    time_full        = ymd_hms(paste("2025-01-01", time)),
    crossed_midnight = c(FALSE, diff(time_full) < 0),     # TRUE each rollover
    time_cum         = time_full + days(cumsum(crossed_midnight))
  ) %>%
  mutate(
    elapsed_sec = as.numeric(difftime(time_cum, first(time_cum), units = "secs")),
    day         = elapsed_sec %/% 86400,                  # whole days
    hour        = (elapsed_sec %% 86400) %/% 3600,        # hours within day
    minute      = (elapsed_sec %% 3600)  %/% 60,          # minutes within hour
    elapsed_h = as.numeric(difftime(time_cum, first(time_cum), units = "hours")),
    time_lbl    = sprintf("Day %d\n%02d:%02d", day, hour, minute)
  ) %>%
  select(-matches("^na"), -crossed_midnight, -time_full, -time_cum, -elapsed_sec, -day, -hour, -minute)
write.csv(growth, paste0("./data/phe_growth/",sample, ".csv"), row.names = TRUE)

# After converting the dt file load the CSV of interest

growth_long <- growth %>%                          
  pivot_longer(
    cols      = matches("^(KS|EHC|NC)"),  # KS01-1, KS07-2, …
    names_to  = "Iso_rep",        # new column = old col-name
    values_to = "OD600"
  ) %>% 
  separate(Iso_rep, into = c("Isolate", "Replicate"),   # KS01-1 → KS01 | 1
           sep = "_", remove = FALSE) %>% 
  mutate(
    Replicate = factor(Replicate),
    Isolate   = factor(Isolate)
  )

# Calculate the sliding median
sliding_window <- 5
growth_smooth <- growth_long %>% 
  group_by(Isolate, Replicate) %>% 
  mutate(OD_smooth = slide_dbl(OD600,
                               median,
                               .before = floor(sliding_window/2), .after = floor(sliding_window/2),
                               .complete = TRUE)) %>% 
  ungroup()

breaks_sel <- growth$point[seq(1, nrow(growth), by = round(nrow(growth)/16))] # total 16 points

growth_plot <- growth_smooth %>%
  mutate(Iso_rep = sub("_", " ", Iso_rep))

combination_colors <- c(
  '#e56562',  # warm red
  '#51c558',  # green
  '#67b4ef',  # light blue
  '#e6ab02',  # mustard yellow
  '#fed996',  # light peach
  '#1b9e77',  # dark green
  '#7570b3',  # purple
  '#fbd3e1',  # pink
  '#a6cee3',  # sky blue
  '#b2df8a',  # pastel green
  '#ff7f00',  # orange
  '#cab2d6',  # lavender
  '#ffff99',  # yellow
  '#fb9a99',  # salmon
  '#66c2a5',  # teal
  '#fc8d62',  # coral
  '#8da0cb',  # steel blue
  '#e78ac3'   # rose
)

ggplot(growth_plot,
       aes(x = point, y = OD_smooth,
           color = Iso_rep          # colour → strain family
           )) +     # group  = interaction(Isolate, Replicate) different symbol for -1 / -2
  geom_line(linewidth = 0.4) +
  geom_point(aes(y = OD600),size = 0.02, alpha = 0.6) +
  scale_x_continuous(
    name = "Datapoint",
    breaks = breaks_sel,
    sec.axis = dup_axis(
      breaks = breaks_sel,
      labels = growth$time_lbl[match(breaks_sel, growth$point)],
      name   = "Day & Time since start"
    )
  ) +
  scale_fill_manual(values = combination_colors) +
  # scale_colour_brewer(palette = "Set1") +
  labs(y = expression(OD[600]),
       colour  = "Isolate",
       title   = "Growth curves for Communities") +
  theme_minimal() +
  coord_cartesian(ylim = c(-0.6, 1)) +
  theme(
    text = element_text(size = 12),
    axis.text.x      = element_text(size = 10),  # bottom axis numbers
    axis.text.x.top  = element_text(size = 8, angle = 45, hjust = 0),
    plot.title = element_text(hjust = 0.5),
    axis.title.x.top = element_text(size = 10, margin = margin(b = 8))
  )

ggsave(paste0("./results/growth_phe/", sample,"-growth-curves.pdf"), width = 8, height = 5)
