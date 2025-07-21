library(gplots)
library(tidyverse)
library(lubridate)
library(slider)
# ---- Samples are
# This part is to modify the raw data in dt files and save it as CSV files with a correct header
samples <- c("20250702-m1.dt1", "20250702-m2.dt2", "20250702-m4.dt4")
samples <- c("20250702-m1.dt1", "20250702-m2.dt2", "20250702-m4.dt4")
sample <- samples[1]

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

ggplot(growth_plot,
       aes(x = point, y = OD_smooth,
           colour = Iso_rep,          # colour → strain family
           group  = interaction(Isolate, Replicate))) +     # different symbol for -1 / -2
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
  scale_colour_brewer(palette = "Set1") +
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

# Do statistical analysis of the data
# 
# library(dplyr)
# library(tidyr)
# library(growthcurver)
# growth_long_stat <- growth_long
# 
# fits <- growth_long_stat %>% 
#   group_by(Isolate, Replicate) %>% 
#   summarize(
#     model = list(
#       SummarizeGrowth(elapsed_h, OD600)   # returns a gcfit object
#     ),
#     .groups = "drop"                      # keep one row per group
#   ) %>%
#   mutate(
#     vals = map(model, ~ .x$vals),          # extract the parameters
#     data = map(model, ~ .x$data)           # extract the data
#   )
# # Extract the data from the fits and make it into long format, ready for plotting
# pred <- fits %>%                                   
#   mutate(pred_df = map(data, \(d)                   
#                        tibble(elapsed_h = d$t,                       
#                               OD_fit    = d$N)                       
#   )) %>%
#   select(Isolate, Replicate, pred_df) %>%          
#   unnest(pred_df) %>%                             
#   mutate(Replicate = factor(Replicate)) 
# 
# # Plot the fitted curves
# shape_vals    <- c(`1` = 16, `2` = 17)             
# linetype_vals <- c(`1` = "dashed", `2` = "dotted")
# 
# ggplot(data = growth_long, 
#        aes(x = elapsed_h, y = OD600,
#            colour = Isolate,
#            shape  = factor(Replicate))) +
#   geom_point(size = 0.5, alpha = 0.8) +
#   geom_line(linewidth = 0.1) +
#   geom_line(data = pred,
#             aes(x = elapsed_h, y = OD_fit,
#                 colour   = Isolate,
#                 linetype = Replicate,                       # dashed vs dotted
#                 group    = interaction(Isolate, Replicate)),
#             linewidth = 0.5) +
#   scale_shape_manual(values = shape_vals,    name = "Replicate") +
#   scale_linetype_manual(values = linetype_vals, name = "Replicate") +
#   scale_colour_brewer(palette = "Set1") +
#   scale_x_continuous(
#     name   = "Datapoint",
#     breaks = breaks_sel,
#     sec.axis = dup_axis(
#       breaks = breaks_sel,
#       labels = growth$time_lbl[match(breaks_sel, growth$point)],
#       name   = "Day & Time since start"
#     )
#   ) +
#   labs(y = expression(OD[600]),
#        colour  = "Isolate",
#        title   = "Growth curves: raw vs model") +
#   theme_minimal(base_size = 12) +
#   theme(
#     axis.text.x.top  = element_text(size = 8, angle = 45, hjust = 0),
#     plot.title       = element_text(hjust = 0.5)
#   )
