library(gplots)
library(tidyverse)
library(lubridate)

growth <- read.table("20250619-isolates.dt3", sep = ",", header = FALSE)
colnames(growth) <- c("point", "time", "na", "KS01-1","KS01-2","KS07-1","KS07-2","KS08-1","KS08-2")
# Remove first datapoint
growth <- growth %>%
  slice(-1) %>%
  mutate(point = point-2) %>%
  mutate(
    time_full = ymd_hms(paste("2025-01-01", time)),
    crossed_midnight = c(FALSE, diff(time_full) < 0),
    day = cumsum(crossed_midnight),
    time_cum = time_full + days(day),
    hour_min = format(time_cum, "%H:%M"),
    time_lbl = paste0("Day ", day, "\n", hour_min)
  )

growth_long <- growth %>%
  pivot_longer(
    cols = c("KS01-1","KS01-2","KS07-1","KS07-2","KS08-1","KS08-2"),
    names_to = "Isolate",
    values_to = "OD600"
  )

breaks_sel <- growth$point[seq(1, nrow(growth), by = 20)]

ggplot(growth_long, aes(x = point, y = OD600, colour = Isolate)) +
  geom_line(size = 0.4) +
  geom_point(size = 0.2) +
  scale_x_continuous(
    name = "Datapoint",
    breaks = breaks_sel,
    sec.axis = dup_axis(
      breaks = breaks_sel,
      labels = growth$time_lbl[match(breaks_sel, growth$point)],
      name   = "Day & Time"
    )
  ) +
  labs(y = expression(OD[600]), title = "Isolate Growth Curves") +
  theme_minimal() +
  theme(
    text = element_text(size = 12),
    axis.text.x      = element_text(size = 10),  # bottom axis numbers
    axis.text.x.top  = element_text(size = 8, angle = 45, hjust = 0),
    plot.title = element_text(hjust = 0.5),
    axis.title.x.top = element_text(size = 10, margin = margin(b = 8))
  )

ggsave("Isolate-growth-curves.pdf", width = 8, height = 5)

# Do statistical analysis of the data

library(growthcurver)
growth_long_stat <- growth_long %>%
  select(-na, -crossed_midnight, -hour_min) %>%
  mutate(
    Replicate = sub(".*-", "", Isolate),
    Isolate =  sub("-.*", "", Isolate),
    OD600 = as.numeric(OD600)
  )

fits <- growth_long_stat %>%                      
  group_by(Isolate, Replicate) %>%                  
  group_modify(~ {
    sg <- SummarizeGrowth(.x$point, .x$OD600)        
    as_tibble(sg$vals)                              
  }) %>%
  ungroup()

param_tbl <- fits %>%
  select(Bacteria, Replicate, mu, lag, K)

anova_mu  <- aov(mu  ~ Bacteria, data = param_tbl)
anova_lag <- aov(lag ~ Bacteria, data = param_tbl)
anova_K   <- aov(K   ~ Bacteria, data = param_tbl)

summary(anova_mu)     # F & p-value for μmax
TukeyHSD(anova_mu)    # post-hoc pairwise


# Overlay fitted curves
df_preds <- fits %>%
  mutate(pred_grid = map2(.x = Isolate, .y = Replicate, ~ {
    tibble(time = seq(min(growth_long_stat$time), max(growth_long_stat$time), length.out = 200))
  })) %>%
  unnest(pred_grid) %>%
  mutate(OD_fit = K / (1 + exp(-r*(time - t_mid))))   # change formula for Gompertz/Baranyi

ggplot(growth_long_stat, aes(time, OD600, colour = Isolate)) +
  geom_point(alpha = 0.6) +
  geom_line(data = df_preds, aes(y = OD_fit), size = 1) +
  facet_wrap(~Isolate) + theme_minimal()

