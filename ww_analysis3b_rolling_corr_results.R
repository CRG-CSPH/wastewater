#############################################################################
##### CHARACTERIZING DYNAMIC POPULATIONS IN WASTEWATER MONITORING AREAS #####
############################## SUMMER 2025 ##################################
################################## EJW ######################################

# ANALYSIS 3B. ROLLING CORRELATIONS BETWEEN DEVICE COUNTS AND FLOW RATES (RESULTS)

# set working directory
setwd("/Users/emwu9912/Documents/CU Anschutz/COVID Wastewater/")

# load libraries
library(pacman)
p_load(grid, gridExtra, tidyverse)

# read in rolling correlation results
corr_df <- readRDS("DataProcessed/corr_df.rds")

# read in sewershed names file
sewershed_names <- read.csv("DataRaw/sewershed_names.csv")

# read in mobile device data
load("DataProcessed/visit_with_cluster_sample.Rda")

# read in combined data
load("DataProcessed/flow_mobile_corr_daily_fill.Rda")

corr_sewersheds <- unique(unlist(flow_mobile_corr_daily_fill$sewershed_analysis_name))

corrs <- list()
i_list <- list()
for(i in seq_along(corr_sewersheds)){
  # subsetting into individual sewersheds
  corrs[[i]] <- flow_mobile_corr_daily_fill %>%
    filter(sewershed_analysis_name == corr_sewersheds[[i]]) %>%
    select(-sewershed_analysis_name)
  i_list[[i]] <- corrs[[i]] %>%
    mutate(time = as.numeric(difftime(date, min(corrs[[i]]$date), units = "days"))) %>%
    select(time, flow_rate_7daymovingavg, total_devices_daily_7daymovingavg) %>%
    as.matrix()
  assign(paste0("input_", corr_sewersheds[[i]]), i_list[[i]], envir = .GlobalEnv)
}

corr_df_long <- bind_rows(corr_df) %>%
  left_join(sewershed_names) %>%
  mutate(cluster = case_when(sewershed_analysis_name %in% visit_with_cluster_sample$sewershed_analysis_name[visit_with_cluster_sample$cluster=="High-Variation"] ~ "High-Variation",
                             TRUE ~ "Low-Variation")) %>%
  rename(date = end_date) %>%
  select(date, sewershed_result_name, cluster, cor_pval_run) %>%
  arrange(sewershed_result_name)

corrs_scaled_list <- list()
result_names <- unique(unlist(corr_df_long$sewershed_result_name))
for(i in 1:length(corr_sewersheds)){
  corrs_scaled_list[[i]] <- corrs[[i]] %>% filter(sewershed_result_name == result_names[i])
}

scale_plots <- list()
for(i in 1:length(corr_sewersheds)){
  scale_plots[[i]] <- ggplot(data = corrs[[i]] %>% filter(sewershed_result_name == result_names[i])) +
    geom_line(aes(x = date, y = scale(total_devices_daily_7daymovingavg), color = "color1")) +
    geom_line(aes(x = date, y = scale(flow_rate_7daymovingavg), color = "color2")) +
    labs(x = NULL, y = "Scaled Population Dynamics, 7-Day Moving Average") +
    scale_x_date(limits = c(min(corrs_scaled_list[[i]]$date)-14, max(corrs_scaled_list[[i]]$date)+14),
                 date_labels="%b %Y", date_breaks = "3 months") +
    scale_color_manual(values = c("color1" = "goldenrod2", "color2" = "darkgreen"),
                       labels = c("color1" = "Device Count", "color2" = "Flow Rate")) +
    guides(color = guide_legend(override.aes = list(linewidth = 8))) +
    theme(plot.margin = margin(10, 20, 0, 0, "pt"),
          legend.title = element_blank(),
          legend.position = "bottom",
          axis.title.y = element_text(margin = margin(20, 10, 0, 10, "pt")),
          axis.text.x = element_text(angle = 90, hjust = 0.95, vjust = 0.4))
}

# plot rolling correlations
corr_plots <- list()
for(i in 1:length(corr_sewersheds)){
  corr_plots[[i]] <- ggplot(data = corr_df_long %>% filter(sewershed_result_name == result_names[i])) +
    geom_line(aes(x = date, y = cor_pval_run), color = "black") +
    geom_hline(yintercept = 0, linetype = "dotted") +
    labs(x = "Rolling Window End Date", y = "30-Day Rolling Correlation") +
    scale_x_date(limits = c(min(corrs_scaled_list[[i]]$date)-14, max(corrs_scaled_list[[i]]$date)+14),
                 date_labels="%b %Y", date_breaks = "3 months") +
    scale_y_continuous(limits = c(-1, 1), labels = scales::number_format(accuracy = 0.01)) +
    theme(plot.margin = margin(10, 20, 0, 0, "pt"),
          axis.title.y = element_text(margin = margin(l = 10, r = 10)),
          axis.text.x = element_text(angle = 90, hjust = 0.95, vjust = 0.4))
}

combined_plots <- list()
for(i in seq_along(corr_sewersheds)){
  combined_plots[[i]] <- grid.arrange(scale_plots[[i]], corr_plots[[i]],
                                      top = textGrob(paste0(result_names[i]), hjust = 0.5,
                                                     gp = gpar(fontsize = 20, font = 8)))
  ggsave(paste0("Figures/Rolling Correlation (Cluster)/", corr_sewersheds[[i]], " combined.png"), height = 10, width = 8, plot = combined_plots[[i]])
}

rolling_corr_summary <- corr_df_long %>%
  #mutate(median = median(cor_pval_run)) %>%
  group_by(sewershed_result_name) %>%
  mutate(strong = ifelse(cor_pval_run > 0.7, cor_pval_run, NA),
         n_strong = sum(!is.na(strong)),
         n_total = sum(!is.na(cor_pval_run)),
         prop = n_strong/n_total,
         mean_corr = mean(cor_pval_run),
         sd_corr = sd(cor_pval_run),
         cov = abs(sd_corr/mean_corr)) %>%
  slice(1) %>%
  mutate(strength = ifelse(prop > 0.1, "High Strength", "Low Strength"),
         stability = ifelse(cov < 1, "High Stability", "Low Stability"))

rolling_corr_summary_high <- rolling_corr_summary %>% filter(cluster == "High-Variation")
rolling_corr_summary_low <- rolling_corr_summary %>% filter(cluster == "Low-Variation")

table(rolling_corr_summary_high$strength, rolling_corr_summary_high$stability)
table(rolling_corr_summary_low$strength, rolling_corr_summary_low$stability)

# ranked sewersheds by strength
rolling_corr_summary_id_strength <- rolling_corr_summary %>%
  ungroup() %>%
  arrange(prop) %>%
  #group_by(sewershed_result_name) %>%
  mutate(id = match(sewershed_result_name, unique(sewershed_result_name)))

plot_rolling_corr_strength <- ggplot(data = rolling_corr_summary_id_strength) +
  geom_bar(aes(x = prop, y = as.factor(id), fill = cluster), stat = "identity") +
  geom_vline(xintercept = 0.1, linetype = "dashed") +
  labs(x = "Percent of Rolling Correlation Values > 0.7", y = NULL) +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.1), labels = scales::percent) +
  scale_fill_manual(values = c("orange2", "darkorchid4"),
                    labels = c("High-Variation", "Low-Variation")) +
  guides(fill = guide_legend(title = "Mobility Cluster")) +
  theme(plot.title = element_text(hjust = 0.5),
        axis.title.x = element_text(margin = margin(t = 5)),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())

ggsave("Figures/plot_rolling_corr_strength.png", width = 9, height = 9, plot = plot_rolling_corr_strength)


# ranked sewersheds by stability
rolling_corr_summary_id_stability <- rolling_corr_summary %>%
  ungroup() %>%
  arrange(cov) %>%
  #group_by(sewershed_result_name) %>%
  mutate(id = match(sewershed_result_name, unique(sewershed_result_name)))

plot_rolling_corr_stability <- ggplot(data = rolling_corr_summary_id_stability) +
  geom_bar(aes(x = cov, y = as.factor(id), fill = cluster), stat = "identity") +
  geom_vline(xintercept = 1, linetype = "dashed") +
  labs(x = "Coefficient of Variation (Absolute Value)", y = NULL) +
  scale_x_continuous(breaks = seq(0, 20, 2)) +
  scale_fill_manual(values = c("orange2", "darkorchid4"),
                    labels = c("High-Variation", "Low-Variation")) +
  guides(fill = guide_legend(title = "Mobility Cluster")) +
  theme(plot.title = element_text(hjust = 0.5),
        axis.title.x = element_text(margin = margin(t = 5)),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())

ggsave("Figures/plot_rolling_corr_stability.png", width = 9, height = 9, plot = plot_rolling_corr_stability)

# spaghetti plot showing rolling correlations

rolling_corr_spaghetti_high <- ggplot(data = corr_df_long %>% filter(cluster == "High-Variation")) +
  geom_line(aes(x = date, y = cor_pval_run, group = sewershed_result_name,
                color = sewershed_result_name), alpha = 0.4) +
  labs(x = "90-Day Rolling Window End Date", y = "Dynamic Spearman Correlation Coefficient") +
  ggtitle("High-Variation Sewersheds (n = 15)") +
  scale_x_date(limits = as.Date(c(min(corr_df_long$date), max(corr_df_long$date))),
               date_breaks = "4 months", date_labels = "%b %Y") +
  guides(color = "none") +
  theme(plot.title = element_text(hjust = 0.5),
        axis.title.x = element_text(margin = margin(t = 5)),
        axis.text.x = element_text(angle = 90, hjust = 0.95, vjust = 0.4))

rolling_corr_spaghetti_low <- ggplot(data = corr_df_long %>% filter(cluster == "Low-Variation")) +
  geom_line(aes(x = date, y = cor_pval_run, group = sewershed_result_name,
                color = sewershed_result_name), alpha = 0.4) +
  labs(x = "90-Day Rolling Window End Date", y = "Dynamic Spearman Correlation Coefficient") +
  ggtitle("Low-Variation Sewersheds (n = 51)") +
  scale_x_date(limits = as.Date(c(min(corr_df_long$date), max(corr_df_long$date))),
               date_breaks = "4 months", date_labels = "%b %Y") +
  guides(color = "none") +
  theme(plot.title = element_text(hjust = 0.5),
        axis.title.x = element_text(margin = margin(t = 5)),
        axis.text.x = element_text(angle = 90, hjust = 0.95, vjust = 0.4))

rolling_corr_spaghetti <- grid.arrange(rolling_corr_spaghetti_high, rolling_corr_spaghetti_low)

ggsave("Figures/rolling_corr_spaghetti.png", height = 9, width = 11, plot = rolling_corr_spaghetti)

########################################################################################################

