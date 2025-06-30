#############################################################################
##### CHARACTERIZING DYNAMIC POPULATIONS IN WASTEWATER MONITORING AREAS #####
############################## SUMMER 2025 ##################################
################################## EJW ######################################

# ANALYSIS 3A. SIMPLE CORRELATIONS BETWEEN DEVICE COUNTS AND FLOW RATES

# set working directory
setwd("/Users/emwu9912/Documents/CU Anschutz/COVID Wastewater/")

# load libraries
library(pacman)
p_load(gridExtra, lubridate, tidyverse)

# read in sewershed names file
sewershed_names <- read.csv("DataRaw/sewershed_names.csv")

# read in mobile device data
load("DataProcessed/visit_with_cluster_sample.Rda")

# read in flow data
load("DataProcessed/cdphe_flow.Rda")

# set up mega dataset for correlation analysis
corr_names <- as.list(unique(cdphe_flow$sewershed_analysis_name))
corr_list <- list()
for (i in seq_along(corr_names)){
  corr_list[[i]] <- cdphe_flow %>% filter(sewershed_analysis_name == corr_names[[i]])
  corr_list_dates <- as.data.frame(seq(min(corr_list[[i]]$date),
                                       max(corr_list[[i]]$date), "days")) %>%
    rename(date = `seq(min(corr_list[[i]]$date), max(corr_list[[i]]$date), "days")`)
  # index to daily timestep with NA's on non-reporting days
  corr_list[[i]] <- merge(corr_list[[i]], corr_list_dates, "date", all = T) %>%
    mutate(sewershed_analysis_name = corr_names[[i]],
           # create extra variable to help with restricting to overlapping date ranges
           daily = "daily")
}

# combine into one big dataset
flow_corr_daily <- bind_rows(corr_list) %>%
  select(date, sewershed_analysis_name, flow_rate, daily)

# merge with mobile device data
flow_mobile_corr_daily <- merge(flow_corr_daily, visit_with_cluster_sample,
                                c("date", "sewershed_analysis_name"), all = T) %>%
  # restrict to overlapping date range
  # if flow data ends earlier than mobility data, end on the last day of the flow data
  # if flow data ends later than mobility data, end on 12/31/2024
  drop_na(daily) %>%
  select(date, sewershed_analysis_name, sewershed_result_name, day_of_week, month, season, year, total_devices_daily, flow_rate, cluster) %>%
  arrange(sewershed_analysis_name)

# save as Rda
save(flow_mobile_corr_daily, file = "DataProcessed/flow_mobile_corr_daily.Rda")

# output to csv
#write.csv(flow_mobile_corr_daily, "DataProcessed/flow_mobile_corr_daily.csv")

# create a single time series for all sewersheds in each cluster with median correlation
# merge with precipitation data
precip <- read.csv("DataRaw/precip.csv") %>%
  rename(cdphe_flow_name = wwtp_name) %>%
  mutate(date = as.Date(date, "%m/%d/%y"),
         month = month(date, label = T, abbr = F)) %>%
  left_join(sewershed_names) %>%
  left_join(visit_with_cluster_sample) %>%
  select(sewershed_result_name, cluster, date, month, melt, rain, total)

# filter to complete observations only 
flow_mobile_corr_complete <- flow_mobile_corr_daily %>%
  drop_na(flow_rate) %>%
  group_by(cluster, month) %>%
  mutate(corr = cor(total_devices_daily, flow_rate, method = c("spearman")))

flow_mobile_corr <- flow_mobile_corr_complete %>%
  add_count() %>%
  slice(1) %>%
  mutate(corr_ci_lower = tanh(atanh(corr) - (1.96/sqrt(n - 3))),
         corr_ci_upper = tanh(atanh(corr) + (1.96/sqrt(n - 3))))

spearman_df1 <- flow_mobile_corr_complete %>%
  select(corr, total_devices_daily) %>%
  arrange(desc(total_devices_daily)) %>%
  slice(1)

spearman_df2 <- flow_mobile_corr_complete %>%
  select(corr, flow_rate) %>%
  arrange(desc(flow_rate)) %>%
  slice(1) %>%
  left_join(spearman_df1) %>%
  as.data.frame()

# scatterplots with correlation coefficients
scatter <- ggplot(data = flow_mobile_corr_complete) +
  geom_point(aes(x = total_devices_daily, y = flow_rate)) +
  labs(x = "Daily Device Count", y = "Flow Rate (mgd)") +
  facet_wrap(cluster ~ month, ncol = 6, scales = "free") +
  scale_x_continuous(labels = scales::comma, expand = expansion(mult = c(0.1, 0.2))) +
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.2))) +
  geom_text(data = spearman_df2, aes(x = total_devices_daily*0.9,
                                     y = flow_rate*1.2,
                                     label = paste("rs =", signif(corr, 3)))); scatter

ggsave("Figures/scatter.png", height = 10, width = 16, plot = scatter)

flow_mobile_corr_summary <- precip %>%
  group_by(cluster, month) %>%
  summarize(mean_melt = mean(melt)/10,
            mean_rain = mean(rain)/10) %>%
  left_join(flow_mobile_corr) %>%
  pivot_longer(cols = c("corr", "mean_melt", "mean_rain"), names_to = "measure",
               values_to = "value") %>%
  mutate(monthnum = month(date, label = F))

corr_precip_high_plot <- ggplot(data = flow_mobile_corr_summary %>% filter(cluster == "High-Variation")) +
  geom_col(aes(x = month, y = value, group = measure, fill = measure),
           position = "dodge", width = 0.5) +
  geom_errorbar(aes(x = monthnum-0.16, ymin = corr_ci_lower, ymax = corr_ci_upper), width = 0.1) +
  labs(x = NULL, y = "Spearman Correlation Coefficient") +
  ggtitle("High-Variation Sewersheds") +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2),
                     sec.axis = sec_axis(transform = ~.*10, name = "Mean Preciciptation (mm)")) +
  scale_fill_manual(values = c("firebrick", "skyblue2", "dodgerblue4"),
                    labels = c("Correlation", "Snowmelt", "Rainfall")) +
  theme(plot.title = element_text(hjust = 0.5),
        plot.margin = margin(t = 20, l = 5),
        axis.text.x = element_text(angle = 90, hjust = 0.95, vjust = 0.4),
        legend.title = element_blank())

corr_precip_low_plot <- ggplot(data = flow_mobile_corr_summary %>% filter(cluster == "Low-Variation")) +
  geom_col(aes(x = month, y = value, group = measure, fill = measure),
           position = "dodge", width = 0.5) +
  geom_errorbar(aes(x = monthnum-0.16, ymin = corr_ci_lower, ymax = corr_ci_upper), width = 0.1) +
  labs(x = NULL, y = "Spearman Correlation Coefficient") +
  ggtitle("Low-Variation Sewersheds") +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2),
                     sec.axis = sec_axis(transform = ~.*10, name = "Mean Preciciptation (mm)")) +
  scale_fill_manual(values = c("firebrick", "skyblue2", "dodgerblue4"),
                    labels = c("Correlation", "Snowmelt", "Rainfall")) +
  theme(plot.title = element_text(hjust = 0.5),
        plot.margin = margin(t = 20, l = 5),
        axis.text.x = element_text(angle = 90, hjust = 0.95, vjust = 0.4),
        legend.title = element_blank())

corr_precip_plot <- grid.arrange(corr_precip_high_plot, corr_precip_low_plot)

ggsave("Figures/corr_precip_plot.png", height = 8, width = 12, plot = corr_precip_plot)

########################################################################################################
