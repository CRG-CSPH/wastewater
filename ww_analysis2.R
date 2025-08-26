#############################################################################
##### CHARACTERIZING DYNAMIC POPULATIONS IN WASTEWATER MONITORING AREAS #####
############################## SUMMER 2025 ##################################
################################## EJW ######################################

# ANALYSIS 2. DESCRIPTIVE STATISTICS AND REGRESSION ON DEVICE COUNTS AND FLOW RATES

# set working directory
setwd("/Users/emwu9912/Documents/CU Anschutz/COVID Wastewater/")

# load libraries
library(pacman)
p_load(distributional, EnvStats, ggdist, gridExtra, Hmisc, lme4, lubridate, tidyverse, viridis)

#turn off scientific notation
options(scipen = 999)

# read in sewershed names file
sewershed_names <- read.csv("DataRaw/sewershed_names.csv")

# read in SS/CBG proportion data from ArcGIS
cbg_proportions <- read.csv("DataRaw/new_overlap_jan2025.csv") %>%
  select(-AREA) %>%
  rename(cdphe_flow_name = wwtp) %>%
  mutate(GEOID = as.character(paste0("0", GEOID))) %>%
  left_join(sewershed_names) %>%
  select(GEOID, sewershed_analysis_name, PERCENTAGE) %>%
  group_by(GEOID) %>%
  # determine in how many sewersheds each CBG shows up
  add_count()

# assign CBGs
cbg_sewershed_assign <- cbg_proportions %>%
  group_by(GEOID, sewershed_analysis_name) %>%
  # winner take all: if a CBG's land mass is occupied by more than one sewershed,
  # assign the CBG with the greatest overlap
  arrange(desc(n), GEOID, desc(PERCENTAGE)) %>%
  group_by(GEOID) %>%
  slice(1) %>%
  # verify that each CBG now shows up only once
  select(-n) %>%
  add_count() %>%
  # 0.2% areal overlap: if a sewershed occupies >= 0.2% of a CBG's land mass,
  # assign that CBG to the sewershed
  group_by(sewershed_analysis_name) %>%
  filter(PERCENTAGE >= 0.2)

# calculate the number of CBG's in each sewershed
cbg_count <- cbg_sewershed_assign %>%
  select(-n) %>%
  group_by(sewershed_analysis_name) %>%
  add_count() %>%
  slice(1)

# create CBG table for appendix
write.csv(cbg_count, "DataProcessed/cbg_count.csv")

# read in monthly device count data
visit_monthly <- read.csv("DataRaw/2025_04_04_monthly_devices.csv") %>%
  rename(cdphe_flow_name = AREA_SEWERSHED,
         internal_devices_monthly = NUMBER_DEVICES_RESIDING,
         external_devices_monthly = SOURCE_AREA_DEVICE_COUNT) %>%
  mutate(total_devices_monthly = internal_devices_monthly + external_devices_monthly) %>%
  mutate(date = as.Date(DATE_RANGE_START, "%Y-%m-%d")) %>%
  select(cdphe_flow_name, date, internal_devices_monthly, external_devices_monthly, total_devices_monthly)

# read in daily device count data
visit_daily <- read.csv("DataRaw/2025_04_04_daily_visits_sum.csv") %>%
  rename(cdphe_flow_name = AREA_SEWERSHED,
         total_devices_daily = STOPS_BY_DAY_L) %>%
  mutate(date = as.Date(DAY, "%Y-%m-%d")) %>%
  select(cdphe_flow_name, date, total_devices_daily)

# merge monthly and daily
visit <- merge(visit_monthly, visit_daily, c("date", "cdphe_flow_name"), all = T) %>%
  filter(date <= "2024-12-31") %>%
  mutate(year = year(date)) %>%
  left_join(sewershed_names) %>%
  select(-sewershed_shapefile_name)

# how many sewersheds in the mobile device data?
unique(unlist(visit$sewershed_analysis_name))

# read in CDPHE flow data (2020-2023)
cdphe_flow1 <- read.csv("DataRaw/ww_flow_pop_data1.csv") %>%
  mutate(date = as.Date(sample_collect_date, "%m/%d/%Y")) %>%
  filter(date <= "2023-12-31") %>%
  select(wwtp_name, date, flow_rate)

# read in CDPHE flow data (2024)
cdphe_flow2 <- read.csv("DataRaw/ww_flow_pop_data2.csv") %>%
  filter(pcr_target == "sars-cov-2") %>%
  mutate(date = as.Date(sample_collect_date, "%m/%d/%y")) %>%
  filter(date <= "2024-12-31") %>%
  select(wwtp_name, date, flow_rate)

# combine datasets
cdphe_flow3 <- bind_rows(cdphe_flow1, cdphe_flow2) %>%
  rename(cdphe_flow_name = wwtp_name) %>%
  group_by(cdphe_flow_name, date) %>%
  slice(1) %>%
  left_join(sewershed_names) %>%
  #left_join(visit_with_cluster) %>%
  mutate(reporting_month = month(date),
         reporting_year = year(date))

# how many sewersheds in the initial flow data?
unique(unlist(cdphe_flow3$cdphe_flow_name))

# drop missing flow data
cdphe_flow4 <- cdphe_flow3 %>%
  drop_na(flow_rate)

# determine how many sewersheds began reporting in August 2020
reporting_start <- cdphe_flow4 %>%
  group_by(sewershed_analysis_name) %>%
  arrange(date) %>%
  slice(1) %>%
  rename(reporting_month_start = reporting_month,
         reporting_year_start = reporting_year) %>%
  filter(reporting_month_start == 8,
         reporting_year_start == 2020)

# determine how many sewersheds reported until December 2024
reporting_end <- cdphe_flow4 %>%
  group_by(sewershed_analysis_name) %>%
  arrange(desc(date)) %>%
  slice(1) %>%
  rename(reporting_month_end = reporting_month,
         reporting_year_end = reporting_year) %>%
  filter(reporting_month_end == 12,
         reporting_year_end == 2024)

# determine how many sewersheds reported for the whole period
reporting_whole <- reporting_end %>%
  filter(cdphe_flow_name %in% reporting_start$cdphe_flow_name) 

cdphe_flow <- cdphe_flow4 %>%
  group_by(sewershed_analysis_name) %>%
  mutate(ndays = as.numeric(difftime(max(date), min(date), "days"))) %>%
  # exclude sewersheds that reported over less than a 365-day period
  filter(ndays >= 365) %>%
  ungroup() %>%
  select(cdphe_flow_name, date, flow_rate, contains("reporting"), ndays)

save(cdphe_flow, file = "DataProcessed/cdphe_flow.Rda")

# merge mobile device and flow data
visit_with_flow <- merge(visit, cdphe_flow, c("date", "cdphe_flow_name"), all= T) %>%
  # restrict to sewersheds appearing in flow data
  filter(cdphe_flow_name %in% cdphe_flow$cdphe_flow_name) %>%
  # impute to check log transformation
  mutate(total_devices_daily = ifelse(total_devices_daily == 0, 1, total_devices_daily),
         flow_rate = ifelse(flow_rate == 0, 0.001, flow_rate))
         
# check distribution of device counts by individual sewershed (using daily data)
sewersheds_analysis <- unique(unlist(visit_with_flow$sewershed_analysis_name))
sewersheds_result <- unique(unlist(visit_with_flow$sewershed_result_name))
sewershed_device_dist <- list()
for (i in seq_along(sewersheds_result)){
  sewershed_device_dist[[i]] <- ggplot(data = visit_with_flow %>% filter(sewershed_analysis_name == sewersheds_analysis[i])) +
    geom_histogram(aes(x = total_devices_daily)) +
    labs(x = "Daily Device Count", y = "Frequency (Sewershed-Days)") +
    ggtitle(sewersheds_result[i]) +
    scale_x_continuous(labels = scales::comma) +
    theme(plot.title = element_text(hjust = 0.5))
  ggsave(paste0("Figures/Device Distribution/", sewersheds_analysis[[i]], "_original.png"), height = 8, width = 8, sewershed_device_dist[[i]])
}

# check log-transformed distributions
sewershed_device_dist_log <- list()
for (i in seq_along(sewersheds_result)){
  sewershed_device_dist_log[[i]] <- ggplot(data = visit_with_flow %>% filter(sewershed_analysis_name == sewersheds_analysis[i])) +
    geom_histogram(aes(x = log(total_devices_daily))) +
    labs(x = "log(Daily Device Count)", y = "Frequency (Sewershed-Days)") +
    ggtitle(paste(sewersheds_result[i], "(Log-Transformed)")) +
    scale_x_continuous(labels = scales::comma) +
    theme(plot.title = element_text(hjust = 0.5))
  ggsave(paste0("Figures/Device Distribution/", sewersheds_analysis[[i]], "_transformed.png"), height = 8, width = 8, sewershed_device_dist_log[[i]])
}

# compare skewness of raw vs. log-transformed device counts
skewness <- visit_with_flow %>%
  group_by(sewershed_result_name) %>%
  mutate(skewness = skewness(total_devices_daily),
         skewness_log = skewness(log(total_devices_daily)),
         normal_skewed = ifelse(abs(skewness) < 0.5, "Normal", "Skewed"),
         normal_skewed_log = ifelse(abs(skewness_log) < 0.5, "Normal", "Skewed"))%>%
  slice(1) %>%
  select(sewershed_result_name, skewness, skewness_log, normal_skewed, normal_skewed_log)

# 2x2 tables describing skew
table(skewness$normal_skewed)
table(skewness$normal_skewed_log)
# better when log-transformed so use log-transformed mobile device data

# check distribution of flow rates by individual sewershed
sewershed_flow_dist <- list()
for (i in seq_along(sewersheds_result)){
  sewershed_flow_dist[[i]] <- ggplot(data = visit_with_flow %>% filter(sewershed_analysis_name == sewersheds_analysis[i])) +
    geom_histogram(aes(x = flow_rate)) +
    labs(x = "Flow Rate", y = "Frequency (Sewershed-Days)") +
    ggtitle(sewersheds_result[i]) +
    #scale_x_continuous(labels = scales::comma) +
    theme(plot.title = element_text(hjust = 0.5))
  ggsave(paste0("Figures/Flow Rate Distribution/", sewersheds_analysis[[i]], "_original.png"), height = 8, width = 8, sewershed_flow_dist[[i]])
}

# check log-transformed flow rate distributions
sewershed_flow_dist_log <- list()
for (i in seq_along(sewersheds_result)){
  sewershed_flow_dist_log[[i]] <- ggplot(data = visit_with_flow %>% filter(sewershed_analysis_name == sewersheds_analysis[i])) +
    geom_histogram(aes(x = log(flow_rate))) +
    labs(x = "log(Flow Rate)", y = "Frequency (Sewershed-Days)") +
    ggtitle(sewersheds_result[i]) +
    theme(plot.title = element_text(hjust = 0.5))
  ggsave(paste0("Figures/Flow Rate Distribution/", sewersheds_analysis[[i]], "_transformed.png"), height = 8, width = 8, sewershed_flow_dist_log[[i]])
}

# compare skewness of raw vs. log-transformed flow rates
skewness_flow <- visit_with_flow %>%
  group_by(sewershed_analysis_name) %>%
  drop_na(flow_rate) %>%
  mutate(skewness_flow = skewness(flow_rate),
         skewness_flow_log = skewness(log(flow_rate)),
         normal_skewed_flow = ifelse(abs(skewness_flow) < 0.5, "Normal", "Skewed"),
         normal_skewed_flow_log = ifelse(abs(skewness_flow_log) < 0.5, "Normal", "Skewed")) %>%
  slice(1) %>%
  select(sewershed_analysis_name, sewershed_result_name, skewness_flow, skewness_flow_log, normal_skewed_flow, normal_skewed_flow_log)

# 2x2 tables describing skew for flow data
table(skewness_flow$normal_skewed_flow)
# log-transforming does not change how many are normal vs. skewed
# (but this is because some turn from skewed to normal and vice versa)
# however, log-transform flow rates anyway to stay consistent with mobile device data
table(skewness_flow$normal_skewed_flow_log)

# log-transform device counts and flow rates
visit_with_flow_log <- visit_with_flow %>%
  # impute daily device counts of zero as 1
  mutate(total_devices_daily_log = log(total_devices_daily),
         total_devices_monthly_log = log(total_devices_monthly),
         internal_devices_monthly_log = log(internal_devices_monthly),
         flow_rate_log = log(flow_rate))

# read in population data to compare device counts by year
populations <- read.csv("DataRaw/co-est2024-pop-08.csv") %>%
  mutate(county = str_remove(county, ", Colorado"),
         est_baseline = as.numeric(str_remove(est_baseline, ",")),
         est_2020 = as.numeric(str_remove(est_2020, ",")),
         est_2021 = as.numeric(str_remove(est_2021, ",")),
         est_2022 = as.numeric(str_remove(est_2022, ",")),
         est_2023 = as.numeric(str_remove(est_2023, ",")),
         est_2024 = as.numeric(str_remove(est_2024, ",")),
         change_5yr = as.numeric(str_remove(change_5yr, "%")),
         change_2020_2022 = as.numeric(str_remove(change_2020_2022, "%")),
         change_2022_2024 = as.numeric(str_remove(change_2022_2024, "%")))

# read in county data
counties <- read.csv("DataRaw/counties.csv") %>%
  group_by(cdphe_flow_name) %>%
  slice(1) %>%
  select(-date) %>%
  left_join(sewershed_names) %>%
  left_join(populations) %>%
  select(cdphe_flow_name, sewershed_result_name, contains("est_2")) %>%
  pivot_longer(!c(cdphe_flow_name, sewershed_result_name), names_to = "year", values_to = "county_pop") %>%
  mutate(year = as.numeric(str_remove(year, "est_"))) %>%
  arrange(sewershed_result_name)

yearly_average <- visit_with_flow_log %>%
  group_by(sewershed_result_name) %>%
  mutate(total_devices_daily_z = scale(total_devices_daily_log)) %>%
  group_by(sewershed_result_name, year) %>%
  mutate(median_year_z = median(total_devices_daily_z)) %>%
  slice(1) %>%
  left_join(counties) %>%
  group_by(sewershed_result_name) %>%
  mutate(county_pop_z = scale(county_pop))

devices_vs_pop <- list()
for (i in seq_along(sewersheds_result)){
  devices_vs_pop[[i]] <- ggplot(data = yearly_average %>% filter(sewershed_result_name == sewersheds_result[i])) +
    geom_point(aes(x = year, y = total_devices_daily_z), color = "dodgerblue3") +
    geom_line(aes(x = year, y = total_devices_daily_z), color = "dodgerblue3", linetype = "dotted") +
    geom_point(aes(x = year, y = county_pop_z), color = "firebrick") +
    geom_line(aes(x = year, y = county_pop_z), color = "firebrick", linetype = "dotted") +
    labs(x = NULL, y = "log(Yearly Median Daily Device Count), Z-Score Normalized") +
    ggtitle(sewersheds_result[i]) +
    scale_y_continuous(labels = scales::comma,
                       sec.axis = sec_axis(~.*1, name = "Yearly Estimated County Population (Z-Score Normalized)", labels = scales::comma)) +
    theme(plot.title = element_text(hjust = 0.5),
          axis.title.y.left = element_text(color = "dodgerblue3"),
          axis.text.y.left = element_text(color = "dodgerblue3"),
          axis.title.y.right = element_text(color = "firebrick"),
          axis.text.y.right = element_text(color = "firebrick"))
  ggsave(paste0("Figures/Devices vs. Population/", sewersheds_analysis[[i]], ".png"), height = 7, width = 10, devices_vs_pop[[i]])
}

# z-score standardize by year
# yearly_average_norm <- visit_with_flow_log %>%
#   group_by(sewershed_result_name, year) %>%
#   mutate(median_year_z_norm = median(scale(total_devices_daily_log))) %>%
#   slice(1)
# 
# yearly_average_combined <- left_join(yearly_average, yearly_average_norm)

  # geom_line(aes(x = year, y = median_year_z, group = sewershed_result_name),
  #           color = "grey10", alpha = 0.3) +
  # geom_line(aes(x = year, y = county_pop_z, group = sewershed_result_name),
  #           color = "firebrick", alpha = 0.2) +
  # labs(x = NULL, y = "median(Daily Device Count), Z-Score Standardized by Sewershed Only") +
  # ggtitle("All Sewersheds") +
  # scale_y_continuous(labels = scales::comma,
  #                    sec.axis = sec_axis(~.*1, name = "Yearly Estimated County Population, Z-Score Standardized")) +
  # theme(plot.title = element_text(hjust = 0.5),
  #       axis.title.y.right = element_text(color = "firebrick", margin = margin(r = 20)),
  #       axis.text.y.right = element_text(color = "firebrick"))

# spaghetti_yearly_average_norm <- ggplot(data = yearly_average_combined) +
#   geom_line(aes(x = year, y = median_year_z_norm, group = sewershed_result_name),
#             color = "grey10", alpha = 0.3) +
#   geom_line(aes(x = year, y = county_pop_z, group = sewershed_result_name),
#             color = "firebrick", alpha = 0.2) +
#   labs(x = NULL, y = "median(Daily Device Count), Z-Score Standardized by Sewershed and Year") +
#   ggtitle("All Sewersheds") +
#   scale_y_continuous(labels = scales::comma,
#                      sec.axis = sec_axis(~.*1, name = "Yearly Estimated County Population, Z-Score Standardized")) +
#   theme(plot.title = element_text(hjust = 0.5),
#         axis.title.y.right = element_text(color = "firebrick"),
#         axis.text.y.right = element_text(color = "firebrick"))
# 
# spaghetti_yearly_average_combined <- grid.arrange(spaghetti_yearly_average, spaghetti_yearly_average_norm, ncol = 2)
# ggsave("Figures/spaghetti_yearly_average_combined.png", height = 7, width = 20, plot = spaghetti_yearly_average_combined)

visit_with_flow_log_norm <- visit_with_flow_log %>%
  group_by(sewershed_analysis_name, year) %>%
  mutate(total_devices_daily_log_norm = scale(total_devices_daily_log),
         total_devices_monthly_log_norm = scale(total_devices_monthly_log),
         internal_devices_monthly_log_norm = scale(internal_devices_monthly_log),
         flow_rate_log_norm = scale(flow_rate_log),
         # remove the two observations from the two sewersheds that had only one observation in 2021
         flow_rate = ifelse(is.nan(flow_rate_log_norm), NA, flow_rate),
         flow_rate_log_norm = ifelse(is.nan(flow_rate_log_norm), NA, flow_rate_log_norm),
         # remove any outliers that are >3 sd's from the mean
         total_devices_daily_log_norm = case_when(total_devices_daily_log_norm >= -3 & total_devices_daily_log_norm <= 3 ~ total_devices_daily_log_norm,
                                                         TRUE ~ NA_integer_),
         total_devices_monthly_log_norm = case_when(total_devices_monthly_log_norm >= -3 & total_devices_monthly_log_norm <= 3 ~ total_devices_monthly_log_norm,
                                                    TRUE ~ NA_integer_),
         internal_devices_monthly_log_norm = case_when(internal_devices_monthly_log_norm >= -3 & internal_devices_monthly_log_norm <= 3 ~ internal_devices_monthly_log_norm,
                                                       TRUE ~ NA_integer_),
         flow_rate_log_norm = case_when(flow_rate_log_norm >= -3 & flow_rate_log_norm <= 3 ~ flow_rate_log_norm,
                                        TRUE ~ NA_integer_))

# perform K-Means clustering on the 66 sewersheds that have both mobile device data and flow data,
# and reported over a period of 365 days or more (in Python)

# read in K-Means cluster data
clusters <- read.csv("DataProcessed/2025_07_10_cluster_results/monthly/2025_06_25_monthly_kmeans_cluster_labels.csv") %>%
  rename(cdphe_flow_name = AREA_SEWERSHED) %>%
  select(cdphe_flow_name, kmeans_k2)

# merge with clusters and define temporal groupings
visit_with_cluster <- left_join(visit_with_flow_log_norm, clusters) %>%
  arrange(sewershed_result_name) %>%
  filter(date >= min(cdphe_flow$date),
         date <= max(cdphe_flow$date)) %>%
  group_by(sewershed_analysis_name) %>%
  mutate(cluster = ifelse(kmeans_k2 == 1, "Cluster 1", "Cluster 2"),
         month = as.factor(month(date, label = T, abbr = F)),
         month_num = as.numeric(month(date, label = F)),
         season = as.factor(case_when(month %in% c("December", "January", "February") ~ "Winter",
                                      month %in% c("March", "April", "May") ~ "Spring",
                                      month %in% c("June", "July", "August") ~ "Summer",
                                      TRUE ~ "Fall")),
         year = as.factor(year),
         day_of_week = wday(date, label = T, abbr = F),
         day_of_week_num = wday(date, label = F, week_start = 1),
         weekday_weekend = as.factor(case_when(day_of_week %in% c("Monday", "Tuesday", "Wednesday",
                                                                  "Thursday", "Friday") ~ "Weekday",
                                               TRUE ~ "Weekend"))) %>%
  select(-contains("reporting"), -ndays, -kmeans_k2)

save(visit_with_cluster, file = "DataProcessed/visit_with_cluster.Rda")

device_counts <- visit_with_cluster %>%
  replace(is.na(.), 0) %>%
  group_by(sewershed_analysis_name, month) %>%
  mutate(sum_monthly_devices = sum(total_devices_monthly),
         sum_daily_devices = sum(total_devices_daily),
         ratio_daily_monthly = sum_daily_devices/sum_monthly_devices) %>%
  slice(1) %>%
  select(sewershed_analysis_name, cluster, month, sum_monthly_devices, sum_daily_devices, ratio_daily_monthly)

device_compare1 <- ggplot(data = device_counts %>% filter(cluster == "Cluster 1")) +
  geom_line(aes(x = month, y = ratio_daily_monthly, group = sewershed_analysis_name), color = "orange3", alpha = 0.6) +
  labs(x = NULL, y = "Ratio of Summed Daily to Monthly Devices") +
  ggtitle("Cluster 1 Sewersheds (n = 13)") +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 90, hjust = 0.95, vjust = 0.4),
        axis.title.y = element_text(margin = margin(l = 10, r = 5)),
        plot.margin = margin(15, 20, 10, 0, "pt"))

device_compare2 <- ggplot(data = device_counts %>% filter(cluster == "Cluster 2")) +
  geom_line(aes(x = month, y = ratio_daily_monthly, group = sewershed_analysis_name), color = "dodgerblue3", alpha = 0.6) +
  labs(x = NULL, y = "Ratio of Summed Daily to Monthly Devices") +
  ggtitle("Cluster 2 Sewersheds (n = 53)") +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 90, hjust = 0.95, vjust = 0.4),
        axis.title.y = element_text(margin = margin(l = 10, r = 5)),
        plot.margin = margin(15, 20, 10, 0, "pt"))

device_compare <- grid.arrange(device_compare1, device_compare2)
ggsave("Figures/device_compare.png", height = 8, width = 10, plot = device_compare)

# get sample size by mobility cluster
visit_with_cluster %>%
  group_by(sewershed_analysis_name) %>%
  slice(1) %>%
  group_by(cluster) %>%
  count()

# create cluster table
cluster_table <- visit_with_cluster %>%
  group_by(sewershed_result_name) %>%
  slice(1) %>%
  arrange(cluster, sewershed_result_name) %>%
  select(sewershed_result_name, cluster)

write.csv(cluster_table, "DataProcessed/cluster_table.csv")

spaghetti1 <- visit_with_cluster %>%
  filter(cluster == "Cluster 1") %>%
  group_by(sewershed_analysis_name, year) %>%
  drop_na(total_devices_monthly_log) %>%
  mutate(spaghetti_value = scale(total_devices_monthly_log))

spaghetti2 <- visit_with_cluster %>%
  filter(cluster == "Cluster 2") %>%
  group_by(sewershed_analysis_name, year) %>%
  drop_na(total_devices_monthly_log) %>%
  mutate(spaghetti_value = scale(total_devices_monthly_log))

spaghetti_plot1 <- ggplot(spaghetti1) +
  geom_line(aes(x = date, y = spaghetti_value, group = sewershed_analysis_name), color = "orange3", alpha = 0.5) +
  labs(x = NULL, y = "Normalized Monthly Device Count") +
  ggtitle("Cluster 1 Sewersheds (n = 13)") +
  scale_x_date(date_breaks = "4 months", date_labels = "%b %Y") +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 90, hjust = 0.95, vjust = 0.4),
        axis.title.y = element_text(margin = margin(l = 10, r = 5)),
        plot.margin = margin(10, 20, 10, 0, "pt"))

spaghetti_plot2 <- ggplot(spaghetti2) +
  geom_line(aes(x = date, y = spaghetti_value, group = sewershed_analysis_name), color = "dodgerblue3", alpha = 0.5) +
  labs(x = NULL, y = "Normalized Monthly Device Count") +
  ggtitle("Cluster 2 Sewersheds (n = 53)") +
  scale_x_date(date_breaks = "4 months", date_labels = "%b %Y") +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 90, hjust = 0.95, vjust = 0.4),
        axis.title.y = element_text(margin = margin(l = 10, r = 5)),
        plot.margin = margin(10, 20, 10, 0, "pt"))

spaghetti_cluster <- grid.arrange(spaghetti_plot1, spaghetti_plot2, ncol = 1)

# save spaghetti plot
ggsave("Figures/spaghetti_cluster.png", height = 8, width = 10, plot = spaghetti_cluster)

# sensitivity analysis using daily data
spaghetti_daily1 <- visit_with_cluster %>%
  filter(cluster == "Cluster 1") %>%
  group_by(sewershed_analysis_name, year) %>%
  mutate(spaghetti_value = scale(total_devices_daily_log))

spaghetti_daily2 <- visit_with_cluster %>%
  filter(cluster == "Cluster 2") %>%
  group_by(sewershed_analysis_name, year) %>%
  mutate(spaghetti_value = scale(total_devices_daily_log))

spaghetti_plot_daily1 <- ggplot(spaghetti_daily1) +
  geom_line(aes(x = date, y = spaghetti_value, group = sewershed_analysis_name), color = "orange3", alpha = 0.5) +
  labs(x = NULL, y = "Normalized Daily Device Count") +
  ggtitle("Cluster 1 Sewersheds (n = 13)") +
  scale_x_date(date_breaks = "4 months", date_labels = "%b %Y") +
  scale_y_continuous(limits = c(-10, 10)) +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 90, hjust = 0.95, vjust = 0.4),
        axis.title.y = element_text(margin = margin(l = 10, r = 5)),
        plot.margin = margin(10, 20, 10, 0, "pt"))

spaghetti_plot_daily2 <- ggplot(spaghetti_daily2) +
  geom_line(aes(x = date, y = spaghetti_value, group = sewershed_analysis_name), color = "dodgerblue3", alpha = 0.5) +
  labs(x = NULL, y = "Normalized Daily Device Count") +
  ggtitle("Cluster 2 Sewersheds (n = 53)") +
  scale_x_date(date_breaks = "4 months", date_labels = "%b %Y") +
  scale_y_continuous(limits = c(-10, 10)) +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 90, hjust = 0.95, vjust = 0.4),
        axis.title.y = element_text(margin = margin(l = 10, r = 5)),
        plot.margin = margin(10, 20, 10, 0, "pt"))

spaghetti_cluster_daily <- grid.arrange(spaghetti_plot_daily1, spaghetti_plot_daily2, ncol = 1)

# save spaghetti plot
ggsave("Figures/spaghetti_cluster_daily.png", height = 8, width = 10, plot = spaghetti_cluster_daily)

# make flow reporting timeline
timeline <- visit_with_cluster %>%
  group_by(cdphe_flow_name) %>%
  # forward fill flow rate to identify start date
  fill(flow_rate, .direction = "down") %>%
  drop_na(flow_rate) %>%
  mutate(start_date = min(date)) %>%
  group_by(cdphe_flow_name, start_date) %>%
  arrange(desc(start_date)) %>%
  ungroup() %>%
  mutate(id = match(cdphe_flow_name, unique(cdphe_flow_name))) %>%
  select(date, start_date, sewershed_result_name, cluster, id)

# plot flow reporting timeline
flow_reporting_timeline <- ggplot(data = timeline) +
  geom_line(aes(x = date, y = sewershed_result_name, group = id, color = cluster), linewidth = 0.8) +
  labs(x = "Reporting Date Range", y = NULL) +
  scale_x_date(date_breaks = "2 months", date_labels = "%b %Y") +
  scale_y_discrete(limits = unique(timeline$sewershed_result_name)) +
  scale_color_manual(values = c("orange3", "dodgerblue3"),
                     labels = c("Cluster 1", "Cluster 2")) +
  guides(color = guide_legend(title = "Mobility Cluster")) +
  theme(plot.title = element_text(hjust = 0.5),
        axis.title.x = element_text(margin = margin(t = 5)),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 0.95, vjust = 0.4)); flow_reporting_timeline

ggsave("Figures/flow_reporting_timeline.png", width = 9, height = 11, plot = flow_reporting_timeline)

# reorder levels for plotting
visit_with_cluster$season <- factor(visit_with_cluster$season, levels = c("Winter", "Spring", "Summer", "Fall"))
visit_with_cluster$day_of_week <- factor(visit_with_cluster$day_of_week, levels = c("Monday", "Tuesday", "Wednesday", "Thursday", "Friday", "Saturday", "Sunday"))

# mean device count by day of week
# visit_daily_stat_dayofweek <- visit_with_cluster %>%
#   group_by(cluster, day_of_week) %>%
#   mutate(mean = mean(total_devices_daily_log_norm),
#          sd = sd(total_devices_daily_log_norm)) %>%
#   slice(1)

density_dayofweek1 <- ggplot(data = visit_with_cluster %>% filter(cluster == "Cluster 1")) +
  stat_slabinterval(aes(x = day_of_week, y = total_devices_daily_log_norm),
               # set point interval (default width is c(.66, .95)
               # that means the thick line contains 66% of the density
               # and the thin line contains 95% of the density
               normalize = "none", point_interval = "mean_qi", .width = c(0.66, 0.66), fill = "orange3",
               point_color = "sienna4", interval_color = "sienna4", alpha = 0.6) +
  labs(x = NULL, y = "Normalized Daily Device Count") +
  ggtitle("Cluster 1 Sewersheds (n = 13)") +
  scale_y_continuous(limits = c(-3, 3)) +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 90, hjust = 0.95, vjust = 0.4),
        axis.title.y = element_text(margin = margin(l = 10, r = 5)),
        plot.margin = margin(15, 20, 10, 0, "pt"))

density_dayofweek2 <- ggplot(data = visit_with_cluster %>% filter(cluster == "Cluster 2")) +
  stat_slabinterval(aes(x = day_of_week, y = total_devices_daily_log_norm),
               # set point interval (default width is c(.66, .95)
               # that means the thick line contains 66% of the density
               # and the thin line contains 95% of the density
               normalize = "none", point_interval = "mean_qi", .width = c(0.66, 0.66), fill = "dodgerblue3",
               point_color = "midnightblue", interval_color = "midnightblue", alpha = 0.6) +
  labs(x = NULL, y = "Normalized Daily Device Count") +
  ggtitle("Cluster 2 Sewersheds (n = 53)") +
  scale_y_continuous(limits = c(-3, 3)) +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 90, hjust = 0.95, vjust = 0.4),
        axis.title.y = element_text(margin = margin(l = 10, r = 5)),
        plot.margin = margin(15, 20, 10, 0, "pt"))

density_dayofweek <- grid.arrange(density_dayofweek1, density_dayofweek2)

ggsave("Figures/density_dayofweek.png", height = 9, width = 9, plot = density_dayofweek)

density_wdwe1 <- ggplot(data = visit_with_cluster %>% filter(cluster == "Cluster 1")) +
  stat_slabinterval(aes(x = weekday_weekend, y = total_devices_daily_log_norm),
                    # set point interval (default width is c(.66, .95)
                    # that means the thick line contains 66% of the density
                    # and the thin line contains 95% of the density
                    normalize = "none", point_interval = "mean_qi", .width = c(0.66, 0.66), fill = "orange3",
                    point_color = "sienna4", interval_color = "sienna4", alpha = 0.6) +
  labs(x = NULL, y = "Normalized Daily Device Count") +
  ggtitle("Cluster 1 Sewersheds (n = 13)") +
  scale_y_continuous(limits = c(-3, 3)) +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 90, hjust = 0.95, vjust = 0.4),
        axis.title.y = element_text(margin = margin(l = 10, r = 5)),
        plot.margin = margin(15, 20, 10, 0, "pt"))

density_wdwe2 <- ggplot(data = visit_with_cluster %>% filter(cluster == "Cluster 2")) +
  stat_slabinterval(aes(x = weekday_weekend, y = total_devices_daily_log_norm),
                    # set point interval (default width is c(.66, .95)
                    # that means the thick line contains 66% of the density
                    # and the thin line contains 95% of the density
                    normalize = "none", point_interval = "mean_qi", .width = c(0.66, 0.66), fill = "dodgerblue3",
                    point_color = "midnightblue", interval_color = "midnightblue", alpha = 0.6) +
  labs(x = NULL, y = "Normalized Daily Device Count") +
  ggtitle("Cluster 2 Sewersheds (n = 53)") +
  scale_y_continuous(limits = c(-3, 3)) +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 90, hjust = 0.95, vjust = 0.4),
        axis.title.y = element_text(margin = margin(l = 10, r = 5)),
        plot.margin = margin(15, 20, 10, 0, "pt"))

density_wdwe <- grid.arrange(density_wdwe1, density_wdwe2)

ggsave("Figures/density_wdwe.png", height = 9, width = 9, plot = density_wdwe)

density_month1 <- ggplot(data = visit_with_cluster %>% filter(cluster == "Cluster 1")) +
  stat_slabinterval(aes(x = month, y = total_devices_monthly_log_norm),
               # set point interval (default width is c(.66, .95)
               # that means the thick line contains 66% of the density
               # and the thin line contains 95% of the density
               # for now, make it .66 to encompass approximately 1 sd
               normalize = "none", point_interval = "mean_qi", .width = c(0.66, 0.66), fill = "orange3",
               point_color = "sienna4", interval_color = "sienna4", alpha = 0.6) +
  stat_slabinterval(aes(x = month, y = flow_rate_log_norm),
                    normalize = "none", point_interval = "mean_qi", .width = c(0.66, 0.66),
                    position = position_nudge(x = 0.15), fill = "grey20", alpha = 0.6) +
  labs(x = NULL, y = "Normalized Monthly Device Count") +
  ggtitle("Cluster 1 Sewersheds (n = 13)") +
  scale_y_continuous(limits = c(-3, 3), sec.axis = sec_axis(~.*1, name = "Normalized Monthly Flow Rate")) +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 90, hjust = 0.95, vjust = 0.4),
        axis.title.y.left = element_text(color = "orange3", margin = margin(l = 10, r = 5)),
        axis.text.y.left = element_text(color = "orange3"),
        plot.margin = margin(15, 20, 10, 0, "pt"))
  
density_month2 <- ggplot(data = visit_with_cluster %>% filter(cluster == "Cluster 2")) +
  stat_slabinterval(aes(x = month, y = total_devices_monthly_log_norm),
               # set point interval (default width is c(.66, .95)
               # that means the thick line contains 66% of the density
               # and the thin line contains 95% of the density
               # for now, make it .66 to encompass approximately 1 sd
               normalize = "none", point_interval = "mean_qi", .width = c(0.66, 0.66), fill = "dodgerblue3",
               point_color = "midnightblue", interval_color = "midnightblue", alpha = 0.6) +
  stat_slabinterval(aes(x = month, y = flow_rate_log_norm),
                    normalize = "none", point_interval = "mean_qi", .width = c(0.66, 0.66),
                    position = position_nudge(x = 0.15), fill = "grey20", alpha = 0.6) +
  labs(x = NULL, y = "Normalized Monthly Device Count") +
  ggtitle("Cluster 2 Sewersheds (n = 53)") +
  scale_y_continuous(limits = c(-3, 3), sec.axis = sec_axis(~.*1, name = "Normalized Monthly Flow Rate")) +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 90, hjust = 0.95, vjust = 0.4),
        axis.title.y.left = element_text(color = "dodgerblue3", margin = margin(l = 10, r = 5)),
        axis.text.y.left = element_text(color = "dodgerblue3"),
        plot.margin = margin(15, 20, 10, 0, "pt"))

density_month <- grid.arrange(density_month1, density_month2)

ggsave("Figures/density_month.png", height = 9, width = 12, plot = density_month)

density_season1 <- ggplot(data = visit_with_cluster %>% filter(cluster == "Cluster 1")) +
  stat_slabinterval(aes(x = season, y = total_devices_monthly_log_norm),
                    # set point interval (default width is c(.66, .95)
                    # that means the thick line contains 66% of the density
                    # and the thin line contains 95% of the density
                    # for now, make it .66 to encompass approximately 1 sd
                    normalize = "none", point_interval = "mean_qi", .width = c(0.66, 0.66), fill = "orange3",
                    point_color = "sienna4", interval_color = "sienna4", alpha = 0.6) +
  stat_slabinterval(aes(x = season, y = flow_rate_log_norm),
                    normalize = "none", point_interval = "mean_qi", .width = c(0.66, 0.66),
                    position = position_nudge(x = 0.15), fill = "grey20", alpha = 0.6) +
  labs(x = NULL, y = "Normalized Monthly Device Count") +
  ggtitle("Cluster 1 Sewersheds (n = 13)") +
  scale_y_continuous(limits = c(-3, 3), sec.axis = sec_axis(~.*1, name = "Normalized Monthly Flow Rate")) +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 90, hjust = 0.95, vjust = 0.4),
        axis.title.y.left = element_text(color = "orange3", margin = margin(l = 10, r = 5)),
        axis.text.y.left = element_text(color = "orange3"),
        plot.margin = margin(15, 20, 10, 0, "pt"))

density_season2 <- ggplot(data = visit_with_cluster %>% filter(cluster == "Cluster 2")) +
  stat_slabinterval(aes(x = season, y = total_devices_monthly_log_norm),
                    # set point interval (default width is c(.66, .95)
                    # that means the thick line contains 66% of the density
                    # and the thin line contains 95% of the density
                    # for now, make it .66 to encompass approximately 1 sd
                    normalize = "none", point_interval = "mean_qi", .width = c(0.66, 0.66), fill = "dodgerblue3",
                    point_color = "midnightblue", interval_color = "midnightblue", alpha = 0.6) +
  stat_slabinterval(aes(x = season, y = flow_rate_log_norm),
                    normalize = "none", point_interval = "mean_qi", .width = c(0.66, 0.66),
                    position = position_nudge(x = 0.15), fill = "grey20", alpha = 0.6) +
  labs(x = NULL, y = "Normalized Monthly Device Count") +
  ggtitle("Cluster 2 Sewersheds (n = 53)") +
  scale_y_continuous(limits = c(-3, 3), sec.axis = sec_axis(~.*1, name = "Normalized Monthly Flow Rate")) +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 90, hjust = 0.95, vjust = 0.4),
        axis.title.y.left = element_text(color = "dodgerblue3", margin = margin(l = 10, r = 5)),
        axis.text.y.left = element_text(color = "dodgerblue3"),
        plot.margin = margin(15, 20, 10, 0, "pt"))

density_season <- grid.arrange(density_season1, density_season2)

ggsave("Figures/density_season.png", height = 9, width = 12, plot = density_season)

  
# calculate mean and sd for device counts
calc_dayofweek <- visit_with_cluster %>%
  group_by(cluster, day_of_week) %>%
  summarize(mean = round(mean(total_devices_daily_log_norm), 3),
            sd = round(sd(total_devices_daily_log_norm), 3)) %>%
  arrange(cluster, factor(day_of_week, levels = c("Monday", "Tuesday", "Wednesday", "Thursday", "Friday", "Saturday", "Sunday")))
write.csv(calc_dayofweek, "DataProcessed/calc_dayofweek.csv")

calc_wdwe <- visit_with_cluster %>%
  group_by(cluster, weekday_weekend) %>%
  summarize(mean = round(mean(total_devices_daily_log_norm), 3),
            sd = round(sd(total_devices_daily_log_norm), 3))
write.csv(calc_wdwe, "DataProcessed/calc_wdwe.csv")

calc_month <- visit_with_cluster %>%
  drop_na(total_devices_monthly) %>%
  group_by(cluster, month) %>%
  summarize(mean = round(mean(total_devices_monthly_log_norm), 3),
            sd = round(sd(total_devices_monthly_log_norm), 3))
write.csv(calc_month, "DataProcessed/calc_month.csv")

calc_season <- visit_with_cluster %>%
  drop_na(total_devices_monthly) %>%
  group_by(cluster, season) %>%
  summarize(mean = round(mean(total_devices_monthly_log_norm), 3),
            sd = round(sd(total_devices_monthly_log_norm), 3))
write.csv(calc_season, "DataProcessed/calc_season.csv")

# linear mixed effect models and LRT on log-transformed device counts
mod_dayofweek1 <- lmer(total_devices_daily_log_norm ~ 1 + (1 | sewershed_analysis_name), data = visit_with_cluster %>% filter(cluster == "Cluster 1"), REML = F)
mod1_dayofweek1 <- lmer(total_devices_daily_log_norm ~ 1 + day_of_week + (1 | sewershed_analysis_name), data = visit_with_cluster %>% filter(cluster == "Cluster 1"), REML = F)
anova(mod_dayofweek1, mod1_dayofweek1, test = "LRT")

mod_dayofweek2 <- lmer(total_devices_daily_log_norm ~ 1 + (1 | sewershed_analysis_name), data = visit_with_cluster %>% filter(cluster == "Cluster 2"), REML = F)
mod1_dayofweek2 <- lmer(total_devices_daily_log_norm ~ 1 + day_of_week + (1 | sewershed_analysis_name), data = visit_with_cluster %>% filter(cluster == "Cluster 2"), REML = F)
anova(mod_dayofweek2, mod1_dayofweek2, test = "LRT")

mod_wdwe1 <- lmer(total_devices_daily_log_norm ~ 1 + (1 | sewershed_analysis_name), data = visit_with_cluster %>% filter(cluster == "Cluster 1"), REML = F)
mod1_wdwe1 <- lmer(total_devices_daily_log_norm ~ 1 + weekday_weekend + (1 | sewershed_analysis_name), data = visit_with_cluster %>% filter(cluster == "Cluster 1"), REML = F)
anova(mod_wdwe1, mod1_wdwe1, test = "LRT")

mod_wdwe2 <- lmer(total_devices_daily_log_norm ~ 1 + (1 | sewershed_analysis_name), data = visit_with_cluster %>% filter(cluster == "Cluster 2"), REML = F)
mod1_wdwe2 <- lmer(total_devices_daily_log_norm ~ 1 + weekday_weekend + (1 | sewershed_analysis_name), data = visit_with_cluster %>% filter(cluster == "Cluster 2"), REML = F)
anova(mod_wdwe2, mod1_wdwe2, test = "LRT")

mod_month1 <- lmer(total_devices_monthly_log_norm ~ 1 + (1 | sewershed_analysis_name), data = visit_with_cluster %>% filter(cluster == "Cluster 1"), REML = F)
mod1_month1 <- lmer(total_devices_monthly_log_norm ~ 1 + month + (1 | sewershed_analysis_name), data = visit_with_cluster %>% filter(cluster == "Cluster 1"), REML = F)
anova(mod_month1, mod1_month1, test = "LRT")

mod_month2 <- lmer(total_devices_monthly_log_norm ~ 1 + (1 | sewershed_analysis_name), data = visit_with_cluster %>% filter(cluster == "Cluster 2"), REML = F)
mod1_month2 <- lmer(total_devices_monthly_log_norm ~ 1 + month + (1 | sewershed_analysis_name), data = visit_with_cluster %>% filter(cluster == "Cluster 2"), REML = F)
anova(mod_month2, mod1_month2, test = "LRT")

mod_season1 <- lmer(total_devices_monthly_log_norm ~ 1 + (1 | sewershed_analysis_name), data = visit_with_cluster %>% filter(cluster == "Cluster 1"), REML = F)
mod1_season1 <- lmer(total_devices_monthly_log_norm ~ 1 + season + (1 | sewershed_analysis_name), data = visit_with_cluster %>% filter(cluster == "Cluster 1"), REML = F)
anova(mod_season1, mod1_season1, test = "LRT")

mod_season2 <- lmer(total_devices_monthly_log_norm ~ 1 + (1 | sewershed_analysis_name), data = visit_with_cluster %>% filter(cluster == "Cluster 2"), REML = F)
mod1_season2 <- lmer(total_devices_monthly_log_norm ~ 1 + season + (1 | sewershed_analysis_name), data = visit_with_cluster %>% filter(cluster == "Cluster 2"), REML = F)
anova(mod_season2, mod1_season2, test = "LRT")

# sensitivity analysis on internal devices
calc_month_internal <- visit_with_cluster %>%
  drop_na(total_devices_monthly) %>%
  group_by(cluster, month) %>%
  summarize(mean = round(mean(internal_devices_monthly_log_norm), 3),
            sd = round(sd(internal_devices_monthly_log_norm), 3))
write.csv(calc_month_internal, "DataProcessed/calc_month_internal.csv")

calc_season_internal <- visit_with_cluster %>%
  drop_na(total_devices_monthly) %>%
  group_by(cluster, season) %>%
  summarize(mean = round(mean(internal_devices_monthly_log_norm), 3),
            sd = round(sd(internal_devices_monthly_log_norm), 3)) %>%
  arrange(cluster, factor(season, levels = c("Winter", "Spring", "Summer", "Fall")))
write.csv(calc_season_internal, "DataProcessed/calc_season_internal.csv")

mod_month_internal1 <- lmer(internal_devices_monthly_log_norm ~ 1 + (1 | sewershed_analysis_name), data = visit_with_cluster %>% filter(cluster == "Cluster 1"), REML = F)
mod1_month_internal1 <- lmer(internal_devices_monthly_log_norm ~ 1 + month + (1 | sewershed_analysis_name), data = visit_with_cluster %>% filter(cluster == "Cluster 1"), REML = F)
anova(mod_month_internal1, mod1_month_internal1, test = "LRT")

mod_month_internal2 <- lmer(internal_devices_monthly_log_norm ~ 1 + (1 | sewershed_analysis_name), data = visit_with_cluster %>% filter(cluster == "Cluster 2"), REML = F)
mod1_month_internal2 <- lmer(internal_devices_monthly_log_norm ~ 1 + month + (1 | sewershed_analysis_name), data = visit_with_cluster %>% filter(cluster == "Cluster 2"), REML = F)
anova(mod_month_internal2, mod1_month_internal2, test = "LRT")

mod_season_internal1 <- lmer(internal_devices_monthly_log_norm ~ 1 + (1 | sewershed_analysis_name), data = visit_with_cluster %>% filter(cluster == "Cluster 1"), REML = F)
mod1_season_internal1 <- lmer(internal_devices_monthly_log_norm ~ 1 + season + (1 | sewershed_analysis_name), data = visit_with_cluster %>% filter(cluster == "Cluster 1"), REML = F)
anova(mod_season_internal1, mod1_season_internal1, test = "LRT")

mod_season_internal2 <- lmer(internal_devices_monthly_log_norm ~ 1 + (1 | sewershed_analysis_name), data = visit_with_cluster %>% filter(cluster == "Cluster 2"), REML = F)
mod1_season_internal2 <- lmer(internal_devices_monthly_log_norm ~ 1 + season + (1 | sewershed_analysis_name), data = visit_with_cluster %>% filter(cluster == "Cluster 2"), REML = F)
anova(mod_season_internal2, mod1_season_internal2, test = "LRT")

# calculate mean and sd for flow rates
calc_month_flow <- visit_with_cluster %>%
  drop_na(flow_rate) %>%
  group_by(cluster, month) %>%
  summarize(mean = signif(mean(flow_rate_log_norm), 3),
            sd = signif(sd(flow_rate_log_norm), 3))
write.csv(calc_month_flow, "DataProcessed/calc_month_flow.csv")

calc_season_flow <- visit_with_cluster %>%
  drop_na(flow_rate) %>%
  group_by(cluster, season) %>%
  summarize(mean = signif(mean(flow_rate_log_norm), 3),
            sd = signif(sd(flow_rate_log_norm), 3)) %>%
  arrange(cluster, factor(season, levels = c("Winter", "Spring", "Summer", "Fall")))
write.csv(calc_season_flow, "DataProcessed/calc_season_flow.csv")

mod_month_flow1 <-  lmer(flow_rate_log_norm ~ 1 +               (1 | sewershed_analysis_name), data = visit_with_cluster %>% filter(cluster == "Cluster 1"), REML = F)
mod1_month_flow1 <- lmer(flow_rate_log_norm ~ 1 + month + (1 | sewershed_analysis_name), data = visit_with_cluster %>% filter(cluster == "Cluster 1"), REML = F)
anova(mod_month_flow1, mod1_month_flow1, test = "LRT")

mod_month_flow2 <- lmer(flow_rate_log_norm ~ 1 + (1 | sewershed_analysis_name), data = visit_with_cluster %>% filter(cluster == "Cluster 2"), REML = F)
mod1_month_flow2 <- lmer(flow_rate_log_norm ~ 1 + month + (1 | sewershed_analysis_name), data = visit_with_cluster %>% filter(cluster == "Cluster 2"), REML = F)
anova(mod_month_flow2, mod1_month_flow2, test = "LRT")

mod_season_flow1 <-  lmer(flow_rate_log_norm ~ 1 +               (1 | sewershed_analysis_name), data = visit_with_cluster %>% filter(cluster == "Cluster 1"), REML = F)
mod1_season_flow1 <- lmer(flow_rate_log_norm ~ 1 + season + (1 | sewershed_analysis_name), data = visit_with_cluster %>% filter(cluster == "Cluster 1"), REML = F)
anova(mod_season_flow1, mod1_season_flow1, test = "LRT")

mod_season_flow2 <- lmer(flow_rate_log_norm ~ 1 + (1 | sewershed_analysis_name), data = visit_with_cluster %>% filter(cluster == "Cluster 2"), REML = F)
mod1_season_flow2 <- lmer(flow_rate_log_norm ~ 1 + season + (1 | sewershed_analysis_name), data = visit_with_cluster %>% filter(cluster == "Cluster 2"), REML = F)
anova(mod_season_flow2, mod1_season_flow2, test = "LRT")

########################################################################################################
