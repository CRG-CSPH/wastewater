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

#read in monthly device count data
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

# z-score normalize by year
yearly_average_norm <- visit_with_flow_log %>%
  group_by(sewershed_result_name, year) %>%
  mutate(median_year_z_norm = median(scale(total_devices_daily_log))) %>%
  slice(1)

yearly_average_combined <- left_join(yearly_average, yearly_average_norm)

spaghetti_yearly_average <- ggplot(data = yearly_average_combined) +
  geom_line(aes(x = year, y = median_year_z, group = sewershed_result_name),
            color = "grey10", alpha = 0.3) +
  geom_line(aes(x = year, y = county_pop_z, group = sewershed_result_name),
            color = "firebrick", alpha = 0.2) +
  labs(x = NULL, y = "median(Daily Device Count), Z-Score Standardized by Sewershed Only") +
  ggtitle("All Sewersheds") +
  scale_y_continuous(labels = scales::comma,
                     sec.axis = sec_axis(~.*1, name = "Yearly Estimated County Population, Z-Score Standardized")) +
  theme(plot.title = element_text(hjust = 0.5),
        axis.title.y.right = element_text(color = "firebrick", margin = margin(r = 20)),
        axis.text.y.right = element_text(color = "firebrick"))

spaghetti_yearly_average_norm <- ggplot(data = yearly_average_combined) +
  geom_line(aes(x = year, y = median_year_z_norm, group = sewershed_result_name),
            color = "grey10", alpha = 0.3) +
  geom_line(aes(x = year, y = county_pop_z, group = sewershed_result_name),
            color = "firebrick", alpha = 0.2) +
  labs(x = NULL, y = "median(Daily Device Count), Z-Score Standardized by Sewershed and Year") +
  ggtitle("All Sewersheds") +
  scale_y_continuous(labels = scales::comma,
                     sec.axis = sec_axis(~.*1, name = "Yearly Estimated County Population, Z-Score Standardized")) +
  theme(plot.title = element_text(hjust = 0.5),
        axis.title.y.right = element_text(color = "firebrick"),
        axis.text.y.right = element_text(color = "firebrick"))

spaghetti_yearly_average_combined <- grid.arrange(spaghetti_yearly_average, spaghetti_yearly_average_norm, ncol = 2)
ggsave("Figures/spaghetti_yearly_average_combined.png", height = 7, width = 20, plot = spaghetti_yearly_average_combined)


# compare yearly median device counts (z-score normalized by sewershed only) to
# yearly median device counts (z-score normalized by both sewershed and year)
# spaghetti_yearly_average <- ggplot(data = yearly_average) +
#   geom_line(aes(x = year, y = median_year_z,
#                   group = sewershed_result_name, color = sewershed_result_name)) +
#   labs(x = NULL, y = "log(Yearly Median Daily Device Count), Z-Score Normalized") +
#   guides(color = "none") +
#   scale_color_viridis(option = "viridis", discrete = T) +
#   theme(plot.title = element_text(hjust = 0.5))
# 
# spaghetti_yearly_average_norm <- ggplot(data = yearly_average_norm) +
#   geom_line(aes(x = year, y = median_year_z_norm,
#                 group = sewershed_result_name, color = sewershed_result_name)) +
#   labs(x = NULL, y = "log(Yearly Median Daily Device Count), Z-Score Normalized by Year") +
#   guides(color = "none") +
#   scale_color_viridis(option = "viridis", discrete = T) +
#   theme(plot.title = element_text(hjust = 0.5))
# 
# spaghetti_yearly_average_combined <- grid.arrange(spaghetti_yearly_average, spaghetti_yearly_average_norm, ncol = 2)
# ggsave("Figures/spaghetti_yearly_average_combined.png", height = 7, width = 20, plot = spaghetti_yearly_average_combined)

# to account for changing opt-in rates and data collection practices, z-score normalize log-transformed device counts by sewershed and year
# for consistency, do the same for flow rates

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

# perform K-Means clustering on the 66 sewersheds that have both mobile device data, flow data,
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
  mutate(cluster = ifelse(kmeans_k2 == 1, "High-Variation", "Low-Variation"),
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

spaghetti_high <- visit_with_cluster %>%
  filter(cluster == "High-Variation") %>%
  group_by(sewershed_analysis_name, year) %>%
  drop_na(total_devices_monthly_log) %>%
  mutate(spaghetti_value = scale(total_devices_monthly_log))

spaghetti_low <- visit_with_cluster %>%
  filter(cluster == "Low-Variation") %>%
  group_by(sewershed_analysis_name, year) %>%
  drop_na(total_devices_monthly_log) %>%
  mutate(spaghetti_value = scale(total_devices_monthly_log))

spaghetti_plot_high <- ggplot(spaghetti_high) +
  geom_line(aes(x = date, y = spaghetti_value, group = sewershed_analysis_name), color = "orange3", alpha = 0.5) +
            #color = "orange3", alpha = 0.6) +
  labs(x = NULL, y = "log(Monthly Device Count), Year-Normalized") +
  ggtitle("High-Variation Sewersheds (n = 13)") +
  scale_x_date(date_breaks = "4 months", date_labels = "%b %Y") +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 90, hjust = 0.95, vjust = 0.4),
        axis.title.y = element_text(margin = margin(l = 10, r = 5)),
        plot.margin = margin(10, 20, 10, 0, "pt"))

spaghetti_plot_low <- ggplot(spaghetti_low) +
  geom_line(aes(x = date, y = spaghetti_value, group = sewershed_analysis_name), color = "dodgerblue3", alpha = 0.5) +
            #color = "darkorchid4", alpha = 0.6) +
  labs(x = NULL, y = "log(Monthly Device Count), Year-Normalized") +
  ggtitle("Low-Variation Sewersheds (n = 53)") +
  scale_x_date(date_breaks = "4 months", date_labels = "%b %Y") +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 90, hjust = 0.95, vjust = 0.4),
        axis.title.y = element_text(margin = margin(l = 10, r = 5)),
        plot.margin = margin(10, 20, 10, 0, "pt"))

spaghetti_cluster <- grid.arrange(spaghetti_plot_high, spaghetti_plot_low, ncol = 1)

# save spaghetti plot
ggsave("Figures/spaghetti_cluster.png", height = 8, width = 10, plot = spaghetti_cluster)

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

# sample size for reporting timeline is 66 sewersheds that have flow data, cluster
# information, and at least one year of data
#unique(unlist(cdphe_flow$sewershed_result_name))

# get sample size by mobility cluster
# cdphe_flow %>%
#   group_by(sewershed_analysis_name) %>%
#   slice(1) %>%
#   group_by(cluster) %>%
#   count()

# plot flow reporting timeline
flow_reporting_timeline <- ggplot(data = timeline) +
  geom_line(aes(x = date, y = sewershed_result_name, group = id, color = cluster), linewidth = 0.8) +
  labs(x = "Reporting Date Range", y = NULL) +
  scale_x_date(date_breaks = "2 months", date_labels = "%b %Y") +
  scale_y_discrete(limits = unique(timeline$sewershed_result_name)) +
  scale_color_manual(values = c("orange3", "dodgerblue3"),
                     labels = c("High-Variation", "Low-Variation")) +
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
visit_daily_stat_dayofweek <- visit_with_cluster %>%
  group_by(cluster, day_of_week) %>%
  mutate(mean = mean(total_devices_daily_log_norm),#median(daily_devices_log_scale),
         sd = sd(total_devices_daily_log_norm)) %>%
  slice(1)

density_high_dayofweek <- ggplot(data = visit_with_cluster %>% filter(cluster == "High-Variation")) +
  stat_slabinterval(aes(x = day_of_week, y = total_devices_daily_log_norm),
               # set point interval (default width is c(.66, .95)
               # that means the thick line contains 66% of the density
               # and the thin line contains 95% of the density
               normalize = "none", point_interval = "mean_qi", .width = c(0.66, 0.66), fill = "orange3",
               point_color = "sienna4", interval_color = "sienna4", alpha = 0.6) +
  labs(x = NULL, y = "log(Daily Device Count), Sewershed- and Year-Standardized") +
  ggtitle("High-Variation Sewersheds (n = 13)") +
  scale_y_continuous(limits = c(-3, 3)) +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 90, hjust = 0.95, vjust = 0.4),
        axis.title.y = element_text(margin = margin(l = 10, r = 5)),
        plot.margin = margin(15, 20, 10, 0, "pt"))

density_low_dayofweek <- ggplot(data = visit_with_cluster %>% filter(cluster == "Low-Variation")) +
  stat_slabinterval(aes(x = day_of_week, y = total_devices_daily_log_norm),
               # set point interval (default width is c(.66, .95)
               # that means the thick line contains 66% of the density
               # and the thin line contains 95% of the density
               normalize = "none", point_interval = "mean_qi", .width = c(0.66, 0.66), fill = "dodgerblue3",
               point_color = "midnightblue", interval_color = "midnightblue", alpha = 0.6) +
  labs(x = NULL, y = "log(Daily Device Count), Sewershed- and Year-Standardized") +
  ggtitle("Low-Variation Sewersheds (n = 53)") +
  scale_y_continuous(limits = c(-3, 3)) +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 90, hjust = 0.95, vjust = 0.4),
        axis.title.y = element_text(margin = margin(l = 10, r = 5)),
        plot.margin = margin(15, 20, 10, 0, "pt"))

density_dayofweek <- grid.arrange(density_high_dayofweek, density_low_dayofweek)

ggsave("Figures/density_dayofweek.png", height = 9, width = 9, plot = density_dayofweek)

density_high_month <- ggplot(data = visit_with_cluster %>% filter(cluster == "High-Variation")) +
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
  labs(x = NULL, y = "log(Daily Device Count), Sewershed- and Year-Standardized") +
  ggtitle("High-Variation Sewersheds (n = 13)") +
  scale_y_continuous(limits = c(-3, 3), sec.axis = sec_axis(~.*1, name = "log(Monthly Flow Rate), Sewershed- and Year-Standardized")) +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 90, hjust = 0.95, vjust = 0.4),
        axis.title.y.left = element_text(color = "orange3", margin = margin(l = 10, r = 5)),
        axis.text.y.left = element_text(color = "orange3"),
        plot.margin = margin(15, 20, 10, 0, "pt"))
  
density_low_month <- ggplot(data = visit_with_cluster %>% filter(cluster == "Low-Variation")) +
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
  labs(x = NULL, y = "log(Daily Device Count), Sewershed- and Year-Standardized") +
  ggtitle("Low-Variation Sewersheds (n = 53)") +
  scale_y_continuous(limits = c(-3, 3), sec.axis = sec_axis(~.*1, name = "log(Monthly Flow Rate), Sewershed- and Year-Standardized")) +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 90, hjust = 0.95, vjust = 0.4),
        axis.title.y.left = element_text(color = "dodgerblue3", margin = margin(l = 10, r = 5)),
        axis.text.y.left = element_text(color = "dodgerblue3"),
        plot.margin = margin(15, 20, 10, 0, "pt"))

density_month <- grid.arrange(density_high_month, density_low_month)

ggsave("Figures/density_month.png", height = 9, width = 12, plot = density_month)
  
# median device counts and flow rates by month
visit_monthly_stat_month <- visit_with_cluster %>%
  drop_na(total_devices_monthly_log_norm) %>%
  group_by(cluster, month) %>%
  mutate(mean_total = mean(total_devices_monthly_log_norm),
         sd_total = sd(total_devices_monthly_log_norm),
         mean_internal = mean(internal_devices_monthly_log_norm)) %>%
  arrange(cluster, month) %>%
  slice(1) %>%
  select(month, month_num, cluster, mean_total, sd_total, mean_internal)

# median flow rate by month
flow_stat_month <- visit_with_cluster %>%#cdphe_flow %>%
  #left_join(visit_with_cluster) %>%
  drop_na(flow_rate) %>%
  drop_na(cluster) %>%
  group_by(cluster, month) %>%
  mutate(mean_flow = median(flow_rate_log_norm),
         sd_flow = sd(flow_rate_log_norm)) %>%
arrange(cluster, month) %>%
  slice(1) %>%
  select(cluster, month, mean_flow, sd_flow)

visit_flow_stat_month <- left_join(visit_monthly_stat_month, flow_stat_month)

plot_month_high <- ggplot(data = visit_flow_stat_month %>% filter(cluster == "High-Variation")) +
  geom_pointrange(aes(x = month, y = mean_total,
                      ymin = mean_total - sd_total,
                      ymax = mean_total + sd_total), color = "orange3") +
  geom_line(aes(x = month_num, y = mean_total),
            color = "orange3", linetype = "dotted") +
  geom_pointrange(aes(x = month_num + 0.1, y = mean_flow,
                      ymin = mean_flow - sd_flow,
                      ymax = mean_flow + sd_flow), color = "grey10") +
  geom_line(aes(x = month_num + 0.1, y = mean_flow),
            color = "grey10", linetype = "dotted") +
  labs(x = NULL, y = "log(Monthly Device Count), Year-Normalized") +
  ggtitle("High-Variation Sewersheds (n = 13)") +
  scale_y_continuous(sec.axis = sec_axis(~.*1, name = "log(Monthly Flow Rate), Year-Normalized", labels = scales::comma)) +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 90, hjust = 0.95, vjust = 0.4),
        axis.title.y.left = element_text(margin = margin(l = 10, r = 5), color = "orange3"),
        axis.title.y.right = element_text(color = "grey10"),
        axis.text.y.left = element_text(color = "orange3"),
        axis.text.y.right = element_text(color = "grey10"),
        plot.margin = margin(15, 20, 10, 0, "pt"))

plot_month_low <- ggplot(data = visit_flow_stat_month %>% filter(cluster == "Low-Variation")) +
  geom_pointrange(aes(x = month, y = mean_total,
                      ymin = mean_total - sd_total,
                      ymax = mean_total + sd_total), color = "dodgerblue3") +
  geom_line(aes(x = month_num, y = mean_total),
            color = "dodgerblue3", linetype = "dotted") +
  geom_pointrange(aes(x = month_num + 0.1, y = mean_flow,
                      ymin = mean_flow - sd_flow,
                      ymax = mean_flow + sd_flow), color = "grey10") +
  geom_line(aes(x = month_num + 0.1, y = mean_flow),
            color = "grey10", linetype = "dotted") +
  labs(x = NULL, y = "log(Monthly Device Count), Year-Normalized") +
  ggtitle("Low-Variation Sewersheds (n = 53)") +
  scale_y_continuous(sec.axis = sec_axis(~.*1, name = "log(Monthly Flow Rate), Year-Normalized", labels = scales::comma)) +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 90, hjust = 0.95, vjust = 0.4),
        axis.title.y.left = element_text(margin = margin(l = 10, r = 5), color = "dodgerblue3"),
        axis.title.y.right = element_text(color = "grey10"),
        axis.text.y.left = element_text(color = "dodgerblue3"),
        axis.text.y.right = element_text(color = "grey10"),
        plot.margin = margin(15, 20, 10, 0, "pt"))

plot_month <- grid.arrange(plot_month_high, plot_month_low, ncol = 1)

ggsave("Figures/plot_month_mobile_flow.png", height = 9, width = 9, plot = plot_month)

violin_season_high <- ggplot(data = visit_with_cluster %>% filter(cluster == "High-Variation"),
                           aes(x = season, y = total_devices_daily_log_norm)) +
  #geom_boxplot(fill = "orange3", coef = Inf, staplewidth = 0.5) +
  geom_violin(fill = "orange3") +
  stat_summary(aes(x = season, y = total_devices_daily_log_norm),
               fun.min = function(z) {quantile(z, 0.25)},
               fun.max = function(z) {quantile(z, 0.75)},
               fun = median) +
  #stat_summary(fun.data = mean_sdl, geom = "pointrange") +
  #geom_hline(yintercept = 0, linetype = "dashed") +
  labs(x = NULL, y = "log(Daily Device Count), Year-Normalized") +
  ggtitle("High-Variation Sewersheds (n = 13)") +
  scale_y_continuous(limits = c(-10, 10)) +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 90, hjust = 0.95, vjust = 0.4),
        axis.title.y = element_text(margin = margin(l = 10, r = 5)),
        plot.margin = margin(15, 20, 10, 0, "pt"))

violin_season_low <- ggplot(data = visit_with_cluster %>% filter(cluster == "Low-Variation"),
                          aes(x = season, y = total_devices_daily_log_norm)) +
  #geom_boxplot(fill = "darkorchid3", coef = Inf, staplewidth = 0.5) +
  geom_violin(fill = "darkorchid3") +
  stat_summary(aes(x = season, y = total_devices_daily_log_norm),
               fun.min = function(z) {quantile(z, 0.25)},
               fun.max = function(z) {quantile(z, 0.75)},
               fun = median) +
  #stat_summary(fun.data = mean_sdl, geom = "pointrange") +
  #geom_hline(yintercept = 0, linetype = "dashed") +
  labs(x = NULL, y = "log(Daily Device Count), Year-Normalized") +
  ggtitle("Low-Variation Sewersheds (n = 53)") +
  scale_y_continuous(limits = c(-10, 10)) +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 90, hjust = 0.95, vjust = 0.4),
        axis.title.y = element_text(margin = margin(l = 10, r = 5)),
        plot.margin = margin(15, 20, 10, 0, "pt"))

#violin_season <- grid.arrange(violin_season_high, violin_season_low, ncol = 2)

summary(aov(visit_with_cluster_high$total_devices_daily_log_norm ~ visit_with_cluster_high$season))
summary(aov(visit_with_cluster_low$total_devices_daily_log_norm ~ visit_with_cluster_low$season))

violin_season_high_flow <- ggplot(data = visit_with_cluster %>% filter(cluster == "High-Variation"),
                             aes(x = season, y = flow_rate_log_norm)) +
  #geom_boxplot(fill = "orange3", coef = Inf, staplewidth = 0.5) +
  geom_violin(fill = "orange3") +
  stat_summary(aes(x = season, y = total_devices_daily_log_norm),
               fun.min = function(z) {quantile(z, 0.25)},
               fun.max = function(z) {quantile(z, 0.75)},
               fun = median) +
  #stat_summary(fun.data = mean_sdl, geom = "pointrange") +
  #geom_hline(yintercept = 0, linetype = "dashed") +
  labs(x = NULL, y = "log(Flow Rate), Year-Normalized") +
  ggtitle("High-Variation Sewersheds (n = 13)") +
  scale_y_continuous(limits = c(-10, 10)) +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 90, hjust = 0.95, vjust = 0.4),
        axis.title.y = element_text(margin = margin(l = 10, r = 5)),
        plot.margin = margin(15, 20, 10, 0, "pt"))

violin_season_low_flow <- ggplot(data = visit_with_cluster %>% filter(cluster == "Low-Variation"),
                            aes(x = season, y = flow_rate_log_norm)) +
  #geom_boxplot(fill = "darkorchid3", coef = Inf, staplewidth = 0.5) +
  geom_violin(fill = "darkorchid3") +
  stat_summary(aes(x = season, y = total_devices_daily_log_norm),
               fun.min = function(z) {quantile(z, 0.25)},
               fun.max = function(z) {quantile(z, 0.75)},
               fun = median) +
  #stat_summary(fun.data = mean_sdl, geom = "pointrange") +
  #geom_hline(yintercept = 0, linetype = "dashed") +
  labs(x = NULL, y = "log(Flow Rate), Year-Normalized") +
  ggtitle("Low-Variation Sewersheds (n = 53)") +
  scale_y_continuous(limits = c(-10, 10)) +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 90, hjust = 0.95, vjust = 0.4),
        axis.title.y = element_text(margin = margin(l = 10, r = 5)),
        plot.margin = margin(15, 20, 10, 0, "pt"))

violin_season <- grid.arrange(violin_season_high, violin_season_low,
  violin_season_high_flow, violin_season_low_flow, ncol = 2)

ggsave("Figures/violin_season.png", height = 9, width = 8, plot = violin_season)

summary(aov(visit_with_cluster_high$flow_rate_log_norm ~ visit_with_cluster_high$season))
summary(aov(visit_with_cluster_low$flow_rate_log_norm ~ visit_with_cluster_low$season))

# median flow rate by month
# flow_stat_month <- visit_with_cluster %>%#cdphe_flow %>%
#   #left_join(visit_with_cluster) %>%
#   drop_na(flow_rate) %>%
#   drop_na(cluster) %>%
#   group_by(cluster, month) %>%
#   mutate(median_flow_monthly_log_scale = median(flow_rate_log_norm),
#             lower_flow_monthly_log_scale = quantile(flow_rate_log_norm, 0.25),
#             upper_flow_monthly_log_scale = quantile(flow_rate_log_norm, 0.75)) %>%
#   slice(1)
# 
# plot_month_flow_high <- ggplot(data = flow_stat_month %>% filter(cluster == "High-Variation")) +
#   geom_pointrange(aes(x = month, y = median_flow_monthly_log_scale,
#                       ymin = lower_flow_monthly_log_scale,
#                       ymax = upper_flow_monthly_log_scale), color = "orange3") +
#   geom_line(aes(x = month_num, y = median_flow_monthly_log_scale),
#             color = "orange3", linetype = "dotted") +
#   labs(x = NULL, y = "log(Flow Rate), Year-Normalized") +
#   ggtitle("High-Variation Sewersheds (n = 13)") +
#   scale_y_continuous(labels = scales::comma) +
#   theme(plot.title = element_text(hjust = 0.5),
#         axis.text.x = element_text(angle = 90, hjust = 0.95, vjust = 0.4),
#         axis.title.y = element_text(margin = margin(l = 10, r = 5)),
#         plot.margin = margin(15, 20, 10, 0, "pt"))
# 
# plot_month_flow_low <- ggplot(data = flow_stat_month %>% filter(cluster == "Low-Variation")) +
#   geom_pointrange(aes(x = month, y = median_flow_monthly_log_scale,
#                       ymin = lower_flow_monthly_log_scale,
#                       ymax = upper_flow_monthly_log_scale), color = "darkorchid4") +
#   geom_line(aes(x = month_num, y = median_flow_monthly_log_scale),
#             color = "darkorchid4", linetype = "dotted") +
#   labs(x = NULL, y = "log(Flow Rate), Year-Normalized") +
#   ggtitle("Low-Variation Sewersheds (n = 53)") +
#   scale_y_continuous(labels = scales::comma) +
#   theme(plot.title = element_text(hjust = 0.5),
#         axis.text.x = element_text(angle = 90, hjust = 0.95, vjust = 0.4),
#         axis.title.y = element_text(margin = margin(l = 10, r = 5)),
#         plot.margin = margin(15, 20, 10, 0, "pt"))
# 
# plot_month_flow <- grid.arrange(plot_month_flow_high, plot_month_flow_low, ncol = 1)
# 
# ggsave("Figures/plot_month_flow.png", height = 9, width = 9, plot = plot_month_flow)

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

# calc_year <- visit_with_cluster %>%
#   drop_na(total_devices_monthly) %>%
#   group_by(cluster, year) %>%
#   summarize(mean = round(mean(total_devices_monthly), 0),
#             sd = round(sd(total_devices_monthly), 0))
# write.csv(calc_year, "DataProcessed/calc_year.csv")

# linear mixed effect models and LRT on log-transformed device counts
mod_dayofweek_high <- lmer(total_devices_daily_log_norm ~ 1 + (1 | sewershed_analysis_name), data = visit_with_cluster %>% filter(cluster == "High-Variation"), REML = F)
mod1_dayofweek_high <- lmer(total_devices_daily_log_norm ~ 1 + day_of_week + (1 | sewershed_analysis_name), data = visit_with_cluster %>% filter(cluster == "High-Variation"), REML = F)
anova(mod_dayofweek_high, mod1_dayofweek_high, test = "LRT")

mod_dayofweek_low <- lmer(total_devices_daily_log_norm ~ 1 + (1 | sewershed_analysis_name), data = visit_with_cluster %>% filter(cluster == "Low-Variation"), REML = F)
mod1_dayofweek_low <- lmer(total_devices_daily_log_norm ~ 1 + day_of_week + (1 | sewershed_analysis_name), data = visit_with_cluster %>% filter(cluster == "Low-Variation"), REML = F)
anova(mod_dayofweek_low, mod1_dayofweek_low, test = "LRT")

mod_wdwe_high <- lmer(total_devices_daily_log_norm ~ 1 + (1 | sewershed_analysis_name), data = visit_with_cluster %>% filter(cluster == "High-Variation"), REML = F)
mod1_wdwe_high <- lmer(total_devices_daily_log_norm ~ 1 + weekday_weekend + (1 | sewershed_analysis_name), data = visit_with_cluster %>% filter(cluster == "High-Variation"), REML = F)
anova(mod_wdwe_high, mod1_wdwe_high, test = "LRT")

mod_wdwe_low <- lmer(total_devices_daily_log_norm ~ 1 + (1 | sewershed_analysis_name), data = visit_with_cluster %>% filter(cluster == "Low-Variation"), REML = F)
mod1_wdwe_low <- lmer(total_devices_daily_log_norm ~ 1 + weekday_weekend + (1 | sewershed_analysis_name), data = visit_with_cluster %>% filter(cluster == "Low-Variation"), REML = F)
anova(mod_wdwe_low, mod1_wdwe_low, test = "LRT")

mod_month_high <- lmer(total_devices_monthly_log_norm ~ 1 + (1 | sewershed_analysis_name), data = visit_with_cluster %>% filter(cluster == "High-Variation"), REML = F)
mod1_month_high <- lmer(total_devices_monthly_log_norm ~ 1 + month + (1 | sewershed_analysis_name), data = visit_with_cluster %>% filter(cluster == "High-Variation"), REML = F)
anova(mod_month_high, mod1_month_high, test = "LRT")

mod_month_low <- lmer(total_devices_monthly_log_norm ~ 1 + (1 | sewershed_analysis_name), data = visit_with_cluster %>% filter(cluster == "Low-Variation"), REML = F)
mod1_month_low <- lmer(total_devices_monthly_log_norm ~ 1 + month + (1 | sewershed_analysis_name), data = visit_with_cluster %>% filter(cluster == "Low-Variation"), REML = F)
anova(mod_month_low, mod1_month_low, test = "LRT")

mod_season_high <- lmer(total_devices_monthly_log_norm ~ 1 + (1 | sewershed_analysis_name), data = visit_with_cluster %>% filter(cluster == "High-Variation"), REML = F)
mod1_season_high <- lmer(total_devices_monthly_log_norm ~ 1 + season + (1 | sewershed_analysis_name), data = visit_with_cluster %>% filter(cluster == "High-Variation"), REML = F)
anova(mod_season_high, mod1_season_high, test = "LRT")

mod_season_low <- lmer(total_devices_monthly_log_norm ~ 1 + (1 | sewershed_analysis_name), data = visit_with_cluster %>% filter(cluster == "Low-Variation"), REML = F)
mod1_season_low <- lmer(total_devices_monthly_log_norm ~ 1 + season + (1 | sewershed_analysis_name), data = visit_with_cluster %>% filter(cluster == "Low-Variation"), REML = F)
anova(mod_season_low, mod1_season_low, test = "LRT")

# mod_year_high <- lmer(total_devices_monthly_log ~ 1 + (1 | sewershed_analysis_name), data = log_df_monthly %>% filter(cluster == "High-Variation"), REML = F)
# mod1_year_high <- lmer(total_devices_monthly_log ~ 1 + year + (1 | sewershed_analysis_name), data = log_df_monthly %>% filter(cluster == "High-Variation"), REML = F)
# anova(mod_year_high, mod1_year_high, test = "LRT")
# 
# mod_year_low <- lmer(total_devices_monthly_log ~ 1 + (1 | sewershed_analysis_name), data = log_df_monthly %>% filter(cluster == "Low-Variation"), REML = F)
# mod1_year_low <- lmer(total_devices_monthly_log ~ 1 + year + (1 | sewershed_analysis_name), data = log_df_monthly %>% filter(cluster == "Low-Variation"), REML = F)
# anova(mod_year_low, mod1_year_low, test = "LRT")

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

# calc_year_internal <- visit_with_cluster %>%
#   drop_na(total_devices_monthly) %>%
#   group_by(cluster, year) %>%
#   summarize(mean = round(mean(internal_devices_monthly), 0),
#             sd = round(sd(internal_devices_monthly), 0))
# write.csv(calc_year_internal, "DataProcessed/calc_year_internal.csv")

mod_month_high_internal <- lmer(internal_devices_monthly_log_norm ~ 1 + (1 | sewershed_analysis_name), data = visit_with_cluster %>% filter(cluster == "High-Variation"), REML = F)
mod1_month_high_internal <- lmer(internal_devices_monthly_log_norm ~ 1 + month + (1 | sewershed_analysis_name), data = visit_with_cluster %>% filter(cluster == "High-Variation"), REML = F)
anova(mod_month_high_internal, mod1_month_high_internal, test = "LRT")

mod_month_low_internal <- lmer(internal_devices_monthly_log_norm ~ 1 + (1 | sewershed_analysis_name), data = visit_with_cluster %>% filter(cluster == "Low-Variation"), REML = F)
mod1_month_low_internal <- lmer(internal_devices_monthly_log_norm ~ 1 + month + (1 | sewershed_analysis_name), data = visit_with_cluster %>% filter(cluster == "Low-Variation"), REML = F)
anova(mod_month_low_internal, mod1_month_low_internal, test = "LRT")

mod_season_high_internal <- lmer(internal_devices_monthly_log_norm ~ 1 + (1 | sewershed_analysis_name), data = visit_with_cluster %>% filter(cluster == "High-Variation"), REML = F)
mod1_season_high_internal <- lmer(internal_devices_monthly_log_norm ~ 1 + season + (1 | sewershed_analysis_name), data = visit_with_cluster %>% filter(cluster == "High-Variation"), REML = F)
anova(mod_season_high_internal, mod1_season_high_internal, test = "LRT")

mod_season_low_internal <- lmer(internal_devices_monthly_log_norm ~ 1 + (1 | sewershed_analysis_name), data = visit_with_cluster %>% filter(cluster == "Low-Variation"), REML = F)
mod1_season_low_internal <- lmer(internal_devices_monthly_log_norm ~ 1 + season + (1 | sewershed_analysis_name), data = visit_with_cluster %>% filter(cluster == "Low-Variation"), REML = F)
anova(mod_season_low_internal, mod1_season_low_internal, test = "LRT")

# mod_year_high_internal <- lmer(internal_devices_monthly_log ~ 1 + (1 | sewershed_analysis_name), data = log_df_monthly %>% filter(cluster == "High-Variation"), REML = F)
# mod1_year_high_internal <- lmer(internal_devices_monthly_log ~ 1 + year + (1 | sewershed_analysis_name), data = log_df_monthly %>% filter(cluster == "High-Variation"), REML = F)
# anova(mod_year_high_internal, mod1_year_high_internal, test = "LRT")
# 
# mod_year_low_internal <- lmer(internal_devices_monthly_log ~ 1 + (1 | sewershed_analysis_name), data = log_df_monthly %>% filter(cluster == "Low-Variation"), REML = F)
# mod1_year_low_internal <- lmer(internal_devices_monthly_log ~ 1 + year + (1 | sewershed_analysis_name), data = log_df_monthly %>% filter(cluster == "Low-Variation"), REML = F)
# anova(mod_year_low_internal, mod1_year_low_internal, test = "LRT")

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

# calc_year_flow <- cdphe_flow %>%
#   group_by(cluster, year) %>%
#   summarize(mean = signif(mean(flow_rate), 3),
#             sd = signif(sd(flow_rate), 3))
# write.csv(calc_year_flow, "DataProcessed/calc_year_flow.csv")

mod_month_flow_high <-  lmer(flow_rate_log_norm ~ 1 +               (1 | sewershed_analysis_name), data = visit_with_cluster %>% filter(cluster == "High-Variation"), REML = F)
mod1_month_flow_high <- lmer(flow_rate_log_norm ~ 1 + month + (1 | sewershed_analysis_name), data = visit_with_cluster %>% filter(cluster == "High-Variation"), REML = F)
anova(mod_month_flow_high, mod1_month_flow_high, test = "LRT")

mod_month_flow_low <- lmer(flow_rate_log_norm ~ 1 + (1 | sewershed_analysis_name), data = visit_with_cluster %>% filter(cluster == "Low-Variation"), REML = F)
mod1_month_flow_low <- lmer(flow_rate_log_norm ~ 1 + month + (1 | sewershed_analysis_name), data = visit_with_cluster %>% filter(cluster == "Low-Variation"), REML = F)
anova(mod_month_flow_low, mod1_month_flow_low, test = "LRT")

mod_season_flow_high <-  lmer(flow_rate_log_norm ~ 1 +               (1 | sewershed_analysis_name), data = visit_with_cluster %>% filter(cluster == "High-Variation"), REML = F)
mod1_season_flow_high <- lmer(flow_rate_log_norm ~ 1 + season + (1 | sewershed_analysis_name), data = visit_with_cluster %>% filter(cluster == "High-Variation"), REML = F)
anova(mod_season_flow_high, mod1_season_flow_high, test = "LRT")

mod_season_flow_low <- lmer(flow_rate_log_norm ~ 1 + (1 | sewershed_analysis_name), data = visit_with_cluster %>% filter(cluster == "Low-Variation"), REML = F)
mod1_season_flow_low <- lmer(flow_rate_log_norm ~ 1 + season + (1 | sewershed_analysis_name), data = visit_with_cluster %>% filter(cluster == "Low-Variation"), REML = F)
anova(mod_season_flow_low, mod1_season_flow_low, test = "LRT")

########################################################################################################


# mod_year_flow_high <-  lmer(flow_rate_log ~ 1 +               (1 | sewershed_analysis_name), data = cdphe_flow_log %>% filter(cluster == "High-Variation"), REML = F)
# mod1_year_flow_high <- lmer(flow_rate_log ~ 1 + year + (1 | sewershed_analysis_name), data = cdphe_flow_log %>% filter(cluster == "High-Variation"), REML = F)
# anova(mod_year_flow_high, mod1_year_flow_high, test = "LRT")
# 
# mod_year_flow_low <- lmer(flow_rate_log ~ 1 + (1 | sewershed_analysis_name), data = cdphe_flow_log %>% filter(cluster == "Low-Variation"), REML = F)
# mod1_year_flow_low <- lmer(flow_rate_log ~ 1 + year + (1 | sewershed_analysis_name), data = cdphe_flow_log %>% filter(cluster == "Low-Variation"), REML = F)
# anova(mod_year_flow_low, mod1_year_flow_low, test = "LRT")


# median device count by year
# visit_monthly_stat_year <- visit_with_cluster %>%
#   drop_na(total_devices_monthly) %>%
#   group_by(cluster, year) %>%
#   summarize(median_total_devices_year = median(total_devices_monthly),
#             lower_total_devices_year = quantile(total_devices_monthly, 0.25),
#             upper_total_devices_year = quantile(total_devices_monthly, 0.75),
#             median_internal_devices_year = median(internal_devices_monthly))
# 
# plot_year_high <- ggplot(data = visit_monthly_stat_year %>% filter(cluster == "High-Variation")) +
#   geom_pointrange(aes(x = year, y = median_total_devices_year,
#                       ymin = lower_total_devices_year,
#                       ymax = upper_total_devices_year)) +
#   labs(x = NULL, y = "Monthly Device Count (Median and IQR)") +
#   ggtitle("High-Variation Sewersheds (n = 15)") +
#   scale_y_continuous(labels = scales::comma) +
#   theme(plot.title = element_text(hjust = 0.5),
#         axis.text.x = element_text(angle = 90, hjust = 0.95, vjust = 0.4),
#         axis.title.y = element_text(margin = margin(l = 10, r = 5)),
#         plot.margin = margin(15, 20, 10, 0, "pt"))
# 
# plot_year_low <- ggplot(data = visit_monthly_stat_year %>% filter(cluster == "Low-Variation")) +
#   geom_pointrange(aes(x = year, y = median_total_devices_year,
#                       ymin = lower_total_devices_year,
#                       ymax = upper_total_devices_year)) +
#   labs(x = NULL, y = "Monthly Device Count (Median and IQR)") +
#   ggtitle("Low-Variation Sewersheds (n = 51)") +
#   scale_y_continuous(labels = scales::comma) +
#   theme(plot.title = element_text(hjust = 0.5),
#         axis.text.x = element_text(angle = 90, hjust = 0.95, vjust = 0.4),
#         axis.title.y = element_text(margin = margin(l = 10, r = 5)),
#         plot.margin = margin(15, 20, 10, 0, "pt"))
# 
# # median flow rate by year
# flow_stat_year <- cdphe_flow %>%
#   left_join(visit_with_cluster) %>%
#   drop_na(flow_rate) %>%
#   drop_na(cluster) %>%
#   group_by(cluster, year) %>%
#   summarize(median_flow_yearly = median(flow_rate),
#             lower_flow_yearly = quantile(flow_rate, 0.25),
#             upper_flow_yearly = quantile(flow_rate, 0.75))
# 
# plot_year_flow_high <- ggplot(data = flow_stat_year %>% filter(cluster == "High-Variation")) +
#   geom_pointrange(aes(x = year, y = median_flow_yearly,
#                       ymin = lower_flow_yearly,
#                       ymax = upper_flow_yearly)) +
#   labs(x = NULL, y = "Yearly Flow Rate (Median and IQR)") +
#   ggtitle("High-Variation Sewersheds (n = 15)") +
#   scale_y_continuous(labels = scales::comma) +
#   theme(plot.title = element_text(hjust = 0.5),
#         axis.text.x = element_text(angle = 90, hjust = 0.95, vjust = 0.4),
#         axis.title.y = element_text(margin = margin(l = 10, r = 5)),
#         plot.margin = margin(15, 20, 10, 0, "pt"))
# 
# plot_year_flow_low <- ggplot(data = flow_stat_year %>% filter(cluster == "Low-Variation")) +
#   geom_pointrange(aes(x = year, y = median_flow_yearly,
#                       ymin = lower_flow_yearly,
#                       ymax = upper_flow_yearly)) +
#   labs(x = NULL, y = "Yearly Flow Rate (Median and IQR)") +
#   ggtitle("Low-Variation Sewersheds (n = 51)") +
#   scale_y_continuous(labels = scales::comma) +
#   theme(plot.title = element_text(hjust = 0.5),
#         axis.text.x = element_text(angle = 90, hjust = 0.95, vjust = 0.4),
#         axis.title.y = element_text(margin = margin(l = 10, r = 5)),
#         plot.margin = margin(15, 20, 10, 0, "pt"))
# 
# plot_year_mobile_flow <- grid.arrange(plot_year_high, plot_year_flow_high,
#                                       plot_year_low, plot_year_flow_low, ncol = 2)
# 
# ggsave("Figures/plot_year_mobile_flow.png", height = 7, width = 14, plot = plot_year_mobile_flow)
# 
# # log-transformed daily device counts
# log_df_daily <- visit_with_cluster %>%
#   mutate(total_devices_daily = ifelse(total_devices_daily == 0, 1, total_devices_daily),
#          total_devices_daily_log = log(total_devices_daily))
# 
# hist_raw_daily_high <- ggplot(data = log_df_daily %>% filter(cluster == "High-Variation")) +
#   geom_histogram(aes(x = total_devices_daily)) +
#   labs(x = "Daily Device Count", y = "Frequency (Sewershed-Days)") +
#   ggtitle("High-Variation Sewersheds") +
#   scale_x_continuous(labels = scales::comma) +
#   scale_y_continuous(labels = scales::comma) +
#   theme(plot.title = element_text(hjust = 0.5))
# 
# hist_log_daily_high <- ggplot(data = log_df_daily %>% filter(cluster == "High-Variation")) +
#   geom_histogram(aes(x = total_devices_daily_log)) +
#   labs(x = "log(Daily Device Count)", y = "Frequency (Sewershed-Days)") +
#   ggtitle("High-Variation Sewersheds") +
#   scale_x_continuous(labels = scales::comma) +
#   scale_y_continuous(labels = scales::comma) +
#   theme(plot.title = element_text(hjust = 0.5))
# 
# qq_daily_high <- ggplot(data = log_df_daily %>% filter(cluster == "High-Variation"),
#                         aes(sample = total_devices_daily_log)) +
#   stat_qq() +
#   stat_qq_line() +
#   labs(x = "Theoretical Quantiles", y = "Sample Quantiles") +
#   ggtitle("High-Variation Sewersheds") +
#   theme(plot.title = element_text(hjust = 0.5))
# 
# hist_raw_daily_low <- ggplot(data = log_df_daily %>% filter(cluster == "Low-Variation")) +
#   geom_histogram(aes(x = total_devices_daily)) +
#   labs(x = "Daily Device Count", y = "Frequency (Sewershed-Days)") +
#   ggtitle("Low-Variation Sewersheds") +
#   scale_x_continuous(labels = scales::comma) +
#   scale_y_continuous(labels = scales::comma) +
#   theme(plot.title = element_text(hjust = 0.5))
# 
# hist_log_daily_low <- ggplot(data = log_df_daily %>% filter(cluster == "Low-Variation")) +
#   geom_histogram(aes(x = total_devices_daily_log)) +
#   labs(x = "log(Daily Device Count)", y = "Frequency (Sewershed-Days)") +
#   ggtitle("Low-Variation Sewersheds") +
#   scale_x_continuous(labels = scales::comma) +
#   scale_y_continuous(labels = scales::comma) +
#   theme(plot.title = element_text(hjust = 0.5))
# 
# qq_daily_low <- ggplot(data = log_df_daily %>% filter(cluster == "Low-Variation"),
#                         aes(sample = total_devices_daily_log)) +
#   stat_qq() +
#   stat_qq_line() +
#   labs(x = "Theoretical Quantiles", y = "Sample Quantiles") +
#   ggtitle("Low-Variation Sewersheds") +
#   theme(plot.title = element_text(hjust = 0.5))
# 
# raw_log_qq_daily <- grid.arrange(hist_raw_daily_high, hist_log_daily_high, qq_daily_high,
#                                  hist_raw_daily_low, hist_log_daily_low, qq_daily_low, ncol = 3)
# 
# ggsave("Figures/raw_log_qq_daily.png", height = 8, width = 14, plot = raw_log_qq_daily)
# 
# # log-transformed monthly device counts
# log_df_monthly <- visit_with_cluster %>%
#   drop_na(total_devices_monthly) %>%
#   mutate(total_devices_monthly_log = log(total_devices_monthly),
#          internal_devices_monthly_log = log(internal_devices_monthly))
# 
# hist_raw_monthly_high <- ggplot(data = log_df_monthly %>% filter(cluster == "High-Variation")) +
#   geom_histogram(aes(x = total_devices_monthly)) +
#   labs(x = "Monthly Device Count", y = "Frequency (Sewershed-Months)") +
#   ggtitle("High-Variation Sewersheds") +
#   scale_x_continuous(labels = scales::comma) +
#   scale_y_continuous(labels = scales::comma) +
#   theme(plot.title = element_text(hjust = 0.5))
# 
# hist_log_monthly_high <- ggplot(data = log_df_monthly %>% filter(cluster == "High-Variation")) +
#   geom_histogram(aes(x = total_devices_monthly_log)) +
#   labs(x = "log(Monthly Device Count)", y = "Frequency (Sewershed-Months)") +
#   ggtitle("High-Variation Sewersheds") +
#   scale_x_continuous(labels = scales::comma) +
#   scale_y_continuous(labels = scales::comma) +
#   theme(plot.title = element_text(hjust = 0.5))
# 
# qq_monthly_high <- ggplot(data = log_df_monthly %>% filter(cluster == "High-Variation"),
#                           aes(sample = total_devices_monthly_log)) +
#   stat_qq() +
#   stat_qq_line() +
#   labs(x = "Theoretical Quantiles", y = "Sample Quantiles") +
#   ggtitle("High-Variation Sewersheds") +
#   theme(plot.title = element_text(hjust = 0.5))
# 
# hist_raw_monthly_low <- ggplot(data = log_df_monthly %>% filter(cluster == "Low-Variation")) +
#   geom_histogram(aes(x = total_devices_monthly)) +
#   labs(x = "Monthly Device Count", y = "Frequency (Sewershed-Months)") +
#   ggtitle("Low-Variation Sewersheds") +
#   scale_x_continuous(labels = scales::comma) +
#   scale_y_continuous(labels = scales::comma) +
#   theme(plot.title = element_text(hjust = 0.5))
# 
# hist_log_monthly_low <- ggplot(data = log_df_monthly %>% filter(cluster == "Low-Variation")) +
#   geom_histogram(aes(x = total_devices_monthly_log)) +
#   labs(x = "log(Monthly Device Count)", y = "Frequency (Sewershed-Months)") +
#   ggtitle("Low-Variation Sewersheds") +
#   scale_x_continuous(labels = scales::comma) +
#   scale_y_continuous(labels = scales::comma) +
#   theme(plot.title = element_text(hjust = 0.5))
# 
# qq_monthly_low <- ggplot(data = log_df_monthly %>% filter(cluster == "Low-Variation"),
#                           aes(sample = total_devices_monthly_log)) +
#   stat_qq() +
#   stat_qq_line() +
#   labs(x = "Theoretical Quantiles", y = "Sample Quantiles") +
#   ggtitle("Low-Variation Sewersheds") +
#   theme(plot.title = element_text(hjust = 0.5))
# 
# raw_log_qq_monthly <- grid.arrange(hist_raw_monthly_high, hist_log_monthly_high, qq_monthly_high,
#                                    hist_raw_monthly_low, hist_log_monthly_low, qq_monthly_low, ncol = 3)
# 
# ggsave("Figures/raw_log_qq_monthly.png", height = 8, width = 14, plot = raw_log_qq_monthly)
# 
# # log-transformed flow rates
# cdphe_flow_log <- cdphe_flow %>%
#   drop_na(flow_rate) %>%
#   mutate(flow_rate = ifelse(flow_rate == 0, 0.001, flow_rate),
#          flow_rate_log = log(flow_rate)) %>%
#   filter(day_of_week %in% c("Monday", "Tuesday", "Wednesday", "Thursday", "Friday"))
# 
# hist_raw_flow_high <- ggplot(data = cdphe_flow_log %>% filter(cluster == "High-Variation")) +
#   geom_histogram(aes(x = flow_rate)) +
#   labs(x = "Flow Rate", y = "Frequency (Sewershed-Days)") +
#   ggtitle("High-Variation Sewersheds") +
#   scale_x_continuous(labels = scales::comma) +
#   scale_y_continuous(labels = scales::comma) +
#   theme(plot.title = element_text(hjust = 0.5)); hist_raw_flow_high
# 
# hist_log_flow_high <- ggplot(data = cdphe_flow_log %>% filter(cluster == "High-Variation")) +
#   geom_histogram(aes(x = flow_rate_log)) +
#   labs(x = "log(Flow Rate)", y = "Frequency (Sewershed-Days)") +
#   ggtitle("High-Variation Sewersheds") +
#   scale_x_continuous(labels = scales::comma) +
#   scale_y_continuous(labels = scales::comma) +
#   theme(plot.title = element_text(hjust = 0.5)); hist_log_flow_high
# 
# qq_flow_high <- ggplot(data = cdphe_flow_log %>% filter(cluster == "High-Variation"),
#                        aes(sample = flow_rate_log)) +
#   stat_qq() +
#   stat_qq_line() +
#   labs(x = "Theoretical Quantiles", y = "Sample Quantiles") +
#   ggtitle("High-Variation Sewersheds") +
#   theme(plot.title = element_text(hjust = 0.5)); qq_flow_high
# 
# hist_raw_flow_low <- ggplot(data = cdphe_flow_log %>% filter(cluster == "Low-Variation")) +
#   geom_histogram(aes(x = flow_rate)) +
#   labs(x = "Flow Rate", y = "Frequency (Sewershed-Days)") +
#   ggtitle("Low-Variation Sewersheds") +
#   scale_x_continuous(labels = scales::comma) +
#   scale_y_continuous(labels = scales::comma) +
#   theme(plot.title = element_text(hjust = 0.5)); hist_raw_flow_low
# 
# hist_log_flow_low <- ggplot(data = cdphe_flow_log %>% filter(cluster == "Low-Variation")) +
#   geom_histogram(aes(x = flow_rate_log)) +
#   labs(x = "log(Flow Rate)", y = "Frequency (Sewershed-Days)") +
#   ggtitle("Low-Variation Sewersheds") +
#   scale_x_continuous(labels = scales::comma) +
#   scale_y_continuous(labels = scales::comma) +
#   theme(plot.title = element_text(hjust = 0.5)); hist_log_flow_low
# 
# qq_flow_low <- ggplot(data = cdphe_flow_log %>% filter(cluster == "Low-Variation"),
#                       aes(sample = flow_rate_log)) +
#   stat_qq() +
#   stat_qq_line() +
#   labs(x = "Theoretical Quantiles", y = "Sample Quantiles") +
#   ggtitle("Low-Variation Sewersheds") +
#   theme(plot.title = element_text(hjust = 0.5)); qq_flow_low
# 
# raw_log_qq_flow <- grid.arrange(hist_raw_flow_high, hist_log_flow_high, qq_flow_high,
#                                 hist_raw_flow_low, hist_log_flow_low, qq_flow_low, ncol = 3)
# 
# ggsave("Figures/raw_log_qq_flow.png", height = 8, width = 14, plot = raw_log_qq_flow)

# read in daily cluster data
# clusters_daily <- read.csv("DataProcessed/2025_07_10_cluster_results/daily/2025_06_25_daily_kmeans_cluster_labels.csv") %>%
#   mutate(cdphe_flow_name = AREA_SEWERSHED) %>%
#   select(cdphe_flow_name, kmeans_k2)
# 
# visit_with_cluster_daily <- left_join(visit, clusters_daily) %>%
#   filter(cdphe_flow_name %in% cdphe_flow$cdphe_flow_name,
#          date >= min(cdphe_flow$date),
#          date <= max(cdphe_flow$date)) %>%
#   mutate(cluster = ifelse(kmeans_k2 == 1, "High-Variation", "Low-Variation"))
# 
# visit_with_cluster_daily %>%
#   group_by(sewershed_analysis_name) %>%
#   slice(1) %>%
#   group_by(cluster) %>%
#   count()
# 
# spaghetti_high_daily <- visit_with_cluster_daily %>%
#   filter(cluster == "High-Variation") %>%
#   group_by(sewershed_analysis_name) %>%
#   drop_na() %>%
#   mutate(spaghetti_value = scale(log(total_devices_daily)))
# 
# spaghetti_low_daily <- visit_with_cluster_daily %>%
#   filter(cluster == "Low-Variation") %>%
#   group_by(sewershed_analysis_name) %>%
#   drop_na() %>%
#   mutate(spaghetti_value = scale(log(total_devices_daily)))
# 
# spaghetti_plot_high_daily <- ggplot(spaghetti_high_daily) +
#   geom_line(aes(x = date, y = spaghetti_value, group = sewershed_analysis_name, color = sewershed_analysis_name),
#             color = "orange2", alpha = 0.6) +
#   labs(x = NULL, y = "log(Daily Device Count), Z-Score Normalized") +
#   ggtitle("High-Variation Sewersheds (n = 14)") +
#   scale_x_date(date_breaks = "4 months", date_labels = "%b %Y") +
#   theme(plot.title = element_text(hjust = 0.5),
#         axis.text.x = element_text(angle = 90, hjust = 0.95, vjust = 0.4),
#         axis.title.y = element_text(margin = margin(l = 10, r = 5)),
#         plot.margin = margin(10, 20, 10, 0, "pt"))
# 
# spaghetti_plot_low_daily <- ggplot(spaghetti_low_daily) +
#   geom_line(aes(x = date, y = spaghetti_value, group = sewershed_analysis_name, color = sewershed_analysis_name),
#             color = "darkorchid4", alpha = 0.6) +
#   labs(x = NULL, y = "log(Daily Device Count), Z-Score Normalized") +
#   ggtitle("Low-Variation Sewersheds (n = 52)") +
#   scale_x_date(date_breaks = "4 months", date_labels = "%b %Y") +
#   theme(plot.title = element_text(hjust = 0.5),
#         axis.text.x = element_text(angle = 90, hjust = 0.95, vjust = 0.4),
#         axis.title.y = element_text(margin = margin(l = 10, r = 5)),
#         plot.margin = margin(10, 20, 10, 0, "pt"))
# 
# spaghetti_cluster_daily <- grid.arrange(spaghetti_plot_high_daily,
#                                         spaghetti_plot_low_daily, ncol = 1)
# 
# # save spaghetti plot
# ggsave("Figures/spaghetti_cluster_daily.png", height = 8, width = 10, plot = spaghetti_cluster_daily)

# get sample size by mobility cluster (flow)
# cdphe_flow %>%
#   group_by(sewershed_analysis_name) %>%
#   slice(1) %>%
#   group_by(cluster) %>%
#   count()

#%>%


# plot_dayofweek_high <- ggplot(data = visit_daily_stat_dayofweek %>% filter(cluster == "High-Variation")) +
#   geom_pointrange(aes(x = day_of_week, y = mean,
#                       ymin = mean - sd,
#                       ymax = mean + sd), color = "orange3") +
#   geom_line(aes(x = day_of_week_num, y = mean),
#             color = "orange3", linetype = "dotted") +
#   labs(x = NULL, y = "log(Daily Device Count), Year-Normalized") +
#   ggtitle("High-Variation Sewersheds (n = 13)") +
#   scale_y_continuous(limits = c(-1.5, 1.5)) +
#   theme(plot.title = element_text(hjust = 0.5),
#         axis.text.x = element_text(angle = 90, hjust = 0.95, vjust = 0.4),
#         axis.title.y = element_text(margin = margin(l = 10, r = 5)),
#         plot.margin = margin(15, 20, 10, 0, "pt"))
# 
# plot_dayofweek_low <- ggplot(data = visit_daily_stat_dayofweek %>% filter(cluster == "Low-Variation")) +
#   geom_pointrange(aes(x = day_of_week, y = mean,
#                       ymin = mean - sd,
#                       ymax = mean + sd), color = "dodgerblue3") +
#   geom_line(aes(x = day_of_week_num, y = mean),
#             color = "dodgerblue3", linetype = "dotted") +
#   labs(x = NULL, y = "Log(Daily Device Count), Year-Normalized") +
#   ggtitle("Low-Variation Sewersheds (n = 53)") +
#   scale_y_continuous(limits = c(-1.5, 1.5)) +
#   theme(plot.title = element_text(hjust = 0.5),
#         axis.text.x = element_text(angle = 90, hjust = 0.95, vjust = 0.4),
#         axis.title.y = element_text(margin = margin(l = 10, r = 5)),
#         plot.margin = margin(15, 20, 10, 0, "pt"))
# 
# plot_dayofweek <- grid.arrange(plot_dayofweek_high, plot_dayofweek_low, ncol = 1)
# 
# ggsave("Figures/plot_dayofweek.png", height = 9, width = 9, plot = plot_dayofweek)
# 
# violin_dayofweek_high <- ggplot(data = visit_with_cluster %>% filter(cluster == "High-Variation"),
#                            aes(x = day_of_week, y = total_devices_daily_log_norm)) +
#   geom_violin(fill = "orange3") +
#   stat_summary(aes(x = day_of_week, y = total_devices_daily_log_norm),
#                fun.min = function(z) {quantile(z, 0.25)},
#                fun.max = function(z) {quantile(z, 0.75)},
#                fun = median) +
#   labs(x = NULL, y = "log(Daily Device Count), Year-Normalized") +
#   ggtitle("High-Variation Sewersheds (n = 13)") +
#   scale_y_continuous(limits = c(-10, 10)) +
#   theme(plot.title = element_text(hjust = 0.5),
#         axis.text.x = element_text(angle = 90, hjust = 0.95, vjust = 0.4),
#         axis.title.y = element_text(margin = margin(l = 10, r = 5)),
#         plot.margin = margin(15, 20, 10, 0, "pt"))
# 
# violin_dayofweek_low <- ggplot(data = visit_with_cluster %>% filter(cluster == "Low-Variation"),
#                           aes(x = day_of_week, y = total_devices_daily_log_norm)) +
#   geom_violin(fill = "darkorchid3") +
#   stat_summary(aes(x = day_of_week, y = total_devices_daily_log_norm),
#                fun.min = function(z) {quantile(z, 0.25)},
#                fun.max = function(z) {quantile(z, 0.75)},
#                fun = median) +
#   labs(x = NULL, y = "log(Daily Device Count), Year-Normalized") +
#   ggtitle("Low-Variation Sewersheds (n = 53)") +
#   scale_y_continuous(limits = c(-10, 10)) +
#   theme(plot.title = element_text(hjust = 0.5),
#         axis.text.x = element_text(angle = 90, hjust = 0.95, vjust = 0.4),
#         axis.title.y = element_text(margin = margin(l = 10, r = 5)),
#         plot.margin = margin(15, 20, 10, 0, "pt")); violin_dayofweek_low
# 
# violin_dayofweek <- grid.arrange(violin_dayofweek_high, violin_dayofweek_low, ncol = 1)
# 
# violin_wdwe_high <- ggplot(data = visit_with_cluster %>% filter(cluster == "High-Variation"),
#                         aes(x = weekday_weekend, y = total_devices_daily_log_norm)) +
#   geom_violin(fill = "orange3") +
#   stat_summary(aes(x = weekday_weekend, y = total_devices_daily_log_norm),
#                fun.min = function(z) {quantile(z, 0.25)},
#                fun.max = function(z) {quantile(z, 0.75)},
#                fun = median) +
#   labs(x = NULL, y = "log(Daily Device Count), Year-Normalized") +
#   ggtitle("High-Variation Sewersheds (n = 13)") +
#   scale_y_continuous(limits = c(-10, 10)) +
#   theme(plot.title = element_text(hjust = 0.5),
#         axis.text.x = element_text(angle = 90, hjust = 0.95, vjust = 0.4),
#         axis.title.y = element_text(margin = margin(l = 10, r = 5)),
#         plot.margin = margin(15, 20, 10, 0, "pt"))
# 
# violin_wdwe_low <- ggplot(data = visit_with_cluster %>% filter(cluster == "Low-Variation"),
#                        aes(x = weekday_weekend, y = total_devices_daily_log_norm)) +
#   geom_violin(fill = "darkorchid3") +
#   stat_summary(aes(x = weekday_weekend, y = total_devices_daily_log_norm),
#                fun.min = function(z) {quantile(z, 0.25)},
#                fun.max = function(z) {quantile(z, 0.75)},
#                fun = median) +
#   labs(x = NULL, y = "log(Daily Device Count), Year-Normalized") +
#   ggtitle("Low-Variation Sewersheds (n = 53)") +
#   scale_y_continuous(limits = c(-10, 10)) +
#   theme(plot.title = element_text(hjust = 0.5),
#         axis.text.x = element_text(angle = 90, hjust = 0.95, vjust = 0.4),
#         axis.title.y = element_text(margin = margin(l = 10, r = 5)),
#         plot.margin = margin(15, 20, 10, 0, "pt"))
# 
# violin_wdwe <- grid.arrange(violin_wdwe_high, violin_wdwe_low, ncol = 2)
# 
# ggsave("Figures/violin_wdwe.png", height = 8, width = 12, plot = violin_wdwe)
# 
# visit_with_cluster_high <- visit_with_cluster %>%
#   filter(cluster == "High-Variation")
# 
# t.test(visit_with_cluster_high$total_devices_daily_log_norm ~ visit_with_cluster_high$weekday_weekend)
# 
# visit_with_cluster_low <- visit_with_cluster %>%
#   filter(cluster == "Low-Variation")
# 
# t.test(visit_with_cluster_low$total_devices_daily_log_norm ~ visit_with_cluster_low$weekday_weekend)