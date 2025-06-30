#############################################################################
##### CHARACTERIZING DYNAMIC POPULATIONS IN WASTEWATER MONITORING AREAS #####
############################## SUMMER 2025 ##################################
################################## EJW ######################################

# ANALYSIS 2. DESCRIPTIVE STATISTICS AND REGRESSION ON DEVICE COUNTS AND FLOW RATES

# set working directory
setwd("/Users/emwu9912/Documents/CU Anschutz/COVID Wastewater/")

# load libraries
library(pacman)
p_load(gridExtra, lme4, lubridate, tidyverse)

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
  # 0.2% areal overlap: if a sewershed occupies >= 0.2% of a CBG's land mass, assign that CBG to the sewershed
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
  left_join(sewershed_names)

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
  filter(cdphe_flow_name %in% reporting_start$cdphe_flow_name) #%>%

cdphe_flow <- cdphe_flow4 %>%
  group_by(sewershed_analysis_name) %>%
  mutate(ndays = as.numeric(difftime(max(date), min(date), "days"))) %>%
  filter(ndays >= 365) #%>%
  # exclude the two values reported by Aspen and Breckenridge in December 2021
  #mutate(date = case_when(cluster == "Low-Variation" ~ date,
                          #cluster == "High-Variation" & date > as.Date("2022-01-01") ~ date,
                          #TRUE ~ NA_Date_)) %>%
  #drop_na(date)

# how many sewersheds included in the analysis?
unique(unlist(cdphe_flow$cdphe_flow_name))

save(cdphe_flow, file = "DataProcessed/cdphe_flow.Rda")

# perform K-Means clustering on the 66 sewersheds that have both mobile device data, flow data,
# and reported over a period of 365 days or more

# read in K-Means cluster data
clusters <- read.csv("DataProcessed/new_cluster_results/monthly/2025_06_25_monthly_kmeans_cluster_labels.csv") %>%
  rename(cdphe_flow_name = AREA_SEWERSHED) %>%
  select(cdphe_flow_name, kmeans_k2)

# merge with clusters and define temporal groupings
visit_with_cluster <- left_join(visit, clusters) %>%
  filter(cdphe_flow_name %in% cdphe_flow$cdphe_flow_name,
         date >= min(cdphe_flow$date),
         date <= max(cdphe_flow$date)) %>%
  mutate(cluster = ifelse(kmeans_k2 == 1, "High-Variation", "Low-Variation"),
         month = as.factor(month(date, label = T, abbr = F)),
         season = as.factor(case_when(month %in% c("December", "January", "February") ~ "Winter",
                                      month %in% c("March", "April", "May") ~ "Spring",
                                      month %in% c("June", "July", "August") ~ "Summer",
                                      TRUE ~ "Fall")),
         year = as.factor(year),
         day_of_week = wday(date, label = T, abbr = F),
         weekday_weekend = as.factor(case_when(day_of_week %in% c("Monday", "Tuesday", "Wednesday", "Thursday", "Friday") ~ "Weekday",
                                               TRUE ~ "Weekend")))

# get sample size by mobility cluster
visit_with_cluster %>%
  group_by(sewershed_analysis_name) %>%
  slice(1) %>%
  group_by(cluster) %>%
  count()

# create cluster table for appendix
cluster_table <- visit_with_cluster %>%
  group_by(sewershed_result_name) %>%
  slice(1) %>%
  arrange(cluster, sewershed_result_name) %>%
  select(sewershed_result_name, cluster)

write.csv(cluster_table, "DataProcessed/cluster_table.csv")

spaghetti_high <- visit_with_cluster %>%
  filter(cluster == "High-Variation") %>%
  group_by(sewershed_analysis_name) %>%
  drop_na() %>%
  mutate(spaghetti_value = scale(log(total_devices_monthly)))

spaghetti_low <- visit_with_cluster %>%
  filter(cluster == "Low-Variation") %>%
  group_by(sewershed_analysis_name) %>%
  drop_na() %>%
  mutate(spaghetti_value = scale(log(total_devices_monthly)))

spaghetti_plot_high <- ggplot(spaghetti_high) +
  geom_line(aes(x = date, y = spaghetti_value, group = sewershed_analysis_name, color = sewershed_analysis_name),
            color = "orange2", alpha = 0.6) +
  labs(x = NULL, y = "log(Monthly Device Count), Z-Score Normalized") +
  ggtitle("High-Variation Sewersheds (n = 15)") +
  scale_x_date(date_breaks = "4 months", date_labels = "%b %Y") +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 90, hjust = 0.95, vjust = 0.4),
        axis.title.y = element_text(margin = margin(l = 10, r = 5)),
        plot.margin = margin(10, 20, 10, 0, "pt"))

spaghetti_plot_low <- ggplot(spaghetti_low) +
  geom_line(aes(x = date, y = spaghetti_value, group = sewershed_analysis_name, color = sewershed_analysis_name),
            color = "darkorchid4", alpha = 0.6) +
  labs(x = NULL, y = "log(Monthly Device Count), Z-Score Normalized") +
  ggtitle("Low-Variation Sewersheds (n = 51)") +
  scale_x_date(date_breaks = "4 months", date_labels = "%b %Y") +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 90, hjust = 0.95, vjust = 0.4),
        axis.title.y = element_text(margin = margin(l = 10, r = 5)),
        plot.margin = margin(10, 20, 10, 0, "pt"))

spaghetti_cluster <- grid.arrange(spaghetti_plot_high, spaghetti_plot_low, ncol = 1)

# save spaghetti plot
ggsave("Figures/spaghetti_cluster.png", height = 8, width = 10, plot = spaghetti_cluster)

# read in daily cluster data
clusters_daily <- read.csv("DataProcessed/new_cluster_results/daily/2025_06_25_daily_kmeans_cluster_labels.csv") %>%
  mutate(cdphe_flow_name = AREA_SEWERSHED) %>%
  select(cdphe_flow_name, kmeans_k2)

visit_with_cluster_daily <- left_join(visit, clusters_daily) %>%
  filter(cdphe_flow_name %in% cdphe_flow$cdphe_flow_name,
         date >= min(cdphe_flow$date),
         date <= max(cdphe_flow$date)) %>%
  mutate(cluster = ifelse(kmeans_k2 == 1, "High-Variation", "Low-Variation"))

visit_with_cluster_daily %>%
  group_by(sewershed_analysis_name) %>%
  slice(1) %>%
  group_by(cluster) %>%
  count()

spaghetti_high_daily <- visit_with_cluster_daily %>%
  filter(cluster == "High-Variation") %>%
  group_by(sewershed_analysis_name) %>%
  drop_na() %>%
  mutate(spaghetti_value = scale(log(total_devices_daily)))

spaghetti_low_daily <- visit_with_cluster_daily %>%
  filter(cluster == "Low-Variation") %>%
  group_by(sewershed_analysis_name) %>%
  drop_na() %>%
  mutate(spaghetti_value = scale(log(total_devices_daily)))

spaghetti_plot_high_daily <- ggplot(spaghetti_high_daily) +
  geom_line(aes(x = date, y = spaghetti_value, group = sewershed_analysis_name, color = sewershed_analysis_name),
            color = "orange2", alpha = 0.6) +
  labs(x = NULL, y = "log(Daily Device Count), Z-Score Normalized") +
  ggtitle("High-Variation Sewersheds (n = 14)") +
  scale_x_date(date_breaks = "4 months", date_labels = "%b %Y") +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 90, hjust = 0.95, vjust = 0.4),
        axis.title.y = element_text(margin = margin(l = 10, r = 5)),
        plot.margin = margin(10, 20, 10, 0, "pt"))

spaghetti_plot_low_daily <- ggplot(spaghetti_low_daily) +
  geom_line(aes(x = date, y = spaghetti_value, group = sewershed_analysis_name, color = sewershed_analysis_name),
            color = "darkorchid4", alpha = 0.6) +
  labs(x = NULL, y = "log(Daily Device Count), Z-Score Normalized") +
  ggtitle("Low-Variation Sewersheds (n = 52)") +
  scale_x_date(date_breaks = "4 months", date_labels = "%b %Y") +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 90, hjust = 0.95, vjust = 0.4),
        axis.title.y = element_text(margin = margin(l = 10, r = 5)),
        plot.margin = margin(10, 20, 10, 0, "pt"))

spaghetti_cluster_daily <- grid.arrange(spaghetti_plot_high_daily,
                                        spaghetti_plot_low_daily, ncol = 1)

# save spaghetti plot
ggsave("Figures/spaghetti_cluster_daily.png", height = 8, width = 10, plot = spaghetti_cluster_daily)


# get sample size by mobility cluster (flow)
cdphe_flow4 %>%
  group_by(sewershed_analysis_name) %>%
  slice(1) %>%
  group_by(cluster) %>%
  count()

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
  filter(cdphe_flow_name %in% reporting_start$cdphe_flow_name) #%>%

cdphe_flow <- cdphe_flow4 %>%
  group_by(sewershed_analysis_name) %>%
  mutate(ndays = as.numeric(difftime(max(date), min(date), "days"))) %>%
  filter(ndays >= 365) %>%
  left_join(visit_with_cluster) %>%
  # exclude the two values reported by the two high-variation sewersheds in December 2021
  mutate(date = case_when(cluster == "Low-Variation" ~ date,
                          cluster == "High-Variation" & date > as.Date("2022-01-01") ~ date,
                          TRUE ~ NA_Date_)) %>%
  drop_na(date)

save(cdphe_flow, file = "DataProcessed/cdphe_flow.Rda")

# make flow reporting timeline
timeline <- cdphe_flow %>%
  #left_join(visit_with_pop) %>%
  group_by(sewershed_analysis_name) %>%
  mutate(start_date = min(date)) %>%
  group_by(sewershed_analysis_name, start_date) %>%
  arrange(desc(start_date)) %>%
  ungroup() %>%
  mutate(id = match(sewershed_analysis_name, unique(sewershed_analysis_name))) %>%
  select(date, sewershed_result_name, cluster, id)

# sample size for reporting timeline is 69 sewersheds that have flow data, cluster
# information, and at least one year of data
unique(unlist(cdphe_flow$sewershed_result_name))

# get sample size by mobility cluster
cdphe_flow %>%
  group_by(sewershed_analysis_name) %>%
  slice(1) %>%
  group_by(cluster) %>%
  count()

# plot flow reporting timeline
flow_reporting_timeline <- ggplot(data = timeline) +
  geom_line(aes(x = date, y = sewershed_result_name, group = id, color = cluster), linewidth = 0.8) +
  labs(x = "Reporting Date Range", y = NULL) +
  scale_x_date(date_breaks = "4 months", date_labels = "%b %Y") +
  scale_y_discrete(limits = unique(timeline$sewershed_result_name)) +
  scale_color_manual(values = c("orange2", "darkorchid4"),
                     labels = c("High-Variation", "Low-Variation")) +
  guides(color = guide_legend(title = "Mobility Cluster")) +
  theme(plot.title = element_text(hjust = 0.5),
        axis.title.x = element_text(margin = margin(t = 5)),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 0.95, vjust = 0.4)); flow_reporting_timeline

ggsave("Figures/flow_reporting_timeline.png", width = 9, height = 11, plot = flow_reporting_timeline)

# plot flow reporting timeline (anonymized)
flow_reporting_timeline_anonymized <- ggplot(data = timeline) +
  geom_line(aes(x = date, y = id, group = id, color = cluster), linewidth = 0.8) +
  labs(x = "Reporting Date Range", y = NULL) +
  scale_x_date(date_breaks = "4 months", date_labels = "%b %Y") +
  scale_y_discrete(limits = factor(unique(timeline$id))) +
  scale_color_manual(values = c("orange2", "darkorchid4"),
                     labels = c("High-Variation", "Low-Variation")) +
  guides(color = guide_legend(title = "Mobility Cluster")) +
  theme(plot.title = element_text(hjust = 0.5),
        axis.title.x = element_text(margin = margin(t = 5)),
        axis.text.x = element_text(angle = 90, hjust = 0.95, vjust = 0.4)); flow_reporting_timeline_anonymized

ggsave("Figures/flow_reporting_timeline_anonymized.png", width = 9, height = 11, plot = flow_reporting_timeline_anonymized)

# exclude sewersheds that were not in the final flow dataset
# restrict start date for mobile device data to start of flow data
visit_with_cluster_sample <- visit_with_cluster %>%
  filter(cdphe_flow_name %in% cdphe_flow$cdphe_flow_name,
         date >= min(cdphe_flow$date),
         date <= max(cdphe_flow$date))

save(visit_with_cluster_sample, file = "DataProcessed/visit_with_cluster_sample.Rda")

# reorder levels for plotting
visit_with_cluster_sample$season <- factor(visit_with_cluster_sample$season, levels = c("Winter", "Spring", "Summer", "Fall"))
visit_with_cluster_sample$day_of_week <- factor(visit_with_cluster_sample$day_of_week, levels = c("Monday", "Tuesday", "Wednesday", "Thursday", "Friday", "Saturday", "Sunday"))

# median device count by day of week
visit_daily_stat_dayofweek <- visit_with_cluster_sample %>%
  drop_na(total_devices_daily) %>%
  group_by(cluster, day_of_week) %>%
  summarize(median_total_devices_daily = median(total_devices_daily),
            lower_total_devices_daily = quantile(total_devices_daily, 0.25),
            upper_total_devices_daily = quantile(total_devices_daily, 0.75))

plot_dayofweek_high <- ggplot(data = visit_daily_stat_dayofweek %>% filter(cluster == "High-Variation")) +
  geom_pointrange(aes(x = day_of_week, y = median_total_devices_daily,
                      ymin = lower_total_devices_daily,
                      ymax = upper_total_devices_daily)) +
  labs(x = NULL, y = "Daily Device Count (Median and IQR)") +
  ggtitle("High-Variation Sewersheds (n = 15)") +
  scale_y_continuous(labels = scales::comma) +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 90, hjust = 0.95, vjust = 0.4),
        axis.title.y = element_text(margin = margin(l = 10, r = 5)),
        plot.margin = margin(15, 20, 10, 0, "pt"))

plot_dayofweek_low <- ggplot(data = visit_daily_stat_dayofweek %>% filter(cluster == "Low-Variation")) +
  geom_pointrange(aes(x = day_of_week, y = median_total_devices_daily,
                      ymin = lower_total_devices_daily,
                      ymax = upper_total_devices_daily)) +
  labs(x = NULL, y = "Daily Device Count (Median and IQR)") +
  ggtitle("Low-Variation Sewersheds (n = 51)") +
  scale_y_continuous(labels = scales::comma, limits = c(0, 30000), breaks = seq(0, 30000, 10000)) +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 90, hjust = 0.95, vjust = 0.4),
        axis.title.y = element_text(margin = margin(l = 10, r = 5)),
        plot.margin = margin(15, 20, 10, 0, "pt"))

plot_dayofweek <- grid.arrange(plot_dayofweek_high, plot_dayofweek_low, ncol = 1)

ggsave("Figures/plot_dayofweek.png", height = 7, width = 7, plot = plot_dayofweek)

# median device count by month
visit_monthly_stat_month <- visit_with_cluster_sample %>%
  drop_na(total_devices_monthly) %>%
  group_by(cluster, month) %>%
  summarize(median_total_devices_monthly = median(total_devices_monthly),
            lower_total_devices_monthly = quantile(total_devices_monthly, 0.25),
            upper_total_devices_monthly = quantile(total_devices_monthly, 0.75),
            median_internal_devices_monthly = median(internal_devices_monthly))

plot_month_high <- ggplot(data = visit_monthly_stat_month %>% filter(cluster == "High-Variation")) +
  geom_pointrange(aes(x = month, y = median_total_devices_monthly,
                      ymin = lower_total_devices_monthly,
                      ymax = upper_total_devices_monthly)) +
  labs(x = NULL, y = "Monthly Device Count (Median and IQR)") +
  ggtitle("High-Variation Sewersheds (n = 15)") +
  scale_y_continuous(labels = scales::comma) +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 90, hjust = 0.95, vjust = 0.4),
        axis.title.y = element_text(margin = margin(l = 10, r = 5)),
        plot.margin = margin(15, 20, 10, 0, "pt"))

plot_month_low <- ggplot(data = visit_monthly_stat_month %>% filter(cluster == "Low-Variation")) +
  geom_pointrange(aes(x = month, y = median_total_devices_monthly,
                      ymin = lower_total_devices_monthly,
                      ymax = upper_total_devices_monthly)) +
  labs(x = NULL, y = "Monthly Device Count (Median and IQR)") +
  ggtitle("Low-Variation Sewersheds (n = 51)") +
  scale_y_continuous(labels = scales::comma) +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 90, hjust = 0.95, vjust = 0.4),
        axis.title.y = element_text(margin = margin(l = 10, r = 5)),
        plot.margin = margin(15, 20, 10, 0, "pt"))

# median flow rate by month
flow_stat_month <- cdphe_flow %>%
  left_join(visit_with_cluster) %>%
  drop_na(flow_rate) %>%
  drop_na(cluster) %>%
  group_by(cluster, month) %>%
  summarize(median_flow_monthly = median(flow_rate),
            lower_flow_monthly = quantile(flow_rate, 0.25),
            upper_flow_monthly = quantile(flow_rate, 0.75))

plot_month_flow_high <- ggplot(data = flow_stat_month %>% filter(cluster == "High-Variation")) +
  geom_pointrange(aes(x = month, y = median_flow_monthly,
                      ymin = lower_flow_monthly,
                      ymax = upper_flow_monthly)) +
  labs(x = NULL, y = "Monthly Flow Rate (Median and IQR)") +
  ggtitle("High-Variation Sewersheds (n = 15)") +
  scale_y_continuous(labels = scales::comma) +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 90, hjust = 0.95, vjust = 0.4),
        axis.title.y = element_text(margin = margin(l = 10, r = 5)),
        plot.margin = margin(15, 20, 10, 0, "pt"))

plot_month_flow_low <- ggplot(data = flow_stat_month %>% filter(cluster == "Low-Variation")) +
  geom_pointrange(aes(x = month, y = median_flow_monthly,
                      ymin = lower_flow_monthly,
                      ymax = upper_flow_monthly)) +
  labs(x = NULL, y = "Monthly Flow Rate (Median and IQR)") +
  ggtitle("Low-Variation Sewersheds (n = 51)") +
  scale_y_continuous(labels = scales::comma) +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 90, hjust = 0.95, vjust = 0.4),
        axis.title.y = element_text(margin = margin(l = 10, r = 5)),
        plot.margin = margin(15, 20, 10, 0, "pt"))

plot_month_mobile_flow <- grid.arrange(plot_month_high, plot_month_flow_high,
                                       plot_month_low, plot_month_flow_low, ncol = 2)

ggsave("Figures/plot_month_mobile_flow.png", height = 8, width = 14, plot = plot_month_mobile_flow)

# median device count by year
visit_monthly_stat_year <- visit_with_cluster_sample %>%
  drop_na(total_devices_monthly) %>%
  group_by(cluster, year) %>%
  summarize(median_total_devices_year = median(total_devices_monthly),
            lower_total_devices_year = quantile(total_devices_monthly, 0.25),
            upper_total_devices_year = quantile(total_devices_monthly, 0.75),
            median_internal_devices_year = median(internal_devices_monthly))

plot_year_high <- ggplot(data = visit_monthly_stat_year %>% filter(cluster == "High-Variation")) +
  geom_pointrange(aes(x = year, y = median_total_devices_year,
                      ymin = lower_total_devices_year,
                      ymax = upper_total_devices_year)) +
  labs(x = NULL, y = "Monthly Device Count (Median and IQR)") +
  ggtitle("High-Variation Sewersheds (n = 15)") +
  scale_y_continuous(labels = scales::comma) +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 90, hjust = 0.95, vjust = 0.4),
        axis.title.y = element_text(margin = margin(l = 10, r = 5)),
        plot.margin = margin(15, 20, 10, 0, "pt"))

plot_year_low <- ggplot(data = visit_monthly_stat_year %>% filter(cluster == "Low-Variation")) +
  geom_pointrange(aes(x = year, y = median_total_devices_year,
                      ymin = lower_total_devices_year,
                      ymax = upper_total_devices_year)) +
  labs(x = NULL, y = "Monthly Device Count (Median and IQR)") +
  ggtitle("Low-Variation Sewersheds (n = 51)") +
  scale_y_continuous(labels = scales::comma) +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 90, hjust = 0.95, vjust = 0.4),
        axis.title.y = element_text(margin = margin(l = 10, r = 5)),
        plot.margin = margin(15, 20, 10, 0, "pt"))

# median flow rate by year
flow_stat_year <- cdphe_flow %>%
  left_join(visit_with_cluster) %>%
  drop_na(flow_rate) %>%
  drop_na(cluster) %>%
  group_by(cluster, year) %>%
  summarize(median_flow_yearly = median(flow_rate),
            lower_flow_yearly = quantile(flow_rate, 0.25),
            upper_flow_yearly = quantile(flow_rate, 0.75))

plot_year_flow_high <- ggplot(data = flow_stat_year %>% filter(cluster == "High-Variation")) +
  geom_pointrange(aes(x = year, y = median_flow_yearly,
                      ymin = lower_flow_yearly,
                      ymax = upper_flow_yearly)) +
  labs(x = NULL, y = "Yearly Flow Rate (Median and IQR)") +
  ggtitle("High-Variation Sewersheds (n = 15)") +
  scale_y_continuous(labels = scales::comma) +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 90, hjust = 0.95, vjust = 0.4),
        axis.title.y = element_text(margin = margin(l = 10, r = 5)),
        plot.margin = margin(15, 20, 10, 0, "pt"))

plot_year_flow_low <- ggplot(data = flow_stat_year %>% filter(cluster == "Low-Variation")) +
  geom_pointrange(aes(x = year, y = median_flow_yearly,
                      ymin = lower_flow_yearly,
                      ymax = upper_flow_yearly)) +
  labs(x = NULL, y = "Yearly Flow Rate (Median and IQR)") +
  ggtitle("Low-Variation Sewersheds (n = 51)") +
  scale_y_continuous(labels = scales::comma) +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 90, hjust = 0.95, vjust = 0.4),
        axis.title.y = element_text(margin = margin(l = 10, r = 5)),
        plot.margin = margin(15, 20, 10, 0, "pt"))

plot_year_mobile_flow <- grid.arrange(plot_year_high, plot_year_flow_high,
                                      plot_year_low, plot_year_flow_low, ncol = 2)

ggsave("Figures/plot_year_mobile_flow.png", height = 7, width = 14, plot = plot_year_mobile_flow)

# log-transformed daily device counts
log_df_daily <- visit_with_cluster_sample %>%
  mutate(total_devices_daily = ifelse(total_devices_daily == 0, 1, total_devices_daily),
         total_devices_daily_log = log(total_devices_daily))

hist_raw_daily_high <- ggplot(data = log_df_daily %>% filter(cluster == "High-Variation")) +
  geom_histogram(aes(x = total_devices_daily)) +
  labs(x = "Daily Device Count", y = "Frequency (Sewershed-Days)") +
  ggtitle("High-Variation Sewersheds") +
  scale_x_continuous(labels = scales::comma) +
  scale_y_continuous(labels = scales::comma) +
  theme(plot.title = element_text(hjust = 0.5))

hist_log_daily_high <- ggplot(data = log_df_daily %>% filter(cluster == "High-Variation")) +
  geom_histogram(aes(x = total_devices_daily_log)) +
  labs(x = "log(Daily Device Count)", y = "Frequency (Sewershed-Days)") +
  ggtitle("High-Variation Sewersheds") +
  scale_x_continuous(labels = scales::comma) +
  scale_y_continuous(labels = scales::comma) +
  theme(plot.title = element_text(hjust = 0.5))

qq_daily_high <- ggplot(data = log_df_daily %>% filter(cluster == "High-Variation"),
                        aes(sample = total_devices_daily_log)) +
  stat_qq() +
  stat_qq_line() +
  labs(x = "Theoretical Quantiles", y = "Sample Quantiles") +
  ggtitle("High-Variation Sewersheds") +
  theme(plot.title = element_text(hjust = 0.5))

hist_raw_daily_low <- ggplot(data = log_df_daily %>% filter(cluster == "Low-Variation")) +
  geom_histogram(aes(x = total_devices_daily)) +
  labs(x = "Daily Device Count", y = "Frequency (Sewershed-Days)") +
  ggtitle("Low-Variation Sewersheds") +
  scale_x_continuous(labels = scales::comma) +
  scale_y_continuous(labels = scales::comma) +
  theme(plot.title = element_text(hjust = 0.5))

hist_log_daily_low <- ggplot(data = log_df_daily %>% filter(cluster == "Low-Variation")) +
  geom_histogram(aes(x = total_devices_daily_log)) +
  labs(x = "log(Daily Device Count)", y = "Frequency (Sewershed-Days)") +
  ggtitle("Low-Variation Sewersheds") +
  scale_x_continuous(labels = scales::comma) +
  scale_y_continuous(labels = scales::comma) +
  theme(plot.title = element_text(hjust = 0.5))

qq_daily_low <- ggplot(data = log_df_daily %>% filter(cluster == "Low-Variation"),
                        aes(sample = total_devices_daily_log)) +
  stat_qq() +
  stat_qq_line() +
  labs(x = "Theoretical Quantiles", y = "Sample Quantiles") +
  ggtitle("Low-Variation Sewersheds") +
  theme(plot.title = element_text(hjust = 0.5))

raw_log_qq_daily <- grid.arrange(hist_raw_daily_high, hist_log_daily_high, qq_daily_high,
                                 hist_raw_daily_low, hist_log_daily_low, qq_daily_low, ncol = 3)

ggsave("Figures/raw_log_qq_daily.png", height = 8, width = 14, plot = raw_log_qq_daily)

# log-transformed monthly device counts
log_df_monthly <- visit_with_cluster_sample %>%
  drop_na(total_devices_monthly) %>%
  mutate(total_devices_monthly_log = log(total_devices_monthly),
         internal_devices_monthly_log = log(internal_devices_monthly))

hist_raw_monthly_high <- ggplot(data = log_df_monthly %>% filter(cluster == "High-Variation")) +
  geom_histogram(aes(x = total_devices_monthly)) +
  labs(x = "Monthly Device Count", y = "Frequency (Sewershed-Months)") +
  ggtitle("High-Variation Sewersheds") +
  scale_x_continuous(labels = scales::comma) +
  scale_y_continuous(labels = scales::comma) +
  theme(plot.title = element_text(hjust = 0.5))

hist_log_monthly_high <- ggplot(data = log_df_monthly %>% filter(cluster == "High-Variation")) +
  geom_histogram(aes(x = total_devices_monthly_log)) +
  labs(x = "log(Monthly Device Count)", y = "Frequency (Sewershed-Months)") +
  ggtitle("High-Variation Sewersheds") +
  scale_x_continuous(labels = scales::comma) +
  scale_y_continuous(labels = scales::comma) +
  theme(plot.title = element_text(hjust = 0.5))

qq_monthly_high <- ggplot(data = log_df_monthly %>% filter(cluster == "High-Variation"),
                          aes(sample = total_devices_monthly_log)) +
  stat_qq() +
  stat_qq_line() +
  labs(x = "Theoretical Quantiles", y = "Sample Quantiles") +
  ggtitle("High-Variation Sewersheds") +
  theme(plot.title = element_text(hjust = 0.5))

hist_raw_monthly_low <- ggplot(data = log_df_monthly %>% filter(cluster == "Low-Variation")) +
  geom_histogram(aes(x = total_devices_monthly)) +
  labs(x = "Monthly Device Count", y = "Frequency (Sewershed-Months)") +
  ggtitle("Low-Variation Sewersheds") +
  scale_x_continuous(labels = scales::comma) +
  scale_y_continuous(labels = scales::comma) +
  theme(plot.title = element_text(hjust = 0.5))

hist_log_monthly_low <- ggplot(data = log_df_monthly %>% filter(cluster == "Low-Variation")) +
  geom_histogram(aes(x = total_devices_monthly_log)) +
  labs(x = "log(Monthly Device Count)", y = "Frequency (Sewershed-Months)") +
  ggtitle("Low-Variation Sewersheds") +
  scale_x_continuous(labels = scales::comma) +
  scale_y_continuous(labels = scales::comma) +
  theme(plot.title = element_text(hjust = 0.5))

qq_monthly_low <- ggplot(data = log_df_monthly %>% filter(cluster == "Low-Variation"),
                          aes(sample = total_devices_monthly_log)) +
  stat_qq() +
  stat_qq_line() +
  labs(x = "Theoretical Quantiles", y = "Sample Quantiles") +
  ggtitle("Low-Variation Sewersheds") +
  theme(plot.title = element_text(hjust = 0.5))

raw_log_qq_monthly <- grid.arrange(hist_raw_monthly_high, hist_log_monthly_high, qq_monthly_high,
                                   hist_raw_monthly_low, hist_log_monthly_low, qq_monthly_low, ncol = 3)

ggsave("Figures/raw_log_qq_monthly.png", height = 8, width = 14, plot = raw_log_qq_monthly)

# log-transformed flow rates
cdphe_flow_log <- cdphe_flow %>%
  drop_na(flow_rate) %>%
  mutate(flow_rate = ifelse(flow_rate == 0, 0.001, flow_rate),
         flow_rate_log = log(flow_rate)) %>%
  filter(day_of_week %in% c("Monday", "Tuesday", "Wednesday", "Thursday", "Friday"))

hist_raw_flow_high <- ggplot(data = cdphe_flow_log %>% filter(cluster == "High-Variation")) +
  geom_histogram(aes(x = flow_rate)) +
  labs(x = "Flow Rate", y = "Frequency (Sewershed-Days)") +
  ggtitle("High-Variation Sewersheds") +
  scale_x_continuous(labels = scales::comma) +
  scale_y_continuous(labels = scales::comma) +
  theme(plot.title = element_text(hjust = 0.5)); hist_raw_flow_high

hist_log_flow_high <- ggplot(data = cdphe_flow_log %>% filter(cluster == "High-Variation")) +
  geom_histogram(aes(x = flow_rate_log)) +
  labs(x = "log(Flow Rate)", y = "Frequency (Sewershed-Days)") +
  ggtitle("High-Variation Sewersheds") +
  scale_x_continuous(labels = scales::comma) +
  scale_y_continuous(labels = scales::comma) +
  theme(plot.title = element_text(hjust = 0.5)); hist_log_flow_high

qq_flow_high <- ggplot(data = cdphe_flow_log %>% filter(cluster == "High-Variation"),
                       aes(sample = flow_rate_log)) +
  stat_qq() +
  stat_qq_line() +
  labs(x = "Theoretical Quantiles", y = "Sample Quantiles") +
  ggtitle("High-Variation Sewersheds") +
  theme(plot.title = element_text(hjust = 0.5)); qq_flow_high

hist_raw_flow_low <- ggplot(data = cdphe_flow_log %>% filter(cluster == "Low-Variation")) +
  geom_histogram(aes(x = flow_rate)) +
  labs(x = "Flow Rate", y = "Frequency (Sewershed-Days)") +
  ggtitle("Low-Variation Sewersheds") +
  scale_x_continuous(labels = scales::comma) +
  scale_y_continuous(labels = scales::comma) +
  theme(plot.title = element_text(hjust = 0.5)); hist_raw_flow_low

hist_log_flow_low <- ggplot(data = cdphe_flow_log %>% filter(cluster == "Low-Variation")) +
  geom_histogram(aes(x = flow_rate_log)) +
  labs(x = "log(Flow Rate)", y = "Frequency (Sewershed-Days)") +
  ggtitle("Low-Variation Sewersheds") +
  scale_x_continuous(labels = scales::comma) +
  scale_y_continuous(labels = scales::comma) +
  theme(plot.title = element_text(hjust = 0.5)); hist_log_flow_low

qq_flow_low <- ggplot(data = cdphe_flow_log %>% filter(cluster == "Low-Variation"),
                      aes(sample = flow_rate_log)) +
  stat_qq() +
  stat_qq_line() +
  labs(x = "Theoretical Quantiles", y = "Sample Quantiles") +
  ggtitle("Low-Variation Sewersheds") +
  theme(plot.title = element_text(hjust = 0.5)); qq_flow_low

raw_log_qq_flow <- grid.arrange(hist_raw_flow_high, hist_log_flow_high, qq_flow_high,
                                hist_raw_flow_low, hist_log_flow_low, qq_flow_low, ncol = 3)

ggsave("Figures/raw_log_qq_flow.png", height = 8, width = 14, plot = raw_log_qq_flow)

# calculate mean and sd for device counts
calc_dayofweek <- visit_with_cluster_sample %>%
  group_by(cluster, day_of_week) %>%
  summarize(mean = round(mean(total_devices_daily), 0),
            sd = round(sd(total_devices_daily), 0)) %>%
  arrange(cluster, factor(day_of_week, levels = c("Monday", "Tuesday", "Wednesday", "Thursday", "Friday", "Saturday", "Sunday")))
write.csv(calc_dayofweek, "DataProcessed/calc_dayofweek.csv")

calc_wdwe <- visit_with_cluster_sample %>%
  group_by(cluster, weekday_weekend) %>%
  summarize(mean = round(mean(total_devices_daily), 0),
            sd = round(sd(total_devices_daily), 0))
write.csv(calc_wdwe, "DataProcessed/calc_wdwe.csv")

calc_month <- visit_with_cluster_sample %>%
  drop_na(total_devices_monthly) %>%
  group_by(cluster, month) %>%
  summarize(mean = round(mean(total_devices_monthly), 0),
            sd = round(sd(total_devices_monthly), 0))
write.csv(calc_month, "DataProcessed/calc_month.csv")

calc_season <- visit_with_cluster_sample %>%
  drop_na(total_devices_monthly) %>%
  group_by(cluster, season) %>%
  summarize(mean = round(mean(total_devices_monthly)),
            sd = round(sd(total_devices_monthly)))
write.csv(calc_season, "DataProcessed/calc_season.csv")

calc_year <- visit_with_cluster_sample %>%
  drop_na(total_devices_monthly) %>%
  group_by(cluster, year) %>%
  summarize(mean = round(mean(total_devices_monthly), 0),
            sd = round(sd(total_devices_monthly), 0))
write.csv(calc_year, "DataProcessed/calc_year.csv")

# linear mixed effect models and LRT on log-transformed device counts
mod_dayofweek_high <- lmer(total_devices_daily_log ~ 1 + (1 | sewershed_analysis_name), data = log_df_daily %>% filter(cluster == "High-Variation"), REML = F)
mod1_dayofweek_high <- lmer(total_devices_daily_log ~ 1 + day_of_week + (1 | sewershed_analysis_name), data = log_df_daily %>% filter(cluster == "High-Variation"), REML = F)
anova(mod_dayofweek_high, mod1_dayofweek_high, test = "LRT")

mod_dayofweek_low <- lmer(total_devices_daily_log ~ 1 + (1 | sewershed_analysis_name), data = log_df_daily %>% filter(cluster == "Low-Variation"), REML = F)
mod1_dayofweek_low <- lmer(total_devices_daily_log ~ 1 + day_of_week + (1 | sewershed_analysis_name), data = log_df_daily %>% filter(cluster == "Low-Variation"), REML = F)
anova(mod_dayofweek_low, mod1_dayofweek_low, test = "LRT")

mod_wdwe_high <- lmer(total_devices_daily_log ~ 1 + (1 | sewershed_analysis_name), data = log_df_daily %>% filter(cluster == "High-Variation"), REML = F)
mod1_wdwe_high <- lmer(total_devices_daily_log ~ 1 + weekday_weekend + (1 | sewershed_analysis_name), data = log_df_daily %>% filter(cluster == "High-Variation"), REML = F)
anova(mod_wdwe_high, mod1_wdwe_high, test = "LRT")

mod_wdwe_low <- lmer(total_devices_daily_log ~ 1 + (1 | sewershed_analysis_name), data = log_df_daily %>% filter(cluster == "Low-Variation"), REML = F)
mod1_wdwe_low <- lmer(total_devices_daily_log ~ 1 + weekday_weekend + (1 | sewershed_analysis_name), data = log_df_daily %>% filter(cluster == "Low-Variation"), REML = F)
anova(mod_wdwe_low, mod1_wdwe_low, test = "LRT")

mod_month_high <- lmer(total_devices_monthly_log ~ 1 + (1 | sewershed_analysis_name), data = log_df_monthly %>% filter(cluster == "High-Variation"), REML = F)
mod1_month_high <- lmer(total_devices_monthly_log ~ 1 + month + (1 | sewershed_analysis_name), data = log_df_monthly %>% filter(cluster == "High-Variation"), REML = F)
anova(mod_month_high, mod1_month_high, test = "LRT")

mod_month_low <- lmer(total_devices_monthly_log ~ 1 + (1 | sewershed_analysis_name), data = log_df_monthly %>% filter(cluster == "Low-Variation"), REML = F)
mod1_month_low <- lmer(total_devices_monthly_log ~ 1 + month + (1 | sewershed_analysis_name), data = log_df_monthly %>% filter(cluster == "Low-Variation"), REML = F)
anova(mod_month_low, mod1_month_low, test = "LRT")

mod_season_high <- lmer(total_devices_monthly_log ~ 1 + (1 | sewershed_analysis_name), data = log_df_monthly %>% filter(cluster == "High-Variation"), REML = F)
mod1_season_high <- lmer(total_devices_monthly_log ~ 1 + season + (1 | sewershed_analysis_name), data = log_df_monthly %>% filter(cluster == "High-Variation"), REML = F)
anova(mod_season_high, mod1_season_high, test = "LRT")

mod_season_low <- lmer(total_devices_monthly_log ~ 1 + (1 | sewershed_analysis_name), data = log_df_monthly %>% filter(cluster == "Low-Variation"), REML = F)
mod1_season_low <- lmer(total_devices_monthly_log ~ 1 + season + (1 | sewershed_analysis_name), data = log_df_monthly %>% filter(cluster == "Low-Variation"), REML = F)
anova(mod_season_low, mod1_season_low, test = "LRT")

mod_year_high <- lmer(total_devices_monthly_log ~ 1 + (1 | sewershed_analysis_name), data = log_df_monthly %>% filter(cluster == "High-Variation"), REML = F)
mod1_year_high <- lmer(total_devices_monthly_log ~ 1 + year + (1 | sewershed_analysis_name), data = log_df_monthly %>% filter(cluster == "High-Variation"), REML = F)
anova(mod_year_high, mod1_year_high, test = "LRT")

mod_year_low <- lmer(total_devices_monthly_log ~ 1 + (1 | sewershed_analysis_name), data = log_df_monthly %>% filter(cluster == "Low-Variation"), REML = F)
mod1_year_low <- lmer(total_devices_monthly_log ~ 1 + year + (1 | sewershed_analysis_name), data = log_df_monthly %>% filter(cluster == "Low-Variation"), REML = F)
anova(mod_year_low, mod1_year_low, test = "LRT")

# sensitivity analysis on internal devices
calc_month_internal <- visit_with_cluster_sample %>%
  drop_na(total_devices_monthly) %>%
  group_by(cluster, month) %>%
  summarize(mean = round(mean(internal_devices_monthly), 0),
            sd = round(sd(internal_devices_monthly), 0))
write.csv(calc_month_internal, "DataProcessed/calc_month_internal.csv")

calc_season_internal <- visit_with_cluster_sample %>%
  drop_na(total_devices_monthly) %>%
  group_by(cluster, season) %>%
  summarize(mean = round(mean(internal_devices_monthly), 0),
            sd = round(sd(internal_devices_monthly), 0)) %>%
  arrange(cluster, factor(season, levels = c("Winter", "Spring", "Summer", "Fall")))
write.csv(calc_season_internal, "DataProcessed/calc_season_internal.csv")

calc_year_internal <- visit_with_cluster_sample %>%
  drop_na(total_devices_monthly) %>%
  group_by(cluster, year) %>%
  summarize(mean = round(mean(internal_devices_monthly), 0),
            sd = round(sd(internal_devices_monthly), 0))
write.csv(calc_year_internal, "DataProcessed/calc_year_internal.csv")

mod_month_high_internal <- lmer(internal_devices_monthly_log ~ 1 + (1 | sewershed_analysis_name), data = log_df_monthly %>% filter(cluster == "High-Variation"), REML = F)
mod1_month_high_internal <- lmer(internal_devices_monthly_log ~ 1 + month + (1 | sewershed_analysis_name), data = log_df_monthly %>% filter(cluster == "High-Variation"), REML = F)
anova(mod_month_high_internal, mod1_month_high_internal, test = "LRT")

mod_month_low_internal <- lmer(internal_devices_monthly_log ~ 1 + (1 | sewershed_analysis_name), data = log_df_monthly %>% filter(cluster == "Low-Variation"), REML = F)
mod1_month_low_internal <- lmer(internal_devices_monthly_log ~ 1 + month + (1 | sewershed_analysis_name), data = log_df_monthly %>% filter(cluster == "Low-Variation"), REML = F)
anova(mod_month_low_internal, mod1_month_low_internal, test = "LRT")

mod_season_high_internal <- lmer(internal_devices_monthly_log ~ 1 + (1 | sewershed_analysis_name), data = log_df_monthly %>% filter(cluster == "High-Variation"), REML = F)
mod1_season_high_internal <- lmer(internal_devices_monthly_log ~ 1 + season + (1 | sewershed_analysis_name), data = log_df_monthly %>% filter(cluster == "High-Variation"), REML = F)
anova(mod_season_high_internal, mod1_season_high_internal, test = "LRT")

mod_season_low_internal <- lmer(internal_devices_monthly_log ~ 1 + (1 | sewershed_analysis_name), data = log_df_monthly %>% filter(cluster == "Low-Variation"), REML = F)
mod1_season_low_internal <- lmer(internal_devices_monthly_log ~ 1 + season + (1 | sewershed_analysis_name), data = log_df_monthly %>% filter(cluster == "Low-Variation"), REML = F)
anova(mod_season_low_internal, mod1_season_low_internal, test = "LRT")

mod_year_high_internal <- lmer(internal_devices_monthly_log ~ 1 + (1 | sewershed_analysis_name), data = log_df_monthly %>% filter(cluster == "High-Variation"), REML = F)
mod1_year_high_internal <- lmer(internal_devices_monthly_log ~ 1 + year + (1 | sewershed_analysis_name), data = log_df_monthly %>% filter(cluster == "High-Variation"), REML = F)
anova(mod_year_high_internal, mod1_year_high_internal, test = "LRT")

mod_year_low_internal <- lmer(internal_devices_monthly_log ~ 1 + (1 | sewershed_analysis_name), data = log_df_monthly %>% filter(cluster == "Low-Variation"), REML = F)
mod1_year_low_internal <- lmer(internal_devices_monthly_log ~ 1 + year + (1 | sewershed_analysis_name), data = log_df_monthly %>% filter(cluster == "Low-Variation"), REML = F)
anova(mod_year_low_internal, mod1_year_low_internal, test = "LRT")

# calculate mean and sd for flow rates
calc_month_flow <- cdphe_flow %>%
  group_by(cluster, month) %>%
  summarize(mean = signif(mean(flow_rate), 3),
            sd = signif(sd(flow_rate), 3))
write.csv(calc_month_flow, "DataProcessed/calc_month_flow.csv")

calc_season_flow <- cdphe_flow %>%
  group_by(cluster, season) %>%
  summarize(mean = signif(mean(flow_rate), 3),
            sd = signif(sd(flow_rate), 3)) %>%
  arrange(cluster, factor(season, levels = c("Winter", "Spring", "Summer", "Fall")))
write.csv(calc_season_flow, "DataProcessed/calc_season_flow.csv")

calc_year_flow <- cdphe_flow %>%
  group_by(cluster, year) %>%
  summarize(mean = signif(mean(flow_rate), 3),
            sd = signif(sd(flow_rate), 3))
write.csv(calc_year_flow, "DataProcessed/calc_year_flow.csv")

mod_month_flow_high <-  lmer(flow_rate_log ~ 1 +               (1 | sewershed_analysis_name), data = cdphe_flow_log %>% filter(cluster == "High-Variation"), REML = F)
mod1_month_flow_high <- lmer(flow_rate_log ~ 1 + month + (1 | sewershed_analysis_name), data = cdphe_flow_log %>% filter(cluster == "High-Variation"), REML = F)
anova(mod_month_flow_high, mod1_month_flow_high, test = "LRT")

mod_month_flow_low <- lmer(flow_rate_log ~ 1 + (1 | sewershed_analysis_name), data = cdphe_flow_log %>% filter(cluster == "Low-Variation"), REML = F)
mod1_month_flow_low <- lmer(flow_rate_log ~ 1 + month + (1 | sewershed_analysis_name), data = cdphe_flow_log %>% filter(cluster == "Low-Variation"), REML = F)
anova(mod_month_flow_low, mod1_month_flow_low, test = "LRT")

mod_season_flow_high <-  lmer(flow_rate_log ~ 1 +               (1 | sewershed_analysis_name), data = cdphe_flow_log %>% filter(cluster == "High-Variation"), REML = F)
mod1_season_flow_high <- lmer(flow_rate_log ~ 1 + season + (1 | sewershed_analysis_name), data = cdphe_flow_log %>% filter(cluster == "High-Variation"), REML = F)
anova(mod_season_flow_high, mod1_season_flow_high, test = "LRT")

mod_season_flow_low <- lmer(flow_rate_log ~ 1 + (1 | sewershed_analysis_name), data = cdphe_flow_log %>% filter(cluster == "Low-Variation"), REML = F)
mod1_season_flow_low <- lmer(flow_rate_log ~ 1 + season + (1 | sewershed_analysis_name), data = cdphe_flow_log %>% filter(cluster == "Low-Variation"), REML = F)
anova(mod_season_flow_low, mod1_season_flow_low, test = "LRT")

mod_year_flow_high <-  lmer(flow_rate_log ~ 1 +               (1 | sewershed_analysis_name), data = cdphe_flow_log %>% filter(cluster == "High-Variation"), REML = F)
mod1_year_flow_high <- lmer(flow_rate_log ~ 1 + year + (1 | sewershed_analysis_name), data = cdphe_flow_log %>% filter(cluster == "High-Variation"), REML = F)
anova(mod_year_flow_high, mod1_year_flow_high, test = "LRT")

mod_year_flow_low <- lmer(flow_rate_log ~ 1 + (1 | sewershed_analysis_name), data = cdphe_flow_log %>% filter(cluster == "Low-Variation"), REML = F)
mod1_year_flow_low <- lmer(flow_rate_log ~ 1 + year + (1 | sewershed_analysis_name), data = cdphe_flow_log %>% filter(cluster == "Low-Variation"), REML = F)
anova(mod_year_flow_low, mod1_year_flow_low, test = "LRT")

########################################################################################################
