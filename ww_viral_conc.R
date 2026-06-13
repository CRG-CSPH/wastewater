# set working directory
setwd("/Users/emwu9912/Documents/CU Anschutz/Wastewater/")

# turn off scientific notation
options(scipen = 999)

library(pacman)
p_load(car, DHARMa, digest, FSA, glmmTMB, lme4, lmerTest, purrr, scales, sf, slider, tidycensus, tidyverse, tigris, zoo)

theme <- theme(plot.title = element_text(hjust = 0.5),
               #axis.title.x = element_text(margin = margin(t = 5)),
               #axis.ticks.y = element_blank(),
               axis.text.x = element_text(angle = 90),
               strip.text = element_text(size = 12))

ww_capacity <- read.csv("DataRaw/ww_capacity.csv") %>%
  rename(utility = wwtp_name)

sentinel <- read.csv("DataRaw/sentinel.csv")

sewershed_names <- read.csv("DataRaw/sewershed_names.csv")
sewershed_populations <- read.csv("DataProcessed/sewershed_populations.csv")


# overlap <- read.csv("DataRaw/new_overlap_jan2025.csv") %>%
#   rename(cdphe_flow_name = wwtp) %>%
#   left_join(sewershed_names) %>%
#   rename(utility = sewershed_sentinel_name)
# 
# assigned_cbgs <- overlap %>%
#   mutate(GEOID = as.character(GEOID)) %>%
#   group_by(GEOID) %>%
#   arrange(desc(PERCENTAGE)) %>%
#   slice(1) %>%
#   filter(PERCENTAGE >= 0.2) %>%
#   group_by(utility) %>%
#   add_count() %>%
#   filter(utility %in% sentinel$utility)
#
# cbg_2019 <- get_acs(geography = "block group",
#                     variables = "B01003_001",
#                     state = "CO",
#                     year = 2019) %>%
#   mutate(GEOID = str_remove(GEOID, "^0+")) %>%
#   rename(pop_2019 = estimate) %>%
#   select(GEOID, pop_2019)
# 
# sewershed_populations <- assigned_cbgs %>%
#   left_join(cbg_2019) %>%
#   group_by(utility) %>%
#   rename(n_cbgs = n) %>%
#   mutate(census_population = sum(pop_2019)) %>%
#   slice(1) %>%
#   select(utility, n_cbgs, census_population)

#write.csv(sewershed_populations, "DataProcessed/sewershed_populations.csv")
sentinel_with_populations <- left_join(sentinel, sewershed_populations) %>%
  mutate(category_id = case_when(category == "Metro" ~ "ME",
                                 category == "Other Large" ~ "OL",
                                 category == "High Tourism" ~ "HT",
                                 TRUE ~ "OS")) %>%
  arrange(match(category, c("ME", "OL", "HT", "OS")))

category_order <- unique(sentinel_with_populations$category)

# set seed for reproducibility
set.seed(42)
sentinel_with_populations_random <- sentinel_with_populations %>%
  mutate(category_ordered = factor(category, levels = category_order)) %>%
  group_by(category_ordered) %>%
  slice_sample(prop = 1) %>%
  mutate(sewershed_id = row_number(),
         sewershed_code = paste0(category_id, sewershed_id)) %>%
  left_join(ww_capacity) %>%
  mutate(TOTAL_DESIGN_FLOW_NMBR_10 = TOTAL_DESIGN_FLOW_NMBR*1.1)

utility_order <- unique(sentinel_with_populations_random$utility)
sewershed_code_order <- unique(sentinel_with_populations_random$sewershed_code)
#category_order <- unique(sentinel_with_populations_random$category)

targets <- c(FLUA = "fluav", RSV = "rsv", COVID = "sars-cov-2")
titles <- c("Influenza A", "RSV", "SARS-CoV-2")
pcr_targets <- c("fluav", "rsv", "sars-cov-2")

lp3_start <- "2024-09-30"

visit <- read.csv("DataProcessed/Mobile Device Data/2026_04_20_daily_visits.csv") %>%
  rename(cdphe_flow_name = DESTINATION_AREA_SEWERSHED,
         total_devices_daily = STOPS_BY_DAY_L) %>%
  mutate(date = as.Date(DAY)) %>%
  left_join(sewershed_names) %>%
  #left_join(sentinel) %>%
  rename(utility = sewershed_sentinel_name) %>%
  left_join(sentinel_with_populations_random) %>%
  filter(utility %in% sentinel_with_populations_random$utility,
         date >= lp3_start)

visit_weekly <- visit %>%
  mutate(week_end = ceiling_date(date, "weeks", week_start = getOption("lubridate.week.start", 6))) %>%
  group_by(utility, week_end) %>%
  mutate(average_devices = signif(mean(total_devices_daily), 4)) %>%
  select(-total_devices_daily) %>%
  left_join(sentinel_with_populations_random) %>%
  mutate(sewershed_code_ordered = factor(sewershed_code, levels = sewershed_code_order))

ww1 <- read.csv("DataRaw/2026/CSPH_WWData_2026-02-26.csv") %>%
  mutate(date = as.Date(sample_collect_date, "%m/%d/%y"),
         dataset = 1) %>%
  group_by(wwtp_name, sample_collect_date, pcr_target) %>%
  # deduplicate
  slice(1)

ww2 <- read.csv("DataRaw/CSPH Data _ April 2026/CSPH_WWData_2026-04-13.csv") %>%
  mutate(date = as.Date(sample_collect_date, "%m/%d/%y"),
         dataset = 2) %>%
  group_by(wwtp_name, sample_collect_date, pcr_target) %>%
  # deduplicate
  slice(1)

ww <- bind_rows(ww1, ww2) %>%
  rename(utility = wwtp_name) %>%
  mutate(lod_sewage = 1200,
         pcr_target = tolower(pcr_target),
         # indicator variable for detect/non-detect
         detect = ifelse(pcr_target_avg_conc <= lod_sewage, 0, 1),
         # impute non-detects with half the LOD
         pcr_target_avg_conc = ifelse(detect == 0, lod_sewage/2, pcr_target_avg_conc)) %>%
  filter(utility %in% sentinel$utility,
         pcr_target %in% c("fluav", "rsv", "sars-cov-2")) %>%
  group_by(utility, pcr_target, date) %>%
  # if repeated values from date overlap, select the more recent dataset
  arrange(desc(dataset)) %>%
  slice(1) %>%
  left_join(sentinel_with_populations)

# proportion of observations with missing or zero flow rates
percent(sum(is.na(ww$flow_rate) | ww$flow_rate == 0, na.rm = TRUE)/nrow(ww), accuracy = 0.01)

# impute missing or zero flow rates with the average of the last two non-missing nonzero flow rates
ww_na_flow <- ww %>%
  group_by(utility, pcr_target) %>%
  mutate(flow_rate = na_if(flow_rate, 0),
         flow_rate_last2avg = slide_dbl(flow_rate,
                                        ~ mean(tail(na.omit(.x), 2)),
                                        .before = Inf,
                                        .complete = F),
         flow_rate = ifelse(is.na(flow_rate), flow_rate_last2avg, flow_rate)) %>%
  select(-flow_rate_last2avg, -sample_collect_date, -population_served) %>%
  left_join(sentinel_with_populations_random)

# proportion of observations where flow rate value exceeds max capacity
percent(length(which(ww_na_flow$flow_rate > ww_na_flow$TOTAL_DESIGN_FLOW_NMBR_10))/nrow(ww_na_flow), accuracy = 0.01)

# by sewershed
above_capacity <- ww_na_flow %>%
  filter(pcr_target == "sars-cov-2") %>%
  mutate(above_capacity_flag = ifelse(flow_rate > TOTAL_DESIGN_FLOW_NMBR_10, "yes", "no")) %>%
  group_by(utility) %>%
  mutate(n_above = sum(above_capacity_flag == "yes"),
         percent_above = round(n_above/n()*100, 1)) %>%
  slice(1) %>%
  select(utility, size, category, TOTAL_DESIGN_FLOW_NMBR_10, percent_above)

ww_impute_flow <- ww_na_flow %>%
  mutate(flow_rate = ifelse(flow_rate > TOTAL_DESIGN_FLOW_NMBR_10, TOTAL_DESIGN_FLOW_NMBR_10, flow_rate),
         week_end = ceiling_date(date, "weeks", week_start = getOption("lubridate.week.start", 6))) %>%
  group_by(utility, pcr_target, week_end) %>%
  add_count() %>%
  rename(samp_freq = n) %>%
  mutate(sewershed_code_ordered = factor(sewershed_code, levels = sewershed_code_order))
  
write.csv(ww_impute_flow, "DataProcessed/ww_impute_flow.csv") # output data for Python WVAL code

sampling_freq <- ww_impute_flow %>%
  group_by(utility, pcr_target, week_end) %>%
  #arrange(category_id, sewershed_id) %>%
  slice(1)

sample_freq_plot_anon <- ggplot(data = sampling_freq %>% filter(pcr_target == "sars-cov-2")) +
  geom_bar(aes(x = date, y = samp_freq), stat = "identity") +
  labs(x = NULL, y = "Number of Samples Per Week") +
  scale_y_continuous(limits = c(0, 3), breaks = seq(0, 3, 1)) +
  scale_x_date(breaks = "4 months", date_labels = "%b %d %Y") +
  facet_wrap(~sewershed_code_ordered, scales = "free_y") +
  theme; sample_freq_plot_anon #+
  #theme(axis.text.x = element_text(angle = 90),
        #strip.text = element_text(size = 12)); sample_freq_plot_anon
ggsave("Figures/sample_freq_plot_anon.png", height = 7, width = 12, plot = sample_freq_plot_anon)

flow_rate_over_time_anon <- ggplot(data = ww_impute_flow %>% filter(pcr_target == "fluav")) +
  geom_point(aes(x = date, y = flow_rate*3.785, color = sample_type)) +
  geom_line(aes(x = date, y = flow_rate*3.785)) +
  labs(x = NULL, y = "Daily Flow Rate (MLD)", color = "Sample Type") +
  scale_y_continuous(labels = scales::comma) +
  scale_x_date(breaks = "4 months", date_labels = "%b %d %Y") +
  facet_wrap(~sewershed_code_ordered, scales = "free_y") +
  theme; flow_rate_over_time_anon
ggsave("Figures/flow_rate_over_time_anon.png", height = 7, width = 15, plot = flow_rate_over_time_anon)

device_count_over_time_anon <- ggplot(data = visit_weekly) +
  geom_line(aes(x = date, y = average_devices/1000), color = "purple3") +
  labs(x = NULL, y = "Weekly Mean of Daily Device Count (Thousands)") +
  scale_x_date(breaks = "4 months", date_labels = "%b %d %Y") +
  scale_y_continuous(labels = scales::comma) +
  facet_wrap(~sewershed_code_ordered, scales = "free_y") +
  theme; device_count_over_time_anon
ggsave("Figures/device_count_over_time_anon.png", height = 7, width = 12, plot = device_count_over_time_anon)

# proportion of total observations that are non-detects
below_lod <- ww_impute_flow %>%
  mutate(below_lod_flag = pcr_target_below_lod == "yes") %>%
  group_by(utility, pcr_target) %>%
  summarise(n_obs = n(),
            pct_below_lod = mean(below_lod_flag, na.rm = TRUE)*100,
            .groups = "drop_last") %>%
  ungroup()

write.csv(below_lod, "DataProcessed/below_lod.csv")

plot_seasonal <- list()
for(i in 1:length(titles)){
  plot_seasonal[[i]] <- ggplot(data = ww_impute_flow %>% filter(pcr_target == pcr_targets[i],
                                                                date >= lp3_start,
                                                                date <= max(visit$date))) +
    geom_line(aes(x = date, y = pcr_target_avg_conc/1000, group = utility), color = "darkblue", alpha = 0.2) +
    labs(x = NULL, y = "Raw Concentration (1000s copies/L") +
    #ggtitle(titles[i]) +
    scale_x_date(limits = as.Date(c("2024-09-15", "2026-03-31")), date_labels = "%b %d %Y", date_breaks = "2 months") +
    scale_y_continuous(labels = scales::comma) +
    theme
  ggsave(paste0("Figures/plot_seasonal_", pcr_targets[i], ".png"), width = 8, height = 5, plot = plot_seasonal[[i]])
}

resp_season_start1 <- "2024-10-01"
resp_season_end1 <- "2025-05-15"
resp_season_start2 <- "2025-10-01"
resp_season_end2 <- "2026-03-31"

# read in Python WVAL output
wval <- imap_dfr(targets , ~ {
  read.csv(glue::glue("DataProcessed/WVAL/2026_04_21_dataset/2026_04_21_{.y}_wval_noagg.csv")) %>%
    mutate(pcr_target = .x,
           date = as.Date(sample_collect_date, "%Y-%m-%d")) %>%
    rename(utility = wwtp_name) %>%
    select(date, pcr_target, utility, WVAL)
})

# generate normalized wastewater concentration dataset
ww_with_wval <- ww_impute_flow %>%
  left_join(wval) %>%
  filter(date >= lp3_start,
         date <= max(visit$date))

# proportion of missing WVAL
percent(sum(is.na(ww_with_wval$WVAL))/nrow(ww_with_wval), accuracy = 0.01)

ww_with_wval_fill <- ww_with_wval %>%
  group_by(utility) %>%
  fill(WVAL, .direction = "down") %>%
  rename(WVAL_conc = WVAL) %>%
  mutate(WVAL_conc = signif(WVAL_conc, 3)/10)
  #select(utility, size, category_id, sewershed_id, sewershed_code, pcr_target, category, sample_type,
         #date, week_end, flow_rate_L, total_devices_daily, detect, detect_seasonal, census_population, contains("conc"))

ww_norm <- ww_with_wval_fill %>%
  left_join(visit_weekly) %>%
  mutate(raw_conc = case_when(pcr_target == "sars-cov-2" ~ pcr_target_avg_conc/10000,
                              TRUE ~ pcr_target_avg_conc/1000),
         avg_ww_use_L = 265,
         flow_rate_L = flow_rate*3.785*1000000,
         flow_norm_conc = signif(raw_conc*flow_rate_L/census_population, 3)/1000,
         mobile_norm_conc = signif(raw_conc/average_devices*avg_ww_use_L, 3)*100,
         combo_norm_conc = signif(raw_conc*(flow_rate_L/average_devices), 3)/100,
         season = case_when(date >= resp_season_start1 & date <= resp_season_end1 ~ "respiratory season",
                            date >= resp_season_start2 & date <= resp_season_end2 ~ "respiratory season",
                            TRUE ~ NA_character_),
         detect_seasonal = case_when(pcr_target %in% c("fluav", "rsv") & season == "respiratory season" ~ detect,
                                     pcr_target == "sars-cov-2" ~ detect,
                                     TRUE ~ NA_integer_)) %>%
           drop_na(detect_seasonal)

ww_norm_weekly <- ww_norm %>%
  group_by(utility, pcr_target, week_end) %>%
  mutate(raw_conc_weekly = mean(raw_conc),
         flow_norm_conc_weekly = mean(flow_norm_conc),
         mobile_norm_conc_weekly = mean(mobile_norm_conc),
         combo_norm_conc_weekly = mean(combo_norm_conc),
         WVAL_conc_weekly = mean(WVAL_conc),
         detect_weekly = as.numeric(any(detect_seasonal == 1))) %>%
  slice(1) %>%
  select(utility, size, category, contains("id"), sewershed_code, census_population, pcr_target,
         contains("week"))

# read in hospital admissions data
hosp_covid1 <- read.csv("DataRaw/2026/COVID_agg.csv") %>%
  mutate(date = as.Date(Event.Onset.Date, "%m/%d/%Y"),
         dataset = 1)

hosp_covid2 <- read.csv("DataRaw/CSPH Data _ April 2026/COVID_agg_April_update.csv") %>%
  mutate(date = as.Date(Event.Onset.Date, "%m/%d/%Y"),
         dataset = 2)

hosp_covid <- bind_rows(hosp_covid1, hosp_covid2) %>%
  mutate(Region.Name = case_when(Region.Level == "County" ~ paste(Region.Name, "County"),
                                 TRUE ~ Region.Name)) %>%
  select(-Case.Count) %>%
  pivot_wider(names_from = Region.Level, values_from = Hospitalized.Count) %>%
  mutate(pcr_target = "sars-cov-2") %>%
  rename(county_hosps = County,
         sewershed_hosps = Sewershed) %>%
  select(-Pathogen.Name, -Event.Onset.Date)

hosp_covid_county <- hosp_covid %>%
  filter(Region.Name %in% sentinel$county) %>%
  select(-sewershed_hosps) %>%
  rename(county = Region.Name)

hosp_covid_sewershed <- hosp_covid %>%
  filter(Region.Name %in% sentinel$utility[sentinel$size == "Small"]) %>%
  select(-county_hosps) %>%
  rename(utility = Region.Name)

hosp_covid_small <- sentinel_with_populations %>%
  filter(size == "Small") %>%
  left_join(hosp_covid_county) %>%
  left_join(hosp_covid_sewershed) %>%
  ungroup() %>%
  select(utility, county, size, pcr_target, date, county_hosps, sewershed_hosps)

hosp_covid_model <- hosp_covid %>%
  filter(Region.Name %in% sentinel$utility[sentinel$size == "Large"]) %>%
  mutate(size = "Large") %>%
  rename(utility = Region.Name) %>%
  select(utility, size, pcr_target, date, county_hosps, sewershed_hosps) %>%
  bind_rows(hosp_covid_small)

hosp_flu1 <- read.csv("DataRaw/2026/Flu_agg.csv") %>%
  filter(Pathogen.Name == "INFLUENZA A") %>%
  mutate(date = as.Date(Event.Onset.Date, "%m/%d/%Y"),
         dataset = 1)

hosp_flu2 <- read.csv("DataRaw/CSPH Data _ April 2026/Flu_agg_April_update.csv") %>%
  filter(Pathogen.Name == "INFLUENZA A") %>%
  mutate(date = as.Date(Event.Onset.Date, "%m/%d/%Y"),
         dataset = 2)

hosp_flu <- bind_rows(hosp_flu1, hosp_flu2) %>%
  mutate(Region.Name = case_when(Region.Level == "County" ~ paste(Region.Name, "County"),
                                 TRUE ~ Region.Name)) %>%
  pivot_wider(names_from = Region.Level, values_from = Hospitalized.Count) %>%
  mutate(pcr_target = "fluav") %>%
  rename(county_hosps = County,
         sewershed_hosps = Sewershed) %>%
  select(-Pathogen.Name, -Event.Onset.Date)

hosp_flu_county <- hosp_flu %>%
  filter(Region.Name %in% sentinel$county) %>%
  select(-sewershed_hosps) %>%
  rename(county = Region.Name)

hosp_flu_sewershed <- hosp_flu %>%
  filter(Region.Name %in% sentinel$utility[sentinel$size == "Small"]) %>%
  select(-county_hosps) %>%
  rename(utility = Region.Name)

hosp_flu_small <- sentinel_with_populations %>%
  filter(size == "Small") %>%
  left_join(hosp_flu_county) %>%
  left_join(hosp_flu_sewershed) %>%
  ungroup() %>%
  add_row(utility = "Telluride", county = "San Miguel County", size = "Small", pcr_target = "fluav",
          date = as.Date("2026-01-31"), county_hosps = 0, sewershed_hosps = 0) %>%
  select(utility, county, size, pcr_target, date, county_hosps, sewershed_hosps)

hosp_flu_model <- hosp_flu %>%
  filter(Region.Name %in% sentinel$utility[sentinel$size == "Large"]) %>%
  mutate(size = "Large") %>%
  rename(utility = Region.Name) %>%
  select(utility, size, pcr_target, date, county_hosps, sewershed_hosps) %>%
  bind_rows(hosp_flu_small)

hosp_rsv1 <- read.csv("DataRaw/2026/RSV_agg.csv") %>%
  group_by(Region.Level, Region.Name, Event.Onset.Date) %>%
  mutate(Hospitalized.Count = sum(Hospitalized.Count)) %>%
  slice(1) %>%
  select(-Age.Group) %>%
  mutate(date = as.Date(Event.Onset.Date, "%m/%d/%Y"),
         dataset = 1)

hosp_rsv2 <- read.csv("DataRaw/CSPH Data _ April 2026/RSV_agg_April_update.csv") %>%
  group_by(Region.Level, Region.Name, Event.Onset.Date) %>%
  mutate(Hospitalized.Count = sum(Hospitalized.Count)) %>%
  slice(1) %>%
  select(-Age.Group) %>%
  mutate(date = as.Date(Event.Onset.Date, "%m/%d/%Y"),
         dataset = 2)

hosp_rsv <- bind_rows(hosp_rsv1, hosp_rsv2) %>%
  mutate(Region.Name = case_when(Region.Level == "County" ~ paste(Region.Name, "County"),
                                 TRUE ~ Region.Name)) %>%
  pivot_wider(names_from = Region.Level, values_from = Hospitalized.Count) %>%
  mutate(pcr_target = "rsv") %>%
  rename(county_hosps = County,
         sewershed_hosps = Sewershed) %>%
  ungroup() %>%
  select(-Pathogen.Name, -Event.Onset.Date)

hosp_rsv_county <- hosp_rsv %>%
  filter(Region.Name %in% sentinel$county) %>%
  select(-sewershed_hosps) %>%
  rename(county = Region.Name)

hosp_rsv_sewershed <- hosp_rsv %>%
  filter(Region.Name %in% sentinel$utility[sentinel$size == "Small"]) %>%
  select(-county_hosps) %>%
  rename(utility = Region.Name)

hosp_rsv_small <- sentinel_with_populations %>%
  filter(size == "Small") %>%
  left_join(hosp_rsv_county) %>%
  left_join(hosp_rsv_sewershed) %>%
  ungroup() %>%
  select(utility, county, size, pcr_target, date, county_hosps, sewershed_hosps)

hosp_rsv_model <- hosp_rsv %>%
  filter(Region.Name %in% sentinel$utility[sentinel$size == "Large"]) %>%
  mutate(size = "Large") %>%
  rename(utility = Region.Name) %>%
  select(utility, size, pcr_target, date, county_hosps, sewershed_hosps) %>%
  bind_rows(hosp_rsv_small)

hosp_small <- bind_rows(hosp_covid_small, hosp_flu_small, hosp_rsv_small) %>%
  drop_na(pcr_target) %>%
  replace(is.na(.), 0) %>%
  filter(date >= lp3_start,
         date <= max(visit$date))

plot_hosp_small_sewershed <- list()
for(i in 1:length(titles)){
  plot_hosp_small_sewershed[[i]] <- ggplot(data = hosp_small %>% filter(pcr_target == pcr_targets[i])) +
    geom_bar(aes(x = date, y = sewershed_hosps), stat = "identity") +
    labs(x = NULL, y = "Daily Hospital Admissions") +
    ggtitle(paste(titles[i], "Hospital Admissions Over Time")) +
    scale_x_date(date_labels = "%b %Y") +
    facet_wrap(~utility, scales = "free_y") +
    theme
  ggsave(paste0("Figures/plot_hosp_small_sewershed_", pcr_targets[i], ".png"), width = 17, height = 9, plot = plot_hosp_small_sewershed[[i]])
}

plot_hosp_small_county <- list()
for(i in 1:length(titles)){
  plot_hosp_small_county[[i]] <- ggplot(data = hosp_small %>% filter(pcr_target == pcr_targets[i])) +
    geom_bar(aes(x = date, y = county_hosps), stat = "identity") +
    labs(x = NULL, y = "Daily Hospital Admissions") +
    ggtitle(paste(titles[i], "Hospital Admissions Over Time")) +
    scale_x_date(date_labels = "%b %Y") +
    facet_wrap(~county, scales = "free_y") +
    theme
  ggsave(paste0("Figures/plot_hosp_small_county_", pcr_targets[i], ".png"), width = 17, height = 9, plot = plot_hosp_small_county[[i]])
}

hosp_small_table <- hosp_small %>%
  group_by(utility, pcr_target) %>%
  summarize(total_sewershed_hosps = sum(sewershed_hosps),
            total_county_hosps = sum(county_hosps))

write.csv(hosp_small_table, "DataProcessed/hosp_small_table.csv")

hosp <- bind_rows(hosp_covid_model, hosp_flu_model, hosp_rsv_model) %>%
  drop_na(pcr_target) %>%
  mutate(admissions = ifelse(is.na(county), sewershed_hosps, county_hosps)) %>%
  select(utility, size, pcr_target, date, admissions) %>%
  filter(date >= lp3_start,
         date <= max(visit$date))

hosp_weekly <- hosp %>%
  mutate(week_end = ceiling_date(date, "weeks", week_start = getOption("lubridate.week.start", 6))) %>%
  group_by(utility, pcr_target, week_end) %>%
  mutate(admissions_weekly = sum(admissions)) %>%
  slice(1) %>%
  select(-date, -admissions)

ww_hosp_weekly <- merge(ww_norm_weekly, hosp_weekly, c("utility", "size", "pcr_target", "week_end"), all = TRUE) %>%
  # drop weeks with a value for hospital admissions but not wastewater
  drop_na(category) %>%
  mutate(admissions_weekly = replace_na(admissions_weekly, 0))

desc <- ww_hosp_weekly %>%
  mutate(sewershed_code_ordered = factor(sewershed_code, levels = sewershed_code_order),
         category_ordered = factor(category, levels = category_order))

norm_names <- c("Raw Concentration", "Flow-Normalized Concentration", "Mobile Device-Normalized Concentration",
                "Combination-Normalized Concentration", "WVAL")

desc_long <- desc %>%
  pivot_longer(cols = c(raw_conc_weekly, flow_norm_conc_weekly, mobile_norm_conc_weekly, combo_norm_conc_weekly, WVAL_conc_weekly),
               names_to = "norm_method", values_to = "value")

norm_methods <- desc_long$norm_method[1:5]

facet_names_pathogens <- c("fluav" = "Influenza A", "rsv" = "RSV", "sars-cov-2" = "SARS-CoV-2")

boxplots_anon <- list()
for(i in 1:length(norm_methods)){
  boxplots_anon[[i]] <- ggplot(data = desc_long %>% filter(norm_method == norm_methods[i])) +
    geom_boxplot(aes(x = log10(value), y = sewershed_code_ordered, group = sewershed_code_ordered, fill = category_ordered)) +
    labs(x = paste0("log10(", norm_names[i], ")"), y = NULL, fill = "Category") +
    scale_y_discrete(limits = rev) +
    #guides(fill = guide_legend(reverse = TRUE)) +
    facet_grid(~pcr_target, scales = "free", labeller = as_labeller(facet_names_pathogens)) +
    theme +
    theme(axis.text.x = element_text(angle = 0))
  ggsave(paste0("Figures/boxplots_anon_", norm_methods[i], ".png"), height = 6, width = 11, plot = boxplots_anon[[i]])
}

norm_order <- c("raw_conc", "flow_norm_conc", "mobile_norm_conc", "combo_norm_conc", "WVAL_conc")

summary_table <- desc %>%
  group_by(utility, pcr_target) %>%
  mutate(detect_freq = percent(mean(detect_weekly, na.rm = TRUE), accuracy = 0.1),
         hosp_incidence = signif(sum(admissions_weekly)/census_population*100000, 3)) %>%
  slice(1) %>%
  ungroup() %>%
  select(sewershed_code_ordered, category, pcr_target, detect_freq, hosp_incidence) %>%
  pivot_wider(names_from = pcr_target, values_from = c(detect_freq, hosp_incidence)) %>%
  select(sewershed_code_ordered, category, contains("fluav"), contains("rsv"), contains("sars-cov-2")) %>%
  arrange(sewershed_code_ordered)

write.csv(summary_table, "DataProcessed/summary_table.csv")

# set up dataset for cross_correlation
cross_corr <- ww_hosp_weekly %>%
  filter(utility %in% c("Metro Wastewater RWHTF - PRC", "Aurora", "CO Springs - Las Vegas")) %>%
  mutate(total_pop = sum(unique(census_population))) %>%
  select(week_end, pcr_target, raw_conc_weekly, admissions_weekly, total_pop) %>%
  group_by(pcr_target, week_end) %>%
  mutate(ww = mean(raw_conc_weekly),
         hosp = sum(admissions_weekly)) %>%
  slice(1) %>%
  mutate(hosp_per_100K = hosp/total_pop*100000) %>%
  select(-hosp) %>%
  pivot_wider(id_cols = week_end, names_from = pcr_target,
              values_from = c(ww, hosp_per_100K)) %>%
  #replace(is.na(.), 0)
  drop_na()

find_best_lag <- function(ww, hosp, max_lag = 4){
  
  df <- tibble(ww = ww, hosp = hosp) %>%
    filter(!is.na(ww), !is.na(hosp))
  
  if(nrow(df) == 0){
    stop("No complete observations remain after filtering.")
  }
  
  ccf_result <- ccf(x = df$ww, y = df$hosp,
                    lag.max = max_lag, plot = FALSE)
  
  tibble(lag = as.numeric(ccf_result$lag),
         correlation = as.numeric(ccf_result$acf)) %>%
    filter(lag >= 0) %>%
    arrange(desc(abs(correlation)))
}

lag_results <- bind_rows(fluav = find_best_lag(cross_corr$ww_fluav, cross_corr$hosp_per_100K_fluav),
                     rsv = find_best_lag(cross_corr$ww_rsv, cross_corr$hosp_per_100K_rsv),
                     `sars-cov-2` = find_best_lag(cross_corr$`ww_sars-cov-2`, cross_corr$`hosp_per_100K_sars-cov-2`),
                     .id = "pcr_target")

write.csv(lag_results, "DataProcessed/lag_results.csv")

# test the data in a negative binomial model to see if log-transforming improves the
# residual structure around the fitted values
library(MASS)
m_orig <- glm.nb(admissions_weekly ~ raw_conc_weekly + detect_weekly,
                 data = ww_hosp_weekly)

m_log <- glm.nb(admissions_weekly ~ log10(raw_conc_weekly) + detect_weekly,
                data = ww_hosp_weekly)

plot_df <- bind_rows(tibble(fitted = fitted(m_orig),
                            residual = residuals(m_orig, type = "pearson"),
                            model = "Original"),
                     tibble(fitted = fitted(m_log),
                            residual = residuals(m_log, type = "pearson"),
                            model = "Log-Transformed"))

residuals <- ggplot(plot_df, aes(fitted, residual)) +
  geom_point(alpha = 0.3) +
  geom_smooth(method = "loess",
              se = FALSE) +
  facet_wrap(~fct_rev(model)) +
  geom_hline(yintercept = 0, linetype = 2); residuals

sim_orig <- simulateResiduals(m_orig)
sim_log <- simulateResiduals(m_log)

testDispersion(sim_orig)
testDispersion(sim_log)

testUniformity(sim_orig)
testUniformity(sim_log)

AIC(m_orig, m_log)

auto <- glm.nb(admissions_weekly ~ raw_conc_weekly + detect_weekly,
               data = ww_hosp_weekly %>% filter(utility == "Metro Wastewater RWHTF - PRC",
                                                pcr_target == "sars-cov-2"))
sim <- simulateResiduals(auto)
plot(sim)
acf(residuals(auto)) # some autocorrelation at 1-5 weeks, include lagged admissions as predictor

detach("package:MASS")

ww_hosp_model <- ww_hosp_weekly %>%
  mutate(utility = as.factor(utility)) %>%
  group_by(utility, pcr_target) %>%
  mutate(admissions_weekly_w_1 = lag(admissions_weekly, 1)) %>%
  mutate(utility_ordered = factor(utility, levels = utility_order),
         sewershed_code_ordered = factor(sewershed_code, levels = sewershed_code_order),
         category_ordered = factor(category, levels = category_order))

plot_hosp_ts_anon <- list()
for(i in 1:length(titles)){
  plot_hosp_ts_anon[[i]] <- ggplot(data = ww_hosp_model %>% filter(pcr_target == pcr_targets[i])) +
    geom_bar(aes(x = week_end, y = admissions_weekly/census_population*100000, fill = size), stat = "identity") +
    labs(x = NULL, y = "Weekly Total Hospital Admissions per 100K", fill = "Size") +
    scale_x_date(limits = c(as.Date("2024-07-15"), as.Date("2026-04-01")),
                 date_breaks = "4 months",
                 date_labels = "%b %d %Y") +
    scale_fill_manual(labels = c("Large", "Small"),
                      values = c("darkblue", "violetred")) +
    facet_wrap(~sewershed_code_ordered) + #scales = "free_y") +
    theme
  ggsave(paste0("Figures/plot_hosp_ts_anon_", pcr_targets[i], ".png"), height = 7, width = 13, plot = plot_hosp_ts_anon[[i]])
}

ww_hosp_plot_example_large <- ggplot(data = ww_hosp_model %>% filter(sewershed_code == "OL3")) +
  geom_bar(aes(x = week_end, y = admissions_weekly*10), stat = "identity", alpha = 0.6) +
  geom_line(aes(x = week_end, y = raw_conc_weekly, color = "color1")) +
  geom_line(aes(x = week_end, y = flow_norm_conc_weekly, color = "color2")) +
  geom_line(aes(x = week_end, y = mobile_norm_conc_weekly, color = "color3")) +
  geom_line(aes(x = week_end, y = combo_norm_conc_weekly, color = "color4")) +
  geom_line(aes(x = week_end, y = WVAL_conc_weekly, color = "color5")) +
  labs(x = "Week End Date", y = "Viral Concentration (copies/L/person)", color = "Normalization Method") +
  scale_x_date(date_breaks = "2 months", date_labels = "%b %d %Y") +
  scale_y_continuous(sec.axis = sec_axis(~. * 0.1, name = "Weekly Total Hospital Admissions")) +
  scale_color_manual(labels = c("color1" = "Raw",
                                "color2" = "Flow-Normalized",
                                "color3" = "Mobile Device-Normalized",
                                "color4" = "Combination-Normalized",
                                "color5" = "WVAL"),
                     values = c("color1" = "darkorchid4",
                                "color2" = "turquoise2",
                                "color3" = "firebrick",
                                "color4" = "violet",
                                "color5" = "orange")) +
  guides(color = guide_legend(override.aes = list(linewidth = 8), nrow = 5)) +
  facet_wrap(~pcr_target, nrow = 3, scales = "free_y", labeller = as_labeller(facet_names_pathogens)) +
  theme; ww_hosp_plot_example_large
ggsave("Figures/ww_hosp_plot_example_large.png", width = 9, height = 6, plot = ww_hosp_plot_example_large)

ww_hosp_plot_example_small <- ggplot(data = ww_hosp_model %>% filter(sewershed_code == "HT2")) +
  geom_bar(aes(x = week_end, y = admissions_weekly*10), stat = "identity", alpha = 0.6) +
  geom_line(aes(x = week_end, y = raw_conc_weekly, color = "color1")) +
  geom_line(aes(x = week_end, y = flow_norm_conc_weekly, color = "color2")) +
  geom_line(aes(x = week_end, y = mobile_norm_conc_weekly, color = "color3")) +
  geom_line(aes(x = week_end, y = combo_norm_conc_weekly, color = "color4")) +
  geom_line(aes(x = week_end, y = WVAL_conc_weekly, color = "color5")) +
  labs(x = "Week End Date", y = "Viral Concentration (copies/L/person)", color = "Normalization Method") +
  scale_x_date(date_breaks = "2 months", date_labels = "%b %d %Y") +
  scale_y_continuous(sec.axis = sec_axis(~. * 0.1, name = "Weekly Total Hospital Admissions per 100K")) +
  scale_color_manual(labels = c("color1" = "Raw",
                                "color2" = "Flow-Normalized",
                                "color3" = "Mobile Device-Normalized",
                                "color4" = "Combination-Normalized",
                                "color5" = "WVAL"),
                     values = c("color1" = "darkorchid4",
                                "color2" = "turquoise2",
                                "color3" = "firebrick",
                                "color4" = "violet",
                                "color5" = "orange")) +
  guides(color = guide_legend(override.aes = list(linewidth = 8), nrow = 5)) +
  facet_wrap(~pcr_target, nrow = 3, scales = "free_y", labeller = as_labeller(facet_names_pathogens)) +
  theme; ww_hosp_plot_example_small
ggsave("Figures/ww_hosp_plot_example_small.png", width = 9, height = 6, plot = ww_hosp_plot_example_small)


#ww_hosp_plot <- list()
# for(i in 1:length(titles)){
#   ww_hosp_plot[[i]] <- ggplot(data = ww_hosp_model %>% filter(pcr_target == pcr_targets[i])) +
#     geom_bar(aes(x = week_end, y = admissions_weekly/census_population*100000), stat = "identity") +
#     geom_line(aes(x = week_end, y = raw_conc_weekly, color = "color1")) +
#     geom_line(aes(x = week_end, y = flow_norm_conc_weekly, color = "color2")) +
#     geom_line(aes(x = week_end, y = mobile_norm_conc_weekly, color = "color3")) +
#     geom_line(aes(x = week_end, y = combo_norm_conc_weekly, color = "color4")) +
#     geom_line(aes(x = week_end, y = WVAL_conc_weekly, color = "color5")) +
#     labs(x = "Week End Date", y = "Value") +
#     scale_color_manual(labels = c("color1" = "Raw",
#                                   "color2" = "Flow-Normalized",
#                                   "color3" = "Mobile Device-Normalized",
#                                   "color4" = "Combination-Normalized",
#                                   "color5" = "WVAL"),
#                        values = c("color1" = "darkorchid4",
#                                   "color2" = "turquoise2",
#                                   "color3" = "firebrick",
#                                   "color4" = "violet",
#                                   "color5" = "orange")) +
#     guides(color = guide_legend(override.aes = list(linewidth = 8), nrow = 5)) +
#     facet_wrap(~utility, scales = "free_y") +
#     theme
#   ggsave(paste0("Figures/ww_hosp_plot_", pcr_targets[i], ".png"), width = 17, height = 9, plot = ww_hosp_plot[[i]])
# }


ww_hosp_model_exclude <- ww_hosp_model %>%
  group_by(utility, pcr_target) %>%
  filter(sum(admissions_weekly) > 5)

# variables to loop over
norm_vars <- c("raw_conc_weekly", "flow_norm_conc_weekly", "mobile_norm_conc_weekly",
               "combo_norm_conc_weekly", "WVAL_conc_weekly")

norm_names <- c("raw", "flow", "mobile", "combo", "wval")

utilities <- unique(ww_hosp_model$utility_ordered)

# create a grid to organize pathogen x sewershed combos
model_grid <- expand_grid(sewershed = utilities,
                          target = pcr_targets,
                          var = norm_vars) %>%
  mutate(row_id = row_number())

library(MASS)
# function to build NB model formula for any normalization method
build_formula <- function(var) {
  as.formula(paste("admissions_weekly ~ poly(", var, ", 2, raw = TRUE) + admissions_weekly_w_1"))
}

fit_nb_model <- function(var, target, sewershed, row_id) {
  
  dat <- ww_hosp_model_exclude %>%
    filter(utility == sewershed,
           pcr_target == target)
  
  if (nrow(dat) == 0) {
    return(list(model = NULL,
                warning = NA_character_,
                error = "no_data",
                row_id = row_id))
  }
  
  warning_msg <- NULL
  error_msg <- NULL
  
  model <- withCallingHandlers(tryCatch(
    glm.nb(formula = build_formula(var),
           data = dat),
    error = function(e) {
      error_msg <<- conditionMessage(e)
      return(NULL)
      }),
    warning = function(w) {
      warning_msg <<- paste(warning_msg, conditionMessage(w), sep = "; ")
      invokeRestart("muffleWarning")
    })
  
  list(model = model,
       warning = warning_msg,
       error = error_msg,
       row_id = row_id)
}

results <- model_grid %>%
  mutate(fit = pmap(list(var, target, sewershed, row_id), fit_nb_model))

results_tidy <- results %>%
  mutate(model_obj = map(fit, "model"),
    #warning = map_chr(fit, "warning"),
    warning = map_chr(fit, ~ {
      w <- .x$warning
      
      if (is.null(w) || length(w) == 0) {
        return(NA_character_)
      }
      w
    }),
    #error = map_chr(fit, "error"),
    error = map_chr(fit, ~ {
      e <- .x$error
      
      if (is.null(e) || length(e) == 0) {
        return(NA_character_)
      }
      e
    }),
    converged = map_lgl(model_obj, ~ {
      if (is.null(.x)) return(FALSE)
      isTRUE(.x$converged)
    }),
    AIC = map_dbl(model_obj, ~ {
      if (is.null(.x)) return(NA_real_)
      AIC(.x)
    })
  )

calc_mae <- function(model) {
  
  if (is.null(model)) return(NA_real_)
  
  y <- model$model$admissions_weekly
  yhat <- predict(model, type = "response")
  
  mean(abs(y - yhat), na.rm = TRUE)
}

detach("package:MASS")

var_order <- c("Raw", "Flow", "Mobile", "Combo", "WVAL")

results_full <- results_tidy %>%
  mutate(var = factor(recode(var, "raw_conc_weekly" = "Raw",
                                "flow_norm_conc_weekly" = "Flow",
                                "mobile_norm_conc_weekly" = "Mobile",
                                "combo_norm_conc_weekly" = "Combo",
                                "WVAL_conc_weekly" = "WVAL"), levels = var_order),
         model = "quadratic",
         MAE = map_dbl(model_obj, calc_mae)) %>%
  select(sewershed, target, var, model, row_id, warning, error, converged, AIC, MAE)

write.csv(results_full, "DataProcessed/results_full.csv")

results_table <- results_full %>%
  group_by(sewershed, target) %>%
  arrange(AIC) %>%
  mutate(delta_AIC_cumul = AIC - min(AIC),
         delta_AIC_overall = format(round(max(AIC) - min(AIC), 2), nsmall = 2),
         delta_AIC_overall = case_when(delta_AIC_overall == "NA" ~ "",
                                       TRUE ~ delta_AIC_overall),
         rank = case_when(is.na(AIC) ~ NA_character_,
                          TRUE ~ as.factor(row_number())),
         distance = as.factor(case_when(delta_AIC_cumul < 2 ~ "<2",
                                        delta_AIC_cumul >= 2 & delta_AIC_cumul < 4 ~ ">=2 and <4",
                                        delta_AIC_cumul >= 4 & delta_AIC_cumul < 6 ~ ">=4 and <6",
                                        delta_AIC_cumul >= 6 & delta_AIC_cumul < 8 ~ ">=6 and <8",
                                        delta_AIC_cumul >= 8 ~ ">=8",
                                        TRUE ~ NA_character_))) %>%
  rename(utility = sewershed) %>%
  left_join(sentinel_with_populations_random) %>%
  mutate(utility_ordered = factor(utility, levels = utility_order),
         sewershed_code_ordered = factor(sewershed_code, levels = sewershed_code_order),
         category_ordered = factor(category, levels = category_order),
         MAE_100K = round(MAE/census_population*100000, 2))%>%
  select(utility, target, var, warning, error, converged, contains("AIC"), contains("MAE"), rank, distance,
         contains("ordered"))

install.packages("ggh4x")
library(ggh4x)
heatmap_rank_anon <- ggplot(data = results_table) +
  geom_tile(aes(x = var, y = sewershed_code_ordered, fill = distance), color = "white", linewidth = 0.1) +
  geom_text(aes(x = var, y = sewershed_code_ordered, label = rank, color = distance)) +
  geom_text(data = results_table %>% distinct(sewershed_code_ordered, delta_AIC_overall),
            aes(x = 5.5, y = sewershed_code_ordered, label = delta_AIC_overall, nudge_x = 0.35)) +
  labs(x = NULL, y = NULL, fill = "Difference in AIC from Leading Model") +
  scale_x_discrete(expand = expansion(mult = c(0, 0.3))) +
  scale_y_discrete(limits = rev) +
  scale_color_manual(values = c("<2" = "white",
                                ">=2 and <4" = "white",
                                ">=4 and <6" = "white",
                                ">=6 and <8" = "black",
                                ">=8" = "black")) +
  scale_fill_manual(values = c("<2" = "midnightblue",
                               ">=2 and <4" = "dodgerblue4",
                               ">=4 and <6" = "steelblue3",
                               ">=6 and <8" = "skyblue2",
                               ">=8" = "lightblue1"),
                    breaks = c("<2", ">=2 and <4", ">=4 and <6", ">=6 and <8", ">=8"),
                    na.value = "grey30") +
  guides(color = "none", fill = guide_legend(title.position = "top", title.hjust = 0.5)) +
  ggh4x::facet_grid2(~target, scales = "free_y", independent = "y", labeller = as_labeller(facet_names_pathogens)) +
  theme +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5),
        legend.position = "bottom",
        legend.key.spacing.x = unit(1, "cm")); heatmap_rank_anon

ggsave("Figures/heatmap_rank_anon.png", width = 15, height = 7, plot = heatmap_rank_anon)

MAE_summary <- results_table %>%
  group_by(utility, target) %>%
  arrange(target, match(sewershed_code_ordered, sewershed_code_order), match(var, var_order)) %>%
  ungroup() %>%
  mutate(Pathogen = recode(target, "fluav" = "Influenza A",
                         "rsv" = "RSV",
                         "sars-cov-2" = "SARS-CoV-2")) %>%
  rename(Category = category_ordered,
         `Sewershed Code` = sewershed_code_ordered) %>%
  select(Pathogen, Category, `Sewershed Code`, var, MAE_100K) %>%
  pivot_wider(names_from = var, values_from = MAE_100K)

write.csv(MAE_summary, "DataProcessed/MAE_summary.csv")

MAE_summary_pathogen <- MAE_summary %>%
  group_by(Pathogen) %>%
  summarize(across(c(Raw, Flow, Mobile, Combination, WVAL), ~ round(mean(.x, na.rm = TRUE), 2)))

write.csv(MAE_summary_pathogen, "DataProcessed/MAE_summary_pathogen.csv")

MAE_summary_category <- MAE_summary %>%
  group_by(Category) %>%
  summarize(across(c(Raw, Flow, Mobile, Combination, WVAL), ~ round(mean(.x, na.rm = TRUE), 2)))

write.csv(MAE_summary_category, "DataProcessed/MAE_summary_category.csv")

MAE_summary_overall <- MAE_summary %>%
  group_by(Pathogen) %>%
  summarize(across(c(Raw, Flow, Mobile, Combo, WVAL), ~ round(mean(.x, na.rm = TRUE), 2))) %>%
  mutate(Pathogen = recode(Pathogen, "Influenza A" = "Influenza A Overall",
                         "RSV" = "RSV Overall",
                         "SARS-CoV-2" = "SARS-CoV-2 Overall"),
         Category = "All") %>%
  bind_rows(MAE_summary) %>%
  arrange(Pathogen, `Sewershed Code`) %>%
  rowwise() %>%
  mutate(`Smallest MAE` = {
    vals <- c_across(Raw:WVAL)
    method_names <- c("Raw", "Flow", "Mobile", "Combo", "WVAL")
    
    if (all(is.na(vals))) {
      NA_character_
    } else {
      paste(
        method_names[vals == min(vals, na.rm = TRUE)],
        collapse = ", ")
    }
  }) %>%
  select(Pathogen, Category, `Sewershed Code`, Raw, Flow, Mobile, Combo, WVAL, `Smallest MAE`)

write.csv(MAE_summary_overall, "DataProcessed/MAE_summary_overall.csv")

# sensitivity analysis with linear
library(MASS)
# function to build NB model formula for any normalization method
build_formula_linear <- function(var) {
  as.formula(paste("admissions_weekly ~", var, "+ detect_weekly + admissions_weekly_w_1"))
}

fit_nb_model_linear <- function(var, target, sewershed, row_id) {
  
  dat <- ww_hosp_model_exclude %>%
    filter(utility == sewershed,
           pcr_target == target)
  
  if (nrow(dat) == 0) {
    return(list(model = NULL,
                warning = NA_character_,
                error = "no_data",
                row_id = row_id))
  }
  
  warning_msg <- NULL
  error_msg <- NULL
  
  model <- withCallingHandlers(tryCatch(
    glm.nb(formula = build_formula_linear(var),
           data = dat),
    error = function(e) {
      error_msg <<- conditionMessage(e)
      return(NULL)
    }),
    warning = function(w) {
      warning_msg <<- paste(warning_msg, conditionMessage(w), sep = "; ")
      invokeRestart("muffleWarning")
    })
  
  list(model = model,
       warning = warning_msg,
       error = error_msg,
       row_id = row_id)
}

results_linear <- model_grid %>%
  mutate(fit = pmap(list(var, target, sewershed, row_id), fit_nb_model_linear))

results_tidy_linear <- results_linear %>%
  mutate(model_obj = map(fit, "model"),
         #warning = map_chr(fit, "warning"),
         warning = map_chr(fit, ~ {
           w <- .x$warning
           
           if (is.null(w) || length(w) == 0) {
             return(NA_character_)
           }
           w
         }),
         #error = map_chr(fit, "error"),
         error = map_chr(fit, ~ {
           e <- .x$error
           
           if (is.null(e) || length(e) == 0) {
             return(NA_character_)
           }
           e
         }),
         converged = map_lgl(model_obj, ~ {
           if (is.null(.x)) return(FALSE)
           isTRUE(.x$converged)
         }),
         AIC = map_dbl(model_obj, ~ {
           if (is.null(.x)) return(NA_real_)
           AIC(.x)
         })
  )

detach("package:MASS")

results_full_linear <- results_tidy_linear %>%
  mutate(var = factor(recode(var, "raw_conc_weekly" = "Raw",
                             "flow_norm_conc_weekly" = "Flow",
                             "mobile_norm_conc_weekly" = "Mobile",
                             "combo_norm_conc_weekly" = "Combination",
                             "WVAL_conc_weekly" = "WVAL"), levels = var_order),
         model = "linear",
         MAE = map_dbl(model_obj, calc_mae)) %>%
  select(sewershed, target, var, model, row_id, warning, error, converged, AIC, MAE)

write.csv(results_full_linear, "DataProcessed/results_full_linear.csv")

results_table_linear <- results_full_linear %>%
  group_by(sewershed, target) %>%
  arrange(AIC) %>%
  mutate(delta_AIC_cumul = AIC - min(AIC),
         delta_AIC_overall = format(round(max(AIC) - min(AIC), 2), nsmall = 2),
         delta_AIC_overall = case_when(delta_AIC_overall == "NA" ~ "",
                                       TRUE ~ delta_AIC_overall),
         rank = case_when(is.na(AIC) ~ NA_character_,
                          TRUE ~ as.factor(row_number())),
         distance = as.factor(case_when(delta_AIC_cumul < 2 ~ "<2",
                                        delta_AIC_cumul >= 2 & delta_AIC_cumul < 4 ~ ">=2 and <4",
                                        delta_AIC_cumul >= 4 & delta_AIC_cumul < 6 ~ ">=4 and <6",
                                        delta_AIC_cumul >= 6 & delta_AIC_cumul < 8 ~ ">=6 and <8",
                                        delta_AIC_cumul >= 8 ~ ">=8",
                                        TRUE ~ NA_character_))) %>%
  rename(utility = sewershed) %>%
  left_join(sentinel_with_populations_random) %>%
  mutate(utility_ordered = factor(utility, levels = utility_order),
         sewershed_code_ordered = factor(sewershed_code, levels = sewershed_code_order),
         category_ordered = factor(category, levels = category_order),
         MAE_100K = round(MAE/census_population*100000, 2))%>%
  select(utility, target, var, warning, error, converged, contains("AIC"), contains("MAE"), rank, distance,
         contains("ordered"))

heatmap_rank_linear_anon <- ggplot(data = results_table_linear) +
  geom_tile(aes(x = var, y = sewershed_code_ordered, fill = distance), color = "white", linewidth = 0.1) +
  geom_text(aes(x = var, y = sewershed_code_ordered, label = rank, color = distance)) +
  geom_text(data = results_table_linear %>% distinct(sewershed_code_ordered, delta_AIC_overall),
            aes(x = 5.5, y = sewershed_code_ordered, label = delta_AIC_overall, nudge_x = 0.35)) +
  labs(x = NULL, y = NULL, fill = "Difference in AIC from Leading Model") +
  scale_x_discrete(expand = expansion(mult = c(0, 0.3))) +
  scale_y_discrete(limits = rev) +
  scale_color_manual(values = c("<2" = "white",
                                ">=2 and <4" = "white",
                                ">=4 and <6" = "white",
                                ">=6 and <8" = "black",
                                ">=8" = "black")) +
  scale_fill_manual(values = c("<2" = "midnightblue",
                               ">=2 and <4" = "dodgerblue4",
                               ">=4 and <6" = "steelblue3",
                               ">=6 and <8" = "skyblue2",
                               ">=8" = "lightblue1"),
                    breaks = c("<2", ">=2 and <4", ">=4 and <6", ">=6 and <8", ">=8"),
                    na.value = "grey30") +
  guides(color = "none", fill = guide_legend(title.position = "top", title.hjust = 0.5)) +
  ggh4x::facet_grid2(~target, scales = "free_y", independent = "y", labeller = as_labeller(facet_names_pathogens)) +
  theme +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5),
        legend.position = "bottom",
        legend.key.spacing.x = unit(1, "cm")); heatmap_rank_linear_anon

ggsave("Figures/heatmap_rank_linear_anon.png", width = 15, height = 7, plot = heatmap_rank_linear_anon)

MAE_summary_linear <- results_table_linear %>%
  group_by(utility, target) %>%
  arrange(target, match(sewershed_code_ordered, sewershed_code_order), match(var, var_order)) %>%
  ungroup() %>%
  mutate(Pathogen = recode(target, "fluav" = "Influenza A",
                           "rsv" = "RSV",
                           "sars-cov-2" = "SARS-CoV-2")) %>%
  rename(Category = category_ordered,
         `Sewershed Code` = sewershed_code_ordered) %>%
  select(Pathogen, Category, `Sewershed Code`, var, MAE_100K) %>%
  pivot_wider(names_from = var, values_from = MAE_100K)

write.csv(MAE_summary_linear, "DataProcessed/MAE_summary_linear.csv")

MAE_summary_pathogen_linear <- MAE_summary_linear %>%
  group_by(Pathogen) %>%
  summarize(across(c(Raw, Flow, Mobile, Combination, WVAL), ~ round(mean(.x, na.rm = TRUE), 2)))

write.csv(MAE_summary_pathogen_linear, "DataProcessed/MAE_summary_pathogen_linear.csv")

MAE_summary_category_linear <- MAE_summary_linear %>%
  group_by(Category) %>%
  summarize(across(c(Raw, Flow, Mobile, Combination, WVAL), ~ round(mean(.x, na.rm = TRUE), 2)))

write.csv(MAE_summary_category_linear, "DataProcessed/MAE_summary_category_linear.csv")


MAE_summary_linear_overall <- MAE_summary_linear %>%
  group_by(Pathogen) %>%
  summarize(across(c(Raw, Flow, Mobile, Combo, WVAL), ~ round(mean(.x, na.rm = TRUE), 2))) %>%
  mutate(Pathogen = recode(Pathogen, "Influenza A" = "Influenza A Overall",
                           "RSV" = "RSV Overall",
                           "SARS-CoV-2" = "SARS-CoV-2 Overall"),
         Category = "All") %>%
  bind_rows(MAE_summary_linear) %>%
  arrange(Pathogen, `Sewershed Code`) %>%
  rowwise() %>%
  mutate(`Smallest MAE` = {
    vals <- c_across(Raw:WVAL)
    method_names <- c("Raw", "Flow", "Mobile", "Combo", "WVAL")
    
    if (all(is.na(vals))) {
      NA_character_
    } else {
      paste(
        method_names[vals == min(vals, na.rm = TRUE)],
        collapse = ", ")
    }
  }) %>%
  select(Pathogen, Category, `Sewershed Code`, Raw, Flow, Mobile, Combo, WVAL, `Smallest MAE`)

write.csv(MAE_summary_linear_overall, "DataProcessed/MAE_summary_linear_overall.csv")

sewershed_names_to_codes <- ww_impute_flow %>%
  ungroup() %>%
  select(utility, sewershed_code_ordered) %>%
  group_by(utility) %>%
  slice(1)

write.csv(sewershed_names_to_codes, "DataProcessed/sewershed_names_to_codes.csv")

################################### GARBAGE BELOW ########################################################







# fit_nb_model <- function(var, target, sewershed) {
#   
#   dat <- ww_hosp_model_exclude %>%
#     filter(utility == sewershed,
#            pcr_target == target)
#   
#   if (nrow(dat) == 0) {
#     return(list(
#       model = NULL,
#       warning = NA_character_,
#       error = NA_character_))
#   }
#   
#   glm.nb(formula = build_formula(var),
#          data = dat)
# }



# safe_fit_nb <- purrr::possibly(
#   fit_nb_model,
#   otherwise = NULL
# )

# safe_AIC <- function(model) {
#   if (is.null(model)) return(NA_real_)
#   AIC(model)
# }
# 
# results <- model_grid %>%
#   mutate(
#     model = pmap(
#       list(var, target, sewershed),
#       fit_nb_model
#     )
#   )
# 
# results <- results %>% mutate(converged = map_lgl(model, ~ {
#   if (is.null(.x)) return(FALSE)
#   isTRUE(.x$converged)
#   }))
# 
# 
# 
# results <- results %>%
#   mutate(AIC = map_dbl(model, safe_AIC))
# 
# 
# 
# 
# results <- model_grid %>%
#   mutate(model = pmap(list(var, target, sewershed), fit_nb_model))
# 
# results %>%
#   mutate(status = map_chr(model, "status"))
#          
         
         
         





# function to fit NB model for any pathogen and sewershed
models <- function(df, conc_var) {
  formula <- as.formula(paste("admissions_weekly ~ poly(", var, ", 2, raw = TRUE) + admissions_weekly_w_1"))
  fit_nb_model <- glm.nb(formula, data = df)
}






fit_nb_model <- function(var, target, sewershed) {
  glmmTMB(formula = build_formula(var),
          family = nbinom2,
          data = ww_hosp_model_exclude %>% filter(utility == sewershed, pcr_target == target),
          control = glmmTMBControl(optimizer = optim, optArgs = list(method = "BFGS")))
}





linear_vs_quadratic <- function(df, conc_var) {
  
  #linear_formula <- as.formula(paste("admissions_weekly ~", conc_var, "+ detect_weekly + admissions_weekly_w_1"))
  quadratic_formula <- as.formula(paste("admissions_weekly ~ poly(", conc_var, ", 2, raw = TRUE)", "+ admissions_weekly_w_1"))
  
  #linear <- glm.nb(linear_formula,
                   #data = df)
  
  quadratic <- glm.nb(quadratic_formula,
                      data = df)
  
  # linear <- glmmTMB(linear_formula,
  #                   family = nbinom2,
  #                   data = df,
  #                   control = glmmTMBControl(optimizer = optim, optArgs = list(method = "BFGS")))
  # 
  # quadratic <- glmmTMB(quadratic_formula,
  #                      family = nbinom2,
  #                      data = df,
  #                      control = glmmTMBControl(optimizer = optim, optArgs = list(method = "BFGS")))
  
  list(linear = linear,
       quadratic = quadratic)
  
}

linear_vs_quadratic_results <- map_dfr(norm_vars, function(norm_var) {
  
  ww_hosp_model_exclude %>%
    group_by(utility, pcr_target) %>%
    group_nest() %>%
    mutate(
      models = map(data, linear_vs_quadratic, conc_var = norm_var),
      aic = map(models, ~ tibble(
        model = c("linear", "quadratic"),
        AIC = c(
          AIC(.x$linear),
          AIC(.x$quadratic)
        )
      ))
    ) %>%
    #select(utility, pcr_target, aic) %>%
    unnest(aic) %>%
    mutate(normalization = norm_var)
  
})

detach("package:MASS")















# test for multicollinearity
mc <- ww_hosp_model_exclude %>%
  group_by(utility, pcr_target) %>%
  summarize(correlation = cor(raw_conc_weekly, detect_weekly,
                              use = "complete.obs"),
            .groups = "drop")

mc_vif <- ww_hosp_model_exclude %>%
  group_by(utility, pcr_target, .add = TRUE) %>%
  group_split() %>%
  map_dfr(function(sub_df) {
    
    model <- lm(admissions_weekly ~ raw_conc_weekly + detect_weekly, data = sub_df)
    
    vif_vals <- vif(model)
    
    tibble(utility = unique(sub_df$utility),
           pcr_target  = unique(sub_df$pcr_target),
           vif = vif_vals["raw_conc_weekly"]
    )
  })

mc_combined <- left_join(mc, mc_vif)

write.csv(mc_combined, "DataProcessed/mc_combined.csv")

## with quadratic term
norm_vars <- c("raw_conc_weekly", "flow_norm_conc_weekly", "mobile_norm_conc_weekly",
               "combo_norm_conc_weekly", "WVAL_conc_weekly")

library(MASS)

linear_vs_quadratic <- function(df, conc_var) {
  
  linear_formula <- as.formula(paste("admissions_weekly ~", conc_var, "+ detect_weekly + admissions_weekly_w_1"))
  quadratic_formula <- as.formula(paste("admissions_weekly ~ poly(", conc_var, ", 2, raw = TRUE)", "+ admissions_weekly_w_1"))
  
  linear <- glm.nb(linear_formula,
                   data = df)
  
  quadratic <- glm.nb(quadratic_formula,
                      data = df)
  
  # linear <- glmmTMB(linear_formula,
  #                   family = nbinom2,
  #                   data = df,
  #                   control = glmmTMBControl(optimizer = optim, optArgs = list(method = "BFGS")))
  # 
  # quadratic <- glmmTMB(quadratic_formula,
  #                      family = nbinom2,
  #                      data = df,
  #                      control = glmmTMBControl(optimizer = optim, optArgs = list(method = "BFGS")))

  list(linear = linear,
    quadratic = quadratic)
  
}

linear_vs_quadratic_results <- map_dfr(norm_vars, function(norm_var) {
  
  ww_hosp_model_exclude %>%
    group_by(utility, pcr_target) %>%
    group_nest() %>%
    mutate(
      models = map(data, linear_vs_quadratic, conc_var = norm_var),
      aic = map(models, ~ tibble(
        model = c("linear", "quadratic"),
        AIC = c(
          AIC(.x$linear),
          AIC(.x$quadratic)
        )
      ))
    ) %>%
    #select(utility, pcr_target, aic) %>%
    unnest(aic) %>%
    mutate(normalization = norm_var)
  
})

detach("package:MASS")

# quadratic is primary analysis
# linear + detect is sensitivity
# all in glm.nb, generate convergence diagnostics and flag the ones that give a warning

linear_vs_quadratic_results_glmnb <- linear_vs_quadratic_results %>%
  select(utility, pcr_target, normalization, model, AIC)

write.csv(linear_vs_quadratic_results_glmnb, "DataProcessed/linear_vs_quadratic_results_glmnb.csv")

library(MASS)
glm.nb(admissions_weekly ~ poly(WVAL_conc_weekly, 2, raw = TRUE) + admissions_weekly_w_1,
       data = ww_hosp_model_exclude %>% filter(utility == "Boulder",
                                               pcr_target == "fluav"))

library(glmmTMB)
library(dplyr)
library(purrr)
library(tibble)

linear_vs_quadratic <- function(df, conc_var) {
  
  # ---- 1. safety checks ----
  if (nrow(df) < 10) {
    return(list(
      linear = NULL,
      quadratic = NULL,
      status = "too_few_rows"
    ))
  }
  
  # ---- 2. scale predictor safely ----
  df <- df %>%
    mutate(
      conc_z = as.numeric(scale(.data[[conc_var]]))
    )
  
  # ---- 3. formulas (stable specification) ----
  linear_formula <- admissions_weekly ~ conc_z + detect_weekly + admissions_weekly_w_1
  quadratic_formula <- admissions_weekly ~ poly(conc_z, 2, raw = TRUE) + admissions_weekly_w_1
  
  # ---- 4. safe model fitting wrapper ----
  fit_safe <- function(formula, data) {
    
    tryCatch(
      {
        model <- glmmTMB(
          formula,
          family = nbinom2,
          data = data,
          control = glmmTMBControl(
            optimizer = optim,
            optArgs = list(method = "BFGS")
          )
        )
        
        list(
          model = model,
          converged = is.null(model$sdr$pdHess) || model$sdr$pdHess,
          warning = NULL,
          error = NULL
        )
        
      },
      warning = function(w) {
        # capture warning but continue
        return(list(
          model = NULL,
          converged = NA,
          warning = conditionMessage(w),
          error = NULL
        ))
      },
      error = function(e) {
        return(list(
          model = NULL,
          converged = FALSE,
          warning = NULL,
          error = conditionMessage(e)
        ))
      }
    )
  }
  
  # ---- 5. fit models ----
  linear <- fit_safe(linear_formula, df)
  quadratic <- fit_safe(quadratic_formula, df)
  
  # ---- 6. return structured output ----
  list(
    linear = linear$model,
    quadratic = quadratic$model,
    diagnostics = tibble(
      model = c("linear", "quadratic"),
      converged = c(linear$converged, quadratic$converged),
      warning = c(linear$warning, quadratic$warning),
      error = c(linear$error, quadratic$error)
    )
  )
}

safe_AIC <- function(model) {
  if (is.null(model)) return(NA_real_)
  tryCatch(AIC(model), error = function(e) NA_real_)
}

linear_vs_quadratic_results <- purrr::map_dfr(norm_vars, function(norm_var) {
  
  ww_hosp_model_exclude %>%
    group_by(utility, pcr_target) %>%
    group_nest() %>%
    mutate(
      models = map(data, linear_vs_quadratic, conc_var = norm_var),
      
      aic = map(models, ~ tibble(
        model = c("linear", "quadratic"),
        AIC = c(
          safe_AIC(.x$linear),
          safe_AIC(.x$quadratic)
        )
      )),
      
      diagnostics = map(models, "diagnostics")
    ) %>%
    
    unnest(aic) %>%
    mutate(normalization = norm_var)
  
})

diag_tbl <- linear_vs_quadratic_results %>%
  select(utility, pcr_target, normalization, diagnostics) %>%
  group_by(utility, pcr_target, normalization) %>%
  slice_tail(n = 1) %>%
  tidyr::unnest(diagnostics) %>%
  select(-warning)

write.csv(diag_tbl, "DataProcessed/diag_tbl.csv")

glmmTMB(admissions_weekly ~ WVAL_conc_weekly + detect_weekly + admissions_weekly_w_1,
        family = nbinom2,
        data = ww_hosp_model_exclude %>% filter(utility == "Aspen",
                                                pcr_target == "fluav"),
        control = glmmTMBControl(
          optimizer = optim,
          optArgs = list(method = "BFGS")))

glmmTMB(admissions_weekly ~ WVAL_conc_weekly + I(WVAL_conc_weekly^2) + admissions_weekly_w_1,
        family = nbinom2,
        data = ww_hosp_model_exclude %>% filter(utility == "Aspen",
                                                pcr_target == "fluav"),
        control = glmmTMBControl(
          optimizer = optim,
          optArgs = list(method = "BFGS")))

linear_vs_quadratic_results <- linear_vs_quadratic_results %>%
  group_by(utility, pcr_target, normalization) %>%
  mutate(AIC_diff = AIC - lag(AIC)) %>%
  select(utility, pcr_target, normalization, model, AIC, AIC_diff)

write.csv(linear_vs_quadratic_results, "DataProcessed/linear_vs_quadratic_results.csv")

test_metrocc_flu_mobile <- linear_vs_quadratic(df = ww_hosp_model_exclude %>% filter(utility == "Metro Wastewater RWHTF - CC",
                                                                                     pcr_target == "fluav"),
                                               conc_var = "mobile_norm_conc_weekly")
test_metrocc_flu_mobile

linear_vs_quadratic_results %>%
  filter(AIC_diff < -2) %>%
  nrow()

# grouped_data <- ww_hosp_model_exclude %>%
#   group_by(utility, pcr_target) %>%
#   group_nest()
# 
# models_df_raw <- ww_hosp_model_exclude %>%
#   group_by(utility, pcr_target) %>%
#   group_nest() %>%
#   mutate(models = map(data, compare_model_type_raw))
# 
# aic_long_raw <- models_df_raw %>%
#   mutate(
#     aic = map(models, ~ tibble(
#       model = c("linear_detect", "quadratic"),
#       aic = AIC(.x$linear_detect, .x$quadratic)
#     ))
#   ) %>%
#   select(utility, pcr_target, aic) %>%
#   unnest(aic) %>%
#   rename(AIC_raw = aic)
# 
# 
# 
# lrt_df <- models_df %>%
#   mutate(
#     with_quad = map(models, ~ anova(.x$with_quad, .x$without_quad, test="LRT")),
#   ) %>%
#   select(utility, pcr_target, with_quad)
# 
# extract_lrt <- function(x, prefix) {
#   tibble(
#     #!!paste0(prefix, "_chisq") := x$Chisq[2],
#     #!!paste0(prefix, "_df")    := x$`Chi Df`[2],
#     !!paste0(prefix, "_lrt_p")     := x$`Pr(>Chisq)`[2]
#   )
# }
# 
# lrt_clean <- lrt_df %>%
#   mutate(
#     quad = map(with_quad, ~ extract_lrt(.x, "with_quad")),
#   ) %>%
#   select(utility, pcr_target, quad) %>%
#   unnest(c(quad))
# 
# write.csv(lrt_clean, "DataProcessed/with_quad_lrt_no_detect.csv")

###


  
#################################### GARBAGE BELOW ###############################################


# # run models for each pathogen and sewershed
# models_covid_aspen_hightourism <- lapply(norm_vars, fit_nb_model, target = "sars-cov-2", sewershed = "Aspen")
# models_covid_alamosa_other <- lapply(norm_vars, fit_nb_model, target = "sars-cov-2", sewershed = "Alamosa")
# models_covid_boulder_other <- lapply(norm_vars, fit_nb_model, target = "sars-cov-2", sewershed = "Boulder")
# 
# models_flu_aspen_hightourism <- lapply(norm_vars, fit_nb_model, target = "fluav", sewershed = "Aspen")
# models_flu_alamosa_other <- lapply(norm_vars, fit_nb_model, target = "fluav", sewershed = "Alamosa")
# models_flu_boulder_other <- lapply(norm_vars, fit_nb_model, target = "fluav", sewershed = "Boulder")
# 
# 
# # function to generate model performance metrics
# compare_models_all <- function(model_list) {
#   data.frame(normalization = norm_names,
#              AIC = sapply(model_list, AIC))
#              #beta = sapply(model_list, \(m) summary(m)$coefficients$cond[2,1]),
#              #p_value = sapply(model_list, \(m) summary(m)$coefficients$cond[2,4]))
# }
# 
# # compare models
# model_comparisons <- mget(ls(pattern = "^models_")) %>%
#   map(~ compare_models_all(.x) %>% arrange(AIC))
# 
# all_models_tbl_primary <- imap_dfr(model_comparisons, ~ {
#   .x %>% 
#     mutate(
#       model_set = .y)
# })
# 
# all_models_tbl_primary <- all_models_tbl_primary %>%
#   separate(model_set, into = c("model", "pathogen", "sewershed", "category"), sep = "_", remove = FALSE) %>%
#   select(pathogen, sewershed, category, normalization, AIC) #%>%
#   # mutate(irr = round(exp(beta), 2),
#   #        p_value = round(p_value, 3),
#   #        significant = ifelse(p_value < 0.05, "yes", "no"))
# 
# write.csv(all_models_tbl_primary, "DataProcessed/model_comparisons_all_primary_20260323.csv", row.names = FALSE)
# 
# 
# 
# 
# ############################################################################
# 
# 
# 
# 
# # function to fit NB model to any PCR target (all categories)
# fit_nb_model_all <- function(var, target) {
#   glmmTMB(formula = build_formula(var),
#           family = nbinom2,
#           data = dplyr::filter(ww_hosp_scale, pcr_target == target),
#           control = glmmTMBControl(optimizer = optim, optArgs = list(method = "BFGS")))
# }
# 
# # function to fit NB model to any PCR target (each category)
# fit_nb_model_cat <- function(var, target, ww_cat) {
#   glmmTMB(formula = build_formula(var),
#           family = nbinom2,
#           data = dplyr::filter(ww_hosp_scale, pcr_target == target, category == ww_cat),
#           control = glmmTMBControl(optimizer = optim, optArgs = list(method = "BFGS")))
# }
# 
# # run models for each pathogen and category
# models_covid_all <- lapply(norm_vars, fit_nb_model_all, target = "sars-cov-2")
# models_covid_hightourism <- lapply(norm_vars, fit_nb_model_cat, target = "sars-cov-2", ww_cat = "High Tourism")
# models_covid_metro <- lapply(norm_vars, fit_nb_model_cat, target = "sars-cov-2", ww_cat = "Metro")
# models_covid_other <- lapply(norm_vars, fit_nb_model_cat, target = "sars-cov-2", ww_cat = "Other")
# 
# models_flu_all <- lapply(norm_vars, fit_nb_model_all, target = "fluav")
# models_flu_hightourism <- lapply(norm_vars, fit_nb_model_cat, target = "fluav", ww_cat = "High Tourism")
# models_flu_metro <- lapply(norm_vars, fit_nb_model_cat, target = "fluav", ww_cat = "Metro")
# models_flu_other <- lapply(norm_vars, fit_nb_model_cat, target = "fluav", ww_cat = "Other")
# 
# models_rsv_all <- lapply(norm_vars, fit_nb_model_all, target = "rsv")
# models_rsv_hightourism <- lapply(norm_vars, fit_nb_model_cat, target = "rsv", ww_cat = "High Tourism")
# models_rsv_metro <- lapply(norm_vars, fit_nb_model_cat, target = "rsv", ww_cat = "Metro")
# models_rsv_other <- lapply(norm_vars, fit_nb_model_cat, target = "rsv", ww_cat = "Other")
# 
# # function to generate model performance metrics
# compare_models_all <- function(model_list) {
#   data.frame(normalization = norm_names,
#              AIC = sapply(model_list, AIC),
#              beta = sapply(model_list, \(m) summary(m)$coefficients$cond[2,1]),
#              p_value = sapply(model_list, \(m) summary(m)$coefficients$cond[2,4]))
# }
# 
# # compare models
# model_comparisons <- mget(ls(pattern = "^models_")) %>%
#   map(~ compare_models_all(.x) %>% arrange(AIC))
# 
# all_models_tbl_primary <- imap_dfr(model_comparisons, ~ {
#   .x %>% 
#     mutate(
#       model_set = .y)
# })
# 
# all_models_tbl_primary <- all_models_tbl_primary %>%
#   separate(model_set, into = c("model", "pathogen", "category"), sep = "_", remove = FALSE) %>%
#   select(pathogen, category, normalization, AIC, beta, p_value) %>%
#   mutate(irr = round(exp(beta), 2),
#          p_value = round(p_value, 3),
#          significant = ifelse(p_value < 0.05, "yes", "no"))
# 
# a <- all_models_tbl_primary %>%
#   filter(category == "all")
# 
# b <- all_models_tbl_primary %>%
#   filter(category != "all")
# 
# write.csv(all_models_tbl_primary, "DataProcessed/model_comparisons_all_primary_20260323.csv", row.names = FALSE)
# 
# all_models_irr_primary <- all_models_tbl_primary %>%
#   arrange(pathogen, category, match(normalization, c("raw", "flow", "mobile", "combo", "wval"))) %>%
#   select(pathogen, category, normalization, irr) %>%
#   pivot_wider(names_from = normalization, values_from = irr)
# 
# write.csv(all_models_irr_primary, "DataProcessed/all_models_irr_primary_20260323.csv")
# 
# all_models_aic_primary <- all_models_tbl_primary %>%
#   group_by(pathogen, category) %>%
#   mutate(rank = ifelse(is.na(AIC), NA, row_number()),
#          delta_AIC_next = AIC - lag(AIC),
#          delta_AIC_cumul = AIC - min(AIC),
#          AIC_weight = round(exp(-0.5*delta_AIC_cumul) / sum(exp(-0.5*delta_AIC_cumul)), 3)) %>%
#   select(pathogen, category, normalization, contains("AIC"), rank)
# 
# write.csv(all_models_aic_primary, "DataProcessed/all_models_aic_primary_20260323.csv")
# 
# all_models_rank_primary <- all_models_aic_primary %>%
#   arrange(pathogen, category, match(normalization, c("raw", "flow", "mobile", "combo", "wval"))) %>%
#   select(pathogen, category, normalization, rank) %>%
#   pivot_wider(names_from = normalization, values_from = rank)
# 
# write.csv(all_models_rank_primary, "DataProcessed/all_models_rank_primary_20260323.csv")
# 
# 
# 
# 
# 
# 
# ###############################################
# 
# 
# 
# 
# 
# 
# 
# 
# plot_norm_methods <- list()
# for(i in 1:length(titles)){
#   plot_norm_methods[[i]] <- ggplot(data = ww_hosp_scale_long %>% filter(pcr_target == pcr_targets[i])) +
#     geom_line(aes(x = week_end, y = value, group = norm_method, color = norm_method)) +
#     labs(x = "Week End Date", y = "Z-Scored Value", color = "Normalization Method") +
#     ggtitle(paste("Normalized", titles[i], "Viral Concentrations")) +
#     scale_x_date(date_labels = "%b %Y") +
#     scale_color_manual(labels = c("Combination-Normalized", "Flow-Normalized",
#                                   "Mobile Device-Normalized", "Raw Concentration", "WVAL"),
#                        values = c("violet", "turquoise2", "firebrick",
#                                   "darkorchid4" , "darkorange")) +
#     guides(color = guide_legend(override.aes = list(linewidth = 8), nrow = 5)) +
#     facet_wrap(~utility, scales = "free_y") +
#     theme
#   ggsave(paste0("Figures/plot_norm_methods_", pcr_targets[i], ".png"), width = 17, height = 9, plot = plot_norm_methods[[i]])
# }
# 
# plot_norm_methods_anon <- list()
# for(i in 1:length(titles)){
#   plot_norm_methods_anon[[i]] <- ggplot(data = ww_hosp_scale_long %>% filter(pcr_target == pcr_targets[i])) +
#     geom_line(aes(x = week_end, y = value, group = norm_method, color = norm_method)) +
#     labs(x = "Week End Date", y = "Z-Scored Value", color = "Normalization Method") +
#     ggtitle(paste("Normalized", titles[i], "Viral Concentrations")) +
#     scale_x_date(date_labels = "%b %Y") +
#     scale_color_manual(labels = c("Combination-Normalized", "Flow-Normalized",
#                                   "Mobile Device-Normalized", "Raw Concentration", "WVAL"),
#                        values = c("violet", "turquoise2", "firebrick",
#                                   "darkorchid4" , "darkorange")) +
#     guides(color = guide_legend(override.aes = list(linewidth = 8), nrow = 5)) +
#     facet_wrap(~sewershed_code_factor, scales = "free_y") +
#     theme
#   ggsave(paste0("Figures/plot_norm_methods_anon_", pcr_targets[i], ".png"), width = 17, height = 9, plot = plot_norm_methods_anon[[i]])
# }
# 
# plot_hosp_ts <- list()
# for(i in 1:length(titles)){
#   plot_hosp_ts[[i]] <- ggplot(data = ww_hosp_model %>% filter(pcr_target == pcr_targets[i])) +
#     geom_bar(aes(x = week_end, y = admissions_weekly), stat = "identity") +
#     labs(x = "Week End Date", y = "Weekly Total Hospital Admissions") +
#     ggtitle(paste(titles[i], "Hospital Admissions Over Time")) +
#     scale_x_date(date_labels = "%b %Y") +
#     facet_wrap(~utility, scales = "free_y") +
#     theme
#   ggsave(paste0("Figures/plot_hosp_ts_", pcr_targets[i], ".png"), width = 17, height = 9, plot = plot_hosp_ts[[i]])
# }
# 
# 
# 
# 
# 
# 
# 
# 
# 
#   
# # hosp_daily <- hosp %>%
# #   group_by(utility, pcr_target) %>%
# #   complete(date = seq(min(date), max(date), by = "day"))
# 
# ww_hosp_daily <- merge(ww_daily, hosp_daily, c("pcr_target", "utility", "date"), all = TRUE) %>%
#   select(-lab_phase, -sample_type) %>%
#   filter(utility %in% sentinel$utility) %>%
#   left_join(visit) %>%
#   mutate(mobile_end_date = max(visit$date)) %>%
#   filter(date <= mobile_end_date) %>%
#   group_by(utility, pcr_target) %>%
#   fill(lod_sewage, category, census_population)
# 
# ww_hosp_daily_fill <- ww_hosp_daily %>%
#   group_by(utility, pcr_target) %>%
#   mutate(pcr_target_avg_conc = na.approx(pcr_target_avg_conc, na.rm = FALSE),
#          WVAL_conc = na.approx(WVAL_conc, na.rm = FALSE),
#          flow_rate = na.approx(flow_rate, na.rm = FALSE),
#          admissions = replace_na(admissions, 0)) %>%
#   fill(pcr_target_avg_conc, .direction = "downup") %>%
#   fill(WVAL_conc, .direction = "downup") %>%
#   fill(flow_rate, .direction = "downup") %>%
#   mutate(detect = ifelse(pcr_target_avg_conc <= lod_sewage, 0, 1)) %>%
#   select(-pcr_target_below_lod)
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# fit_zinb <- glmmTMB(admissions ~ detect * log_raw_conc_w_1 + offset(log(census_population)) + (1 | utility),
#                     family = nbinom2(),
#                     data = ww_hosp_scale %>% filter(pcr_target == "sars-cov-2"))
# 
# summary(fit_zinb)
# 
# formula <- bf(admissions ~ detect + detect:log_raw_conc_w_1 + offset(log(census_population)) + (1 | utility),
#               hu ~ detect + (1 | utility))
# 
# fit <- brm(formula = formula,
#   data = ww_hosp_scale,
#   family = hurdle_negbinomial(),
#   chains = 2,
#   cores = 2,
#   iter = 2000,
#   warmup = 1000,
#   threads = threading(2),
#   control = list(adapt_delta = 0.99,
#                  max_treedepth = 12),
#   seed = 123)
# 
# 
# 
# 
# 
# 
# 
# 
# # formula builder
# build_formula <- function(var) {
#   as.formula(paste("admissions ~", var,
#                    "+ detect + offset(log(census_population)) + (1 | utility)"))
# }
# 
# fit_raw <- brm(
#   formula = formula_simple,   # or formula_simple
#   data = ww_hosp_scale,
#   family = hurdle_negbinomial(),
#   chains = 4,
#   cores = 4,
#   iter = 4000,
#   warmup = 2000,
#   control = list(
#     adapt_delta = 0.99,
#     max_treedepth = 12
#   ),
#   seed = 123
# )
# 
# 
# 
# 
# 
# 
# 
# fit <- brm(bf(admissions ~ detect + log_raw_conc_w_1 + offset(log(census_population)) + (1 | utility)),
#   data = ww_hosp_scale %>% filter(pcr_target == "sars-cov-2"),
#   family = negbinomial(),
#   prior = c(prior(normal(0, 0.5), class = "b"),
#             prior(normal(0, 1), class = "Intercept"),
#             prior(exponential(1), class = "shape")),
#   iter = 4000,
#   warmup = 2000,
#   chains = 4,
#   cores = 4,
#   control = list(adapt_delta = 0.995, max_treedepth = 15))
# 
# 
# 
# fit1 <- brm(bf(admissions ~ detect + flow_norm_conc_w_1 + offset(log(census_population)) + (1 | utility)),
#             data = ww_hosp_scale %>% filter(pcr_target == "sars-cov-2"),
#             family = negbinomial(),
#             chains = 4, cores = 4, iter = 4000,
#             control = list(adapt_delta = 0.95)
# )
# 
# fit2 <- brm(bf(admissions ~ detect + mobile_norm_conc_w_1 + offset(log(census_population)) + (1 | utility)),
#             data = ww_hosp_scale %>% filter(pcr_target == "sars-cov-2"),
#             family = negbinomial(),
#             chains = 4, cores = 4, iter = 4000,
#             control = list(adapt_delta = 0.95)
# )
# 
# fit3 <- brm(bf(admissions ~ detect + raw_conc_w_1 + offset(log(census_population)) + (1 | utility)),
#            data = ww_hosp_scale %>% filter(pcr_target == "fluav"),
#            family = negbinomial(),
#            chains = 4, cores = 4, iter = 4000,
#            control = list(adapt_delta = 0.95)
# )
# 
# fit4 <- brm(bf(admissions ~ detect + flow_norm_conc_w_1 + offset(log(census_population)) + (1 | utility)),
#             data = ww_hosp_scale %>% filter(pcr_target == "fluav"),
#             family = negbinomial(),
#             chains = 4, cores = 4, iter = 4000,
#             control = list(adapt_delta = 0.95)
# )
# 
# fit5 <- brm(bf(admissions ~ detect + mobile_norm_conc_w_1 + offset(log(census_population)) + (1 | utility)),
#             data = ww_hosp_scale %>% filter(pcr_target == "fluav"),
#             family = negbinomial(),
#             chains = 4, cores = 4, iter = 4000,
#             control = list(adapt_delta = 0.95)
# )
# 
# 
# # organize normalization variables
# norm_vars <- c("raw_conc_w_1", "flow_norm_conc_w_1", "mobile_norm_conc_w_1", "combo_norm_conc_w_1")
# 
# targets <- c("sars-cov-2", "fluav", "rsv")
# 
# # formula builder
# build_formula <- function(var) {
#   as.formula(paste("admissions ~", var,
#                    "+ detect + offset(log(census_population)) + (1 | utility)"))
# }
# 
# # model fitter
# fit_bayes_model_all <- function(var, target) {
#   brm(formula = build_formula(var),
#       data = filter(ww_hosp_scale, pcr_target == target),
#       family = negbinomial(),
#       chains = 4,
#       cores = 4,
#       iter = 4000,
#       warmup = 1000,
#       backend = "cmdstanr",
#       control = list(adapt_delta = 0.95),
#       seed = 686)
# }
# 
# models_covid_all <- lapply(norm_vars, fit_bayes_model_all, target = "sars-cov-2")
# 
# 
# ww_hosp_norm <- ww_hosp %>%
#   mutate(detect = if_else(pcr_target_below_lod == "yes", 0, 1),
#          is_below = pcr_target_below_lod == "yes",
#          # base LODs
#          lod_flow = lod_sewage * flow_rate,
#          lod_mobile = lod_sewage / total_devices_daily,
#          lod_combo = lod_sewage * (flow_rate / total_devices_daily),
#          # base concentrations
#          raw_conc = pcr_target_avg_conc,
#          flow_conc = pcr_target_avg_conc * flow_rate,
#          mobile_conc = pcr_target_avg_conc / total_devices_daily,
#          combo_conc = pcr_target_avg_conc * (flow_rate / total_devices_daily)) %>%
#   # generate all substitution variants
#   {
#     factors <- c(sub1 = 1, sub2 = 1/2, sub3 = 1/sqrt(2))
#     
#     reduce(names(factors), function(df, nm) {
#       f <- factors[[nm]]
#       
#       df %>%
#         mutate(!!paste0("raw_conc_", nm) := sub_fun(raw_conc, lod_sewage, is_below, f),
#                !!paste0("flow_norm_conc_", nm) := sub_fun(flow_conc, lod_flow, is_below, f),
#                !!paste0("mobile_norm_conc_", nm) := sub_fun(mobile_conc, lod_mobile, is_below, f),
#                !!paste0("combo_norm_conc_", nm) := sub_fun(combo_conc, lod_combo, is_below, f))
#     }, .init = .)
#   } %>%
# 
# 
# 
# 
# 
# # organize normalization variables
# norm_vars <- list(sub1 = c("raw_conc_sub1_w_1", "flow_norm_conc_sub1_w_1", "mobile_norm_conc_sub1_w_1", "combo_norm_conc_sub1_w_1"),
#                   sub2 = c("raw_conc_sub2_w_1", "flow_norm_conc_sub2_w_1", "mobile_norm_conc_sub2_w_1", "combo_norm_conc_sub2_w_1"),
#                   sub3 = c("raw_conc_sub3_w_1", "flow_norm_conc_sub3_w_1", "mobile_norm_conc_sub3_w_1", "combo_norm_conc_sub3_w_1"))
# 
# targets <- c("sars-cov-2", "fluav", "rsv")
# 
# # formula builder
# build_formula <- function(var) {
#   as.formula(paste("admissions ~", var,
#                    "+ detect + offset(log(census_population)) + (1 | utility)"))
# }
# 
# # model fitter
# fit_bayes_model_all <- function(var, target) {
#   brm(formula = build_formula(var),
#       data = filter(ww_hosp_scale, pcr_target == target),
#       family = negbinomial(),
#       chains = 4,
#       cores = 4,
#       iter = 4000,
#       warmup = 1000,
#       backend = "cmdstanr",
#       control = list(adapt_delta = 0.9),
#       seed = 686)
# }
# 
# # FULL pipeline across targets + substitution methods
# results <- map(targets, function(target) {
#   map(norm_vars, function(vars) {
#     
#     models <- map(vars, fit_bayes_model_all, target = target)
#     names(models) <- vars
#     
#     loos <- map(models, loo)
#     
#     list(
#       models = models,
#       loos = loos,
#       comparison = loo_compare(loos)
#     )
#   })
# })
# # started at 8pm Wednesday 3/18, done by midnight
# names(results) <- targets
# names(results[[1]]) <- names(norm_vars)  # assign sub1/sub2/sub3 names
# 
# saveRDS(results, "DataProcessed/ww_bayes_results.rds")
# 
# results$`sars-cov-2`$sub1$comparison
# results$`sars-cov-2`$sub2$comparison
# results$`sars-cov-2`$sub3$comparison
# 
# results$`fluav`$sub1$comparison
# results$`fluav`$sub2$comparison
# results$`fluav`$sub3$comparison
# 
# results$`rsv`$sub1$comparison
# results$`rsv`$sub2$comparison
# results$`rsv`$sub3$comparison
# 
# ####################################################################################
# 
# 
# 
# 
# fit_bayes_model_all <- function(var, target) {
#   brm(formula = build_formula(var),
#       data = dplyr::filter(ww_hosp_scale, pcr_target == target),
#       family = negbinomial(),
#       chains = 2,
#       cores = 4,
#       iter = 2000,
#       control = list(adapt_delta = 0.95),
#       seed = 686)
# }
# 
# models_sub1_covid_all <- lapply(norm_vars_sub1, fit_bayes_model_all, target = "sars-cov-2")
# names(models_sub1_covid_all) <- norm_vars_sub1
# loos_sub1_covid_all <- lapply(models_sub1_covid_all, loo)
# loo_compare(loos_sub1_covid_all)
# 
# models_sub2_covid_all <- lapply(norm_vars_sub2, fit_bayes_model_all, target = "sars-cov-2")
# names(models_sub2_covid_all) <- norm_vars_sub2
# loos_sub2_covid_all <- lapply(models_sub2_covid_all, loo)
# loo_compare(loos_sub2_covid_all)
# 
# models_sub3_covid_all <- lapply(norm_vars_sub3, fit_bayes_model_all, target = "sars-cov-2")
# names(models_sub3_covid_all) <- norm_vars_sub3
# loos_sub3_covid_all <- lapply(models_sub3_covid_all, loo)
# loo_compare(loos_sub3_covid_all)
# 
# models_sub1_flu_all <- lapply(norm_vars_sub1, fit_bayes_model_all, target = "fluav")
# names(models_sub1_flu_all) <- norm_vars_sub1
# loos_sub1_flu_all <- lapply(models_sub1_flu_all, loo)
# loo_compare(loos_sub1_flu_all)
# 
# models_sub2_flu_all <- lapply(norm_vars_sub2, fit_bayes_model_all, target = "fluav")
# names(models_sub2_flu_all) <- norm_vars_sub2
# loos_sub2_flu_all <- lapply(models_sub2_flu_all, loo)
# loo_compare(loos_sub2_flu_all)
# 
# models_sub3_flu_all <- lapply(norm_vars_sub3, fit_bayes_model_all, target = "fluav")
# names(models_sub3_flu_all) <- norm_vars_sub3
# loos_sub3_flu_all <- lapply(models_sub3_flu_all, loo)
# loo_compare(loos_sub3_flu_all)
# 
# models_sub1_rsv_all <- lapply(norm_vars_sub1, fit_bayes_model_all, target = "rsv")
# names(models_sub1_rsv_all) <- norm_vars_sub1
# loos_sub1_rsv_all <- lapply(models_sub1_rsv_all, loo)
# loo_compare(loos_sub1_rsv_all)
# 
# models_sub2_rsv_all <- lapply(norm_vars_sub2, fit_bayes_model_all, target = "rsv")
# names(models_sub2_rsv_all) <- norm_vars_sub2
# loos_sub2_rsv_all <- lapply(models_sub2_rsv_all, loo)
# loo_compare(loos_sub2_rsv_all)
# 
# models_sub3_rsv_all <- lapply(norm_vars_sub3, fit_bayes_model_all, target = "rsv")
# names(models_sub3_rsv_all) <- norm_vars_sub3
# loos_sub3_rsv_all <- lapply(models_sub3_rsv_all, loo)
# loo_compare(loos_sub3_rsv_all)
# 
# 
# 
# 
# fit_bayes_model_cat <- function(var, target, ww_cat) {
#   brm(formula = build_formula(var),
#       data = dplyr::filter(ww_hosp_scale, pcr_target == target, category == ww_cat),
#       family = negbinomial(),
#       chains = 2,
#       cores = 4,
#       iter = 2000,
#       control = list(adapt_delta = 0.95),
#       seed = 686)
# }
# 
# 
# 
# 
# models_primary_covid_hightourism <- lapply(norm_vars, fit_bayes_model_cat, target = "sars-cov-2", ww_cat = "High Tourism")
# models_primary_covid_metro <- lapply(norm_vars, fit_bayes_model_cat, target = "sars-cov-2", ww_cat = "Metro")
# models_primary_covid_other <- lapply(norm_vars, fit_bayes_model_cat, target = "sars-cov-2", ww_cat = "Other")
# 
# models_primary_flu_all <- lapply(norm_vars, fit_bayes_model_all, target = "fluav")
# models_primary_flu_hightourism <- lapply(norm_vars, fit_bayes_model_cat, target = "fluav", ww_cat = "High Tourism")
# models_primary_flu_metro <- lapply(norm_vars, fit_bayes_model_cat, target = "fluav", ww_cat = "Metro")
# models_primary_flu_other <- lapply(norm_vars, fit_bayes_model_cat, target = "fluav", ww_cat = "Other")
# 
# models_primary_rsv_all <- lapply(norm_vars, fit_bayes_model_all, target = "rsv")
# models_primary_rsv_hightourism <- lapply(norm_vars, fit_bayes_model_cat, target = "rsv", ww_cat = "High Tourism")
# models_primary_rsv_metro <- lapply(norm_vars, fit_bayes_model_cat, target = "rsv", ww_cat = "Metro")
# models_primary_rsv_other <- lapply(norm_vars, fit_bayes_model_cat, target = "rsv", ww_cat = "Other")
# 
# # function to generate model performance metrics
# compare_models_all <- function(model_list) {
#   data.frame(normalization = norm_names,
#              AIC = sapply(model_list, AIC),
#              beta = sapply(model_list, \(m) summary(m)$coefficients$cond[2,1]),
#              p_value = sapply(model_list, \(m) summary(m)$coefficients$cond[2,4]))
# }
# 
# # compare models
# model_comparisons <- mget(ls(pattern = "^models_primary_")) %>%
#   map(~ compare_models_all(.x) %>% arrange(AIC))
# 
# all_models_tbl_primary <- imap_dfr(model_comparisons, ~ {
#   .x %>% 
#     mutate(
#       model_set = .y)
# })
# 
# all_models_tbl_primary <- all_models_tbl_primary %>%
#   separate(model_set, into = c("model", "analysis", "pathogen", "category"), sep = "_", remove = FALSE) %>%
#   select(analysis, pathogen, category, normalization, AIC, beta, p_value) %>%
#   mutate(irr = round(exp(beta), 2),
#          p_value = round(p_value, 3),
#          significant = ifelse(p_value < 0.05, "yes", "no"))
# 
# write.csv(all_models_tbl_primary, "DataProcessed/model_comparisons_all_primary.csv", row.names = FALSE)
# 
# 
# 
# 
# models_primary_covid_hightourism <- lapply(norm_vars, fit_nb_model_cat, target = "sars-cov-2", ww_cat = "High Tourism")
# models_primary_covid_metro <- lapply(norm_vars, fit_nb_model_cat, target = "sars-cov-2", ww_cat = "Metro")
# models_primary_covid_other <- lapply(norm_vars, fit_nb_model_cat, target = "sars-cov-2", ww_cat = "Other")
# 
# 
# 
# # function to fit NB model to any PCR target (all categories)
# fit_nb_model_all <- function(var, target) {
#   glmmTMB(formula = build_formula(var),
#           family = nbinom2,
#           data = dplyr::filter(ww_hosp_weekly_scale, pcr_target == target))
# }
# 
# # function to fit NB model to any PCR target (each category)
# fit_nb_model_cat <- function(var, target, ww_cat) {
#   glmmTMB(formula = build_formula(var),
#           family = nbinom2,
#           data = dplyr::filter(ww_hosp_weekly_scale, pcr_target == target, category == ww_cat))
# }
# 
# fits <- lapply(predictors, function(pred, target, ww_cat) {
#   brm(formula = as.formula(paste0("admissions ~ ", pred, " + (1 | utility) + offset(log(census_population))")),
#       data = dplyr::filter(ww_hosp_scale, pcr_target == target, category == ww_cat),
#       family = negbinomial(),
#       chains = 4,
#       cores = 4,
#       iter = 4000,
#       control = list(adapt_delta = 0.95),
#       seed = 686)
# })
# 
# names(fits) <- predictors
# 
# loos <- lapply(fits, loo)
# loo_compare(loos)
# 
# 
# models_primary_covid_hightourism <- lapply(norm_vars, fit_nb_model_cat, target = "sars-cov-2", ww_cat = "High Tourism")
# models_primary_covid_metro <- lapply(norm_vars, fit_nb_model_cat, target = "sars-cov-2", ww_cat = "Metro")
# models_primary_covid_other <- lapply(norm_vars, fit_nb_model_cat, target = "sars-cov-2", ww_cat = "Other")
# 
# 
# fit <- brm(formula = admissions ~ raw_conc + (1 | utility) + offset(log(census_population)),
#            data = ww_hosp_scale,
#            family = negbinomial(),
#            chains = 4,
#            cores = 4,
#            iter = 4000,
#            control = list(adapt_delta = 0.95))
# 
# fit1 <- brm(formula = admissions ~ flow_norm_conc + (1 | utility) + offset(log(census_population)),
#             data = ww_hosp_scale,
#             family = negbinomial(),
#             chains = 4,
#             cores = 4,
#             iter = 4000,
#             control = list(adapt_delta = 0.95))
# 
# fit2 <- brm(formula = admissions ~ mobile_norm_conc + (1 | utility) + offset(log(census_population)),
#             data = ww_hosp_scale,
#             family = negbinomial(),
#             chains = 4,
#             cores = 4,
#             iter = 4000,
#             control = list(adapt_delta = 0.95))
# 
# loo(fit, fit1)
# 
# norm_vars <- c("raw_conc", "flow_norm_conc", "mobile_norm_conc", "combo_norm_conc")
# 
# fit_brms_model <- function(var_name, data) {
#   data_model <- ww_hosp_scale %>%
#     mutate(conc = ww_hosp_scale[[var_name]])
#   
#   brm(formula = admissions ~ conc + offset(log(census_population)) + (1 | utility),
#       data = data_model,
#       family = negbinomial(),
#       chains = 4,
#       cores = 4,
#       iter = 4000,
#       control = list(adapt_delta = 0.95),
#       seed = 686)
# }
# 
# fits <- map(norm_vars, ~ fit_brms_model(.x, ww_hosp_scale))
# 
# 
# 
# 
# 
# ww_remove_outliers <- ww1 %>%
#   group_by(utility, pcr_target) %>%
#   arrange(date) %>%
#   mutate(roll_mean = cummean(pcr_target_avg_conc),
#          n = row_number(),
#          roll_var = cumsum((pcr_target_avg_conc - roll_mean)^2) / pmax(n - 1, 1),
#          roll_sd = sqrt(roll_var),
#          z = (pcr_target_avg_conc - roll_mean) / roll_sd) %>%
#   filter(abs(z) <= 4 | is.na(z)) %>%
#   select(-(c(roll_mean, n, roll_var, roll_sd, z)))
# 
# write.csv(ww_remove_outliers, "DataProcessed/ww_remove_outliers_20260305.csv")
#   
# 
# # ww <- read.csv("DataRaw/2_CDPHE_Wastewater_Data_ 2025-08-04 .csv") %>%
# #   # cleaning and formatting
# #   mutate(utility = wwtp_name,
# #          wwtp_name = tolower(wwtp_name),
# #          pcr_target = tolower(pcr_target),
# #          date = as.Date(sample_collect_date, "%Y-%m-%d"),
# #          test_result_date = as.Date(test_result_date, "%Y-%m-%d")) %>%
# #   filter(utility %in% sentinel$utility,
# #          pcr_target %in% c("fluav", "rsv", "sars-cov-2")) %>%
# #          #date >= "2024-09-30") %>%
# #   left_join(flow_missing_impute) %>%
# #   # impute anything BLOD with half the LOD
# #   mutate(pcr_target_avg_conc = case_when(pcr_target_avg_conc <= lod_sewage ~ lod_sewage/2,
# #                                          TRUE ~ pcr_target_avg_conc),
# #          # impute missing and zero flow rates with recovered data
# #          flow_rate = case_when(flow_rate == 0 ~ flow_rate_impute,
# #   # set remaining zero flow rates to missing
# #   TRUE ~ flow_rate)) %>%
# #   select(-flow_rate_impute) %>%
# #   group_by(utility, pcr_target) %>%
# #   # impute the remaining missing flow rates with the average of the last two non-missing flow rates
# #   mutate(flow_rate_last2avg = slide_dbl(flow_rate,
# #                                         ~ mean(tail(na.omit(.x), 2)),
# #                                         .before = Inf,
# #                                         .complete = F),
# #          flow_rate = ifelse(is.na(flow_rate), flow_rate_last2avg, flow_rate)) %>%
# #   select(-flow_rate_last2avg) %>%
# #   left_join(sentinel)
# 
# 
# 
# # sampling frequency dataset
# samp_freq <- desc %>% slice(1)
# 
# 
# sample_type_plot <- list()
# for(i in 1:length(titles)){
#   sample_type_plot[[i]] <- ggplot(data = desc %>% filter(pcr_target == pcr_targets[i])) +
#                                   #aes(x = date, y = flow_rate_per_capita, color = sample_type)) +
#     geom_point(aes(x = date, y = flow_rate_per_capita, color = sample_type)) +
#     geom_line(aes(x = date, y = flow_rate_per_capita)) +
#     labs(x = NULL, y = "Daily Flow Rate per Capita (gal)", color = "Sample Type") +
#     ggtitle(paste("Flow Rate by Sampling Technique,", titles[i])) +
#     facet_wrap(~utility, scales = "free") +
#     theme +
#     theme(axis.text.x = element_text(angle = 0))
#   ggsave(paste0("Figures/sample_type_plot_", pcr_targets[i], ".png"), width = 17, height = 9, plot = sample_type_plot[[i]])
# }
# 
# 
# 
# # sample_ndays_plot <- list()
# # for(i in 1:length(titles)){
# #   sample_ndays_plot[[i]] <- ggplot(data = desc %>% filter(pcr_target == pcr_targets[i])) +
# #     geom_bar(aes(x = date, y = ndays), stat = "identity") +
# #     geom_line(aes(x = date, y = mean_ndays), color = "red") +
# #     labs(x = NULL, y = "Number of Days") +
# #     scale_y_continuous(limits = c(0, 21), breaks = seq(0, 21, 3)) +
# #     ggtitle(paste("Number of Days Between Sample Collection and Test Result,", titles[i])) +
# #     facet_wrap(~utility, scales = "free") +
# #     theme +
# #     theme(axis.text.x = element_text(angle = 0))
# #   ggsave(paste0("Figures/sample_ndays_plot_", pcr_targets[i], ".png"), width = 15, height = 9, plot = sample_ndays_plot[[i]])
# # }
# 
# targets <- c(COVID = "sars-cov-2",
#              FLUA = "fluav",
#              RSV = "rsv")
# 
# wval <- imap_dfr(targets , ~ {
#   read.csv(glue::glue("DataProcessed/2026_03_12_dataset/2026_03_12_{.y}_wval.csv")) %>%
#     mutate(pcr_target = .x,
#            week_end = as.Date(week_end, "%Y-%m-%d")) %>%
#     rename(WVAL_conc_weekly = WVAL,
#            utility = wwtp_name) %>%
#     select(week_end, pcr_target, utility, WVAL_conc_weekly) %>%
#     filter(week_end <= "2026-01-31")
# })
# 
# ww_weekly <- desc %>%
#   group_by(utility, pcr_target) %>%
#   rename(raw_conc = conc_obs) %>%
#   # create normalized concentration values
#   mutate(flow_norm_conc = raw_conc*flow_rate,
#          mobile_norm_conc = raw_conc/total_devices_daily,
#          combo_norm_conc = raw_conc*(flow_rate/total_devices_daily)) %>%
#   # aggregate to weekly
#   group_by(utility, pcr_target, week_end) %>%
#   mutate(raw_conc_weekly = median(raw_conc),
#          flow_norm_conc_weekly = median(flow_norm_conc),
#          mobile_norm_conc_weekly = median(mobile_norm_conc),
#          combo_norm_conc_weekly = median(combo_norm_conc)) %>%
#   slice(1) %>%
#   select(utility, category, pcr_target, week_end, contains("conc_weekly"), census_population) %>%
#   left_join(wval)
# 
# 
# hosp_weekly <- hosp %>%
#   #mutate(week_end = ceiling_date(date, "weeks", week_start = getOption("lubridate.week.start", 6))) %>%
#   group_by(utility, pcr_target, week_end) %>%
#   mutate(admissions_weekly = sum(admissions)) %>%
#   slice(1) %>%
#   filter(week_end <= "2026-01-31",
#          utility %in% sentinel$utility) %>%
#   select(-c(admissions, date))
# 
# ww_hosp_weekly <- merge(ww_weekly, hosp_weekly, c("pcr_target", "utility", "week_end"), all = T) %>%
#   mutate(admissions_weekly = replace_na(admissions_weekly, 0)) %>%
#   group_by(utility, pcr_target) %>%
#   fill(category, census_population, raw_conc_weekly, flow_norm_conc_weekly,
#        mobile_norm_conc_weekly, combo_norm_conc_weekly, WVAL_conc_weekly) %>%
#   # create lagged variables
#   mutate(admissions_weekly_w_1 = lag(admissions_weekly, 1)) %>%
#   mutate(across(c(raw_conc_weekly, flow_norm_conc_weekly, mobile_norm_conc_weekly,
#                   combo_norm_conc_weekly, WVAL_conc_weekly),
#                 list(w_1 = ~lag(.x, 1)),
#                 .names = "{.col}_{.fn}")) %>%
#   select(utility, category, census_population, pcr_target, contains("week"))
# 
# write.csv(ww_hosp_weekly, "DataProcessed/ww_hosp_weekly.csv")
# 
# scatter <- list()
# for(i in 1:length(titles)){
#   scatter[[i]] <- ggplot(data = ww_hosp_weekly %>% filter(pcr_target == pcr_targets[i]),
#                          aes(x = raw_conc_weekly, y = admissions_weekly)) +
#     geom_point() +
#     geom_smooth(method = "glm.nb", se = TRUE, color = "blue") +
#     labs(x = "Weekly Raw Concentration (copies/L)", y = "Weekly Total Hospital Admissions (n)") +
#     ggtitle(paste(titles[i], "Hospital Admissions vs. Raw Viral Concentrations")) +
#     scale_x_continuous(labels = scales::comma) +
#     facet_wrap(~utility, scales = "free") +
#     theme +
#     theme(axis.text.x = element_text(angle = 0))
#   ggsave(paste0("Figures/scatter_", pcr_targets[i], ".png"), width = 15, height = 9, plot = scatter[[i]])
# }
# 
# 
# # do we need zero inflation?
# m_nb <- glmmTMB(admissions_weekly ~ raw_conc_weekly +
#                   offset(log(census_population)) + (1 | utility),
#                 family = nbinom2,
#                 data = ww_hosp_weekly_scale %>% filter(pcr_target == "sars-cov-2"))
# 
# sim <- simulateResiduals(m_nb)
# plot(sim)
# testZeroInflation(sim) # p = 0.944, no evidence of excess/structural zeros, do not use ZI model
# 
# # check ACF of residuals
# acf(residuals(m_nb)) # some autocorrelation at 1-8 weeks, include lagged admissions as predictor
# 
# # variables to loop over
# norm_vars <- c("raw_conc_weekly", "flow_norm_conc_weekly", "mobile_norm_conc_weekly",
#                "combo_norm_conc_weekly", "WVAL_conc_weekly")
# 
# norm_names <- c("raw", "flow", "mobile", "combo", "wval")
# 
# # function to build NB model formula for any normalization method
# build_formula <- function(var) {
#   as.formula(paste("admissions_weekly ~", var,
#                    "+ admissions_weekly_w_1 + offset(log(census_population)) + (1 | utility)"))
# }
# 
# # function to fit NB model to any PCR target (all categories)
# fit_nb_model_all <- function(var, target) {
#   glmmTMB(formula = build_formula(var),
#           family = nbinom2,
#           data = dplyr::filter(ww_hosp_weekly_scale, pcr_target == target))
# }
# 
# # function to fit NB model to any PCR target (each category)
# fit_nb_model_cat <- function(var, target, ww_cat) {
#   glmmTMB(formula = build_formula(var),
#           family = nbinom2,
#           data = dplyr::filter(ww_hosp_weekly_scale, pcr_target == target, category == ww_cat))
# }
# 
# # run models for each pathogen and category
# models_primary_covid_all <- lapply(norm_vars, fit_nb_model_all, target = "sars-cov-2")
# models_primary_covid_hightourism <- lapply(norm_vars, fit_nb_model_cat, target = "sars-cov-2", ww_cat = "High Tourism")
# models_primary_covid_metro <- lapply(norm_vars, fit_nb_model_cat, target = "sars-cov-2", ww_cat = "Metro")
# models_primary_covid_other <- lapply(norm_vars, fit_nb_model_cat, target = "sars-cov-2", ww_cat = "Other")
# 
# models_primary_flu_all <- lapply(norm_vars, fit_nb_model_all, target = "fluav")
# models_primary_flu_hightourism <- lapply(norm_vars, fit_nb_model_cat, target = "fluav", ww_cat = "High Tourism")
# models_primary_flu_metro <- lapply(norm_vars, fit_nb_model_cat, target = "fluav", ww_cat = "Metro")
# models_primary_flu_other <- lapply(norm_vars, fit_nb_model_cat, target = "fluav", ww_cat = "Other")
# 
# models_primary_rsv_all <- lapply(norm_vars, fit_nb_model_all, target = "rsv")
# models_primary_rsv_hightourism <- lapply(norm_vars, fit_nb_model_cat, target = "rsv", ww_cat = "High Tourism")
# models_primary_rsv_metro <- lapply(norm_vars, fit_nb_model_cat, target = "rsv", ww_cat = "Metro")
# models_primary_rsv_other <- lapply(norm_vars, fit_nb_model_cat, target = "rsv", ww_cat = "Other")
# 
# # function to generate model performance metrics
# compare_models_all <- function(model_list) {
#   data.frame(normalization = norm_names,
#              AIC = sapply(model_list, AIC),
#              beta = sapply(model_list, \(m) summary(m)$coefficients$cond[2,1]),
#              p_value = sapply(model_list, \(m) summary(m)$coefficients$cond[2,4]))
# }
# 
# # compare models
# model_comparisons <- mget(ls(pattern = "^models_primary_")) %>%
#   map(~ compare_models_all(.x) %>% arrange(AIC))
# 
# all_models_tbl_primary <- imap_dfr(model_comparisons, ~ {
#   .x %>% 
#     mutate(
#       model_set = .y)
# })
# 
# all_models_tbl_primary <- all_models_tbl_primary %>%
#   separate(model_set, into = c("model", "analysis", "pathogen", "category"), sep = "_", remove = FALSE) %>%
#   select(analysis, pathogen, category, normalization, AIC, beta, p_value) %>%
#   mutate(irr = round(exp(beta), 2),
#          p_value = round(p_value, 3),
#          significant = ifelse(p_value < 0.05, "yes", "no"))
# 
# write.csv(all_models_tbl_primary, "DataProcessed/model_comparisons_all_primary.csv", row.names = FALSE)
# 
# all_models_irr_primary <- all_models_tbl_primary %>%
#   arrange(pathogen, category, match(normalization, c("raw", "flow", "mobile", "combo", "wval"))) %>%
#   select(pathogen, category, normalization, irr) %>%
#   pivot_wider(names_from = normalization, values_from = irr)
# 
# write.csv(all_models_irr_primary, "DataProcessed/all_models_irr_primary.csv")
# 
# all_models_aic_primary <- all_models_tbl_primary %>%
#   group_by(pathogen, category) %>%
#   mutate(rank = ifelse(is.na(AIC), NA, row_number()),
#          delta_AIC_next = AIC - lag(AIC),
#          delta_AIC_cumul = AIC - min(AIC),
#          AIC_weight = round(exp(-0.5*delta_AIC_cumul) / sum(exp(-0.5*delta_AIC_cumul)), 3)) %>%
#   select(pathogen, category, normalization, contains("AIC"), rank)
# 
# write.csv(all_models_aic_primary, "DataProcessed/all_models_aic_primary.csv")
# 
# all_models_rank_primary <- all_models_aic_primary %>%
#   arrange(pathogen, category, match(normalization, c("raw", "flow", "mobile", "combo", "wval"))) %>%
#   select(pathogen, category, normalization, rank) %>%
#   pivot_wider(names_from = normalization, values_from = rank)
# 
# write.csv(all_models_rank_primary, "DataProcessed/all_models_rank_primary.csv")
# 
# # SENSITIVITY ANALYSIS
# # variables to loop over
# norm_vars_sens <- c("raw_conc_weekly_w_1", "flow_norm_conc_weekly_w_1", "mobile_norm_conc_weekly_w_1",
#                "combo_norm_conc_weekly_w_1", "WVAL_conc_weekly_w_1")
# 
# norm_names_sens <- c("raw_w_1", "flow_w_1", "mobile_w_1", "combo_w_1", "wval_w_1")
# 
# # run models for each pathogen and category (sensitivity analysis)
# models_sensitivity_covid_all <- lapply(norm_vars_sens, fit_nb_model_all, target = "sars-cov-2")
# models_sensitivity_covid_hightourism <- lapply(norm_vars_sens, fit_nb_model_cat, target = "sars-cov-2", ww_cat = "High Tourism")
# models_sensitivity_covid_metro <- lapply(norm_vars_sens, fit_nb_model_cat, target = "sars-cov-2", ww_cat = "Metro")
# models_sensitivity_covid_other <- lapply(norm_vars_sens, fit_nb_model_cat, target = "sars-cov-2", ww_cat = "Other")
# 
# models_sensitivity_flu_all <- lapply(norm_vars_sens, fit_nb_model_all, target = "fluav")
# models_sensitivity_flu_hightourism <- lapply(norm_vars_sens, fit_nb_model_cat, target = "fluav", ww_cat = "High Tourism")
# models_sensitivity_flu_metro <- lapply(norm_vars_sens, fit_nb_model_cat, target = "fluav", ww_cat = "Metro")
# models_sensitivity_flu_other <- lapply(norm_vars_sens, fit_nb_model_cat, target = "fluav", ww_cat = "Other")
# 
# models_sensitivity_rsv_all <- lapply(norm_vars_sens, fit_nb_model_all, target = "rsv")
# models_sensitivity_rsv_hightourism <- lapply(norm_vars_sens, fit_nb_model_cat, target = "rsv", ww_cat = "High Tourism")
# models_sensitivity_rsv_metro <- lapply(norm_vars_sens, fit_nb_model_cat, target = "rsv", ww_cat = "Metro")
# models_sensitivity_rsv_other <- lapply(norm_vars_sens, fit_nb_model_cat, target = "rsv", ww_cat = "Other")
# 
# # compare models
# model_comparisons_sens <- mget(ls(pattern = "^models_sensitivity_")) %>%
#   map(~ compare_models_all(.x) %>% arrange(AIC))
# 
# all_models_tbl_sensitivity <- imap_dfr(model_comparisons_sens, ~ {
#   .x %>% 
#     mutate(
#       model_set = .y)
# })
# 
# all_models_tbl_sensitivity <- all_models_tbl_sensitivity %>%
#   separate(model_set, into = c("model", "analysis", "pathogen", "category"), sep = "_", remove = FALSE) %>%
#   select(analysis, pathogen, category, normalization, AIC, beta, p_value) %>%
#   mutate(irr = round(exp(beta), 2),
#          p_value = round(p_value, 3),
#          significant = ifelse(p_value < 0.05, "yes", "no"))
# 
# write.csv(all_models_tbl_sensitivity, "DataProcessed/model_comparisons_all_sensitivity.csv", row.names = FALSE)
# 
# all_models_irr_sensitivity <- all_models_tbl_sensitivity %>%
#   arrange(pathogen, category, match(normalization, c("raw", "flow", "mobile", "combo", "wval"))) %>%
#   select(pathogen, category, normalization, irr) %>%
#   pivot_wider(names_from = normalization, values_from = irr)
# 
# write.csv(all_models_irr_sensitivity, "DataProcessed/all_models_irr_sensitivity.csv")
# 
# all_models_aic_sensitivity <- all_models_tbl_sensitivity %>%
#   group_by(pathogen, category) %>%
#   mutate(rank = ifelse(is.na(AIC), NA, row_number()),
#          delta_AIC_next = AIC - lag(AIC),
#          delta_AIC_cumul = AIC - min(AIC),
#          AIC_weight = round(exp(-0.5*delta_AIC_cumul) / sum(exp(-0.5*delta_AIC_cumul)), 3)) %>%
#   select(pathogen, category, normalization, contains("AIC"), rank)
# 
# write.csv(all_models_aic_sensitivity, "DataProcessed/all_models_aic_sensitivity.csv")
# 
# all_models_rank_sensitivity <- all_models_aic_sensitivity %>%
#   arrange(pathogen, category, match(normalization, c("raw", "flow", "mobile", "combo", "wval"))) %>%
#   select(pathogen, category, normalization, rank) %>%
#   pivot_wider(names_from = normalization, values_from = rank)
# 
# write.csv(all_models_rank_sensitivity, "DataProcessed/all_models_rank_sensitivity.csv")
# 
# # ww_hosp_norm <- ww_hosp %>%
# #   mutate(detect = ifelse(pcr_target_below_lod == "yes", 0, 1),
# #          raw_conc_sub1 = ifelse(pcr_target_below_lod == "yes", lod_sewage, pcr_target_avg_conc),
# #          lod_flow = lod_sewage*flow_rate,
# #          lod_mobile = lod_sewage/total_devices_daily,
# #          lod_combo = lod_sewage*(flow_rate/total_devices_daily),
# #          flow_norm_conc_sub1 = ifelse(pcr_target_below_lod == "yes", lod_flow, pcr_target_avg_conc*flow_rate),
# #          mobile_norm_conc_sub1 = ifelse(pcr_target_below_lod == "yes", lod_mobile, pcr_target_avg_conc/total_devices_daily),
# #          combo_norm_conc_sub1 = ifelse(pcr_target_below_lod == "yes", lod_combo, pcr_target_avg_conc*(flow_rate/total_devices_daily)),
# #          raw_conc_sub2 = ifelse(pcr_target_below_lod == "yes", lod_sewage/2, pcr_target_avg_conc),
# #          flow_norm_conc_sub2 = ifelse(pcr_target_below_lod == "yes", lod_flow/2, pcr_target_avg_conc*flow_rate),
# #          mobile_norm_conc_sub2 = ifelse(pcr_target_below_lod == "yes", lod_mobile/2, pcr_target_avg_conc/total_devices_daily),
# #          combo_norm_conc_sub2 = ifelse(pcr_target_below_lod == "yes", lod_combo/2, pcr_target_avg_conc*(flow_rate/total_devices_daily)),
# #          raw_conc_sub3 = ifelse(pcr_target_below_lod == "yes", lod_sewage/sqrt(2), pcr_target_avg_conc),
# #          flow_norm_conc_sub3 = ifelse(pcr_target_below_lod == "yes", lod_flow/sqrt(2), pcr_target_avg_conc*flow_rate),
# #          mobile_norm_conc_sub3 = ifelse(pcr_target_below_lod == "yes", lod_mobile/sqrt(2), pcr_target_avg_conc/total_devices_daily),
# #          combo_norm_conc_sub3 = ifelse(pcr_target_below_lod == "yes", lod_combo/sqrt(2), pcr_target_avg_conc*(flow_rate/total_devices_daily))) %>%
# #   mutate(across(c(contains("conc")),
# #                 list(w_1 = ~lag(.x, 2)),
# #                 .names = "{.col}_{.fn}")) %>%
# #   select(utility, category, pcr_target, date, detect, contains("conc"), -contains("lod"),
# #          -contains("pcr_target_avg_conc"), admissions, census_population)
# 
# # 
# # 
# # 
# # cbg_2010 <- get_decennial(geography = "block group",
# #                           variables = "P001001",
# #                           state = "CO",
# #                           year = 2010,
# #                           geometry = TRUE) %>%
# #   rename(pop_2010 = value) %>%
# #   select(GEOID, pop_2010)
# # 
# # cbg_2019 <- get_acs(geography = "block group",
# #                     variables = "B01003_001",
# #                     state = "CO",
# #                     year = 2019,
# #                     geometry = FALSE) %>%
# #   rename(pop_2019 = estimate) %>%
# #   select(GEOID, pop_2019)
# # 
# # cbg_growth = merge(cbg_2010, cbg_2019, "GEOID", all = T) %>%
# #   mutate(growth_factor = pop_2019/pop_2010)
# # 
# # co_counties <- counties(state = "CO", cb = TRUE, class = "sf")
# # 
# # cbg_growth <- ggplot(cbg_growth, aes(fill = growth_factor)) +
# #   geom_sf(color = NA) +
# #   geom_sf(data = co_counties, fill = NA, color = "black", size = 0.1) +
# #   scale_fill_viridis_c(option = "magma", direction = -1); cbg_growth
# # 
# # ggsave("Figures/cbg_growth.png", width = 5, height = 4, plot = cbg_growth)
# 
# 
# 
# sample_type_plot <- list()
# for(i in 1:length(titles)){
#   sample_type_plot[[i]] <- ggplot(data = desc %>% filter(pcr_target == pcr_targets[i])) +
#     #aes(x = date, y = flow_rate_per_capita, color = sample_type)) +
#     geom_point(aes(x = date, y = flow_rate_per_capita, color = sample_type)) +
#     geom_line(aes(x = date, y = flow_rate_per_capita)) +
#     labs(x = NULL, y = "Daily Flow Rate per Capita (gal)", color = "Sample Type") +
#     ggtitle(paste("Flow Rate by Sampling Technique,", titles[i])) +
#     facet_wrap(~utility, scales = "free") +
#     theme +
#     theme(axis.text.x = element_text(angle = 0))
#   ggsave(paste0("Figures/sample_type_plot_", pcr_targets[i], ".png"), width = 17, height = 9, plot = sample_type_plot[[i]])
# }
# 
# hosp_covid <- bind_rows(hosp_covid1, hosp_covid2) %>%
#   mutate(pcr_target = "sars-cov-2") %>%
#   rename(utility = Region.Name,
#          admissions = Hospitalized.Count) %>%
#   filter(utility %in% sentinel$utility) %>%
#   select(utility, pcr_target, date, admissions, dataset) %>%
#   group_by(utility, pcr_target, date) %>%
#   # identify any duplicate observations (resulting from overlapping dates in the two datasets)
#   add_count() %>%
#   group_by(utility, pcr_target, date) %>%
#   slice(1)
# 
# hosp_flu1 <- read.csv("DataRaw/2026/Flu_agg.csv") %>%
#   filter(Region.Level == "Sewershed",
#          Pathogen.Name == "INFLUENZA A") %>%
#   mutate(date = as.Date(Event.Onset.Date, "%m/%d/%Y"),
#          dataset = 1)
# 
# hosp_flu2 <- read.csv("DataRaw/CSPH Data _ April 2026/Flu_agg_April_update.csv") %>%
#   filter(Region.Level == "Sewershed",
#          Pathogen.Name == "INFLUENZA A") %>%
#   mutate(date = as.Date(Event.Onset.Date, "%m/%d/%Y"),
#          dataset = 2)
# 
# hosp_flu <- bind_rows(hosp_flu1, hosp_flu2) %>%
#   mutate(pcr_target = "fluav") %>%
#   rename(utility = Region.Name,
#          admissions = Hospitalized.Count) %>%
#   filter(utility %in% sentinel$utility) %>%
#   select(utility, pcr_target, date, admissions, dataset) %>%
#   group_by(utility, pcr_target, date) %>%
#   # identify any duplicate observations (resulting from overlapping dates in the two datasets)
#   add_count() %>%
#   group_by(utility, pcr_target, date) %>%
#   slice(1)
# 
# hosp_rsv1 <- read.csv("DataRaw/2026/RSV_agg.csv") %>%
#   filter(Region.Level == "Sewershed") %>%
#   group_by(Region.Name, Event.Onset.Date) %>%
#   mutate(Hospitalized.Count = sum(Hospitalized.Count)) %>%
#   slice(1) %>%
#   select(-Age.Group) %>%
#   mutate(date = as.Date(Event.Onset.Date, "%m/%d/%Y"),
#          dataset = 1)
# 
# hosp_rsv2 <- read.csv("DataRaw/CSPH Data _ April 2026/RSV_agg_April_update.csv") %>%
#   filter(Region.Level == "Sewershed") %>%
#   group_by(Region.Name, Event.Onset.Date) %>%
#   mutate(Hospitalized.Count = sum(Hospitalized.Count)) %>%
#   slice(1) %>%
#   select(-Age.Group) %>%
#   mutate(date = as.Date(Event.Onset.Date, "%m/%d/%Y"),
#          dataset = 2)
# 
# hosp_rsv <- bind_rows(hosp_rsv1, hosp_rsv2) %>%
#   ungroup() %>%
#   mutate(pcr_target = "rsv") %>%
#   rename(utility = Region.Name,
#          admissions = Hospitalized.Count) %>%
#   filter(utility %in% sentinel$utility) %>%
#   select(utility, pcr_target, date, admissions, dataset) %>%
#   group_by(utility, pcr_target, date) %>%
#   # identify any duplicate observations (resulting from overlapping dates in the two datasets)
#   add_count() %>%
#   group_by(utility, pcr_target, date) %>%
#   slice(1)

# mc_gt10 <- mc_combined %>%
#   filter(vif > 10)
# 
# ww_hosp_model_exclude_high_vif <- ww_hosp_model_exclude %>%
#   semi_join(mc_gt10)
# 
# # fit models on the most problematic sewersheds with VIF > 10
# fit_models <- function(df) {
#   
#   detect_only <- glmmTMB(admissions_weekly ~ detect_weekly + admissions_weekly_w_1,
#                          family = nbinom2, data = df)
#   
#   conc_only <- glmmTMB(admissions_weekly ~ log_raw_conc_weekly + admissions_weekly_w_1,
#                        family = nbinom2, data = df)
#   
#   both <- glmmTMB(admissions_weekly ~ log_raw_conc_weekly + detect_weekly + admissions_weekly_w_1,
#                   family = nbinom2, data = df)
#   
#   list(detect_only = detect_only,
#        conc_only = conc_only,
#        both = both)
# }
# 
# grouped_data <- ww_hosp_model_exclude_high_vif %>%
#   group_by(utility, pcr_target) %>%
#   group_nest()
# 
# models_df <- ww_hosp_model_exclude_high_vif %>%
#   group_by(utility, pcr_target) %>%
#   group_nest() %>%
#   mutate(models = map(data, fit_models))
# 
# aic_long <- models_df %>%
#   mutate(
#     aic = map(models, ~ tibble(
#       model = c("detect_only", "conc_only", "both"),
#       aic = AIC(.x$detect_only, .x$conc_only, .x$both)
#     ))
#   ) %>%
#   select(utility, pcr_target, aic) %>%
#   unnest(aic)
# 
# lrt_df <- models_df %>%
#   mutate(
#     lrt_detect = map(models, ~ anova(.x$detect_only, .x$both, test="LRT")),
#     lrt_conc   = map(models, ~ anova(.x$conc_only, .x$both, test="LRT"))
#   ) %>%
#   select(utility, pcr_target, lrt_detect, lrt_conc)
# 
# extract_lrt <- function(x, prefix) {
#   tibble(
#     !!paste0(prefix, "_chisq") := x$Chisq[2],
#     !!paste0(prefix, "_df")    := x$`Chi Df`[2],
#     !!paste0(prefix, "_p")     := x$`Pr(>Chisq)`[2]
#   )
# }
# 
# lrt_clean <- lrt_df %>%
#   mutate(
#     detect = map(lrt_detect, ~ extract_lrt(.x, "detect")),
#     conc   = map(lrt_conc, ~ extract_lrt(.x, "conc"))
#   ) %>%
#   select(utility, pcr_target, detect, conc) %>%
#   unnest(c(detect, conc))
# 
# ##########################################################################################
# 
