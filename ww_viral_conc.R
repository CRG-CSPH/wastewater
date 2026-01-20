# set working directory
setwd("/Users/emwu9912/Documents/CU Anschutz/Wastewater/")

# turn off scientific notation
options(scipen = 999)

library(pacman)
p_load(DHARMa, FSA, glmmTMB, lme4, lubridate, MASS, performance, reformulas, slider, tidyverse, TMB)

theme <- theme(plot.title = element_text(hjust = 0.5),
               axis.title.x = element_text(margin = margin(t = 5)),
               axis.ticks.y = element_blank(),
               axis.text.x = element_text(angle = 90, hjust = 0.95, vjust = 0.4))

sentinel <- read.csv("DataRaw/sentinel.csv")
sewershed_names <- read.csv("DataRaw/sewershed_names.csv")

visit <- read.csv("DataProcessed/2025_10_21_daily_visits.csv") %>%
  rename(cdphe_flow_name = DESTINATION_AREA_SEWERSHED,
         total_devices_daily = STOPS_BY_DAY_L) %>%
  mutate(date = as.Date(DAY)) %>%
  left_join(sewershed_names) %>%
  rename(utility = sewershed_sentinel_name) %>%
  filter(#date >= "2024-09-30",
         utility %in% sentinel$utility) %>%
  select(utility, date, total_devices_daily)

# recover data from utilities with missing or zero flow
flow_missing_impute <- read.csv("DataRaw/20251106_MissingFlowRates_CSPH.csv") %>%
  mutate(date = as.Date(sample_collect_date, "%m/%d/%y")) %>%
  rename(flow_rate_impute = flow_rate) %>%
  select(utility, pcr_target, date, flow_rate_impute)

ww <- read.csv("DataRaw/2_CDPHE_Wastewater_Data_ 2025-08-04 .csv") %>%
  # cleaning and formatting
  mutate(utility = wwtp_name,
         wwtp_name = tolower(wwtp_name),
         pcr_target = tolower(pcr_target),
         date = as.Date(sample_collect_date, "%Y-%m-%d"),
         test_result_date = as.Date(test_result_date, "%Y-%m-%d")) %>%
  filter(utility %in% sentinel$utility,
         pcr_target %in% c("fluav", "rsv", "sars-cov-2")) %>%
         #date >= "2024-09-30") %>%
  left_join(flow_missing_impute) %>%
  # impute anything BLOD with half the LOD
  mutate(pcr_target_avg_conc = case_when(pcr_target_avg_conc <= lod_sewage ~ lod_sewage/2,
                                         TRUE ~ pcr_target_avg_conc),
         # impute missing and zero flow rates with recovered data
         flow_rate = case_when(flow_rate == 0 ~ flow_rate_impute,
  # set remaining zero flow rates to missing
  TRUE ~ flow_rate)) %>%
  select(-flow_rate_impute) %>%
  group_by(utility, pcr_target) %>%
  # impute the remaining missing flow rates with the average of the last two non-missing flow rates
  mutate(flow_rate_last2avg = slide_dbl(flow_rate,
                                        ~ mean(tail(na.omit(.x), 2)),
                                        .before = Inf,
                                        .complete = F),
         flow_rate = ifelse(is.na(flow_rate), flow_rate_last2avg, flow_rate)) %>%
  select(-flow_rate_last2avg) %>%
  left_join(sentinel)

ww_remove_outliers <- ww %>%
  group_by(utility, pcr_target) %>%
  arrange(date) %>%
  mutate(roll_mean = cummean(pcr_target_avg_conc),
         n = row_number(),
         roll_var = cumsum((pcr_target_avg_conc - roll_mean)^2) / pmax(n - 1, 1),
         roll_sd = sqrt(roll_var),
         z = (pcr_target_avg_conc - roll_mean) / roll_sd) %>%
  filter(abs(z) <= 4 | is.na(z))

write.csv(ww_remove_outliers, "DataProcessed/ww_remove_outliers.csv")

desc <- ww_remove_outliers %>%
  rename(capacity_based_population = population_served) %>%
  select(lab_phase, wwtp_name, utility, category, capacity_based_population, census_population, pcr_target, sample_type, pcr_target_avg_conc, flow_rate, date, test_result_date, category) %>%
  left_join(visit) %>%
  mutate(flow_rate_per_capita = flow_rate*1000000/census_population,
         device_count_per_capita = total_devices_daily/census_population,
         ndays = as.numeric(difftime(test_result_date, date, units = "days")),
         week_end = ceiling_date(date, "weeks", week_start = getOption("lubridate.week.start", 6))) %>%
  group_by(utility, pcr_target) %>%
  mutate(mean_ndays = mean(ndays)) %>%
  filter(date >= "2024-09-30") %>%
  ungroup() %>%
  arrange(desc(category), desc(utility)) %>%
  mutate(utility_factor = factor(utility, unique(utility))) %>%
  group_by(utility, pcr_target, week_end) %>%
  add_count()

desc_covid <- desc %>%
  filter(pcr_target == "sars-cov-2")

boxplot_flow <- ggplot(data = desc_covid) +
  geom_boxplot(aes(x = flow_rate_per_capita, y = utility_factor, group = utility, fill = category)) +
  labs(x = "Daily Flow Rate Per Capita (gal)", y = NULL, fill = "Category") +
  theme_bw() +
  theme(panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank()); boxplot_flow

ggsave("Figures/boxplot_flow.png", height = 8, width = 14, plot = boxplot_flow)

boxplot_mobile <- ggplot(data = desc_covid) +
  geom_boxplot(aes(x = device_count_per_capita, y = utility_factor, group = utility, fill = category)) +
  labs(x = "Daily Device Count per Capita", y = NULL, fill = "Category") +
  theme_bw() +
  theme(panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank()); boxplot_mobile

ggsave("Figures/boxplot_mobile.png", height = 8, width = 14, plot = boxplot_mobile)

# Kruskal-Wallis test to determine if there are significant differences in
# flow rate per capita or device count per capita between sewershed categories
kruskal.test(flow_rate_per_capita ~ category, data = desc_covid) # p < 0.001 (significant differences)
kruskal.test(device_count_per_capita ~ category, data = desc_covid) # p = 0.20 (no significant differences)

# Dunn's test for pairwise comparisons (flow)
dunnTest(flow_rate_per_capita ~ category, data = desc_covid,
         method = "bonferroni")
# High Tourism - Metro: p < 0.001
# High Tourism - Other: p < 0.001
# Metro - Other: p = 0.019
# all significantly different from one another

# sampling frequency dataset
samp_freq <- desc %>% slice(1)

titles <- c("SARS-CoV-2", "Flu A", "RSV")
pcr_targets <- c("sars-cov-2", "fluav", "rsv")
sample_type_plot <- list()
for(i in 1:length(titles)){
  sample_type_plot[[i]] <- ggplot(data = desc %>% filter(pcr_target == pcr_targets[i])) +
                                  #aes(x = date, y = flow_rate_per_capita, color = sample_type)) +
    geom_point(aes(x = date, y = flow_rate_per_capita, color = sample_type)) +
    geom_line(aes(x = date, y = flow_rate_per_capita)) +
    labs(x = NULL, y = "Daily Flow Rate per Capita (gal)", color = "Sample Type") +
    ggtitle(paste("Flow Rate by Sampling Technique,", titles[i])) +
    facet_wrap(~utility, scales = "free") +
    theme +
    theme(axis.text.x = element_text(angle = 0))
  ggsave(paste0("Figures/sample_type_plot_", pcr_targets[i], ".png"), width = 17, height = 9, plot = sample_type_plot[[i]])
}

sample_freq_plot <- list()
for(i in 1:length(titles)){
  sample_freq_plot[[i]] <- ggplot(data = samp_freq %>% filter(pcr_target == pcr_targets[i])) +
    geom_bar(aes(x = week_end, y = n), stat = "identity") +
    labs(x = NULL, y = "Number of Times Sampled Per Week") +
    scale_y_continuous(limits = c(0, 3), breaks = seq(0, 3, 1)) +
    ggtitle(paste("Sampling Frequency,", titles[i])) +
    facet_wrap(~utility, scales = "free") +
    theme +
    theme(axis.text.x = element_text(angle = 0))
  ggsave(paste0("Figures/sample_freq_plot_", pcr_targets[i], ".png"), width = 15, height = 9, plot = sample_freq_plot[[i]])
}

sample_ndays_plot <- list()
for(i in 1:length(titles)){
  sample_ndays_plot[[i]] <- ggplot(data = desc %>% filter(pcr_target == pcr_targets[i])) +
    geom_bar(aes(x = date, y = ndays), stat = "identity") +
    geom_line(aes(x = date, y = mean_ndays), color = "red") +
    labs(x = NULL, y = "Number of Days") +
    scale_y_continuous(limits = c(0, 21), breaks = seq(0, 21, 3)) +
    ggtitle(paste("Number of Days Between Sample Collection and Test Result,", titles[i])) +
    facet_wrap(~utility, scales = "free") +
    theme +
    theme(axis.text.x = element_text(angle = 0))
  ggsave(paste0("Figures/sample_ndays_plot_", pcr_targets[i], ".png"), width = 15, height = 9, plot = sample_ndays_plot[[i]])
}

targets <- c(COVID = "sars-cov-2",
             FLUA = "fluav",
             RSV = "rsv")

wval <- imap_dfr(targets , ~ {
  read.csv(glue::glue("DataProcessed/2025_12_04_dataset/2025_12_04_{.y}_wval.csv")) %>%
    mutate(pcr_target = .x,
           week_end = as.Date(week_end, "%Y-%m-%d")) %>%
    rename(WVAL_conc_weekly = WVAL) %>%
    select(week_end, wwtp_name, pcr_target, WVAL_conc_weekly)
})

ww_weekly <- desc %>%
  group_by(utility, pcr_target) %>%
  rename(raw_conc = pcr_target_avg_conc) %>%
  mutate(flow_norm_conc = raw_conc*flow_rate,
         mobile_norm_conc = raw_conc/total_devices_daily,
         combo_norm_conc = raw_conc*(flow_rate/total_devices_daily)) %>%
  # aggregate to weekly
  group_by(utility, pcr_target, week_end) %>%
  mutate(raw_conc_weekly = median(raw_conc),
         flow_norm_conc_weekly = median(flow_norm_conc),
         mobile_norm_conc_weekly = median(mobile_norm_conc),
         combo_norm_conc_weekly = median(combo_norm_conc)) %>%
  slice(1) %>%
  left_join(wval)

hosp_covid <- read.csv("DataRaw/CSPH/20251006_COVID_agg.csv") %>%
  filter(Region.Level == "Sewershed")

hosp_flu <- read.csv("DataRaw/CSPH/20250926_Flu_agg.csv") %>%
  filter(Region.Level == "Sewershed",
         Pathogen.Name == "INFLUENZA A")

hosp_rsv <- read.csv("DataRaw/CSPH/20250922_RSV_agg.csv") %>%
  filter(Region.Level == "Sewershed") %>%
  group_by(Region.Name, Event.Onset.Date) %>%
  mutate(Hospitalized.Count = sum(Hospitalized.Count)) %>%
  slice(1) %>%
  select(-Age.Group)

hosp <- bind_rows(hosp_covid, hosp_flu, hosp_rsv) %>%
  mutate(date = as.Date(Event.Onset.Date, "%m/%d/%Y"),
         pcr_target = case_when(Pathogen.Name == "INFLUENZA A" ~ "fluav",
                                Pathogen.Name == "RESPIRATORY SYNCYTIAL VIRUS" ~ "rsv",
                                TRUE ~ "sars-cov-2"),
         week_end = ceiling_date(date, "weeks", week_start = getOption("lubridate.week.start", 6))) %>%
  rename(utility = Region.Name,
         admissions = Hospitalized.Count) %>%
  select(utility, pcr_target, date, week_end, admissions) %>%
  filter(date >= "2024-09-30")

hosp_weekly <- hosp %>%
  mutate(week_end = ceiling_date(date, "weeks", week_start = getOption("lubridate.week.start", 6))) %>%
  group_by(utility, pcr_target, week_end) %>%
  mutate(admissions_weekly = sum(admissions)) %>%
  slice(1) %>%
  select(-c(admissions, date))

ww_hosp_weekly <- merge(ww_weekly, hosp_weekly, c("pcr_target", "utility", "week_end"), all = T) %>%
  filter(utility %in% sentinel$utility) %>%
  mutate(admissions_weekly = replace_na(admissions_weekly, 0)) %>%
  group_by(utility, pcr_target) %>%
  fill(category, census_population, raw_conc_weekly, flow_norm_conc_weekly,
       mobile_norm_conc_weekly, combo_norm_conc_weekly, WVAL_conc_weekly) %>%
  # create lagged variables
  mutate(admissions_weekly_w_1 = lag(admissions_weekly, 1)) %>%
  mutate(across(c(raw_conc_weekly, flow_norm_conc_weekly, mobile_norm_conc_weekly,
                  combo_norm_conc_weekly, WVAL_conc_weekly),
                list(w_1 = ~lag(.x, 1)),
                .names = "{.col}_{.fn}")) %>%
  select(utility, category, census_population, pcr_target, date, contains("week"))

write.csv(ww_hosp_weekly, "DataProcessed/ww_hosp_weekly.csv")

scatter <- list()
for(i in 1:length(titles)){
  scatter[[i]] <- ggplot(data = ww_hosp_weekly %>% filter(pcr_target == pcr_targets[i]),
                         aes(x = raw_conc_weekly, y = admissions_weekly)) +
    geom_point() +
    geom_smooth(method = "glm.nb", se = TRUE, color = "blue") +
    labs(x = "Weekly Raw Concentration (copies/L)", y = "Weekly Total Hospital Admissions (n)") +
    ggtitle(paste(titles[i], "Hospital Admissions vs. Raw Viral Concentrations")) +
    scale_x_continuous(labels = scales::comma) +
    facet_wrap(~utility, scales = "free") +
    theme +
    theme(axis.text.x = element_text(angle = 0))
  ggsave(paste0("Figures/scatter_", pcr_targets[i], ".png"), width = 15, height = 9, plot = scatter[[i]])
}

ww_hosp_weekly_scale <- ww_hosp_weekly %>%
  group_by(utility, pcr_target) %>%
  mutate(across(contains("conc"), scale)) %>%
  select(-date) %>%
  drop_na()

write.csv(ww_hosp_weekly_scale, "DataProcessed/ww_hosp_weekly_scale.csv")

ww_hosp_weekly_scale_long <- ww_hosp_weekly_scale %>%
  select(!contains("_w_"), -c(category, census_population, admissions_weekly)) %>%
  pivot_longer(c(contains("weekly")), names_to = "norm_method", values_to = "value")

plot_norm_methods <- list()
for(i in 1:length(titles)){
  plot_norm_methods[[i]] <- ggplot(data = ww_hosp_weekly_scale_long %>% filter(pcr_target == pcr_targets[i])) +
    geom_line(aes(x = week_end, y = value, group = norm_method, color = norm_method)) +
    labs(x = "Week End Date", y = "Z-Scored Value", color = "Normalization Method") +
    ggtitle(paste("Normalized", titles[i], "Viral Concentrations")) +
    scale_x_date(date_labels = "%b %Y") +
    scale_color_manual(labels = c("Combination-Normalized", "Flow-Normalized",
                                  "Mobile Device-Normalized", "Raw Concentration", "WVAL"),
                       values = c("violet", "turquoise2", "firebrick",
                                  "darkorchid4", "darkorange")) +
    guides(color = guide_legend(override.aes = list(linewidth = 8), nrow = 5)) +
    facet_wrap(~utility, scales = "free_y") +
    theme
  ggsave(paste0("Figures/plot_norm_methods_", pcr_targets[i], ".png"), width = 17, height = 9, plot = plot_norm_methods[[i]])
}

plot_hosp_ts <- list()
for(i in 1:length(titles)){
  plot_hosp_ts[[i]] <- ggplot(data = ww_hosp_weekly %>% filter(pcr_target == pcr_targets[i])) +
    geom_bar(aes(x = week_end, y = admissions_weekly), stat = "identity") +
    labs(x = "Week End Date", y = "Weekly Hospital Admissions") +
    ggtitle(paste(titles[i], "Hospital Admissions Over Time")) +
    scale_x_date(date_labels = "%b %Y") +
    facet_wrap(~utility, scales = "free_y") +
    theme
  ggsave(paste0("Figures/plot_hosp_ts_", pcr_targets[i], ".png"), width = 17, height = 9, plot = plot_hosp_ts[[i]])
}

# do we need zero inflation?
m_nb <- glmmTMB(admissions_weekly ~ raw_conc_weekly +
                  offset(log(census_population)) + (1 | utility),
                family = nbinom2,
                data = ww_hosp_weekly_scale %>% filter(pcr_target == "sars-cov-2"))

sim <- simulateResiduals(m_nb)
plot(sim)
testZeroInflation(sim) # p = 0.528, no evidence of excess/structural zeros, do not use ZI model

# check ACF of residuals
acf(residuals(m_nb)) # some autocorrelation at 1-8 weeks, include lagged admissions as predictor

# variables to loop over
norm_vars <- c("raw_conc_weekly", "flow_norm_conc_weekly", "mobile_norm_conc_weekly",
               "combo_norm_conc_weekly", "WVAL_conc_weekly")

norm_names <- c("raw", "flow", "mobile", "combo", "wval")

# function to build NB model formula for any normalization method
build_formula <- function(var) {
  as.formula(paste("admissions_weekly ~", var,
                   "+ admissions_weekly_w_1 + offset(log(census_population)) + (1 | utility)"))
}

# function to fit NB model to any PCR target (all categories)
fit_nb_model_all <- function(var, target) {
  glmmTMB(formula = build_formula(var),
          family = nbinom2,
          data = dplyr::filter(ww_hosp_weekly_scale, pcr_target == target))
}

# function to fit NB model to any PCR target (each category)
fit_nb_model_cat <- function(var, target, ww_cat) {
  glmmTMB(formula = build_formula(var),
          family = nbinom2,
          data = dplyr::filter(ww_hosp_weekly_scale, pcr_target == target, category == ww_cat))
}

# run models for each pathogen and category
models_primary_covid_all <- lapply(norm_vars, fit_nb_model_all, target = "sars-cov-2")
models_primary_covid_hightourism <- lapply(norm_vars, fit_nb_model_cat, target = "sars-cov-2", ww_cat = "High Tourism")
models_primary_covid_metro <- lapply(norm_vars, fit_nb_model_cat, target = "sars-cov-2", ww_cat = "Metro")
models_primary_covid_other <- lapply(norm_vars, fit_nb_model_cat, target = "sars-cov-2", ww_cat = "Other")

models_primary_flu_all <- lapply(norm_vars, fit_nb_model_all, target = "fluav")
models_primary_flu_hightourism <- lapply(norm_vars, fit_nb_model_cat, target = "fluav", ww_cat = "High Tourism")
models_primary_flu_metro <- lapply(norm_vars, fit_nb_model_cat, target = "fluav", ww_cat = "Metro")
models_primary_flu_other <- lapply(norm_vars, fit_nb_model_cat, target = "fluav", ww_cat = "Other")

models_primary_rsv_all <- lapply(norm_vars, fit_nb_model_all, target = "rsv")
models_primary_rsv_hightourism <- lapply(norm_vars, fit_nb_model_cat, target = "rsv", ww_cat = "High Tourism")
models_primary_rsv_metro <- lapply(norm_vars, fit_nb_model_cat, target = "rsv", ww_cat = "Metro")
models_primary_rsv_other <- lapply(norm_vars, fit_nb_model_cat, target = "rsv", ww_cat = "Other")

# function to generate model performance metrics
compare_models_all <- function(model_list) {
  data.frame(normalization = norm_names,
             AIC = sapply(model_list, AIC),
             beta = sapply(model_list, \(m) summary(m)$coefficients$cond[2,1]),
             p_value = sapply(model_list, \(m) summary(m)$coefficients$cond[2,4]))
}

# compare models
model_comparisons <- mget(ls(pattern = "^models_primary_")) %>%
  map(~ compare_models_all(.x) %>% arrange(AIC))

all_models_tbl_primary <- imap_dfr(model_comparisons, ~ {
  .x %>% 
    mutate(
      model_set = .y)
})

all_models_tbl_primary <- all_models_tbl_primary %>%
  separate(model_set, into = c("model", "analysis", "pathogen", "category"), sep = "_", remove = FALSE) %>%
  select(analysis, pathogen, category, normalization, AIC, beta, p_value) %>%
  mutate(irr = round(exp(beta), 2),
         p_value = round(p_value, 3),
         significant = ifelse(p_value < 0.05, "yes", "no"))

write.csv(all_models_tbl_primary, "DataProcessed/model_comparisons_all_primary.csv", row.names = FALSE)

all_models_irr_primary <- all_models_tbl_primary %>%
  arrange(pathogen, category, match(normalization, c("raw", "flow", "mobile", "combo", "wval"))) %>%
  select(pathogen, category, normalization, irr) %>%
  pivot_wider(names_from = normalization, values_from = irr)

write.csv(all_models_irr_primary, "DataProcessed/all_models_irr_primary.csv")

all_models_aic_primary <- all_models_tbl_primary %>%
  group_by(pathogen, category) %>%
  mutate(rank = ifelse(is.na(AIC), NA, row_number()),
         delta_AIC_next = AIC - lag(AIC),
         delta_AIC_cumul = AIC - min(AIC),
         AIC_weight = round(exp(-0.5*delta_AIC_cumul) / sum(exp(-0.5*delta_AIC_cumul)), 3)) %>%
  select(pathogen, category, normalization, contains("AIC"), rank)

write.csv(all_models_aic_primary, "DataProcessed/all_models_aic_primary.csv")

all_models_rank_primary <- all_models_aic_primary %>%
  arrange(pathogen, category, match(normalization, c("raw", "flow", "mobile", "combo", "wval"))) %>%
  select(pathogen, category, normalization, rank) %>%
  pivot_wider(names_from = normalization, values_from = rank)

write.csv(all_models_rank_primary, "DataProcessed/all_models_rank_primary.csv")

# SENSITIVITY ANALYSIS
# variables to loop over
norm_vars_sens <- c("raw_conc_weekly_w_1", "flow_norm_conc_weekly_w_1", "mobile_norm_conc_weekly_w_1",
               "combo_norm_conc_weekly_w_1", "WVAL_conc_weekly_w_1")

norm_names_sens <- c("raw_w_1", "flow_w_1", "mobile_w_1", "combo_w_1", "wval_w_1")

# run models for each pathogen and category (sensitivity analysis)
models_sensitivity_covid_all <- lapply(norm_vars_sens, fit_nb_model_all, target = "sars-cov-2")
models_sensitivity_covid_hightourism <- lapply(norm_vars_sens, fit_nb_model_cat, target = "sars-cov-2", ww_cat = "High Tourism")
models_sensitivity_covid_metro <- lapply(norm_vars_sens, fit_nb_model_cat, target = "sars-cov-2", ww_cat = "Metro")
models_sensitivity_covid_other <- lapply(norm_vars_sens, fit_nb_model_cat, target = "sars-cov-2", ww_cat = "Other")

models_sensitivity_flu_all <- lapply(norm_vars_sens, fit_nb_model_all, target = "fluav")
models_sensitivity_flu_hightourism <- lapply(norm_vars_sens, fit_nb_model_cat, target = "fluav", ww_cat = "High Tourism")
models_sensitivity_flu_metro <- lapply(norm_vars_sens, fit_nb_model_cat, target = "fluav", ww_cat = "Metro")
models_sensitivity_flu_other <- lapply(norm_vars_sens, fit_nb_model_cat, target = "fluav", ww_cat = "Other")

models_sensitivity_rsv_all <- lapply(norm_vars_sens, fit_nb_model_all, target = "rsv")
models_sensitivity_rsv_hightourism <- lapply(norm_vars_sens, fit_nb_model_cat, target = "rsv", ww_cat = "High Tourism")
models_sensitivity_rsv_metro <- lapply(norm_vars_sens, fit_nb_model_cat, target = "rsv", ww_cat = "Metro")
models_sensitivity_rsv_other <- lapply(norm_vars_sens, fit_nb_model_cat, target = "rsv", ww_cat = "Other")

# compare models
model_comparisons_sens <- mget(ls(pattern = "^models_sensitivity_")) %>%
  map(~ compare_models_all(.x) %>% arrange(AIC))

all_models_tbl_sensitivity <- imap_dfr(model_comparisons_sens, ~ {
  .x %>% 
    mutate(
      model_set = .y)
})

all_models_tbl_sensitivity <- all_models_tbl_sensitivity %>%
  separate(model_set, into = c("model", "analysis", "pathogen", "category"), sep = "_", remove = FALSE) %>%
  select(analysis, pathogen, category, normalization, AIC, beta, p_value) %>%
  mutate(irr = round(exp(beta), 2),
         p_value = round(p_value, 3),
         significant = ifelse(p_value < 0.05, "yes", "no"))

write.csv(all_models_tbl_sensitivity, "DataProcessed/model_comparisons_all_sensitivity.csv", row.names = FALSE)

all_models_irr_sensitivity <- all_models_tbl_sensitivity %>%
  arrange(pathogen, category, match(normalization, c("raw", "flow", "mobile", "combo", "wval"))) %>%
  select(pathogen, category, normalization, irr) %>%
  pivot_wider(names_from = normalization, values_from = irr)

write.csv(all_models_irr_sensitivity, "DataProcessed/all_models_irr_sensitivity.csv")

all_models_aic_sensitivity <- all_models_tbl_sensitivity %>%
  group_by(pathogen, category) %>%
  mutate(rank = ifelse(is.na(AIC), NA, row_number()),
         delta_AIC_next = AIC - lag(AIC),
         delta_AIC_cumul = AIC - min(AIC),
         AIC_weight = round(exp(-0.5*delta_AIC_cumul) / sum(exp(-0.5*delta_AIC_cumul)), 3)) %>%
  select(pathogen, category, normalization, contains("AIC"), rank)

write.csv(all_models_aic_sensitivity, "DataProcessed/all_models_aic_sensitivity.csv")

all_models_rank_sensitivity <- all_models_aic_sensitivity %>%
  arrange(pathogen, category, match(normalization, c("raw", "flow", "mobile", "combo", "wval"))) %>%
  select(pathogen, category, normalization, rank) %>%
  pivot_wider(names_from = normalization, values_from = rank)

write.csv(all_models_rank_sensitivity, "DataProcessed/all_models_rank_sensitivity.csv")

##########################################################################################

