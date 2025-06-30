#############################################################################
##### CHARACTERIZING DYNAMIC POPULATIONS IN WASTEWATER MONITORING AREAS #####
############################## SUMMER 2025 ##################################
################################## EJW ######################################

# ANALYSIS 3B. ROLLING CORRELATIONS BETWEEN DEVICE COUNTS AND FLOW RATES

# set working directory
setwd("/Users/emwu9912/Documents/CU Anschutz/COVID Wastewater/")

# load libraries
library(pacman)
p_load(gridExtra, NonParRolCor, tidyverse)

load("DataProcessed/flow_mobile_corr_daily.Rda")

# fill daily timestep so there are no missing values
flow_mobile_corr_daily_fill <- flow_mobile_corr_daily %>%
  group_by(sewershed_analysis_name) %>%
  # last observation carried forward
  fill(c(flow_rate), .direction = "down") %>%
  # create 7-day moving average for device count and flow
  mutate(total_devices_daily_7daymovingavg = zoo::rollmean(total_devices_daily, k = 7, fill = NA, align = "right"),
         flow_rate_7daymovingavg = zoo::rollmean(flow_rate, k = 7, fill = NA, align = "right")) %>%
  ungroup() %>%
  drop_na()

save(flow_mobile_corr_daily_fill, file = "DataProcessed/flow_mobile_corr_daily_fill.Rda")

# check how many sewersheds remain (which also creates the vector for the for-loop)
corr_sewersheds <- unique(unlist(flow_mobile_corr_daily_fill$sewershed_analysis_name))
corr_sewersheds

# get sample sizes
flow_mobile_corr_daily_fill %>%
  group_by(sewershed_analysis_name) %>%
  slice(1) %>%
  group_by(cluster) %>%
  count()

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

# calculate rolling correlations for all sewersheds
corr_estim <- list()
corr_list <- mget(ls(pattern = "input_"))
for (i in seq_along(corr_list)){
  corr_estim[[i]] <- rolcor_estim_1win(corr_list[[i]], CorMethod = "spearman", widthwin = 90,
                                       Align = "right",rmltrd = T, Scale = T, MCSim = 1000,
                                       Np = 2, prob = 0.95)
}

corr_df <- list()
for(i in 1:length(corr_estim)){
  corr_df[[i]] <- as.data.frame(corr_estim[[i]]$Correlation_coefficients) %>%
    mutate(start_day_index = seq(0, nrow(corr_estim[[i]]$Correlation_coefficients)-1, 1),
           end_date = min(corrs[[i]]$date) + 89 + start_day_index,
           sewershed_analysis_name = corr_sewersheds[[i]])
}

# save correlation results
saveRDS(corr_df, "DataProcessed/corr_df.rds")

########################################################################################################
