###################################
### Code for: 
### Määttä et al. (2026) A cross-site comparison of ecosystem- and plot-scale methane fluxes across multiple sites
### Code created by Tiia Määttä, with parts written with GPT 4, 4o and 5.4

##################################
###----STATISTICAL ANALYSES----### 
##################################

# Remove all variables
rm(list=ls(all = TRUE))
# Clear the console window
cat("\014")
# Close all graphics windows
graphics.off()

# open libraries

library("ggplot2")
library("dplyr")
library("lubridate")
library("viridis")
library("ggpubr")
library("moments")
library("tidyverse")
library("stringr")
library("ggpmisc")
library("car")
library("MASS")
library("mgcv")
library("nlme")
library("patchwork")
library("scales")
library("see")
library("ggokabeito")
library("conover.test")
library("multcompView")
library("rcompanion")
library("colorspace")
library("corrplot")
library("bestNormalize")
library("MuMIn")

# IMPORT DATA SETS #

# ALL SITE DATA SETS (incl. predictors)

# available for download at https://doi.org/10.5281/zenodo.17312404

# half-hourly aggregation
all_hh_aggr <- read.csv("path/halfhourly_aggregation.csv")

# hourly aggregation
all_hr_aggr <- read.csv("path/hourly_aggregation.csv")

# daily aggregation
all_d_aggr <- read.csv("path/daily_aggregation.csv")

# weekly aggregation
all_week_aggr <- read.csv("path/weekly_aggregation.csv")

# monthly aggregation
all_month_aggr <- read.csv("path/monthly_aggregation.csv")

# annual aggregation
all_yr_aggr <- read.csv("path/monthly_aggregation.csv")

# the following datasets were not published due to the raw chamber data, some of which have not been published by the data providers. 
# Chamber datasets can be requested from data providers.

# all sites (raw data, i.e., chamber and EC data have not been aggregated)

# with duplicated EC
all_sites_pred_dupl <- read.csv("path/allsites_raw_with_duplicates.csv")
 
# without duplicated EC
all_sites_pred_no_dupl <- read.csv("path/allsites_raw_no_duplicates.csv")

################################################################################################################################

# calculate delta_FCH4 (difference between ecosystem-scale FCH4 and plot-scale FCH4) based on EC and chamber FCH4 medians

all_d_aggr$DELTA_FCH4 <- all_d_aggr$EC_FCH4_MEDIAN - all_d_aggr$CH_FCH4_MEDIAN
all_hh_aggr$DELTA_FCH4 <- all_hh_aggr$EC_FCH4 - all_hh_aggr$CH_FCH4_MEDIAN

all_sites_pred_dupl$DELTA_FCH4 <- all_sites_pred_dupl$EC_FCH4_comb - all_sites_pred_dupl$ch_FCH4_nmolCH4m2s1
all_week_aggr$DELTA_FCH4 <- all_week_aggr$EC_FCH4_MEDIAN - all_week_aggr$CH_FCH4_MEDIAN
all_month_aggr$DELTA_FCH4 <- all_month_aggr$EC_FCH4_MEDIAN - all_month_aggr$CH_FCH4_MEDIAN
all_yr_aggr$DELTA_FCH4 <- all_yr_aggr$EC_FCH4_MEDIAN - all_yr_aggr$CH_FCH4_MEDIAN
all_hr_aggr$DELTA_FCH4 <- all_hr_aggr$EC_FCH4_median - all_hr_aggr$CH_FCH4_MEDIAN
all_sites_pred_dupl$DELTA_FCH4 <- all_sites_pred_dupl$EC_FCH4_comb - all_sites_pred_dupl$ch_FCH4_nmolCH4m2s1

# calculate delta_FCH4 (difference between ecosystem-scale FCH4 and plot-scale FCH4) based on EC and chamber FCH4 means

all_d_aggr$DELTA_FCH4_means <- all_d_aggr$EC_FCH4_MEAN - all_d_aggr$CH_FCH4_MEAN
all_hh_aggr$DELTA_FCH4_means <- all_hh_aggr$EC_FCH4 - all_hh_aggr$CH_FCH4_MEAN
all_week_aggr$DELTA_FCH4_means <- all_week_aggr$EC_FCH4_MEAN - all_week_aggr$CH_FCH4_MEAN
all_month_aggr$DELTA_FCH4_means <- all_month_aggr$EC_FCH4_MEAN - all_month_aggr$CH_FCH4_MEAN
all_yr_aggr$DELTA_FCH4_means <- all_yr_aggr$EC_FCH4_MEAN - all_yr_aggr$CH_FCH4_MEAN
all_hr_aggr$DELTA_FCH4_means <- all_hr_aggr$EC_FCH4_mean - all_hr_aggr$CH_FCH4_MEAN

# check the distributions of the predictors for the d and dlcaggr data

hist(all_sites_pred_no_dupl$EC_LE_F_ANNOPTLM)
hist(all_sites_pred_no_dupl$NEE_F) 
hist(all_sites_pred_no_dupl$EC_PPFD_IN_F)
hist(all_sites_pred_no_dupl$EC_USTAR)
hist(all_sites_pred_no_dupl$EC_WS_F)
hist(all_sites_pred_no_dupl$EC_TA_F)
hist(all_sites_pred_no_dupl$EC_P_F)
hist(all_sites_pred_no_dupl$EC_RECO_DT)
hist(all_sites_pred_no_dupl$EC_GPP_DT)
hist(all_sites_pred_no_dupl$EC_VPD_F)
hist(all_sites_pred_no_dupl$EC_PA_F)
hist(all_d_aggr$U_WIND_F_MEAN) 
hist(all_d_aggr$V_WIND_F_MEAN)



#################################################################################################################################

# BASIC STATISTICS: MEAN AND MEDIAN CH4 FLUX PER SITE, delta FCH4, WILCOXON-MANN-WHITNEY TESTS

### daily

# per site
all_d_summary <- all_d_aggr %>%
  group_by(SITE) %>%
  summarise_at(c("CH_FCH4_MEDIAN", "EC_FCH4_MEDIAN", "DELTA_FCH4"),
               list(mean = ~mean(.x, na.rm = TRUE),
                    sd = ~sd(.x, na.rm = TRUE),
                    median = ~median(.x, na.rm = TRUE),
                    IQR = ~IQR(.x, na.rm = TRUE),
                    CV = ~raster::cv(.x, na.rm = TRUE)))
all_d_summary <- as.data.frame(all_d_summary)
all_d_summary

# across sites
all_d_summary <- all_d_aggr %>%
  summarise_at(c("CH_FCH4_MEDIAN", "EC_FCH4_MEDIAN", "DELTA_FCH4"),
               list(mean = ~mean(.x, na.rm = TRUE),
                    sd = ~sd(.x, na.rm = TRUE),
                    median = ~median(.x, na.rm = TRUE),
                    IQR = ~IQR(.x, na.rm = TRUE),
                    CV = ~raster::cv(.x, na.rm = TRUE)))
all_d_summary <- as.data.frame(all_d_summary)
all_d_summary

# percentage increase in EC in comparison to chamber (CH)

# daily
all_d_aggr <- all_d_aggr %>%
  mutate(
    pct_ECCH_mean   = ifelse(CH_FCH4_MEAN  != 0,
                             100*(EC_FCH4_MEAN   - CH_FCH4_MEAN)  /
                               CH_FCH4_MEAN, NA_real_),
    pct_ECCH_median = ifelse(CH_FCH4_MEDIAN != 0,
                             100*(EC_FCH4_MEDIAN - CH_FCH4_MEDIAN) /
                               CH_FCH4_MEDIAN, NA_real_),
    # Optional symmetric %
    pct_sym_mean   = ifelse(EC_FCH4_MEAN + CH_FCH4_MEAN > 0,
                            200*(EC_FCH4_MEAN - CH_FCH4_MEAN) /
                              (EC_FCH4_MEAN + CH_FCH4_MEAN), NA_real_),
    pct_sym_median = ifelse(EC_FCH4_MEDIAN + CH_FCH4_MEDIAN > 0,
                            200*(EC_FCH4_MEDIAN - CH_FCH4_MEDIAN) /
                              (EC_FCH4_MEDIAN + CH_FCH4_MEDIAN), NA_real_)
  )

# half-hourly
all_hh_aggr <- all_hh_aggr %>%
  mutate(
    pct_ECCH_mean   = ifelse(CH_FCH4_MEAN  != 0,
                             100*(EC_FCH4   - CH_FCH4_MEAN)  /
                               CH_FCH4_MEAN, NA_real_),
    pct_ECCH_median = ifelse(CH_FCH4_MEDIAN != 0,
                             100*(EC_FCH4 - CH_FCH4_MEDIAN) /
                               CH_FCH4_MEDIAN, NA_real_),
    # Optional symmetric %
    pct_sym_mean   = ifelse(EC_FCH4 + CH_FCH4_MEAN > 0,
                            200*(EC_FCH4 - CH_FCH4_MEAN) /
                              (EC_FCH4 + CH_FCH4_MEAN), NA_real_),
    pct_sym_median = ifelse(EC_FCH4 + CH_FCH4_MEDIAN > 0,
                            200*(EC_FCH4 - CH_FCH4_MEDIAN) /
                              (EC_FCH4 + CH_FCH4_MEDIAN), NA_real_)
  )

# hourly
all_hr_aggr <- all_hr_aggr %>%
  mutate(
    pct_ECCH_mean   = ifelse(CH_FCH4_MEAN  != 0,
                             100*(EC_FCH4_mean   - CH_FCH4_MEAN)  /
                               CH_FCH4_MEAN, NA_real_),
    pct_ECCH_median = ifelse(CH_FCH4_MEDIAN != 0,
                             100*(EC_FCH4_median - CH_FCH4_MEDIAN) /
                               CH_FCH4_MEDIAN, NA_real_),
    # Optional symmetric %
    pct_sym_mean   = ifelse(EC_FCH4_mean + CH_FCH4_MEAN > 0,
                            200*(EC_FCH4_mean - CH_FCH4_MEAN) /
                              (EC_FCH4_mean + CH_FCH4_MEAN), NA_real_),
    pct_sym_median = ifelse(EC_FCH4_median + CH_FCH4_MEDIAN > 0,
                            200*(EC_FCH4_median - CH_FCH4_MEDIAN) /
                              (EC_FCH4_median + CH_FCH4_MEDIAN), NA_real_)
  )

# weekly
all_week_aggr <- all_week_aggr %>%
  mutate(
    pct_ECCH_mean   = ifelse(CH_FCH4_MEAN  != 0,
                             100*(EC_FCH4_MEAN   - CH_FCH4_MEAN)  /
                               CH_FCH4_MEAN, NA_real_),
    pct_ECCH_median = ifelse(CH_FCH4_MEDIAN != 0,
                             100*(EC_FCH4_MEDIAN - CH_FCH4_MEDIAN) /
                               CH_FCH4_MEDIAN, NA_real_),
    # Optional symmetric %
    pct_sym_mean   = ifelse(EC_FCH4_MEAN + CH_FCH4_MEAN > 0,
                            200*(EC_FCH4_MEAN - CH_FCH4_MEAN) /
                              (EC_FCH4_MEAN + CH_FCH4_MEAN), NA_real_),
    pct_sym_median = ifelse(EC_FCH4_MEDIAN + CH_FCH4_MEDIAN > 0,
                            200*(EC_FCH4_MEDIAN - CH_FCH4_MEDIAN) /
                              (EC_FCH4_MEDIAN + CH_FCH4_MEDIAN), NA_real_)
  )

# monthly
all_month_aggr <- all_month_aggr %>%
  mutate(
    pct_ECCH_mean   = ifelse(CH_FCH4_MEAN  != 0,
                             100*(EC_FCH4_MEAN   - CH_FCH4_MEAN)  /
                               CH_FCH4_MEAN, NA_real_),
    pct_ECCH_median = ifelse(CH_FCH4_MEDIAN != 0,
                             100*(EC_FCH4_MEDIAN - CH_FCH4_MEDIAN) /
                               CH_FCH4_MEDIAN, NA_real_),
    # Optional symmetric %
    pct_sym_mean   = ifelse(EC_FCH4_MEAN + CH_FCH4_MEAN > 0,
                            200*(EC_FCH4_MEAN - CH_FCH4_MEAN) /
                              (EC_FCH4_MEAN + CH_FCH4_MEAN), NA_real_),
    pct_sym_median = ifelse(EC_FCH4_MEDIAN + CH_FCH4_MEDIAN > 0,
                            200*(EC_FCH4_MEDIAN - CH_FCH4_MEDIAN) /
                              (EC_FCH4_MEDIAN + CH_FCH4_MEDIAN), NA_real_)
  )

# annual
all_yr_aggr <- all_yr_aggr %>%
  mutate(
    pct_ECCH_mean   = ifelse(CH_FCH4_MEAN  != 0,
                             100*(EC_FCH4_MEAN   - CH_FCH4_MEAN)  /
                               CH_FCH4_MEAN, NA_real_),
    pct_ECCH_median = ifelse(CH_FCH4_MEDIAN != 0,
                             100*(EC_FCH4_MEDIAN - CH_FCH4_MEDIAN) /
                               CH_FCH4_MEDIAN, NA_real_),
    # Optional symmetric %
    pct_sym_mean   = ifelse(EC_FCH4_MEAN + CH_FCH4_MEAN > 0,
                            200*(EC_FCH4_MEAN - CH_FCH4_MEAN) /
                              (EC_FCH4_MEAN + CH_FCH4_MEAN), NA_real_),
    pct_sym_median = ifelse(EC_FCH4_MEDIAN + CH_FCH4_MEDIAN > 0,
                            200*(EC_FCH4_MEDIAN - CH_FCH4_MEDIAN) /
                              (EC_FCH4_MEDIAN + CH_FCH4_MEDIAN), NA_real_)
  )


### mean version

# per site
all_d_summary <- all_d_aggr %>%
  group_by(SITE) %>%
  summarise_at(c("CH_FCH4_MEAN", "EC_FCH4_MEAN", "DELTA_FCH4_means"),
               list(mean = ~mean(.x, na.rm = TRUE),
                    sd = ~sd(.x, na.rm = TRUE),
                    median = ~median(.x, na.rm = TRUE),
                    IQR = ~IQR(.x, na.rm = TRUE),
                    CV = ~raster::cv(.x, na.rm = TRUE)))
all_d_summary <- as.data.frame(all_d_summary)
all_d_summary

# all sites together
all_d_summary <- all_d_aggr %>%
  summarise_at(c("CH_FCH4_MEAN", "EC_FCH4_mean", "DELTA_FCH4_means"),
               list(mean = ~mean(.x, na.rm = TRUE),
                    sd = ~sd(.x, na.rm = TRUE),
                    median = ~median(.x, na.rm = TRUE),
                    IQR = ~IQR(.x, na.rm = TRUE),
                    CV = ~raster::cv(.x, na.rm = TRUE)))
all_d_summary <- as.data.frame(all_d_summary)
all_d_summary 

# half hourly

# per site
hh_summary <- all_hh_aggr %>% 
  group_by(SITE) %>%
  summarise_at(c("CH_FCH4_MEDIAN", "EC_FCH4", "DELTA_FCH4"), 
               list(mean = ~mean(.x, na.rm = TRUE),
                    sd = ~sd(.x, na.rm = TRUE),
                    median = ~median(.x, na.rm = TRUE),
                    IQR = ~IQR(.x, na.rm = TRUE),
                    CV = ~raster::cv(.x, na.rm = TRUE)))
hh_summary <- as.data.frame(hh_summary)
hh_summary

# all sites together:
hh_summary <- all_hh_aggr %>% 
  summarise_at(c("CH_FCH4_MEDIAN", "EC_FCH4", "DELTA_FCH4"), 
               list(mean = ~mean(.x, na.rm = TRUE),
                    sd = ~sd(.x, na.rm = TRUE),
                    median = ~median(.x, na.rm = TRUE),
                    IQR = ~IQR(.x, na.rm = TRUE),
                    CV = ~raster::cv(.x, na.rm = TRUE)))
hh_summary <- as.data.frame(hh_summary)
hh_summary

## mean version

# all sites together:
hh_summary <- all_hh_aggr %>% 
  summarise_at(c("CH_FCH4_MEAN", "EC_FCH4", "DELTA_FCH4_means"), 
               list(mean = ~mean(.x, na.rm = TRUE),
                    sd = ~sd(.x, na.rm = TRUE),
                    median = ~median(.x, na.rm = TRUE),
                    IQR = ~IQR(.x, na.rm = TRUE),
                    CV = ~raster::cv(.x, na.rm = TRUE)))
hh_summary <- as.data.frame(hh_summary)
hh_summary

# hourly 
hr_summary <- all_hr_aggr %>% 
  group_by(SITE) %>%
  summarise_at(c("CH_FCH4_MEDIAN", "EC_FCH4_median", "DELTA_FCH4"), 
               list(mean = ~mean(.x, na.rm = TRUE),
                    sd = ~sd(.x, na.rm = TRUE),
                    median = ~median(.x, na.rm = TRUE),
                    IQR = ~IQR(.x, na.rm = TRUE),
                    CV = ~raster::cv(.x, na.rm = TRUE)))
hr_summary <- as.data.frame(hr_summary)
hr_summary

# all sites together: 

hr_summary <- all_hr_aggr %>% 
  summarise_at(c("CH_FCH4_MEDIAN", "EC_FCH4_median", "DELTA_FCH4"), 
               list(mean = ~mean(.x, na.rm = TRUE),
                    sd = ~sd(.x, na.rm = TRUE),
                    median = ~median(.x, na.rm = TRUE),
                    IQR = ~IQR(.x, na.rm = TRUE),
                    CV = ~raster::cv(.x, na.rm = TRUE)))
hr_summary <- as.data.frame(hr_summary)
hr_summary

## mean version

# all sites together: 

hr_summary <- all_hr_aggr %>% 
  summarise_at(c("CH_FCH4_MEAN", "EC_FCH4_mean", "DELTA_FCH4_means"), 
               list(mean = ~mean(.x, na.rm = TRUE),
                    sd = ~sd(.x, na.rm = TRUE),
                    median = ~median(.x, na.rm = TRUE),
                    IQR = ~IQR(.x, na.rm = TRUE),
                    CV = ~raster::cv(.x, na.rm = TRUE)))
hr_summary <- as.data.frame(hr_summary)
hr_summary


# weekly

# per site

week_summary <- all_week_aggr %>% 
  group_by(SITE) %>%
  summarise_at(c("CH_FCH4_MEDIAN", "EC_FCH4_MEDIAN", "DELTA_FCH4"), 
               list(mean = ~mean(.x, na.rm = TRUE),
                    sd = ~sd(.x, na.rm = TRUE),
                    median = ~median(.x, na.rm = TRUE),
                    IQR = ~IQR(.x, na.rm = TRUE),
                    CV = ~raster::cv(.x, na.rm = TRUE)))
week_summary <- as.data.frame(week_summary)
week_summary

# all sites together: 

week_summary <- all_week_aggr %>% 
  summarise_at(c("CH_FCH4_MEDIAN", "EC_FCH4_MEDIAN", "DELTA_FCH4"), 
               list(mean = ~mean(.x, na.rm = TRUE),
                    sd = ~sd(.x, na.rm = TRUE),
                    median = ~median(.x, na.rm = TRUE),
                    IQR = ~IQR(.x, na.rm = TRUE),
                    CV = ~raster::cv(.x, na.rm = TRUE)))

week_summary <- as.data.frame(week_summary)

week_summary

## mean version

# all sites together: 

week_summary <- all_week_aggr %>% 
  summarise_at(c("CH_FCH4_MEAN", "EC_FCH4_MEAN", "DELTA_FCH4_means"), 
               list(mean = ~mean(.x, na.rm = TRUE),
                    sd = ~sd(.x, na.rm = TRUE),
                    median = ~median(.x, na.rm = TRUE),
                    IQR = ~IQR(.x, na.rm = TRUE),
                    CV = ~raster::cv(.x, na.rm = TRUE)))
week_summary <- as.data.frame(week_summary)
week_summary

# monthly

# per site:
month_summary <- all_month_aggr %>% 
  group_by(SITE) %>%
  summarise_at(c("CH_FCH4_MEDIAN", "EC_FCH4_MEDIAN", "DELTA_FCH4"), 
               list(mean = ~mean(.x, na.rm = TRUE),
                    sd = ~sd(.x, na.rm = TRUE),
                    median = ~median(.x, na.rm = TRUE),
                    IQR = ~IQR(.x, na.rm = TRUE),
                    CV = ~raster::cv(.x, na.rm = TRUE)))
month_summary <- as.data.frame(month_summary)
month_summary

# all sites together: 

month_summary <- all_month_aggr %>% 
  summarise_at(c("CH_FCH4_MEDIAN", "EC_FCH4_MEDIAN", "DELTA_FCH4"), 
               list(mean = ~mean(.x, na.rm = TRUE),
                    sd = ~sd(.x, na.rm = TRUE),
                    median = ~median(.x, na.rm = TRUE),
                    IQR = ~IQR(.x, na.rm = TRUE),
                    CV = ~raster::cv(.x, na.rm = TRUE)))
month_summary <- as.data.frame(month_summary)
month_summary

## mean version

# all sites together: 

month_summary <- all_month_aggr %>% 
  summarise_at(c("CH_FCH4_MEAN", "EC_FCH4_MEAN", "DELTA_FCH4_means"), 
               list(mean = ~mean(.x, na.rm = TRUE),
                    sd = ~sd(.x, na.rm = TRUE),
                    median = ~median(.x, na.rm = TRUE),
                    IQR = ~IQR(.x, na.rm = TRUE),
                    CV = ~raster::cv(.x, na.rm = TRUE)))
month_summary <- as.data.frame(month_summary)
month_summary

# annual

# per site:
yr_summary <- all_yr_aggr %>% 
  group_by(SITE) %>%
  summarise_at(c("CH_FCH4_MEDIAN", "EC_FCH4_MEDIAN", "DELTA_FCH4"), 
               list(mean = ~mean(.x, na.rm = TRUE),
                    sd = ~sd(.x, na.rm = TRUE),
                    median = ~median(.x, na.rm = TRUE),
                    IQR = ~IQR(.x, na.rm = TRUE),
                    CV = ~raster::cv(.x, na.rm = TRUE)))
yr_summary <- as.data.frame(yr_summary)
yr_summary

# all sites together: 

yr_summary <- all_yr_aggr %>% 
  summarise_at(c("CH_FCH4_MEDIAN", "EC_FCH4_MEDIAN", "DELTA_FCH4"), 
               list(mean = ~mean(.x, na.rm = TRUE),
                    sd = ~sd(.x, na.rm = TRUE),
                    median = ~median(.x, na.rm = TRUE),
                    IQR = ~IQR(.x, na.rm = TRUE),
                    CV = ~raster::cv(.x, na.rm = TRUE)))
yr_summary <- as.data.frame(yr_summary)
yr_summary

## mean version

# all sites together: 

yr_summary <- all_yr_aggr %>% 
  summarise_at(c("CH_FCH4_MEAN", "EC_FCH4_MEAN", "DELTA_FCH4_means"), 
               list(mean = ~mean(.x, na.rm = TRUE),
                    sd = ~sd(.x, na.rm = TRUE),
                    median = ~median(.x, na.rm = TRUE),
                    IQR = ~IQR(.x, na.rm = TRUE),
                    CV = ~raster::cv(.x, na.rm = TRUE)))
yr_summary <- as.data.frame(yr_summary)
yr_summary


### WILCOXON-MANN-WHITNEY TESTS BETWEEN ECOSYSTEM- AND PLOT-SCALE FCH4

### prepare the data ###

# daily
all_d_ch <- subset(all_d_aggr, select = c(TIMESTAMP, CH_FCH4_MEDIAN))
all_d_ch$method <- "chamber"
all_d_ec <- subset(all_d_aggr, select = c(TIMESTAMP, EC_FCH4_MEDIAN))
all_d_ec$method <- "EC"
colnames(all_d_ch)[2] <- "flux"
colnames(all_d_ec)[2] <- "flux"
all_d <- rbind(all_d_ch, all_d_ec)
all_d <- all_d %>% drop_na(flux)

# mean version

all_d_ch <- subset(all_d_aggr, select = c(TIMESTAMP, CH_FCH4_MEAN))
all_d_ch$method <- "chamber"
all_d_ec <- subset(all_d_aggr, select = c(TIMESTAMP, EC_FCH4_MEAN))
all_d_ec$method <- "EC"
colnames(all_d_ch)[2] <- "flux"
colnames(all_d_ec)[2] <- "flux"
all_d <- rbind(all_d_ch, all_d_ec)
all_d <- all_d %>% drop_na(flux)


# half-hourly
hh_aggr_ch <- subset(all_hh_aggr, select = c(TIMESTAMP_START, CH_FCH4_MEDIAN))
hh_aggr_ch$method <- "chamber"
hh_aggr_ec <- subset(all_hh_aggr, select = c(TIMESTAMP_START, EC_FCH4))
hh_aggr_ec$method <- "EC"
colnames(hh_aggr_ch)[2] <- "flux"
colnames(hh_aggr_ec)[2] <- "flux"
hh_aggr <- rbind(hh_aggr_ch, hh_aggr_ec)
hh_aggr <- hh_aggr %>% drop_na(flux)

# mean version

hh_aggr_ch <- subset(all_hh_aggr, select = c(TIMESTAMP_START, CH_FCH4_MEAN))
hh_aggr_ch$method <- "chamber"
hh_aggr_ec <- subset(all_hh_aggr, select = c(TIMESTAMP_START, EC_FCH4))
hh_aggr_ec$method <- "EC"
colnames(hh_aggr_ch)[2] <- "flux"
colnames(hh_aggr_ec)[2] <- "flux"
hh_aggr <- rbind(hh_aggr_ch, hh_aggr_ec)
hh_aggr <- hh_aggr %>% drop_na(flux)

# hourly
all_hr_ch <- subset(all_hr_aggr, select = c(TIMESTAMP, CH_FCH4_MEDIAN))
all_hr_ch$method <- "chamber"
all_hr_ec <- subset(all_hr_aggr, select = c(TIMESTAMP, EC_FCH4_median))
all_hr_ec$method <- "EC"
colnames(all_hr_ch)[2] <- "flux"
colnames(all_hr_ec)[2] <- "flux"
all_hr <- rbind(all_hr_ch, all_hr_ec)
all_hr <- all_hr %>% drop_na(flux)

# mean version

all_hr_ch <- subset(all_hr_aggr, select = c(TIMESTAMP, CH_FCH4_MEAN))
all_hr_ch$method <- "chamber"
all_hr_ec <- subset(all_hr_aggr, select = c(TIMESTAMP, EC_FCH4_mean))
all_hr_ec$method <- "EC"
colnames(all_hr_ch)[2] <- "flux"
colnames(all_hr_ec)[2] <- "flux"
all_hr <- rbind(all_hr_ch, all_hr_ec)
all_hr <- all_hr %>% drop_na(flux)

# weekly
all_week_ch <- subset(all_week_aggr, select = c(YEAR, Week_of_year, CH_FCH4_MEDIAN))
all_week_ch$method <- "chamber"
all_week_ec <- subset(all_week_aggr, select = c(YEAR, Week_of_year, EC_FCH4_MEDIAN))
all_week_ec$method <- "EC"
colnames(all_week_ch)[3] <- "flux"
colnames(all_week_ec)[3] <- "flux"
all_week <- rbind(all_week_ch, all_week_ec)
all_week <- all_week %>% drop_na(flux)

# mean version

all_week_ch <- subset(all_week_aggr, select = c(YEAR, Week_of_year, CH_FCH4_MEAN))
all_week_ch$method <- "chamber"
all_week_ec <- subset(all_week_aggr, select = c(YEAR, Week_of_year, EC_FCH4_MEAN))
all_week_ec$method <- "EC"
colnames(all_week_ch)[3] <- "flux"
colnames(all_week_ec)[3] <- "flux"
all_week <- rbind(all_week_ch, all_week_ec)
all_week <- all_week %>% drop_na(flux)

# monthly

all_month_ch <- subset(all_month_aggr, select = c(YEAR, MONTH, CH_FCH4_MEDIAN))
all_month_ch$method <- "chamber"
all_month_ec <- subset(all_month_aggr, select = c(YEAR, MONTH, EC_FCH4_MEDIAN))
all_month_ec$method <- "EC"
colnames(all_month_ch)[3] <- "flux"
colnames(all_month_ec)[3] <- "flux"
all_month <- rbind(all_month_ch, all_month_ec)
all_month <- all_month %>% drop_na(flux)

## mean version

all_month_ch <- subset(all_month_aggr, select = c(YEAR, MONTH, CH_FCH4_MEAN))
all_month_ch$method <- "chamber"
all_month_ec <- subset(all_month_aggr, select = c(YEAR, MONTH, EC_FCH4_MEAN))
all_month_ec$method <- "EC"
colnames(all_month_ch)[3] <- "flux"
colnames(all_month_ec)[3] <- "flux"
all_month <- rbind(all_month_ch, all_month_ec)
all_month <- all_month %>% drop_na(flux)

# monthly

all_yr_ch <- subset(all_yr_aggr, select = c(YEAR, CH_FCH4_MEDIAN))
all_yr_ch$method <- "chamber"
all_yr_ec <- subset(all_yr_aggr, select = c(YEAR, EC_FCH4_MEDIAN))
all_yr_ec$method <- "EC"
colnames(all_yr_ch)[2] <- "flux"
colnames(all_yr_ec)[2] <- "flux"
all_yr <- rbind(all_yr_ch, all_yr_ec)
all_yr <- all_yr %>% drop_na(flux)


## mean version

all_yr_ch <- subset(all_yr_aggr, select = c(YEAR, CH_FCH4_MEAN))
all_yr_ch$method <- "chamber"
all_yr_ec <- subset(all_yr_aggr, select = c(YEAR, EC_FCH4_MEAN))
all_yr_ec$method <- "EC"
colnames(all_yr_ch)[2] <- "flux"
colnames(all_yr_ec)[2] <- "flux"
all_yr <- rbind(all_yr_ch, all_yr_ec)
all_yr <- all_yr %>% drop_na(flux)

# daily
wilcox.test(flux ~ method, data = all_d)

# number of observations per EC and chamber
all_d %>%
  filter(!is.na(flux) & method %in% c("EC", "chamber")) %>%
  group_by(method) %>%
  summarise(count = n())

# half-hourly
wilcox.test(flux ~ method, data = hh_aggr)

# number of observations per EC and chamber
hhh_aggr %>%
  filter(!is.na(flux) & method %in% c("EC", "chamber")) %>%
  group_by(method) %>%
  summarise(count = n())

# hourly
wilcox.test(flux ~ method, data = all_hr)

# number of observations per EC and chamber
all_hr %>%
  filter(!is.na(flux) & method %in% c("EC", "chamber")) %>%
  group_by(method) %>%
  summarise(count = n())

# weekly
wilcox.test(flux ~ method, data = all_week)

# number of observations per EC and chamber
all_week %>%
  filter(!is.na(flux) & method %in% c("EC", "chamber")) %>%
  group_by(method) %>%
  summarise(count = n())

# monthly
wilcox.test(flux ~ method, data = all_month) 

# number of observations per EC and chamber
all_month %>%
  filter(!is.na(flux) & method %in% c("EC", "chamber")) %>%
  group_by(method) %>%
  summarise(count = n())

# annual
wilcox.test(flux ~ method, data = all_yr) 

# number of observations per EC and chamber
all_yr %>%
  filter(!is.na(flux) & method %in% c("EC", "chamber")) %>%
  group_by(method) %>%
  summarise(count = n())

### run Wilcoxon test for the difference between ecosystem- and plot-scale FCH4 at each site

# define the function to process each SITE group
process_site_d <- function(df) {
  # Reshape the data to long format
  long_df <- df %>%
    gather(key = "variable", value = "value", EC_FCH4_MEDIAN, CH_FCH4_MEDIAN) %>%
    filter(!is.na(value))
  
  # Perform the Wilcoxon test
  test_result <- if (length(unique(long_df$variable)) == 2) {
    wilcox.test(value ~ variable, data = long_df)
  } else {
    list(p.value = NA, statistic = NA)
  }
  
  # Extract p-value and W value
  p_value <- test_result$p.value
  W_value <- test_result$statistic
  
  # Count non-NA observations
  non_na_counts <- df %>%
    summarise(non_na_EC = sum(!is.na(EC_FCH4_MEDIAN)),
              non_na_CH = sum(!is.na(CH_FCH4_MEDIAN)))
  
  return(data.frame(p_value = p_value, W_value, non_na_counts))
}

# apply the function to each SITE group
results <- all_d_aggr %>%
  group_by(SITE) %>%
  do(process_site_d(.))

results

# hourly

results <- all_hr_aggr %>%
  group_by(SITE) %>%
  do(process_site_d(.))

results

# weekly

results <- all_week_aggr %>%
  group_by(SITE) %>%
  do(process_site_d(.))

results

# monthly

results <- all_month_aggr %>%
  group_by(SITE) %>%
  do(process_site_d(.))

results

# annual

results <- all_yr_aggr %>%
  group_by(SITE) %>%
  do(process_site_d(.))

results

# function for half-hourly aggregation

process_site_auto <- function(df) {
  # Reshape the data to long format
  long_df <- df %>%
    gather(key = "variable", value = "value", EC_FCH4, CH_FCH4_MEDIAN) %>%
    filter(!is.na(value))
  
  # Perform the Wilcoxon test
  test_result <- if (length(unique(long_df$variable)) == 2) {
    wilcox.test(value ~ variable, data = long_df)
  } else {
    list(p.value = NA, statistic = NA)
  }
  
  # Extract p-value and W value
  p_value <- test_result$p.value
  W_value <- test_result$statistic
  
  # Count non-NA observations
  non_na_counts <- df %>%
    summarise(non_na_EC = sum(!is.na(EC_FCH4)),
              non_na_CH = sum(!is.na(CH_FCH4_MEDIAN)))
  
  return(data.frame(p_value = p_value, W_value, non_na_counts))
}

results <- all_hh_aggr %>%
  group_by(SITE) %>%
  do(process_site_auto(.))

results <- as.data.frame(results)
results

# function for hourly aggregation
process_site_hr <- function(df) {
  # Reshape the data to long format
  long_df <- df %>%
    gather(key = "variable", value = "value", EC_FCH4_median, CH_FCH4_MEDIAN) %>%
    filter(!is.na(value))
  
  # Perform the Wilcoxon test
  test_result <- if (length(unique(long_df$variable)) == 2) {
    wilcox.test(value ~ variable, data = long_df)
  } else {
    list(p.value = NA, statistic = NA)
  }
  
  # Extract p-value and W value
  p_value <- test_result$p.value
  W_value <- test_result$statistic
  
  # Count non-NA observations
  non_na_counts <- df %>%
    summarise(non_na_EC = sum(!is.na(EC_FCH4_median)),
              non_na_CH = sum(!is.na(CH_FCH4_MEDIAN)))
  
  return(data.frame(p_value = p_value, W_value, non_na_counts))
}

results <- all_hr_aggr %>%
  group_by(SITE) %>%
  do(process_site_hr(.))

results <- as.data.frame(results)
results

# daily percent change in delta_FCH4:
# overall medians of delta_FCH4
delta_median <- median(all_d_aggr$DELTA_FCH4, na.rm = TRUE)
delta_mean   <- median(all_d_aggr$DELTA_FCH4_means, na.rm = TRUE)

# percent change
((delta_mean - delta_median) / abs(delta_median)) * 100

# calculate percent change of delta_FCH4 for the rest of the aggregations

((median(all_hh_aggr$DELTA_FCH4_means, na.rm = TRUE) - median(all_hh_aggr$DELTA_FCH4, na.rm = TRUE)) / abs(median(all_hh_aggr$DELTA_FCH4, na.rm = T))) * 100

((median(all_hr_aggr$DELTA_FCH4_means, na.rm = TRUE) - median(all_hr_aggr$DELTA_FCH4, na.rm = TRUE)) / abs(median(all_hr_aggr$DELTA_FCH4, na.rm = T))) * 100

((median(all_week_aggr$DELTA_FCH4_means, na.rm = TRUE) - median(all_week_aggr$DELTA_FCH4, na.rm = TRUE)) / abs(median(all_week_aggr$DELTA_FCH4, na.rm = T))) * 100

((median(all_month_aggr$DELTA_FCH4_means, na.rm = TRUE) - median(all_month_aggr$DELTA_FCH4, na.rm = TRUE)) / abs(median(all_month_aggr$DELTA_FCH4, na.rm = T))) * 100

((median(all_yr_aggr$DELTA_FCH4_means, na.rm = TRUE) - median(all_yr_aggr$DELTA_FCH4, na.rm = TRUE)) / abs(median(all_yr_aggr$DELTA_FCH4, na.rm = T))) * 100


################################################################################################################################################################################

### Exploring the contribution of very high CH4 emissions observed per chamber and EC per site on site annual CH4 emission

### calculate 90th percentiles in chamber and EC FCH4 and % of fluxes above that value per aggregation

## Column holding mean-based chamber fluxes
ch_col <- "CH_FCH4_MEAN"

## half-hourly
x_hh <- all_hh_aggr[["EC_FCH4"]]
x_hh <- x_hh[is.finite(x_hh) & x_hh > 0]   # keep only positive fluxes
n_hh <- length(x_hh); total_hh <- sum(x_hh)
cut90_hh <- as.numeric(quantile(x_hh, 0.90, na.rm = TRUE, type = 7))
cut95_hh <- as.numeric(quantile(x_hh, 0.95, na.rm = TRUE, type = 7))
hot90_hh <- sum(x_hh[x_hh > cut90_hh]); hot95_hh <- sum(x_hh[x_hh > cut95_hh])
pctval90_hh <- 100 * mean(x_hh > cut90_hh); pctval95_hh <- 100 * mean(x_hh > cut95_hh)
pctflux90_hh <- 100 * hot90_hh / total_hh; pctflux95_hh <- 100 * hot95_hh / total_hh
df_hh <- rbind(
  data.frame(aggregation="half-hourly", percentile=0.90, cutoff=cut90_hh, n=n_hh,
             pct_values_above=pctval90_hh, pct_flux_from_above=pctflux90_hh,
             total_flux=total_hh, hot_flux=hot90_hh),
  data.frame(aggregation="half-hourly", percentile=0.95, cutoff=cut95_hh, n=n_hh,
             pct_values_above=pctval95_hh, pct_flux_from_above=pctflux95_hh,
             total_flux=total_hh, hot_flux=hot95_hh)
)

## hourly
x_hr <- all_hr_aggr[["EC_FCH4_mean"]]
x_hr <- x_hr[is.finite(x_hr) & x_hr > 0] 
n_hr <- length(x_hr); total_hr <- sum(x_hr)
cut90_hr <- as.numeric(quantile(x_hr, 0.90, na.rm = TRUE, type = 7))
cut95_hr <- as.numeric(quantile(x_hr, 0.95, na.rm = TRUE, type = 7))
hot90_hr <- sum(x_hr[x_hr > cut90_hr]); hot95_hr <- sum(x_hr[x_hr > cut95_hr])
pctval90_hr <- 100 * mean(x_hr > cut90_hr); pctval95_hr <- 100 * mean(x_hr > cut95_hr)
pctflux90_hr <- 100 * hot90_hr / total_hr; pctflux95_hr <- 100 * hot95_hr / total_hr
df_hr <- rbind(
  data.frame(aggregation="hourly", percentile=0.90, cutoff=cut90_hr, n=n_hr,
             pct_values_above=pctval90_hr, pct_flux_from_above=pctflux90_hr,
             total_flux=total_hr, hot_flux=hot90_hr),
  data.frame(aggregation="hourly", percentile=0.95, cutoff=cut95_hr, n=n_hr,
             pct_values_above=pctval95_hr, pct_flux_from_above=pctflux95_hr,
             total_flux=total_hr, hot_flux=hot95_hr)
)

## daily
x_dy <- all_d_aggr[["EC_FCH4_MEAN"]]
x_dy <- x_dy[is.finite(x_dy) & x_dy > 0] 
n_dy <- length(x_dy); total_dy <- sum(x_dy)
cut90_dy <- as.numeric(quantile(x_dy, 0.90, na.rm = TRUE, type = 7))
cut95_dy <- as.numeric(quantile(x_dy, 0.95, na.rm = TRUE, type = 7))
hot90_dy <- sum(x_dy[x_dy > cut90_dy]); hot95_dy <- sum(x_dy[x_dy > cut95_dy])
pctval90_dy <- 100 * mean(x_dy > cut90_dy); pctval95_dy <- 100 * mean(x_dy > cut95_dy)
pctflux90_dy <- 100 * hot90_dy / total_dy; pctflux95_dy <- 100 * hot95_dy / total_dy
df_dy <- rbind(
  data.frame(aggregation="daily", percentile=0.90, cutoff=cut90_dy, n=n_dy,
             pct_values_above=pctval90_dy, pct_flux_from_above=pctflux90_dy,
             total_flux=total_dy, hot_flux=hot90_dy),
  data.frame(aggregation="daily", percentile=0.95, cutoff=cut95_dy, n=n_dy,
             pct_values_above=pctval95_dy, pct_flux_from_above=pctflux95_dy,
             total_flux=total_dy, hot_flux=hot95_dy)
)

## weekly
x_wk <- all_week_aggr[["EC_FCH4_MEAN"]]
x_wk <- x_wk[is.finite(x_wk) & x_wk > 0] 
n_wk <- length(x_wk); total_wk <- sum(x_wk)
cut90_wk <- as.numeric(quantile(x_wk, 0.90, na.rm = TRUE, type = 7))
cut95_wk <- as.numeric(quantile(x_wk, 0.95, na.rm = TRUE, type = 7))
hot90_wk <- sum(x_wk[x_wk > cut90_wk]); hot95_wk <- sum(x_wk[x_wk > cut95_wk])
pctval90_wk <- 100 * mean(x_wk > cut90_wk); pctval95_wk <- 100 * mean(x_wk > cut95_wk)
pctflux90_wk <- 100 * hot90_wk / total_wk; pctflux95_wk <- 100 * hot95_wk / total_wk
df_wk <- rbind(
  data.frame(aggregation="weekly", percentile=0.90, cutoff=cut90_wk, n=n_wk,
             pct_values_above=pctval90_wk, pct_flux_from_above=pctflux90_wk,
             total_flux=total_wk, hot_flux=hot90_wk),
  data.frame(aggregation="weekly", percentile=0.95, cutoff=cut95_wk, n=n_wk,
             pct_values_above=pctval95_wk, pct_flux_from_above=pctflux95_wk,
             total_flux=total_wk, hot_flux=hot95_wk)
)

## monthly
x_mo <- all_month_aggr[["EC_FCH4_MEAN"]]
x_mo <- x_mo[is.finite(x_mo) & x_mo > 0] 
n_mo <- length(x_mo); total_mo <- sum(x_mo)
cut90_mo <- as.numeric(quantile(x_mo, 0.90, na.rm = TRUE, type = 7))
cut95_mo <- as.numeric(quantile(x_mo, 0.95, na.rm = TRUE, type = 7))
hot90_mo <- sum(x_mo[x_mo > cut90_mo]); hot95_mo <- sum(x_mo[x_mo > cut95_mo])
pctval90_mo <- 100 * mean(x_mo > cut90_mo); pctval95_mo <- 100 * mean(x_mo > cut95_mo)
pctflux90_mo <- 100 * hot90_mo / total_mo; pctflux95_mo <- 100 * hot95_mo / total_mo
df_mo <- rbind(
  data.frame(aggregation="monthly", percentile=0.90, cutoff=cut90_mo, n=n_mo,
             pct_values_above=pctval90_mo, pct_flux_from_above=pctflux90_mo,
             total_flux=total_mo, hot_flux=hot90_mo),
  data.frame(aggregation="monthly", percentile=0.95, cutoff=cut95_mo, n=n_mo,
             pct_values_above=pctval95_mo, pct_flux_from_above=pctflux95_mo,
             total_flux=total_mo, hot_flux=hot95_mo)
)

## annual
x_yr <- all_yr_aggr[["EC_FCH4_MEAN"]]
x_yr <- x_yr[is.finite(x_yr) & x_yr > 0] 
n_yr <- length(x_yr); total_yr <- sum(x_yr)
cut90_yr <- as.numeric(quantile(x_yr, 0.90, na.rm = TRUE, type = 7))
cut95_yr <- as.numeric(quantile(x_yr, 0.95, na.rm = TRUE, type = 7))
hot90_yr <- sum(x_yr[x_yr > cut90_yr]); hot95_yr <- sum(x_yr[x_yr > cut95_yr])
pctval90_yr <- 100 * mean(x_yr > cut90_yr); pctval95_yr <- 100 * mean(x_yr > cut95_yr)
pctflux90_yr <- 100 * hot90_yr / total_yr; pctflux95_yr <- 100 * hot95_yr / total_yr
df_yr <- rbind(
  data.frame(aggregation="annual", percentile=0.90, cutoff=cut90_yr, n=n_yr,
             pct_values_above=pctval90_yr, pct_flux_from_above=pctflux90_yr,
             total_flux=total_yr, hot_flux=hot90_yr),
  data.frame(aggregation="annual", percentile=0.95, cutoff=cut95_yr, n=n_yr,
             pct_values_above=pctval95_yr, pct_flux_from_above=pctflux95_yr,
             total_flux=total_yr, hot_flux=hot95_yr)
)

## combine
tail_share_table_means <- rbind(df_hh, df_hr, df_dy, df_wk, df_mo, df_yr)
tail_share_table_means

#### calculate p90s and percentages based on unaggregated data and annual CH4 flux

# organize site levels
site_levels <- sort(unique(all_sites_pred_no_dupl$SITE))

# create YEAR column
all_sites_pred_no_dupl <- all_sites_pred_no_dupl %>% mutate(YEAR = year(TIMESTAMP_START))

# chamber: % of annual emissions from half-hourly > site-year p90 (positives only)
ch_siteyear_flux90 <- all_sites_pred_no_dupl %>%
  group_by(SITE, YEAR) %>%
  # keep only positive (emission) fluxes for percentile + totals
  summarise(
    cutoff90 = as.numeric(quantile(ch_FCH4_nmolCH4m2s1[ch_FCH4_nmolCH4m2s1 > 0], 0.90,
                                   na.rm = TRUE, type = 7)),
    total_flux = sum(ch_FCH4_nmolCH4m2s1[ch_FCH4_nmolCH4m2s1 > 0], na.rm = TRUE),
    hot_flux   = {
      x <- ch_FCH4_nmolCH4m2s1
      pos <- x > 0
      sum(x[pos & x > cutoff90], na.rm = TRUE)
    },
    n_pos = sum(ch_FCH4_nmolCH4m2s1 > 0, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    pct_flux_above_p90 = ifelse(is.finite(total_flux) & total_flux > 0,
                                100 * hot_flux / total_flux, NA_real_),
    method = "Chamber"
  )

# EC: % of annual emissions from half-hourlies > site-year p90 (positives only)
ec_siteyear_flux90 <- all_sites_pred_no_dupl %>%
  group_by(SITE, YEAR) %>%
  summarise(
    cutoff90 = as.numeric(quantile(EC_FCH4_comb[EC_FCH4_comb > 0], 0.90,
                                   na.rm = TRUE, type = 7)),
    total_flux = sum(EC_FCH4_comb[EC_FCH4_comb > 0], na.rm = TRUE),
    hot_flux   = {
      x <- EC_FCH4_comb
      pos <- x > 0
      sum(x[pos & x > cutoff90], na.rm = TRUE)
    },
    n_pos = sum(EC_FCH4_comb > 0, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    pct_flux_above_p90 = ifelse(is.finite(total_flux) & total_flux > 0,
                                100 * hot_flux / total_flux, NA_real_),
    method = "EC"
  )

# combine and plot
site_levels <- sort(unique(all_sites_pred_no_dupl$SITE))

plot_df <- bind_rows(ch_siteyear_flux90, ec_siteyear_flux90) %>%
  mutate(SITE = factor(SITE, levels = site_levels)) %>%
  arrange(SITE, YEAR, method)

ggplot(plot_df, aes(x = SITE, y = pct_flux_above_p90, shape = method, fill = method)) +
  geom_point(position = position_jitterdodge(jitter.width = 0.08, dodge.width = 0.5),
             size = 4, alpha = 0.8, color = "black", stroke = 0.5) +
  scale_shape_manual(values = c("Chamber" = 21, "EC" = 24)) + 
  scale_fill_manual(values = c("Chamber" = "#d8031c",
                               "EC" = "#6497bf")) +
  labs(
    x = "",
    y = "% of annual CH4 emission > site–year p90",
    shape = NULL
  ) +
  theme_bw() +
  scale_y_continuous(limits = c(0, 100), expand = c(0, 0)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 16),
        axis.text.y = element_text(size=16),
        axis.title.y = element_text(size=16),
        legend.text = element_text(size=16),
        legend.position = "top",
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank())

#####################################################################################################

#### ENVIRONMENTAL VARIABLES FOR TABLE 1 ####

# WTL MEANS

# use the raw data
# ensure that there are no duplicates in the data
all_sites_pred_no_dupl <- all_sites_pred_no_dupl[, !duplicated(names(all_sites_pred_no_dupl)) | !seq_along(names(all_sites_pred_no_dupl)) == match(names(all_sites_pred_no_dupl), names(all_sites_pred_no_dupl))]

# calculate site WTL based on the available WTL values for each site
all_sites_pred_no_dupl <- all_sites_pred_no_dupl %>%
  mutate(
    EC_WTL_F = if_else(SITE == "US-STJ", EC_WTL_F, EC_WTL_F * 100)
  ) %>%
  mutate(
    SITE_WTL = case_when(
      SITE %in% c("FI-SI2", "US-LA1", "US-LA2") ~ ch_WTL,
      SITE %in% c("US-HO1", "US-OWC", "US-UAF", "US-STJ") ~ EC_WTL_F, # include US-STJ if you want EC here
      SITE == "CN-HGU" ~ NA_real_,
      TRUE ~ rowMeans(cbind(ch_WTL, EC_WTL_F), na.rm = TRUE)
    )
  )

# calculate the annual mean water table depth (WTL) for each site and year based on the EC and chamber overlap period
annual_mean_wtd <- all_sites_pred_no_dupl %>%
  filter(SITE != "CN-HGU") %>%
  group_by(SITE, YEAR = year(TIMESTAMP_START)) %>%
  summarize(Annual_Mean_WTL = mean(SITE_WTL, na.rm = TRUE),
            annual_CV = raster::cv(SITE_WTL, na.rm=TRUE)) %>%
  ungroup()

# calculate the mean of annual means for each site
mean_annual_wtd_per_site <- annual_mean_wtd %>%
  group_by(SITE) %>%
  summarize(Mean_Annual_WTL = mean(Annual_Mean_WTL, na.rm = TRUE),
            mean_annual_CV = mean(annual_CV, na.rm = TRUE))

mean_annual_wtd_per_site

# get min and max WTL per site
all_sites_pred_no_dupl %>%
  group_by(SITE) %>%
  summarise(
    min_SITE_WTL = min(SITE_WTL, na.rm = TRUE),
    max_SITE_WTL = max(SITE_WTL, na.rm = TRUE)
  )

# same for air temperature (available from FLUXNET-CH4 EC data for each site)
annual_mean_ta <- all_sites_pred_no_dupl %>%
  group_by(SITE, YEAR = year(TIMESTAMP_START)) %>%
  summarize(Annual_Mean_TA = mean(EC_TA_F, na.rm = TRUE)) %>%
  ungroup()

# calculate the mean of annual means for each site
annual_mean_ta %>%
  group_by(SITE) %>%
  summarize(Mean_Annual_TA = mean(Annual_Mean_TA, na.rm = TRUE))

# same for precipitation
annual_precip <- all_sites_pred_no_dupl %>%
  group_by(SITE, YEAR = year(TIMESTAMP_START)) %>%
  summarize(
    Annual_Precip = sum(EC_P_F, na.rm = TRUE),   # sum instead of mean
    .groups = "drop"
  )

# calculate the mean of annual sums for each site
p <- annual_precip %>%
  group_by(SITE) %>%
  summarize(
    Mean_Annual_Precip = mean(Annual_Precip, na.rm = TRUE),
    .groups = "drop"
  )
p <- as.data.frame(p)
p


#########################################################################################################################################################

######## RELATIONSHIPS BETWEEN ECOSYSTEM- AND PLOT-SCALE FCH4 ##########

# SCATTER PLOT + SPEARMAN

# daily
# untransformed
daily_scatter <- ggplot(all_d_aggr, 
                        aes(x = CH_FCH4_MEDIAN, y = EC_FCH4_MEDIAN)) +
  geom_point(size = 4, shape = 21, color = "black", fill = "darkgrey", stroke = 0.6) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black", linewidth = 0.5) +  
  stat_cor(method = "spearman", label.x.npc = "left", label.y.npc = "top", 
           cor.coef.name = "rho", size = 6) +  
  theme_bw() +
  labs(x = expression(paste("Chamber FCH"[4]," (nmol m"^-2, "s"^-1, ")")),
       y = expression(paste("EC FCH"[4]," (nmol m"^-2, "s"^-1, ")")),
       title = "Daily") +
  theme(axis.text = element_text(size = 18),
        axis.title = element_text(size = 18),
        plot.title = element_text(size = 18, hjust = 0.5)) + 
  scale_y_continuous(limits = c(0, 3300)) + 
  scale_x_continuous(limits = c(0, 3300))
daily_scatter

# asinh-transformed version

# data range checks
range(all_d_aggr$CH_FCH4_MEDIAN, na.rm = TRUE)
range(all_d_aggr$EC_FCH4_median, na.rm = TRUE)

# choose limits
lims <- c(-4, 3220) 
# major/minor ticks in original units; ggplot will place them on the asinh scale
major_breaks <- c(-3, -2, -1, 0, 1, 2, 3, 5, 10, 20, 50, 100, 200, 300, 500, 1000, 2000, 3000)

# plot
daily_scatter_asinh <- ggplot(
  all_d_aggr,
  aes(x = CH_FCH4_MEDIAN, y = EC_FCH4_MEDIAN)
) +
  geom_point(size = 4, shape = 21, color = "black", fill = "darkgrey", stroke = 0.6) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black", linewidth = 0.5) +
  stat_cor(method = "spearman", label.x.npc = "left", label.y.npc = "top",
           cor.coef.name = "rho", size = 6) +
  theme_bw() +
  labs(
    x = expression(paste("Chamber FCH"[4], " (nmol m"^-2, " s"^-1, ")")),
    y = expression(paste("EC FCH"[4], " (nmol m"^-2, " s"^-1, ")")),
    title = "Daily"
  ) +
  theme(
    axis.text.x  = element_text(size = 14, angle = 45, hjust = 1),
    axis.text.y  = element_text(size = 14),
    axis.title   = element_text(size = 18),
    plot.title   = element_text(size = 18, hjust = 0.5),
    legend.title = element_text(size = 16),
    legend.text  = element_text(size = 16),
    panel.grid.major = element_line(colour = "grey70", linewidth = 0.2),
    panel.grid.minor = element_blank()   # hide minor grid
  ) +
  scale_x_continuous(
    trans = asinh_trans(),
    limits = lims,
    breaks = major_breaks[major_breaks >= lims[1] & major_breaks <= lims[2]],
    minor_breaks = NULL                  # no minor ticks/grid
  ) +
  scale_y_continuous(
    trans = asinh_trans(),
    limits = lims,
    breaks = major_breaks[major_breaks >= lims[1] & major_breaks <= lims[2]],
    minor_breaks = NULL
  ) +
  coord_equal()

daily_scatter_asinh

# half-hourly

# untransformed version
halfhourly_scatter <- ggplot(all_hh_aggr, 
                             aes(x = CH_FCH4_MEDIAN, 
                                 y = EC_FCH4)) +
  geom_hex(bins=20) +  # Use hex binning
  scale_fill_viridis_c(option = "magma",
                       trans = "log",
                       breaks = c(1, 10, 100, 1000, 10000),
                       labels = c("1", "10", "100", "1000", "10000"),
                       name = "Count") +  # magma color scale + log count
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", 
              color = "black", linewidth = 0.5) +  # 1:1 line
  stat_cor(method = "spearman", label.x.npc = "left", label.y.npc = "top", 
           cor.coef.name = "rho", size = 6) +  # Spearman correlation
  theme_bw() +
  labs(x = expression(paste("Chamber FCH"[4]," (nmol m"^-2, "s"^-1, ")")),
       y = expression(paste("EC FCH"[4]," (nmol m"^-2, "s"^-1, ")")),
       title = "Half-hourly") +
  theme(axis.text = element_text(size = 18),
        axis.title = element_text(size = 18),
        plot.title = element_text(size = 18, hjust = 0.5),
        legend.title=element_text(size=16),
        legend.text = element_text(size=16),
        panel.background = element_rect(fill = "transparent", color = NA),
        plot.background  = element_rect(fill = "transparent", color = NA),
        legend.background = element_rect(fill = "transparent", color = NA),
        legend.box.background = element_rect(fill = "transparent", color = NA)) +
  scale_y_continuous(limits = c(-125, 650)) +  # Force y-axis max limit to 3300
  scale_x_continuous(limits = c(-125, 650))

halfhourly_scatter

# asinh-transformed version

# data range checks
range(all_hh_aggr$CH_FCH4_MEDIAN, na.rm = TRUE)
range(all_hh_aggr$EC_FCH4,         na.rm = TRUE)

# choose limits
lims <- c(-110, 530)

# major/minor ticks in original units; ggplot will place them on the asinh scale
major_breaks <- c(-200, -100, -50, -20, -10, -5, -2, -1, 0, 1, 2,  5, 10, 20, 50, 100, 200, 500)

# plot
hh_scatter_asinh <- ggplot(
  all_hh_aggr,
  aes(x = CH_FCH4_MEDIAN, y = EC_FCH4)
) +
  geom_hex(bins = 60) +
  scale_fill_viridis_c(trans = "log10", name = "Obs.\ncount", option = "magma") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black", linewidth = 0.5) +
  stat_cor(method = "spearman", label.x.npc = "left", label.y.npc = "top",
           cor.coef.name = "rho", size = 6) +
  theme_bw() +
  labs(
    x = expression(paste("Plot FCH"[4], " (nmol m"^-2, " s"^-1, ")")),
    y = expression(paste("Ecosystem FCH"[4], " (nmol m"^-2, " s"^-1, ")")),
    title = "Half-hourly"
  ) +
  theme(
    axis.text.x  = element_text(size = 14, angle = 45, hjust = 1),
    axis.text.y  = element_text(size = 14),
    axis.title   = element_text(size = 18),
    plot.title   = element_text(size = 18, hjust = 0.5),
    legend.title = element_text(size = 16),
    legend.text  = element_text(size = 16),
    panel.grid.major = element_line(colour = "grey70", linewidth = 0.2),
    panel.grid.minor = element_blank()   # ← hide minor grid
  ) +
  scale_x_continuous(
    trans = asinh_trans(),
    limits = lims,
    breaks = major_breaks[major_breaks >= lims[1] & major_breaks <= lims[2]],
    minor_breaks = NULL                  # ← no minor ticks/grid
  ) +
  scale_y_continuous(
    trans = asinh_trans(),
    limits = lims,
    breaks = major_breaks[major_breaks >= lims[1] & major_breaks <= lims[2]],
    minor_breaks = NULL
  ) +
  coord_equal()

hh_scatter_asinh

# hourly

# untransformed version

hourly_scatter <- ggplot(all_hr_aggr, 
                         aes(x = CH_FCH4_MEDIAN, 
                             y = EC_FCH4_median)) +
  geom_hex(bins=20) +  # Use hex binning
  scale_fill_viridis_c(option = "magma",
                       trans = "log",
                       breaks = c(1, 10, 100, 1000, 10000),
                       labels = c("1", "10", "100", "1000", "10000"),
                       name = "Count") +  # magma color scale + log count
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", 
              color = "black", linewidth = 0.5) +  # 1:1 line
  stat_cor(method = "spearman", label.x.npc = "left", label.y.npc = "top", 
           cor.coef.name = "rho", size = 6) +  # Spearman correlation
  theme_bw() +
  labs(x = expression(paste("Chamber FCH"[4]," (nmol m"^-2, "s"^-1, ")")),
       y = expression(paste("EC FCH"[4]," (nmol m"^-2, "s"^-1, ")")),
       title = "Hourly") +
  theme(axis.text = element_text(size = 18),
        axis.title = element_text(size = 18),
        plot.title = element_text(size = 18, hjust = 0.5),
        legend.text = element_text(size=16),
        legend.title= element_text(size=16),
        panel.background = element_rect(fill = "transparent", color = NA),
        plot.background  = element_rect(fill = "transparent", color = NA),
        legend.background = element_rect(fill = "transparent", color = NA),
        legend.box.background = element_rect(fill = "transparent", color = NA)) +
  scale_y_continuous(limits = c(-10, 500)) +  # Force y-axis max limit to 3300
  scale_x_continuous(limits = c(-10, 500))

hourly_scatter


# asinh-transformed version

# data range checks
range(all_hr_aggr$CH_FCH4_MEDIAN, na.rm = TRUE)
range(all_hr_aggr$EC_FCH4_median,         na.rm = TRUE)

# choose limits: 
lims <- c(-110, 460) 

# major/minor ticks in original units; ggplot will place them on the asinh scale
major_breaks <- c(-200, -100, -50, -20, -10, -5, -2, -1, 0, 1, 2,  5, 10, 20, 50, 100, 200, 300, 500)

# plot
hr_scatter_asinh <- ggplot(
  all_hr_aggr,
  aes(x = CH_FCH4_MEDIAN, y = EC_FCH4_median)
) +
  geom_hex(bins = 60) +
  scale_fill_viridis_c(trans = "log10", name = "Obs.\ncount", option = "magma") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black", linewidth = 0.5) +
  stat_cor(method = "spearman", label.x.npc = "left", label.y.npc = "top",
           cor.coef.name = "rho", size = 6) +
  theme_bw() +
  labs(
    x = expression(paste("Plot FCH"[4], " (nmol m"^-2, " s"^-1, ")")),
    y = expression(paste("Ecosystem FCH"[4], " (nmol m"^-2, " s"^-1, ")")),
    title = "Hourly"
  ) +
  theme(
    axis.text.x  = element_text(size = 14, angle = 45, hjust = 1),
    axis.text.y  = element_text(size = 14),
    axis.title   = element_text(size = 18),
    plot.title   = element_text(size = 18, hjust = 0.5),
    legend.title = element_text(size = 16),
    legend.text  = element_text(size = 16),
    panel.grid.major = element_line(colour = "grey70", linewidth = 0.2),
    panel.grid.minor = element_blank()   # ← hide minor grid
  ) +
  scale_x_continuous(
    trans = asinh_trans(),
    limits = lims,
    breaks = major_breaks[major_breaks >= lims[1] & major_breaks <= lims[2]],
    minor_breaks = NULL                  # ← no minor ticks/grid
  ) +
  scale_y_continuous(
    trans = asinh_trans(),
    limits = lims,
    breaks = major_breaks[major_breaks >= lims[1] & major_breaks <= lims[2]],
    minor_breaks = NULL
  ) +
  coord_equal()

hr_scatter_asinh

## weekly

# untransformed version

weekly_scatter <- ggplot(all_week_aggr, 
                         aes(x = CH_FCH4_MEDIAN, y = EC_FCH4_MEDIAN)) +
  geom_point(size = 4, shape = 21, color = "black", fill = "darkgrey", stroke = 0.6) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black", linewidth = 0.5) +  # 1:1 line
  #geom_smooth(method = "lm", se = FALSE) +  # Optional regression line
  stat_cor(method = "spearman", label.x.npc = "left", label.y.npc = "top", 
           cor.coef.name = "rho", size = 6) +  # Use rho (ρ) instead of R
  # stat_poly_eq(aes(label = paste(..eq.label.., ..rr.label.., ..p.value.label.., sep = "~~~")), 
  #              formula = y ~ x, parse = TRUE, size = 6) +  # Linear regression equation & R²
  theme_bw() +
  labs(x = expression(paste("Chamber FCH"[4]," (nmol m"^-2, "s"^-1, ")")),
       y = expression(paste("EC FCH"[4]," (nmol m"^-2, "s"^-1, ")")),
       title = "Weekly") +
  theme(axis.text = element_text(size = 18),
        axis.title = element_text(size = 18),
        plot.title = element_text(size = 18, hjust = 0.5)) +  # Centered title with size 18
  scale_y_continuous(limits = c(-10, 1670)) +  # Force y-axis max limit to 3300
  scale_x_continuous(limits = c(-10, 1670))    # Force x-axis max limit to 3300
weekly_scatter


# data range checks
range(all_week_aggr$CH_FCH4_MEDIAN, na.rm = TRUE)
range(all_week_aggr$EC_FCH4_MEDIAN, na.rm = TRUE)

# choose limits
lims <- c(-3, 1800) 

# major/minor ticks in original units; ggplot will place them on the asinh scale
major_breaks <- c(-2, -1, 0, 1, 2, 3, 5, 10, 20, 50, 100, 200, 300, 500, 1000, 1500)

# plot
week_scatter_asinh <- ggplot(
  all_week_aggr,
  aes(x = CH_FCH4_MEDIAN, y = EC_FCH4_MEDIAN)
) +
  geom_point(size = 4, shape = 21, color = "black", fill = "darkgrey", stroke = 0.6) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black", linewidth = 0.5) +
  stat_cor(method = "spearman", label.x.npc = "left", label.y.npc = "top",
           cor.coef.name = "rho", size = 6) +
  theme_bw() +
  labs(
    x = expression(paste("Plot FCH"[4], " (nmol m"^-2, " s"^-1, ")")),
    y = expression(paste("Ecosystem FCH"[4], " (nmol m"^-2, " s"^-1, ")")),
    title = "Weekly"
  ) +
  theme(
    axis.text.x  = element_text(size = 14, angle = 45, hjust = 1),
    axis.text.y  = element_text(size = 14),
    axis.title   = element_text(size = 18),
    plot.title   = element_text(size = 18, hjust = 0.5),
    legend.title = element_text(size = 16),
    legend.text  = element_text(size = 16),
    panel.grid.major = element_line(colour = "grey70", linewidth = 0.2),
    panel.grid.minor = element_blank()   # ← hide minor grid
  ) +
  scale_x_continuous(
    trans = asinh_trans(),
    limits = lims,
    breaks = major_breaks[major_breaks >= lims[1] & major_breaks <= lims[2]],
    minor_breaks = NULL                  # ← no minor ticks/grid
  ) +
  scale_y_continuous(
    trans = asinh_trans(),
    limits = lims,
    breaks = major_breaks[major_breaks >= lims[1] & major_breaks <= lims[2]],
    minor_breaks = NULL
  ) +
  coord_equal()

week_scatter_asinh

# monthly

# untransformed version

monthly_scatter <- ggplot(all_month_aggr, 
                          aes(x = CH_FCH4_MEDIAN, y = EC_FCH4_MEDIAN)) +
  geom_point(size = 4, shape = 21, color = "black", fill = "darkgrey", stroke = 0.6) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black", linewidth = 0.5) +
  stat_cor(method = "spearman", label.x.npc = "left", label.y.npc = "top", 
           cor.coef.name = "rho", size = 6) +
  theme_bw() +
  labs(x = expression(paste("Chamber FCH"[4]," (nmol m"^-2, "s"^-1, ")")),
       y = expression(paste("EC FCH"[4]," (nmol m"^-2, "s"^-1, ")")),
       title = "Monthly") +
  theme(axis.text = element_text(size = 18),
        axis.title = element_text(size = 18),
        plot.title = element_text(size = 18, hjust = 0.5)) +
  scale_y_continuous(limits = c(-10, 1670)) + 
  scale_x_continuous(limits = c(-10, 1670))
monthly_scatter

# data range checks
range(all_month_aggr$CH_FCH4_MEDIAN, na.rm = TRUE)
range(all_month_aggr$EC_FCH4_MEDIAN,         na.rm = TRUE)


# asinh-transformed version

# choose limits
lims <- c(-2, 1660)

# major/minor ticks in original units; ggplot will place them on the asinh scale
major_breaks <- c(-2, -1, 0, 1, 2, 3, 5, 10, 20, 50, 100, 200, 300, 500, 1000, 1500)

# plot
month_scatter_asinh <- ggplot(
  all_month_aggr,
  aes(x = CH_FCH4_MEDIAN, y = EC_FCH4_MEDIAN)
) +
  geom_point(size = 4, shape = 21, color = "black", fill = "darkgrey", stroke = 0.6) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black", linewidth = 0.5) +
  stat_cor(method = "spearman", label.x.npc = "left", label.y.npc = "top",
           cor.coef.name = "rho", size = 6) +
  theme_bw() +
  labs(
    x = expression(paste("Plot FCH"[4], " (nmol m"^-2, " s"^-1, ")")),
    y = expression(paste("Ecosystem FCH"[4], " (nmol m"^-2, " s"^-1, ")")),
    title = "Monthly"
  ) +
  theme(
    axis.text.x  = element_text(size = 14, angle = 45, hjust = 1),
    axis.text.y  = element_text(size = 14),
    axis.title   = element_text(size = 18),
    plot.title   = element_text(size = 18, hjust = 0.5),
    legend.title = element_text(size = 16),
    legend.text  = element_text(size = 16),
    panel.grid.major = element_line(colour = "grey70", linewidth = 0.2),
    panel.grid.minor = element_blank()   # ← hide minor grid
  ) +
  scale_x_continuous(
    trans = asinh_trans(),
    limits = lims,
    breaks = major_breaks[major_breaks >= lims[1] & major_breaks <= lims[2]],
    minor_breaks = NULL                  # ← no minor ticks/grid
  ) +
  scale_y_continuous(
    trans = asinh_trans(),
    limits = lims,
    breaks = major_breaks[major_breaks >= lims[1] & major_breaks <= lims[2]],
    minor_breaks = NULL
  ) +
  coord_equal()

month_scatter_asinh

## annual

# untransformed version

annual_scatter <- ggplot(all_yr_aggr, 
                         aes(x = CH_FCH4_MEDIAN, y = EC_FCH4_MEDIAN)) +
  geom_point(size = 4, shape = 21, color = "black", fill = "darkgrey", stroke = 0.6) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black", linewidth = 0.5) +  
  stat_cor(method = "spearman", label.x.npc = "left", label.y.npc = "top", 
           cor.coef.name = "rho", size = 6) + 
  theme_bw() +
  labs(x = expression(paste("Chamber FCH"[4]," (nmol m"^-2, "s"^-1, ")")),
       y = expression(paste("EC FCH"[4]," (nmol m"^-2, "s"^-1, ")")),
       title = "Annual") +
  theme(axis.text = element_text(size = 18),
        axis.title = element_text(size = 18),
        plot.title = element_text(size = 18, hjust = 0.5)) + 
  scale_y_continuous(limits = c(-2, 880)) + 
  scale_x_continuous(limits = c(-2, 880)) 
annual_scatter

# asinh-transformed version

# data range checks 
range(all_yr_aggr$CH_FCH4_MEDIAN, na.rm = TRUE)
range(all_yr_aggr$EC_FCH4_MEDIAN,         na.rm = TRUE)

# choose limits

lims <- c(-2, 860) 
# major/minor ticks in original units; ggplot will place them on the asinh scale
major_breaks <- c(-3, -2, -1, 0, 1, 2, 3, 5, 10, 20, 50, 100, 200, 300, 500, 800)

# plot
year_scatter_asinh <- ggplot(
  all_yr_aggr,
  aes(x = CH_FCH4_MEDIAN, y = EC_FCH4_MEDIAN)
) +
  geom_point(size = 4, shape = 21, color = "black", fill = "darkgrey", stroke = 0.6) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black", linewidth = 0.5) +
  stat_cor(method = "spearman", label.x.npc = "left", label.y.npc = "top",
           cor.coef.name = "rho", size = 6) +
  theme_bw() +
  labs(
    x = expression(paste("Plot FCH"[4], " (nmol m"^-2, " s"^-1, ")")),
    y = expression(paste("Ecosystem FCH"[4], " (nmol m"^-2, " s"^-1, ")")),
    title = "Annual"
  ) +
  theme(
    axis.text.x  = element_text(size = 14, angle = 45, hjust = 1),
    axis.text.y  = element_text(size = 14),
    axis.title   = element_text(size = 18),
    plot.title   = element_text(size = 18, hjust = 0.5),
    legend.title = element_text(size = 16),
    legend.text  = element_text(size = 16),
    panel.grid.major = element_line(colour = "grey70", linewidth = 0.2),
    panel.grid.minor = element_blank()   # ← hide minor grid
  ) +
  scale_x_continuous(
    trans = asinh_trans(),
    limits = lims,
    breaks = major_breaks[major_breaks >= lims[1] & major_breaks <= lims[2]],
    minor_breaks = NULL                  # ← no minor ticks/grid
  ) +
  scale_y_continuous(
    trans = asinh_trans(),
    limits = lims,
    breaks = major_breaks[major_breaks >= lims[1] & major_breaks <= lims[2]],
    minor_breaks = NULL
  ) +
  coord_equal()

year_scatter_asinh


# combine into one figure
require(grid)

# untransformed
all_scatter_spearman <- ggarrange(
  halfhourly_scatter + rremove("ylab") + rremove("xlab") + theme(plot.margin = unit(c(0.1, 0.05, 0.1, 0.05), "cm")),
  hourly_scatter + rremove("ylab") + rremove("xlab") + theme(plot.margin = unit(c(0.1, 0.05, 0.1, 0.05), "cm")), 
  daily_scatter + rremove("ylab") + rremove("xlab") + theme(plot.margin = unit(c(0.1, 0.05, 0.1, 0.05), "cm")), 
  weekly_scatter + rremove("ylab") + rremove("xlab") + theme(plot.margin = unit(c(0.1, 0.05, 0.1, 0.05), "cm")), 
  monthly_scatter + rremove("ylab") + rremove("xlab")+ theme(plot.margin = unit(c(0.1, 0.05, 0.1, 0.05), "cm")),
  annual_scatter + rremove("ylab") + rremove("xlab") + theme(plot.margin = unit(c(0.1, 0.05, 0.1, 0.05), "cm")),
  ncol = 3, 
  nrow = 2,
  align = "hv"
)

all_scatter_spearman <- annotate_figure(
  all_scatter_spearman, 
  left = textGrob(expression(paste("EC nmol m"^-2, " s"^-1)), rot = 90, vjust = 1, hjust = 0.3, gp = gpar(fontsize = 20)),
  bottom = textGrob(expression(paste("Chamber nmol m"^-2, " s"^-1)), gp = gpar(fontsize = 20))
)

all_scatter_spearman

# asinh-transformed

### common legend for half-hourly and hourly plots: 

## compute a shared max for the hex counts
b_hh <- ggplot_build(hh_scatter_asinh)
b_hr <- ggplot_build(hr_scatter_asinh)

max_count <- max(
  max(b_hh$data[[1]]$count, na.rm = TRUE),
  max(b_hr$data[[1]]$count, na.rm = TRUE)
)

shared_fill <- scale_fill_viridis_c(
  trans  = "log10",
  name   = "Obs.\ncount",
  option = "magma",
  limits = c(1, max_count),
  breaks = breaks_log(n = 4),
  oob    = squish
)

## apply the exact same fill scale to hh + hr
hh_scatter_asinh2 <- hh_scatter_asinh + shared_fill
hr_scatter_asinh2 <- hr_scatter_asinh + shared_fill

## remove legends from the other plots
daily2 <- daily_scatter_asinh + theme(legend.position = "none")
week2  <- week_scatter_asinh  + theme(legend.position = "none")
month2 <- month_scatter_asinh + theme(legend.position = "none")
year2  <- year_scatter_asinh  + theme(legend.position = "none")

## ggarrange with a common legend
all_scatter_spearman_asinh <- ggarrange(
  hh_scatter_asinh2 + rremove("ylab") + rremove("xlab") + theme(plot.margin = unit(c(0.1, 0.05, 0.1, 0.05), "cm")),
  hr_scatter_asinh2 + rremove("ylab") + rremove("xlab") + theme(plot.margin = unit(c(0.1, 0.05, 0.1, 0.05), "cm")),
  daily2 + rremove("ylab") + rremove("xlab") + theme(plot.margin = unit(c(0.1, 0.05, 0.1, 0.05), "cm")),
  week2  + rremove("ylab") + rremove("xlab") + theme(plot.margin = unit(c(0.1, 0.05, 0.1, 0.05), "cm")),
  month2 + rremove("ylab") + rremove("xlab") + theme(plot.margin = unit(c(0.1, 0.05, 0.1, 0.05), "cm")),
  year2  + rremove("ylab") + rremove("xlab") + theme(plot.margin = unit(c(0.1, 0.05, 0.1, 0.05), "cm")),
  ncol = 3, nrow = 2, align = "hv",
  common.legend = TRUE, legend = "right"
)

all_scatter_spearman_asinh <- annotate_figure(
  all_scatter_spearman_asinh,
  left   = textGrob(expression(paste("Ecosystem nmol m"^-2, " s"^-1)), rot = 90, vjust = 1, hjust = 0.3, gp = gpar(fontsize = 20)),
  bottom = textGrob(expression(paste("Plot nmol m"^-2, " s"^-1)), gp = gpar(fontsize = 20))
)

all_scatter_spearman_asinh


### calculate NRMSE for the scatter plots

# daily

sqrt(
  mean(
    (all_d_aggr$EC_FCH4_MEDIAN -
       all_d_aggr$CH_FCH4_MEDIAN)^2,
    na.rm = TRUE
  )
) / sd(
  c(all_d_aggr$EC_FCH4_MEDIAN,
    all_d_aggr$CH_FCH4_MEDIAN),
  na.rm = TRUE
)

# weekly

sqrt(
  mean(
    (all_week_aggr$EC_FCH4_MEDIAN -
       all_week_aggr$CH_FCH4_MEDIAN)^2,
    na.rm = TRUE
  )
) / sd(
  c(all_week_aggr$EC_FCH4_MEDIAN,
    all_week_aggr$CH_FCH4_MEDIAN),
  na.rm = TRUE
)


# monthly

sqrt(
  mean(
    (all_month_aggr$EC_FCH4_MEDIAN -
       all_month_aggr$CH_FCH4_MEDIAN)^2,
    na.rm = TRUE
  )
) / sd(
  c(all_month_aggr$EC_FCH4_MEDIAN,
    all_month_aggr$CH_FCH4_MEDIAN),
  na.rm = TRUE
)

# annual

sqrt(
  mean(
    (all_yr_aggr$EC_FCH4_MEDIAN -
       all_yr_aggr$CH_FCH4_MEDIAN)^2,
    na.rm = TRUE
  )
) / sd(
  c(all_yr_aggr$EC_FCH4_MEDIAN,
    all_yr_aggr$CH_FCH4_MEDIAN),
  na.rm = TRUE
)


# hourly

sqrt(
  mean(
    (all_hr_aggr$EC_FCH4_median -
       all_hr_aggr$CH_FCH4_MEDIAN)^2,
    na.rm = TRUE
  )
) / sd(
  c(all_hr_aggr$EC_FCH4_median,
    all_hr_aggr$CH_FCH4_MEDIAN),
  na.rm = TRUE
)

# half-hourly

sqrt(
  mean(
    (all_hh_aggr$EC_FCH4 -
       all_hh_aggr$CH_FCH4_MEDIAN)^2,
    na.rm = TRUE
  )
) / sd(
  c(all_hh_aggr$EC_FCH4,
    all_hh_aggr$CH_FCH4_MEDIAN),
  na.rm = TRUE
)


##########################################################################################################################################

##### RELATIONSHIP BETWEEN ABSOLUTE DELTA_FCH4 AND FCH4 MAGNITUDE (ROW-WISE MEAN BETWEEN EC AND CHAMBER FCH4), FIG. 4

### calculate the row-wise mean between EC and chamber FCH4

all_d_aggr <- mutate(all_d_aggr, FCH4_mean = rowMeans(dplyr::select(all_d_aggr, 
                                                                                          c(EC_FCH4_MEDIAN, CH_FCH4_MEDIAN)), na.rm = TRUE)) 

all_hr_aggr <- mutate(all_hr_aggr, FCH4_mean = rowMeans(dplyr::select(all_hr_aggr, 
                                                                      c(EC_FCH4_median, CH_FCH4_MEDIAN)), na.rm = TRUE)) 

all_week_aggr <- mutate(all_week_aggr, FCH4_mean = rowMeans(dplyr::select(all_week_aggr, 
                                                                          c(EC_FCH4_MEDIAN, CH_FCH4_MEDIAN)), na.rm = TRUE)) 

all_month_aggr <- mutate(all_month_aggr, FCH4_mean = rowMeans(dplyr::select(all_month_aggr, 
                                                                            c(EC_FCH4_MEDIAN, CH_FCH4_MEDIAN)), na.rm = TRUE)) 

all_yr_aggr <- mutate(all_yr_aggr, FCH4_mean = rowMeans(dplyr::select(all_yr_aggr, 
                                                                      c(EC_FCH4_MEDIAN, CH_FCH4_MEDIAN)), na.rm = TRUE)) 

all_hh_aggr <- mutate(all_hh_aggr, FCH4_mean = rowMeans(dplyr::select(all_hh_aggr, 
                                                                                              c(EC_FCH4, CH_FCH4_MEDIAN)), na.rm = TRUE)) 

## calculate % of observations that belong to CN-HGU where DELTA_FCH4 < 0

# hourly
all_hr_aggr %>%
  # Filter for FCH4_mean < 0
  filter(FCH4_mean < 0) %>%
  # Summarize the counts
  summarize(
    # Count how many observations belong to SITE = "CN-HGU"
    cn_hgu_count = sum(SITE == "CN-HGU"),
    # Total count of observations where FCH4_mean < 0 (all sites)
    total_count = n(),
    # Calculate the percentage of CN-HGU observations
    percentage_cn_hgu = (cn_hgu_count / total_count) * 100
  )

# half-hourly
all_hh_aggr %>%
  # Filter for FCH4_mean < 0
  filter(FCH4_mean < 0) %>%
  # Summarize the counts
  summarize(
    # Count how many observations belong to SITE = "CN-HGU"
    cn_hgu_count = sum(SITE == "CN-HGU"),
    # Total count of observations where FCH4_mean < 0 (all sites)
    total_count = n(),
    # Calculate the percentage of CN-HGU observations
    percentage_cn_hgu = (cn_hgu_count / total_count) * 100
  )


# calculate percentage of how many delta_FCH4 are negative and positive when FCH4_mean < 0 or FCH4_mean > 0 across sites

# create a new column for the sign of delta_FCH4 for all aggregations

all_d_aggr$DELTA_FCH4_orig_sign <- sign(all_d_aggr$DELTA_FCH4)
all_hh_aggr$DELTA_FCH4_orig_sign <- sign(all_hh_aggr$DELTA_FCH4)
all_hr_aggr$DELTA_FCH4_orig_sign <- sign(all_hr_aggr$DELTA_FCH4)
all_week_aggr$DELTA_FCH4_orig_sign <- sign(all_week_aggr$DELTA_FCH4)
all_month_aggr$DELTA_FCH4_orig_sign <- sign(all_month_aggr$DELTA_FCH4)
all_yr_aggr$DELTA_FCH4_orig_sign <- sign(all_yr_aggr$DELTA_FCH4)

# Replace -1 with "negative", 1 with "positive", and 0 with "zero"
all_d_aggr$DELTA_FCH4_orig_sign <- ifelse(all_d_aggr$DELTA_FCH4_orig_sign == -1, "negative",
                                                    ifelse(all_d_aggr$DELTA_FCH4_orig_sign == 1, "positive", "zero"))

all_hh_aggr$DELTA_FCH4_orig_sign <- ifelse(all_hh_aggr$DELTA_FCH4_orig_sign == -1, "negative",
                                                      ifelse(all_hh_aggr$DELTA_FCH4_orig_sign == 1, "positive", "zero"))

all_hr_aggr$DELTA_FCH4_orig_sign <- ifelse(all_hr_aggr$DELTA_FCH4_orig_sign == -1, "negative",
                                          ifelse(all_hr_aggr$DELTA_FCH4_orig_sign == 1, "positive", "zero"))

all_week_aggr$DELTA_FCH4_orig_sign <- ifelse(all_week_aggr$DELTA_FCH4_orig_sign == -1, "negative",
                                            ifelse(all_week_aggr$DELTA_FCH4_orig_sign == 1, "positive", "zero"))

all_month_aggr$DELTA_FCH4_orig_sign <- ifelse(all_month_aggr$DELTA_FCH4_orig_sign == -1, "negative",
                                             ifelse(all_month_aggr$DELTA_FCH4_orig_sign == 1, "positive", "zero"))

all_yr_aggr$DELTA_FCH4_orig_sign <- ifelse(all_yr_aggr$DELTA_FCH4_orig_sign == -1, "negative",
                                          ifelse(all_yr_aggr$DELTA_FCH4_orig_sign == 1, "positive", "zero"))

# hourly and half-hourly

negative_ecchmean_data_hr <- subset(all_hr_aggr, FCH4_mean < 0)
negative_ecchmean_data_hh <- subset(all_hh_aggr, FCH4_mean < 0)

# Calculate the number of negative delta_FCH4 observations
sign_counts_hr <- table(negative_ecchmean_data_hr$DELTA_FCH4_orig_sign)
sign_counts_hh <- table(negative_ecchmean_data_hh$DELTA_FCH4_orig_sign)

# Calculate the total number of negative FCH4_mean observations
total_negative_hr <- nrow(negative_ecchmean_data_hr)
total_negative_hh <- nrow(negative_ecchmean_data_hh)

# Calculate percentages
(sign_counts_hr / total_negative_hr) * 100
(sign_counts_hh / total_negative_hh) * 100

## calculate percentage of how many delta_FCH4 are negative when ECCHmean > 0

pos_ecchmean_data_hr <- subset(all_hr_aggr, FCH4_mean > 0)
pos_ecchmean_data_hh <- subset(all_hh_aggr, FCH4_mean > 0)

# Calculate the number of negative delta_FCH4 observations
sign_counts_hr <- table(pos_ecchmean_data_hr$DELTA_FCH4_orig_sign)
sign_counts_hh <- table(pos_ecchmean_data_hh$DELTA_FCH4_orig_sign)

# Calculate the total number of negative FCH4_mean observations
total_pos_hr <- nrow(pos_ecchmean_data_hr)
total_pos_hh <- nrow(pos_ecchmean_data_hh)

# Calculate percentages
(sign_counts_hr / total_pos_hr) * 100
(sign_counts_hh / total_pos_hh) * 100


# per site

# Calculate the percentage of negative delta_FCH4 for each site where FCH4_mean < 0

all_hr_aggr %>%
  # Filter for FCH4_mean < 0
  filter(FCH4_mean < 0) %>%
  # Group by SITE
  group_by(SITE) %>%
  # Summarize to calculate the percentage of negative observations within FCH4_mean < 0
  summarize(
    # Count how many of the filtered observations have DELTA_FCH4_orig_sign == "negative"
    negative_count = sum(DELTA_FCH4_orig_sign == "negative"),
    # Total count of observations where FCH4_mean < 0
    total_count = n(),
    # Calculate the percentage
    percentage = (negative_count / total_count) * 100
  ) %>%
  # Handle division by zero by replacing NA with 0
  mutate(percentage = replace_na(percentage, 0))

# Calculate the percentage of positive delta_FCH4 for each site where FCH4_mean < 0

all_hr_aggr %>%
  # Filter for FCH4_mean < 0
  filter(FCH4_mean < 0) %>%
  # Group by SITE
  group_by(SITE) %>%
  # Summarize to calculate the percentage of negative observations within FCH4_mean < 0
  summarize(
    # Count how many of the filtered observations have DELTA_FCH4_orig_sign == "negative"
    positive_count = sum(DELTA_FCH4_orig_sign == "positive"),
    # Total count of observations where FCH4_mean < 0
    total_count = n(),
    # Calculate the percentage
    percentage = (positive_count / total_count) * 100
  ) %>%
  # Handle division by zero by replacing NA with 0
  mutate(percentage = replace_na(percentage, 0))


# half-hourly

all_hh_aggr %>%
  # Filter for FCH4_mean < 0
  filter(FCH4_mean < 0) %>%
  # Group by SITE
  group_by(SITE) %>%
  # Summarize to calculate the percentage of negative observations within FCH4_mean < 0
  summarize(
    # Count how many of the filtered observations have DELTA_FCH4_orig_sign == "negative"
    positive_count = sum(DELTA_FCH4_orig_sign == "negative"),
    # Total count of observations where FCH4_mean < 0
    total_count = n(),
    # Calculate the percentage
    percentage = (positive_count / total_count) * 100
  ) %>%
  # Handle division by zero by replacing NA with 0
  mutate(percentage = replace_na(percentage, 0))

all_hh_aggr %>%
  # Filter for FCH4_mean < 0
  filter(FCH4_mean < 0) %>%
  # Group by SITE
  group_by(SITE) %>%
  # Summarize to calculate the percentage of negative observations within FCH4_mean < 0
  summarize(
    # Count how many of the filtered observations have DELTA_FCH4_orig_sign == "negative"
    positive_count = sum(DELTA_FCH4_orig_sign == "positive"),
    # Total count of observations where FCH4_mean < 0
    total_count = n(),
    # Calculate the percentage
    percentage = (positive_count / total_count) * 100
  ) %>%
  # Handle division by zero by replacing NA with 0
  mutate(percentage = replace_na(percentage, 0))


# calculate percentage of how many delta_FCH4 are negative and positive when FCH4_mean < 0 in daily aggregation

negative_data <- subset(all_d_aggr, DELTA_FCH4_orig_sign == "negative")

# Calculate the number of sites where negative delta_FCH4 observations
site_counts <- table(negative_data$SITE)

# Calculate the total number of negative FCH4_mean observations
total_negative <- nrow(negative_data)

# Calculate percentages for each site
percentages <- (site_counts / total_negative) * 100
percentages

# calculate percentage of how many delta_FCH4 are negative and positive when FCH4_mean < 0

negative_ecchmean_data <- subset(all_d_aggr, FCH4_mean < 0)

# Calculate the number of negative and positive delta_FCH4 observations
sign_counts <- table(negative_ecchmean_data$DELTA_FCH4_orig_sign)

# Calculate the total number of negative FCH4_mean observations
total_negative <- nrow(negative_ecchmean_data)

# Calculate percentages
(sign_counts / total_negative) * 100

# calculate percentage of how many delta_FCH4 are negative and positive when FCH4_mean > 0

pos_ecchmean_data <- subset(all_d_aggr, FCH4_mean > 0)

# Calculate the number of negative and positive delta_FCH4 observations
sign_counts <- table(pos_ecchmean_data$DELTA_FCH4_orig_sign)

# Calculate the total number of positive FCH4_mean observations
total_pos <- nrow(pos_ecchmean_data)

# Calculate percentages
(sign_counts / total_pos) * 100

# calculate percentage of how many delta_FCH4 are positive and negative when ECCHmean > 0

pos_ecchmean_data <- subset(all_d_aggr, FCH4_mean > 0)

# Calculate the number of positive and negative delta_FCH4 observations per site
sign_counts <- table(pos_ecchmean_data$DELTA_FCH4_orig_sign)

# Calculate the total number of positive FCH4_mean observations
total_pos <- nrow(pos_ecchmean_data)

# Calculate percentages
(sign_counts / total_pos) * 100


# weekly

# calculate percentage of how many delta_FCH4 are negative and positive when ECCHmean < 0

negative_ecchmean_data <- subset(all_week_aggr, FCH4_mean < 0)

# Calculate the number of negative and positive delta_FCH4 observations
sign_counts <- table(negative_ecchmean_data$DELTA_FCH4_orig_sign)

# Calculate the total number of negative FCH4_mean observations
total_negative <- nrow(negative_ecchmean_data)

# Calculate percentages
(sign_counts / total_negative) * 100

# calculate percentage of how many delta_FCH4 are negative and positive when FCH4_mean > 0

pos_ecchmean_data <- subset(all_week_aggr, FCH4_mean > 0)

# Calculate the number of negative and positive delta_FCH4 observations
sign_counts <- table(pos_ecchmean_data$DELTA_FCH4_orig_sign)

# Calculate the total number of positive FCH4_mean observations
total_pos <- nrow(pos_ecchmean_data)

# Calculate percentages
(sign_counts / total_pos) * 100


# monthly

# calculate percentage of how many delta_FCH4 are negative and positive when FCH4_mean < 0

negative_ecchmean_data <- subset(all_month_aggr, FCH4_mean < 0)

# Calculate the number of negative and positive delta_FCH4 observations
sign_counts <- table(negative_ecchmean_data$DELTA_FCH4_orig_sign)

# Calculate the total number of negative FCH4_mean observations
total_negative <- nrow(negative_ecchmean_data)

# Calculate percentages
(sign_counts / total_negative) * 100

# calculate percentage of how many delta_FCH4 are negative when FCH4_mean > 0

pos_ecchmean_data <- subset(all_month_aggr, FCH4_mean > 0)

# Calculate the number of negative and positive delta_FCH4 observations
sign_counts <- table(pos_ecchmean_data$DELTA_FCH4_orig_sign)

# Calculate the total number of positive FCH4_mean observations
total_pos <- nrow(pos_ecchmean_data)

# Calculate percentages
(sign_counts / total_pos) * 100


# annual

# calculate percentage of how many delta_FCH4 are negative and positive when FCH4_mean > 0

pos_ecchmean_data <- subset(all_yr_aggr, FCH4_mean > 0)

# Calculate the number of negative and positive delta_FCH4 observations
sign_counts <- table(pos_ecchmean_data$DELTA_FCH4_orig_sign)

# Calculate the total number of positive FCH4_mean observations
total_pos <- nrow(pos_ecchmean_data)

# Calculate percentages
(sign_counts / total_pos) * 100


# calculate % of obs where delta_FCH4 sign is negative and site is cn-hgu

# Filter for negative DELTA_FCH4_orig_sign and FCH4_mean < 0 for CN-HGU
negative_CNHGU_data <- subset(all_hr_aggr, 
                               DELTA_FCH4_orig_sign == "negative" & 
                                 FCH4_mean < 0 & 
                                 SITE == "CN-HGU")

# Calculate the number of negative observations for CN-HGU
count_negative_CNHGU <- nrow(negative_CNHGU_data)

# Calculate the total number of negative observations for CN-HGU
total_negative_CNHGU <- nrow(subset(all_hr_aggr, 
                                     DELTA_FCH4_orig_sign == "negative" & 
                                       SITE == "CN-HGU"))

# Calculate the percentage for CN-HGU
percentage_negative_CNHGU <- (count_negative_CNHGU / total_negative_CNHGU) * 100
percentage_negative_CNHGU


# hourly plot
## density cloud

blandaltman_hourly_abs <- ggplot(all_hr_aggr, 
                                 aes(x = FCH4_mean, 
                                     y = abs(DELTA_FCH4))) +
  geom_hex(bins=20) +
  scale_fill_viridis_c(option = "magma",
                       trans = "log",
                       breaks = c(1, 10, 100, 1000, 10000),
                       labels = c("1", "10", "100", "1000", "10000"),
                       name = "Count") +
  geom_hline(yintercept = 0, color = "blue", linetype = "dashed", linewidth = 1) +
  theme_bw() +
  ylab(expression(paste("EC - CH difference (nmol m"^-2, "s"^-1, ")"))) +
  xlab(expression(paste("EC - CH mean (nmol m"^-2, "s"^-1, ")"))) +
  ggtitle("Hourly") +
  theme(axis.text = element_text(size = 18),
        axis.title = element_text(size = 18),
        plot.title = element_text(size = 18, hjust = 0.5),
        legend.title=element_text(size=16),
        legend.text = element_text(size=16),
        panel.background = element_rect(fill = "transparent", color = NA),
        plot.background  = element_rect(fill = "transparent", color = NA),
        legend.background = element_rect(fill = "transparent", color = NA),
        legend.box.background = element_rect(fill = "transparent", color = NA)) +
  scale_y_continuous(limits = c(-1, 770))
blandaltman_hourly_abs


# plots without outliers (for Fig. 4) and with outliers (Fig. B17)

# remove outliers

# daily

# Sort the dataframe by FCH4_mean in descending order and exclude the top two values
all_d_aggr$abs_DELTA_FCH4 <- abs(all_d_aggr$DELTA_FCH4)

# Sort the dataframe by the absolute values of DELTA_FCH4 in descending order and exclude the top three values
filtered_daily <- all_d_aggr[order(all_d_aggr$abs_DELTA_FCH4, decreasing = TRUE), ]
filtered_daily<- filtered_daily[-c(1, 2, 3), ]

# Remove the temporary column
all_d_aggr$abs_DELTA_FCH4 <- NULL

### ---> 3 outliers from abs(DELTA_FCH4) removed

# plot without outliers (Fig. 4)

blandaltman_daily_abs_filtered <- ggplot(filtered_daily, aes(x = FCH4_mean, y = abs(DELTA_FCH4))) +
  geom_point(size=4, shape = 21, stroke = 0.6, color = "black", fill = "darkgrey") +
  geom_hline(yintercept = 0, color = "blue", linetype = "dashed", linewidth = 1) +
  ylab(expression(paste("Absolute "*Delta*"FCH"[4]," (nmol m"^-2, "s"^-1, ")"))) +
  xlab(expression(paste("FCH"[4*"_"*mean]," (nmol m"^-2, "s"^-1, ")"))) +
  ggtitle("Daily") +
  theme_bw() +
  theme(
    axis.title.x = element_text(size = 20),     
    axis.title.y = element_text(size = 20),     
    axis.text = element_text(size = 20),
    plot.title = element_text(size = 20, hjust = 0.5),
    legend.title = element_text(size = 20),  
    legend.text = element_text(size = 20)      
  ) +
  scale_y_continuous(limits = c(0, 770))
blandaltman_daily_abs_filtered

# with outliers (Fig. B17)

blandaltman_daily_abs <- ggplot(all_d_aggr, aes(x = FCH4_mean, y = abs(DELTA_FCH4))) +
  geom_point(size=4, shape = 21, stroke = 0.6, color = "black", fill = "darkgrey") +
  geom_hline(yintercept = 0, color = "blue", linetype = "dashed", linewidth = 1) +
  ylab(expression(paste("Absolute "*Delta*"FCH"[4]," (nmol m"^-2, "s"^-1, ")"))) +
  xlab(expression(paste("FCH"[4*"_"*mean]," (nmol m"^-2, "s"^-1, ")"))) +
  ggtitle("Daily") +
  theme_bw() +
  theme(
    axis.title.x = element_text(size = 20),     
    axis.title.y = element_text(size = 20),     
    axis.text = element_text(size = 20),
    plot.title = element_text(size = 20, hjust = 0.5),
    legend.title = element_text(size = 20),  
    legend.text = element_text(size = 20)      
  ) +
  scale_y_continuous(limits = c(0, 1000))
blandaltman_daily_abs


# weekly

# Sort the dataframe by FCH4_mean in descending order and exclude the top three values
filtered_weekly <- all_week_aggr[order(all_week_aggr$FCH4_mean, decreasing = TRUE), ]
filtered_weekly<- filtered_weekly[-c(1, 2, 3, 4, 5, 6, 7), ]

# ---> 7 outlier points removed from FCH4_mean

# and then the highest 2 from abs DELTA_FCH4
# Sort the dataframe by the absolute values of DELTA_FCH4 in descending order and exclude the top three values
filtered_weekly$DELTA_FCH4_abs <- abs(filtered_weekly$DELTA_FCH4)
filtered_weekly <- filtered_weekly[order(filtered_weekly$DELTA_FCH4_abs, decreasing = TRUE), ]
filtered_weekly<- filtered_weekly[-c(1, 2, 3), ]

# ----> 3 outliers removed from DELTA_FCH4_abs

# plot without outliers (Fig. 4)

blandaltman_weekly_abs_filtered <- ggplot(filtered_weekly, aes(x = FCH4_mean, y = abs(DELTA_FCH4))) +
  geom_point(size=4, shape = 21, stroke = 0.6, color = "black",  fill = "darkgrey") +
  geom_hline(yintercept = 0, color = "blue", linetype = "dashed", linewidth = 1) +
  ylab(expression(paste("Absolute "*Delta*"FCH"[4]," (nmol m"^-2, "s"^-1, ")"))) +
  xlab(expression(paste("FCH"[4*"_"*mean]," (nmol m"^-2, "s"^-1, ")"))) +
  ggtitle("Weekly") +
  theme_bw() +
  theme(
    axis.title.x = element_text(size = 20),     
    axis.title.y = element_text(size = 20),     
    axis.text = element_text(size = 20),
    plot.title = element_text(size = 20, hjust = 0.5),
    legend.title = element_text(size = 20),  
    legend.text = element_text(size = 20)      
  ) +
  scale_y_continuous(limits = c(0, 250))

blandaltman_weekly_abs_filtered

# plot with outliers (Fig. B17)

blandaltman_weekly_abs <- ggplot(all_week_aggr, aes(x = FCH4_mean, y = abs(DELTA_FCH4))) +
  geom_point(size=4, shape = 21, stroke = 0.6, color = "black", fill = "darkgrey") +
  geom_hline(yintercept = 0, color = "blue", linetype = "dashed", linewidth = 1) +
  ylab(expression(paste("Absolute "*Delta*"FCH"[4]," (nmol m"^-2, "s"^-1, ")"))) +
  xlab(expression(paste("FCH"[4*"_"*mean]," (nmol m"^-2, "s"^-1, ")"))) +
  ggtitle("Weekly") +
  theme_bw() +
  theme(
    axis.title.x = element_text(size = 20),     
    axis.title.y = element_text(size = 20),     
    axis.text = element_text(size = 20),
    plot.title = element_text(size = 20, hjust = 0.5),
    legend.title = element_text(size = 20),  
    legend.text = element_text(size = 20)      
  ) +
  scale_y_continuous(limits = c(0, 250))

blandaltman_weekly_abs


# monthly

# Sort the dataframe by FCH4_mean in descending order and exclude the top three values
filtered_monthly <- all_month_aggr[order(all_month_aggr$FCH4_mean, decreasing = TRUE), ]
filtered_monthly<- filtered_monthly[-c(1, 2, 3, 4, 5, 6), ]

# ---> 6 outlier points removed from FCH4_mean

# and then the highest 2 from abs DELTA_FCH4
# Sort the dataframe by the absolute values of DELTA_FCH4 in descending order and exclude the top three values
filtered_monthly$DELTA_FCH4_abs <- abs(filtered_monthly$DELTA_FCH4)
filtered_monthly <- filtered_monthly[order(filtered_monthly$DELTA_FCH4_abs, decreasing = TRUE), ]
filtered_monthly<- filtered_monthly[-c(1, 2), ]

# ----> 2 outliers removed from DELTA_FCH4_abs

# plot without outliers (Fig. 4)
blandaltman_monthly_abs_filtered <- ggplot(filtered_monthly, aes(x = FCH4_mean, y = abs(DELTA_FCH4))) +
  geom_point(size=4, shape = 21, stroke = 0.6, color = "black", fill = "darkgrey") +
  geom_hline(yintercept = 0, color = "blue", linetype = "dashed", linewidth = 1) +
  ylab(expression(paste("Absolute "*Delta*"FCH"[4]," (nmol m"^-2, "s"^-1, ")"))) +
  xlab(expression(paste("FCH"[4*"_"*mean]," (nmol m"^-2, "s"^-1, ")"))) +
  ggtitle("Monthly") +
  theme_bw() +
  theme(
    axis.title.x = element_text(size = 20),     
    axis.title.y = element_text(size = 20),     
    axis.text = element_text(size = 20),
    plot.title = element_text(size = 20, hjust = 0.5),
    legend.title = element_text(size = 20),  
    legend.text = element_text(size = 20)      
  ) +
  scale_y_continuous(limits = c(0, 250))

blandaltman_monthly_abs_filtered

# plot without outliers (Fig. B17)

blandaltman_monthly_abs <- ggplot(all_month_aggr, aes(x = FCH4_mean, y = abs(DELTA_FCH4))) +
  geom_point(size=4, shape = 21, stroke = 0.6, color = "black", fill = "darkgrey") +
  geom_hline(yintercept = 0, color = "blue", linetype = "dashed", linewidth = 1) +
  ylab(expression(paste("Absolute "*Delta*"FCH"[4]," (nmol m"^-2, "s"^-1, ")"))) +
  xlab(expression(paste("FCH"[4*"_"*mean]," (nmol m"^-2, "s"^-1, ")"))) +
  ggtitle("Monthly") +
  theme_bw() +
  theme(
    axis.title.x = element_text(size = 20),     
    axis.title.y = element_text(size = 20),     
    axis.text = element_text(size = 20),
    plot.title = element_text(size = 20, hjust = 0.5),
    legend.title = element_text(size = 20),  
    legend.text = element_text(size = 20)      
  ) +
  scale_y_continuous(limits = c(0, 250))

blandaltman_monthly_abs

# annual

# remove outlier
filtered_annual$DELTA_FCH4_abs <- abs(filtered_annual$DELTA_FCH4)
filtered_annual <- filtered_annual[order(filtered_annual$DELTA_FCH4_abs, decreasing = TRUE), ]
filtered_annual<- filtered_annual[-c(1), ]

# plot without outlier (Fig. 4)

blandaltman_annual_abs_filtered <- ggplot(filtered_annual, aes(x = FCH4_mean, y = abs(DELTA_FCH4))) +
  geom_point(size=4, shape = 21, stroke = 0.6, color = "black", fill = "darkgrey") +
  geom_hline(yintercept = 0, color = "blue", linetype = "dashed", linewidth = 1) +
  ylab(expression(paste("Absolute "*Delta*"FCH"[4]," (nmol m"^-2, "s"^-1, ")"))) +
  xlab(expression(paste("FCH"[4*"_"*mean]," (nmol m"^-2, "s"^-1, ")"))) +
  ggtitle("Annual") +
  theme_bw() +
  theme(
    axis.title.x = element_text(size = 20),     
    axis.title.y = element_text(size = 20),     
    axis.text = element_text(size = 20),
    plot.title = element_text(size = 20, hjust = 0.5),
    legend.title = element_text(size = 20),  
    legend.text = element_text(size = 20)      
  ) +
  scale_y_continuous(limits = c(0, 250))
blandaltman_annual_abs_filtered 

# with outlier (Fig. B17)

blandaltman_annual_abs <- ggplot(all_yr_aggr, aes(x = FCH4_mean, y = abs(DELTA_FCH4))) +
  geom_point(size=4, shape = 21, stroke = 0.6, color = "black", fill = "darkgrey") +
  geom_hline(yintercept = 0, color = "blue", linetype = "dashed", linewidth = 1) +
  ylab(expression(paste("Absolute "*Delta*"FCH"[4]," (nmol m"^-2, "s"^-1, ")"))) +
  xlab(expression(paste("FCH"[4*"_"*mean]," (nmol m"^-2, "s"^-1, ")"))) +
  ggtitle("Annual") +
  theme_bw() +
  theme(
    axis.title.x = element_text(size = 20),     
    axis.title.y = element_text(size = 20), 
    plot.title = element_text(size = 20, hjust = 0.5),
    axis.text = element_text(size = 20),
    legend.title = element_text(size = 20),  
    legend.text = element_text(size = 20)      
  ) +
  scale_y_continuous(limits = c(0, 250)) +
  scale_x_continuous(limits = c(0, 160))
blandaltman_annual_abs 


# combine plots into one figure

require(grid)

# Fig. 4
blandaltman_all_temporal_abs_no_outliers <- ggarrange(
  blandaltman_halfhourly_abs + rremove("ylab") + rremove("xlab") + theme(plot.margin = unit(c(0.1, 0.05, 0.1, 0.05), "cm")),
  blandaltman_hourly_abs + rremove("ylab") + rremove("xlab")+ theme(plot.margin = unit(c(0.1, 0.05, 0.1, 0.05), "cm")), 
  blandaltman_daily_abs_filtered + rremove("ylab") + rremove("xlab")+ theme(plot.margin = unit(c(0.1, 0.05, 0.1, 0.05), "cm")), 
  blandaltman_weekly_abs_filtered + rremove("ylab") + rremove("xlab") + theme(plot.margin = unit(c(0.1, 0.05, 0.1, 0.05), "cm")), 
  blandaltman_monthly_abs_filtered + rremove("ylab") + rremove("xlab")+ theme(plot.margin = unit(c(0.1, 0.05, 0.1, 0.05), "cm")),
  blandaltman_annual_abs + rremove("ylab") + rremove("xlab")+ theme(plot.margin = unit(c(0.1, 0.05, 0.1, 0.05), "cm")),
  ncol = 3, 
  nrow = 2, 
  align = "hv"
)

blandaltman_all_temporal_abs_no_outliers <- annotate_figure(
  blandaltman_all_temporal_abs_no_outliers, 
  left = textGrob(expression(paste("Absolute "*Delta*"FCH"[4]," (nmol m"^-2, "s"^-1, ")")), rot = 90, vjust = 1, hjust = 0.3, gp = gpar(fontsize = 20)),
  bottom = textGrob(expression(paste("FCH"[4*"_"*mean]," (nmol m"^-2, "s"^-1, ")")), gp = gpar(fontsize = 20))
)

blandaltman_all_temporal_abs_no_outliers

# correlations
cor.test(all_hh_aggr$FCH4_mean, abs(all_hh_aggr$DELTA_FCH4), method = "spearman")
cor.test(all_hr_aggr$FCH4_mean, abs(all_hr_aggr$DELTA_FCH4), method = "spearman")
cor.test(all_d_aggr$FCH4_mean, abs(all_d_aggr$DELTA_FCH4), method = "spearman")
cor.test(all_week_aggr$FCH4_mean, abs(all_week_aggr$DELTA_FCH4), method = "spearman")
cor.test(all_month_aggr$FCH4_mean, abs(all_month_aggr$DELTA_FCH4), method = "spearman")
cor.test(all_yr_aggr$FCH4_mean, abs(all_yr_aggr$DELTA_FCH4), method = "spearman")

#### combine the legend for half-hourly and hourly plots

## compute a shared max for the hex counts
b_hh <- ggplot_build(blandaltman_halfhourly_abs)
b_hr <- ggplot_build(blandaltman_hourly_abs)

max_count <- max(
  max(b_hh$data[[1]]$count, na.rm = TRUE),
  max(b_hr$data[[1]]$count, na.rm = TRUE)
)

shared_fill <- scale_fill_viridis_c(
  trans  = "log10",
  name   = "Count",
  option = "magma",
  limits = c(1, max_count),
  breaks = breaks_log(n = 4),
  oob    = squish
)

## apply the exact same fill scale to hh + hr
blandaltman_halfhourly_abs2 <- blandaltman_halfhourly_abs + shared_fill
blandaltman_hourly_abs2     <- blandaltman_hourly_abs     + shared_fill

## remove legends from the other plots
daily2  <- blandaltman_daily_abs_filtered   + theme(legend.position = "none")
week2   <- blandaltman_weekly_abs_filtered  + theme(legend.position = "none")
month2  <- blandaltman_monthly_abs_filtered + theme(legend.position = "none")
annual2 <- blandaltman_annual_abs           + theme(legend.position = "none")

## ggarrange with a common legend
blandaltman_all_temporal_abs_no_outliers <- ggarrange(
  blandaltman_halfhourly_abs2 + rremove("ylab") + rremove("xlab") + theme(plot.margin = unit(c(0.1, 0.05, 0.1, 0.05), "cm")),
  blandaltman_hourly_abs2     + rremove("ylab") + rremove("xlab") + theme(plot.margin = unit(c(0.1, 0.05, 0.1, 0.05), "cm")),
  daily2  + rremove("ylab") + rremove("xlab") + theme(plot.margin = unit(c(0.1, 0.05, 0.1, 0.05), "cm")),
  week2   + rremove("ylab") + rremove("xlab") + theme(plot.margin = unit(c(0.1, 0.05, 0.1, 0.05), "cm")),
  month2  + rremove("ylab") + rremove("xlab") + theme(plot.margin = unit(c(0.1, 0.05, 0.1, 0.05), "cm")),
  annual2 + rremove("ylab") + rremove("xlab") + theme(plot.margin = unit(c(0.1, 0.05, 0.1, 0.05), "cm")),
  ncol = 3, nrow = 2, align = "hv",
  common.legend = TRUE, legend = "right"
)

## add shared axis labels
blandaltman_all_temporal_abs_no_outliers <- annotate_figure(
  blandaltman_all_temporal_abs_no_outliers,
  left   = textGrob(expression(paste("Absolute "*Delta*"FCH"[4], " (nmol m"^-2, " s"^-1, ")")),
                    rot = 90, vjust = 1, hjust = 0.3, gp = gpar(fontsize = 20)),
  bottom = textGrob(expression(paste("FCH"[4*"_"*mean], " (nmol m"^-2, " s"^-1, ")")),
                    gp = gpar(fontsize = 20))
)

blandaltman_all_temporal_abs_no_outliers


#### calculate NRMSE for each aggregation

# half-hourly

sqrt(
  mean(all_hh_aggr$DELTA_FCH4^2, na.rm = TRUE)
) / sd(
  c(all_hh_aggr$EC_FCH4,
    all_hh_aggr$CH_FCH4_MEDIAN),
  na.rm = TRUE
)

# hourly
sqrt(
  mean(all_hr_aggr$DELTA_FCH4^2, na.rm = TRUE)
) / sd(
  c(all_hr_aggr$EC_FCH4_median,
    all_hh_aggr$CH_FCH4_MEDIAN),
  na.rm = TRUE
)

# daily
sqrt(
  mean(all_d_aggr$DELTA_FCH4^2, na.rm = TRUE)
) / sd(
  c(all_d_aggr$EC_FCH4_MEDIAN,
    all_d_aggr$CH_FCH4_MEDIAN),
  na.rm = TRUE
)

# weekly
sqrt(
  mean(all_week_aggr$DELTA_FCH4^2, na.rm = TRUE)
) / sd(
  c(all_week_aggr$EC_FCH4_MEDIAN,
    all_week_aggr$CH_FCH4_MEDIAN),
  na.rm = TRUE
)

# monthly

sqrt(
  mean(all_month_aggr$DELTA_FCH4^2, na.rm = TRUE)
) / sd(
  c(all_month_aggr$EC_FCH4_MEDIAN,
    all_month_aggr$CH_FCH4_MEDIAN),
  na.rm = TRUE
)

# annual

sqrt(
  mean(all_yr_aggr$DELTA_FCH4^2, na.rm = TRUE)
) / sd(
  c(all_yr_aggr$EC_FCH4_MEDIAN,
    all_yr_aggr$CH_FCH4_MEDIAN),
  na.rm = TRUE
)


# with outliers (Fig. B17)

blandaltman_all_temporal_with_outliers <- ggarrange(
  blandaltman_daily_abs + rremove("ylab") + rremove("xlab")  + theme(plot.margin = unit(c(0.1, 0.05, 0.3, 0.05), "cm")), 
  blandaltman_weekly_abs + rremove("ylab") + rremove("xlab")+ theme(plot.margin = unit(c(0.1, 0.05, 0.3, 0.05), "cm")), 
  blandaltman_monthly_abs + rremove("ylab") + rremove("xlab") + theme(plot.margin = unit(c(0.1, 0.05, 0.3, 0.05), "cm")),
  blandaltman_annual_abs + rremove("ylab") + rremove("xlab") + theme(plot.margin = unit(c(0.1, 0.05, 0.3, 0.05), "cm")),
  ncol = 2, 
  nrow = 2, 
  common.legend = TRUE, 
  legend = "right",
  align = "hv"
)

blandaltman_all_temporal_with_outliers <- annotate_figure(
  blandaltman_all_temporal_with_outliers, 
  left = textGrob(expression(paste("Absolute "*Delta*"FCH"[4]," (nmol m"^-2, "s"^-1, ")")), rot = 90, vjust = 1, hjust = 0.3, gp = gpar(fontsize = 20)),
  bottom = textGrob(expression(paste("FCH"[4*"_"*mean]," (nmol m"^-2, "s"^-1, ")")), gp = gpar(fontsize = 20))
)

blandaltman_all_temporal_with_outliers

###########################################################################################################################################

### RELATIONSHIP BETWEEN FCH4 VARIATION BETWEEN CHAMBERS OR EC FCH4 MEASUREMENT TIMESTAMPS AND DELTA_FCH4 
### (Fig. 5)

# calculate variation (variance, standard deviation, IQR and CV) per site per temporal aggregation unit (e.g. day, week, month, year)
# based on the unaggregated dataset

# chamber
# per date and site
all_sites_pred_no_dupl <- all_sites_pred_no_dupl %>%
  group_by(SITE, DATE) %>%
  mutate(
    ch_var_perdatesite = var(ch_FCH4_nmolCH4m2s1, na.rm = TRUE),
    ch_sd_perdatesite = sd(ch_FCH4_nmolCH4m2s1, na.rm = TRUE), 
    ch_iqr_perdatesite = IQR(ch_FCH4_nmolCH4m2s1, na.rm= TRUE),
    ch_CV_perdatesite = raster::cv(ch_FCH4_nmolCH4m2s1, na.rm=TRUE) 
  )

# per week, year and site
all_sites_pred_no_dupl <- all_sites_pred_no_dupl %>%
  group_by(SITE, YEAR, isoweek(DATE)) %>%
  mutate(
    ch_var_perweeksite = var(ch_FCH4_nmolCH4m2s1, na.rm = TRUE),
    ch_sd_perweeksite = sd(ch_FCH4_nmolCH4m2s1, na.rm = TRUE), 
    ch_iqr_perweeksite = IQR(ch_FCH4_nmolCH4m2s1, na.rm = TRUE),
    ch_CV_perweeksite = raster::cv(ch_FCH4_nmolCH4m2s1, na.rm=TRUE)  
  )

# per month, year and site
all_sites_pred_no_dupl <- all_sites_pred_no_dupl %>%
  group_by(SITE, YEAR, MONTH) %>%
  mutate(
    ch_var_permonthsite = var(ch_FCH4_nmolCH4m2s1, na.rm = TRUE),
    ch_sd_permonthsite = sd(ch_FCH4_nmolCH4m2s1, na.rm = TRUE), 
    ch_iqr_permonthsite = IQR(ch_FCH4_nmolCH4m2s1, na.rm = TRUE), 
    ch_CV_permonthsite = raster::cv(ch_FCH4_nmolCH4m2s1, na.rm=TRUE)  
  )

# per year and site
all_sites_pred_no_dupl <- all_sites_pred_no_dupl %>%
  group_by(SITE, YEAR) %>%
  mutate(
    ch_var_peryearsite = var(ch_FCH4_nmolCH4m2s1, na.rm = TRUE),
    ch_sd_peryearsite = sd(ch_FCH4_nmolCH4m2s1, na.rm = TRUE), 
    ch_iqr_peryearsite = IQR(ch_FCH4_nmolCH4m2s1, na.rm = TRUE),
    ch_CV_peryearsite = raster::cv(ch_FCH4_nmolCH4m2s1, na.rm=TRUE)  
  )

all_sites_pred_no_dupl <- as.data.frame(all_sites_pred_no_dupl)


# EC

# per date and site

all_sites_pred_no_dupl <- all_sites_pred_no_dupl %>%
  group_by(SITE, DATE) %>%
  mutate(
    EC_var_perdatesite = var(EC_FCH4_F_comb, na.rm = TRUE), 
    EC_CV_perdatesite = raster::cv(EC_FCH4_F_comb, na.rm=TRUE),
    EC_IQR_perdatesite = IQR(EC_FCH4_F_comb, na.rm = TRUE)
  )

# per week, year and site

all_sites_pred_no_dupl <- all_sites_pred_no_dupl %>%
  group_by(SITE, YEAR, isoweek(DATE)) %>%
  mutate(
    EC_var_perweeksite = var(EC_FCH4_F_comb, na.rm = TRUE), 
    EC_IQR_perweeksite = IQR(EC_FCH4_F_comb, na.rm = TRUE),
    EC_CV_perweeksite = raster::cv(EC_FCH4_F_comb, na.rm=TRUE) 
  )

# per month, year and site

all_sites_pred_no_dupl <- all_sites_pred_no_dupl %>%
  group_by(SITE, YEAR, MONTH) %>%
  mutate(
    EC_var_permonthsite = var(EC_FCH4_F_comb, na.rm = TRUE),
    EC_IQR_permonthsite = IQR(EC_FCH4_F_comb, na.rm = TRUE),
    EC_CV_permonthsite = raster::cv(EC_FCH4_F_comb, na.rm=TRUE)
  )

# per year and site

all_sites_pred_no_dupl <- all_sites_pred_no_dupl %>%
  group_by(SITE, YEAR) %>%
  mutate(
    EC_var_peryearsite = var(EC_FCH4_F_comb, na.rm = TRUE),
    EC_IQR_peryearsite = IQR(EC_FCH4_F_comb, na.rm = TRUE),
    EC_CV_peryearsite = raster::cv(EC_FCH4_F_comb, na.rm=TRUE)
  )

all_sites_pred_no_dupl <- as.data.frame(all_sites_pred_no_dupl)

## subset the new columns with time info

daily_var <- subset(all_sites_pred_no_dupl, select = c(EC_var_perdatesite, EC_CV_perdatesite, EC_IQR_perdatesite, 
                                                       ch_var_perdatesite, ch_iqr_perdatesite, ch_CV_perdatesite, 
                                                       DATE, SITE))

week_var <- subset(all_sites_pred_no_dupl, select = c(EC_var_perweeksite, EC_CV_perweeksite, EC_IQR_perweeksite,
                                                      ch_var_perweeksite, ch_CV_perweeksite, ch_iqr_perweeksite,
                                                      YEAR, Week_of_year, SITE))

month_var <- subset(all_sites_pred_no_dupl, select = c(EC_var_permonthsite, EC_CV_permonthsite, EC_IQR_permonthsite,
                                                       ch_var_permonthsite, ch_CV_permonthsite, ch_iqr_permonthsite,
                                                       YEAR, MONTH, SITE))

year_var <- subset(all_sites_pred_no_dupl, select = c(EC_var_peryearsite, EC_CV_peryearsite, EC_IQR_peryearsite,
                                                      ch_var_peryearsite, ch_CV_peryearsite, ch_iqr_peryearsite,
                                                      YEAR, SITE))

daily_var <- unique(daily_var)

week_var <- unique(week_var)

month_var <- unique(month_var)

year_var <- unique(year_var)

# add these to the aggregations

# daily

colnames(daily_var)[7] <- "TIMESTAMP"

all_d_aggr <- merge(all_d_aggr,daily_var,by=c("SITE","TIMESTAMP"), all = TRUE)
all_d_aggr <- merge(all_d_aggr,daily_var,by=c("SITE","TIMESTAMP"), all = TRUE)

# weekly

all_week_aggr$Week_of_year <- as.character(all_week_aggr$Week_of_year)
week_var$Week_of_year <- as.character(week_var$Week_of_year)

all_week_aggr <- merge(all_week_aggr,week_var,by=c("SITE","Week_of_year", "YEAR"), all = TRUE)

# month

all_month_aggr <- merge(all_month_aggr,month_var,by=c("SITE","YEAR", "MONTH"), all = TRUE)

# year

all_yr_aggr <- merge(all_yr_aggr,year_var,by=c("SITE","YEAR"), all = TRUE)


# plot (use IQR)

# daily

# untransformed (Fig. B18)
# chamber IQR

d_ch_IQR_untransformed <- ggplot(all_d_aggr, 
                                 aes(x = ch_iqr_perdatesite, 
                                     y = abs(DELTA_FCH4), 
                                     color = SITE)) +
  geom_point(size = 3) +  # all sites in color
  # extra layer: only US-UAF, styled manually
  geom_point(data = filter(all_d_aggr, SITE == "US-UAF"),
             aes(x = ch_iqr_perdatesite, y = abs(DELTA_FCH4)),
             shape = 21, fill = "white", color = "black", size = 3, stroke = 1) +
  labs(x = expression("Chamber FCH" [4] * " IQR"), 
       y = expression(paste("Absolute "*Delta*"FCH"[4], 
                            " nmol m"^-2, " s"^-1))) +
  stat_cor(method = "spearman", show.legend = FALSE, size = 5) +
  scale_color_okabeito() +
  theme_bw() +
  theme(axis.title = element_text(size = 20),
        axis.text = element_text(size = 20))

d_ch_IQR_untransformed

# EC IQR

d_ec_IQR_untransformed <- ggplot(all_d_aggr, 
                                 aes(x = EC_IQR_perdatesite, 
                                     y = abs(DELTA_FCH4), 
                                     color = SITE)) +
  geom_point(size = 3) +  # all sites in color
  # extra layer: only US-UAF, styled manually
  geom_point(data = filter(all_d_aggr, SITE == "US-UAF"),
             aes(x = EC_IQR_perdatesite, y = abs(DELTA_FCH4)),
             shape = 21, fill = "white", color = "black", size = 3, stroke = 1) +
  labs(x = expression("EC FCH" [4] * " IQR"), 
       y = expression(paste("Absolute "*Delta*"FCH"[4], 
                            " nmol m"^-2, " s"^-1))) +
  stat_cor(method = "spearman", show.legend = FALSE, size = 5) +
  scale_color_okabeito() +
  theme_bw() +
  theme(axis.title = element_text(size = 20),
        axis.text = element_text(size = 20))


d_ec_IQR_untransformed

# spearman correlations
cor.test(all_d_aggr$ch_iqr_perdatesite, abs(all_d_aggr$DELTA_FCH4), method = "spearman")
cor.test(all_d_aggr$EC_IQR_perdatesite, abs(all_d_aggr$DELTA_FCH4), method = "spearman")

# correlations per site
# chamber
all_d_aggr %>%
  group_by(SITE) %>%
  summarize(correlation = cor.test(ch_iqr_perdatesite, abs(DELTA_FCH4), method = "spearman")$estimate,
            p_value = cor.test(ch_iqr_perdatesite, abs(DELTA_FCH4), method = "spearman")$p.value,
            count = n())

# EC
all_d_aggr %>%
  group_by(SITE) %>%
  summarize(correlation = cor.test(EC_IQR_perdatesite, abs(DELTA_FCH4), method = "spearman")$estimate,
            p_value = cor.test(EC_IQR_perdatesite, abs(DELTA_FCH4), method = "spearman")$p.value,
            count = n())

# log-transformed chamber IQR and absolute delta_FCH4
d_IQR_ch <- ggplot(all_d_aggr, aes(x = log(ch_iqr_perdatesite + 0.01), y = log(abs(DELTA_FCH4)), color = SITE, shape = SITE)) +
  geom_point(size=2) +
  labs(x = expression("Chamber FCH" [4] * " IQR (nmol m"^-2, " s"^-1, ", log +0.01)"), y = expression(paste("Absolute "*Delta*"FCH"[4], " nmol m"^-2, " s"^-1, ", log"))) +
  geom_smooth(aes(x = log(ch_iqr_perdatesite + 0.01), y = log(abs(DELTA_FCH4))), , alpha=0.2, method = "lm", color = "blue", se = T) +  # General regression line
  # Add regression equations for each group (SITE)
  stat_poly_eq(aes(label = paste(..eq.label.., ..rr.label.., ..p.value.label.., sep = "~~~"), color = SITE),
               formula = y ~ x,
               parse = TRUE,
               label.x.npc = "left",
               label.y.npc = "top",
               size = 5) +
  scale_color_okabeito() +
  scale_shape_manual(values = c(0, 1, 2, 3, 4, 5, 6, 7, 8)) +
  theme_bw() +
  theme(axis.title = element_text(size=20),
        axis.text = element_text(size = 20),
        legend.text = element_text(size = 15))

d_IQR_ch



# Okabe–Ito (9 colors) for separating sites
ok_cols <- c(
  "#E69F00", "#56B4E9", "#009E73", "#F0E442",
  "#0072B2", "#D55E00", "#CC79A7", "#999999", "#000000"
)

# Factor-safe site levels
sites <- if (is.factor(all_d_aggr$SITE)) {
  levels(all_d_aggr$SITE)
} else {
  sort(unique(all_d_aggr$SITE))
}

stopifnot("US-UAF" %in% sites)

# Split out US-UAF
other_sites <- setdiff(sites, "US-UAF")

if (length(other_sites) > length(ok_cols)) {
  stop("You have ", length(other_sites), " non-US-UAF sites but only 9 Okabe–Ito colors.")
}

# Colors: Okabe–Ito for others, white for US-UAF
my_cols <- setNames(rep(NA_character_, length(sites)), sites)
my_cols[other_sites] <- ok_cols[seq_along(other_sites)]
my_cols["US-UAF"] <- "white"   # keeps US-UAF white in the legend and base layer

# Shapes: 9 distinct for others; 21 for US-UAF (base layer will be invisible; special layer shows outline)
shape_pool <- c(15,16,17,18,19,20, 0,2,5)
my_shapes <- setNames(rep(NA_integer_, length(sites)), sites)
my_shapes[other_sites] <- shape_pool[seq_along(other_sites)]
my_shapes["US-UAF"] <- 21

# plot log-transformed chamber IQR vs log-transformed absolute delta_FCH4
ggplot(all_d_aggr, 
       aes(x = log(ch_iqr_perdatesite + 0.01), 
           y = log(abs(DELTA_FCH4)), 
           color = SITE, shape = SITE)) +
  geom_point(size = 2) +
  # US-UAF points
  geom_point(data = filter(all_d_aggr, SITE == "US-UAF"),
             aes(x = log(ch_iqr_perdatesite + 0.01), 
                 y = log(abs(DELTA_FCH4))),
             shape = 21, fill = "white", color = "black", size = 1, stroke = 0.8) +
  labs(x = expression("Chamber FCH" [4] * " IQR (nmol m"^-2, " s"^-1, ", log +0.01)"), 
       y = expression(paste("Absolute "*Delta*"FCH"[4], 
                            " nmol m"^-2, " s"^-1, ", log"))) +
  geom_smooth(method = "lm", se = TRUE, alpha = 0.2) +
  geom_smooth(data = filter(all_d_aggr, SITE == "US-UAF"),
              aes(x = log(ch_iqr_perdatesite + 0.01), 
                  y = log(abs(DELTA_FCH4))),
              method = "lm", se = TRUE, alpha = 0.2, color = "black") +
  stat_poly_eq(aes(label = paste(..eq.label.., ..rr.label.., ..p.value.label.., sep = "~~~"), color = SITE),
               formula = y ~ x,
               parse = TRUE,
               label.x.npc = "left",
               label.y.npc = "top",
               size = 5) +
  scale_color_manual(values = my_cols) +
  scale_shape_manual(values = my_shapes) +
  theme_bw() +
  theme(axis.title = element_text(size=20),
        axis.text = element_text(size = 20),
        legend.text = element_text(size = 15))
d_IQR_ch


# per-site regression stats
df <- all_d_aggr %>%
  transmute(
    SITE,
    x = log(ch_iqr_perdatesite + 0.01),
    y = log(abs(DELTA_FCH4))
  ) %>%
  filter(is.finite(x), is.finite(y))

site_stats <- df %>%
  group_by(SITE) %>%
  group_modify(~{
    d <- .x
    n <- nrow(d)
    if (n < 3 || var(d$x) == 0 || var(d$y) == 0) {
      return(tibble(
        n = n, r2 = NA_real_, adj_r2 = NA_real_,
        slope = NA_real_, slope_p = NA_real_, model_p = NA_real_,
        note = "Insufficient data or zero variance"
      ))
    }
    fit <- lm(y ~ x, data = d)
    s   <- summary(fit)
    slope_p <- coef(s)["x","Pr(>|t|)"]
    model_p <- if (!is.null(s$fstatistic)) {
      pf(s$fstatistic[1], s$fstatistic[2], s$fstatistic[3], lower.tail = FALSE)
    } else NA_real_
    tibble(
      n = n,
      r2 = s$r.squared,            
      adj_r2 = s$adj.r.squared, 
      slope = coef(fit)[["x"]],
      slope_p = slope_p, 
      model_p = model_p,
      note = NA_character_
    )
  }) %>%
  ungroup()

# See the table
site_stats %>%
  mutate(across(c(r2, adj_r2, slope), ~round(.x, 3)),
         across(c(slope_p, model_p), ~signif(.x, 3))) %>%
  arrange(desc(r2)) %>%
  print(n = Inf)


# based on these stats, color only significant predictor sites

d_IQR_ch <- ggplot(all_d_aggr, 
                   aes(x = log(ch_iqr_perdatesite + 0.01), 
                       y = log(abs(DELTA_FCH4)))) +
  geom_hline(yintercept = 0, color = "blue", linetype = "dashed", linewidth = 1) +
  # Set color and shape mappings in geom_point to avoid grouping by SITE in geom_smooth
  geom_point(aes(color = site_color, shape = site_color), size = 2, stroke=1) +  # Add points with shapes and colors
  labs(x = expression("Chamber FCH" [4] * " IQR (nmol m"^-2, " s"^-1, ", log +0.01)"), 
       y = expression("Absolute "*Delta*"FCH"[4] * "(nmol m"^{-2}*" s"^{-1}*"), log")) +  # Label axes
  
  # Specify colors for chosen sites and light gray for "Other"
  scale_color_manual(values = c("SE-DEG" = "#009E73", 
                                "CN-HGU" = "#E69F00", 
                                "US-HO1" = "#F5C710", 
                                "US-UAF" = "black", 
                                "US-LA1" = "#0072B2",
                                "Other" = "lightgray")) +
  
  # Specify shapes for chosen sites and circle (16) for "Other"
  scale_shape_manual(values = c("SE-DEG" = 0, 
                                "CN-HGU" = 1, 
                                "US-HO1" = 2, 
                                "US-UAF" = 3, 
                                "US-LA1" = 4,
                                "Other" = 16)) +
  
  # Display regression equations for each group without drawing the lines
  stat_poly_eq(aes(label = paste(..eq.label.., ..rr.label.., ..p.value.label.., sep = "~~~"), 
                   color = site_color),
               formula = y ~ x,
               parse = TRUE,
               label.x.npc = "left",
               label.y.npc = "top",
               size = 5) +
  theme_bw() +
  theme(axis.title = element_text(size = 20),
        axis.text = element_text(size = 20),
        legend.text = element_text(size = 15),
        legend.title = element_blank())

d_IQR_ch

all_d_aggr <- all_d_aggr %>%
  arrange(ifelse(SITE == "SE-DEG", 1, 0))

## log-transformed EC IQR (first show all sites)

d_IQR_EC <- ggplot(all_d_aggr, aes(x = log(EC_IQR_perdatesite + 0.01), y = log(abs(DELTA_FCH4)), color = SITE)) +
  geom_point(size=2) +
  labs(x = expression("EC FCH" [4] * " IQR (nmol m"^-2, " s"^-1, ", log +0.01)"), y = expression(paste("Log absolute "*Delta*"FCH"[4], " nmol m"^-2, " s"^-1))) + 
  stat_poly_eq(aes(label = paste(..eq.label.., ..rr.label.., ..p.value.label.., sep = "~~~"), color = SITE),
               formula = y ~ x,
               parse = TRUE,
               label.x.npc = "left",
               label.y.npc = "top",
               size = 5) +
  scale_color_okabeito()+
  theme_bw() +
  theme(axis.title = element_text(size=20),
        axis.text = element_text(size = 20),
        legend.text = element_text(size = 15))
d_IQR_EC


# per-site regression stats
df <- all_d_aggr %>%
  transmute(
    SITE,
    x = log(EC_IQR_perdatesite + 0.01),
    y = log(abs(DELTA_FCH4))
  ) %>%
  filter(is.finite(x), is.finite(y))

site_stats <- df %>%
  group_by(SITE) %>%
  group_modify(~{
    d <- .x
    n <- nrow(d)
    if (n < 3 || var(d$x) == 0 || var(d$y) == 0) {
      return(tibble(
        n = n, r2 = NA_real_, adj_r2 = NA_real_,
        slope = NA_real_, slope_p = NA_real_, model_p = NA_real_,
        note = "Insufficient data or zero variance"
      ))
    }
    fit <- lm(y ~ x, data = d)
    s   <- summary(fit)
    slope_p <- coef(s)["x","Pr(>|t|)"]
    model_p <- if (!is.null(s$fstatistic)) {
      pf(s$fstatistic[1], s$fstatistic[2], s$fstatistic[3], lower.tail = FALSE)
    } else NA_real_
    tibble(
      n = n,
      r2 = s$r.squared,           
      adj_r2 = s$adj.r.squared,   
      slope = coef(fit)[["x"]],
      slope_p = slope_p,      
      model_p = model_p,
      note = NA_character_
    )
  }) %>%
  ungroup()

# See the table
site_stats %>%
  mutate(across(c(r2, adj_r2, slope), ~round(.x, 3)),
         across(c(slope_p, model_p), ~signif(.x, 3))) %>%
  arrange(desc(r2)) %>%
  print(n = Inf)


# color only significant predictors
all_d_aggr$site_color_ec <- ifelse(all_d_aggr$SITE %in% c("SE-DEG"), all_d_aggr$SITE, "Other")

# Reorder data so that US-HO1 is second last and SE-DEG is last
all_d_aggr <- all_d_aggr %>%
  arrange(ifelse(site_color_ec == "SE-DEG", 2, ifelse(site_color_ec == "US-HO1", 1, 0)))

# plot with only significant sites
d_IQR_ec <- ggplot(all_d_aggr, 
                   aes(x = log(ec_iqr_perdatesite_plus001), 
                       y = log(abs(DELTA_FCH4)))) +
  # Set color and shape mappings in geom_point and add stroke for larger borders
  geom_point(aes(color = site_color_ec, shape = site_color_ec), size = 2, stroke = 1) + 
  labs(x = expression("EC FCH" [4] * " IQR (nmol m"^-2, " s"^-1, ", log +0.01)"), 
       y = expression(paste("Absolute "*Delta*"FCH"[4], " nmol m"^-2, " s"^-1, ", log"))) + 
  
  # General regression line only for the whole dataset
  geom_smooth(method = "lm", color = "blue", se = TRUE) +
  
  # Display regression equations for each group without drawing the lines
  stat_poly_eq(aes(label = paste(..eq.label.., ..rr.label.., ..p.value.label.., sep = "~~~"), 
                   color = site_color_ec),
               formula = y ~ x,
               parse = TRUE,
               label.x.npc = "left",
               label.y.npc = "top",
               size = 5) +
  
  # Display the general regression equation
  stat_poly_eq(aes(label = paste(..eq.label.., ..rr.label.., ..p.value.label.., sep = "~~~")),
               formula = y ~ x,
               parse = TRUE,
               label.x.npc = "right",
               label.y.npc = "bottom",
               size = 5,
               color = "blue") +
  
  # Define colors for specific sites and light gray for "Other"
  scale_color_manual(values = c("SE-DEG" = "#009E73", 
                                "US-HO1" = "#F5C710", 
                                "Other" = "lightgray")) +
  
  # Define shapes for specific sites and circle (16) for "Other"
  scale_shape_manual(values = c("SE-DEG" = 0, 
                                "US-HO1" = 2, 
                                "Other" = 16)) +
  
  theme_bw() +
  theme(axis.title = element_text(size = 20),
        axis.text = element_text(size = 20),
        legend.text = element_text(size = 15))

d_IQR_ec


## run linear mixed models

# ADD NONA DF VERSION CODE HERE

# daily
# chamber IQR
m_d_iqr <- lme(
  log(abs(DELTA_FCH4)) ~ log(ch_iqr_perdatesite + 0.01), 
  random = ~1 | SITE, 
  data = d_noNA
)

summary(m_d_iqr)

r.squaredGLMM(m_d_iqr)

# chamber and EC IQR
m_d_iqr_2 <- lme(
  log(abs(DELTA_FCH4)) ~ log(ch_iqr_perdatesite + 0.01) + log(EC_IQR_perdatesite + 0.01), 
  random = ~1 | SITE, 
  data = d_noNA
)

summary(m_d_iqr)
summary(m_d_iqr_2)

r.squaredGLMM(m_d_iqr)
r.squaredGLMM(m_d_iqr_2)



# weekly

# chamber IQR
m_week_iqr <- lme(
  log(abs(DELTA_FCH4)) ~ log(ch_iqr_perweeksite + 0.01), 
  random = ~1 | SITE, 
  data = week_noNA
)

summary(m_week_iqr)
r.squaredGLMM(m_week_iqr)

# EC IQR
m_week_iqr <- lme(
  log(abs(DELTA_FCH4)) ~ log(EC_IQR_perweeksite + 0.01), 
  random = ~1 | SITE, 
  data = week_noNA
)

summary(m_week_iqr)
r.squaredGLMM(m_week_iqr)

# chamber and EC IQR
m_week_iqr_2 <- lme(
  log(abs(DELTA_FCH4)) ~ log(ch_iqr_perweeksite + 0.01) + log(EC_IQR_perweeksite + 0.01), 
  random = ~1 | SITE, 
  data = week_noNA
)
summary(m_week_iqr_2)
r.squaredGLMM(m_week_iqr_2)


# monthly

# chamber IQR
m_month_iqr <- lme(
  log(abs(DELTA_FCH4)) ~ log(ch_iqr_permonthsite + 0.01), 
  random = ~1 | SITE, 
  data = month_noNA
)

summary(m_month_iqr)
r.squaredGLMM(m_month_iqr)

# EC IQR
m_month_iqr <- lme(
  log(abs(DELTA_FCH4)) ~ log(EC_IQR_permonthsite + 0.01), 
  random = ~1 | SITE, 
  data = month_noNA
)

summary(m_month_iqr)
r.squaredGLMM(m_month_iqr)


# chamber and EC IQR
m_month_iqr_2 <- lme(
  log(abs(DELTA_FCH4)) ~ log(ch_iqr_permonthsite + 0.01) + log(EC_IQR_permonthsite + 0.01), 
  random = ~1 | SITE, 
  data = month_noNA
)
summary(m_month_iqr_2)
r.squaredGLMM(m_month_iqr_2)

# annual

# chamber IQR

m_yr_iqr <- lme(
  log(abs(DELTA_FCH4)) ~ log(ch_iqr_peryearsite + 0.01), 
  random = ~1 | SITE, 
  data = all_yr_aggr
)

summary(m_yr_iqr)
r.squaredGLMM(m_yr_iqr)

# chamber and EC IQR

m_yr_iqr_2 <- lme(
  log(abs(DELTA_FCH4)) ~ log(ch_iqr_peryearsite + 0.01) + log(EC_IQR_peryearsite + 0.01), 
  random = ~1 | SITE, 
  data = yr_noNA
)
summary(m_yr_iqr_2)
r.squaredGLMM(m_yr_iqr_2)

# EC IQR

m_yr_iqr <- lme(
  log(abs(DELTA_FCH4)) ~ log(EC_IQR_peryearsite + 0.01), 
  random = ~1 | SITE, 
  data = all_yr_aggr
)

summary(m_yr_iqr)
r.squaredGLMM(m_yr_iqr)


# weekly plots

# chamber IQR
# first with all sites, untransformed

week_ch_IQR_untransformed <- ggplot(all_week_aggr, 
                                    aes(x = ch_iqr_perweeksite, 
                                        y = abs(DELTA_FCH4), 
                                        color = SITE)) +
  geom_point(size = 3) + 
  # extra layer: US-UAF
  geom_point(data = filter(all_week_aggr, SITE == "US-UAF"),
             aes(x = ch_iqr_perweeksite, y = abs(DELTA_FCH4)),
             shape = 21, fill = "white", color = "black", size = 3, stroke = 1) +
  labs(x = expression("Chamber FCH" [4] * " IQR"), 
       y = expression(paste("Absolute "*Delta*"FCH"[4], 
                            " nmol m"^-2, " s"^-1))) +
  stat_cor(method = "spearman", show.legend = FALSE, size = 5) +
  scale_color_okabeito() +
  theme_bw() +
  theme(axis.title = element_text(size = 20),
        axis.text = element_text(size = 20))

week_ch_IQR_untransformed

# EC IQR
# first all sites, untransformed

week_ec_IQR_untransformed <- ggplot(all_week_aggr, 
                                    aes(x = EC_IQR_perweeksite, 
                                        y = abs(DELTA_FCH4), 
                                        color = SITE)) +
  geom_point(size = 3) +
  # extra layer: US-UAF
  geom_point(data = filter(all_week_aggr, SITE == "US-UAF"),
             aes(x = EC_IQR_perweeksite, y = abs(DELTA_FCH4)),
             shape = 21, fill = "white", color = "black", size = 3, stroke = 1) +
  labs(x = expression("EC FCH" [4] * " IQR (nmol m"^-2, " s"^-1, ")"), 
       y = expression(paste("Absolute "*Delta*"FCH"[4], 
                            " nmol m"^-2, " s"^-1))) +
  stat_cor(method = "spearman", show.legend = FALSE, size = 5) +
  scale_color_okabeito() +
  theme_bw() +
  theme(axis.title = element_text(size = 20),
        axis.text = element_text(size = 20))

week_ec_IQR_untransformed


# check spearman correlations

cor.test(all_week_aggr$ch_iqr_perweeksite, abs(all_week_aggr$DELTA_FCH4), method = "spearman")
cor.test(all_week_aggr$EC_IQR_perweeksite, abs(all_week_aggr$DELTA_FCH4), method = "spearman")

# chamber IQR
all_week_aggr %>%
  group_by(SITE) %>%
  summarize(correlation = cor.test(ch_iqr_perweeksite.x, abs(DELTA_FCH4), method = "spearman")$estimate,
            p_value = cor.test(ch_iqr_perweeksite.x, abs(DELTA_FCH4), method = "spearman")$p.value,
            count = n())

# EC IQR
all_week_aggr %>%
  group_by(SITE) %>%
  summarize(correlation = cor.test(EC_IQR_perweeksite, abs(DELTA_FCH4), method = "spearman")$estimate,
            p_value = cor.test(EC_IQR_perweeksite, abs(DELTA_FCH4), method = "spearman")$p.value,
            count = n())


# log-transformed

# chamber IQR

# sites separately, log-transformed

ggplot(all_week_aggr, 
       aes(x = log(ch_iqr_perweeksite + 0.01), 
           y = log(abs(DELTA_FCH4)), 
           color = SITE, shape = SITE)) +
  geom_point(size = 2) +
  # US-UAF:
  geom_point(data = filter(all_week_aggr, SITE == "US-UAF"),
             aes(x = log(ch_iqr_perweeksite + 0.01), 
                 y = log(abs(DELTA_FCH4))),
             shape = 21, fill = "white", color = "black", size = 1, stroke = 0.8) +
  labs(x = expression("Chamber FCH" [4] * " IQR (nmol m"^-2, " s"^-1, ", log +0.01)"), 
       y = expression(paste("Absolute "*Delta*"FCH"[4], 
                            " nmol m"^-2, " s"^-1, ", log"))) +
  geom_smooth(method = "lm", se = TRUE, alpha = 0.2) +
  geom_smooth(data = filter(all_week_aggr, SITE == "US-UAF"),
              aes(x = log(ch_iqr_perweeksite + 0.01), 
                  y = log(abs(DELTA_FCH4))),
              method = "lm", se = TRUE, alpha = 0.2, color = "black") +
  stat_poly_eq(aes(label = paste(..eq.label.., ..rr.label.., ..p.value.label.., sep = "~~~"), color = SITE),
               formula = y ~ x,
               parse = TRUE,
               label.x.npc = "left",
               label.y.npc = "top",
               size = 5) +
  scale_color_manual(values = my_cols) +
  scale_shape_manual(values = my_shapes) +
  theme_bw() +
  theme(axis.title = element_text(size=20),
        axis.text = element_text(size = 20),
        legend.text = element_text(size = 15))



# per-site regression stats
df <- all_week_aggr %>%
  transmute(
    SITE,
    x = log(ch_iqr_perweeksite + 0.01),
    y = log(abs(DELTA_FCH4))
  ) %>%
  filter(is.finite(x), is.finite(y))

site_stats <- df %>%
  group_by(SITE) %>%
  group_modify(~{
    d <- .x
    n <- nrow(d)
    if (n < 3 || var(d$x) == 0 || var(d$y) == 0) {
      return(tibble(
        n = n, r2 = NA_real_, adj_r2 = NA_real_,
        slope = NA_real_, slope_p = NA_real_, model_p = NA_real_,
        note = "Insufficient data or zero variance"
      ))
    }
    fit <- lm(y ~ x, data = d)
    s   <- summary(fit)
    slope_p <- coef(s)["x","Pr(>|t|)"]
    model_p <- if (!is.null(s$fstatistic)) {
      pf(s$fstatistic[1], s$fstatistic[2], s$fstatistic[3], lower.tail = FALSE)
    } else NA_real_
    tibble(
      n = n,
      r2 = s$r.squared,          
      adj_r2 = s$adj.r.squared,  
      slope = coef(fit)[["x"]],
      slope_p = slope_p, 
      model_p = model_p,
      note = NA_character_
    )
  }) %>%
  ungroup()

# See the table
site_stats %>%
  mutate(across(c(r2, adj_r2, slope), ~round(.x, 3)),
         across(c(slope_p, model_p), ~signif(.x, 3))) %>%
  arrange(desc(r2)) %>%
  print(n = Inf)

# color only significant predictors
all_week_aggr$site_color_ch <- ifelse(all_week_aggr$SITE %in% c("US-LA1","US-OWC"), all_week_aggr$SITE, "Other")

# Reorder data
all_week_aggr <- all_week_aggr %>%
  arrange(
    ifelse(site_color_ch == "SE-DEG", 3, 
           ifelse(site_color_ch == "US-OWC", 2, 
                  ifelse(site_color_ch == "US-LA1", 1, 0)))
  )


# plot with only significant sites
week_IQR_ch <- ggplot(all_week_aggr, 
                      aes(x = log(ch_iqr_perweeksite + 0.01), 
                          y = log(abs(DELTA_FCH4)))) +
  geom_point(aes(color = site_color_ch, shape = site_color_ch), size = 2, stroke = 1) +
  geom_hline(yintercept = 0, color = "blue", linetype = "dashed", linewidth = 1) +
  labs(x = expression("Chamber FCH" [4] * " IQR (nmol m"^{-2}*" s"^{-1}*", log +0.01)"), 
       y = expression("Absolute "*Delta*"FCH"[4]*" (nmol m"^{-2}*" s"^{-1}*", log)")) +
  geom_smooth(aes(x = log(ch_iqr_perweeksite + 0.01), y = log(abs(DELTA_FCH4)), color = site_color_ch), 
              method = "lm",  
              se = F) +  # Add overall regression line
  stat_poly_eq(aes(label = paste(..eq.label.., ..rr.label.., ..p.value.label.., sep = "~~~"), color = site_color_ch), 
               formula = y ~ x, 
               parse = TRUE, 
               label.x.npc = "left", 
               label.y.npc = "top", 
               size = 5) + 
  scale_color_manual(values = c("US-LA1" = "#0072B2", "US-OWC" = "#CC79A7", "Other" = "lightgray")) +
  scale_shape_manual(values = c("US-LA1" = 4,
                                "US-OWC" = 6,
                                "Other" = 16)) +
  theme_bw() +
  theme(axis.title = element_text(size = 20),
        axis.text = element_text(size = 20),
        legend.text = element_text(size = 15),
        panel.background = element_rect(fill = "transparent", color = NA),
        plot.background  = element_rect(fill = "transparent", color = NA),
        legend.background = element_rect(fill = "transparent", color = NA),
        legend.box.background = element_rect(fill = "transparent", color = NA))

week_IQR_ch


# EC IQR, log-transformed, all sites

ggplot(all_week_aggr, 
       aes(x = log(EC_IQR_perweeksite + 0.01), 
           y = log(abs(DELTA_FCH4)), 
           color = SITE, shape = SITE)) +
  geom_point(size = 2) +
  # US-UAF:
  geom_point(data = filter(all_week_aggr, SITE == "US-UAF"),
             aes(x = log(EC_IQR_perweeksite + 0.01), 
                 y = log(abs(DELTA_FCH4))),
             shape = 21, fill = "white", color = "black", size = 1, stroke = 0.8) +
  labs(x = expression("EC FCH" [4] * " IQR (nmol m"^-2, " s"^-1, ", log +0.01)"), 
       y = expression(paste("Log absolute FCH"[4], 
                            " nmol m"^-2, " s"^-1, ", log"))) +
  geom_smooth(method = "lm", se = TRUE, alpha = 0.2) +
  geom_smooth(data = filter(all_week_aggr, SITE == "US-UAF"),
              aes(x = log(EC_IQR_perweeksite + 0.01), 
                  y = log(abs(DELTA_FCH4))),
              method = "lm", se = TRUE, alpha = 0.2, color = "black") +
  stat_poly_eq(aes(label = paste(..eq.label.., ..rr.label.., ..p.value.label.., sep = "~~~"), color = SITE),
               formula = y ~ x,
               parse = TRUE,
               label.x.npc = "left",
               label.y.npc = "top",
               size = 5) +
  scale_color_manual(values = my_cols) +
  scale_shape_manual(values = my_shapes) +
  theme_bw() +
  theme(axis.title = element_text(size=20),
        axis.text = element_text(size = 20),
        legend.text = element_text(size = 15))



# per-site regression stats
df <- all_week_aggr %>%
  transmute(
    SITE,
    x = log(EC_IQR_perweeksite + 0.01),
    y = log(abs(DELTA_FCH4))
  ) %>%
  filter(is.finite(x), is.finite(y))

site_stats <- df %>%
  group_by(SITE) %>%
  group_modify(~{
    d <- .x
    n <- nrow(d)
    if (n < 3 || var(d$x) == 0 || var(d$y) == 0) {
      return(tibble(
        n = n, r2 = NA_real_, adj_r2 = NA_real_,
        slope = NA_real_, slope_p = NA_real_, model_p = NA_real_,
        note = "Insufficient data or zero variance"
      ))
    }
    fit <- lm(y ~ x, data = d)
    s   <- summary(fit)
    slope_p <- coef(s)["x","Pr(>|t|)"]
    model_p <- if (!is.null(s$fstatistic)) {
      pf(s$fstatistic[1], s$fstatistic[2], s$fstatistic[3], lower.tail = FALSE)
    } else NA_real_
    tibble(
      n = n,
      r2 = s$r.squared, 
      adj_r2 = s$adj.r.squared,
      slope = coef(fit)[["x"]],
      slope_p = slope_p,   
      model_p = model_p,
      note = NA_character_
    )
  }) %>%
  ungroup()

# See the table
site_stats %>%
  mutate(across(c(r2, adj_r2, slope), ~round(.x, 3)),
         across(c(slope_p, model_p), ~signif(.x, 3))) %>%
  arrange(desc(r2)) %>%
  print(n = Inf)


# color only significant predictors
all_week_aggr$site_color_ec <- ifelse(all_week_aggr$SITE %in% c("SE-DEG"), all_week_aggr$SITE, "Other")

# Reorder data so that SE-DEG is last
all_week_aggr <- all_week_aggr %>%
  arrange(ifelse(site_color_ec == "SE-DEG", 1, 0))

# EC IQR log-transformed, with only SE-DEG colored
week_IQR_ec <- ggplot(all_week_aggr, 
                      aes(x = log(EC_IQR_perweeksite + 0.01), 
                          y = log(abs(DELTA_FCH4)))) +
  geom_point(aes(color = site_color_ec, shape = site_color_ec), size = 2, stroke = 1) +
  geom_hline(yintercept = 0, color = "blue", linetype = "dashed", linewidth = 1) +
  labs(x = expression("EC FCH" [4] * " IQR (nmol m"^{-2}*" s"^{-1}*", log +0.01) "),  
       y = expression("Log absolute "*Delta*"FCH"[4]*" (nmol m"^{-2}*" s"^{-1}*", log)")) +
  
  # regression line only for SE-DEG
  geom_smooth(data = dplyr::filter(all_week_aggr, site_color_ec == "SE-DEG"),
              aes(color = site_color_ec),
              method = "lm", se = FALSE, linewidth = 1) +
  
  stat_poly_eq(aes(label = paste(..eq.label.., ..rr.label.., ..p.value.label.., sep = "~~~"), 
                   color = site_color_ec), 
               formula = y ~ x, 
               parse = TRUE, 
               label.x.npc = "left", 
               label.y.npc = "top", 
               size = 5) +
  scale_color_manual(values = c("SE-DEG" = "#009E73", "Other" = "lightgray")) +
  scale_shape_manual(values = c("SE-DEG" = 0, "Other" = 16)) +
  theme_bw() +
  theme(axis.title = element_text(size = 20),
        axis.text = element_text(size = 20),
        legend.text = element_text(size = 15),
        panel.background = element_rect(fill = "transparent", color = NA),
        plot.background  = element_rect(fill = "transparent", color = NA),
        legend.background = element_rect(fill = "transparent", color = NA),
        legend.box.background = element_rect(fill = "transparent", color = NA))
week_IQR_ec


# monthly plots

# spearman correlations
cor.test(all_month_aggr$ch_iqr_permonthsite, abs(all_month_aggr$DELTA_FCH4), method = "spearman")
cor.test(all_month_aggr$EC_IQR_permonthsite, abs(all_month_aggr$DELTA_FCH4), method = "spearman")

# per site
all_month_aggr %>%
  group_by(SITE) %>%
  summarize(correlation = cor.test(ch_iqr_permonthsite, abs(DELTA_FCH4), method = "spearman")$estimate,
            p_value = cor.test(ch_iqr_permonthsite, abs(DELTA_FCH4), method = "spearman")$p.value,
            count = n())

all_month_aggr %>%
  group_by(SITE) %>%
  summarize(correlation = cor.test(EC_IQR_permonthsite, abs(DELTA_FCH4), method = "spearman")$estimate,
            p_value = cor.test(EC_IQR_permonthsite, abs(DELTA_FCH4), method = "spearman")$p.value,
            count = n())

# arrange sites for plotting
all_month_aggr <- all_month_aggr %>%
  arrange(ifelse(site_color_ch == "US-LA1", 1, 0))
all_month_aggr <- all_month_aggr %>%
  arrange(ifelse(site_color_ec == "SE-DEG", 1, 0))


# EC IQR all sites, untransformed
month_ec_IQR_untransformed <- ggplot(all_month_aggr, 
                                     aes(x = EC_IQR_permonthsite, 
                                         y = abs(DELTA_FCH4), 
                                         color = SITE)) +
  geom_point(size = 3) +
  # US-UAF:
  geom_point(data = filter(all_month_aggr, SITE == "US-UAF"),
             aes(x = EC_IQR_permonthsite, y = abs(DELTA_FCH4)),
             shape = 21, fill = "white", color = "black", size = 3, stroke = 1) +
  labs(x = expression("Chamber FCH" [4] * " IQR (nmol m"^-2, " s"^-1, ")"), 
       y = expression(paste("Absolute "*Delta*"FCH"[4], 
                            " nmol m"^-2, " s"^-1))) +
  stat_cor(method = "spearman", show.legend = FALSE, size = 5) +
  scale_color_okabeito() +
  theme_bw() +
  theme(axis.title = element_text(size = 20),
        axis.text = element_text(size = 20))
month_ec_IQR_untransformed


# chamber IQR, all sites, log-transformed
ggplot(all_month_aggr, 
       aes(x = log(ch_iqr_permonthsite + 0.01), 
           y = log(abs(DELTA_FCH4)), 
           color = SITE, shape = SITE)) +
  geom_point(size = 2) +
  # US-UAF:
  geom_point(data = filter(all_month_aggr, SITE == "US-UAF"),
             aes(x = log(ch_iqr_permonthsite + 0.01), 
                 y = log(abs(DELTA_FCH4))),
             shape = 21, fill = "white", color = "black", size = 1, stroke = 0.8) +
  labs(x = expression("Chamber FCH" [4] * " IQR (nmol m"^-2, " s"^-1, ", log +0.01)"), 
       y = expression(paste("Log absolute "*Delta*"FCH"[4], 
                            " nmol m"^-2, " s"^-1))) +
  geom_smooth(method = "lm", se = TRUE, alpha = 0.2) +
  geom_smooth(data = filter(all_month_aggr, SITE == "US-UAF"),
              aes(x = log(ch_iqr_permonthsite + 0.01), 
                  y = log(abs(DELTA_FCH4))),
              method = "lm", se = TRUE, alpha = 0.2, color = "black") +
  stat_poly_eq(aes(label = paste(..eq.label.., ..rr.label.., ..p.value.label.., sep = "~~~"), color = SITE),
               formula = y ~ x,
               parse = TRUE,
               label.x.npc = "left",
               label.y.npc = "top",
               size = 5) +
  scale_color_manual(values = my_cols) +
  scale_shape_manual(values = my_shapes) +
  theme_bw() +
  theme(axis.title = element_text(size=20),
        axis.text = element_text(size = 20),
        legend.text = element_text(size = 15))

# per-site regression stats
df <- all_month_aggr %>%
  transmute(
    SITE,
    x = log(ch_iqr_permonthsite + 0.01),
    y = log(abs(DELTA_FCH4))
  ) %>%
  filter(is.finite(x), is.finite(y))

site_stats <- df %>%
  group_by(SITE) %>%
  group_modify(~{
    d <- .x
    n <- nrow(d)
    if (n < 3 || var(d$x) == 0 || var(d$y) == 0) {
      return(tibble(
        n = n, r2 = NA_real_, adj_r2 = NA_real_,
        slope = NA_real_, slope_p = NA_real_, model_p = NA_real_,
        note = "Insufficient data or zero variance"
      ))
    }
    fit <- lm(y ~ x, data = d)
    s   <- summary(fit)
    slope_p <- coef(s)["x","Pr(>|t|)"]
    model_p <- if (!is.null(s$fstatistic)) {
      pf(s$fstatistic[1], s$fstatistic[2], s$fstatistic[3], lower.tail = FALSE)
    } else NA_real_
    tibble(
      n = n,
      r2 = s$r.squared,        
      adj_r2 = s$adj.r.squared,
      slope = coef(fit)[["x"]],
      slope_p = slope_p, 
      model_p = model_p,
      note = NA_character_
    )
  }) %>%
  ungroup()

# See the table
site_stats %>%
  mutate(across(c(r2, adj_r2, slope), ~round(.x, 3)),
         across(c(slope_p, model_p), ~signif(.x, 3))) %>%
  arrange(desc(r2)) %>%
  print(n = Inf)

# color only significant predictors
all_month_aggr$site_color_ch <- ifelse(all_month_aggr$SITE %in% c("US-LA1", "US-OWC"), all_month_aggr$SITE, "Other")

# plot
month_IQR_ch <- ggplot(all_month_aggr, 
                       aes(x = log(ch_iqr_permonthsite + 0.01), 
                           y = log(abs(DELTA_FCH4)))) +
  geom_point(aes(color = site_color_ch, shape = site_color_ch), size = 2, stroke = 1) +
  geom_hline(yintercept = 0, color = "blue", linetype = "dashed", linewidth = 1) +
  labs(x = expression("Chamber FCH" [4] * " IQR (nmol m"^-2, " s"^-1, ", log +0.01)"), 
       y = expression("Absolute "*Delta*"FCH"[4] * "(nmol m"^{-2}*" s"^{-1}*", log)")) + 
  geom_smooth(data = dplyr::filter(all_month_aggr, site_color_ch %in% c("US-LA1", "US-LA2")),
              aes(color = site_color_ch),
              method = "lm", se = FALSE, linewidth = 1) +
  # Group-specific regression equations without lines
  stat_poly_eq(aes(label = paste(..eq.label.., ..rr.label.., ..p.value.label.., sep = "~~~"), color = site_color_ch), 
               formula = y ~ x, 
               parse = TRUE, 
               label.x.npc = "left", 
               label.y.npc = "top", 
               size = 5) + 
  scale_color_manual(values = c("US-LA1" = "#0072B2", "US-OWC" = "#CC79A7", "Other" = "lightgray")) +
  # Specify shapes for chosen sites and circle (16) for "Other"
  scale_shape_manual(values = c("US-LA1" = 4,
                                "US-OWC" = 6,
                                "Other" = 16)) +
  theme_bw() +
  theme(axis.title = element_text(size = 20),
        axis.text = element_text(size = 20),
        legend.text = element_text(size = 15),
        panel.background = element_rect(fill = "transparent", color = NA),
        plot.background  = element_rect(fill = "transparent", color = NA),
        legend.background = element_rect(fill = "transparent", color = NA),
        legend.box.background = element_rect(fill = "transparent", color = NA))

month_IQR_ch


# EC IQR log-transformed, all sites

# Since some sites do not have enough data points, pick groups that can fit a line (>=3 finite pts, non-zero variance)
df_labels <- all_month_aggr %>%
  mutate(
    x = log(EC_IQR_permonthsite + 0.01),
    y = log(abs(DELTA_FCH4))
  ) %>%
  group_by(SITE) %>%
  filter(
    sum(is.finite(x) & is.finite(y)) >= 3,
    sd(x[is.finite(x)], na.rm = TRUE) > 0,
    sd(y[is.finite(y)], na.rm = TRUE) > 0
  ) %>%
  ungroup()

# plot
ggplot(all_month_aggr, 
       aes(x = log(EC_IQR_permonthsite + 0.01), 
           y = log(abs(DELTA_FCH4)), 
           color = SITE, shape = SITE)) +
  geom_point(size = 2) +
  # US-UAF:
  geom_point(data = filter(all_month_aggr, SITE == "US-UAF"),
             aes(x = log(EC_IQR_permonthsite + 0.01), 
                 y = log(abs(DELTA_FCH4))),
             shape = 21, fill = "white", color = "black", size = 1, stroke = 0.8) +
  labs(x = expression("EC FCH" [4] * " IQR (nmol m"^-2, " s"^-1, ", log +0.01)"), 
       y = expression(paste("Log absolute "*Delta*"FCH"[4], 
                            " nmol m"^-2, " s"^-1))) +
  geom_smooth(method = "lm", se = TRUE, alpha = 0.2) +
  geom_smooth(data = filter(all_month_aggr, SITE == "US-UAF"),
              aes(x = log(EC_IQR_permonthsite + 0.01), 
                  y = log(abs(DELTA_FCH4))),
              method = "lm", se = TRUE, alpha = 0.2, color = "black") +
  ggpmisc::stat_poly_eq(
    data = df_labels,
    aes(x = x, y = y, group = SITE,
        label = paste(..eq.label.., ..rr.label.., ..p.value.label.., sep = "~~~")),
    formula = y ~ x,
    parse = TRUE,
    label.x.npc = "left",
    label.y.npc = "top",
    size = 5,
    color = "black",
    na.rm = TRUE
  ) +
  scale_color_manual(values = my_cols) +
  scale_shape_manual(values = my_shapes) +
  theme_bw() +
  theme(axis.title = element_text(size=20),
        axis.text = element_text(size = 20),
        legend.text = element_text(size = 15))


# per-site regression stats
df <- all_month_aggr %>%
  transmute(
    SITE,
    x = log(EC_IQR_permonthsite + 0.01),
    y = log(abs(DELTA_FCH4))
  ) %>%
  filter(is.finite(x), is.finite(y))

site_stats <- df %>%
  group_by(SITE) %>%
  group_modify(~{
    d <- .x
    n <- nrow(d)
    if (n < 3 || var(d$x) == 0 || var(d$y) == 0) {
      return(tibble(
        n = n, r2 = NA_real_, adj_r2 = NA_real_,
        slope = NA_real_, slope_p = NA_real_, model_p = NA_real_,
        note = "Insufficient data or zero variance"
      ))
    }
    fit <- lm(y ~ x, data = d)
    s   <- summary(fit)
    slope_p <- coef(s)["x","Pr(>|t|)"]
    model_p <- if (!is.null(s$fstatistic)) {
      pf(s$fstatistic[1], s$fstatistic[2], s$fstatistic[3], lower.tail = FALSE)
    } else NA_real_
    tibble(
      n = n,
      r2 = s$r.squared,           
      adj_r2 = s$adj.r.squared, 
      slope = coef(fit)[["x"]],
      slope_p = slope_p, 
      model_p = model_p,
      note = NA_character_
    )
  }) %>%
  ungroup()

# See the table
site_stats %>%
  mutate(across(c(r2, adj_r2, slope), ~round(.x, 3)),
         across(c(slope_p, model_p), ~signif(.x, 3))) %>%
  arrange(desc(r2)) %>%
  print(n = Inf)

# color only significant sites
all_month_aggr$site_color_ec <- ifelse(all_month_aggr$SITE %in% c("SE-DEG"), all_month_aggr$SITE, "Other")

all_month_aggr <- all_month_aggr %>%
  arrange(ifelse(site_color_ec == "SE-DEG", 1, 0))

# plot only with significant sites
month_IQR_ec <- ggplot(all_month_aggr, 
                       aes(x = log(EC_IQR_permonthsite + 0.01), 
                           y = log(abs(DELTA_FCH4)))) +
  geom_point(aes(color = site_color_ec, shape = site_color_ec), size = 2, stroke = 1) + 
  geom_hline(yintercept = 0, color = "blue", linetype = "dashed", linewidth = 1) +
  labs(x = expression("EC FCH" [4] * " IQR (nmol m"^{-2}*" s"^{-1}*", log +0.01) "), 
       y = expression("Absolute "*Delta*"FCH"[4] * "(nmol m"^{-2}*" s"^{-1}*", log)")) +
  geom_smooth(aes(x = log(EC_IQR_permonthsite + 0.01), y = log(abs(DELTA_FCH4)), color = site_color_ec), 
              method = "lm", se = F) +
  stat_poly_eq(aes(label = paste(..eq.label.., ..rr.label.., ..p.value.label.., sep = "~~~"), color = site_color_ec), 
               formula = y ~ x, 
               parse = TRUE, 
               label.x.npc = "left", 
               label.y.npc = "bottom", 
               size = 5) + 
  scale_color_manual(values = c("SE-DEG" = "#009E73", "Other" = "lightgray")) +
  scale_shape_manual(values = c("SE-DEG" = 0,
                                "Other" = 16)) +
  theme_bw() +
  theme(axis.title = element_text(size = 20),
        axis.text = element_text(size = 20),
        legend.text = element_text(size = 15),
        panel.background = element_rect(fill = "transparent", color = NA),
        plot.background  = element_rect(fill = "transparent", color = NA),
        legend.background = element_rect(fill = "transparent", color = NA),
        legend.box.background = element_rect(fill = "transparent", color = NA))

month_IQR_ec

# annual plots

# all sites, untransformed
# chamber IQR
year_ch_IQR_untransformed <- ggplot(all_yr_aggr, 
                                    aes(x = ch_iqr_peryearsite, 
                                        y = abs(DELTA_FCH4), 
                                        color = SITE)) +
  geom_point(size = 3) +
  # US-UAF:
  geom_point(data = filter(all_yr_aggr, SITE == "US-UAF"),
             aes(x = ch_iqr_peryearsite, y = abs(DELTA_FCH4)),
             shape = 21, fill = "white", color = "black", size = 3, stroke = 1) +
  labs(x = expression("Chamber FCH" [4] * " IQR (nmol m"^-2, " s"^-1, ")"), 
       y = expression(paste("Absolute "*Delta*"FCH"[4], 
                            " nmol m"^-2, " s"^-1))) +
  stat_cor(method = "spearman", show.legend = FALSE, size = 5) +
  scale_color_okabeito() +
  theme_bw() +
  theme(axis.title = element_text(size = 20),
        axis.text = element_text(size = 20))

year_ch_IQR_untransformed

# all sites, untransformed
# EC IQR

year_ec_IQR_untransformed <- ggplot(all_yr_aggr, 
                                    aes(x = EC_IQR_peryearsite, 
                                        y = abs(DELTA_FCH4), 
                                        color = SITE)) +
  geom_point(size = 3) + 
  # US-UAF:
  geom_point(data = filter(all_yr_aggr, SITE == "US-UAF"),
             aes(x = EC_IQR_peryearsite, y = abs(DELTA_FCH4)),
             shape = 21, fill = "white", color = "black", size = 3, stroke = 1) +
  labs(x = expression("EC FCH" [4] * " IQR (nmol m"^-2, " s"^-1, ")"), 
       y = expression(paste("Absolute "*Delta*"FCH"[4], 
                            " nmol m"^-2, " s"^-1))) +
  stat_cor(method = "spearman", show.legend = FALSE, size = 5) +
  scale_color_okabeito() +
  theme_bw() +
  theme(axis.title = element_text(size = 20),
        axis.text = element_text(size = 20))
year_ec_IQR_untransformed

# spearman correlations
cor.test(all_yr_aggr$ch_iqr_peryearsite, abs(all_yr_aggr$DELTA_FCH4), method = "spearman")
cor.test(all_yr_aggr$EC_IQR_peryearsite, abs(all_yr_aggr$DELTA_FCH4), method = "spearman")

# per site
all_yr_aggr %>%
  filter(SITE != c("US-LA1", "US-LOS")) %>% # not enough points
  group_by(SITE) %>%
  summarize(correlation = cor.test(ch_iqr_peryearsite, abs(DELTA_FCH4), method = "spearman")$estimate,
            p_value = cor.test(ch_iqr_peryearsite, abs(DELTA_FCH4), method = "spearman")$p.value,
            count = n())

all_yr_aggr %>%
  filter(SITE != c("US-LA1", "US-LOS")) %>% # not enough points
  group_by(SITE) %>%
  summarize(correlation = cor.test(EC_IQR_peryearsite, abs(DELTA_FCH4), method = "spearman")$estimate,
            p_value = cor.test(EC_IQR_peryearsite, abs(DELTA_FCH4), method = "spearman")$p.value,
            count = n())

# chamber IQR, log-transformed

ggplot(all_yr_aggr, 
       aes(x = log(ch_iqr_peryearsite + 0.01), 
           y = log(abs(DELTA_FCH4)), 
           color = SITE, shape = SITE)) +
  geom_point(size = 2) +
  # US-UAF:
  geom_point(data = filter(all_yr_aggr, SITE == "US-UAF"),
             aes(x = log(ch_iqr_peryearsite + 0.01), 
                 y = log(abs(DELTA_FCH4))),
             shape = 21, fill = "white", color = "black", size = 1, stroke = 0.8) +
  labs(x = expression("Chamber FCH" [4] * " IQR (nmol m"^-2, " s"^-1, "log + 0.01) "), 
       y = expression(paste("Absolute FCH"[4] *  
                            " nmol m"^-2, " s"^-1, ", log)"))) +
  geom_smooth(method = "lm", se = TRUE, alpha = 0.2) +
  geom_smooth(data = filter(all_yr_aggr, SITE == "US-UAF"),
              aes(x = log(ch_iqr_peryearsite + 0.01), 
                  y = log(abs(DELTA_FCH4))),
              method = "lm", se = TRUE, alpha = 0.2, color = "black") +
  stat_poly_eq(aes(label = paste(..eq.label.., ..rr.label.., ..p.value.label.., sep = "~~~"), color = SITE),
               formula = y ~ x,
               parse = TRUE,
               label.x.npc = "left",
               label.y.npc = "top",
               size = 5) +
  scale_color_manual(values = my_cols) +
  scale_shape_manual(values = my_shapes) +
  theme_bw() +
  theme(axis.title = element_text(size=20),
        axis.text = element_text(size = 20),
        legend.text = element_text(size = 15))


# per-site regression stats
df <- all_yr_aggr %>%
  transmute(
    SITE,
    x = log(ch_iqr_peryearsite + 0.01),
    y = log(abs(DELTA_FCH4))
  ) %>%
  filter(is.finite(x), is.finite(y))

site_stats <- df %>%
  group_by(SITE) %>%
  group_modify(~{
    d <- .x
    n <- nrow(d)
    if (n < 3 || var(d$x) == 0 || var(d$y) == 0) {
      return(tibble(
        n = n, r2 = NA_real_, adj_r2 = NA_real_,
        slope = NA_real_, slope_p = NA_real_, model_p = NA_real_,
        note = "Insufficient data or zero variance"
      ))
    }
    fit <- lm(y ~ x, data = d)
    s   <- summary(fit)
    slope_p <- coef(s)["x","Pr(>|t|)"]
    model_p <- if (!is.null(s$fstatistic)) {
      pf(s$fstatistic[1], s$fstatistic[2], s$fstatistic[3], lower.tail = FALSE)
    } else NA_real_
    tibble(
      n = n,
      r2 = s$r.squared, 
      adj_r2 = s$adj.r.squared,
      slope = coef(fit)[["x"]],
      slope_p = slope_p,  
      model_p = model_p,
      note = NA_character_
    )
  }) %>%
  ungroup()

# See the table
site_stats %>%
  mutate(across(c(r2, adj_r2, slope), ~round(.x, 3)),
         across(c(slope_p, model_p), ~signif(.x, 3))) %>%
  arrange(desc(r2)) %>%
  print(n = Inf)

# none of the sites were significant
# plot without highlighting individual sites
yr_IQR_ch <- ggplot(all_yr_aggr, aes(x = log(ch_iqr_peryearsite + 0.01), y = log(abs(DELTA_FCH4)))) +
  geom_point(size=2, color = "lightgrey") + 
  geom_hline(yintercept = 0, color = "blue", linetype = "dashed", linewidth = 1) +
  labs(x = expression("Chamber FCH" [4] * " IQR (nmol m"^{-2}*" s"^{-1}*", log +0.01)"), y = expression("Absolute "*Delta*"FCH" [4] * " (nmol m"^{-2}*" s"^{-1}*")")) +
  theme_bw()  +
  theme(axis.title = element_text(size=20),
        axis.text = element_text(size = 20),
        legend.text = element_text(size = 15),
        panel.background = element_rect(fill = "transparent", color = NA),
        plot.background  = element_rect(fill = "transparent", color = NA),
        legend.background = element_rect(fill = "transparent", color = NA),
        legend.box.background = element_rect(fill = "transparent", color = NA))
yr_IQR_ch

# EC IQR log-transformed

ggplot(all_yr_aggr, 
       aes(x = log(EC_IQR_peryearsite + 0.01), 
           y = log(abs(DELTA_FCH4)), 
           color = SITE, shape = SITE)) +
  geom_point(size = 2) +
  # US-UAF:
  geom_point(data = filter(all_yr_aggr, SITE == "US-UAF"),
             aes(x = log(EC_IQR_peryearsite + 0.01), 
                 y = log(abs(DELTA_FCH4))),
             shape = 21, fill = "white", color = "black", size = 1, stroke = 0.8) +
  labs(x = expression("EC FCH" [4] * " IQR (nmol m"^-2, " s"^-1, ", log +0.01)"), 
       y = expression(paste("Absolute "*Delta*"FCH"[4], 
                            " nmol m"^-2, " s"^-1, ", log"))) +
  geom_smooth(method = "lm", se = TRUE, alpha = 0.2) +
  geom_smooth(data = filter(all_yr_aggr, SITE == "US-UAF"),
              aes(x = log(EC_IQR_peryearsite + 0.01), 
                  y = log(abs(DELTA_FCH4))),
              method = "lm", se = TRUE, alpha = 0.2, color = "black") +
  stat_poly_eq(aes(label = paste(..eq.label.., ..rr.label.., ..p.value.label.., sep = "~~~"), color = SITE),
               formula = y ~ x,
               parse = TRUE,
               label.x.npc = "left",
               label.y.npc = "top",
               size = 5) +
  scale_color_manual(values = my_cols) +
  scale_shape_manual(values = my_shapes) +
  theme_bw() +
  theme(axis.title = element_text(size=20),
        axis.text = element_text(size = 20),
        legend.text = element_text(size = 15))

# none of the sites were significant --> plot without highlighting individual sites
yr_IQR_ec <- ggplot(all_yr_aggr, aes(x = log(EC_IQR_peryearsite + 0.01), y = log(abs(DELTA_FCH4)))) +
  geom_point(size=2, color="lightgrey") +  
  geom_hline(yintercept = 0, color = "blue", linetype = "dashed", linewidth = 1) +
  labs(x = expression("EC FCH" [4] * " IQR (nmol m"^{-2}*" s"^{-1}*", log +0.01)"), y = expression("Absolute "*Delta*"FCH" [4] * " (nmol m"^{-2}*" s"^{-1}*")")) + 
  theme_bw()  +
  theme(axis.title = element_text(size=20),
        axis.text = element_text(size = 20),
        legend.text = element_text(size = 15),
        legend.title = element_blank(),
        panel.background = element_rect(fill = "transparent", color = NA),
        plot.background  = element_rect(fill = "transparent", color = NA),
        legend.background = element_rect(fill = "transparent", color = NA),
        legend.box.background = element_rect(fill = "transparent", color = NA))
yr_IQR_ec

### combine all to same figure (Fig. 5) (daily plots excluded)

ch_EC_IQR_reg <- ggarrange(#d_IQR_ch + rremove("ylab") ,
  #   d_IQR_ec + rremove("ylab"), 
  week_IQR_ch + rremove("ylab"), 
  month_IQR_ch + rremove("ylab"),
  yr_IQR_ch + rremove("ylab"),
  week_IQR_ec + rremove("ylab"),
  month_IQR_ec + rremove("ylab") ,
  yr_IQR_ec + rremove("ylab"),
  ncol = 3, nrow = 2,
  align = "hv",
  common.legend = TRUE)

ch_EC_IQR_reg <- annotate_figure(ch_EC_IQR_reg, left = textGrob(expression("Log absolute "*Delta*"FCH"[4] * "(nmol m"^{-2}*" s"^{-1}*")"), rot = 90, vjust = 1, gp = gpar(fontsize = 20)),
                                 bottom = textGrob(""))
ch_EC_IQR_reg


# combine the untransformed plots (Fig. B18) (daily plots included)

ch_EC_IQR_untransformed <- ggarrange(d_ch_IQR_untransformed + rremove("ylab") + rremove("xlab"),
                                     d_ec_IQR_untransformed + rremove("ylab")+ rremove("xlab"), 
                                     week_ch_IQR_untransformed + rremove("ylab")+ rremove("xlab"), 
                                     week_ec_IQR_untransformed + rremove("ylab")+ rremove("xlab"), 
                                     month_ch_IQR_untransformed + rremove("ylab")+ rremove("xlab"),
                                     month_ec_IQR_untransformed + rremove("ylab")+ rremove("xlab"),
                                     year_ch_IQR_untransformed + rremove("ylab"),
                                     year_ec_IQR_untransformed + rremove("ylab"),
                                     ncol = 2, nrow = 4,
                                     align = "hv",
                                     common.legend = TRUE)

ch_EC_IQR_untransformed <- annotate_figure(ch_EC_IQR_untransformed, left = textGrob(expression(paste("Absolute "*Delta*"FCH"[4]," (nmol m"^-2, "s"^-1, ")")), rot = 90, vjust = 1, gp = gpar(fontsize = 20)),
                                           bottom = textGrob(""))
ch_EC_IQR_untransformed


######################################################################################################################################################################

# make plots showing both chamber and EC flux at the same time (Figs. B8-B12)

# first for each site separately

# daily aggregation

all_d_aggr <- all_d_aggr %>%
  mutate(ECCH_median = pmap_dbl(dplyr::select(., EC_FCH4_MEDIAN, CH_FCH4_MEDIAN), ~ median(c(...), na.rm = TRUE)))

# Reshape the data to long format

data_long <- all_d_aggr %>%
  pivot_longer(
    cols = c(CH_FCH4_MEDIAN, EC_FCH4_MEDIAN), 
    names_to = "flux_type", 
    values_to = "FCH4"
  ) %>%
  mutate(flux_type = ifelse(flux_type == "CH_FCH4_MEDIAN", "Chamber", "EC"))

data_long <- as.data.frame(data_long)

# half-hourly dataset
data_long_auto <- all_hh_aggr %>%
  pivot_longer(
    cols = c(CH_FCH4_MEDIAN, EC_FCH4), 
    names_to = "flux_type", 
    values_to = "FCH4"
  ) %>%
  mutate(flux_type = ifelse(flux_type == "CH_FCH4_MEDIAN", "Chamber", "EC"))

data_long_auto <- as.data.frame(data_long_auto)

# half-hourly, with hour, month and per site

data_long_hh_summary_site <- all_hh_aggr %>%
  pivot_longer(
    cols = c(CH_FCH4_MEDIAN, EC_FCH4), 
    names_to = "flux_type", 
    values_to = "FCH4"
  ) %>%
  mutate(flux_type = ifelse(flux_type == "CH_FCH4_MEDIAN", "Chamber", "EC")) %>%
  group_by(HOUR, MONTH, flux_type, SITE) %>%
  summarize(
    median_FCH4 = median(FCH4, na.rm = TRUE),
    IQR_lower = quantile(FCH4, 0.25, na.rm = TRUE),
    IQR_upper = quantile(FCH4, 0.75, na.rm = TRUE),
    median_DELTA_FCH4 = median(DELTA_FCH4, na.rm=TRUE),
    .groups = "drop"
  ) %>%
  mutate(IQR = IQR_upper - IQR_lower)

# half-hourly with hour, month, across sites
data_long_hh_summary <- all_hh_aggr %>%
  pivot_longer(
    cols = c(CH_FCH4_MEDIAN, EC_FCH4), 
    names_to = "flux_type", 
    values_to = "FCH4"
  ) %>%
  mutate(flux_type = ifelse(flux_type == "CH_FCH4_MEDIAN", "Chamber", "EC")) %>%
  group_by(HOUR, MONTH, flux_type) %>%
  summarize(
    median_FCH4 = median(FCH4, na.rm = TRUE),
    IQR_lower = quantile(FCH4, 0.25, na.rm = TRUE),
    IQR_upper = quantile(FCH4, 0.75, na.rm = TRUE),
    median_DELTA_FCH4 = median(DELTA_FCH4, na.rm=TRUE),
    .groups = "drop"
  ) %>%
  mutate(IQR = IQR_upper - IQR_lower)

# exclude outliers from us-ho1
data_long_filtered <- data_long %>%
  filter(!(SITE == "US-HO1" & flux_type == "Chamber" & FCH4 > 3.5),  # Remove rows where Chamber FCH4 > 15 for US-HO1
         !(SITE == "US-HO1" & DELTA_FCH4 < -3.5)) 

data_long_filtered <- data_long %>%
  filter(!(SITE == "US-HO1" & flux_type == "Chamber" & FCH4 > 7.5)) 

# plot site trends (from daily aggregation) (Fig. B14)

ggplot(data_long_filtered, aes(x = TIMESTAMP, y = FCH4, color = flux_type)) +
  
  geom_point(aes(size = case_when(
    SITE == "US-HO1" ~ 1.5,
    SITE %in% c("US-UAF", "SE-DEG", "CN-HGU") ~ 3,
    TRUE ~ 4
  ))) +
  
  # selected sites (alpha = 0.9)
geom_point(
  data = subset(data_long_filtered, 
                SITE %in% c("US-HO1", "CN-HGU", "SE-DEG", "US-UAF")),
  aes(
    y = DELTA_FCH4,
    size = case_when(
      SITE == "US-HO1" ~ 1.5,
      SITE %in% c("US-UAF", "SE-DEG", "CN-HGU") ~ 2.5,
      TRUE ~ 4
    )
  ),
  alpha = 0.4,
  color = "#555555",
  shape = 2,
  stroke = 0.7
) +
  
  # all other sites (alpha = 1)
geom_point(
  data = subset(data_long_filtered, 
                !SITE %in% c("US-HO1", "CN-HGU", "SE-DEG", "US-UAF")),
  aes(
    y = DELTA_FCH4,
    size = case_when(
      SITE == "US-HO1" ~ 1.5,
      SITE %in% c("US-UAF", "SE-DEG", "CN-HGU") ~ 2.5,
      TRUE ~ 4
    )
  ),
  alpha = 1,
  color = "#555555",
  shape = 2,
  stroke = 0.7
) +
  
  labs(x = " ", 
       y = expression("FCH"[4] * " (nmol m"^-2 * " s"^-1 * ")")) +
  
  scale_color_manual(values = c("Chamber" = "#CD4071FF", "EC" = "#0096FF")) +
  scale_shape_manual(values = c("Chamber" = 16, "EC" = 17)) +
  scale_size_identity() +
  scale_x_date(date_labels = "%Y-%m", date_breaks = "1 month") +
  
  facet_wrap(~ SITE, scales = "free", ncol = 2) +
  
  theme_bw() +
  theme(axis.title = element_text(size = 20),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text = element_text(size = 20),
        legend.title = element_blank(),
        legend.text = element_text(size = 20),
        strip.background = element_blank(),
        strip.text = element_text(color = "black", size = 15))

##COOLER MONTH PLOTS FOR CN-HGU (Fig. B15)

data_subset <- data_long %>%
  filter(SITE == "CN-HGU",
         Month_letter %in% c("Feb", "Mar", "Apr", "May", "Jun")) %>%
  mutate(
    seasonal_date = as.DATE(paste0("2000-", month(TIMESTAMP), "-", day(TIMESTAMP)))
  )

ggplot(data_subset,
       aes(x = seasonal_date, y = FCH4, color = flux_type)) +
  geom_point(aes(size = 2)) +
  geom_point(aes(y = DELTA_FCH4, size = 2),
             color = "#555555",
             shape = 2,
             stroke = 0.7,
             alpha = 0.4) +
  labs(x = " ",
       y = expression("FCH"[4] * " (nmol m"^-2 * " s"^-1 * ")")) +
  scale_color_manual(values = c("Chamber" = "#CD4071FF", "EC" = "#0096FF")) +
  scale_size_identity() +
  scale_x_date(date_labels = "%b", date_breaks = "1 month") +
  theme_bw() +
  theme(axis.title = element_text(size = 20), 
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text = element_text(size = 20),
        legend.title = element_blank(),
        legend.text = element_text(size = 20))

#### half-hourly plots with month facets

# Filter data for SITE = "US-UAF"
data_long_hh_summary_USUAF <- data_long_hh_summary_site %>%
  filter(SITE == "US-UAF")

data_long_hh_summary_USUAF <- data_long_hh_summary_USUAF %>% mutate(Month_letter =
                                                                          case_when(data_long_hh_summary_USUAF$MONTH == 06 ~ "Jun", 
                                                                                    data_long_hh_summary_USUAF$MONTH == 07 ~ "Jul",
                                                                                    data_long_hh_summary_USUAF$MONTH == 08 ~ "Aug",
                                                                                    data_long_hh_summary_USUAF$MONTH == 09 ~ "Sep",
                                                                                    data_long_hh_summary_USUAF$MONTH == 10 ~ "Oct",
                                                                                    data_long_hh_summary_USUAF$MONTH == 05 ~ "May")
)


month_order <- c("Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")

# Reorder the Month_letter column in the dataframe
data_long_hh_summary_USUAF$Month_letter <- factor(data_long_hh_summary_USUAF$Month_letter, 
                                                    levels = month_order, 
                                                    ordered = TRUE)



nud <- 0.18   # horizontal shift for plotting

data_long_hh_summary_USUAF$HOUR <- as.numeric(data_hh_hh_summary_USUAF$HOUR)

# plot (Fig. B12)

ggplot(data_long_hh_summary_USUAF, aes(x = HOUR, y = median_FCH4)) +
  
  ## EC nudged left
  geom_errorbar(data = subset(data_long_hh_summary_USUAF, flux_type == "EC"),
                aes(ymin = IQR_lower, ymax = IQR_upper, color = flux_type),
                width = 0.2, position = position_nudge(x = -nud)) +
  geom_line(data = subset(data_long_hh_summary_USUAF, flux_type == "EC"),
            aes(color = flux_type), position = position_nudge(x = -nud)) +
  geom_point(data = subset(data_long_hh_summary_USUAF, flux_type == "EC"),
             aes(color = flux_type), position = position_nudge(x = -nud), size = 2) +
  
  ## Chamber nudged right
  geom_errorbar(data = subset(data_long_hh_summary_USUAF, flux_type == "Chamber"),
                aes(ymin = IQR_lower, ymax = IQR_upper, color = flux_type),
                width = 0.2, position = position_nudge(x = +nud)) +
  geom_line(data = subset(data_long_hh_summary_USUAF, flux_type == "Chamber"),
            aes(color = flux_type), position = position_nudge(x = +nud)) +
  geom_point(data = subset(data_long_hh_summary_USUAF, flux_type == "Chamber"),
             aes(color = flux_type), position = position_nudge(x = +nud), size = 2) +
  
  labs(
    x = "HOUR (0-23)",
    y = expression("FCH"[4] * " (nmol m"^-2 * " s"^-1 * ")")
  ) +
  scale_color_manual(values = c("Chamber" = "#CD4071FF", "EC" = "#0096FF"),
                     breaks = c("Chamber","EC")) +
  scale_x_continuous(breaks = 0:23, expand = expansion(mult = 0, add = 0.25)) +
  facet_wrap(~ Month_letter, scales = "free") +
  theme_bw() +
  theme(
    axis.title = element_text(size = 20),
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text = element_text(size = 20),
    legend.title = element_blank(),
    legend.text = element_text(size = 20),
    strip.background = element_rect(fill = "white", color = NA), 
    strip.text = element_text(color = "black", size = 15)
  )



# SITE = US-Ho1
data_long_hh_summary_USHO1 <- data_long_hh_summary_site %>%
  filter(SITE == "US-HO1")


data_long_hh_summary_USHO1 <- data_long_hh_summary_USHO1 %>% mutate(Month_letter =
                                                                          case_when(data_long_hh_summary_USHO1$MONTH == 06 ~ "Jun", 
                                                                                    data_long_hh_summary_USHO1$MONTH == 07 ~ "Jul",
                                                                                    data_long_hh_summary_USHO1$MONTH == 08 ~ "Aug",
                                                                                    data_long_hh_summary_USHO1$MONTH == 09 ~ "Sep",
                                                                                    data_long_hh_summary_USHO1$MONTH == 10 ~ "Oct",
                                                                                    data_long_hh_summary_USHO1$MONTH == 11 ~ "Nov",
                                                                                    data_long_hh_summary_USHO1$MONTH == 03 ~ "Mar",
                                                                                    data_long_hh_summary_USHO1$MONTH == 04 ~ "Apr",
                                                                                    data_long_hh_summary_USHO1$MONTH == 05 ~ "May")
)


month_order <- c("Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov")

# Reorder the Month_letter column in the dataframe
data_long_hh_summary_USHO1$Month_letter <- factor(data_long_hh_summary_USHO1$Month_letter, 
                                                    levels = month_order, 
                                                    ordered = TRUE)




# Create the ggplot faceted by MONTH (Fig. B11)
nud <- 0.18   # horizontal shift for plotting

data_long_hh_summary_USHO1$HOUR <- as.numeric(data_long_hh_summary_USHO1$HOUR)

# plot
ggplot(data_long_hh_summary_USHO1, aes(x = HOUR, y = median_FCH4)) +
  
  ## EC nudged left
  geom_errorbar(data = subset(data_long_hh_summary_USHO1, flux_type == "EC"),
                aes(ymin = IQR_lower, ymax = IQR_upper, color = flux_type),
                width = 0.2, position = position_nudge(x = -nud)) +
  geom_line(data = subset(data_long_hh_summary_USHO1, flux_type == "EC"),
            aes(color = flux_type), position = position_nudge(x = -nud)) +
  geom_point(data = subset(data_long_hh_summary_USHO1, flux_type == "EC"),
             aes(color = flux_type), position = position_nudge(x = -nud), size = 2) +
  
  ## Chamber nudged right
  geom_errorbar(data = subset(data_long_hh_summary_USHO1, flux_type == "Chamber"),
                aes(ymin = IQR_lower, ymax = IQR_upper, color = flux_type),
                width = 0.2, position = position_nudge(x = +nud)) +
  geom_line(data = subset(data_long_hh_summary_USHO1, flux_type == "Chamber"),
            aes(color = flux_type), position = position_nudge(x = +nud)) +
  geom_point(data = subset(data_long_hh_summary_USHO1, flux_type == "Chamber"),
             aes(color = flux_type), position = position_nudge(x = +nud), size = 2) +
  
  labs(
    x = "HOUR (0-23)",
    y = expression("FCH"[4] * " (nmol m"^-2 * " s"^-1 * ")")
  ) +
  scale_color_manual(values = c("Chamber" = "#CD4071FF", "EC" = "#0096FF"),
                     breaks = c("Chamber","EC")) +
  scale_x_continuous(breaks = 0:23, expand = expansion(mult = 0, add = 0.25)) +
  facet_wrap(~ Month_letter, scales = "free") +
  theme_bw() +
  theme(
    axis.title = element_text(size = 20),
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text = element_text(size = 20),
    legend.title = element_blank(),
    legend.text = element_text(size = 20),
    strip.background = element_rect(fill = "white", color = NA), 
    strip.text = element_text(color = "black", size = 15)
  )



# SITE = SE-DEG
data_long_hh_summary_SEDEG <- data_long_hh_summary_site %>%
  filter(SITE == "SE-DEG")


data_long_hh_summary_SEDEG <- data_long_hh_summary_SEDEG %>% mutate(Month_letter =
                                                                          case_when(data_long_hh_summary_SEDEG$MONTH == 06 ~ "Jun", 
                                                                                    data_long_hh_summary_SEDEG$MONTH == 07 ~ "Jul",
                                                                                    data_long_hh_summary_SEDEG$MONTH == 08 ~ "Aug",
                                                                                    data_long_hh_summary_SEDEG$MONTH == 09 ~ "Sep",
                                                                                    data_long_hh_summary_SEDEG$MONTH == 10 ~ "Oct",
                                                                                    data_long_hh_summary_SEDEG$MONTH == 11 ~ "Nov",
                                                                                    data_long_hh_summary_SEDEG$MONTH == 05 ~ "May")
)


month_order <- c("May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov")

# Reorder the Month_letter column in the dataframe
data_long_hh_summary_SEDEG$Month_letter <- factor(data_long_hh_summary_SEDEG$Month_letter, 
                                                    levels = month_order, 
                                                    ordered = TRUE)


# Create the ggplot faceted by MONTH
nud <- 0.18   # horizontal shift for plotting

data_long_hh_summary_SEDEG$HOUR <- as.numeric(data_long_hh_summary_SEDEG$HOUR)

# plot (Fig. B10)
ggplot(data_long_hh_summary_SEDEG, aes(x = HOUR, y = median_FCH4)) +
  
  ## EC nudged left
  geom_errorbar(data = subset(data_long_hh_summary_SEDEG, flux_type == "EC"),
                aes(ymin = IQR_lower, ymax = IQR_upper, color = flux_type),
                width = 0.2, position = position_nudge(x = -nud)) +
  geom_line(data = subset(data_long_hh_summary_SEDEG, flux_type == "EC"),
            aes(color = flux_type), position = position_nudge(x = -nud)) +
  geom_point(data = subset(data_long_hh_summary_SEDEG, flux_type == "EC"),
             aes(color = flux_type), position = position_nudge(x = -nud), size = 2) +
  
  ## Chamber nudged right
  geom_errorbar(data = subset(data_long_hh_summary_SEDEG, flux_type == "Chamber"),
                aes(ymin = IQR_lower, ymax = IQR_upper, color = flux_type),
                width = 0.2, position = position_nudge(x = +nud)) +
  geom_line(data = subset(data_long_hh_summary_SEDEG, flux_type == "Chamber"),
            aes(color = flux_type), position = position_nudge(x = +nud)) +
  geom_point(data = subset(data_long_hh_summary_SEDEG, flux_type == "Chamber"),
             aes(color = flux_type), position = position_nudge(x = +nud), size = 2) +
  
  labs(
    x = "HOUR (0-23)",
    y = expression("FCH"[4] * " (nmol m"^-2 * " s"^-1 * ")")
  ) +
  scale_color_manual(values = c("Chamber" = "#CD4071FF", "EC" = "#0096FF"),
                     breaks = c("Chamber","EC")) +
  scale_x_continuous(breaks = 0:23, expand = expansion(mult = 0, add = 0.25)) +
  facet_wrap(~ Month_letter, scales = "free") +
  theme_bw() +
  theme(
    axis.title = element_text(size = 20),
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text = element_text(size = 20),
    legend.title = element_blank(),
    legend.text = element_text(size = 20),
    strip.background = element_rect(fill = "white", color = NA),  # remove facet box outline
    strip.text = element_text(color = "black", size = 15)
  )


# SITE = CN-HGU
data_long_hh_summary_CNHGU <- data_long_hh_summary_site %>%
  filter(SITE == "CN-HGU")

data_long_hh_summary_CNHGU <- data_long_hh_summary_CNHGU %>% mutate(Month_letter =
                                                                          case_when(data_long_hh_summary_CNHGU$MONTH == 06 ~ "Jun", 
                                                                                    data_long_hh_summary_CNHGU$MONTH == 07 ~ "Jul",
                                                                                    data_long_hh_summary_CNHGU$MONTH == 08 ~ "Aug",
                                                                                    data_long_hh_summary_CNHGU$MONTH == 09 ~ "Sep",
                                                                                    data_long_hh_summary_CNHGU$MONTH == 10 ~ "Oct",
                                                                                    data_long_hh_summary_CNHGU$MONTH == 11 ~ "Nov",
                                                                                    data_long_hh_summary_CNHGU$MONTH == 12 ~ "Dec",
                                                                                    data_long_hh_summary_CNHGU$MONTH == 02 ~ "Feb",
                                                                                    data_long_hh_summary_CNHGU$MONTH == 03 ~ "Mar",
                                                                                    data_long_hh_summary_CNHGU$MONTH == 04 ~ "Apr",
                                                                                    data_long_hh_summary_CNHGU$MONTH == 05 ~ "May")
)


month_order <- c("Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")

# Reorder the Month_letter column in the dataframe
data_long_hh_summary_CNHGU$Month_letter <- factor(data_long_hh_summary_CNHGU$Month_letter, 
                                                    levels = month_order, 
                                                    ordered = TRUE)


# Create the ggplot faceted by MONTH
nud <- 0.18   # horizontal shift for plotting

data_long_hh_summary_CNHGU$HOUR <- as.numeric(data_long_hh_summary_CNHGU$HOUR)

# plot (Fig. B9)
ggplot(data_long_hh_summary_CNHGU, aes(x = HOUR, y = median_FCH4)) +
  
  ## EC nudged left
  geom_errorbar(data = subset(data_long_hh_summary_CNHGU, flux_type == "EC"),
                aes(ymin = IQR_lower, ymax = IQR_upper, color = flux_type),
                width = 0.2, position = position_nudge(x = -nud)) +
  geom_line(data = subset(data_long_hh_summary_CNHGU, flux_type == "EC"),
            aes(color = flux_type), position = position_nudge(x = -nud)) +
  geom_point(data = subset(data_long_hh_summary_CNHGU, flux_type == "EC"),
             aes(color = flux_type), position = position_nudge(x = -nud), size = 2) +
  
  ## Chamber nudged right
  geom_errorbar(data = subset(data_long_hh_summary_CNHGU, flux_type == "Chamber"),
                aes(ymin = IQR_lower, ymax = IQR_upper, color = flux_type),
                width = 0.2, position = position_nudge(x = +nud)) +
  geom_line(data = subset(data_long_hh_summary_CNHGU, flux_type == "Chamber"),
            aes(color = flux_type), position = position_nudge(x = +nud)) +
  geom_point(data = subset(data_long_hh_summary_CNHGU, flux_type == "Chamber"),
             aes(color = flux_type), position = position_nudge(x = +nud), size = 2) +
  
  labs(
    x = "HOUR (0-23)",
    y = expression("FCH"[4] * " (nmol m"^-2 * " s"^-1 * ")")
  ) +
  scale_color_manual(values = c("Chamber" = "#CD4071FF", "EC" = "#0096FF"),
                     breaks = c("Chamber","EC")) +
  scale_x_continuous(breaks = 0:23, expand = expansion(mult = 0, add = 0.25)) +
  facet_wrap(~ Month_letter, scales = "free") +
  theme_bw() +
  theme(
    axis.title = element_text(size = 20),
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text = element_text(size = 20),
    legend.title = element_blank(),
    legend.text = element_text(size = 20),
    strip.background = element_rect(fill = "white", color = NA),  # remove facet box outline
    strip.text = element_text(color = "black", size = 15)
  )


# all sites together (n=4 sites)

data_long_hh_summary <- data_long_hh_summary %>% mutate(Month_letter =
                                                              case_when(data_long_hh_summary$MONTH == 06 ~ "Jun", 
                                                                        data_long_hh_summary$MONTH == 07 ~ "Jul",
                                                                        data_long_hh_summary$MONTH == 08 ~ "Aug",
                                                                        data_long_hh_summary$MONTH == 09 ~ "Sep",
                                                                        data_long_hh_summary$MONTH == 10 ~ "Oct",
                                                                        data_long_hh_summary$MONTH == 11 ~ "Nov",
                                                                        data_long_hh_summary$MONTH == 12 ~ "Dec",
                                                                        data_long_hh_summary$MONTH == 02 ~ "Feb",
                                                                        data_long_hh_summary$MONTH == 03 ~ "Mar",
                                                                        data_long_hh_summary$MONTH == 04 ~ "Apr",
                                                                        data_long_hh_summary$MONTH == 05 ~ "May")
)


month_order <- c("Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")

# Reorder the Month_letter column in the dataframe
data_long_hh_summary$Month_letter <- factor(data_long_hh_summary$Month_letter, 
                                              levels = month_order, 
                                              ordered = TRUE)

data_long_hh_summary$HOUR <- as.numeric(data_long_hh_summary$HOUR)

# exclude months where not all sites have observations

# Filter the data to include only the months from May to October
filtered_data <- data_long_hh_summary %>%
  filter(Month_letter %in% c("May", "Jun", "Jul", "Aug", "Sep", "Oct"))

ec <- filtered_data %>% filter(flux_type == "EC")
ch <- filtered_data %>% filter(flux_type == "Chamber")

pos <- position_dodge(width = 0.5)

nud <- 0.18  # adjust spacing for plotting

# plot all sites together (Fig. B8)
ggplot(filtered_data, aes(x = HOUR, y = median_FCH4)) +
  # EC nudged left
  geom_errorbar(data = subset(filtered_data, flux_type == "EC"),
                aes(ymin = IQR_lower, ymax = IQR_upper, color = flux_type),
                width = 0.2, position = position_nudge(x = -nud)) +
  geom_line(data = subset(filtered_data, flux_type == "EC"),
            aes(color = flux_type), position = position_nudge(x = -nud)) +
  geom_point(data = subset(filtered_data, flux_type == "EC"),
             aes(color = flux_type), position = position_nudge(x = -nud), size = 2) +
  
  # Chamber nudged right
  geom_errorbar(data = subset(filtered_data, flux_type == "Chamber"),
                aes(ymin = IQR_lower, ymax = IQR_upper, color = flux_type),
                width = 0.2, position = position_nudge(x = +nud)) +
  geom_line(data = subset(filtered_data, flux_type == "Chamber"),
            aes(color = flux_type), position = position_nudge(x = +nud)) +
  geom_point(data = subset(filtered_data, flux_type == "Chamber"),
             aes(color = flux_type), position = position_nudge(x = +nud), size = 2) +
  
  scale_color_manual(values = c("Chamber" = "#CD4071FF", "EC" = "#0096FF"),
                     breaks = c("Chamber", "EC")) +
  scale_x_continuous(breaks = 0:23) +
  facet_wrap(~ Month_letter, scales = "free") +
  theme_bw() +
  theme(
    strip.background = element_rect(fill = "white", color = NA),  # remove facet box outline
    axis.title = element_text(size = 20),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 15),
    axis.text = element_text(size = 20),
    axis.title.x = element_blank(),
    legend.title = element_blank(),
    legend.text = element_text(size = 20),
    strip.text = element_text(color = "black", size = 15)
  )

############################################################################################################################################################
### KRUSKAL-WALLIS tests for how delta_FCH4 differs between HOUR, MONTH AND YEAR ###


all_d_aggr <- all_d_aggr %>% mutate(Month_letter =
                                                            case_when(month(all_d_aggr$TIMESTAMP) == 06 ~ "Jun", 
                                                                      month(all_d_aggr$TIMESTAMP) == 07 ~ "Jul",
                                                                      month(all_d_aggr$TIMESTAMP) == 08 ~ "Aug",
                                                                      month(all_d_aggr$TIMESTAMP) == 09 ~ "Sep",
                                                                      month(all_d_aggr$TIMESTAMP) == 10 ~ "Oct",
                                                                      month(all_d_aggr$TIMESTAMP) == 11 ~ "Nov",
                                                                      month(all_d_aggr$TIMESTAMP) == 12 ~ "Dec",
                                                                      month(all_d_aggr$TIMESTAMP) == 01 ~ "Jan",
                                                                      month(all_d_aggr$TIMESTAMP) == 02 ~ "Feb",
                                                                      month(all_d_aggr$TIMESTAMP) == 03 ~ "Mar",
                                                                      month(all_d_aggr$TIMESTAMP) == 04 ~ "Apr",
                                                                      month(all_d_aggr$TIMESTAMP) == 05 ~ "May")
)

all_week_aggr <- all_week_aggr %>% mutate(Month_letter =
                                            case_when(all_week_aggr$MONTH == "6" ~ "Jun", 
                                                      all_week_aggr$MONTH == "7" ~ "Jul",
                                                      all_week_aggr$MONTH == "8" ~ "Aug",
                                                      all_week_aggr$MONTH == "9" ~ "Sep",
                                                      all_week_aggr$MONTH == "10" ~ "Oct",
                                                      all_week_aggr$MONTH == "11" ~ "Nov",
                                                      all_week_aggr$MONTH == "12" ~ "Dec",
                                                      all_week_aggr$MONTH == "1" ~ "Jan",
                                                      all_week_aggr$MONTH == "2" ~ "Feb",
                                                      all_week_aggr$MONTH == "3" ~ "Mar",
                                                      all_week_aggr$MONTH == "4" ~ "Apr",
                                                      all_week_aggr$MONTH == "5" ~ "May")
)

all_month_aggr <- all_month_aggr %>% mutate(Month_letter =
                                              case_when(all_month_aggr$MONTH == "6" ~ "Jun", 
                                                        all_month_aggr$MONTH == "7" ~ "Jul",
                                                        all_month_aggr$MONTH == "8" ~ "Aug",
                                                        all_month_aggr$MONTH == "9" ~ "Sep",
                                                        all_month_aggr$MONTH == "10" ~ "Oct",
                                                        all_month_aggr$MONTH == "11" ~ "Nov",
                                                        all_month_aggr$MONTH == "12" ~ "Dec",
                                                        all_month_aggr$MONTH == "1" ~ "Jan",
                                                        all_month_aggr$MONTH == "2" ~ "Feb",
                                                        all_month_aggr$MONTH == "3" ~ "Mar",
                                                        all_month_aggr$MONTH == "4" ~ "Apr",
                                                        all_month_aggr$MONTH == "5" ~ "May")
)


all_hr_aggr <- all_hr_aggr %>% mutate(Month_letter =
                                        case_when(all_hr_aggr$MONTH == "6" ~ "Jun", 
                                                  all_hr_aggr$MONTH == "7" ~ "Jul",
                                                  all_hr_aggr$MONTH == "8" ~ "Aug",
                                                  all_hr_aggr$MONTH == "9" ~ "Sep",
                                                  all_hr_aggr$MONTH == "10" ~ "Oct",
                                                  all_hr_aggr$MONTH == "11" ~ "Nov",
                                                  all_hr_aggr$MONTH == "12" ~ "Dec",
                                                  all_hr_aggr$MONTH == "1" ~ "Jan",
                                                  all_hr_aggr$MONTH == "2" ~ "Feb",
                                                  all_hr_aggr$MONTH == "3" ~ "Mar",
                                                  all_hr_aggr$MONTH == "4" ~ "Apr",
                                                  all_hr_aggr$MONTH == "5" ~ "May")
)


all_hh_aggr <- all_hh_aggr %>% mutate(Month_letter =
                                                                case_when(all_hh_aggr$MONTH == "6" ~ "Jun", 
                                                                          all_hh_aggr$MONTH == "7" ~ "Jul",
                                                                          all_hh_aggr$MONTH == "8" ~ "Aug",
                                                                          all_hh_aggr$MONTH == "9" ~ "Sep",
                                                                          all_hh_aggr$MONTH == "10" ~ "Oct",
                                                                          all_hh_aggr$MONTH == "11" ~ "Nov",
                                                                          all_hh_aggr$MONTH == "12" ~ "Dec",
                                                                          all_hh_aggr$MONTH == "1" ~ "Jan",
                                                                          all_hh_aggr$MONTH == "2" ~ "Feb",
                                                                          all_hh_aggr$MONTH == "3" ~ "Mar",
                                                                          all_hh_aggr$MONTH == "4" ~ "Apr",
                                                                          all_hh_aggr$MONTH == "5" ~ "May")
)


# Hourly
all_hr_aggr$TIMESTAMP <- as_datetime(all_hr_aggr$TIMESTAMP)

# Daily

all_d_aggr$TIMESTAMP <- as_datetime(all_d_aggr$TIMESTAMP)
all_d_aggr$YEAR <- year(all_d_aggr$TIMESTAMP)
all_d_aggr$MONTH <- month(all_d_aggr$TIMESTAMP)
all_d_aggr$DOY <- yday(all_d_aggr$TIMESTAMP)


# Kruskal-Wallis and Conover tests

# Half-hourly

kruskal.test(DELTA_FCH4 ~ as.factor(HOUR), data = all_hh_aggr)

conover.test(all_hh_aggr$DELTA_FCH4, all_hh_aggr$HOUR, method = "bonferroni")

kruskal.test(DELTA_FCH4 ~ as.factor(MONTH), data = all_hh_aggr)

conover.test(all_hh_aggr$DELTA_FCH4, all_hh_aggr$MONTH, method = "bonferroni")

kruskal.test(DELTA_FCH4 ~ YEAR, data = all_hh_aggr)

kruskal.test(DELTA_FCH4 ~ SITE, data = all_hh_aggr)

# Hourly

kruskal.test(DELTA_FCH4 ~ HOUR, data = all_hr_aggr)

kruskal.test(DELTA_FCH4 ~ YEAR, data = all_hr_aggr)

kruskal.test(DELTA_FCH4 ~ as.factor(MONTH), data = all_hr_aggr)

conover.test(all_hr_aggr$DELTA_FCH4, all_hr_aggr$MONTH, method = "bonferroni")

kruskal.test(DELTA_FCH4 ~ SITE, data = all_hr_aggr)

# make a heatmap for comparing delta_FCH4 between hours

# half-hourly (Fig. B6)

result <- conover.test(all_hh_aggr$DELTA_FCH4, all_hh_aggr$HOUR, method = "bonferroni")
result <- as.data.frame(result)

# Extract p-values from the result
p_values <- result$P.adjusted

# Get the pairwise comparisons
comparisons <- result$comparisons

# helper to pretty-print p-values
pretty_p <- function(x, dec_digits = 4, sci_digits = 0, sci_cut = 1e-3) {
  vapply(x, function(p) {
    if (is.na(p)) return(NA_character_)
    if (p == 0)  return("0")
    if (p == 1)  return("1")
    if (p < sci_cut) return(sprintf(paste0("%.", sci_digits, "e"), p))  # e.g., 6e-04
    # fixed format, then trim trailing zeros and trailing dot
    fx <- sprintf(paste0("%.", dec_digits, "f"), p)
    sub("\\.?0+$", "", fx)
  }, character(1))
}


# Tidy the conover.test output
pairs_df <- tibble(
  comp = result$comparisons,
  p    = result$P.adjusted
) %>%
  separate(comp, into = c("a", "b"), sep = " - ", convert = TRUE) %>%
  mutate(x = pmax(a, b), y = pmin(a, b)) %>%
  distinct(x, y, .keep_all = TRUE) %>%
  mutate(
    Var1 = factor(x, levels = 0:23),
    Var2 = factor(y, levels = 0:23),
    color_category = case_when(
      p >= 0.05  ~ "p >= 0.05",
      p >= 0.01  ~ "0.01 <= p < 0.05",
      p >= 0.001 ~ "0.001 <= p < 0.01",
      TRUE       ~ "p < 0.001"
    ),
    # mixed formatting
    label_num = pretty_p(p, dec_digits = 4, sci_digits = 0, sci_cut = 1e-3),
    label     = if_else(p < 0.05, paste0(label_num, "*"), label_num),
    text_col  = if_else(p < 0.001, "white", "black"),
    fontface  = if_else(p < 0.05, "bold", "plain")
  )

# Diagonal labels (no tiles on the diagonal)
diag_df <- tibble(
  Var1 = factor(0:23, levels = 0:23),
  Var2 = factor(0:23, levels = 0:23),
  lab  = 0:23
)

# Plot only the chosen triangle (corner at bottom-right)
ggplot(pairs_df, aes(Var1, Var2)) +
  geom_tile(aes(fill = color_category), color = "grey85", linewidth = 0.3) +
  geom_text(aes(label = label, color = text_col, fontface = fontface), size = 3) +
  geom_text(data = diag_df, aes(Var1, Var2, label = lab), inherit.aes = FALSE, size = 6) +
  scale_x_discrete(name = "HOUR", expand = expansion(0)) +
  scale_y_discrete(name = "HOUR", expand = expansion(0)) +
  scale_fill_manual(
    name   = expression(italic(p) ~ "-value"),
    values = c(
      "p >= 0.05" = "white",
      "0.01 <= p < 0.05" = "#FEC98DFF",
      "0.001 <= p < 0.01" = "#F1605DFF",
      "p < 0.001" = "#9F2F7FFF"
    ),
    breaks = c("p >= 0.05", "0.01 <= p < 0.05", "0.001 <= p < 0.01", "p < 0.001"),
    labels = c(
      expression(italic(p) >= 0.05),
      expression("0.01" <= italic(p) ~ "< 0.05"),
      expression("0.001" <= italic(p) ~ "< 0.01"),
      expression(italic(p) < 0.001)
    ),
    na.value = "white"
  ) +
  scale_color_identity() +
  theme_minimal() +
  theme(
    axis.title = element_blank(),
    axis.text  = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank(),
    legend.position = c(0.02, 0.98),
    legend.justification = c(0, 1),
    legend.text = element_text(size = 15),
    legend.title = element_text(size = 15)
  )


##### hourly aggregation (Fig. B7)

result <- conover.test(all_hr_aggr$DELTA_FCH4, all_hr_aggr$HOUR, method = "bonferroni")
result <- as.data.frame(result)

# Extract p-values from the result
p_values <- result$P.adjusted

# Get the pairwise comparisons
comparisons <- result$comparisons

# Tidy the conover.test output
pairs_df <- tibble(
  comp = result$comparisons,
  p    = result$P.adjusted
) %>%
  separate(comp, into = c("a", "b"), sep = " - ", convert = TRUE) %>%
  mutate(x = pmax(a, b), y = pmin(a, b)) %>%
  distinct(x, y, .keep_all = TRUE) %>%
  mutate(
    Var1 = factor(x, levels = 0:23),
    Var2 = factor(y, levels = 0:23),
    color_category = case_when(
      p >= 0.05  ~ "p >= 0.05",
      p >= 0.01  ~ "0.01 <= p < 0.05",
      p >= 0.001 ~ "0.001 <= p < 0.01",
      TRUE       ~ "p < 0.001"
    ),
    # mixed formatting
    label_num = pretty_p(p, dec_digits = 4, sci_digits = 0, sci_cut = 1e-3),
    label     = if_else(p < 0.05, paste0(label_num, "*"), label_num),
    text_col  = if_else(p < 0.001, "white", "black"),
    fontface  = if_else(p < 0.05, "bold", "plain")
  )

# Diagonal labels (no tiles on the diagonal)
diag_df <- tibble(
  Var1 = factor(0:23, levels = 0:23),
  Var2 = factor(0:23, levels = 0:23),
  lab  = 0:23
)

# Plot only the chosen triangle (corner at bottom-right)
ggplot(pairs_df, aes(Var1, Var2)) +
  geom_tile(aes(fill = color_category), color = "grey85", linewidth = 0.3) +
  geom_text(aes(label = label, color = text_col, fontface = fontface), size = 3) +
  geom_text(data = diag_df, aes(Var1, Var2, label = lab), inherit.aes = FALSE, size = 6) +
  scale_x_discrete(name = "HOUR", expand = expansion(0)) +
  scale_y_discrete(name = "HOUR", expand = expansion(0)) +
  scale_fill_manual(
    name   = expression(italic(p) ~ "-value"),
    values = c(
      "p >= 0.05" = "white",
      "0.01 <= p < 0.05" = "#FEC98DFF",
      "0.001 <= p < 0.01" = "#F1605DFF",
      "p < 0.001" = "#9F2F7FFF"
    ),
    breaks = c("p >= 0.05", "0.01 <= p < 0.05", "0.001 <= p < 0.01", "p < 0.001"),
    labels = c(
      expression(italic(p) >= 0.05),
      expression("0.01" <= italic(p) ~ "< 0.05"),
      expression("0.001" <= italic(p) ~ "< 0.01"),
      expression(italic(p) < 0.001)
    ),
    na.value = "white"
  ) +
  scale_color_identity() +
  theme_minimal() +
  theme(
    axis.title = element_blank(),
    axis.text  = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank(),
    legend.position = c(0.02, 0.98),
    legend.justification = c(0, 1),
    legend.text = element_text(size = 15),
    legend.title = element_text(size = 15)
  )


## Kruskal-Wallis and Conover tests for daily-annual aggregations

kruskal.test(DELTA_FCH4 ~ as.factor(MONTH), data = all_d_aggr)

kruskal.test(DELTA_FCH4 ~ YEAR, data = all_d_aggr)

kruskal.test(DELTA_FCH4 ~ SITE, data = all_d_aggr)

kruskal.test(DELTA_FCH4 ~ as.factor(DOMINANT_VEGETATION), data = all_d_aggr)

conover.test(all_d_aggr$DELTA_FCH4, all_d_aggr$DOMINANT_VEGETATION)

conover.test(abs(all_d_aggr$DELTA_FCH4), all_d_aggr$MONTH, method = "bonferroni")

# weekly

kruskal.test(DELTA_FCH4 ~ as.factor(DOMINANT_VEGETATION), data = all_week_aggr)
conover.test(all_week_aggr$DELTA_FCH4, all_week_aggr$DOMINANT_VEGETATION)


kruskal.test(abs(DELTA_FCH4) ~ as.factor(MONTH), data = all_week_aggr)

kruskal.test(DELTA_FCH4 ~ YEAR, data = all_week_aggr)

kruskal.test(DELTA_FCH4 ~ SITE, data = all_week_aggr)

# monthly

kruskal.test(DELTA_FCH4 ~ as.factor(DOMINANT_VEGETATION), data = all_month_aggr)
conover.test(all_month_aggr$DELTA_FCH4, all_month_aggr$DOMINANT_VEGETATION)


kruskal.test(abs(DELTA_FCH4) ~ as.factor(MONTH), data = all_month_aggr)

kruskal.test(DELTA_FCH4 ~ YEAR, data = all_month_aggr)

kruskal.test(DELTA_FCH4 ~ SITE, data = all_month_aggr)

# annual

kruskal.test(DELTA_FCH4 ~ YEAR, data = all_yr_aggr)

kruskal.test(DELTA_FCH4 ~ SITE, data = all_yr_aggr)

kruskal.test(DELTA_FCH4 ~ as.factor(DOMINANT_VEGETATION), data = all_yr_aggr)


### plot monthly delta_FCH4 (Fig. B13) based on conover tests significant differences

# modify dataframes 

# half-hourly

# Identify and remove outliers
filtered_data_hh <- all_hh_aggr %>%
  group_by(Month_letter) %>%
  mutate(
    Q1 = quantile(DELTA_FCH4, 0.25, na.rm = TRUE),
    Q3 = quantile(DELTA_FCH4, 0.75, na.rm = TRUE),
    IQR = Q3 - Q1,
    lower_bound = Q1 - 1.5 * IQR,
    upper_bound = Q3 + 1.5 * IQR
  ) %>%
  filter(DELTA_FCH4 >= lower_bound & DELTA_FCH4 <= upper_bound) %>%
  dplyr::select(-Q1, -Q3, -IQR, -lower_bound, -upper_bound)  # Remove temporary columns

filtered_data_hh$Month_letter <- factor(filtered_data_hh$MONTH, levels = as.character(2:11))

filtered_data_hh$Month_letter <- factor(
  filtered_data_hh$Month_letter,
  levels = c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")
)

filtered_halfhourly <- filtered_data_hh %>%
  mutate(Aggregation = "Half-hourly")

# absolute delta_FCH4
filtered_data_hh_abs <- all_hh_aggr %>%
  group_by(Month_letter) %>%
  mutate(
    Q1 = quantile(DELTA_FCH4_abs, 0.25, na.rm = TRUE),
    Q3 = quantile(DELTA_FCH4_abs, 0.75, na.rm = TRUE),
    IQR = Q3 - Q1,
    lower_bound = Q1 - 1.5 * IQR,
    upper_bound = Q3 + 1.5 * IQR
  ) %>%
  filter(DELTA_FCH4_abs >= lower_bound & DELTA_FCH4_abs <= upper_bound) %>%
  dplyr::select(-Q1, -Q3, -IQR, -lower_bound, -upper_bound)  # Remove temporary columns


filtered_data_hh_abs$Month_letter <- factor(
  filtered_data_hh_abs$Month_letter,
  levels = c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")
)

filtered_halfhourly_abs <- filtered_data_hh_abs %>%
  mutate(Aggregation = "Half-hourly")



# Hourly
# Identify and remove outliers
filtered_data_hr <- all_hr_aggr %>%
  group_by(Month_letter) %>%
  mutate(
    Q1 = quantile(DELTA_FCH4, 0.25, na.rm = TRUE),
    Q3 = quantile(DELTA_FCH4, 0.75, na.rm = TRUE),
    IQR = Q3 - Q1,
    lower_bound = Q1 - 1.5 * IQR,
    upper_bound = Q3 + 1.5 * IQR
  ) %>%
  filter(DELTA_FCH4 >= lower_bound & DELTA_FCH4 <= upper_bound) %>%
  dplyr::select(-Q1, -Q3, -IQR, -lower_bound, -upper_bound)  # Remove temporary columns

filtered_data_hr$Month_letter <- factor(filtered_data_hr$Month_letter, levels = as.character(2:11))

filtered_data_hr$Month_letter <- factor(
  filtered_data_hr$Month_letter,
  levels = c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")
)

filtered_hourly <- filtered_data_hr %>%
  mutate(Aggregation = "Hourly")

# absolute delta_FCH4

all_hr_aggr$DELTA_FCH4_abs <- abs(all_hr_aggr$DELTA_FCH4)

filtered_data_hr_abs <- all_hr_aggr %>%
  group_by(Month_letter) %>%
  mutate(
    Q1 = quantile(DELTA_FCH4_abs, 0.25, na.rm = TRUE),
    Q3 = quantile(DELTA_FCH4_abs, 0.75, na.rm = TRUE),
    IQR = Q3 - Q1,
    lower_bound = Q1 - 1.5 * IQR,
    upper_bound = Q3 + 1.5 * IQR
  ) %>%
  filter(DELTA_FCH4_abs >= lower_bound & DELTA_FCH4_abs <= upper_bound) %>%
  dplyr::select(-Q1, -Q3, -IQR, -lower_bound, -upper_bound)  # Remove temporary columns


filtered_data_hr_abs$Month_letter <- factor(
  filtered_data_hr_abs$Month_letter,
  levels = c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")
)

filtered_hourly_abs <- filtered_data_hr_abs %>%
  mutate(Aggregation = "Hourly")




# daily
# Identify and remove outliers
# absolute delta_FCH4

all_d_aggr$DELTA_FCH4_abs <- abs(all_d_aggr$DELTA_FCH4)

filtered_data_d_abs <- all_d_aggr %>%
  group_by(Month_letter) %>%
  mutate(
    Q1 = quantile(DELTA_FCH4_abs, 0.25, na.rm = TRUE),
    Q3 = quantile(DELTA_FCH4_abs, 0.75, na.rm = TRUE),
    IQR = Q3 - Q1,
    lower_bound = Q1 - 1.5 * IQR,
    upper_bound = Q3 + 1.5 * IQR
  ) %>%
  filter(DELTA_FCH4_abs >= lower_bound & DELTA_FCH4_abs <= upper_bound) %>%
  dplyr::select(-Q1, -Q3, -IQR, -lower_bound, -upper_bound)  # Remove temporary columns


filtered_data_d_abs$Month_letter <- factor(
  filtered_data_d_abs$Month_letter,
  levels = c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")
)

filtered_daily_abs <- filtered_data_d_abs %>%
  mutate(Aggregation = "Daily")

# subset the dataframes

filtered_halfhourly_2_abs <- subset(filtered_halfhourly_abs, select = c(Aggregation, Month_letter, DELTA_FCH4_abs))
filtered_hourly_2_abs <- subset(filtered_hourly_abs, select = c(Aggregation, Month_letter, DELTA_FCH4_abs))
filtered_daily_2_abs <- subset(filtered_daily_abs, select = c(Aggregation, Month_letter, DELTA_FCH4_abs))

# combine dataframes into one
combined_filtered_abs <- bind_rows(
  filtered_halfhourly_2_abs,
  filtered_hourly_2_abs,
  filtered_daily_2_abs
)

# add aggregation info
combined_filtered_abs$Aggregation <- factor(
  combined_filtered_abs$Aggregation,
  levels = c("Half-hourly", "Hourly", "Daily")
)

# check which months are included in each site

all_hh_aggr %>%
  distinct(SITE, Month_letter) %>%
  mutate(present = TRUE) %>%
  tidyr::pivot_wider(names_from = Month_letter, values_from = present, values_fill = FALSE)
# May-Oct

all_hr_aggr %>%
  distinct(SITE, Month_letter) %>%
  mutate(present = TRUE) %>%
  tidyr::pivot_wider(names_from = Month_letter, values_from = present, values_fill = FALSE)
# May-Oct

all_d_aggr %>%
  distinct(SITE, Month_letter) %>%
  mutate(present = TRUE) %>%
  tidyr::pivot_wider(names_from = Month_letter, values_from = present, values_fill = FALSE)

# do May-Oct but mention that in May all except 2 and Jun all except 1

# include only May-Oct

plot_data <- combined_filtered_abs %>%
  filter(Month_letter %in% c("May", "Jun", "Jul", "Aug", "Sep", "Oct"))

# get comparisons

month_order <- c("May", "Jun", "Jul", "Aug", "Sep", "Oct")

plot_data_filtered <- all_d_aggr %>%
  filter(Month_letter %in% month_order)

plot_data_filtered <- all_hr_aggr %>%
  filter(Month_letter %in% month_order)

plot_data_filtered <- all_hh_aggr %>%
  filter(Month_letter %in% month_order)

# Run test
conover_res <- conover.test(plot_data_filtered$DELTA_FCH4_abs, plot_data_filtered$Month_letter)
conover.test(plot_data_filtered$DELTA_FCH4_abs, plot_data_filtered$Month_letter)

# Clean and match comparisons to p-values
comparisons_raw <- conover_res$comparisons
comparisons_clean <- sapply(strsplit(comparisons_raw, " - "), function(x) {
  paste(trimws(x), collapse = " - ")
})

# Parse comparisons and filter to only desired months
group_pairs <- strsplit(comparisons_clean, " - ")
valid <- sapply(group_pairs, function(x) all(x %in% month_order))

# Use only valid p-values and comparisons
pvals <- conover_res$P[valid]
comparisons <- comparisons_clean[valid]

# Generate CLD
cld <- cldList(comparison = comparisons,
               p.value = pvals,
               threshold = 0.05)

# Rename column
names(cld)[1] <- "Month_letter"

label_y <- plot_data %>%
  group_by(Month_letter) %>%
  summarise(y = max(abs(DELTA_FCH4_abs), na.rm = TRUE) * 1.05)

cld_annot <- left_join(cld, label_y, by = "Month_letter")

conover_res

# plot (Fig. B13, cleaned afterward with Inkscape)
ggplot(plot_data, aes(x = Month_letter, y = DELTA_FCH4_abs)) +
  geom_jitter(
    aes(fill = Aggregation),
    position = position_jitterdodge(jitter.width = 0.3, dodge.width = 1),
    size = 1,
    color = "lightgrey",  
    alpha = 0.3,
    show.legend = FALSE
  ) +
  geom_boxplot(
    aes(fill = Aggregation),
    outlier.shape = NA,
    position = position_dodge(width = 1),
    color = "#636363"
  ) +
  geom_text(
    data = cld_annot,
    aes(x = Month_letter, y = y, label = Letter),
    position = position_dodge(width = 1),
    fontface = "bold",
    size = 6,
    vjust = 0
  ) +
  geom_hline(yintercept = 0, color = "blue", linetype = "dashed", size = 1) +
  labs(
    y = expression("Absolute" ~ Delta * FCH[4] ~ (nmol ~ m^{-2} ~ s^{-1})),
    fill = "Temporal aggregation",
    color = "Temporal aggregation"
  ) +
  theme_bw() +
  scale_fill_viridis_d(option = "magma", begin = 0.5, end = 0.85) +
  scale_color_viridis_d(option = "magma", begin = 0.5, end = 0.85) +
  theme(
    plot.title = element_text(size = 20, hjust = 0.5),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 20),
    axis.text = element_text(size = 20),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.text = element_text(size = 20),
    legend.title = element_blank()
  )


# check monthly delta_FCH4 descriptive statistics

# monthly 
month_summary <- all_month_aggr %>%
  group_by(MONTH) %>%
  summarise(
    Median = median(abs(DELTA_FCH4), na.rm = TRUE),
    Mean = mean(abs(DELTA_FCH4), na.rm = TRUE),
    IQR = IQR(abs(DELTA_FCH4), na.rm = TRUE),
    SD = sd(abs(DELTA_FCH4), na.rm = TRUE),
    CV = raster::cv(abs(DELTA_FCH4), na.rm=TRUE)
  )

# weekly
week_summary <- all_week_aggr %>%
  group_by(MONTH) %>%
  summarise(
    Median = median(DELTA_FCH4_abs, na.rm = TRUE),
    Mean = mean(DELTA_FCH4_abs, na.rm = TRUE),
    IQR = IQR(DELTA_FCH4_abs, na.rm = TRUE),
    SD = sd(DELTA_FCH4_abs, na.rm = TRUE),
    CV = raster::cv(DELTA_FCH4_abs, na.rm=TRUE)
  )

# daily
d_summary <- all_d_aggr %>%
  group_by(MONTH) %>%
  summarise(
    Median = median(DELTA_FCH4_abs, na.rm = TRUE),
    Mean = mean(DELTA_FCH4_abs, na.rm = TRUE),
    IQR = IQR(DELTA_FCH4_abs, na.rm = TRUE),
    SD = sd(DELTA_FCH4_abs, na.rm = TRUE),
    CV = raster::cv(DELTA_FCH4_abs, na.rm=TRUE)
  )

# hourly
hr_summary <- all_hr_aggr %>%
  group_by(MONTH) %>%
  summarise(
    Median = median(DELTA_FCH4_abs, na.rm = TRUE),
    Mean = mean(DELTA_FCH4_abs, na.rm = TRUE),
    IQR = IQR(DELTA_FCH4_abs, na.rm = TRUE),
    SD = sd(DELTA_FCH4_abs, na.rm = TRUE),
    CV = raster::cv(DELTA_FCH4_abs, na.rm=TRUE)
  )

# half-hourly
hh_summary <- all_hh_aggr %>%
  group_by(MONTH) %>%
  summarise(
    Median = median(DELTA_FCH4_abs, na.rm = TRUE),
    Mean = mean(DELTA_FCH4_abs, na.rm = TRUE),
    IQR = IQR(DELTA_FCH4_abs, na.rm = TRUE),
    SD = sd(DELTA_FCH4_abs, na.rm = TRUE),
    CV = raster::cv(DELTA_FCH4_abs, na.rm=TRUE)
  )


### diel trends

# Calculate medians and IQRs for each hour
medians_IQRs <- all_hh_aggr %>%
  group_by(HOUR) %>% 
  summarize(
    median_DELTA_FCH4 = median(DELTA_FCH4, na.rm = TRUE),
    IQR = IQR(DELTA_FCH4, na.rm = TRUE) 
  )

# Identify the hours with the highest and lowest medians
highest_median_info <- medians_IQRs %>%
  arrange(desc(median_DELTA_FCH4)) %>%
  slice(1)

lowest_median_info <- medians_IQRs %>%
  arrange(median_DELTA_FCH4) %>%
  slice(1)

highest_median_info$HOUR
highest_median_info$median_DELTA_FCH4
highest_median_info$IQR

lowest_median_info$HOUR
lowest_median_info$median_DELTA_FCH4
lowest_median_info$IQR


## differences between chamber methods (automated vs manual)

# daily
all_d_aggr$CH_METHOD <- as.factor(all_d_aggr$CH_METHOD)

wilcox.test(DELTA_FCH4 ~ CH_METHOD, data = all_d_aggr)

all_d_aggr %>%
  group_by(CH_METHOD) %>%
  summarize(non_na_count = sum(!is.na(DELTA_FCH4)))

# weekly
all_week_aggr$CH_METHOD <- as.factor(all_week_aggr$CH_METHOD)
wilcox.test(DELTA_FCH4 ~ CH_METHOD, data = all_week_aggr)

all_week_aggr %>%
  group_by(CH_METHOD) %>%
  summarize(non_na_count = sum(!is.na(DELTA_FCH4)))

# monthly
wilcox.test(DELTA_FCH4 ~ CH_METHOD, data = all_month_aggr)

all_month_aggr %>%
  group_by(CH_METHOD) %>%
  summarize(non_na_count = sum(!is.na(DELTA_FCH4)))

# annual
wilcox.test(DELTA_FCH4 ~ as.factor(CH_METHOD), data = all_yr_aggr)

all_yr_aggr %>%
  group_by(CH_METHOD) %>%
  summarize(non_na_count = sum(!is.na(DELTA_FCH4)))

## delta_FCH4 medians grouped by chamber method

# weekly

all_week_aggr %>%
  group_by(CH_METHOD) %>%
  summarise(
    median_DELTA_FCH4 = median(DELTA_FCH4, na.rm = TRUE),
    CV_DELTA_FCH4 = raster::cv(DELTA_FCH4, na.rm = TRUE)
  )

# monthly
all_month_aggr %>%
  group_by(CH_METHOD) %>%
  summarise(
    median_DELTA_FCH4 = median(DELTA_FCH4, na.rm = TRUE),
    CV_DELTA_FCH4 = raster::cv(DELTA_FCH4, na.rm = TRUE)
  )

#annual
all_yr_aggr %>%
  group_by(CH_METHOD) %>%
  summarise(
    median_DELTA_FCH4 = median(DELTA_FCH4, na.rm = TRUE),
    CV_DELTA_FCH4 = raster::cv(DELTA_FCH4, na.rm = TRUE)
  )

# daily
all_d_aggr %>%
  group_by(CH_METHOD) %>%
  summarise(
    median_DELTA_FCH4 = median(DELTA_FCH4, na.rm = TRUE),
    CV_DELTA_FCH4 = raster::cv(DELTA_FCH4, na.rm = TRUE)
  )

####################################################################################################################################################################

## Heatmaps of diel delta_FCH4 trends (Figs. B4-B5)

# add Month_letter column
all_hh_aggr <- all_hh_aggr %>% mutate(Month_letter =
                                                                case_when(month(all_hh_aggr$TIMESTAMP_START) == 06 ~ "Jun", 
                                                                          month(all_hh_aggr$TIMESTAMP_START) == 07 ~ "Jul",
                                                                          month(all_hh_aggr$TIMESTAMP_START) == 08 ~ "Aug",
                                                                          month(all_hh_aggr$TIMESTAMP_START) == 09 ~ "Sep",
                                                                          month(all_hh_aggr$TIMESTAMP_START) == 10 ~ "Oct",
                                                                          month(all_hh_aggr$TIMESTAMP_START) == 11 ~ "Nov",
                                                                          month(all_hh_aggr$TIMESTAMP_START) == 12 ~ "Dec",
                                                                          month(all_hh_aggr$TIMESTAMP_START) == 01 ~ "Jan",
                                                                          month(all_hh_aggr$TIMESTAMP_START) == 02 ~ "Feb",
                                                                          month(all_hh_aggr$TIMESTAMP_START) == 03 ~ "Mar",
                                                                          month(all_hh_aggr$TIMESTAMP_START) == 04 ~ "Apr",
                                                                          month(all_hh_aggr$TIMESTAMP_START) == 05 ~ "May")
)


all_hr_aggr <- all_hr_aggr %>% mutate(Month_letter =
                                        case_when(MONTH == 6 ~ "Jun", 
                                                  MONTH == 7 ~ "Jul",
                                                  MONTH == 8 ~ "Aug",
                                                  MONTH == 9 ~ "Sep",
                                                  MONTH == 10 ~ "Oct",
                                                  MONTH == 11 ~ "Nov",
                                                  MONTH == 12 ~ "Dec",
                                                  MONTH == 1 ~ "Jan",
                                                  MONTH == 2 ~ "Feb",
                                                  MONTH == 3 ~ "Mar",
                                                  MONTH == 4 ~ "Apr",
                                                  MONTH == 5 ~ "May")
)

all_hh_aggr$HOUR <- factor(all_hh_aggr$HOUR, levels = as.character(0:23))
all_hr_aggr$HOUR <- factor(all_hr_aggr$HOUR, levels = as.character(0:23))

# all sites, half-hourly
all_hh_aggr$Month_letter <- with(all_hh_aggr,factor(Month_letter,levels = c("Feb", "Mar", "Apr", "May", "Jun",
                                                                                                    "Jul", "Aug", "Sep", "Oct", "Nov")))

# Summarize data
data_summary_auto <- all_hh_aggr %>%
  group_by(HOUR, Month_letter) %>%
  summarize(median_DELTA_FCH4 = median(DELTA_FCH4), .groups = 'drop')

# filter out months that are not in every site
data_filtered <- data_summary_auto %>%
  dplyr::filter(!Month_letter %in% c("Feb", "Mar", "Apr", "Nov"))

# half-hourly plot
heatmap_auto <- ggplot(data_filtered, aes(as.factor(HOUR), as.factor(Month_letter))) +
  geom_tile(aes(fill = median_DELTA_FCH4)) +
  scale_fill_continuous_divergingx(palette = 'RdBu',
                                   name = expression(paste(Delta*"FCH"[4]," (nmol m"^-2, "s"^-1, ")"))) + 
  ggtitle("All sites") +
  theme_bw() + xlab("HOUR") + ylab("MONTH")+
  theme(plot.title = element_text(hjust=0.5, size = 20),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size=20),
        axis.text.x = element_text(size=20, angle=90, hjust=1),
        axis.text.y = element_text(size=20),
        axis.text = element_text(size = 20),
        legend.title = element_blank(),
        legend.text = element_text(size=20))
heatmap_auto  

# Kruskal-Wallis test for each month
kruskal_results <- all_hh_aggr %>%
  filter(!is.na(DELTA_FCH4)) %>%         # Remove rows with NA in DELTA_FCH4
  group_by(Month_letter) %>%            # Group by month
  nest() %>%                            # Nest the data for each month
  mutate(kruskal_test = purrr::map(data, ~kruskal.test(DELTA_FCH4 ~ HOUR, data = .x))) %>%  # Apply Kruskal-Wallis test
  mutate(p_value = purrr::map_dbl(kruskal_test, ~.x$p.value)) %>%  # Extract p-values
  ungroup() %>%                         # Ungroup to avoid issues with select()
  dplyr::select(Month_letter, p_value)         # Select the relevant columns

# View the results
kruskal_results


# heatmap for all sites, hourly aggregation

all_hr_aggr$Month_letter <- with(all_hr_aggr,factor(Month_letter,levels = c("Feb", "Mar", "Apr", "May", "Jun",
                                                                            "Jul", "Aug", "Sep", "Oct", "Nov")))
all_hr_aggr$Month_letter <- as.factor(all_hr_aggr$Month_letter)

# summarize df
data_summary_hr <- all_hr_aggr %>%
  group_by(HOUR, Month_letter) %>%
  summarize(median_DELTA_FCH4 = median(DELTA_FCH4), .groups = 'drop')

# Filter out months that are not included in all sites
data_filtered_hr <- data_summary_hr %>%
  dplyr::filter(!Month_letter %in% c("Feb", "Mar", "Apr", "Nov"))

# hourly heatmap
heatmap_hr <- ggplot(data_filtered_hr, aes(as.factor(HOUR), as.factor(Month_letter))) +    
  geom_tile(aes(fill = median_DELTA_FCH4)) +
  scale_fill_continuous_divergingx(palette = 'RdBu', 
                                   name = expression(paste(Delta*"FCH"[4]," (nmol m"^-2, "s"^-1, ")"))) + 
  ggtitle("All sites") +
  theme_bw() + xlab("HOUR") + ylab("MONTH")+
  theme(plot.title = element_text(hjust=0.5, size = 20),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size=20),
        axis.text.x = element_text(size=20, angle=90, hjust=1),
        axis.text.y = element_text(size=20),
        axis.text = element_text(size = 20),
        legend.title = element_blank(),
        legend.text = element_text(size=20))
heatmap_hr  

# heatmaps for individual sites

# half-hourly
hh_SEDEG <- subset(all_hh_aggr[all_hh_aggr$SITE == "SE-DEG", ])
# remove NA from DELTA_FCH4 because otherwise the heatmap includes NAs as grey cells
hh_SEDEG <- hh_SEDEG %>% drop_na(DELTA_FCH4)
hh_SEDEG$Month_letter <- with(hh_SEDEG,factor(Month_letter,levels = c("Feb", "Mar", "Apr", "May", "Jun",
                                                                          "Jul", "Aug", "Sep", "Oct", "Nov")))

# summarize df
data_summary_hh_SEDEG <- hh_SEDEG %>%
  group_by(HOUR, Month_letter) %>%
  summarize(median_DELTA_FCH4 = median(DELTA_FCH4), .groups = 'drop')

# plot heatmap for SE-Deg
heatmap_hh_SEDEG <- ggplot(data_summary_hh_SEDEG, aes(as.factor(HOUR), as.factor(Month_letter))) +               
  geom_tile(aes(fill = median_DELTA_FCH4)) +
  scale_fill_continuous_divergingx(palette = 'RdBu', 
                                   name = expression(paste(Delta*"FCH"[4]," (nmol m"^-2, "s"^-1, ")"))) + ggtitle("SE-Deg") +
  theme_linedraw() + xlab("HOUR") + ylab("MONTH") +
  theme(plot.title = element_text(hjust=0.5, size = 20),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size=20),
        axis.text.x = element_text(size=20, angle=90, hjust=1),
        axis.text.y = element_text(size=20),
        axis.text = element_text(size = 20),
        legend.title = element_blank(),
        legend.text = element_text(size=20))
heatmap_hh_SEDEG  


# Kruskal-Wallis test for each month
kruskal_results <-hh_SEDEG %>%
  filter(!is.na(DELTA_FCH4)) %>%         # Remove rows with NA in DELTA_FCH4
  group_by(Month_letter) %>%            # Group by month
  nest() %>%                            # Nest the data for each month
  mutate(kruskal_test = purrr::map(data, ~kruskal.test(DELTA_FCH4 ~ HOUR, data = .x))) %>%  # Apply Kruskal-Wallis test
  mutate(p_value = purrr::map_dbl(kruskal_test, ~.x$p.value)) %>%  # Extract p-values
  ungroup() %>%                         # Ungroup to avoid issues with select()
  dplyr::select(Month_letter, p_value)         # Select the relevant columns

# View the results
kruskal_results


# hourly aggregation - SE-Deg

hr_SEDEG <- subset(all_hr_aggr[all_hr_aggr$SITE == "SE-DEG", ])
# remove NA from DELTA_FCH4 because otherwise the heatmap includes NAs as grey cells
hr_SEDEG <- hr_SEDEG %>% drop_na(DELTA_FCH4)
hr_SEDEG$Month_letter <- with(hr_SEDEG,factor(Month_letter,levels = c("Feb", "Mar", "Apr", "May", "Jun",
                                                                      "Jul", "Aug", "Sep", "Oct", "Nov")))
data_summary_hr_SEDEG <- hr_SEDEG %>%
  group_by(HOUR, Month_letter) %>%
  summarize(median_DELTA_FCH4 = median(DELTA_FCH4), .groups = 'drop')

# plot SE-Deg hourly heatmap
heatmap_hr_SEDEG <- ggplot(data_summary_hr_SEDEG, aes(as.factor(HOUR), as.factor(Month_letter))) +  
  geom_tile(aes(fill = median_DELTA_FCH4)) +
  scale_fill_continuous_divergingx(palette = 'RdBu', 
                                   name = expression(paste(Delta*"FCH"[4]," (nmol m"^-2, "s"^-1, ")"))) + ggtitle("SE-Deg") +
  theme_linedraw() + xlab(" ") + ylab(" ") +
  theme(plot.title = element_text(hjust=0.5, size = 20),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size=20),
        axis.text.x = element_text(size=20, angle=90, hjust=1),
        axis.text.y = element_text(size=20),
        axis.text = element_text(size = 20),
        legend.title = element_blank(),
        legend.text = element_text(size=20))
heatmap_hr_SEDEG  


# half-hourly heatmap for CN-Hgu

hh_CNHGU <- subset(all_hh_aggr[all_hh_aggr$SITE == "CN-HGU", ])
# remove NA from DELTA_FCH4 because otherwise the heatmap includes NAs as grey cells
hh_CNHGU <- hh_CNHGU %>% drop_na(DELTA_FCH4)
hh_CNHGU$Month_letter <- with(hh_CNHGU,factor(Month_letter,levels = c("Feb", "Mar", "Apr", "May", "Jun",
                                                                          "Jul", "Aug", "Sep", "Oct", "Nov")))

data_summary_hh_CNHGU <- hh_CNHGU %>%
  group_by(HOUR, Month_letter) %>%
  summarize(median_DELTA_FCH4 = median(DELTA_FCH4), .groups = 'drop')

# plot
heatmap_hh_CNHGU <- ggplot(data_summary_hh_CNHGU, aes(as.factor(HOUR), as.factor(Month_letter))) +          
  geom_tile(aes(fill = median_DELTA_FCH4)) +
  scale_fill_continuous_divergingx(palette = 'RdBu', 
                                   name = expression(paste(Delta*"FCH"[4]," (nmol m"^-2, "s"^-1, ")")))+
  xlab("HOUR") + ylab("MONTH") + ggtitle("CN-Hgu") + theme_linedraw() +
  theme(plot.title = element_text(hjust=0.5, size = 20),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size=20),
        axis.text.x = element_text(size=20, angle=90, hjust=1),
        axis.text.y = element_text(size=20),
        axis.text = element_text(size = 20),
        legend.title = element_blank(),
        legend.text = element_text(size=20))
heatmap_hh_CNHGU


# hourly heatmap for CN-Hgu

hr_CNHGU <- subset(all_hr_aggr[all_hr_aggr$SITE == "CN-HGU", ])
# remove NA from DELTA_FCH4 because otherwise the heatmap includes NAs as grey cells
hr_CNHGU <- hr_CNHGU %>% drop_na(DELTA_FCH4)
hr_CNHGU$Month_letter <- with(hr_CNHGU,factor(Month_letter,levels = c("Feb", "Mar", "Apr", "May", "Jun",
                                                                      "Jul", "Aug", "Sep", "Oct", "Nov")))

data_summary_hr_CNHGU <- hr_CNHGU %>%
  group_by(HOUR, Month_letter) %>%
  summarize(median_DELTA_FCH4 = median(DELTA_FCH4), .groups = 'drop')

# Perform Kruskal-Wallis test for each month
hh_CNHGU$HOUR <- as.factor(hh_CNHGU$HOUR)
hh_CNHGU$Month_letter <- as.factor(hh_CNHGU$Month_letter)

kruskal_results <- hh_CNHGU %>%
  filter(!is.na(DELTA_FCH4)) %>%         # Remove rows with NA in DELTA_FCH4
  group_by(Month_letter) %>%            # Group by month
  nest() %>%                            # Nest the data for each month
  mutate(kruskal_test = purrr::map(data, ~kruskal.test(DELTA_FCH4 ~ HOUR, data = .x))) %>%  # Apply Kruskal-Wallis test
  mutate(p_value = purrr::map_dbl(kruskal_test, ~.x$p.value)) %>%  # Extract p-values
  ungroup() %>%                         # Ungroup to avoid issues with select()
  dplyr::select(Month_letter, p_value)         # Select the relevant columns

# View the results
kruskal_results

# plot hourly CN-Hgu heatmap

heatmap_hr_CNHGU <- ggplot(data_summary_hr_CNHGU, aes(as.factor(HOUR), as.factor(Month_letter))) +    
  geom_tile(aes(fill = median_DELTA_FCH4)) +
  scale_fill_continuous_divergingx(palette = 'RdBu', 
                                   name = expression(paste(Delta*"FCH"[4]," (nmol m"^-2, "s"^-1, ")")))+
  xlab("HOUR") + ylab("MONTH") + ggtitle("CN-Hgu") + theme_linedraw() +
  theme(plot.title = element_text(hjust=0.5, size = 20),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size=20),
        axis.text.x = element_text(size=20, angle=90, hjust=1),
        axis.text.y = element_text(size=20),
        axis.text = element_text(size = 20),
        legend.title = element_blank(),
        legend.text = element_text(size=20))
heatmap_hr_CNHGU


# half-hourly heatmap for US-Ho1

hh_USHO1 <- subset(all_hh_aggr[all_hh_aggr$SITE == "US-HO1", ])
# remove NA from DELTA_FCH4 because otherwise the heatmap includes NAs as grey cells
hh_USHO1$Month_letter <- with(hh_USHO1,factor(Month_letter,levels = c("Feb", "Mar", "Apr", "May", "Jun",
                                                                          "Jul", "Aug", "Sep", "Oct", "Nov")))

data_summary_hh_USHO1 <- hh_USHO1 %>%
  group_by(HOUR, Month_letter) %>%
  summarize(median_DELTA_FCH4 = median(DELTA_FCH4), .groups = 'drop')


# plot
heatmap_hh_USHO1 <- ggplot(data_summary_hh_USHO1, aes(as.factor(HOUR), as.factor(Month_letter))) +            
  geom_tile(aes(fill = median_DELTA_FCH4)) +
  scale_fill_continuous_divergingx(palette = 'RdBu',
                                   name = expression(paste(Delta*"FCH"[4]," (nmol m"^-2, "s"^-1, ")"))) + ggtitle("US-Ho1")+
  theme_bw() + xlab("HOUR") + ylab("MONTH") +
  theme(plot.title = element_text(hjust=0.5, size = 20),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size=20),
        axis.text.x = element_text(size=20, angle=90, hjust=1),
        axis.text.y = element_text(size=20),
        axis.text = element_text(size = 20),
        legend.title = element_blank(),
        legend.text = element_text(size=20))
heatmap_hh_USHO1

# Kruskal-Wallis test for each month
kruskal_results <- hh_USHO1 %>%
  filter(!is.na(DELTA_FCH4)) %>%         # Remove rows with NA in DELTA_FCH4
  group_by(Month_letter) %>%            # Group by month
  nest() %>%                            # Nest the data for each month
  mutate(kruskal_test = purrr::map(data, ~kruskal.test(DELTA_FCH4 ~ HOUR, data = .x))) %>%  # Apply Kruskal-Wallis test
  mutate(p_value = purrr::map_dbl(kruskal_test, ~.x$p.value)) %>%  # Extract p-values
  ungroup() %>%                         # Ungroup to avoid issues with select()
  dplyr::select(Month_letter, p_value)         # Select the relevant columns

# View the results
kruskal_results

# hourly aggregation US-Ho1 heatmap

hr_USHO1 <- subset(all_hr_aggr[all_hr_aggr$SITE == "US-HO1", ])
# remove NA from DELTA_FCH4 because otherwise the heatmap includes NAs as grey cells
hr_USHO1$Month_letter <- with(hr_USHO1,factor(Month_letter,levels = c("Feb", "Mar", "Apr", "May", "Jun",
                                                                      "Jul", "Aug", "Sep", "Oct", "Nov")))

data_summary_hr_USHO1 <- hr_USHO1 %>%
  group_by(HOUR, Month_letter) %>%
  summarize(median_DELTA_FCH4 = median(DELTA_FCH4), .groups = 'drop')

# plot US-Ho1 hourly aggregation heatmap
heatmap_hr_USHO1 <- ggplot(data_summary_hr_USHO1, aes(as.factor(HOUR), as.factor(Month_letter))) +              
  geom_tile(aes(fill = median_DELTA_FCH4)) +
  scale_fill_continuous_divergingx(palette = 'RdBu', 
                                   name = expression(paste(Delta*"FCH"[4]," (nmol m"^-2, "s"^-1, ")"))) + ggtitle("US-Ho1")+
  theme_bw() + xlab("HOUR") + ylab("MONTH")+
  theme(plot.title = element_text(hjust=0.5, size = 20),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size=20),
        axis.text.x = element_text(size=20, angle=90, hjust=1),
        axis.text.y = element_text(size=20),
        axis.text = element_text(size = 20),
        legend.title = element_blank(),
        legend.text = element_text(size=20))
heatmap_hr_USHO1


# US-Uaf half-hourly heatmap

hh_USUAF <- subset(all_hh_aggr[all_hh_aggr$SITE == "US-UAF", ])
# remove NA from DELTA_FCH4 because otherwise the heatmap includes NAs as grey cells
hh_USUAF <- hh_USUAF %>% drop_na(DELTA_FCH4)
hh_USUAF$Month_letter <- with(hh_USUAF,factor(Month_letter,levels = c("Feb", "Mar", "Apr", "May", "Jun",
                                                                          "Jul", "Aug", "Sep", "Oct", "Nov")))

data_summary_hh_USUAF <- hh_USUAF %>%
  group_by(HOUR, Month_letter) %>%
  summarize(median_DELTA_FCH4 = median(DELTA_FCH4), .groups = 'drop')


# plot half-hourly heatmap for US-Uaf
heatmap_hh_USUAF <- ggplot(data_summary_hh_USUAF, aes(as.factor(HOUR), as.factor(Month_letter))) +              
  geom_tile(aes(fill = median_DELTA_FCH4)) +
  scale_fill_continuous_divergingx(palette = 'RdBu', 
                                   name=expression(paste(Delta*"FCH"[4]," (nmol m"^-2, "s"^-1, ")"))) +
  ggtitle("US-Uaf") + theme_bw() + xlab(" ") + ylab(" ")+
  theme(plot.title = element_text(hjust=0.5, size = 20),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size=20),
        axis.text.x = element_text(size=20, angle=90, hjust=1),
        axis.text.y = element_text(size=20),
        axis.text = element_text(size = 20),
        legend.title = element_blank(),
        legend.text = element_text(size=20))

heatmap_hh_USUAF

# Kruskal-Wallis test for each month
kruskal_results <- hh_USUAF %>%
  filter(!is.na(DELTA_FCH4)) %>%         # Remove rows with NA in DELTA_FCH4
  group_by(Month_letter) %>%            # Group by month
  nest() %>%                            # Nest the data for each month
  mutate(kruskal_test = purrr::map(data, ~kruskal.test(DELTA_FCH4 ~ HOUR, data = .x))) %>%  # Apply Kruskal-Wallis test
  mutate(p_value = purrr::map_dbl(kruskal_test, ~.x$p.value)) %>%  # Extract p-values
  ungroup() %>%                         # Ungroup to avoid issues with select()
  dplyr::select(Month_letter, p_value)         # Select the relevant columns

# View the results
kruskal_results

# Hourly aggregation heatmap for US-Uaf

hr_USUAF <- subset(all_hr_aggr[all_hr_aggr$SITE == "US-UAF", ])
# remove NA from DELTA_FCH4 because otherwise the heatmap includes NAs as grey cells
hr_USUAF <- hr_USUAF %>% drop_na(DELTA_FCH4)
hr_USUAF$Month_letter <- with(hr_USUAF,factor(Month_letter,levels = c("Feb", "Mar", "Apr", "May", "Jun",
                                                                      "Jul", "Aug", "Sep", "Oct", "Nov")))

data_summary_hr_USUAF <- hr_USUAF %>%
  group_by(HOUR, Month_letter) %>%
  summarize(median_DELTA_FCH4 = median(DELTA_FCH4), .groups = 'drop')


# plot hourly aggregation heatmap for US-Uaf

heatmap_hr_USUAF <- ggplot(data_summary_hr_USUAF, aes(as.factor(HOUR), as.factor(Month_letter))) +       
  geom_tile(aes(fill = median_DELTA_FCH4)) +
  scale_fill_continuous_divergingx(palette = 'RdBu',# mid = 0, p3 = 0.9, p4=0.9, 
                                   name=expression(paste(Delta*"FCH"[4]," (nmol m"^-2, "s"^-1, ")"))) +
  ggtitle("US-Uaf") + theme_linedraw() + xlab("HOUR") + ylab("MONTH")+
  theme(plot.title = element_text(hjust=0.5, size = 20),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size=20),
        axis.text.x = element_text(size=20, angle=90, hjust=1),
        axis.text.y = element_text(size=20),
        axis.text = element_text(size = 20),
        legend.title = element_blank(),
        legend.text = element_text(size=20))

heatmap_hr_USUAF


# combine heatmaps together into one figure (the figure was further modified in Inkscape)

# half-hourly (Fig. B4)
require(grid)
heatmap_hh_comb <- ggarrange(heatmap_hh + rremove("ylab") + rremove("xlab"), heatmap_hh_CNHGU + rremove("ylab") + rremove("xlab"), heatmap_hh_SEDEG + rremove("ylab") + rremove("xlab"), 
                             heatmap_hh_USHO1 + rremove("ylab") + rremove("xlab"), heatmap_hh_USUAF + rremove("ylab") + rremove("xlab"),
                             ncol = 2, nrow = 3, common.legend = F,
                             align="hv")

heatmap_hh_comb <- annotate_figure(heatmap_hh_comb, left = textGrob("MONTH", rot = 90, vjust = 1),
                                   bottom = textGrob("HOUR"))
heatmap_hh_comb

# hourly (Fig. B5)

heatmap_hr_comb <- ggarrange(heatmap_hr + rremove("ylab") + rremove("xlab"), 
                             heatmap_hr_CNHGU + rremove("ylab") + rremove("xlab"), 
                             heatmap_hr_SEDEG + rremove("ylab") + rremove("xlab"), 
                             heatmap_hr_USHO1 + rremove("ylab") + rremove("xlab"), 
                             heatmap_hr_USUAF + rremove("ylab") + rremove("xlab"),
                             ncol = 2, nrow = 3, common.legend = F,
                             align = "hv")

heatmap_hr_comb <- annotate_figure(heatmap_hr_comb, left = textGrob("MONTH", rot = 90, vjust = 1),
                                   bottom = textGrob("HOUR"))
heatmap_hr_comb


################################################################################################################################################################
##### LINEAR MIXED MODELS #####

# create a DELTA_FCH4_abs column

all_d_aggr$DELTA_FCH4_abs <- abs(all_d_aggr$DELTA_FCH4)
all_hh_aggr$DELTA_FCH4_abs <- abs(all_hh_aggr$DELTA_FCH4)
all_week_aggr$DELTA_FCH4_abs <- abs(all_week_aggr$DELTA_FCH4)
all_month_aggr$DELTA_FCH4_abs <- abs(all_month_aggr$DELTA_FCH4)
all_yr_aggr$DELTA_FCH4_abs <- abs(all_yr_aggr$DELTA_FCH4)
all_hr_aggr$DELTA_FCH4_abs <- abs(all_hr_aggr$DELTA_FCH4)

#####

# for models, exclude CN-HGU (which does not have water table depth WTL data) and remove rows with NAs

# remove CN-HGU

all_sites_preds_d_noCNHGU <- all_d_aggr %>%
  filter(SITE != "CN-HGU")

hh_sites_preds_chaggr_noCNHGU <- all_hh_aggr %>%
  filter(SITE != "CN-HGU")
 
all_hr_aggr_noCNHGU <- all_hr_aggr %>%
  filter(!SITE %in% c("CN-HGU", "US-OWC", "FI-SI2", "US-LA1", "US-LA2", "US-LOS"))

all_week_aggr_noCNHGU <- all_week_aggr %>%
  filter(SITE != "CN-HGU")
 
all_month_aggr_noCNHGU <- all_month_aggr %>%
  filter(SITE != "CN-HGU")

### remove NAs

preds_aggr <- c("SITE_TS_TOP", 
                "SITE_WTL", 
                "USTAR_MEDIAN", 
                "U_WIND_F_MEAN", 
                "V_WIND_F_MEAN",
                "PA_F_MEDIAN",  
                "NEE_F_MEAN", 
                "VPD_F_MEDIAN"
)

preds_hh <- c("SITE_TS_TOP", 
                "SITE_WTL",
                "EC_USTAR", 
                "U_WIND_F", 
                "V_WIND_F",
                "EC_PA_F",  
                "NEE_F", 
                "EC_VPD_F"
)


hh_noNA_noCNHGU <- hh_sites_preds_chaggr_noCNHGU[complete.cases(hh_sites_preds_chaggr_noCNHGU[, preds_hh]), ]
d_noNA_noCNHGU <- all_sites_preds_d_noCNHGU[complete.cases(all_sites_preds_d_noCNHGU[, preds_aggr]), ]
hr_noNA_noCNHGU <- all_hr_aggr_noCNHGU[complete.cases(all_hr_aggr_noCNHGU[, preds_aggr]), ]
week_noNA_noCNHGU <- all_week_aggr_noCNHGU[complete.cases(all_week_aggr_noCNHGU[, preds_aggr]), ]
month_noNA_noCNHGU <- all_month_aggr_noCNHGU[complete.cases(all_month_aggr_noCNHGU[, preds_aggr]), ]

# Apply centering and scaling
d_noNA_noCNHGU <- d_noNA_noCNHGU %>%
  mutate(across(all_of(preds_aggr), 
                list(ctr = ~ (. - mean(.)) / sd(.)),
                .names = "{col}_ctr"))

hh_noNA_noCNHGU <- hh_noNA_noCNHGU %>%
  mutate(across(all_of(preds_hh), 
                list(ctr = ~ (. - mean(.)) / sd(.)),
                .names = "{col}_ctr"))

hr_noNA_noCNHGU <- hr_noNA_noCNHGU %>%
  mutate(across(all_of(preds_aggr), 
                list(ctr = ~ (. - mean(.)) / sd(.)),
                .names = "{col}_ctr"))

week_noNA_noCNHGU <- week_noNA_noCNHGU %>%
  mutate(across(all_of(preds_aggr), 
                list(ctr = ~ (. - mean(.)) / sd(.)),
                .names = "{col}_ctr"))
month_noNA_noCNHGU <- month_noNA_noCNHGU %>%
  mutate(across(all_of(preds_aggr), 
                list(ctr = ~ (. - mean(.)) / sd(.)),
                .names = "{col}_ctr"))


# get the predictors

preds_d <- d_noNA_noCNHGU[, c("SITE_TS_TOP_ctr", 
                                      "SITE_WTL_ctr", 
                                      "USTAR_MEDIAN_ctr", 
                                      "U_WIND_F_MEAN_ctr", 
                                      "V_WIND_F_MEAN_ctr",
                                      "PA_F_MEDIAN_ctr",  
                                      "NEE_F_MEAN_ctr", 
                                      "VPD_F_MEDIAN_ctr"
)]
colnames(preds_d) <- c("TS", "WTL",
                           "USTAR", 
                           "u.wind",
                           "v.wind",
                           "PA", 
                           "NEE", 
                           "VPD"
)


preds_hh <- hh_noNA_noCNHGU[, c("SITE_TS_TOP_ctr", 
                                    "SITE_WTL_ctr", 
                                    "EC_USTAR_ctr", 
                                    "U_WIND_F_ctr", 
                                    "V_WIND_F_ctr",
                                    "EC_PA_F_ctr",  
                                    "NEE_F_ctr", 
                                    "EC_VPD_F_ctr"
)]
colnames(preds_hh) <- c("TS", 
                          "WTL",
                          "USTAR", 
                          "u.wind",
                          "v.wind",
                          "PA", 
                          "NEE", 
                          "VPD"
)
correlation_matrix <- cor(preds_d, use = "pairwise.complete.obs")
correlation_matrix <- cor(preds_hh, use = "pairwise.complete.obs")

# Round correlation matrix for clarity
correlation_matrix <- round(correlation_matrix, 2)

# Visualize correlation matrix using corrplot
corrplot(correlation_matrix, method = "color", type = "upper", order = "hclust", 
         tl.col = "black", tl.srt = 45, diag = FALSE, addCoef.col = "black")

## check histograms

hist(d_noNA_noCNHGU$DELTA_FCH4_abs)
hist(log(d_noNA_noCNHGU$DELTA_FCH4_abs))

all_sites_preds_d_noCNHGU$DELTA_FCH4_abs <- abs(all_sites_preds_d_noCNHGU$DELTA_FCH4)

# Yeo-Johnson transformations for absolute delta_FCH4
yeo_johnson_transformation <- yeojohnson(abs(d_noNA_noCNHGU$DELTA_FCH4))

# Transform the data
d_noNA_noCNHGU$DELTA_FCH4_abs_yjtrans <- predict(yeo_johnson_transformation)

hist(d_noNA_noCNHGU$DELTA_FCH4_abs_yjtrans)
qqPlot(log(d_noNA_noCNHGU$DELTA_FCH4_abs))

d_noNA_noCNHGU$DOMINANT_VEGETATION <- as.factor(d_noNA_noCNHGU$DOMINANT_VEGETATION)
d_noNA_noCNHGU$MONTH <- as.factor(d_noNA_noCNHGU$MONTH)
d_noNA_noCNHGU$SITE <- as.factor(d_noNA_noCNHGU$SITE)

# center and scale
d_noNA_noCNHGU$DELTA_FCH4_abs_yjtrans_ctr <- scale(d_noNA_noCNHGU$DELTA_FCH4_abs_yjtrans)

# first model version with YJ-transformed delta_FCH4
m_yj <- lme(
  DELTA_FCH4_abs_yjtrans_ctr ~ SITE_TS_TOP_ctr +
    SITE_WTL_ctr +
    USTAR_MEDIAN_ctr +
    U_WIND_F_MEAN_ctr +
    V_WIND_F_MEAN_ctr +
    PA_F_MEDIAN_ctr +
    NEE_F_MEAN_ctr +
    VPD_F_MEDIAN_ctr +
    DOMINANT_VEGETATION + 
    MONTH, 
  random = ~1 | SITE,  
  data = d_noNA_noCNHGU
)

# model residuals
residuals_yj <- residuals(m_yj, type = "pearson")

hist(residuals_yj, main = "Histogram of YJTransformed Residuals")

qqnorm(residuals_yj, main = "Q-Q Plot for YJ-Transformed Model")
qqline(residuals_yj)

plot(fitted(m_yj), residuals_yj, main = "Residuals vs Fitted (YJ Model)")
abline(h = 0, col = "red")

AIC(m_yj)

# half-hourly

hh_noNA_noCNHGU$HOUR <- as.factor(hh_noNA_noCNHGU$HOUR)
hh_noNA_noCNHGU$MONTH <- as.factor(hh_noNA_noCNHGU$MONTH)

# YJ transformation for delta_FCH4

yeo_johnson_transformation <- yeojohnson(hh_noNA_noCNHGU$DELTA_FCH4_abs)

# Transform the data
hh_noNA_noCNHGU$DELTA_FCH4_abs_yjtrans <- predict(yeo_johnson_transformation)

hist(hh_noNA_noCNHGU$DELTA_FCH4_abs_yjtrans)
qqPlot(hh_noNA_noCNHGU$DELTA_FCH4_abs_yjtrans)
qqPlot(hh_noNA_noCNHGU$DELTA_FCH4_abs)

# model residuals
residuals_yj <- residuals(m_yj, type = "pearson")

hist(residuals_yj, main = "Histogram of YJTransformed Residuals")

qqline(residuals_log)

qqnorm(residuals_yj, main = "Q-Q Plot for YJ-Transformed Model")
qqline(residuals_yj)

plot(fitted(m_yj), residuals_yj, main = "Residuals vs Fitted (YJ Model)")
abline(h = 0, col = "red")

AIC(m_yj)

# center and scale
hh_noNA_noCNHGU$DELTA_FCH4_abs_yjtrans_ctr <- scale(hh_noNA_noCNHGU$DELTA_FCH4_abs_yjtrans)


# first model version with YJ transformed delta_FCH4 without accounting for temporal autocorrelation and residual heterogeneity
m_auto <- lme(
  DELTA_FCH4_abs_yjtrans_ctr ~ SITE_TS_TOP_ctr +
    SITE_WTL_ctr +
    EC_USTAR_ctr +
    U_WIND_F_ctr +
    V_WIND_F_ctr +
    EC_PA_F_ctr +
    NEE_F_ctr +
    EC_VPD_F_ctr +
    DOMINANT_VEGETATION + 
    HOUR, 
  random = ~1 | SITE, 
  data = hh_noNA_noCNHGU
)

vif(m_auto)


### The following code tries out different residual heterogeneity and temporal autocorrelation structures
### and compares them with AIC and residual diagnostics

## for models, create a new HalfHour variable

# Extract the half-hour index (1–48) from a timestamp
hh_noNA_noCNHGU$HalfHour <- substr(as.character(hh_noNA_noCNHGU$TIMESTAMP_START), 15, 16)
hh_noNA_noCNHGU$HalfHour[1] <- "00"
hh_noNA_noCNHGU$HalfHour[hh_noNA_noCNHGU$HalfHour == ""] <- "00"

hh_noNA_noCNHGU$HalfHour2 <- ifelse(hh_noNA_noCNHGU$HalfHour == "00", 1, 2)
hh_noNA_noCNHGU$HalfHour2 <- as.numeric(hh_noNA_noCNHGU$HalfHour2)

hh_noNA_noCNHGU$DATE <- as_date(hh_noNA_noCNHGU$DATE)

hh_noNA_noCNHGU <- hh_noNA_noCNHGU %>%
  mutate(DateHour = paste(DATE, HOUR, sep = "_"))

m_hh_datehour <- lme(
  DELTA_FCH4_abs_yjtrans_ctr ~ SITE_TS_TOP_ctr +
    SITE_WTL_ctr +
    EC_USTAR_ctr +
    U_WIND_F_ctr +
    V_WIND_F_ctr +
    EC_PA_F_ctr +
    NEE_F_ctr +
    EC_VPD_F_ctr +
    DOMINANT_VEGETATION + 
    HOUR, 
  random = ~1 | DateHour/SITE, 
  weights = varIdent(form = ~1 | MONTH), #doesn't converge
  data = hh_noNA_noCNHGU
)

AIC(m_hh_datehour, m_hh_AR1, m_auto)


m_hh_exp <- lme(
  DELTA_FCH4_abs_yjtrans ~ SITE_TS_TOP_ctr +
    SITE_WTL_ctr +
    EC_USTAR_ctr +
    U_WIND_F_ctr +
    V_WIND_F_ctr +
    EC_PA_F_ctr +
    NEE_F_ctr +
    EC_VPD_F_ctr +
    DOMINANT_VEGETATION + 
    HOUR, 
  random = ~1 | SITE, 
  weights = varComb(
    varExp(form = ~ EC_USTAR_ctr),       # Allow variance to vary by MONTH (categorical)
    varExp(form = ~ EC_PA_F_ctr)        # Exponential variance for EC_PA_F_ctr (continuous)
  ), 
  data = hh_noNA_noCNHGU
)


AIC(m_auto, m_hh_exp, m_hh_AR1, m_hh_sitedatehourre)

m_hh_expPA_datehour <- lme(
  DELTA_FCH4_abs_yjtrans ~ SITE_TS_TOP_ctr +
    SITE_WTL_ctr +
    EC_USTAR_ctr +
    U_WIND_F_ctr +
    V_WIND_F_ctr +
    EC_PA_F_ctr +
    NEE_F_ctr +
    EC_VPD_F_ctr +
    DOMINANT_VEGETATION + 
    HOUR, 
  random = ~1 | SITE/DateHour,  
  weights = varComb(     
    varExp(form = ~ EC_PA_F_ctr)        
  ), 
  data = hh_noNA_noCNHGU
)


m_hh_exp2_datehour <- lme(
  DELTA_FCH4_abs_yjtrans ~ SITE_TS_TOP_ctr +
    SITE_WTL_ctr +
    EC_USTAR_ctr +
    U_WIND_F_ctr +
    V_WIND_F_ctr +
    EC_PA_F_ctr +
    NEE_F_ctr +
    EC_VPD_F_ctr +
    DOMINANT_VEGETATION + 
    HOUR, 
  random = ~1 | SITE/DateHour,  
  weights = varComb(
    varExp(form = ~ EC_USTAR_ctr),       
    varExp(form = ~ EC_PA_F_ctr)        
  ), 
  data = hh_noNA_noCNHGU
)


m_hh_exp3_datehour <- lme(
  DELTA_FCH4_abs_yjtrans ~ SITE_TS_TOP_ctr +
    SITE_WTL_ctr +
    EC_USTAR_ctr +
    U_WIND_F_ctr +
    V_WIND_F_ctr +
    EC_PA_F_ctr +
    NEE_F_ctr +
    EC_VPD_F_ctr +
    DOMINANT_VEGETATION + 
    HOUR, 
  random = ~1 | SITE/DateHour, 
  weights = varComb(
    varExp(form = ~ EC_VPD_F_ctr),
    varExp(form = ~ EC_USTAR_ctr),       
    varExp(form = ~ EC_PA_F_ctr)        
  ), 
  data = hh_noNA_noCNHGU
)

AIC(m_auto, m_hh_exp_date, m_hh_exp, m_hh_exp2_datehour, m_hh_exp3_datehour)

m_hh_exp3_date <- lme(
  DELTA_FCH4_abs_yjtrans ~ SITE_TS_TOP_ctr +
    SITE_WTL_ctr +
    EC_USTAR_ctr +
    U_WIND_F_ctr +
    V_WIND_F_ctr +
    EC_PA_F_ctr +
    NEE_F_ctr +
    EC_VPD_F_ctr +
    DOMINANT_VEGETATION + 
    HOUR, 
  random = ~1 | SITE/DATE, 
  weights = varComb(
    varExp(form = ~ EC_VPD_F_ctr),
    varExp(form = ~ EC_USTAR_ctr),       
    varExp(form = ~ EC_PA_F_ctr)        
  ), 
  data = hh_noNA_noCNHGU
)

m_hh_date <- lme(
  DELTA_FCH4_abs_yjtrans ~ SITE_TS_TOP_ctr +
    SITE_WTL_ctr +
    EC_USTAR_ctr +
    U_WIND_F_ctr +
    V_WIND_F_ctr +
    EC_PA_F_ctr +
    NEE_F_ctr +
    EC_VPD_F_ctr +
    DOMINANT_VEGETATION + 
    HOUR, 
  random = ~1 | SITE/DATE, 
  data = hh_noNA_noCNHGU
)


AIC(m_auto, m_hh_exp3_date)


m_hh_exp3_datehour_AR1 <- lme(
  DELTA_FCH4_abs_yjtrans ~ SITE_TS_TOP_ctr +
    SITE_WTL_ctr +
    EC_USTAR_ctr +
    U_WIND_F_ctr +
    V_WIND_F_ctr +
    EC_PA_F_ctr +
    NEE_F_ctr +
    EC_VPD_F_ctr +
    DOMINANT_VEGETATION + 
    HOUR, 
  random = ~1 | SITE/DateHour, 
  correlation = corAR1(form = ~ HalfHour2 | SITE/DateHour),
  weights = varComb(
    varExp(form = ~ EC_VPD_F_ctr),
    varExp(form = ~ EC_USTAR_ctr),       
    varExp(form = ~ EC_PA_F_ctr)        
  ), 
  data = hh_noNA_noCNHGU
)


m_hh_expPAUSTAR_datehour_AR1 <- lme(
  DELTA_FCH4_abs_yjtrans ~ SITE_TS_TOP_ctr +
    SITE_WTL_ctr +
    EC_USTAR_ctr +
    U_WIND_F_ctr +
    V_WIND_F_ctr +
    EC_PA_F_ctr +
    NEE_F_ctr +
    EC_VPD_F_ctr +
    DOMINANT_VEGETATION + 
    HOUR, 
  random = ~1 | SITE/DateHour, 
  correlation = corAR1(form = ~ HalfHour2 | SITE/DateHour),
  weights = varComb(
    varExp(form = ~ EC_USTAR_ctr),       
    varExp(form = ~ EC_PA_F_ctr)        
  ), 
  data = hh_noNA_noCNHGU
)

m_hh_expPAVPD_datehour_AR1 <- lme(
  DELTA_FCH4_abs_yjtrans ~ SITE_TS_TOP_ctr +
    SITE_WTL_ctr +
    EC_USTAR_ctr +
    U_WIND_F_ctr +
    V_WIND_F_ctr +
    EC_PA_F_ctr +
    NEE_F_ctr +
    EC_VPD_F_ctr +
    DOMINANT_VEGETATION + 
    HOUR, 
  random = ~1 | SITE/DateHour, 
  correlation = corAR1(form = ~ HalfHour2 | SITE/DateHour),
  weights = varComb(
    varExp(form = ~ EC_VPD_F_ctr),     
    varExp(form = ~ EC_PA_F_ctr)        
  ), 
  data = hh_noNA_noCNHGU
)

m_hh_expUSTARVPD_datehour_AR1 <- lme(
  DELTA_FCH4_abs_yjtrans ~ SITE_TS_TOP_ctr +
    SITE_WTL_ctr +
    EC_USTAR_ctr +
    U_WIND_F_ctr +
    V_WIND_F_ctr +
    EC_PA_F_ctr +
    NEE_F_ctr +
    EC_VPD_F_ctr +
    DOMINANT_VEGETATION + 
    HOUR, 
  random = ~1 | SITE/DateHour, 
  correlation = corAR1(form = ~ HalfHour2 | SITE/DateHour),
  weights = varComb(
    varExp(form = ~ EC_VPD_F_ctr),
    varExp(form = ~ EC_USTAR_ctr)        
  ), 
  data = hh_noNA_noCNHGU
)


AIC(m_auto, m_hh_AR1, m_hh_expPA_datehour, m_hh_exp2_datehour, m_hh_exp3_datehour,
    m_hh_expPAVPD_datehour_AR1, m_hh_expPAUSTAR_datehour_AR1, m_hh_expUSTARVPD_datehour_AR1,
    m_hh_exp3_datehour_AR1)
# exp3 datehour AR1 is the best

# some residual het gen remains in VPD, PA, USTAR

# Extract residuals and fitted values (these were used one by one to look at residuals for each model)
E1 <- residuals(m_auto, type = "pearson")
F1 <- fitted(m_auto)

E1 <- residuals(m_hh_AR1, type = "pearson")
F1 <- fitted(m_hh_AR1)

E1 <- residuals(m_hh_exp, type = "pearson")
F1 <- fitted(m_hh_exp)

E1 <- residuals(m_hh_sitedatehourre, type = "pearson")
F1 <- fitted(m_hh_sitedatehourre)

E1 <- residuals(m_hh_exp3_datehour_AR1, type = "pearson")
F1 <- fitted(m_hh_exp3_datehour_AR1)

E1 <- residuals(m_hh_exp3_datehour, type = "pearson")
F1 <- fitted(m_hh_exp3_datehour) 

E1 <- residuals(m_hh_exp3_date, type = "pearson")
F1 <- fitted(m_hh_exp3_date)

E1 <- residuals(m_hh_expPAVPD_datehour_AR1, type = "pearson")
F1 <- fitted(m_hh_expPAVPD_datehour_AR1)


# Diagnostic plots
par(mfrow = c(2, 2), mar = c(5, 5, 2, 2))

# Residuals vs Fitted
plot(F1, E1, xlab = "Fitted values", ylab = "Residuals", main = "Residuals vs Fitted")
abline(h = 0, lty = 2)

# Residuals by SITE
plot(as.factor(hh_noNA_noCNHGU$SITE), E1, xlab = "SITE", ylab = "Residuals", main = "Residuals by SITE")
abline(h = 0, lty = 2)

# Residuals by SITE_WTL
plot(hh_noNA_noCNHGU$SITE_WTL_ctr, E1, xlab = "SITE_WTL", ylab = "Residuals", main = "Residuals by SITE_WTL")
abline(h = 0, lty = 2)

# res by site ts top
plot(hh_noNA_noCNHGU$SITE_TS_TOP_ctr, E1, xlab = "SITE_TS_TOP", ylab = "Residuals", main = "Residuals by SITE_TS_TOP")
abline(h = 0, lty = 2)

# res by VPD
plot(hh_noNA_noCNHGU$EC_VPD_F_ctr, E1, xlab = "VPD", ylab = "Residuals", main = "Residuals by VPD")
abline(h = 0, lty = 2)

# res by NEE
plot(hh_noNA_noCNHGU$NEE_F_ctr, E1, xlab = "NEE", ylab = "Residuals", main = "Residuals by NEE")
abline(h = 0, lty = 2)

# res by USTAR
plot(hh_noNA_noCNHGU$EC_USTAR_ctr, E1, xlab = "USTAR", ylab = "Residuals", main = "Residuals by USTAR")
abline(h = 0, lty = 2)

# res by u.wind
plot(hh_noNA_noCNHGU$U_WIND_F_ctr, E1, xlab = "U_WIND_F", ylab = "Residuals", main = "Residuals by u.wind")
abline(h = 0, lty = 2)

# res by v.wind
plot(hh_noNA_noCNHGU$V_WIND_F_ctr, E1, xlab = "V_WIND_F", ylab = "Residuals", main = "Residuals by v.wind")
abline(h = 0, lty = 2)

# res by PA
plot(hh_noNA_noCNHGU$EC_PA_F_ctr, E1, xlab = "PA", ylab = "Residuals", main = "Residuals by PA")
abline(h = 0, lty = 2)

# res by DOMINANT_VEGETATION
plot(as.factor(hh_noNA_noCNHGU$DOMINANT_VEGETATION), E1, xlab = "DOMINANT_VEGETATION", ylab = "Residuals", main = "Residuals by DOMINANT_VEGETATION")
abline(h = 0, lty = 2)

# Residuals by MONTH
boxplot(E1 ~ hh_noNA_noCNHGU$MONTH, xlab = "MONTH", ylab = "Residuals", main = "Residuals by MONTH")
abline(h = 0, lty = 2)

# Residuals by DATE
boxplot(E1 ~ hh_noNA_noCNHGU$DATE, xlab = "DATE", ylab = "Residuals")
abline(h = 0, lty = 2)

# Reset plotting layout
par(mfrow = c(1, 1))


##### check temporal autocorrelation 

E1 <- residuals(m_auto, type = "pearson")
E1 <- residuals(m_hh_AR1, type = "pearson")
E1 <- residuals(m_hh_exp3_datehour, type = "pearson")
E1 <- residuals(m_hh_exp3_date, type = "pearson")

# check ACF plot
I1 <- !is.na(hh_noNA_noCNHGU$DELTA_FCH4_abs)
Efull <- vector(length = length(hh_noNA_noCNHGU$DELTA_FCH4_abs))
Efull <- NA
Efull[I1] <- E1
acf(Efull, na.action = na.pass,
    main = "Auto correlation plot for residuals", lag.max = 48)


# corAR1

m_hh_AR1 <- lme(
  DELTA_FCH4_abs_yjtrans_ctr ~ SITE_TS_TOP_ctr +
    SITE_WTL_ctr +
    EC_USTAR_ctr +
    U_WIND_F_ctr +
    V_WIND_F_ctr +
    EC_PA_F_ctr +
    NEE_F_ctr +
    EC_VPD_F_ctr +
    DOMINANT_VEGETATION +
    HOUR,
  random = ~1 | SITE/DateHour,
  correlation = corAR1(form = ~ HalfHour2 | SITE/DateHour),
  data = hh_noNA_noCNHGU
)

AIC(m_auto, m_hh_AR1)


# no corAR1, just site and DateHour as nested random effects
m_hh_sitedatehourre <- lme(
  DELTA_FCH4_abs_yjtrans ~ SITE_TS_TOP_ctr +
    SITE_WTL_ctr +
    EC_USTAR_ctr +
    U_WIND_F_ctr +
    V_WIND_F_ctr +
    EC_PA_F_ctr +
    NEE_F_ctr +
    EC_VPD_F_ctr +
    DOMINANT_VEGETATION +
    HOUR,
  random = ~1 | SITE/DateHour,
  data = hh_noNA_noCNHGU
)

AIC(m_auto, m_hh_AR1, m_hh_sitedatehourre)


##### ----> half-hourly model: 

m_hh_exp3_date <- lme(
  DELTA_FCH4_abs_yjtrans ~ SITE_TS_TOP_ctr +
    SITE_WTL_ctr +
    EC_USTAR_ctr +
    U_WIND_F_ctr +
    V_WIND_F_ctr +
    EC_PA_F_ctr +
    NEE_F_ctr +
    EC_VPD_F_ctr +
    DOMINANT_VEGETATION + 
    HOUR, 
  random = ~1 | SITE/DATE, 
  weights = varComb(
    varExp(form = ~ EC_VPD_F_ctr),
    varExp(form = ~ EC_USTAR_ctr),       
    varExp(form = ~ EC_PA_F_ctr)        
  ), 
  data = hh_noNA_noCNHGU
)


#### daily model

# Simple model without temporal autocorrelation or residual het gen structures
m_d <- lme(
  DELTA_FCH4_abs_yjtrans_ctr ~ SITE_TS_TOP_ctr +
    SITE_WTL_ctr +
    USTAR_MEDIAN_ctr +
    U_WIND_F_MEAN_ctr +
    V_WIND_F_MEAN_ctr +
    PA_F_MEDIAN_ctr +
    NEE_F_MEAN_ctr +
    VPD_F_MEDIAN_ctr +
    DOMINANT_VEGETATION +
    MONTH,
  random = ~1 | SITE,
  data = d_noNA_noCNHGU
)

# the following code explores different models with different residual heterogeneity and temporal autocorrelation structures

m_d_exp <- lme(
  DELTA_FCH4_abs_yjtrans ~ SITE_TS_TOP_ctr +
    SITE_WTL_ctr +
    USTAR_MEDIAN_ctr +
    U_WIND_F_MEAN_ctr +
    V_WIND_F_MEAN_ctr +
    PA_F_MEDIAN_ctr +
    NEE_F_mean_ctr +
    VPD_F_MEDIAN_ctr +
    DOMINANT_VEGETATION +
    MONTH,
  random = ~1 | SITE,
  #correlation = corAR1(form = ~ time | group),
  weights = varExp(form = ~PA_F_MEDIAN_ctr),
  data = d_noNA_noCNHGU
)


m_d_expident <- lme(
  DELTA_FCH4_abs_yjtrans ~ SITE_TS_TOP_ctr +
    SITE_WTL_ctr +
    USTAR_MEDIAN_ctr +
    U_WIND_F_MEAN_ctr +
    V_WIND_F_MEAN_ctr +
    PA_F_MEDIAN_ctr +
    NEE_F_mean_ctr +
    VPD_F_MEDIAN_ctr +
    DOMINANT_VEGETATION +
    MONTH,
  random = ~1 | SITE,
  #correlation = corAR1(form = ~ time | group),
  weights = varComb(
    varIdent(form = ~ 1 | MONTH),
    varExp(form = ~ PA_F_MEDIAN_ctr)
  ),
  data = d_noNA_noCNHGU
)


AIC(m_d, m_d_exp, m_d_novar, m_d_expident)

#### temporal autocorrelation

# first use the model without the variance structures
E1 <- residuals(m_d, type = "pearson")

I1 <- !is.na(d_noNA_noCNHGU$DELTA_FCH4_abs)
Efull <- vector(length = length(d_noNA_noCNHGU$DELTA_FCH4_abs))
Efull <- NA
Efull[I1] <- E1
acf(Efull, na.action = na.pass,
    main = "Auto-correlation plot for residuals")
# this shows that the days are correlated with each other

# corAR1

d_noNA_noCNHGU$TIMESTAMP <- as_date(d_noNA_noCNHGU$TIMESTAMP)

m_d_corAR1 <- lme(
  DELTA_FCH4_abs_yjtrans ~ SITE_TS_TOP_ctr +
    SITE_WTL_ctr +
    USTAR_MEDIAN_ctr +
    U_WIND_F_MEAN_ctr +
    V_WIND_F_MEAN_ctr +
    PA_F_MEDIAN_ctr +
    NEE_F_mean_ctr +
    VPD_F_MEDIAN_ctr +
    DOMINANT_VEGETATION +
    MONTH,
  random = ~1 | SITE,
  correlation = corAR1(form = ~ TIMESTAMP | SITE), # SITE ensures that the correlation is different for all sites- correlation is not the same for all sites
  # weights = varComb(
  #   varIdent(form = ~ 1 | MONTH),
  #   varExp(form = ~ PA_F_MEDIAN_ctr)
  #),
  data = d_noNA_noCNHGU
)

# ACF plots for different models
E1 <- residuals(m_d_corAR1, type = "pearson")
E1 <- residuals(m_d_test, type = "pearson")
E1 <- residuals(m_d_corAR1day, type = "pearson")

I1 <- !is.na(d_noNA_noCNHGU$DELTA_FCH4_abs)
Efull <- vector(length = length(d_noNA_noCNHGU$DELTA_FCH4_abs))
Efull <- NA
Efull[I1] <- E1
acf(Efull, na.action = na.pass,
    main = "Auto-correlation plot for residuals")

AIC(m_d, m_d_corAR1day, m_d_test)

# testing different temporal autocorrelation structures:

d_noNA_noCNHGU$DAY <- as.numeric(d_noNA_noCNHGU$DAY)

d_noNA_noCNHGU$YearMonth <- format(d_noNA_noCNHGU$TIMESTAMP, "%Y-%m")


m_d_corAR1day <- lme(
  DELTA_FCH4_abs_yjtrans_ctr ~ SITE_TS_TOP_ctr +
    SITE_WTL_ctr +
    USTAR_MEDIAN_ctr +
    U_WIND_F_MEAN_ctr +
    V_WIND_F_MEAN_ctr +
    PA_F_MEDIAN_ctr +
    NEE_F_MEAN_ctr +
    VPD_F_MEDIAN_ctr +
    DOMINANT_VEGETATION +
    MONTH,
  random = ~1 | SITE/YearMonth,
  correlation = corAR1(form = ~ DAY | SITE/YearMonth), # SITE ensures that the correlation is different for all sites the correlation is not the same for all sites
  data = d_noNA_noCNHGU
)

# Checking residual diagnostics

E1 <- residuals(m_d_test, type = "pearson")
F1 <- fitted(m_d_test)

E1 <- residuals(m_d, type = "pearson")
F1 <- fitted(m_d)

E1 <- residuals(m_d_corAR1day, type = "pearson")
F1 <- fitted(m_d_corAR1day)

E1 <- residuals(m_d_corAR1day_exp, type = "pearson")
F1 <- fitted(m_d_corAR1day_exp)

E1 <- residuals(m_d_corAR1day_exp2, type = "pearson")
F1 <- fitted(m_d_corAR1day_exp2)

# Diagnostic plots
par(mfrow = c(2, 2), mar = c(5, 5, 2, 2))

# Residuals vs Fitted
plot(F1, E1, xlab = "Fitted values", ylab = "Residuals", main = "Residuals vs Fitted")
abline(h = 0, lty = 2)

# Residuals by SITE
plot(as.factor(d_noNA_noCNHGU$SITE), E1, xlab = "SITE", ylab = "Residuals", main = "Residuals by SITE")
abline(h = 0, lty = 2)

# Residuals by SITE_WTL
plot(d_noNA_noCNHGU$SITE_WTL_ctr, E1, xlab = "SITE_WTL", ylab = "Residuals", main = "Residuals by SITE_WTL")
abline(h = 0, lty = 2)

# res by site ts top
plot(d_noNA_noCNHGU$SITE_TS_TOP_ctr, E1, xlab = "SITE_TS_TOP", ylab = "Residuals", main = "Residuals by SITE_TS_TOP")
abline(h = 0, lty = 2)
# fan shape

# res by VPD
plot(d_noNA_noCNHGU$VPD_F_MEDIAN_ctr, E1, xlab = "VPD", ylab = "Residuals", main = "Residuals by VPD")
abline(h = 0, lty = 2)

# res by NEE
plot(d_noNA_noCNHGU$NEE_F_MEAN_ctr, E1, xlab = "NEE", ylab = "Residuals", main = "Residuals by NEE")
abline(h = 0, lty = 2)

# res by USTAR
plot(d_noNA_noCNHGU$USTAR_MEDIAN_ctr, E1, xlab = "USTAR", ylab = "Residuals", main = "Residuals by USTAR")
abline(h = 0, lty = 2)

# res by u.wind
plot(d_noNA_noCNHGU$U_WIND_F_MEAN_ctr, E1, xlab = "U_WIND_F", ylab = "Residuals", main = "Residuals by u.wind")
abline(h = 0, lty = 2)

# res by v.wind
plot(d_noNA_noCNHGU$V_WIND_F_MEAN_ctr, E1, xlab = "V_WIND_F", ylab = "Residuals", main = "Residuals by v.wind")
abline(h = 0, lty = 2)

# res by PA
plot(d_noNA_noCNHGU$PA_F_MEDIAN_ctr, E1, xlab = "PA", ylab = "Residuals", main = "Residuals by PA")
abline(h = 0, lty = 2)

# DOMINANT_VEGETATION
plot(as.factor(d_noNA_noCNHGU$DOMINANT_VEGETATION), E1, xlab = "DOMINANT_VEGETATION", ylab = "Residuals", main = "Residuals by DOMINANT_VEGETATION")
abline(h = 0, lty = 2)

# Residuals by MONTH
boxplot(E1 ~ d_noNA_noCNHGU$MONTH, xlab = "MONTH", ylab = "Residuals", main = "Residuals by MONTH")
abline(h = 0, lty = 2)

# Residuals by date
boxplot(E1 ~ d_noNA_noCNHGU$TIMESTAMP, xlab = "MONTH", ylab = "Residuals")
abline(h = 0, lty = 2)

# Residuals by year
boxplot(E1 ~ d_noNA_noCNHGU$YEAR, xlab = "YEAR", ylab = "Residuals", main = "Residuals by YEAR")
abline(h = 0, lty = 2)

# Reset plotting layout
par(mfrow = c(1, 1))

# testing more models with temporal autocorrelation + residual heterogeneity structures

m_d_corAR1day_exp <- lme(
  DELTA_FCH4_abs_yjtrans_ctr ~ SITE_TS_TOP_ctr +
    SITE_WTL_ctr +
    USTAR_MEDIAN_ctr +
    U_WIND_F_MEAN_ctr +
    V_WIND_F_MEAN_ctr +
    PA_F_MEDIAN_ctr +
    NEE_F_MEAN_ctr +
    VPD_F_MEDIAN_ctr +
    DOMINANT_VEGETATION + 
    MONTH, 
  random = ~1 | SITE/YearMonth, # yearmonth is nested within SITE
  correlation = corAR1(form = ~ DAY | SITE/YearMonth), # SITE ensures that the correlation is different for all sites, the correlation is not the same for all sites
  weights = varExp(form = ~PA_F_MEDIAN_ctr),
  data = d_noNA_noCNHGU
)


m_d_corAR1day_exp2 <- lme(
  DELTA_FCH4_abs_yjtrans ~ SITE_TS_TOP_ctr +
    SITE_WTL_ctr +
    USTAR_MEDIAN_ctr +
    U_WIND_F_MEAN_ctr +
    V_WIND_F_MEAN_ctr +
    PA_F_MEDIAN_ctr +
    NEE_F_MEAN_ctr +
    VPD_F_MEDIAN_ctr +
    DOMINANT_VEGETATION + 
    MONTH, 
  random = ~1 | SITE/YearMonth, # yearmonth is nested within SITE
  correlation = corAR1(form = ~ DAY | SITE/YearMonth), # SITE ensures that the correlation is different for all sites, correlation is not the same for all sites
  weights = varComb(varExp(form = ~PA_F_MEDIAN_ctr),
                    varExp(form = ~SITE_TS_TOP_ctr)),
  data = d_noNA_noCNHGU
)


# m_d_corAR1day_noexp <- lme(
#   DELTA_FCH4_abs_yjtrans ~ SITE_TS_TOP_ctr +
#     SITE_WTL_ctr +
#     USTAR_MEDIAN_ctr +
#     U_WIND_F_MEAN_ctr +
#     V_WIND_F_MEAN_ctr +
#     PA_F_MEDIAN_ctr +
#     NEE_F_mean_ctr +
#     VPD_F_MEDIAN_ctr +
#     DOMINANT_VEGETATION + 
#     MONTH, 
#   random = ~1 | YearMonth/SITE, # yearmonth is nested within SITE
#   correlation = corAR1(form = ~ DAY | YearMonth/SITE), # SITE ensures that the correlation is different for all sites, we cannot assume that it is the same for all sites
#   #weights = varExp(form = ~PA_F_MEDIAN_ctr),
#   data = d_noNA_noCNHGU
# )

# m_d_corAR1day_varfixed <- lme(
#   DELTA_FCH4_abs_yjtrans ~ SITE_TS_TOP_ctr +
#     SITE_WTL_ctr +
#     USTAR_MEDIAN_ctr +
#     U_WIND_F_MEAN_ctr +
#     V_WIND_F_MEAN_ctr +
#     PA_F_MEDIAN_ctr +
#     NEE_F_mean_ctr +
#     VPD_F_MEDIAN_ctr +
#     DOMINANT_VEGETATION + 
#     MONTH, 
#   random = ~1 | YearMonth/SITE, # yearmonth is nested within SITE
#   correlation = corAR1(form = ~ DAY | YearMonth/SITE), # SITE ensures that the correlation is different for all sites, we cannot assume that it is the same for all sites
#   weights = varFixed(~PA_F_MEDIAN_ctr),
#   data = d_noNA_noCNHGU
# )
# 
# m_d_corAR1day_varpower <- lme(
#   DELTA_FCH4_abs_yjtrans ~ SITE_TS_TOP_ctr +
#     SITE_WTL_ctr +
#     USTAR_MEDIAN_ctr +
#     U_WIND_F_MEAN_ctr +
#     V_WIND_F_MEAN_ctr +
#     PA_F_MEDIAN_ctr +
#     NEE_F_mean_ctr +
#     VPD_F_MEDIAN_ctr +
#     DOMINANT_VEGETATION + 
#     MONTH, 
#   random = ~1 | YearMonth/SITE, # yearmonth is nested within SITE
#   correlation = corAR1(form = ~ DAY | YearMonth/SITE), # SITE ensures that the correlation is different for all sites, we cannot assume that it is the same for all sites
#   weights = varPower(form = ~PA_F_MEDIAN_ctr),
#   data = d_noNA_noCNHGU
# )

AIC(m_d, m_d_corAR1day, m_d_corAR1day_exp,
    m_d_corAR1day_exp2)

AIC(m_d_novar, m_d_corAR1, m_d_corAR1day, m_d_corAR1day_exp, m_d_corAR1day_noexp,
    m_d_corAR1day_varfixed, m_d_corAR1day_varpower)
# corAR1 + varExp for PA seems best
# varfixed is really bad
# exp is best

#varcomb with ident month and exp PA would be best for PA het gen but with corAR1 the model doesn't converge

#### --> d final model: 

m_d_corAR1day_exp2 <- lme(
  DELTA_FCH4_abs_yjtrans_ctr ~ SITE_TS_TOP_ctr +
    SITE_WTL_ctr +
    USTAR_MEDIAN_ctr +
    U_WIND_F_MEAN_ctr +
    V_WIND_F_MEAN_ctr +
    PA_F_MEDIAN_ctr +
    NEE_F_MEAN_ctr +
    VPD_F_MEDIAN_ctr +
    DOMINANT_VEGETATION + 
    MONTH, 
  random = ~1 | SITE/YearMonth, # yearmonth is nested within SITE
  correlation = corAR1(form = ~ DAY | SITE/YearMonth), # SITE ensures that the correlation is different for all sites, the correlation is not the same for all sites
  weights = varComb(varExp(form = ~PA_F_MEDIAN_ctr),
                    varExp(form = ~SITE_TS_TOP_ctr)),
  data = d_noNA_noCNHGU
)

### hourly

# Yeo-Johnson transformation for delta_FCH4
hr_noNA_noCNHGU$DELTA_FCH4_abs <- abs(hr_noNA_noCNHGU$DELTA_FCH4)
yeo_johnson_transformation <- yeojohnson(hr_noNA_noCNHGU$DELTA_FCH4_abs)

# Transform the data
hr_noNA_noCNHGU$DELTA_FCH4_abs_yjtrans <- predict(yeo_johnson_transformation)

hist(hr_noNA_noCNHGU$DELTA_FCH4_abs_yjtrans)
qqPlot(hr_noNA_noCNHGU$DELTA_FCH4_abs_yjtrans)

m_yj <- lmer(DELTA_FCH4_abs_yjtrans ~ SITE_TS_TOP_ctr +
               SITE_WTL_ctr +
               USTAR_MEDIAN_ctr +
               U_WIND_F_MEAN_ctr +
               V_WIND_F_MEAN_ctr +
               PA_F_MEDIAN_ctr +
               NEE_F_MEAN_ctr +
               VPD_F_MEDIAN_ctr +
               DOMINANT_VEGETATION + 
               MONTH +
               HOUR +
               (1 | SITE), data = hr_noNA_noCNHGU)

# model residuals

residuals_yj <- residuals(m_yj, type = "pearson")

hist(residuals_yj, main = "Histogram of YJTransformed Residuals")

qqnorm(residuals_yj, main = "Q-Q Plot for YJ-Transformed Model")
qqline(residuals_yj)

plot(fitted(m_yj), residuals_yj, main = "Residuals vs Fitted (YJ Model)")
abline(h = 0, col = "red")

AIC(m_yj)

# center and scale
hr_noNA_noCNHGU$DELTA_FCH4_abs_yjtrans_ctr <- scale(hr_noNA_noCNHGU$DELTA_FCH4_abs_yjtrans)


##### choose variance and correlation structure

hr_noNA_noCNHGU$MONTH <- as.factor(hr_noNA_noCNHGU$MONTH)
hr_noNA_noCNHGU$HOUR <- as.factor(hr_noNA_noCNHGU$HOUR)

# simple model without variance and correlation structures
m_hr <- lme(
  DELTA_FCH4_abs_yjtrans_ctr ~ SITE_TS_TOP_ctr +
    SITE_WTL_ctr +
    USTAR_MEDIAN_ctr +
    U_WIND_F_MEAN_ctr +
    V_WIND_F_MEAN_ctr +
    PA_F_MEDIAN_ctr +
    NEE_F_MEAN_ctr +
    VPD_F_MEDIAN_ctr +
    DOMINANT_VEGETATION +
    MONTH +
    HOUR,
  random = ~1 | SITE,
  data = hr_noNA_noCNHGU
)

vif(m_hr)

# check residuals, and start from temporal autocorrelation etc

E1 <- residuals(m_hr, type = "pearson")
F1 <- fitted(m_hr)

# E1 <- residuals(m_hr_AR1, type = "pearson")
# F1 <- fitted(m_hr_AR1)
# 
# E1 <- residuals(m_hr_AR1_exp, type = "pearson")
# F1 <- fitted(m_hr_AR1_exp)
# 
# E1 <- residuals(m_hr_AR1_exp2, type = "pearson")
# F1 <- fitted(m_hr_AR1_exp2)
# 
# E1 <- residuals(m_hr_AR1_exp3, type = "pearson")
# F1 <- fitted(m_hr_AR1_exp3)

E1 <- residuals(m_hr_exp3_sitedate, type = "pearson")
F1 <- fitted(m_hr_exp3_sitedate)

# Diagnostic plots
par(mfrow = c(2, 2), mar = c(5, 5, 2, 2))

# Residuals vs Fitted
plot(F1, E1, xlab = "Fitted values", ylab = "Residuals", main = "Residuals vs Fitted")
abline(h = 0, lty = 2)

# Residuals by SITE
plot(as.factor(hr_noNA_noCNHGU$SITE), E1, xlab = "SITE", ylab = "Residuals", main = "Residuals by SITE")
abline(h = 0, lty = 2)

# Residuals by SITE_WTL
plot(hr_noNA_noCNHGU$SITE_WTL_ctr, E1, xlab = "SITE_WTL", ylab = "Residuals", main = "Residuals by SITE_WTL")
abline(h = 0, lty = 2)

# res by site ts top
plot(hr_noNA_noCNHGU$SITE_TS_TOP_ctr, E1, xlab = "SITE_TS_TOP", ylab = "Residuals", main = "Residuals by SITE_TS_TOP")
abline(h = 0, lty = 2)

# res by VPD
plot(hr_noNA_noCNHGU$VPD_F_MEDIAN_ctr, E1, xlab = "VPD", ylab = "Residuals", main = "Residuals by VPD")
abline(h = 0, lty = 2)
# looks a bit worrying, higher variance in smaller VPD


# res by NEE
plot(hr_noNA_noCNHGU$NEE_F_MEAN_ctr, E1, xlab = "NEE", ylab = "Residuals", main = "Residuals by NEE")
abline(h = 0, lty = 2)

# res by USTAR
plot(hr_noNA_noCNHGU$USTAR_MEDIAN_ctr, E1, xlab = "USTAR", ylab = "Residuals", main = "Residuals by USTAR")
abline(h = 0, lty = 2)

# res by u.wind
plot(hr_noNA_noCNHGU$U_WIND_F_MEAN_ctr, E1, xlab = "U_WIND_F", ylab = "Residuals", main = "Residuals by u.wind")
abline(h = 0, lty = 2)

# res by v.wind
plot(hr_noNA_noCNHGU$V_WIND_F_MEAN_ctr, E1, xlab = "V_WIND_F", ylab = "Residuals", main = "Residuals by v.wind")
abline(h = 0, lty = 2)

# res by PA
plot(hr_noNA_noCNHGU$PA_F_MEDIAN_ctr, E1, xlab = "PA", ylab = "Residuals", main = "Residuals by PA")
abline(h = 0, lty = 2)

plot(as.factor(hr_noNA_noCNHGU$DOMINANT_VEGETATION), E1, xlab = "DOMINANT_VEGETATION", ylab = "Residuals", main = "Residuals by DOMINANT_VEGETATION")
abline(h = 0, lty = 2)

# Residuals by MONTH
boxplot(E1 ~ hr_noNA_noCNHGU$MONTH, xlab = "MONTH", ylab = "Residuals", main = "Residuals by MONTH")
abline(h = 0, lty = 2)

# Residuals by DATE
boxplot(E1 ~ hr_noNA_noCNHGU$DATE, xlab = "MONTH", ylab = "Residuals")
abline(h = 0, lty = 2)

# Residuals by YEAR
boxplot(E1 ~ hr_noNA_noCNHGU$YEAR, xlab = "YEAR", ylab = "Residuals", main = "Residuals by YEAR")
abline(h = 0, lty = 2)

# Reset plotting layout
par(mfrow = c(1, 1))

### temporal autocorrelation structures

E1 <- residuals(m_hr, type = "pearson")
# E1 <- residuals(m_hr_AR1, type = "pearson")
# E1 <- residuals(m_hr_AR1_exp, type = "pearson")
# E1 <- residuals(m_hr_AR1_exp2, type = "pearson")
# E1 <- residuals(m_hr_AR1_exp3, type = "pearson")
E1 <- residuals(m_hr_exp3_sitedate, type = "pearson")


I1 <- !is.na(hr_noNA_noCNHGU$DELTA_FCH4_abs_yjtrans)
Efull <- vector(length = length(hr_noNA_noCNHGU$DELTA_FCH4_abs_yjtrans))
Efull <- NA
Efull[I1] <- E1
acf(Efull, na.action = na.pass,
    main = "Auto-correlation plot for residuals")

hr_noNA_noCNHGU$TIMESTAMP <- as_datetime(hr_noNA_noCNHGU$TIMESTAMP)
hr_noNA_noCNHGU$YearMonth <- format(hr_noNA_noCNHGU$TIMESTAMP, "%Y-%m")
hr_noNA_noCNHGU$DATE <- date(hr_noNA_noCNHGU$TIMESTAMP)


# m_hr_AR1 <- lme(
#   DELTA_FCH4_abs_yjtrans ~ SITE_TS_TOP_ctr +
#     SITE_WTL_ctr +
#     USTAR_MEDIAN_ctr +
#     U_WIND_F_MEAN_ctr +
#     V_WIND_F_MEAN_ctr +
#     PA_F_MEDIAN_ctr +
#     NEE_F_MEAN_ctr +
#     VPD_F_MEDIAN_ctr +
#     DOMINANT_VEGETATION +
#     MONTH +
#     HOUR,
#   random = ~1 | SITE/DATE,
#   correlation = corAR1(form = ~ HOUR | SITE/DATE), # SITE ensures that the correlation is different for all sites, we cannot assume that it is the same for all sites
#   #weights = varExp(form = ~PA_F_MEDIAN_ctr),
#   data = hr_noNA_noCNHGU
# )
# # 
# AIC(m_hr, m_hr_AR1)
# # AR1 is much better
# # ACF also improved
# 
# 
# m_hr_AR1_exp <- lme(
#   DELTA_FCH4_abs_yjtrans ~ SITE_TS_TOP_ctr +
#     SITE_WTL_ctr +
#     USTAR_MEDIAN_ctr +
#     U_WIND_F_MEAN_ctr +
#     V_WIND_F_MEAN_ctr +
#     PA_F_MEDIAN_ctr +
#     NEE_F_MEAN_ctr +
#     VPD_F_MEDIAN_ctr +
#     DOMINANT_VEGETATION +
#     MONTH +
#     HOUR,
#   random = ~1 | SITE/DATE,
#   correlation = corAR1(form = ~ HOUR | SITE/DATE), # SITE ensures that the correlation is different for all sites, we cannot assume that it is the same for all sites
#   weights = varExp(form = ~PA_F_MEDIAN_ctr),
#   data = hr_noNA_noCNHGU
# )
# 
# AIC(m_hr, m_hr_AR1, m_hr_AR1_exp)
# # exp is best
# 
# m_hr_AR1_exp2 <- lme(
#   DELTA_FCH4_abs_yjtrans ~ SITE_TS_TOP_ctr +
#     SITE_WTL_ctr +
#     USTAR_MEDIAN_ctr +
#     U_WIND_F_MEAN_ctr +
#     V_WIND_F_MEAN_ctr +
#     PA_F_MEDIAN_ctr +
#     NEE_F_MEAN_ctr +
#     VPD_F_MEDIAN_ctr +
#     DOMINANT_VEGETATION + 
#     MONTH +
#     HOUR, 
#   random = ~1 | SITE/DATE, 
#   correlation = corAR1(form = ~ HOUR | SITE/DATE), # SITE ensures that the correlation is different for all sites, we cannot assume that it is the same for all sites
#   weights = varComb(varExp(form = ~PA_F_MEDIAN_ctr),
#                     varExp(form = ~VPD_F_MEDIAN_ctr)),
#   data = hr_noNA_noCNHGU
# )
# 
# AIC(m_hr, m_hr_AR1, m_hr_AR1_exp, m_hr_AR1_exp2)
# # exp2 is best
# 
# m_hr_AR1_exp3 <- lme(
#   DELTA_FCH4_abs_yjtrans ~ SITE_TS_TOP_ctr +
#     SITE_WTL_ctr +
#     USTAR_MEDIAN_ctr +
#     U_WIND_F_MEAN_ctr +
#     V_WIND_F_MEAN_ctr +
#     PA_F_MEDIAN_ctr +
#     NEE_F_MEAN_ctr +
#     VPD_F_MEDIAN_ctr +
#     DOMINANT_VEGETATION +
#     MONTH +
#     HOUR,
#   random = ~1 | SITE/DATE,
#   correlation = corAR1(form = ~ HOUR | SITE/DATE), # SITE ensures that the correlation is different for all sites, we cannot assume that it is the same for all sites
#   weights = varComb(varExp(form = ~PA_F_MEDIAN_ctr),
#                     varExp(form = ~VPD_F_MEDIAN_ctr),
#                     varExp(form = ~USTAR_MEDIAN_ctr)),
#   data = hr_noNA_noCNHGU
# )
# 
# AIC(m_hr, m_hr_AR1, m_hr_AR1_exp, m_hr_AR1_exp2, m_hr_AR1_exp3)

### check if the AR1 leads to non-convergence and include only nested random effect

hr_noNA_noCNHGU$DATE <- as.factor(hr_noNA_noCNHGU$DATE)
hr_noNA_noCNHGU$SITE <- as.factor(hr_noNA_noCNHGU$SITE)

m_hr_exp3_sitedate <- lme(
  DELTA_FCH4_abs_yjtrans ~ SITE_TS_TOP_ctr +
    SITE_WTL_ctr +
    USTAR_MEDIAN_ctr +
    U_WIND_F_MEAN_ctr +
    V_WIND_F_MEAN_ctr +
    PA_F_MEDIAN_ctr +
    NEE_F_MEAN_ctr +
    VPD_F_MEDIAN_ctr +
    DOMINANT_VEGETATION +
    MONTH +
    HOUR,
  random = ~1 | SITE/DATE,
  #correlation = corAR1(form = ~ HOUR | SITE/DATE), # SITE ensures that the correlation is different for all sites, we cannot assume that it is the same for all sites
  weights = varComb(varExp(form = ~PA_F_MEDIAN_ctr),
                    varExp(form = ~VPD_F_MEDIAN_ctr),
                    varExp(form = ~USTAR_MEDIAN_ctr)),
  data = hr_noNA_noCNHGU
)

AIC(m_hr, m_hr_exp3_sitedate)


# m_hr_exp3_sitedate_noVPD <- lme(
#   DELTA_FCH4_abs_yjtrans ~ SITE_TS_TOP_ctr +
#     SITE_WTL_ctr +
#     USTAR_MEDIAN_ctr +
#     U_WIND_F_MEAN_ctr +
#     V_WIND_F_MEAN_ctr +
#     PA_F_MEDIAN_ctr +
#     NEE_F_MEAN_ctr +
#     VPD_F_MEDIAN_ctr +
#     DOMINANT_VEGETATION +
#     MONTH +
#     HOUR,
#   random = ~1 | SITE/DATE,
#   #correlation = corAR1(form = ~ HOUR | SITE/DATE), # SITE ensures that the correlation is different for all sites, we cannot assume that it is the same for all sites
#   weights = varComb(varExp(form = ~PA_F_MEDIAN_ctr),
#                     #varExp(form = ~VPD_F_MEDIAN_ctr),
#                     varExp(form = ~USTAR_MEDIAN_ctr)),
#   data = hr_noNA_noCNHGU
# )
# 
# 
# m_hr_exp3_sitedate_noUSTAR <- lme(
#   DELTA_FCH4_abs_yjtrans ~ SITE_TS_TOP_ctr +
#     SITE_WTL_ctr +
#     USTAR_MEDIAN_ctr +
#     U_WIND_F_MEAN_ctr +
#     V_WIND_F_MEAN_ctr +
#     PA_F_MEDIAN_ctr +
#     NEE_F_MEAN_ctr +
#     VPD_F_MEDIAN_ctr +
#     DOMINANT_VEGETATION +
#     MONTH +
#     HOUR,
#   random = ~1 | SITE/DATE,
#   #correlation = corAR1(form = ~ HOUR | SITE/DATE), # SITE ensures that the correlation is different for all sites, we cannot assume that it is the same for all sites
#   weights = varComb(varExp(form = ~PA_F_MEDIAN_ctr),
#                     varExp(form = ~VPD_F_MEDIAN_ctr)),
#                     #varExp(form = ~USTAR_MEDIAN_ctr)),
#   data = hr_noNA_noCNHGU
# )
# 
# m_hr_exp3_sitedate_noPA <- lme(
#   DELTA_FCH4_abs_yjtrans ~ SITE_TS_TOP_ctr +
#     SITE_WTL_ctr +
#     USTAR_MEDIAN_ctr +
#     U_WIND_F_MEAN_ctr +
#     V_WIND_F_MEAN_ctr +
#     PA_F_MEDIAN_ctr +
#     NEE_F_MEAN_ctr +
#     VPD_F_MEDIAN_ctr +
#     DOMINANT_VEGETATION +
#     MONTH +
#     HOUR,
#   random = ~1 | SITE/DATE,
#   #correlation = corAR1(form = ~ HOUR | SITE/DATE), # SITE ensures that the correlation is different for all sites, we cannot assume that it is the same for all sites
#   weights = varComb(#varExp(form = ~PA_F_MEDIAN_ctr),
#                     varExp(form = ~VPD_F_MEDIAN_ctr),
#                     varExp(form = ~USTAR_MEDIAN_ctr)),
#   data = hr_noNA_noCNHGU
# )
# 
# AIC(m_hr, m_hr_AR1, m_hr_AR1_exp, m_hr_AR1_exp2, m_hr_AR1_exp3, m_hr_exp3_sitedate,
#     m_hr_exp3_sitedate_noUSTAR, m_hr_exp3_sitedate_noVPD, m_hr_exp3_sitedate_noPA)
# 


# m_hr_AR1_exp2_ustarvpd <- lme(
#   DELTA_FCH4_abs_yjtrans ~ SITE_TS_TOP_ctr +
#     SITE_WTL_ctr +
#     USTAR_MEDIAN_ctr +
#     U_WIND_F_MEAN_ctr +
#     V_WIND_F_MEAN_ctr +
#     PA_F_MEDIAN_ctr +
#     NEE_F_MEAN_ctr +
#     VPD_F_MEDIAN_ctr +
#     DOMINANT_VEGETATION + 
#     MONTH +
#     HOUR, 
#   random = ~1 | DATE/SITE, 
#   correlation = corAR1(form = ~ HOUR | DATE/SITE), # SITE ensures that the correlation is different for all sites, we cannot assume that it is the same for all sites
#   weights = varComb(#varExp(form = ~PA_F_MEDIAN_ctr),
#                     varExp(form = ~VPD_F_MEDIAN_ctr),
#                     varExp(form = ~USTAR_MEDIAN_ctr)),
#   data = hr_noNA_noCNHGU
# )
# 
# AIC(m_hr, m_hr_AR1_exp, m_hr_AR1_exp2, m_hr_AR1_exp2_ustarvpd)
# ustarvpd is worse than just AR1 exp


# m_hr_AR1_exp2_ustarpa <- lme(
#   DELTA_FCH4_abs_yjtrans ~ SITE_TS_TOP_ctr +
#     SITE_WTL_ctr +
#     USTAR_MEDIAN_ctr +
#     U_WIND_F_MEAN_ctr +
#     V_WIND_F_MEAN_ctr +
#     PA_F_MEDIAN_ctr +
#     NEE_F_MEAN_ctr +
#     VPD_F_MEDIAN_ctr +
#     DOMINANT_VEGETATION + 
#     MONTH +
#     HOUR, 
#   random = ~1 | DATE/SITE, 
#   correlation = corAR1(form = ~ HOUR | DATE/SITE), # SITE ensures that the correlation is different for all sites, we cannot assume that it is the same for all sites
#   weights = varComb(varExp(form = ~PA_F_MEDIAN_ctr),
#     #varExp(form = ~VPD_F_MEDIAN_ctr),
#     varExp(form = ~USTAR_MEDIAN_ctr)),
#   data = hr_noNA_noCNHGU
# )
# 
# AIC(m_hr, m_hr_AR1_exp, m_hr_AR1_exp2, m_hr_AR1_exp2_ustarvpd, m_hr_AR1_exp2_ustarpa)
# ustarpa is worse than exp2 (pa vpd)


#### --> hourly model: 

m_hr_exp3_sitedate <- lme(
  DELTA_FCH4_abs_yjtrans_ctr ~ SITE_TS_TOP_ctr +
    SITE_WTL_ctr +
    USTAR_MEDIAN_ctr +
    U_WIND_F_MEAN_ctr +
    V_WIND_F_MEAN_ctr +
    PA_F_MEDIAN_ctr +
    NEE_F_MEAN_ctr +
    VPD_F_MEDIAN_ctr +
    DOMINANT_VEGETATION +
    MONTH +
    HOUR,
  random = ~1 | SITE/DATE,
  weights = varComb(varExp(form = ~PA_F_MEDIAN_ctr),
                    varExp(form = ~VPD_F_MEDIAN_ctr),
                    varExp(form = ~USTAR_MEDIAN_ctr)),
  data = hr_noNA_noCNHGU
)

## check multicollinearity for hourly model (with lmer, because GVIF does not work with lme)

preds_hr <- hr_noNA_noCNHGU[, c("SITE_TS_TOP_ctr", 
                                "SITE_WTL_ctr", 
                                "USTAR_MEDIAN_ctr", 
                                "U_WIND_F_MEAN_ctr", 
                                "V_WIND_F_MEAN_ctr",
                                "PA_F_MEDIAN_ctr",  
                                "NEE_F_MEAN_ctr", 
                                "VPD_F_MEDIAN_ctr"
)]


colnames(preds_hr) <- c("TS", 
                        "WTL",
                        "USTAR", 
                        "u.wind",
                        "v.wind",
                        "PA", 
                        "NEE", 
                        "VPD"
)


# Compute correlation matrix
correlation_matrix <- cor(preds_hr, use = "pairwise.complete.obs")

# Round correlation matrix for clarity
correlation_matrix <- round(correlation_matrix, 2)
corrplot(correlation_matrix, method = "color", type = "upper", order = "hclust", 
         tl.col = "black", tl.srt = 45, diag = FALSE, addCoef.col = "black")

# lmer does not support complex variance and correlation structures so they are not included here
m_hr <- lmer(DELTA_FCH4_abs_yjtrans ~ SITE_TS_TOP_ctr +
               SITE_WTL_ctr +
               USTAR_MEDIAN_ctr +
               U_WIND_F_MEAN_ctr +
               V_WIND_F_MEAN_ctr +
               PA_F_MEDIAN_ctr +
               NEE_F_MEAN_ctr +
               VPD_F_MEDIAN_ctr +
               DOMINANT_VEGETATION + 
               MONTH +
               HOUR +
               (1 | SITE), data = hr_noNA_noCNHGU)

# Calculate the GVIF
vif(m_hr)
# all gvifs are <3 (incl. categorical squares)


### weekly model

# correlations between predictors

preds_week <- week_noNA_noCNHGU[, c("SITE_TS_TOP_ctr", 
                                    "SITE_WTL_ctr", 
                                    "USTAR_MEDIAN_ctr", 
                                    "U_WIND_F_MEAN_ctr", 
                                    "V_WIND_F_MEAN_ctr",
                                    "PA_F_MEDIAN_ctr",  
                                    "NEE_F_MEAN_ctr", 
                                    "VPD_F_MEDIAN_ctr"
)]
colnames(preds_week) <- c("TS", 
                          "WTL",
                          "USTAR", 
                          "u.wind",
                          "v.wind",
                          "PA", 
                          "NEE", 
                          "VPD"
)


# Compute correlation matrix
correlation_matrix <- cor(preds_week, use = "pairwise.complete.obs")

# Round correlation matrix for clarity
correlation_matrix <- round(correlation_matrix, 2)

# Visualize correlation matrix using corrplot
corrplot(correlation_matrix, method = "color", type = "upper", order = "hclust", 
         tl.col = "black", tl.srt = 45, diag = FALSE, addCoef.col = "black")

## Yeo-Johnson transform absolute delta_FCH4

week_noNA_noCNHGU$DELTA_FCH4_abs <- abs(week_noNA_noCNHGU$DELTA_FCH4)

yeo_johnson_transformation <- yeojohnson(week_noNA_noCNHGU$DELTA_FCH4_abs)

# Transform the data
week_noNA_noCNHGU$DELTA_FCH4_abs_yjtrans <- predict(yeo_johnson_transformation)

hist(week_noNA_noCNHGU$DELTA_FCH4_abs_yjtrans)
qqPlot(week_noNA_noCNHGU$DELTA_FCH4_abs_yjtrans)

# center and scale
week_noNA_noCNHGU$DELTA_FCH4_abs_yjtrans_ctr <- scale(week_noNA_noCNHGU$DELTA_FCH4_abs_yjtrans)

m_yj <- lme(
  DELTA_FCH4_abs_yjtrans_ctr ~ SITE_TS_TOP_ctr +
    SITE_WTL_ctr +
    USTAR_MEDIAN_ctr +
    U_WIND_F_MEAN_ctr +
    V_WIND_F_MEAN_ctr +
    PA_F_MEDIAN_ctr +
    NEE_F_MEAN_ctr +
    VPD_F_MEDIAN_ctr +
    DOMINANT_VEGETATION + 
    MONTH, 
  random = ~1 | SITE, 
  data = week_noNA_noCNHGU
)

# model residuals
residuals_yj <- residuals(m_yj, type = "pearson")

hist(residuals_yj, main = "Histogram of YJTransformed Residuals")

qqnorm(residuals_yj, main = "Q-Q Plot for YJ-Transformed Model")
qqline(residuals_yj)

plot(fitted(m_yj), residuals_yj, main = "Residuals vs Fitted (YJ Model)")
abline(h = 0, col = "red")

AIC(m_yj)

# separate weekly models to without VPD and without NEE (collinear)

m_week_novpd <- lme(
  DELTA_FCH4_abs_yjtrans_ctr ~ SITE_TS_TOP_ctr +
    SITE_WTL_ctr +
    USTAR_MEDIAN_ctr +
    U_WIND_F_MEAN_ctr +
    V_WIND_F_MEAN_ctr +
    PA_F_MEDIAN_ctr +
    NEE_F_MEAN_ctr +
    #VPD_F_MEDIAN_ctr +
    DOMINANT_VEGETATION + 
    MONTH, 
  random = ~1 | SITE, 
  data = week_noNA_noCNHGU
)

m_week_nonee <- lme(
  DELTA_FCH4_abs_yjtrans_ctr ~ SITE_TS_TOP_ctr +
    SITE_WTL_ctr +
    USTAR_MEDIAN_ctr +
    U_WIND_F_MEAN_ctr +
    V_WIND_F_MEAN_ctr +
    PA_F_MEDIAN_ctr +
    #NEE_F_MEAN_ctr +
    VPD_F_MEDIAN_ctr +
    DOMINANT_VEGETATION + 
    MONTH, 
  random = ~1 | SITE,  
  data = week_noNA_noCNHGU
)

# Calculate the GVIF
vif(m_week_novpd)
vif(m_week_nonee)

#### model residual diagnostics

E1 <- residuals(m_week_novpd, type = "pearson")
F1 <- fitted(m_week_novpd)

E1 <- residuals(m_week_nonee, type = "pearson")
F1 <- fitted(m_week_nonee)

E1 <- residuals(m_week_AR1_expPA, type = "pearson")
F1 <- fitted(m_week_AR1_expPA)

E1 <- residuals(m_week_monthident, type = "pearson")
F1 <- fitted(m_week_monthident)

E1 <- residuals(m_week_yearident_expPA, type = "pearson")
F1 <- fitted(m_week_yearident_expPA)

E1 <- residuals(m_week_AR1_yearident_expPA, type = "pearson")
F1 <- fitted(m_week_AR1_yearident_expPA)

# Diagnostic plots
par(mfrow = c(2, 2), mar = c(5, 5, 2, 2))

# Residuals vs Fitted
plot(F1, E1, xlab = "Fitted values", ylab = "Residuals", main = "Residuals vs Fitted")
abline(h = 0, lty = 2)

# Residuals by SITE
plot(as.factor(week_noNA_noCNHGU$SITE), E1, xlab = "SITE", ylab = "Residuals", main = "Residuals by SITE")
abline(h = 0, lty = 2)

# Residuals by SITE_WTL
plot(week_noNA_noCNHGU$SITE_WTL_ctr, E1, xlab = "SITE_WTL", ylab = "Residuals", main = "Residuals by SITE_WTL")
abline(h = 0, lty = 2)

# res by site ts top
plot(week_noNA_noCNHGU$SITE_TS_TOP_ctr, E1, xlab = "SITE_TS_TOP", ylab = "Residuals", main = "Residuals by SITE_TS_TOP")
abline(h = 0, lty = 2)

# res by VPD
plot(week_noNA_noCNHGU$VPD_F_MEDIAN_ctr, E1, xlab = "VPD", ylab = "Residuals", main = "Residuals by VPD")
abline(h = 0, lty = 2)

# res by NEE
plot(week_noNA_noCNHGU$NEE_F_MEAN_ctr, E1, xlab = "NEE", ylab = "Residuals", main = "Residuals by NEE")
abline(h = 0, lty = 2)

# res by USTAR
plot(week_noNA_noCNHGU$USTAR_MEDIAN_ctr, E1, xlab = "USTAR", ylab = "Residuals", main = "Residuals by USTAR")
abline(h = 0, lty = 2)

# res by u.wind
plot(week_noNA_noCNHGU$U_WIND_F_MEAN_ctr, E1, xlab = "U_WIND_F", ylab = "Residuals", main = "Residuals by u.wind")
abline(h = 0, lty = 2)

# res by v.wind
plot(week_noNA_noCNHGU$V_WIND_F_MEAN_ctr, E1, xlab = "V_WIND_F", ylab = "Residuals", main = "Residuals by v.wind")
abline(h = 0, lty = 2)

# res by PA
plot(week_noNA_noCNHGU$PA_F_MEDIAN_ctr, E1, xlab = "PA", ylab = "Residuals", main = "Residuals by PA")
abline(h = 0, lty = 2)

# DOMINANT_VEGETATION
plot(as.factor(week_noNA_noCNHGU$DOMINANT_VEGETATION), E1, xlab = "DOMINANT_VEGETATION", ylab = "Residuals", main = "Residuals by DOMINANT_VEGETATION")
abline(h = 0, lty = 2)

# Residuals by MONTH
boxplot(E1 ~ week_noNA_noCNHGU$MONTH, xlab = "MONTH", ylab = "Residuals", main = "Residuals by MONTH")
abline(h = 0, lty = 2)

# Residuals by year
boxplot(E1 ~ week_noNA_noCNHGU$YEAR, xlab = "year", ylab = "Residuals", main = "Residuals by year")
abline(h = 0, lty = 2)

# Reset plotting layout
par(mfrow = c(1, 1))

### temporal autocorrelation

E1 <- residuals(m_week_novpd, type = "pearson")
E1 <- residuals(m_week_nonee, type = "pearson")
E1 <- residuals(m_week_AR1_expPA, type = "pearson")
E1 <- residuals(m_week_AR1_yearident_expPA, type = "pearson")

I1 <- !is.na(week_noNA_noCNHGU$DELTA_FCH4_abs_yjtrans)
Efull <- vector(length = length(week_noNA_noCNHGU$DELTA_FCH4_abs_yjtrans))
Efull <- NA
Efull[I1] <- E1
acf(Efull, na.action = na.pass,
    main = "Auto-correlation plot for residuals")

# create a new YearMonth variable for nested random effects and temporal autocorrelation structures
week_noNA_noCNHGU$YearMonth <- paste(week_noNA_noCNHGU$YEAR, week_noNA_noCNHGU$MONTH, sep = "-")

AIC(m_week, m_week_AR1)

# the following code tests different models with different variance and correlation structures
# by commenting out either NEE or VPD (for this version, NEE was left commented out)

m_week_AR1_expPA <- lme(
  DELTA_FCH4_abs_yjtrans ~ SITE_TS_TOP_ctr +
    SITE_WTL_ctr +
    USTAR_MEDIAN_ctr +
    U_WIND_F_MEAN_ctr +
    V_WIND_F_MEAN_ctr +
    PA_F_MEDIAN_ctr +
    #NEE_F_MEAN_ctr +
    VPD_F_MEDIAN_ctr +
    DOMINANT_VEGETATION + 
    MONTH, 
  random = ~1 | SITE/YearMonth, 
  correlation = corAR1(form = ~ Week_of_year | SITE/YearMonth), # SITE ensures that the correlation is different for all sites, correlation is not the same for all sites
  weights = varExp(form = ~PA_F_MEDIAN_ctr),
  data = week_noNA_noCNHGU
)

AIC(m_week, m_week_AR1_expPA)
# expPA is best

m_week_AR1_yearident_expPA <- lme(
  DELTA_FCH4_abs_yjtrans_ctr ~ SITE_TS_TOP_ctr +
    SITE_WTL_ctr +
    USTAR_MEDIAN_ctr +
    U_WIND_F_MEAN_ctr +
    V_WIND_F_MEAN_ctr +
    PA_F_MEDIAN_ctr +
    #NEE_F_MEAN_ctr +
    VPD_F_MEDIAN_ctr +
    DOMINANT_VEGETATION +
    MONTH,
  random = ~1 | SITE/YearMonth,
  correlation = corAR1(form = ~ Week_of_year | SITE/YearMonth), # SITE ensures that the correlation is different for all sites, correlation is not the same for all sites
  weights = varComb(varIdent(form = ~YEAR),
                    varExp(form= ~PA_F_MEDIAN_ctr)),
  data = week_noNA_noCNHGU
)
 
AIC(m_week_nonee, m_week_AR1_expPA , m_week_AR1_yearident_expPA)


week_noNA_noCNHGU$Week_of_year <- as.numeric(week_noNA_noCNHGU$Week_of_year)

m_week_AR1_expPA <- lme(
  DELTA_FCH4_abs_yjtrans_ctr ~ SITE_TS_TOP_ctr +
    SITE_WTL_ctr +
    USTAR_MEDIAN_ctr +
    U_WIND_F_MEAN_ctr +
    V_WIND_F_MEAN_ctr +
    PA_F_MEDIAN_ctr +
    #NEE_F_MEAN_ctr +
    VPD_F_MEDIAN_ctr +
    DOMINANT_VEGETATION + 
    MONTH, 
  random = ~1 | SITE/YearMonth, 
  correlation = corAR1(form = ~ Week_of_year | SITE/YearMonth), # SITE ensures that the correlation is different for all sites, correlation is not the same for all sites
  weights = varExp(form = ~PA_F_MEDIAN_ctr),
  data = week_noNA_noCNHGU
)


m_week_yearident_expWTL <- lme(
  DELTA_FCH4_abs_yjtrans_ctr ~ SITE_TS_TOP_ctr +
    SITE_WTL_ctr +
    USTAR_MEDIAN_ctr +
    U_WIND_F_MEAN_ctr +
    V_WIND_F_MEAN_ctr +
    PA_F_MEDIAN_ctr +
    #NEE_F_mean_ctr +
    VPD_F_MEDIAN_ctr +
    DOMINANT_VEGETATION +
    MONTH,
  random = ~1 | SITE/YearMonth,
  weights = varComb(varIdent(form = ~YEAR),
                    varExp(form= ~SITE_WTL_ctr)),
  data = week_noNA_noCNHGU
)

AIC(m_week_yearident_expPA, m_week_yearident_expWTL, m_week_AR1_yearident_expWTL,
    m_week_AR1_yearident_expPA)


### --> weekly models: 

m_week_nonee <- lme(
  DELTA_FCH4_abs_yjtrans_ctr ~ SITE_TS_TOP_ctr +
    SITE_WTL_ctr +
    USTAR_MEDIAN_ctr +
    U_WIND_F_MEAN_ctr +
    V_WIND_F_MEAN_ctr +
    PA_F_MEDIAN_ctr +
    #NEE_F_MEAN_ctr +
    VPD_F_MEDIAN_ctr +
    DOMINANT_VEGETATION + 
    MONTH, 
  random = ~1 | SITE/YearMonth, 
  correlation = corAR1(form = ~ Week_of_year | SITE/YearMonth), 
  weights = varExp(form = ~PA_F_MEDIAN_ctr),
  data = week_noNA_noCNHGU
)

m_week_novpd <- lme(
  DELTA_FCH4_abs_yjtrans_ctr ~ SITE_TS_TOP_ctr +
    SITE_WTL_ctr +
    USTAR_MEDIAN_ctr +
    U_WIND_F_MEAN_ctr +
    V_WIND_F_MEAN_ctr +
    PA_F_MEDIAN_ctr +
    NEE_F_MEAN_ctr +
    #VPD_F_MEDIAN_ctr +
    DOMINANT_VEGETATION + 
    MONTH, 
  random = ~1 | SITE/YearMonth, 
  correlation = corAR1(form = ~ Week_of_year | SITE/YearMonth), 
  weights = varExp(form = ~PA_F_MEDIAN_ctr),
  data = week_noNA_noCNHGU
)



### monthly


preds_month <- month_noNA_noCNHGU[, c("SITE_TS_TOP_ctr", 
                                      "SITE_WTL_ctr", 
                                      "USTAR_MEDIAN_ctr", 
                                      "U_WIND_F_MEAN_ctr", 
                                      "V_WIND_F_MEAN_ctr",
                                      "PA_F_MEDIAN_ctr",  
                                      "NEE_F_MEAN_ctr", 
                                      "VPD_F_MEDIAN_ctr"
)]
colnames(preds_month) <- c("TS", 
                           "WTL",
                           "USTAR", 
                           "u.wind",
                           "v.wind",
                           "PA", 
                           "NEE", 
                           "VPD"
)


# Compute correlation matrix
correlation_matrix <- cor(preds_month, use = "pairwise.complete.obs")

# Round correlation matrix for clarity
correlation_matrix <- round(correlation_matrix, 2)

# Visualize correlation matrix using corrplot
corrplot(correlation_matrix, method = "color", type = "upper", order = "hclust", 
         tl.col = "black", tl.srt = 45, diag = FALSE, addCoef.col = "black")


# Yeo-Johnson transform delta_FCH4

month_noNA_noCNHGU$DELTA_FCH4_abs <- abs(month_noNA_noCNHGU$DELTA_FCH4)
yeo_johnson_transformation <- yeojohnson(month_noNA_noCNHGU$DELTA_FCH4_abs)

month_noNA_noCNHGU$DELTA_FCH4_abs_yjtrans <- predict(yeo_johnson_transformation)

hist(month_noNA_noCNHGU$DELTA_FCH4_abs_yjtrans)
qqPlot(month_noNA_noCNHGU$DELTA_FCH4_abs_yjtrans)

# scale and center
month_noNA_noCNHGU$DELTA_FCH4_abs_yjtrans_ctr <- scale(month_noNA_noCNHGU$DELTA_FCH4_abs_yjtrans)

m_yj <- lme(
  DELTA_FCH4_abs_yjtrans_ctr ~ SITE_TS_TOP_ctr +
    SITE_WTL_ctr +
    USTAR_MEDIAN_ctr +
    U_WIND_F_MEAN_ctr +
    V_WIND_F_MEAN_ctr +
    PA_F_MEDIAN_ctr +
    NEE_F_MEAN_ctr +
    VPD_F_MEDIAN_ctr +
    DOMINANT_VEGETATION + 
    MONTH, 
  random = ~1 | SITE,  
  data = month_noNA_noCNHGU
)

# model residuals
residuals_yj <- residuals(m_yj, type = "pearson")

hist(residuals_yj, main = "Histogram of YJTransformed Residuals")

qqnorm(residuals_yj, main = "Q-Q Plot for YJ-Transformed Model")
qqline(residuals_yj)

plot(fitted(m_yj), residuals_yj, main = "Residuals vs Fitted (YJ Model)")
abline(h = 0, col = "red")

AIC(m_yj)


# simple model
m_month <- lme(
  DELTA_FCH4_abs_yjtrans_ctr ~ SITE_TS_TOP_ctr +
    SITE_WTL_ctr +
    USTAR_MEDIAN_ctr +
    U_WIND_F_MEAN_ctr +
    V_WIND_F_MEAN_ctr +
    PA_F_MEDIAN_ctr +
    NEE_F_MEAN_ctr +
    DOMINANT_VEGETATION + 
    MONTH, 
  random = ~1 | SITE, 
  data = month_noNA_noCNHGU
)

# Calculate the GVIF
vif(m_month)

# VPD needs to be removed

### the following code compares models with different variance and correlation structures

m_month <- lme(
  DELTA_FCH4_abs_yjtrans_ctr ~ SITE_TS_TOP_ctr +
    SITE_WTL_ctr +
    USTAR_MEDIAN_ctr +
    U_WIND_F_MEAN_ctr +
    V_WIND_F_MEAN_ctr +
    PA_F_MEDIAN_ctr +
    NEE_F_MEAN_ctr +
    DOMINANT_VEGETATION + 
    MONTH, 
  random = ~1 | SITE/YEAR, 
  data = month_noNA_noCNHGU
)


E1 <- residuals(m_month, type = "pearson")
F1 <- fitted(m_month)

E1 <- residuals(m_month_varident_exp_AR1, type = "pearson")
F1 <- fitted(m_month_varident_exp_AR1)

E1 <- residuals(m_month_varident_exp_AR1, type = "pearson")
F1 <- fitted(m_month_varident_exp_AR1)

E1 <- residuals(m_month_varident_exp, type = "pearson")
F1 <- fitted(m_month_varident_exp)

E1 <- residuals(m_month_exp_AR1, type = "pearson")
F1 <- fitted(m_month_exp_AR1)

E1 <- residuals(m_month_exp, type = "pearson")
F1 <- fitted(m_month_exp)

E1 <- residuals(m_month_yearvarident, type = "pearson")
F1 <- fitted(m_month_yearvarident)

E1 <- residuals(m_month_yearvaridentre, type = "pearson")
F1 <- fitted(m_month_yearvaridentre) 

E1 <- residuals(m_month_yearvarident_yearmonthre, type = "pearson")
F1 <- fitted(m_month_yearvarident_yearmonthre) 

E1 <- residuals(m_month_exp2_2, type = "pearson")
F1 <- fitted(m_month_exp2_2) 

E1 <- residuals(m_month_exp2_AR1_identmonth, type = "pearson")
F1 <- fitted(m_month_exp2_AR1_identmonth) 

E1 <- residuals(m_month_exp2_AR1_identyearmonth, type = "pearson")
F1 <- fitted(m_month_exp2_AR1_identyearmonth) 

E1 <- residuals(m_month_season_exp2, type = "pearson")
F1 <- fitted(m_month_season_exp2) 

E1 <- residuals(m_month_exp3_AR1, type = "pearson")
F1 <- fitted(m_month_exp3_AR1) 

E1 <- residuals(m_month_nomonth_exp2_AR1, type = "pearson")
F1 <- fitted(m_month_nomonth_exp2_AR1) 

E1 <- residuals(m_month_nomonth_exp2, type = "pearson")
F1 <- fitted(m_month_nomonth_exp2) 

E1 <- residuals(m_month_nomonth_exp2_identseason, type = "pearson")
F1 <- fitted(m_month_nomonth_exp2_identseason) 

E1 <- residuals(m_month_nomonth_exp2_identyear, type = "pearson")
F1 <- fitted(m_month_nomonth_exp2_identyear) 

E1 <- residuals(m_month_season_exp3, type = "pearson")
F1 <- fitted(m_month_season_exp3) 

E1 <- residuals(m_month_season_exp3_identmonthyear, type = "pearson")
F1 <- fitted(m_month_season_exp3_identmonthyear)

E1 <- residuals(m_month_season_exp2_identmonthyear, type = "pearson")
F1 <- fitted(m_month_season_exp2_identmonthyear) 

E1 <- residuals(m_month_season_exp2, type = "pearson")
F1 <- fitted(m_month_season_exp2) 

# Diagnostic plots
par(mfrow = c(4, 4), mar = c(5, 5, 2, 2))

# Residuals vs Fitted
plot(F1, E1, xlab = "Fitted values", ylab = "Residuals", main = "Residuals vs Fitted")
abline(h = 0, lty = 2)
# het gen, a bit better with AR1 and exp for u wind

# Residuals by SITE
plot(as.factor(month_noNA_noCNHGU$SITE), E1, xlab = "SITE", ylab = "Residuals", main = "Residuals by SITE")
abline(h = 0, lty = 2)
# het gen, improved with AR1 exp2

# Residuals by SITE_WTL
plot(month_noNA_noCNHGU$SITE_WTL_ctr, E1, xlab = "SITE_WTL", ylab = "Residuals", main = "Residuals by SITE_WTL")
abline(h = 0, lty = 2)

# res by site ts top
plot(month_noNA_noCNHGU$SITE_TS_TOP_ctr, E1, xlab = "SITE_TS_TOP", ylab = "Residuals", main = "Residuals by SITE_TS_TOP")
abline(h = 0, lty = 2)

# res by VPD
plot(month_noNA_noCNHGU$VPD_F_MEDIAN_ctr, E1, xlab = "VPD", ylab = "Residuals", main = "Residuals by VPD")
abline(h = 0, lty = 2)

# res by NEE
plot(month_noNA_noCNHGU$NEE_F_MEAN_ctr, E1, xlab = "NEE", ylab = "Residuals", main = "Residuals by NEE")
abline(h = 0, lty = 2)

# res by USTAR
plot(month_noNA_noCNHGU$USTAR_MEDIAN_ctr, E1, xlab = "USTAR", ylab = "Residuals", main = "Residuals by USTAR")
abline(h = 0, lty = 2)

# res by u.wind
plot(month_noNA_noCNHGU$U_WIND_F_MEAN_ctr, E1, xlab = "U_WIND_F", ylab = "Residuals", main = "Residuals by u.wind")
abline(h = 0, lty = 2)

# res by v.wind
plot(month_noNA_noCNHGU$V_WIND_F_MEAN_ctr, E1, xlab = "V_WIND_F", ylab = "Residuals", main = "Residuals by v.wind")
abline(h = 0, lty = 2)

# res by PA
plot(month_noNA_noCNHGU$PA_F_MEDIAN_ctr, E1, xlab = "PA", ylab = "Residuals", main = "Residuals by PA")
abline(h = 0, lty = 2)

# DOMINANT_VEGETATION
plot(as.factor(month_noNA_noCNHGU$DOMINANT_VEGETATION), E1, xlab = "DOMINANT_VEGETATION", ylab = "Residuals", main = "Residuals by DOMINANT_VEGETATION")
abline(h = 0, lty = 2)

# Residuals by MONTH
boxplot(E1 ~ month_noNA_noCNHGU$MONTH, xlab = "MONTH", ylab = "Residuals", main = "Residuals by MONTH")
abline(h = 0, lty = 2)

boxplot(E1 ~ month_noNA_noCNHGU$YEAR, xlab = "year", ylab = "Residuals", main = "Residuals by year")
abline(h = 0, lty = 2)

# Reset plotting layout
par(mfrow = c(1, 1))

# temporal autocorrelation

E1 <- residuals(m_month, type = "pearson")
E1 <- residuals(m_month_varident_exp_AR1, type = "pearson")
E1 <- residuals(m_month_exp2_2, type = "pearson")

# ACF plot
I1 <- !is.na(month_noNA_noCNHGU$DELTA_FCH4_abs_yjtrans)
Efull <- vector(length = length(month_noNA_noCNHGU$DELTA_FCH4_abs_yjtrans))
Efull <- NA
Efull[I1] <- E1
acf(Efull, na.action = na.pass,
    main = "Auto-correlation plot for residuals")

# compare temporal autocorrelation structures
month_noNA_noCNHGU$Month_factor <- as.factor(month_noNA_noCNHGU$MONTH)
month_noNA_noCNHGU$Month_num <- as.numeric(month_noNA_noCNHGU$MONTH)


m_month_exp2_AR1 <- lme(
  DELTA_FCH4_abs_yjtrans_ctr ~ SITE_TS_TOP_ctr +
    SITE_WTL_ctr +
    USTAR_MEDIAN_ctr +
    U_WIND_F_MEAN_ctr +
    V_WIND_F_MEAN_ctr +
    PA_F_MEDIAN_ctr +
    NEE_F_MEAN_ctr +
    DOMINANT_VEGETATION + 
    Month_factor, 
  random = ~1 | SITE/YEAR, 
  correlation = corAR1(form = ~ Month_num | SITE/YEAR), # SITE ensures that the correlation is different for all sites- the correlation is not the same for all sites
  weights = varComb(varExp(form = ~U_WIND_F_MEAN_ctr),
                    varExp(form = ~USTAR_MEDIAN_ctr)),
  data = month_noNA_noCNHGU
)

m_month_exp2_AR1_2 <- lme(
  DELTA_FCH4_abs_yjtrans_ctr ~ SITE_TS_TOP_ctr +
    SITE_WTL_ctr +
    USTAR_MEDIAN_ctr +
    U_WIND_F_MEAN_ctr +
    V_WIND_F_MEAN_ctr +
    PA_F_MEDIAN_ctr +
    NEE_F_MEAN_ctr +
    DOMINANT_VEGETATION + 
    Month_factor, 
  random = ~1 | SITE/YEAR, 
  correlation = corAR1(form = ~ Month_num | SITE/YEAR), # SITE ensures that the correlation is different for all sites- the correlation is not the same for all sites
  weights = varComb(varExp(form = ~U_WIND_F_MEAN_ctr),
                    varExp(form = ~PA_F_MEDIAN_ctr)),
  data = month_noNA_noCNHGU
)

m_month_exp2_2 <- lme(
  DELTA_FCH4_abs_yjtrans_ctr ~ SITE_TS_TOP_ctr +
    SITE_WTL_ctr +
    USTAR_MEDIAN_ctr +
    U_WIND_F_MEAN_ctr +
    V_WIND_F_MEAN_ctr +
    PA_F_MEDIAN_ctr +
    NEE_F_MEAN_ctr +
    DOMINANT_VEGETATION + 
    Month_factor, 
  random = ~1 | SITE, 
  weights = varComb(varExp(form = ~U_WIND_F_MEAN_ctr),
                    varExp(form = ~PA_F_MEDIAN_ctr)),
  data = month_noNA_noCNHGU
)


AIC(m_month, m_month_exp2_AR1, m_month_exp2_AR1_2, m_month_exp2_2)

# exp2_2 improved AIC from 177 to 170 - it is the best

### --> monthly model:

m_month_exp2_2 <- lme(
  DELTA_FCH4_abs_yjtrans_ctr ~ SITE_TS_TOP_ctr +
    SITE_WTL_ctr +
    USTAR_MEDIAN_ctr +
    U_WIND_F_MEAN_ctr +
    V_WIND_F_MEAN_ctr +
    PA_F_MEDIAN_ctr +
    NEE_F_MEAN_ctr +
    DOMINANT_VEGETATION + 
    MONTH, 
  random = ~1 | SITE, 
  weights = varComb(varExp(form = ~U_WIND_F_MEAN_ctr),
                    varExp(form = ~PA_F_MEDIAN_ctr)),
  data = month_noNA_noCNHGU
)


###### backward stepwise variable selection ######

# for selecting the reference level for month

hr_noNA_noCNHGU %>%
  group_by(MONTH) %>%
  summarize(Sites_with_Data = n_distinct(SITE)) %>%
  arrange(desc(Sites_with_Data))

# 5

hh_noNA_noCNHGU %>%
  group_by(MONTH) %>%
  summarize(Sites_with_Data = n_distinct(SITE)) %>%
  arrange(desc(Sites_with_Data))

# 5

d_noNA_noCNHGU %>%
  group_by(MONTH) %>%
  summarize(Sites_with_Data = n_distinct(SITE)) %>%
  arrange(desc(Sites_with_Data))

# 5

week_noNA_noCNHGU %>%
  group_by(MONTH) %>%
  summarize(Sites_with_Data = n_distinct(SITE)) %>%
  arrange(desc(Sites_with_Data))

# 5

month_noNA_noCNHGU %>%
  group_by(MONTH) %>%
  summarize(Sites_with_Data = n_distinct(SITE)) %>%
  arrange(desc(Sites_with_Data))

# 5

hr_noNA_noCNHGU$MONTH <- as.factor(hr_noNA_noCNHGU$MONTH)
# Change the reference level for MONTH
hr_noNA_noCNHGU$MONTH <- relevel(hr_noNA_noCNHGU$MONTH, ref = "5")
hh_noNA_noCNHGU$MONTH <- relevel(as.factor(hh_noNA_noCNHGU$MONTH), ref = "5")
d_noNA_noCNHGU$MONTH <- relevel(as.factor(d_noNA_noCNHGU$MONTH), ref = "5")
week_noNA_noCNHGU$MONTH <- relevel(as.factor(week_noNA_noCNHGU$MONTH), ref = "5")
month_noNA_noCNHGU$MONTH <- relevel(as.factor(month_noNA_noCNHGU$MONTH), ref = "5")

# change ref level to moss sphagnum
week_noNA_noCNHGU$DOMINANT_VEGETATION <- as.factor(week_noNA_noCNHGU$DOMINANT_VEGETATION)
hr_noNA_noCNHGU$DOMINANT_VEGETATION <- as.factor(hr_noNA_noCNHGU$DOMINANT_VEGETATION)
month_noNA_noCNHGU$DOMINANT_VEGETATION <- as.factor(month_noNA_noCNHGU$DOMINANT_VEGETATION)

hr_noNA_noCNHGU$DOMINANT_VEGETATION <- relevel(hr_noNA_noCNHGU$DOMINANT_VEGETATION, ref = "MOSS_SPHAGNUM")
hh_noNA_noCNHGU$DOMINANT_VEGETATION <- relevel(as.factor(hh_noNA_noCNHGU$DOMINANT_VEGETATION), ref = "MOSS_SPHAGNUM")
d_noNA_noCNHGU$DOMINANT_VEGETATION <- relevel(d_noNA_noCNHGU$DOMINANT_VEGETATION, ref = "MOSS_SPHAGNUM")
week_noNA_noCNHGU$DOMINANT_VEGETATION <- relevel(week_noNA_noCNHGU$DOMINANT_VEGETATION, ref = "MOSS_SPHAGNUM")
month_noNA_noCNHGU$DOMINANT_VEGETATION <- relevel(month_noNA_noCNHGU$DOMINANT_VEGETATION, ref = "MOSS_SPHAGNUM")

##### half-hourly

hh_noNA_noCNHGU$HOUR <- as.factor(hh_noNA_noCNHGU$HOUR)
hh_noNA_noCNHGU$MONTH <- as.factor(hh_noNA_noCNHGU$MONTH)
hh_noNA_noCNHGU$SITE <- as.factor(hh_noNA_noCNHGU$SITE)
hh_noNA_noCNHGU$DATE <- as.factor(hh_noNA_noCNHGU$DATE)

# without TS
m_auto <- lme(
  DELTA_FCH4_abs_yjtrans_ctr ~ #SITE_TS_TOP_ctr +
    SITE_WTL_ctr +
    EC_USTAR_ctr +
    U_WIND_F_ctr +
    V_WIND_F_ctr +
    EC_PA_F_ctr +
    NEE_F_ctr +
    EC_VPD_F_ctr +
    HOUR +
    DOMINANT_VEGETATION +
    MONTH, 
  random = ~1 | SITE/DATE, 
  weights = varComb(
    varExp(form = ~ EC_VPD_F_ctr),
    varExp(form = ~ EC_USTAR_ctr),       
    varExp(form = ~ EC_PA_F_ctr)        
  ), 
  data = hh_noNA_noCNHGU,
  method = "ML"
)

summary(m_auto)
anova(m_auto)

# remove dom veg

m_hh_2 <- lme(
  DELTA_FCH4_abs_yjtrans_ctr ~ #SITE_TS_TOP_ctr +
    SITE_WTL_ctr +
    EC_USTAR_ctr +
    U_WIND_F_ctr +
    V_WIND_F_ctr +
    EC_PA_F_ctr +
    NEE_F_ctr +
    EC_VPD_F_ctr +
    HOUR +
    #DOMINANT_VEGETATION +
    MONTH, 
  random = ~1 | SITE/DATE, 
  weights = varComb(
    varExp(form = ~ EC_VPD_F_ctr),
    varExp(form = ~ EC_USTAR_ctr),       
    varExp(form = ~ EC_PA_F_ctr)        
  ), 
  data = hh_noNA_noCNHGU,
  method = "ML"
)

anova(m_auto, m_hh_2)
# p = 0.2164 and AIC decreased --> remove dom veg

anova(m_hh_2)
# remove uwind

m_hh_3 <- lme(
  DELTA_FCH4_abs_yjtrans_ctr ~ #SITE_TS_TOP_ctr +
    SITE_WTL_ctr +
    EC_USTAR_ctr +
    #U_WIND_F_ctr +
    V_WIND_F_ctr +
    EC_PA_F_ctr +
    NEE_F_ctr +
    EC_VPD_F_ctr +
    HOUR +
    #DOMINANT_VEGETATION +
    MONTH, 
  random = ~1 | SITE/DATE, 
  weights = varComb(
    varExp(form = ~ EC_VPD_F_ctr),
    varExp(form = ~ EC_USTAR_ctr),       
    varExp(form = ~ EC_PA_F_ctr)        
  ), 
  data = hh_noNA_noCNHGU,
  method = "ML"
)

anova(m_hh_2, m_hh_3)
# AIC decreased, p = 0.6706 --> remove uwind

anova(m_hh_3)
# remove wtd

m_hh_4 <- lme(
  DELTA_FCH4_abs_yjtrans_ctr ~ #SITE_TS_TOP_ctr +
    #SITE_WTL_ctr +
    EC_USTAR_ctr +
    #U_WIND_F_ctr +
    V_WIND_F_ctr +
    EC_PA_F_ctr +
    NEE_F_ctr +
    EC_VPD_F_ctr +
    HOUR +
    #DOMINANT_VEGETATION +
    MONTH, 
  random = ~1 | SITE/DATE, 
  weights = varComb(
    varExp(form = ~ EC_VPD_F_ctr),
    varExp(form = ~ EC_USTAR_ctr),       
    varExp(form = ~ EC_PA_F_ctr)        
  ), 
  data = hh_noNA_noCNHGU,
  method = "ML"
)

anova(m_hh_3, m_hh_4)
# AIC increased and p<0.0001 --> keep wtd

### --> final model: 

m_hh_3 <- lme(
  DELTA_FCH4_abs_yjtrans_ctr ~ #SITE_TS_TOP_ctr +
    SITE_WTL_ctr +
    EC_USTAR_ctr +
    #U_WIND_F_ctr +
    V_WIND_F_ctr +
    EC_PA_F_ctr +
    NEE_F_ctr +
    EC_VPD_F_ctr +
    HOUR +
    #DOMINANT_VEGETATION +
    MONTH, 
  random = ~1 | SITE/DATE, 
  weights = varComb(
    varExp(form = ~ EC_VPD_F_ctr),
    varExp(form = ~ EC_USTAR_ctr),       
    varExp(form = ~ EC_PA_F_ctr)        
  ), 
  data = hh_noNA_noCNHGU,
  method = "REML"
)

# model with TS instead of month (for SI)

m_auto <- lme(
  DELTA_FCH4_abs_yjtrans_ctr ~ SITE_TS_TOP_ctr +
    SITE_WTL_ctr +
    EC_USTAR_ctr +
    U_WIND_F_ctr +
    V_WIND_F_ctr +
    EC_PA_F_ctr +
    NEE_F_ctr +
    EC_VPD_F_ctr + 
    HOUR + 
    DOMINANT_VEGETATION, 
  random = ~1 | SITE/DATE, 
  weights = varComb(
    varExp(form = ~ EC_VPD_F_ctr),
    varExp(form = ~ EC_USTAR_ctr),       
    varExp(form = ~ EC_PA_F_ctr)        
  ), 
  data = hh_noNA_noCNHGU,
  method="ML"
)

vif(m_auto) 
summary(m_auto)
anova(m_auto) 
# remove u wind

m_hh_2 <- lme(
  DELTA_FCH4_abs_yjtrans_ctr ~ SITE_TS_TOP_ctr +
    SITE_WTL_ctr +
    EC_USTAR_ctr +
    #U_WIND_F_ctr +
    V_WIND_F_ctr +
    EC_PA_F_ctr +
    NEE_F_ctr +
    EC_VPD_F_ctr + 
    HOUR + 
    DOMINANT_VEGETATION, 
  random = ~1 | SITE/DATE, 
  weights = varComb(
    varExp(form = ~ EC_VPD_F_ctr),
    varExp(form = ~ EC_USTAR_ctr),       
    varExp(form = ~ EC_PA_F_ctr)        
  ), 
  data = hh_noNA_noCNHGU,
  method="ML"
)

anova(m_auto, m_hh_2)
# AIC decreased, p=0.8497


anova(m_hh_2)
# remove WTL

m_hh_3 <- lme(
  DELTA_FCH4_abs_yjtrans_ctr ~ SITE_TS_TOP_ctr +
    #SITE_WTL_ctr +
    EC_USTAR_ctr +
    #U_WIND_F_ctr +
    V_WIND_F_ctr +
    EC_PA_F_ctr +
    NEE_F_ctr +
    EC_VPD_F_ctr + 
    HOUR + 
    DOMINANT_VEGETATION, 
  random = ~1 | SITE/DATE, 
  weights = varComb(
    varExp(form = ~ EC_VPD_F_ctr),
    varExp(form = ~ EC_USTAR_ctr),       
    varExp(form = ~ EC_PA_F_ctr)        
  ), 
  data = hh_noNA_noCNHGU,
  method="ML"
)

anova(m_hh_2, m_hh_3)
# AIC increased a bit, p = 0.1551 --> remove WTL

anova(m_hh_3)
# remove dom veg

m_hh_4 <- lme(
  DELTA_FCH4_abs_yjtrans_ctr ~ SITE_TS_TOP_ctr +
    #SITE_WTL_ctr +
    EC_USTAR_ctr +
    #U_WIND_F_ctr +
    V_WIND_F_ctr +
    EC_PA_F_ctr +
    NEE_F_ctr +
    EC_VPD_F_ctr + 
    HOUR,  
  #DOMINANT_VEGETATION, 
  random = ~1 | SITE/DATE, 
  weights = varComb(
    varExp(form = ~ EC_VPD_F_ctr),
    varExp(form = ~ EC_USTAR_ctr),       
    varExp(form = ~ EC_PA_F_ctr)        
  ), 
  data = hh_noNA_noCNHGU,
  method="ML"
)

anova(m_hh_3, m_hh_4)
# AIC decreased, p = 0.2043 --> remove dom veg


anova(m_hh_4)
# all are significant --> stop

# final model (with TS, no month): 

m_hh_4 <- lme(
  DELTA_FCH4_abs_yjtrans_ctr ~ SITE_TS_TOP_ctr +
    #SITE_WTL_ctr +
    EC_USTAR_ctr +
    #U_WIND_F_ctr +
    V_WIND_F_ctr +
    EC_PA_F_ctr +
    NEE_F_ctr +
    EC_VPD_F_ctr + 
    HOUR,  
  #DOMINANT_VEGETATION, 
  random = ~1 | SITE/DATE, 
  weights = varComb(
    varExp(form = ~ EC_VPD_F_ctr),
    varExp(form = ~ EC_USTAR_ctr),       
    varExp(form = ~ EC_PA_F_ctr)        
  ), 
  data = hh_noNA_noCNHGU,
  method="REML"
)

# refit the full models with REML

# without TS, with month
m_auto <- lme(
  DELTA_FCH4_abs_yjtrans_ctr ~ #SITE_TS_TOP_ctr +
    SITE_WTL_ctr +
    EC_USTAR_ctr +
    U_WIND_F_ctr +
    V_WIND_F_ctr +
    EC_PA_F_ctr +
    NEE_F_ctr +
    EC_VPD_F_ctr +
    DOMINANT_VEGETATION +
    HOUR +
    MONTH,
  random = ~1 | SITE/DATE,
  weights = varComb(
    varExp(form = ~ EC_VPD_F_ctr),
    varExp(form = ~ EC_USTAR_ctr),
    varExp(form = ~ EC_PA_F_ctr)
  ),
  data = hh_noNA_noCNHGU,
  method="REML"
)

# without month, with TS
m_auto <- lme(
  DELTA_FCH4_abs_yjtrans_ctr ~ SITE_TS_TOP_ctr +
    SITE_WTL_ctr +
    EC_USTAR_ctr +
    U_WIND_F_ctr +
    V_WIND_F_ctr +
    EC_PA_F_ctr +
    NEE_F_ctr +
    EC_VPD_F_ctr + 
    HOUR, 
  random = ~1 | SITE/DATE, 
  weights = varComb(
    varExp(form = ~ EC_VPD_F_ctr),
    varExp(form = ~ EC_USTAR_ctr),       
    varExp(form = ~ EC_PA_F_ctr)        
  ), 
  data = hh_noNA_noCNHGU,
  method="REML"
)



##### effect sizes for the fixed effects

# Extract raw fixed-effect coefficients

raw_coefficients <- fixed.effects(m_hh_3) # with month, no TS
raw_coefficients <- fixed.effects(m_hh_4) # with TS, no month

raw_coefficients <- fixed.effects(m_auto) # full model without TS or month

# Sort beta coefficients by absolute value
sorted_beta_absolute <- raw_coefficients[order(abs(raw_coefficients), decreasing = TRUE)]

# View sorted coefficients
sorted_beta_absolute

anova(m_hh_3) # without TS, with month
summary(m_hh_3)

anova(m_auto) # full model with TS or month
summary(m_auto) # full model with TS or month

anova(m_hh_4) # with TS, no month
summary(m_hh_4)

# Extract variance components
var_components <- VarCorr(m_hh_4)
var_components <- VarCorr(m_hh_3)
var_components <- VarCorr(m_auto)
VarCorr(m_auto)
var_components

# Total variance = random var + reisdual var
# calculate manually based on var_components

total_var <- 0.82785795 + 0.07128008 + 0.24449454

total_var <- 0.82983216 + 0.07127437 +  0.24449507

total_var <- 0.99689906 + 0.07128932 + 0.24449828

total_var <- 1.16554210 + 0.09286666 + 0.22623930

# site:
(0.82785795 / total_var) * 100
(0.89538139 / total_var) * 100 
(0.82983216 / total_var) * 100
(0.99689906 / total_var) * 100

# date

(0.07128008 / total_var) * 100
(0.09907175 / total_var) * 100 
(0.07127437 / total_var) * 100
(0.07128932 / total_var) * 100
r.squaredGLMM(m_hh_4)
r.squaredGLMM(m_hh_3)
r.squaredGLMM(m_auto)

#### check the residuals again

E1 <- residuals(m_hh_3, type = "pearson")
F1 <- fitted(m_hh_3) 

# Diagnostic plots
par(mfrow = c(4, 4), mar = c(5, 5, 2, 2))

# Residuals vs Fitted
plot(F1, E1, xlab = "Fitted values", ylab = "Residuals", main = "Residuals vs Fitted")
abline(h = 0, lty = 2)
# het gen, but pretty much the same as the one we started with

# Residuals by SITE
plot(as.factor(hh_noNA_noCNHGU$SITE), E1, xlab = "SITE", ylab = "Residuals", main = "Residuals by SITE")
abline(h = 0, lty = 2)
# okay

# Residuals by SITE_WTL
plot(hh_noNA_noCNHGU$SITE_WTL_ctr, E1, xlab = "SITE_WTL", ylab = "Residuals", main = "Residuals by SITE_WTL")
abline(h = 0, lty = 2)
# some het gen but probs okay

# res by site ts top
plot(hh_noNA_noCNHGU$SITE_TS_TOP_ctr, E1, xlab = "SITE_TS_TOP", ylab = "Residuals", main = "Residuals by SITE_TS_TOP")
abline(h = 0, lty = 2)

# res by site ts top
plot(hh_noNA_noCNHGU$EC_VPD_F_ctr, E1, xlab = "VPD", ylab = "Residuals", main = "Residuals by VPD")
abline(h = 0, lty = 2)

# res by site ts top
plot(hh_noNA_noCNHGU$NEE_F_ctr, E1, xlab = "NEE", ylab = "Residuals", main = "Residuals by SITE_TS_TOP")
abline(h = 0, lty = 2)

# res by USTAR
plot(hh_noNA_noCNHGU$EC_USTAR_ctr, E1, xlab = "USTAR", ylab = "Residuals", main = "Residuals by USTAR")
abline(h = 0, lty = 2)

# res by u.wind
plot(hh_noNA_noCNHGU$U_WIND_F_ctr, E1, xlab = "U_WIND_F", ylab = "Residuals", main = "Residuals by u.wind")
abline(h = 0, lty = 2)

# res by v.wind
plot(hh_noNA_noCNHGU$V_WIND_F_ctr, E1, xlab = "V_WIND_F", ylab = "Residuals", main = "Residuals by v.wind")
abline(h = 0, lty = 2)

# res by PA
plot(hh_noNA_noCNHGU$EC_PA_F_ctr, E1, xlab = "PA", ylab = "Residuals", main = "Residuals by PA")
abline(h = 0, lty = 2)

# DOMINANT_VEGETATION
plot(as.factor(hh_noNA_noCNHGU$DOMINANT_VEGETATION), E1, xlab = "DOMINANT_VEGETATION", ylab = "Residuals", main = "Residuals by DOMINANT_VEGETATION")
abline(h = 0, lty = 2)

# Residuals by MONTH
boxplot(E1 ~ hh_noNA_noCNHGU$MONTH, xlab = "MONTH", ylab = "Residuals", main = "Residuals by MONTH")
abline(h = 0, lty = 2)

# Reset plotting layout
par(mfrow = c(1, 1))

# test normality
residuals <- residuals(m_hh_3)

qqPlot(residuals)
hist(residuals)

#### hourly model

hr_noNA_noCNHGU$HOUR <- as.factor(hr_noNA_noCNHGU$HOUR)
hr_noNA_noCNHGU$TIMESTAMP <- as_datetime(hr_noNA_noCNHGU$TIMESTAMP)

hr_noNA_noCNHGU$DOMINANT_VEGETATION <- as.factor(hr_noNA_noCNHGU$DOMINANT_VEGETATION)
hr_noNA_noCNHGU$MONTH <- as.factor(hr_noNA_noCNHGU$MONTH)
hr_noNA_noCNHGU$Hour_factor <- as.factor(hr_noNA_noCNHGU$HOUR)
hr_noNA_noCNHGU$Hour_hour <- hour(as_datetime(hr_noNA_noCNHGU$TIMESTAMP))
hr_noNA_noCNHGU$DATE <- as_date(hr_noNA_noCNHGU$DATE)
hr_noNA_noCNHGU$SITE <- as.factor(hr_noNA_noCNHGU$SITE)

# full model

m_hr <- lme(
  DELTA_FCH4_abs_yjtrans_ctr ~ SITE_TS_TOP_ctr +
    SITE_WTL_ctr +
    USTAR_MEDIAN_ctr +
    U_WIND_F_MEAN_ctr +
    V_WIND_F_MEAN_ctr +
    PA_F_MEDIAN_ctr +
    NEE_F_mean_ctr +
    VPD_F_MEDIAN_ctr +
    HOUR +
    DOMINANT_VEGETATION +
    MONTH,
  random = ~1 | SITE/DATE,
  weights = varComb(varExp(form = ~PA_F_MEDIAN_ctr),
                    varExp(form = ~VPD_F_MEDIAN_ctr),
                    varExp(form = ~USTAR_MEDIAN_ctr)),
  data = hr_noNA_noCNHGU,
  method = "ML",
  control = lmeControl(opt = "nlminb", maxIter = 100, msMaxIter = 100)
)

vif(m_hr)

# check anova and p values

anova(m_hr)
summary(m_hr)

# remove dom veg

m_hr_2 <- lme(
  DELTA_FCH4_abs_yjtrans_ctr ~ SITE_TS_TOP_ctr +
    SITE_WTL_ctr +
    USTAR_MEDIAN_ctr +
    U_WIND_F_MEAN_ctr +
    V_WIND_F_MEAN_ctr +
    PA_F_MEDIAN_ctr +
    NEE_F_mean_ctr +
    VPD_F_MEDIAN_ctr +
    HOUR +
    #DOMINANT_VEGETATION +
    MONTH,
  random = ~1 | SITE/DATE,
  weights = varComb(varExp(form = ~PA_F_MEDIAN_ctr),
                    varExp(form = ~VPD_F_MEDIAN_ctr),
                    varExp(form = ~USTAR_MEDIAN_ctr)),
  data = hr_noNA_noCNHGU,
  method = "ML"
)
# compare models
anova(m_hr, m_hr_2)
# AIC decreased, p=0.2554 --> remove dom veg

anova(m_hr_2)
# remove u wind

m_hr_3 <- lme(
  DELTA_FCH4_abs_yjtrans_ctr ~ SITE_TS_TOP_ctr +
    SITE_WTL_ctr +
    USTAR_MEDIAN_ctr +
    #U_WIND_F_MEAN_ctr +
    V_WIND_F_MEAN_ctr +
    PA_F_MEDIAN_ctr +
    NEE_F_mean_ctr +
    VPD_F_MEDIAN_ctr +
    HOUR +
    #DOMINANT_VEGETATION +
    MONTH,
  random = ~1 | SITE/DATE,
  weights = varComb(varExp(form = ~PA_F_MEDIAN_ctr),
                    varExp(form = ~VPD_F_MEDIAN_ctr),
                    varExp(form = ~USTAR_MEDIAN_ctr)),
  data = hr_noNA_noCNHGU,
  method = "ML"
)
# compare models
anova(m_hr_2, m_hr_3)
# AIC decreased, p=0.8351 --> remove u wind

anova(m_hr_3)
# all are significant --> stop

# final model: 

m_hr_3 <- lme(
  DELTA_FCH4_abs_yjtrans_ctr ~ SITE_TS_TOP_ctr +
    SITE_WTL_ctr +
    USTAR_MEDIAN_ctr +
    #U_WIND_F_MEAN_ctr +
    V_WIND_F_MEAN_ctr +
    PA_F_MEDIAN_ctr +
    NEE_F_mean_ctr +
    VPD_F_MEDIAN_ctr +
    HOUR +
    #DOMINANT_VEGETATION +
    MONTH,
  random = ~1 | SITE/DATE,
  weights = varComb(varExp(form = ~PA_F_MEDIAN_ctr),
                    varExp(form = ~VPD_F_MEDIAN_ctr),
                    varExp(form = ~USTAR_MEDIAN_ctr)),
  data = hr_noNA_noCNHGU,
  method = "REML"
)

### refit the full model with REML

m_hr <- lme(
  DELTA_FCH4_abs_yjtrans_ctr ~ SITE_TS_TOP_ctr +
    SITE_WTL_ctr +
    USTAR_MEDIAN_ctr +
    U_WIND_F_MEAN_ctr +
    V_WIND_F_MEAN_ctr +
    PA_F_MEDIAN_ctr +
    NEE_F_mean_ctr +
    VPD_F_MEDIAN_ctr +
    HOUR +
    DOMINANT_VEGETATION +
    MONTH,
  random = ~1 | SITE/DATE,
  weights = varComb(varExp(form = ~PA_F_MEDIAN_ctr),
                    varExp(form = ~VPD_F_MEDIAN_ctr),
                    varExp(form = ~USTAR_MEDIAN_ctr)),
  data = hr_noNA_noCNHGU,
  method = "REML",
  control = lmeControl(opt = "nlminb", maxIter = 100, msMaxIter = 100)
)

##### effect sizes for the fixed effects

# Extract raw fixed-effect coefficients
raw_coefficients <- fixed.effects(m_hr_3)
raw_coefficients <- fixed.effects(m_hr)

# Sort beta coefficients by absolute value
sorted_beta_absolute <- raw_coefficients[order(abs(raw_coefficients), decreasing = TRUE)]

# View sorted coefficients
sorted_beta_absolute

# random effect variation explained

# Extract variance components
var_components <- VarCorr(m_hr)
var_components <- VarCorr(m_hr_3)
var_components

# Total variance = random var + residual var
total_var <- 0.86123359 + 0.09855655 + 0.22877931
total_var <- 1.11919297 + 0.09856878 + 0.22878657
total_var <- 0.92001789 + 0.08266157 + 0.20014971

# site:
(0.86123359 / total_var) * 100 # 72.45968
(1.11919297 / total_var) * 100 # 77.3699

# date

(0.09855655 / total_var) * 100 # 8.292031
(0.09856878  / total_var) * 100 # 6.814068

r.squaredGLMM(m_hr_3)
r.squaredGLMM(m_hr)

summary(m_hr_3)
anova(m_hr_3)

summary(m_hr)

#### check the residuals again

E1 <- residuals(m_hr_3, type = "pearson")
F1 <- fitted(m_hr_3) 

# Diagnostic plots
par(mfrow = c(4, 4), mar = c(5, 5, 2, 2))

# Residuals vs Fitted
plot(F1, E1, xlab = "Fitted values", ylab = "Residuals", main = "Residuals vs Fitted")
abline(h = 0, lty = 2)
# het gen, but pretty much the same as the one we started with

# Residuals by SITE
plot(as.factor(hr_noNA_noCNHGU$SITE), E1, xlab = "SITE", ylab = "Residuals", main = "Residuals by SITE")
abline(h = 0, lty = 2)
# okay

# Residuals by SITE_WTL
plot(hr_noNA_noCNHGU$SITE_WTL_ctr, E1, xlab = "SITE_WTL", ylab = "Residuals", main = "Residuals by SITE_WTL")
abline(h = 0, lty = 2)
# some het gen but probs okay

# res by site ts top
plot(hr_noNA_noCNHGU$SITE_TS_TOP_ctr, E1, xlab = "SITE_TS_TOP", ylab = "Residuals", main = "Residuals by SITE_TS_TOP")
abline(h = 0, lty = 2)

# res by VPD
plot(hr_noNA_noCNHGU$VPD_F_MEDIAN_ctr, E1, xlab = "VPD", ylab = "Residuals", main = "Residuals by VPD")
abline(h = 0, lty = 2)

# res by NEE
plot(hr_noNA_noCNHGU$NEE_F_mean_ctr, E1, xlab = "NEE", ylab = "Residuals", main = "Residuals by NEE")
abline(h = 0, lty = 2)

# res by USTAR
plot(hr_noNA_noCNHGU$USTAR_MEDIAN_ctr, E1, xlab = "USTAR", ylab = "Residuals", main = "Residuals by USTAR")
abline(h = 0, lty = 2)

# res by u.wind
plot(hr_noNA_noCNHGU$U_WIND_F_MEAN_ctr, E1, xlab = "U_WIND_F", ylab = "Residuals", main = "Residuals by u.wind")
abline(h = 0, lty = 2)

# res by v.wind
plot(hr_noNA_noCNHGU$V_WIND_F_MEAN_ctr, E1, xlab = "V_WIND_F", ylab = "Residuals", main = "Residuals by v.wind")
abline(h = 0, lty = 2)

# res by PA
plot(hr_noNA_noCNHGU$PA_F_MEDIAN_ctr, E1, xlab = "PA", ylab = "Residuals", main = "Residuals by PA")
abline(h = 0, lty = 2)

# res by DOMINANT_VEGETATION
plot(as.factor(hr_noNA_noCNHGU$DOMINANT_VEGETATION), E1, xlab = "DOMINANT_VEGETATION", ylab = "Residuals", main = "Residuals by DOMINANT_VEGETATION")
abline(h = 0, lty = 2)

# Residuals by MONTH
boxplot(E1 ~ hr_noNA_noCNHGU$MONTH, xlab = "MONTH", ylab = "Residuals", main = "Residuals by MONTH")
abline(h = 0, lty = 2)

# res by year
boxplot(E1 ~ hr_noNA_noCNHGU$YEAR, xlab = "year", ylab = "Residuals", main = "Residuals by year")
abline(h = 0, lty = 2)

# Reset plotting layout
par(mfrow = c(1, 1))

# test normality
residuals <- residuals(m_hr_3)

qqPlot(residuals)
hist(residuals)
# looks normal

###### daily model

d_noNA_noCNHGU$MONTH <- as.factor(d_noNA_noCNHGU$MONTH)
d_noNA_noCNHGU$DOMINANT_VEGETATION <- as.factor(d_noNA_noCNHGU$DOMINANT_VEGETATION)
d_noNA_noCNHGU$SITE <- as.factor(d_noNA_noCNHGU$SITE)
d_noNA_noCNHGU$TIMESTAMP <- as_date(d_noNA_noCNHGU$TIMESTAMP)

d_noNA_noCNHGU$DAY <- as.numeric(d_noNA_noCNHGU$DAY)
d_noNA_noCNHGU$YearMonth <- format(d_noNA_noCNHGU$TIMESTAMP, "%Y-%m")


# full model

m_d <- lme(
  DELTA_FCH4_abs_yjtrans_ctr ~ SITE_TS_TOP_ctr +
    SITE_WTL_ctr +
    USTAR_MEDIAN_ctr +
    U_WIND_F_MEAN_ctr +
    V_WIND_F_MEAN_ctr +
    PA_F_MEDIAN_ctr +
    NEE_F_MEAN_ctr +
    VPD_F_MEDIAN_ctr +
    MONTH +
    DOMINANT_VEGETATION, 
  random = ~1 | SITE/YearMonth,
  correlation = corAR1(form = ~ DAY | SITE/YearMonth), 
  weights = varComb(varExp(form = ~PA_F_MEDIAN_ctr),
                    varExp(form = ~SITE_TS_TOP_ctr)),
  data = d_noNA_noCNHGU,
  method = "ML"
)

anova(m_d)
# remove u wind

m_d_2 <- lme(
  DELTA_FCH4_abs_yjtrans_ctr ~ SITE_TS_TOP_ctr +
    SITE_WTL_ctr +
    USTAR_MEDIAN_ctr +
    #U_WIND_F_MEAN_ctr +
    V_WIND_F_MEAN_ctr +
    PA_F_MEDIAN_ctr +
    NEE_F_MEAN_ctr +
    VPD_F_MEDIAN_ctr +
    MONTH, 
  random = ~1 | SITE/YearMonth,
  correlation = corAR1(form = ~ DAY | SITE/YearMonth), 
  weights = varComb(varExp(form = ~PA_F_MEDIAN_ctr),
                    varExp(form = ~SITE_TS_TOP_ctr)),
  data = d_noNA_noCNHGU,
  method = "ML"
)

# compare models

anova(m_d, m_d_2)
# p value = 0.0111 and AIC increased --> do not remove u.wind

#### ---> final model (refit with REML): 

m_d <- lme(
  DELTA_FCH4_abs_yjtrans_ctr ~ SITE_TS_TOP_ctr +
    SITE_WTL_ctr +
    USTAR_MEDIAN_ctr +
    U_WIND_F_MEAN_ctr +
    V_WIND_F_MEAN_ctr +
    PA_F_MEDIAN_ctr +
    NEE_F_MEAN_ctr +
    VPD_F_MEDIAN_ctr +
    MONTH +
    DOMINANT_VEGETATION, 
  random = ~1 | SITE/YearMonth,
  correlation = corAR1(form = ~ DAY | SITE/YearMonth), 
  weights = varComb(varExp(form = ~PA_F_MEDIAN_ctr),
                    varExp(form = ~SITE_TS_TOP_ctr)),
  data = d_noNA_noCNHGU,
  method = "REML"
)

##### effect sizes for the fixed effects

# Extract raw fixed-effect coefficients
raw_coefficients <- fixed.effects(m_d)

# Sort beta coefficients by absolute value
sorted_beta_absolute <- raw_coefficients[order(abs(raw_coefficients), decreasing = TRUE)]

# View sorted coefficients
sorted_beta_absolute

summary(m_d)
# random effect variation explained

# Extract variance components
var_components <- VarCorr(m_d)
var_components

# Total variance = random var + residual var
total_var <- 0.33383460 + 0.13793784 + 0.08845633
total_var <- 0.85255865 + 0.13187448 + 0.08399246
total_var <- 0.39844113 + 0.13226987 + 0.08388384

# site:
(0.33383460 / total_var) * 100 # 59.58898

# yearmonth

(0.13793784 / total_var) * 100 # 24.6217

r.squaredGLMM(m_d)

#### check the residuals again

E1 <- residuals(m_d_4, type = "pearson")
F1 <- fitted(m_d_4) 

# Diagnostic plots
par(mfrow = c(4, 4), mar = c(5, 5, 2, 2))

# Residuals vs Fitted
plot(F1, E1, xlab = "Fitted values", ylab = "Residuals", main = "Residuals vs Fitted")
abline(h = 0, lty = 2)
# het gen, but pretty much the same as the one we started with

# Residuals by SITE
plot(as.factor(d_noNA_noCNHGU$SITE), E1, xlab = "SITE", ylab = "Residuals", main = "Residuals by SITE")
abline(h = 0, lty = 2)
# okay

# Residuals by SITE_WTL
plot(d_noNA_noCNHGU$SITE_WTL_ctr, E1, xlab = "SITE_WTL", ylab = "Residuals", main = "Residuals by SITE_WTL")
abline(h = 0, lty = 2)

# res by site ts top
plot(d_noNA_noCNHGU$SITE_TS_TOP_ctr, E1, xlab = "SITE_TS_TOP", ylab = "Residuals", main = "Residuals by SITE_TS_TOP")
abline(h = 0, lty = 2)

# res by VPD
plot(d_noNA_noCNHGU$VPD_F_MEDIAN_ctr, E1, xlab = "VPD", ylab = "Residuals", main = "Residuals by VPD")
abline(h = 0, lty = 2)

# res by NEE
plot(d_noNA_noCNHGU$NEE_F_mean_ctr, E1, xlab = "NEE", ylab = "Residuals", main = "Residuals by NEE")
abline(h = 0, lty = 2)

# res by USTAR
plot(d_noNA_noCNHGU$USTAR_MEDIAN_ctr, E1, xlab = "USTAR", ylab = "Residuals", main = "Residuals by USTAR")
abline(h = 0, lty = 2)

# res by u.wind
plot(d_noNA_noCNHGU$U_WIND_F_MEAN_ctr, E1, xlab = "U_WIND_F", ylab = "Residuals", main = "Residuals by u.wind")
abline(h = 0, lty = 2)

# res by v.wind
plot(d_noNA_noCNHGU$V_WIND_F_MEAN_ctr, E1, xlab = "V_WIND_F", ylab = "Residuals", main = "Residuals by v.wind")
abline(h = 0, lty = 2)

# res by PA
plot(d_noNA_noCNHGU$PA_F_MEDIAN_ctr, E1, xlab = "PA", ylab = "Residuals", main = "Residuals by PA")
abline(h = 0, lty = 2)

# res by DOMINANT_VEGETATION
plot(as.factor(d_noNA_noCNHGU$DOMINANT_VEGETATION), E1, xlab = "DOMINANT_VEGETATION", ylab = "Residuals", main = "Residuals by DOMINANT_VEGETATION")
abline(h = 0, lty = 2)

# Residuals by MONTH
boxplot(E1 ~ d_noNA_noCNHGU$MONTH, xlab = "MONTH", ylab = "Residuals", main = "Residuals by MONTH")
abline(h = 0, lty = 2)

# res by year
boxplot(E1 ~ d_noNA_noCNHGU$YEAR, xlab = "year", ylab = "Residuals", main = "Residuals by year")
abline(h = 0, lty = 2)

# Reset plotting layout
par(mfrow = c(1, 1))

# test normality
residuals <- residuals(m_d)
shapiro.test(residuals)

qqPlot(residuals)
hist(residuals)

vif(m_d)

###### weekly

week_noNA_noCNHGU$MONTH <- as.factor(week_noNA_noCNHGU$MONTH)
week_noNA_noCNHGU$DOMINANT_VEGETATION <- as.factor(week_noNA_noCNHGU$DOMINANT_VEGETATION)
week_noNA_noCNHGU$SITE <- as.factor(week_noNA_noCNHGU$SITE)

week_noNA_noCNHGU$YearMonth <- paste(week_noNA_noCNHGU$YEAR, week_noNA_noCNHGU$MONTH, sep = "-")


# full model (no NEE)
m_week_nonee <- lme(
  DELTA_FCH4_abs_yjtrans_ctr ~ SITE_TS_TOP_ctr +
    SITE_WTL_ctr +
    USTAR_MEDIAN_ctr +
    U_WIND_F_MEAN_ctr +
    V_WIND_F_MEAN_ctr +
    PA_F_MEDIAN_ctr +
    #NEE_F_MEAN_ctr +
    VPD_F_MEDIAN_ctr +
    DOMINANT_VEGETATION + 
    MONTH, 
  random = ~1 | SITE/YearMonth, 
  correlation = corAR1(form = ~ Week_of_year | SITE/YearMonth), 
  weights = varExp(form = ~PA_F_MEDIAN_ctr),
  data = week_noNA_noCNHGU,
  method = "ML"
)

summary(m_week_nonee)

car::vif(m_week_nonee)

anova(m_week_nonee)
# remove WTL

m_week_nonee_2 <- lme(
  DELTA_FCH4_abs_yjtrans_ctr ~ SITE_TS_TOP_ctr +
    #SITE_WTL_ctr +
    USTAR_MEDIAN_ctr +
    U_WIND_F_MEAN_ctr +
    V_WIND_F_MEAN_ctr +
    PA_F_MEDIAN_ctr +
    #NEE_F_MEAN_ctr +
    VPD_F_MEDIAN_ctr +
    DOMINANT_VEGETATION + 
    MONTH, 
  random = ~1 | SITE/YearMonth, 
  correlation = corAR1(form = ~ Week_of_year | SITE/YearMonth), 
  weights = varExp(form = ~PA_F_MEDIAN_ctr),
  data = week_noNA_noCNHGU,
  method = "ML"
)

anova(m_week_nonee, m_week_nonee_2)
# AIC decreased and p > 0.05 --> remove wtd

anova(m_week_nonee_2)
# remove u.wind

m_week_nonee_3 <- lme(
  DELTA_FCH4_abs_yjtrans_ctr ~ SITE_TS_TOP_ctr +
    #SITE_WTL_ctr +
    USTAR_MEDIAN_ctr +
    #U_WIND_F_MEAN_ctr +
    V_WIND_F_MEAN_ctr +
    PA_F_MEDIAN_ctr +
    #NEE_F_MEAN_ctr +
    VPD_F_MEDIAN_ctr +
    DOMINANT_VEGETATION + 
    MONTH, 
  random = ~1 | SITE/YearMonth, 
  correlation = corAR1(form = ~ Week_of_year | SITE/YearMonth), 
  weights = varExp(form = ~PA_F_MEDIAN_ctr),
  data = week_noNA_noCNHGU,
  method = "ML"
)

# compare models
anova(m_week_nonee_2, m_week_nonee_3)
# p value 0.73 and AIC decreased --> remove u wind

anova(m_week_nonee_3)
# remove ustar

m_week_nonee_4 <- lme(
  DELTA_FCH4_abs_yjtrans_ctr ~ SITE_TS_TOP_ctr +
    #SITE_WTL_ctr +
    #USTAR_MEDIAN_ctr +
    #U_WIND_F_MEAN_ctr +
    V_WIND_F_MEAN_ctr +
    PA_F_MEDIAN_ctr +
    #NEE_F_MEAN_ctr +
    VPD_F_MEDIAN_ctr +
    DOMINANT_VEGETATION + 
    MONTH, 
  random = ~1 | SITE/YearMonth, 
  correlation = corAR1(form = ~ Week_of_year | SITE/YearMonth), 
  weights = varExp(form = ~PA_F_MEDIAN_ctr),
  data = week_noNA_noCNHGU,
  method = "ML"
)

anova(m_week_nonee_3, m_week_nonee_4)
# p value 0.3067, AIC decreased --> drop ustar

anova(m_week_nonee_4)
# remove v wind

m_week_nonee_5 <- lme(
  DELTA_FCH4_abs_yjtrans_ctr ~ SITE_TS_TOP_ctr +
    #SITE_WTL_ctr +
    #USTAR_MEDIAN_ctr +
    #U_WIND_F_MEAN_ctr +
    #V_WIND_F_MEAN_ctr +
    PA_F_MEDIAN_ctr +
    #NEE_F_MEAN_ctr +
    VPD_F_MEDIAN_ctr +
    DOMINANT_VEGETATION + 
    MONTH, 
  random = ~1 | SITE/YearMonth, 
  correlation = corAR1(form = ~ Week_of_year | SITE/YearMonth), 
  weights = varExp(form = ~PA_F_MEDIAN_ctr),
  data = week_noNA_noCNHGU,
  method = "ML"
)

# compare models

anova(m_week_nonee_4, m_week_nonee_5)
# p value 0.4 AIC decreased --> remove v wind

anova(m_week_nonee_5)
# remove month

m_week_nonee_6 <- lme(
  DELTA_FCH4_abs_yjtrans_ctr ~ SITE_TS_TOP_ctr +
    #SITE_WTL_ctr +
    #USTAR_MEDIAN_ctr +
    #U_WIND_F_MEAN_ctr +
    #V_WIND_F_MEAN_ctr +
    PA_F_MEDIAN_ctr +
    #NEE_F_MEAN_ctr +
    VPD_F_MEDIAN_ctr +
    DOMINANT_VEGETATION, 
  #MONTH, 
  random = ~1 | SITE/YearMonth, 
  correlation = corAR1(form = ~ Week_of_year | SITE/YearMonth), 
  weights = varExp(form = ~PA_F_MEDIAN_ctr),
  data = week_noNA_noCNHGU,
  method = "ML"
)

anova(m_week_nonee_5, m_week_nonee_6)
# p=0.6644, AIc decreased --> drop month

anova(m_week_nonee_6)
# remove vpd

m_week_nonee_7 <- lme(
  DELTA_FCH4_abs_yjtrans_ctr ~ SITE_TS_TOP_ctr +
    #SITE_WTL_ctr +
    #USTAR_MEDIAN_ctr +
    #U_WIND_F_MEAN_ctr +
    #V_WIND_F_MEAN_ctr +
    PA_F_MEDIAN_ctr +
    #NEE_F_MEAN_ctr +
    #VPD_F_MEDIAN_ctr +
    DOMINANT_VEGETATION, 
  #MONTH, 
  random = ~1 | SITE/YearMonth, 
  correlation = corAR1(form = ~ Week_of_year | SITE/YearMonth), 
  weights = varExp(form = ~PA_F_MEDIAN_ctr),
  data = week_noNA_noCNHGU,
  method = "ML"
)

anova(m_week_nonee_6, m_week_nonee_7)
# p=0.1744, AIC decreased --> remove vpd

anova(m_week_nonee_7)
# remove ts

m_week_nonee_8 <- lme(
  DELTA_FCH4_abs_yjtrans_ctr ~ #SITE_TS_TOP_ctr +
    #SITE_WTL_ctr +
    #USTAR_MEDIAN_ctr +
    #U_WIND_F_MEAN_ctr +
    #V_WIND_F_MEAN_ctr +
    PA_F_MEDIAN_ctr +
    #NEE_F_MEAN_ctr +
    #VPD_F_MEDIAN_ctr +
    DOMINANT_VEGETATION, 
  #MONTH, 
  random = ~1 | SITE/YearMonth, 
  correlation = corAR1(form = ~ Week_of_year | SITE/YearMonth), 
  weights = varExp(form = ~PA_F_MEDIAN_ctr),
  data = week_noNA_noCNHGU,
  method = "ML"
)

anova(m_week_nonee_7, m_week_nonee_8)
# p=0.8, AIC decreased --> remove ts

anova(m_week_nonee_8)
# remove PA

m_week_nonee_9 <- lme(
  DELTA_FCH4_abs_yjtrans_ctr ~ #SITE_TS_TOP_ctr +
    #SITE_WTL_ctr +
    #USTAR_MEDIAN_ctr +
    #U_WIND_F_MEAN_ctr +
    #V_WIND_F_MEAN_ctr +
    #PA_F_MEDIAN_ctr +
    #NEE_F_MEAN_ctr +
    #VPD_F_MEDIAN_ctr +
    DOMINANT_VEGETATION, 
  #MONTH, 
  random = ~1 | SITE/YearMonth, 
  correlation = corAR1(form = ~ Week_of_year | SITE/YearMonth), 
  weights = varExp(form = ~PA_F_MEDIAN_ctr),
  data = week_noNA_noCNHGU,
  method = "ML"
)

anova(m_week_nonee_8, m_week_nonee_9)
# p=0.051, AIC increased --> keep PA


### ---> final model (refit with REML: 

m_week_nonee_8 <- lme(
  DELTA_FCH4_abs_yjtrans_ctr ~ #SITE_TS_TOP_ctr +
    #SITE_WTL_ctr +
    #USTAR_MEDIAN_ctr +
    #U_WIND_F_MEAN_ctr +
    #V_WIND_F_MEAN_ctr +
    PA_F_MEDIAN_ctr +
    #NEE_F_MEAN_ctr +
    #VPD_F_MEDIAN_ctr +
    DOMINANT_VEGETATION, 
  #MONTH, 
  random = ~1 | SITE/YearMonth, 
  correlation = corAR1(form = ~ Week_of_year | SITE/YearMonth), 
  weights = varExp(form = ~PA_F_MEDIAN_ctr),
  data = week_noNA_noCNHGU,
  method = "REML"
)

anova(m_week_nonee_8)

#### full no VPD model

m_week_novpd <- lme(
  DELTA_FCH4_abs_yjtrans_ctr ~ SITE_TS_TOP_ctr +
    SITE_WTL_ctr +
    USTAR_MEDIAN_ctr +
    U_WIND_F_MEAN_ctr +
    V_WIND_F_MEAN_ctr +
    PA_F_MEDIAN_ctr +
    NEE_F_MEAN_ctr +
    #VPD_F_MEDIAN_ctr +
    DOMINANT_VEGETATION + 
    MONTH, 
  random = ~1 | SITE/YearMonth, 
  correlation = corAR1(form = ~ Week_of_year | SITE/YearMonth), 
  weights = varExp(form = ~PA_F_MEDIAN_ctr),
  data = week_noNA_noCNHGU,
  method = "REML"
)
vif(m_week_novpd)

anova(m_week_novpd)
# remove nee

m_week_novpd_2 <- lme(
  DELTA_FCH4_abs_yjtrans_ctr ~ SITE_TS_TOP_ctr +
    SITE_WTL_ctr +
    USTAR_MEDIAN_ctr +
    U_WIND_F_MEAN_ctr +
    V_WIND_F_MEAN_ctr +
    PA_F_MEDIAN_ctr +
    #NEE_F_MEAN_ctr +
    #VPD_F_MEDIAN_ctr +
    DOMINANT_VEGETATION + 
    MONTH, 
  random = ~1 | SITE/YearMonth, 
  correlation = corAR1(form = ~ Week_of_year | SITE/YearMonth), 
  weights = varExp(form = ~PA_F_MEDIAN_ctr),
  data = week_noNA_noCNHGU,
  method = "ML"
)

# compare models
anova(m_week_novpd, m_week_novpd_2)
# p value 0.7 and AIC decreased --> remove nee

anova(m_week_novpd_2)
# remove WTL

m_week_novpd_3 <- lme(
  DELTA_FCH4_abs_yjtrans_ctr ~ SITE_TS_TOP_ctr +
    #SITE_WTL_ctr +
    USTAR_MEDIAN_ctr +
    U_WIND_F_MEAN_ctr +
    V_WIND_F_MEAN_ctr +
    PA_F_MEDIAN_ctr +
    #NEE_F_MEAN_ctr +
    #VPD_F_MEDIAN_ctr +
    DOMINANT_VEGETATION + 
    MONTH, 
  random = ~1 | SITE/YearMonth, 
  correlation = corAR1(form = ~ Week_of_year | SITE/YearMonth), 
  weights = varExp(form = ~PA_F_MEDIAN_ctr),
  data = week_noNA_noCNHGU,
  method = "ML"
)
# compare models
anova(m_week_novpd_2, m_week_novpd_3)
# p value 0.3 and AIC decreased --> remove WTL

anova(m_week_novpd_3)
# remove u wind

m_week_novpd_4 <- lme(
  DELTA_FCH4_abs_yjtrans_ctr ~ SITE_TS_TOP_ctr +
    #SITE_WTL_ctr +
    USTAR_MEDIAN_ctr +
    #U_WIND_F_MEAN_ctr +
    V_WIND_F_MEAN_ctr +
    PA_F_MEDIAN_ctr +
    #NEE_F_MEAN_ctr +
    #VPD_F_MEDIAN_ctr +
    DOMINANT_VEGETATION + 
    MONTH, 
  random = ~1 | SITE/YearMonth, 
  correlation = corAR1(form = ~ Week_of_year | SITE/YearMonth), 
  weights = varExp(form = ~PA_F_MEDIAN_ctr),
  data = week_noNA_noCNHGU,
  method = "ML"
)

# compare models
anova(m_week_novpd_3, m_week_novpd_4)
# p value 0.38 and AIC decreased --> remove u wind

anova(m_week_novpd_4)
# remove v wind

m_week_novpd_5 <- lme(
  DELTA_FCH4_abs_yjtrans_ctr ~ SITE_TS_TOP_ctr +
    #SITE_WTL_ctr +
    USTAR_MEDIAN_ctr +
    #U_WIND_F_MEAN_ctr +
    #V_WIND_F_MEAN_ctr +
    PA_F_MEDIAN_ctr +
    #NEE_F_MEAN_ctr +
    #VPD_F_MEDIAN_ctr +
    DOMINANT_VEGETATION + 
    MONTH, 
  random = ~1 | SITE/YearMonth, 
  correlation = corAR1(form = ~ Week_of_year | SITE/YearMonth), 
  weights = varExp(form = ~PA_F_MEDIAN_ctr),
  data = week_noNA_noCNHGU,
  method = "ML"
)
# compare models
anova(m_week_novpd_4, m_week_novpd_5)
# p value 0.36 and AIC decreased --> remove v wind

anova(m_week_novpd_5)
# remove ustar

m_week_novpd_6 <- lme(
  DELTA_FCH4_abs_yjtrans_ctr ~ SITE_TS_TOP_ctr +
    #SITE_WTL_ctr +
    #USTAR_MEDIAN_ctr +
    #U_WIND_F_MEAN_ctr +
    #V_WIND_F_MEAN_ctr +
    PA_F_MEDIAN_ctr +
    #NEE_F_MEAN_ctr +
    #VPD_F_MEDIAN_ctr +
    DOMINANT_VEGETATION + 
    MONTH, 
  random = ~1 | SITE/YearMonth, 
  correlation = corAR1(form = ~ Week_of_year | SITE/YearMonth), 
  weights = varExp(form = ~PA_F_MEDIAN_ctr),
  data = week_noNA_noCNHGU,
  method = "ML"
)

# compare models
anova(m_week_novpd_5, m_week_novpd_6)
# p value 0.39 and AIC decreased --> remove ustar

anova(m_week_novpd_6)
# remove month

m_week_novpd_7 <- lme(
  DELTA_FCH4_abs_yjtrans_ctr ~ SITE_TS_TOP_ctr +
    #SITE_WTL_ctr +
    #USTAR_MEDIAN_ctr +
    #U_WIND_F_MEAN_ctr +
    #V_WIND_F_MEAN_ctr +
    PA_F_MEDIAN_ctr +
    #NEE_F_MEAN_ctr +
    #VPD_F_MEDIAN_ctr +
    DOMINANT_VEGETATION, 
  #MONTH, 
  random = ~1 | SITE/YearMonth, 
  correlation = corAR1(form = ~ Week_of_year | SITE/YearMonth), 
  weights = varExp(form = ~PA_F_MEDIAN_ctr),
  data = week_noNA_noCNHGU,
  method = "ML"
)

# compare models
anova(m_week_novpd_6, m_week_novpd_7)
# p value 0.65 and AIC decreased --> remove month

anova(m_week_novpd_7)
# remove ts

m_week_novpd_8 <- lme(
  DELTA_FCH4_abs_yjtrans_ctr ~ #SITE_TS_TOP_ctr +
    #SITE_WTL_ctr +
    #USTAR_MEDIAN_ctr +
    #U_WIND_F_MEAN_ctr +
    #V_WIND_F_MEAN_ctr +
    PA_F_MEDIAN_ctr +
    #NEE_F_MEAN_ctr +
    #VPD_F_MEDIAN_ctr +
    DOMINANT_VEGETATION, 
  #MONTH, 
  random = ~1 | SITE/YearMonth, 
  correlation = corAR1(form = ~ Week_of_year | SITE/YearMonth), 
  weights = varExp(form = ~PA_F_MEDIAN_ctr),
  data = week_noNA_noCNHGU,
  method = "ML"
)

# compare models
anova(m_week_novpd_7, m_week_novpd_8)
# p value 0.35 and AIC decreased --> remove ts

anova(m_week_novpd_8)
# remove PA

m_week_novpd_9 <- lme(
  DELTA_FCH4_abs_yjtrans_ctr ~ #SITE_TS_TOP_ctr +
    #SITE_WTL_ctr +
    #USTAR_MEDIAN_ctr +
    #U_WIND_F_MEAN_ctr +
    #V_WIND_F_MEAN_ctr +
    #PA_F_MEDIAN_ctr +
    #NEE_F_MEAN_ctr +
    #VPD_F_MEDIAN_ctr +
    DOMINANT_VEGETATION, 
  #MONTH, 
  random = ~1 | SITE/YearMonth, 
  correlation = corAR1(form = ~ Week_of_year | SITE/YearMonth), 
  weights = varExp(form = ~PA_F_MEDIAN_ctr),
  data = week_noNA_noCNHGU,
  method = "ML"
)
# compare models
anova(m_week_novpd_8, m_week_novpd_9)
# p value 0.051 and AIC increased --> dont remove PA

# final model (refit with REML):

m_week_novpd_8 <- lme(
  DELTA_FCH4_abs_yjtrans_ctr ~ #SITE_TS_TOP_ctr +
    #SITE_WTL_ctr +
    #USTAR_MEDIAN_ctr +
    #U_WIND_F_MEAN_ctr +
    #V_WIND_F_MEAN_ctr +
    PA_F_MEDIAN_ctr +
    #NEE_F_MEAN_ctr +
    #VPD_F_MEDIAN_ctr +
    DOMINANT_VEGETATION, 
  #MONTH, 
  random = ~1 | SITE/YearMonth, 
  correlation = corAR1(form = ~ Week_of_year | SITE/YearMonth), 
  weights = varExp(form = ~PA_F_MEDIAN_ctr),
  data = week_noNA_noCNHGU,
  method = "REML"
)

##### effect sizes for the fixed effects

# Extract raw fixed-effect coefficients
raw_coefficients <- fixed.effects(m_week_nonee_8)
raw_coefficients <- fixed.effects(m_week_novpd_8)
raw_coefficients <- fixed.effects(m_week_nonee)
raw_coefficients <- fixed.effects(m_week_novpd)

# Sort beta coefficients by absolute value
sorted_beta_absolute <- raw_coefficients[order(abs(raw_coefficients), decreasing = TRUE)]

# View sorted coefficients
sorted_beta_absolute

# random effect variation explained
summary(m_week_novpd_8)
summary(m_week_nonee_8)
summary(m_week_novpd)
summary(m_week_nonee)

# Extract variance components
VarCorr(m_week_nonee_8)

VarCorr(m_week_novpd_8)

VarCorr(m_week_novpd)

VarCorr(m_week_nonee)

# Total variance = random var + reisdual var
total_var <- 3.682265e-01 + 2.812642e-09 + 2.046916e-01 # novpd
total_var <- 3.569289e-01 + 3.571155e-09 + 2.030489e-01 # nonee
total_var <- 3.355515e-01 + 5.490914e-08 + 1.978534e-01 # nonee_8


# site:
(3.355515e-01 / total_var) * 100 # 62.90746
( 3.355515e-01/total_var) * 100 # 62.90746
(3.569289e-01/total_var) * 100 # 63.73983
(3.682265e-01/total_var) * 100 # 63.73983 novpd

# yearmonth

(5.490914e-08 / total_var) * 100 # 1.029408e-05 nonee_8
(3.571155e-09 / total_var) * 100 # 6.377315e-07 nonee
(2.812642e-09 / total_var) * 100 # 4.909326e-07 novpd

r.squaredGLMM(m_week_novpd_8)
r.squaredGLMM(m_week_nonee)
r.squaredGLMM(m_week_novpd)

#### check the residuals again


E1 <- residuals(m_week_nonee_8, type = "pearson")
F1 <- fitted(m_week_nonee_8) 

# Diagnostic plots
par(mfrow = c(4, 4), mar = c(5, 5, 2, 2))

# Residuals vs Fitted
plot(F1, E1, xlab = "Fitted values", ylab = "Residuals", main = "Residuals vs Fitted")
abline(h = 0, lty = 2)
# het gen, but pretty much the same as the one we started with

# Residuals by SITE
plot(as.factor(week_noNA_noCNHGU$SITE), E1, xlab = "SITE", ylab = "Residuals", main = "Residuals by SITE")
abline(h = 0, lty = 2)
# okay

# Residuals by SITE_WTL
plot(week_noNA_noCNHGU$SITE_WTL_ctr, E1, xlab = "SITE_WTL", ylab = "Residuals", main = "Residuals by SITE_WTL")
abline(h = 0, lty = 2)

# res by site ts top
plot(week_noNA_noCNHGU$SITE_TS_TOP_ctr, E1, xlab = "SITE_TS_TOP", ylab = "Residuals", main = "Residuals by SITE_TS_TOP")
abline(h = 0, lty = 2)

# res by VPD
plot(week_noNA_noCNHGU$VPD_F_MEDIAN_ctr, E1, xlab = "VPD", ylab = "Residuals", main = "Residuals by VPD")
abline(h = 0, lty = 2)

# res by NEE
plot(week_noNA_noCNHGU$NEE_F_mean_ctr, E1, xlab = "NEE", ylab = "Residuals", main = "Residuals by NEE")
abline(h = 0, lty = 2)

# res by USTAR
plot(week_noNA_noCNHGU$USTAR_MEDIAN_ctr, E1, xlab = "USTAR", ylab = "Residuals", main = "Residuals by USTAR")
abline(h = 0, lty = 2)

# res by u.wind
plot(week_noNA_noCNHGU$U_WIND_F_MEAN_ctr, E1, xlab = "U_WIND_F", ylab = "Residuals", main = "Residuals by u.wind")
abline(h = 0, lty = 2)

# res by v.wind
plot(week_noNA_noCNHGU$V_WIND_F_MEAN_ctr, E1, xlab = "V_WIND_F", ylab = "Residuals", main = "Residuals by v.wind")
abline(h = 0, lty = 2)

# res by PA
plot(week_noNA_noCNHGU$PA_F_MEDIAN_ctr, E1, xlab = "PA", ylab = "Residuals", main = "Residuals by PA")
abline(h = 0, lty = 2)

# DOMINANT_VEGETATION
plot(as.factor(week_noNA_noCNHGU$DOMINANT_VEGETATION), E1, xlab = "DOMINANT_VEGETATION", ylab = "Residuals", main = "Residuals by DOMINANT_VEGETATION")
abline(h = 0, lty = 2)

# Residuals by MONTH
boxplot(E1 ~ week_noNA_noCNHGU$MONTH, xlab = "MONTH", ylab = "Residuals", main = "Residuals by MONTH")
abline(h = 0, lty = 2)

# year
boxplot(E1 ~ week_noNA_noCNHGU$YEAR, xlab = "year", ylab = "Residuals", main = "Residuals by year")
abline(h = 0, lty = 2)

# Reset plotting layout
par(mfrow = c(1, 1))

# test normality
residuals <- residuals(m_week_nonee_8)

qqPlot(residuals)
hist(residuals)

#### monthly model

month_noNA_noCNHGU$Month_factor <- as.factor(month_noNA_noCNHGU$MONTH)
month_noNA_noCNHGU$Month_num <- as.numeric(month_noNA_noCNHGU$MONTH)

month_noNA_noCNHGU$DOMINANT_VEGETATION <- as.factor(month_noNA_noCNHGU$DOMINANT_VEGETATION)
month_noNA_noCNHGU$SITE <- as.factor(month_noNA_noCNHGU$SITE)

# full model
m_month <- lme(
  DELTA_FCH4_abs_yjtrans_ctr ~ SITE_TS_TOP_ctr +
    SITE_WTL_ctr +
    USTAR_MEDIAN_ctr +
    U_WIND_F_MEAN_ctr +
    V_WIND_F_MEAN_ctr +
    PA_F_MEDIAN_ctr +
    NEE_F_mean_ctr +
    VPD_F_MEDIAN_ctr +
    DOMINANT_VEGETATION + 
    Month_factor, 
  random = ~1 | SITE/YEAR, 
  correlation = corAR1(form = ~ Month_num | SITE/YEAR), 
  weights = varComb(varExp(form = ~U_WIND_F_MEAN_ctr),
                    varExp(form = ~USTAR_MEDIAN_ctr)),
  data = month_noNA_noCNHGU
)

vif(m_month)

# ts, vpd and ustar are problematic

# --> no vpd, ts and ustar

# ---> full model:

m_month <- lme(
  DELTA_FCH4_abs_yjtrans_ctr ~ 
    SITE_WTL_ctr +
    USTAR_MEDIAN_ctr +
    U_WIND_F_MEAN_ctr +
    V_WIND_F_MEAN_ctr +
    PA_F_MEDIAN_ctr +
    NEE_F_MEAN_ctr +
    DOMINANT_VEGETATION + 
    MONTH, 
  random = ~1 | SITE, 
  weights = varComb(varExp(form = ~U_WIND_F_MEAN_ctr),
                    varExp(form = ~PA_F_MEDIAN_ctr)),
  data = month_noNA_noCNHGU,
  method = "ML"
)

anova(m_month)
# remove ustar

m_month_2 <- lme(
  DELTA_FCH4_abs_yjtrans_ctr ~ 
    SITE_WTL_ctr +
    #USTAR_MEDIAN_ctr +
    U_WIND_F_MEAN_ctr +
    V_WIND_F_MEAN_ctr +
    PA_F_MEDIAN_ctr +
    NEE_F_MEAN_ctr +
    DOMINANT_VEGETATION + 
    MONTH, 
  random = ~1 | SITE, 
  weights = varComb(varExp(form = ~U_WIND_F_MEAN_ctr),
                    varExp(form = ~PA_F_MEDIAN_ctr)),
  data = month_noNA_noCNHGU,
  method = "ML"
)

# compare models

anova(m_month, m_month_2)
# p = 0.0087 and AIC increased --> do not remove ustar

# final model (refit with REML): 

m_month <- lme(
  DELTA_FCH4_abs_yjtrans_ctr ~ 
    SITE_WTL_ctr +
    USTAR_MEDIAN_ctr +
    U_WIND_F_MEAN_ctr +
    V_WIND_F_MEAN_ctr +
    PA_F_MEDIAN_ctr +
    NEE_F_MEAN_ctr +
    DOMINANT_VEGETATION + 
    MONTH, 
  random = ~1 | SITE, 
  weights = varComb(varExp(form = ~U_WIND_F_MEAN_ctr),
                    varExp(form = ~PA_F_MEDIAN_ctr)),
  data = month_noNA_noCNHGU,
  method = "REML"
)


##### effect sizes for the fixed effects

# Extract raw fixed-effect coefficients
raw_coefficients <- fixed.effects(m_month)

# Sort beta coefficients by absolute value
sorted_beta_absolute <- raw_coefficients[order(abs(raw_coefficients), decreasing = TRUE)]

# View sorted coefficients
sorted_beta_absolute

# random effect variation explained

# Extract variance components
VarCorr(m_month)

# Total variance = random var + reisdual var
total_var <- 3.899725e-01 + 6.683517e-10 + 1.255412e-01 
total_var <- 0.35601591 + 0.01783855 + 0.09551632 

# site:
(0.2327032 / total_var) * 100
(3.563880e-01 / total_var) * 100
(0.35601591 / total_var) * 100

# year

(6.683517e-10 / total_var) * 100
(1.778174e-09 / total_var) * 100
(0.1289042  / total_var) * 100

r.squaredGLMM(m_month)

anova(m_month)
summary(m_month)

#### check the residuals again

E1 <- residuals(m_month, type = "pearson")
F1 <- fitted(m_month) 

# Residuals vs Fitted
plot(F1, E1, xlab = "Fitted values", ylab = "Residuals", main = "Residuals vs Fitted")
abline(h = 0, lty = 2)

# Residuals by SITE
plot(as.factor(month_noNA_noCNHGU$SITE), E1, xlab = "SITE", ylab = "Residuals", main = "Residuals by SITE")
abline(h = 0, lty = 2)

# Residuals by SITE_WTL
plot(month_noNA_noCNHGU$SITE_WTL_ctr, E1, xlab = "SITE_WTL", ylab = "Residuals", main = "Residuals by SITE_WTL")
abline(h = 0, lty = 2)

# res by site ts top
plot(month_noNA_noCNHGU$SITE_TS_TOP_ctr, E1, xlab = "SITE_TS_TOP", ylab = "Residuals", main = "Residuals by SITE_TS_TOP")
abline(h = 0, lty = 2)

# res by VPD
plot(month_noNA_noCNHGU$VPD_F_MEDIAN_ctr, E1, xlab = "VPD", ylab = "Residuals", main = "Residuals by VPD")
abline(h = 0, lty = 2)

# res by NEE
plot(month_noNA_noCNHGU$NEE_F_mean_ctr, E1, xlab = "NEE", ylab = "Residuals", main = "Residuals by NEE")
abline(h = 0, lty = 2)

# res by USTAR
plot(month_noNA_noCNHGU$USTAR_MEDIAN_ctr, E1, xlab = "USTAR", ylab = "Residuals", main = "Residuals by USTAR")
abline(h = 0, lty = 2)

# res by u.wind
plot(month_noNA_noCNHGU$U_WIND_F_MEAN_ctr, E1, xlab = "U_WIND_F", ylab = "Residuals", main = "Residuals by u.wind")
abline(h = 0, lty = 2)

# res by v.wind
plot(month_noNA_noCNHGU$V_WIND_F_MEAN_ctr, E1, xlab = "V_WIND_F", ylab = "Residuals", main = "Residuals by v.wind")
abline(h = 0, lty = 2)

# res by PA
plot(month_noNA_noCNHGU$PA_F_MEDIAN_ctr, E1, xlab = "PA", ylab = "Residuals", main = "Residuals by PA")
abline(h = 0, lty = 2)

# DOMINANT_VEGETATION
plot(as.factor(month_noNA_noCNHGU$DOMINANT_VEGETATION), E1, xlab = "DOMINANT_VEGETATION", ylab = "Residuals", main = "Residuals by DOMINANT_VEGETATION")
abline(h = 0, lty = 2)

# Residuals by MONTH
boxplot(E1 ~ month_noNA_noCNHGU$MONTH, xlab = "MONTH", ylab = "Residuals", main = "Residuals by MONTH")
abline(h = 0, lty = 2)

# by year
boxplot(E1 ~ month_noNA_noCNHGU$YEAR, xlab = "year", ylab = "Residuals", main = "Residuals by year")
abline(h = 0, lty = 2)

# Reset plotting layout
par(mfrow = c(1, 1))

# test normality
residuals <- residuals(m_month)

qqPlot(residuals)
hist(residuals)


################ LOOCV ##################

### leave one site out

### monthly model

m_month <- lme(
  DELTA_FCH4_abs_yjtrans_ctr ~ 
    SITE_WTL_ctr +
    USTAR_MEDIAN_ctr +
    U_WIND_F_MEAN_ctr +
    V_WIND_F_MEAN_ctr +
    PA_F_MEDIAN_ctr +
    NEE_F_MEAN_ctr +
    DOMINANT_VEGETATION + 
    MONTH, 
  random = ~1 | SITE, 
  weights = varComb(varExp(form = ~U_WIND_F_MEAN_ctr),
                    varExp(form = ~PA_F_MEDIAN_ctr)),
  data = month_noNA_noCNHGU,
  method = "REML"
)

### need to remove DOMINANT_VEGETATION from the model because CV can't predict new sites because DOMINANT_VEGETATION is perfectly site-dependent

month_noNA_noCNHGU <- month_noNA_noCNHGU %>%
  mutate(SITE = factor(SITE),
         MONTH = factor(MONTH))

sites <- levels(month_noNA_noCNHGU$SITE)

cv_df <- data.frame()

for (s in sites) {
  
  train_data <- month_noNA_noCNHGU %>% filter(SITE != s)
  test_data  <- month_noNA_noCNHGU %>% filter(SITE == s)
  
  # IMPORTANT: drop unused factor levels in training set
  train_data <- droplevels(train_data)
  
  # Use only months that are actually present in training data
  train_month_levels <- levels(train_data$MONTH)
  
  # Re-factor MONTH in test set to training months only
  test_data$MONTH <- factor(test_data$MONTH, levels = train_month_levels)
  
  # Drop test rows with months not seen in training
  test_data_cv <- test_data[!is.na(test_data$MONTH), ]
  if (nrow(test_data_cv) == 0) next
  
  m_cv <- lme(
    DELTA_FCH4_abs_yjtrans_ctr ~ 
      SITE_WTL_ctr +
      USTAR_MEDIAN_ctr +
      U_WIND_F_MEAN_ctr +
      V_WIND_F_MEAN_ctr +
      PA_F_MEDIAN_ctr +
      NEE_F_MEAN_ctr +
      MONTH,
    random = ~1 | SITE,
    weights = varComb(
      varExp(form = ~U_WIND_F_MEAN_ctr),
      varExp(form = ~PA_F_MEDIAN_ctr)
    ),
    data = train_data,
    method = "REML"
  )
  
  preds <- predict(m_cv, newdata = test_data_cv, level = 0)
  
  tmp <- data.frame(
    SITE = s,
    MONTH = test_data_cv$MONTH,
    observed = test_data_cv$DELTA_FCH4_abs_yjtrans_ctr,
    predicted = as.numeric(preds)
  )
  
  cv_df <- rbind(cv_df, tmp)
}

# Metrics
# RMSE

sqrt(mean((cv_df$observed - cv_df$predicted)^2, na.rm = TRUE))
#MAE
mean(abs(cv_df$observed - cv_df$predicted), na.rm = TRUE)
#R2
1 - sum((cv_df$observed - cv_df$predicted)^2, na.rm = TRUE) /
  sum((cv_df$observed - mean(cv_df$observed, na.rm = TRUE))^2, na.rm = TRUE)


##### weekly model

# Make sure grouping variables are factors
week_noNA_noCNHGU <- week_noNA_noCNHGU %>%
  mutate(
    SITE = factor(SITE),
    YearMonth = factor(YearMonth)
  )

sites <- levels(week_noNA_noCNHGU$SITE)
cv_df <- data.frame()

for (s in sites) {
  
  train_data <- week_noNA_noCNHGU %>% filter(SITE != s)
  test_data  <- week_noNA_noCNHGU %>% filter(SITE == s)
  
  # Drop unused factor levels in the training data
  train_data <- droplevels(train_data)
  
  # Align YearMonth levels in test data to training (important for model.matrix)
  test_data$YearMonth <- factor(test_data$YearMonth, levels = levels(train_data$YearMonth))
  
  # Drop test rows with YearMonth not seen in training (rare, but can happen)
  test_data_cv <- test_data[!is.na(test_data$YearMonth), ]
  if (nrow(test_data_cv) == 0) next
  
  # Fit CV model (fixed effects only include PA_F_MEDIAN_ctr)
  # NOTE: DOMINANT_VEGETATION excluded here to avoid "new level" when a veg class is unique to held-out site
  m_cv <- lme(
    DELTA_FCH4_abs_yjtrans_ctr ~ PA_F_MEDIAN_ctr,
    random = ~1 | SITE/YearMonth,
    correlation = corAR1(form = ~ Week_of_year | SITE/YearMonth),
    weights = varExp(form = ~ PA_F_MEDIAN_ctr),
    data = train_data,
    method = "REML"
  )
  
  # Predict on held-out site using fixed effects only (new SITE has no random effects)
  preds <- predict(m_cv, newdata = test_data_cv, level = 0)
  
  tmp <- data.frame(
    SITE = s,
    YearMonth = test_data_cv$YearMonth,
    Week_of_year = test_data_cv$Week_of_year,
    observed = test_data_cv$DELTA_FCH4_abs_yjtrans_ctr,
    predicted = as.numeric(preds)
  )
  
  cv_df <- rbind(cv_df, tmp)
}

# Metrics
# RMSE
sqrt(mean((cv_df$observed - cv_df$predicted)^2, na.rm = TRUE))
# MAE
mean(abs(cv_df$observed - cv_df$predicted), na.rm = TRUE)
# R2
1 - sum((cv_df$observed - cv_df$predicted)^2, na.rm = TRUE) /
  sum((cv_df$observed - mean(cv_df$observed, na.rm = TRUE))^2, na.rm = TRUE)


### daily model

# Ensure grouping variables are factors
d_noNA_noCNHGU <- d_noNA_noCNHGU %>%
  mutate(
    SITE = factor(SITE),
    YearMonth = factor(YearMonth),
    MONTH = factor(MONTH)
  )

sites <- levels(d_noNA_noCNHGU$SITE)
cv_df <- data.frame()

for (s in sites) {
  
  # 1) Split data
  train_data <- d_noNA_noCNHGU %>% filter(SITE != s)
  test_data  <- d_noNA_noCNHGU %>% filter(SITE == s)
  
  # 2) Drop unused factor levels in training set
  train_data <- droplevels(train_data)
  
  # 3) Align factor levels in test set to what exists in training set
  # Any unseen levels will become NA and be dropped
  test_data$MONTH    <- factor(test_data$MONTH,    levels = levels(train_data$MONTH))
  test_data$YearMonth<- factor(test_data$YearMonth,levels = levels(train_data$YearMonth))
  
  # 4) Drop test rows that contain factor levels unseen in training
  test_data_cv <- test_data[!is.na(test_data$MONTH) & !is.na(test_data$YearMonth), ]
  if (nrow(test_data_cv) == 0) next
  
  # (Optional but often helps corAR1 behave consistently)
  train_data <- train_data[order(train_data$SITE, train_data$YearMonth, train_data$DAY), ]
  test_data_cv <- test_data_cv[order(test_data_cv$SITE, test_data_cv$YearMonth, test_data_cv$DAY), ]
  
  # 5) Fit CV model
  m_cv <- lme(
    DELTA_FCH4_abs_yjtrans_ctr ~
      SITE_TS_TOP_ctr +
      SITE_WTL_ctr +
      USTAR_MEDIAN_ctr +
      U_WIND_F_MEAN_ctr +
      V_WIND_F_MEAN_ctr +
      PA_F_MEDIAN_ctr +
      NEE_F_MEAN_ctr +
      VPD_F_MEDIAN_ctr +
      MONTH,
    random = ~1 | SITE/YearMonth,
    correlation = corAR1(form = ~ DAY | SITE/YearMonth),
    weights = varComb(
      varExp(form = ~PA_F_MEDIAN_ctr),
      varExp(form = ~SITE_TS_TOP_ctr)
    ),
    data = train_data,
    method = "REML"
  )
  
  # 6) Predict held-out SITE using fixed effects only
  preds <- predict(m_cv, newdata = test_data_cv, level = 0)
  
  # 7) Store results
  tmp <- data.frame(
    SITE = s,
    YearMonth = test_data_cv$YearMonth,
    MONTH = test_data_cv$MONTH,
    DAY = test_data_cv$DAY,
    observed = test_data_cv$DELTA_FCH4_abs_yjtrans_ctr,
    predicted = as.numeric(preds)
  )
  
  cv_df <- rbind(cv_df, tmp)
}

# Metrics
# RMSE
sqrt(mean((cv_df$observed - cv_df$predicted)^2, na.rm = TRUE))
# MAE
mean(abs(cv_df$observed - cv_df$predicted), na.rm = TRUE)
# R2
1 - sum((cv_df$observed - cv_df$predicted)^2, na.rm = TRUE) /
  sum((cv_df$observed - mean(cv_df$observed, na.rm = TRUE))^2, na.rm = TRUE)


#### half-hourly model with TS, no month

# Make sure grouping vars are factors
hh_noNA_noCNHGU <- hh_noNA_noCNHGU %>%
  mutate(
    SITE = factor(SITE),
    DATE = factor(DATE),
    HOUR = factor(HOUR)
  )

sites <- levels(hh_noNA_noCNHGU$SITE)
cv_df <- data.frame()

for (s in sites) {
  
  # 1) Split: leave one SITE out
  train_data <- hh_noNA_noCNHGU %>% filter(SITE != s)
  test_data  <- hh_noNA_noCNHGU %>% filter(SITE == s)
  
  # 2) Drop unused factor levels in training data (critical)
  train_data <- droplevels(train_data)
  
  # 3) Align factor levels in test data to training data
  #    (new levels become NA and will be dropped)
  test_data$HOUR <- factor(test_data$HOUR, levels = levels(train_data$HOUR))
  test_data$DATE <- factor(test_data$DATE, levels = levels(train_data$DATE))
  
  # 4) Drop test rows with unseen factor levels (common for DATE)
  test_data_cv <- test_data[!is.na(test_data$HOUR) & !is.na(test_data$DATE), ]
  if (nrow(test_data_cv) == 0) next
  
  # 5) Fit CV model (DOMINANT_VEGETATION already excluded)
  m_cv <- lme(
    DELTA_FCH4_abs_yjtrans_ctr ~
      SITE_TS_TOP_ctr +
      EC_USTAR_ctr +
      V_WIND_F_ctr +
      EC_PA_F_ctr +
      NEE_F_ctr +
      EC_VPD_F_ctr +
      HOUR,
    random = ~1 | SITE/DATE,
    weights = varComb(
      varExp(form = ~ EC_VPD_F_ctr),
      varExp(form = ~ EC_USTAR_ctr),
      varExp(form = ~ EC_PA_F_ctr)
    ),
    data = train_data,
    method = "REML"
  )
  
  # 6) Predict held-out site using fixed effects only
  preds <- predict(m_cv, newdata = test_data_cv, level = 0)
  
  # 7) Store fold predictions
  tmp <- data.frame(
    SITE = s,
    DATE = test_data_cv$DATE,
    HOUR = test_data_cv$HOUR,
    observed = test_data_cv$DELTA_FCH4_abs_yjtrans_ctr,
    predicted = as.numeric(preds)
  )
  
  cv_df <- rbind(cv_df, tmp)
}

# Metrics
# RMSE
sqrt(mean((cv_df$observed - cv_df$predicted)^2, na.rm = TRUE))
# MAE
mean(abs(cv_df$observed - cv_df$predicted), na.rm = TRUE)
# R2
1 - sum((cv_df$observed - cv_df$predicted)^2, na.rm = TRUE) /
  sum((cv_df$observed - mean(cv_df$observed, na.rm = TRUE))^2, na.rm = TRUE)

## half-hourly without ts, with month

# Make sure grouping vars are factors
hh_noNA_noCNHGU <- hh_noNA_noCNHGU %>%
  mutate(
    SITE  = factor(SITE),
    DATE  = factor(DATE),
    MONTH = factor(MONTH),
    HOUR  = factor(HOUR)
  )

sites <- levels(hh_noNA_noCNHGU$SITE)
cv_df <- data.frame()

for (s in sites) {
  
  # 1) Leave one SITE out
  train_data <- hh_noNA_noCNHGU %>% filter(SITE != s)
  test_data  <- hh_noNA_noCNHGU %>% filter(SITE == s)
  
  # 2) Drop unused factor levels in training set
  train_data <- droplevels(train_data)
  
  # 3) Align factor levels in test set to training set
  #    Any unseen levels become NA
  test_data$MONTH <- factor(test_data$MONTH, levels = levels(train_data$MONTH))
  test_data$HOUR  <- factor(test_data$HOUR,  levels = levels(train_data$HOUR))
  test_data$DATE  <- factor(test_data$DATE,  levels = levels(train_data$DATE))
  
  # 4) Drop rows with unseen factor levels
  test_data_cv <- test_data[!is.na(test_data$MONTH) & !is.na(test_data$HOUR) & !is.na(test_data$DATE), ]
  if (nrow(test_data_cv) == 0) next
  
  # 5) Fit CV model (DOMINANT_VEGETATION excluded)
  m_cv <- lme(
    DELTA_FCH4_abs_yjtrans_ctr ~ 
      SITE_WTL_ctr +
      EC_USTAR_ctr +
      V_WIND_F_ctr +
      EC_PA_F_ctr +
      NEE_F_ctr +
      EC_VPD_F_ctr + 
      HOUR +
      MONTH,
    random = ~1 | SITE/DATE,
    weights = varComb(
      varExp(form = ~ EC_VPD_F_ctr),
      varExp(form = ~ EC_USTAR_ctr),
      varExp(form = ~ EC_PA_F_ctr)
    ),
    method = "REML",
    data = train_data
  )
  
  # 6) Predict held-out SITE using fixed effects only
  preds <- predict(m_cv, newdata = test_data_cv, level = 0)
  
  # 7) Store results
  tmp <- data.frame(
    SITE = s,
    DATE = test_data_cv$DATE,
    MONTH = test_data_cv$MONTH,
    HOUR = test_data_cv$HOUR,
    observed = test_data_cv$DELTA_FCH4_abs_yjtrans_ctr,
    predicted = as.numeric(preds)
  )
  
  cv_df <- rbind(cv_df, tmp)
}

# Metrics
# RMSE
sqrt(mean((cv_df$observed - cv_df$predicted)^2, na.rm = TRUE))
# MAE
mean(abs(cv_df$observed - cv_df$predicted), na.rm = TRUE)
# R2
1 - sum((cv_df$observed - cv_df$predicted)^2, na.rm = TRUE) /
  sum((cv_df$observed - mean(cv_df$observed, na.rm = TRUE))^2, na.rm = TRUE)

## hourly model

# Ensure grouping vars are factors
hr_noNA_noCNHGU <- hr_noNA_noCNHGU %>%
  mutate(
    SITE  = factor(SITE),
    DATE  = factor(DATE),
    MONTH = factor(MONTH),
    HOUR  = factor(HOUR)
  )

sites <- levels(hr_noNA_noCNHGU$SITE)
cv_df <- data.frame()

for (s in sites) {
  
  # 1) Leave one SITE out
  train_data <- hr_noNA_noCNHGU %>% filter(SITE != s)
  test_data  <- hr_noNA_noCNHGU %>% filter(SITE == s)
  
  # 2) Drop unused factor levels in training set
  train_data <- droplevels(train_data)
  
  # 3) Align factor levels in test set to training set
  test_data$MONTH <- factor(test_data$MONTH, levels = levels(train_data$MONTH))
  test_data$HOUR  <- factor(test_data$HOUR,  levels = levels(train_data$HOUR))
  test_data$DATE  <- factor(test_data$DATE,  levels = levels(train_data$DATE))
  
  # 4) Drop test rows with unseen factor levels
  test_data_cv <- test_data[!is.na(test_data$MONTH) & !is.na(test_data$HOUR) & !is.na(test_data$DATE), ]
  if (nrow(test_data_cv) == 0) next
  
  # 5) Fit CV model (DOMINANT_VEGETATION excluded)
  m_cv <- lme(
    DELTA_FCH4_abs_yjtrans_ctr ~
      SITE_TS_TOP_ctr +
      SITE_WTL_ctr +
      USTAR_MEDIAN_ctr +
      V_WIND_F_MEAN_ctr +
      PA_F_MEDIAN_ctr +
      NEE_F_mean_ctr +
      VPD_F_MEDIAN_ctr +
      HOUR +
      MONTH,
    random = ~1 | SITE/DATE,
    weights = varComb(
      varExp(form = ~PA_F_MEDIAN_ctr),
      varExp(form = ~VPD_F_MEDIAN_ctr),
      varExp(form = ~USTAR_MEDIAN_ctr)
    ),
    data = train_data,
    method = "REML"
  )
  
  # 6) Predict held-out SITE using fixed effects only
  preds <- predict(m_cv, newdata = test_data_cv, level = 0)
  
  # 7) Store results
  tmp <- data.frame(
    SITE = s,
    DATE = test_data_cv$DATE,
    MONTH = test_data_cv$MONTH,
    HOUR = test_data_cv$HOUR,
    observed = test_data_cv$DELTA_FCH4_abs_yjtrans_ctr,
    predicted = as.numeric(preds)
  )
  
  cv_df <- rbind(cv_df, tmp)
}

# Metrics
# RMSE
sqrt(mean((cv_df$observed - cv_df$predicted)^2, na.rm = TRUE))
# MAE
mean(abs(cv_df$observed - cv_df$predicted), na.rm = TRUE)
# R2
1 - sum((cv_df$observed - cv_df$predicted)^2, na.rm = TRUE) /
  sum((cv_df$observed - mean(cv_df$observed, na.rm = TRUE))^2, na.rm = TRUE)

##### half-hourly full model

hh_noNA_noCNHGU <- hh_noNA_noCNHGU %>%
  mutate(
    SITE    = factor(SITE),
    DATE    = factor(DATE),
    HOUR    = factor(HOUR),
    MONTH   = factor(MONTH),
    DOMINANT_VEGETATION = factor(DOMINANT_VEGETATION)
  )

sites <- levels(hh_noNA_noCNHGU$SITE)
cv_df <- data.frame()

for (s in sites) {
  
  # 1) Leave one SITE out
  train_data <- hh_noNA_noCNHGU %>% filter(SITE != s)
  test_data  <- hh_noNA_noCNHGU %>% filter(SITE == s)
  
  # 2) Drop unused factor levels in training set
  train_data <- droplevels(train_data)
  
  # 3) Align factor levels in test set to training set
  test_data$HOUR    <- factor(test_data$HOUR,    levels = levels(train_data$HOUR))
  test_data$MONTH   <- factor(test_data$MONTH,   levels = levels(train_data$MONTH))
  test_data$DOMINANT_VEGETATION <- factor(test_data$DOMINANT_VEGETATION, levels = levels(train_data$DOMINANT_VEGETATION))
  test_data$DATE    <- factor(test_data$DATE,    levels = levels(train_data$DATE))
  
  # 4) Drop test rows with unseen factor levels
  test_data_cv <- test_data[
    !is.na(test_data$HOUR) &
      !is.na(test_data$MONTH) &
      !is.na(test_data$DOMINANT_VEGETATION) &
      !is.na(test_data$DATE),
  ]
  if (nrow(test_data_cv) == 0) next
  
  # 5) Fit CV model on training data
  m_cv <- lme(
    DELTA_FCH4_abs_yjtrans_ctr ~
      SITE_WTL_ctr +
      EC_USTAR_ctr +
      U_WIND_F_ctr +
      V_WIND_F_ctr +
      EC_PA_F_ctr +
      NEE_F_ctr +
      EC_VPD_F_ctr +
      DOMINANT_VEGETATION +
      HOUR +
      MONTH,
    random = ~1 | SITE/DATE,
    weights = varComb(
      varExp(form = ~ EC_VPD_F_ctr),
      varExp(form = ~ EC_USTAR_ctr),
      varExp(form = ~ EC_PA_F_ctr)
    ),
    data = train_data,
    method = "REML"
  )
  
  # 6) Predict for held-out SITE using fixed effects only
  preds <- predict(m_cv, newdata = test_data_cv, level = 0)
  
  # 7) Store results
  tmp <- data.frame(
    SITE = s,
    DATE = test_data_cv$DATE,
    MONTH = test_data_cv$MONTH,
    HOUR = test_data_cv$HOUR,
    DOMINANT_VEGETATION = test_data_cv$DOMINANT_VEGETATION,
    observed = test_data_cv$DELTA_FCH4_abs_yjtrans_ctr,
    predicted = as.numeric(preds)
  )
  
  cv_df <- rbind(cv_df, tmp)
}

# Metrics
# RMSE
sqrt(mean((cv_df$observed - cv_df$predicted)^2, na.rm = TRUE))
# MAE
mean(abs(cv_df$observed - cv_df$predicted), na.rm = TRUE)
# R2
1 - sum((cv_df$observed - cv_df$predicted)^2, na.rm = TRUE) /
  sum((cv_df$observed - mean(cv_df$observed, na.rm = TRUE))^2, na.rm = TRUE)


# hourly full model

hr_noNA_noCNHGU <- hr_noNA_noCNHGU %>%
  mutate(
    SITE    = factor(SITE),
    DATE    = factor(DATE),
    HOUR    = factor(HOUR),
    MONTH   = factor(MONTH),
    DOMINANT_VEGETATION = factor(DOMINANT_VEGETATION)
  )

sites <- levels(hr_noNA_noCNHGU$SITE)
cv_df <- data.frame()

for (s in sites) {
  
  # 1) Leave one SITE out
  train_data <- hr_noNA_noCNHGU %>% filter(SITE != s)
  test_data  <- hr_noNA_noCNHGU %>% filter(SITE == s)
  
  # 2) Drop unused factor levels in training set
  train_data <- droplevels(train_data)
  
  # 3) Align factor levels in test set to training set
  test_data$HOUR    <- factor(test_data$HOUR,    levels = levels(train_data$HOUR))
  test_data$MONTH   <- factor(test_data$MONTH,   levels = levels(train_data$MONTH))
  test_data$DOMINANT_VEGETATION <- factor(test_data$DOMINANT_VEGETATION, levels = levels(train_data$DOMINANT_VEGETATION))
  test_data$DATE    <- factor(test_data$DATE,    levels = levels(train_data$DATE))
  
  # 4) Drop test rows with unseen factor levels
  test_data_cv <- test_data[
    !is.na(test_data$HOUR) &
      !is.na(test_data$MONTH) &
      !is.na(test_data$DOMINANT_VEGETATION) &
      !is.na(test_data$DATE),
  ]
  if (nrow(test_data_cv) == 0) next
  
  # 5) Fit CV model on training data
  m_cv <- lme(
    DELTA_FCH4_abs_yjtrans_ctr ~
      SITE_TS_TOP_ctr +
      SITE_WTL_ctr +
      USTAR_MEDIAN_ctr +
      U_WIND_F_MEAN_ctr +
      V_WIND_F_MEAN_ctr +
      PA_F_MEDIAN_ctr +
      NEE_F_mean_ctr +
      VPD_F_MEDIAN_ctr +
      HOUR +
      DOMINANT_VEGETATION +
      MONTH,
    random = ~1 | SITE/DATE,
    weights = varComb(
      varExp(form = ~PA_F_MEDIAN_ctr),
      varExp(form = ~VPD_F_MEDIAN_ctr),
      varExp(form = ~USTAR_MEDIAN_ctr)
    ),
    data = train_data,
    method = "REML",
    control = lmeControl(opt = "nlminb", maxIter = 100, msMaxIter = 100)
  )
  
  # 6) Predict for held-out SITE using fixed effects only
  preds <- predict(m_cv, newdata = test_data_cv, level = 0)
  
  # 7) Store results
  tmp <- data.frame(
    SITE = s,
    DATE = test_data_cv$DATE,
    MONTH = test_data_cv$MONTH,
    HOUR = test_data_cv$HOUR,
    DOMINANT_VEGETATION = test_data_cv$DOMINANT_VEGETATION,
    observed = test_data_cv$DELTA_FCH4_abs_yjtrans_ctr,
    predicted = as.numeric(preds)
  )
  
  cv_df <- rbind(cv_df, tmp)
}

# Metrics
# RMSE
sqrt(mean((cv_df$observed - cv_df$predicted)^2, na.rm = TRUE))
# MAE
mean(abs(cv_df$observed - cv_df$predicted), na.rm = TRUE)
# R2
1 - sum((cv_df$observed - cv_df$predicted)^2, na.rm = TRUE) /
  sum((cv_df$observed - mean(cv_df$observed, na.rm = TRUE))^2, na.rm = TRUE)


### weekly model, full, no nee

# Ensure grouping vars are factors
week_noNA_noCNHGU <- week_noNA_noCNHGU %>%
  mutate(
    SITE       = factor(SITE),
    YearMonth  = factor(YearMonth),
    MONTH      = factor(MONTH)
  )

sites <- levels(week_noNA_noCNHGU$SITE)
cv_df <- data.frame()

for (s in sites) {
  
  # 1) Leave one SITE out
  train_data <- week_noNA_noCNHGU %>% filter(SITE != s)
  test_data  <- week_noNA_noCNHGU %>% filter(SITE == s)
  
  # 2) Drop unused factor levels in training set
  train_data <- droplevels(train_data)
  
  # 3) Align factor levels in test set to training set
  test_data$YearMonth <- factor(test_data$YearMonth, levels = levels(train_data$YearMonth))
  test_data$MONTH     <- factor(test_data$MONTH,     levels = levels(train_data$MONTH))
  
  # 4) Drop test rows with unseen factor levels
  test_data_cv <- test_data[
    !is.na(test_data$YearMonth) &
      !is.na(test_data$MONTH),
  ]
  if (nrow(test_data_cv) == 0) next
  
  # (Optional but helps corAR1 behave consistently)
  train_data <- train_data[order(train_data$SITE, train_data$YearMonth, train_data$Week_of_year), ]
  test_data_cv <- test_data_cv[order(test_data_cv$SITE, test_data_cv$YearMonth, test_data_cv$Week_of_year), ]
  
  # 5) Fit CV model on training data
  m_cv <- lme(
    DELTA_FCH4_abs_yjtrans_ctr ~
      SITE_TS_TOP_ctr +
      SITE_WTL_ctr +
      USTAR_MEDIAN_ctr +
      U_WIND_F_MEAN_ctr +
      V_WIND_F_MEAN_ctr +
      PA_F_MEDIAN_ctr +
      VPD_F_MEDIAN_ctr +
      MONTH,
    random = ~1 | SITE/YearMonth,
    correlation = corAR1(form = ~ Week_of_year | SITE/YearMonth),
    weights = varExp(form = ~PA_F_MEDIAN_ctr),
    data = train_data,
    method = "REML"
  )
  
  # 6) Predict held-out SITE using fixed effects only
  preds <- predict(m_cv, newdata = test_data_cv, level = 0)
  
  # 7) Store results
  tmp <- data.frame(
    SITE = s,
    YearMonth = test_data_cv$YearMonth,
    MONTH = test_data_cv$MONTH,
    Week_of_year = test_data_cv$Week_of_year,
    observed = test_data_cv$DELTA_FCH4_abs_yjtrans_ctr,
    predicted = as.numeric(preds)
  )
  
  cv_df <- rbind(cv_df, tmp)
}

# Metrics
# RMSE
sqrt(mean((cv_df$observed - cv_df$predicted)^2, na.rm = TRUE))
# MAE
mean(abs(cv_df$observed - cv_df$predicted), na.rm = TRUE)
# R2
1 - sum((cv_df$observed - cv_df$predicted)^2, na.rm = TRUE) /
  sum((cv_df$observed - mean(cv_df$observed, na.rm = TRUE))^2, na.rm = TRUE)


# weekly model, full, no vpd

# Ensure grouping vars are factors
week_noNA_noCNHGU <- week_noNA_noCNHGU %>%
  mutate(
    SITE      = factor(SITE),
    YearMonth = factor(YearMonth),
    MONTH     = factor(MONTH)
  )

sites <- levels(week_noNA_noCNHGU$SITE)
cv_df <- data.frame()

for (s in sites) {
  
  # 1) Leave one SITE out
  train_data <- week_noNA_noCNHGU %>% filter(SITE != s)
  test_data  <- week_noNA_noCNHGU %>% filter(SITE == s)
  
  # 2) Drop unused factor levels in training set
  train_data <- droplevels(train_data)
  
  # 3) Align factor levels in test set to training set
  test_data$YearMonth <- factor(test_data$YearMonth, levels = levels(train_data$YearMonth))
  test_data$MONTH     <- factor(test_data$MONTH,     levels = levels(train_data$MONTH))
  
  # 4) Drop test rows with unseen factor levels
  test_data_cv <- test_data[!is.na(test_data$YearMonth) & !is.na(test_data$MONTH), ]
  if (nrow(test_data_cv) == 0) next
  
  # (Optional but helps corAR1 behave consistently)
  train_data   <- train_data[order(train_data$SITE, train_data$YearMonth, train_data$Week_of_year), ]
  test_data_cv <- test_data_cv[order(test_data_cv$SITE, test_data_cv$YearMonth, test_data_cv$Week_of_year), ]
  
  # 5) Fit CV model on training data
  m_cv <- lme(
    DELTA_FCH4_abs_yjtrans_ctr ~
      SITE_TS_TOP_ctr +
      SITE_WTL_ctr +
      USTAR_MEDIAN_ctr +
      U_WIND_F_MEAN_ctr +
      V_WIND_F_MEAN_ctr +
      PA_F_MEDIAN_ctr +
      NEE_F_MEAN_ctr +
      MONTH,
    random = ~1 | SITE/YearMonth,
    correlation = corAR1(form = ~ Week_of_year | SITE/YearMonth),
    weights = varExp(form = ~PA_F_MEDIAN_ctr),
    data = train_data,
    method = "REML"
  )
  
  # 6) Predict held-out SITE using fixed effects only
  preds <- predict(m_cv, newdata = test_data_cv, level = 0)
  
  # 7) Store results
  tmp <- data.frame(
    SITE = s,
    YearMonth = test_data_cv$YearMonth,
    MONTH = test_data_cv$MONTH,
    Week_of_year = test_data_cv$Week_of_year,
    observed = test_data_cv$DELTA_FCH4_abs_yjtrans_ctr,
    predicted = as.numeric(preds)
  )
  
  cv_df <- rbind(cv_df, tmp)
}

# Metrics
sqrt(mean((cv_df$observed - cv_df$predicted)^2, na.rm = TRUE))
mean(abs(cv_df$observed - cv_df$predicted), na.rm = TRUE)
1 - sum((cv_df$observed - cv_df$predicted)^2, na.rm = TRUE) /
  sum((cv_df$observed - mean(cv_df$observed, na.rm = TRUE))^2, na.rm = TRUE)


#################################################################################################################################################
#### MOVE THIS TO OTHER SCRIPTS? 
# Create the new column SITE_TS_TOP
all_d_aggr <- all_d_aggr %>%
  mutate(SITE_TS_TOP = case_when(
    SITE %in% c("FI-SI2", "US-HO1") ~ ch_TS_TOP,
    SITE %in% c("US-OWC", "CN-HGU") ~ EC_TS_TOP,
    TRUE ~ rowMeans(cbind(ch_TS_TOP, EC_TS_TOP), na.rm = TRUE)
  ))


# EC_WTL is m so convert EC_WTL to cm

all_d_aggr$EC_WTL_F_mean <- all_d_aggr$EC_WTL_F_mean * 100

# Create the new column SITE_WTL
all_d_aggr <- all_d_aggr %>%
  mutate(SITE_WTL = case_when(
    SITE %in% c("FI-SI2", "US-LA1", "US-LA2") ~ ch_WTL_mean,
    SITE %in% c("US-HO1", "US-OWC", "US-UAF", "US-LOS") ~ EC_WTL_F_mean,
    SITE == "CN-HGU" ~ NA,
    TRUE ~ rowMeans(cbind(ch_WTL_mean, EC_WTL_F_mean), na.rm = TRUE)
  ))


###############################################################################################################################

# Fig. B1

all_sites_pred_dupl <- all_sites_pred_dupl %>% mutate(Month_letter =
                                                        case_when(all_sites_pred_dupl$MONTH == 06 ~ "Jun", 
                                                                  all_sites_pred_dupl$MONTH == 07 ~ "Jul",
                                                                  all_sites_pred_dupl$MONTH == 08 ~ "Aug",
                                                                  all_sites_pred_dupl$MONTH == 09 ~ "Sep",
                                                                  all_sites_pred_dupl$MONTH == 10 ~ "Oct",
                                                                  all_sites_pred_dupl$MONTH == 11 ~ "Nov",
                                                                  all_sites_pred_dupl$MONTH == 12 ~ "Dec",
                                                                  all_sites_pred_dupl$MONTH == 01 ~ "Jan",
                                                                  all_sites_pred_dupl$MONTH == 02 ~ "Feb",
                                                                  all_sites_pred_dupl$MONTH == 03 ~ "Mar",
                                                                  all_sites_pred_dupl$MONTH == 04 ~ "Apr",
                                                                  all_sites_pred_dupl$MONTH == 05 ~ "May")
)

aggregated_data <- all_sites_pred_dupl %>%
  group_by(SITE, Month_letter) %>%
  summarize(
    unique_years = n_distinct(YEAR),       # Count unique "YEAR"
    unique_ch_IDs = n_distinct(ch_ID),     # Count unique "ch_ID"
    .groups = "drop"
  )

aggregated_data <- all_sites_pred_dupl %>%
  group_by(SITE, Month_letter, YEAR) %>%
  summarize(unique_ch_IDs_year = n_distinct(ch_ID), .groups = "drop") %>%
  group_by(SITE, Month_letter) %>%
  summarize(
    avg_unique_ch_IDs = ceiling(mean(unique_ch_IDs_year)), # Round to 1 decimal place
    unique_years = n_distinct(YEAR),                       # Count unique years
    .groups = "drop"
  )

aggregated_data

aggregated_data$avg_unique_ch_IDs <- as.factor(aggregated_data$avg_unique_ch_IDs)

colors <- magma(11, end = 0.9)
print(colors)

safe_palette <- c(
  "#88CCEE", # Light blue
  "#CC6677", # Light red
  "#DDCC77", # Yellow
  "#117733", # Dark green
  "#332288", # Dark blue
  "#AA4499", # Purple
  "#44AA99", # Teal
  "#999933", # Olive
  "#882255", # Wine
  "#661100"  # Brown
)
# Add two additional colors
extended_safe_palette <- c(safe_palette, "#E69F00", "darkgrey") # Orange and Sky Blue

# Use as many colors as needed
safe_colors <- extended_safe_palette[1:n_distinct(aggregated_data$avg_unique_ch_IDs)]

aggregated_data$Month_letter <- factor(
  aggregated_data$Month_letter,
  levels = c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")
)



aggregated_data <- as.data.frame(aggregated_data)

aggregated_data$SITE <- factor(
  aggregated_data$SITE,
  levels = c("CN-HGU", "FI-SI2", "SE-DEG", "US-HO1", "US-LA1", "US-LA2", "US-LOS", "US-OWC", "US-STJ","US-UAF")
)

ggplot(aggregated_data, aes(x = Month_letter, y = SITE)) +
  geom_point(aes(size = unique_years, color = avg_unique_ch_IDs)) +
  scale_color_manual(
    name = "Average no. chambers",
    values = safe_colors
  ) +
  scale_size_continuous(name = "No. years") +
  theme_bw() +
  guides(
    color = guide_legend(
      override.aes = list(size = 5) # Increase legend point size
    )
  ) +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(size=20, angle = 45, hjust = 1),
        axis.text.y = element_text(size=20),
        legend.text = element_text(size=20),
        legend.title = element_text(size = 20)
  )

#########################################################################################################################################

##### manual vs auto plot (Fig. B19)

df_auto <- all_d_aggr %>% filter(CH_METHOD == "auto")
df_manual <- all_d_aggr %>% filter(CH_METHOD == "manual")

# Compute Spearman correlation separately for each group
cor.test(df_auto$EC_FCH4_MEDIAN, df_auto$CH_FCH4_MEDIAN, method = "spearman")
cor.test(df_manual$EC_FCH4_MEDIAN, df_manual$CH_FCH4_MEDIAN, method = "spearman")

# plot automated chamber vs EC
auto_EC <- ggplot(all_d_aggr %>% filter(CH_METHOD == "auto"), 
                  aes(x = CH_FCH4_MEDIAN, y = EC_FCH4_MEDIAN)) +
  geom_point(size = 3, shape = 21, color = "black", fill = "grey", stroke = 0.6) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black", linewidth = 0.5) +  # 1:1 line
  stat_cor(method = "spearman", label.x.npc = "left", label.y.npc = "top", , 
           cor.coef.name = "rho", size = 6) +
  theme_bw() +
  labs(x = expression(paste("Chamber FCH"[4]," (nmol m"^-2, "s"^-1, ")")),
       y = expression(paste("EC FCH"[4]," (nmol m"^-2, "s"^-1, ")")),
       title = "Automated") +
  theme(axis.text = element_text(size = 18),
        axis.title = element_text(size = 18),
        plot.title = element_text(size = 18, hjust = 0.5)) + 
  scale_x_continuous(limits = c(-2, 100)) +
  scale_y_continuous(limits=c(-2, 100))
auto_EC



# manual chamber FCH4 vs EC FCH4
manual_EC <- ggplot(all_d_aggr %>% filter(CH_METHOD == "manual"), 
                    aes(x = CH_FCH4_MEDIAN, y = EC_FCH4_MEDIAN)) +
  geom_point(size = 3, shape = 21, color = "black", fill = "grey", stroke = 0.6) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black", linewidth = 0.5) +
  stat_cor(method = "spearman", label.x.npc = "center", label.y.npc = "top", 
           cor.coef.name = "rho", size = 6) +
  theme_bw() +
  labs(x = expression(paste("Chamber FCH"[4]," (nmol m"^-2, "s"^-1, ")")),
       y = expression(paste("EC FCH"[4]," (nmol m"^-2, "s"^-1, ")")),
       title = "Manual") +
  theme(axis.text = element_text(size = 18),
        axis.title = element_text(size = 18),
        plot.title = element_text(size = 18, hjust = 0.5)) +
  scale_y_continuous(limits = c(-2, 3300)) +
  scale_x_continuous(limits = c(-2, 3300))
manual_EC

# weekly

# automated
auto_EC_weekly <- ggplot(all_week_aggr %>% filter(CH_METHOD == "auto"), 
                         aes(x = CH_FCH4_MEDIAN, y = EC_FCH4_MEDIAN)) +
  geom_point(size = 3, shape = 21, color = "black", fill = "grey", stroke = 0.6) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black", linewidth = 0.5) +
  stat_cor(method = "spearman", label.x.npc = "left", label.y.npc = "top", , 
           cor.coef.name = "rho", size = 6) +  
  theme_bw() +
  labs(x = expression(paste("Chamber FCH"[4]," (nmol m"^-2, "s"^-1, ")")),
       y = expression(paste("EC FCH"[4]," (nmol m"^-2, "s"^-1, ")")),
       title = "Automated") +
  theme(axis.text = element_text(size = 18),
        axis.title = element_text(size = 18),
        plot.title = element_text(size = 18, hjust = 0.5)) + 
  scale_x_continuous(limits = c(-2, 80)) +
  scale_y_continuous(limits=c(-2, 80))
auto_EC_weekly

# manual
manual_EC_weekly <- ggplot(all_week_aggr %>% filter(CH_METHOD == "manual"), 
                           aes(x = CH_FCH4_MEDIAN, y = EC_FCH4_MEDIAN)) +
  geom_point(size = 3, shape = 21, color = "black", fill = "grey", stroke = 0.6) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black", linewidth = 0.5) +
  stat_cor(method = "spearman", label.x.npc = "center", label.y.npc = "top", 
           cor.coef.name = "rho", size = 6) + 
  theme_bw() +
  labs(x = expression(paste("Chamber FCH"[4]," (nmol m"^-2, "s"^-1, ")")),
       y = expression(paste("EC FCH"[4]," (nmol m"^-2, "s"^-1, ")")),
       title = "Manual") +
  theme(axis.text = element_text(size = 18),
        axis.title = element_text(size = 18),
        plot.title = element_text(size = 18, hjust = 0.5)) +
  scale_y_continuous(limits = c(-2, 2000)) + 
  scale_x_continuous(limits = c(-2, 2000)) 
manual_EC_weekly

# monthly

# automated
auto_EC_monthly <- ggplot(all_month_aggr %>% filter(CH_METHOD == "auto"), 
                          aes(x = CH_FCH4_MEDIAN, y = EC_FCH4_MEDIAN)) +
  geom_point(size = 3, shape = 21, color = "black", fill = "grey", stroke = 0.6) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black", linewidth = 0.5) +
  stat_cor(method = "spearman", label.x.npc = "center", label.y.npc = "top", 
           cor.coef.name = "rho", size = 6) + 
  theme_bw() +
  labs(x = expression(paste("Chamber FCH"[4]," (nmol m"^-2, "s"^-1, ")")),
       y = expression(paste("EC FCH"[4]," (nmol m"^-2, "s"^-1, ")")),
       title = "Automated") +
  theme(axis.text = element_text(size = 18),
        axis.title = element_text(size = 18),
        plot.title = element_text(size = 18, hjust = 0.5)) +
  scale_y_continuous(limits = c(-2, 20)) + 
  scale_x_continuous(limits = c(-2, 20))
auto_EC_monthly

# manual
manual_EC_monthly <- ggplot(all_month_aggr %>% filter(CH_METHOD == "manual"), 
                            aes(x = CH_FCH4_MEDIAN, y = EC_FCH4_MEDIAN)) +
  geom_point(size = 3, shape = 21, color = "black", fill = "grey", stroke = 0.6) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black", linewidth = 0.5) +
  stat_cor(method = "spearman", label.x.npc = "center", label.y.npc = "top", 
           cor.coef.name = "rho", size = 6) + 
  theme_bw() +
  labs(x = expression(paste("Chamber FCH"[4]," (nmol m"^-2, "s"^-1, ")")),
       y = expression(paste("EC FCH"[4]," (nmol m"^-2, "s"^-1, ")")),
       title = "Manual") +
  theme(axis.text = element_text(size = 18),
        axis.title = element_text(size = 18),
        plot.title = element_text(size = 18, hjust = 0.5)) +
  scale_y_continuous(limits = c(-2, 2000)) + 
  scale_x_continuous(limits = c(-2, 2000)) 
manual_EC_monthly

# annual
# manual
manual_EC_annual <- ggplot(all_yr_aggr %>% filter(CH_METHOD == "manual"), 
                           aes(x = CH_FCH4_MEDIAN, y = EC_FCH4_MEDIAN)) +
  geom_point(size = 3, shape = 21, color = "black", fill = "grey", stroke = 0.6) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black", linewidth = 0.5) +
  stat_cor(method = "spearman", label.x.npc = "center", label.y.npc = "top", 
           cor.coef.name = "rho", size = 6) + 
  theme_bw() +
  labs(x = expression(paste("Chamber FCH"[4]," (nmol m"^-2, "s"^-1, ")")),
       y = expression(paste("EC FCH"[4]," (nmol m"^-2, "s"^-1, ")")),
       title = "Manual") +
  theme(axis.text = element_text(size = 18),
        axis.title = element_text(size = 18),
        plot.title = element_text(size = 18, hjust = 0.5)) +  
  scale_y_continuous(limits = c(-2, 1000)) + 
  scale_x_continuous(limits = c(-2, 1000)) 
manual_EC_annual

# automated
auto_EC_annual <- ggplot(all_yr_aggr %>% filter(CH_METHOD == "auto"), 
                         aes(x = CH_FCH4_MEDIAN, y = EC_FCH4_MEDIAN)) +
  geom_point(size = 3, shape = 21, color = "black", fill = "grey", stroke = 0.6) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black", linewidth = 0.5) + 
  stat_cor(method = "spearman", label.x.npc = "center", label.y.npc = "top", 
           cor.coef.name = "rho", size = 6) + 
  theme_bw() +
  labs(x = expression(paste("Chamber FCH"[4]," (nmol m"^-2, "s"^-1, ")")),
       y = expression(paste("EC FCH"[4]," (nmol m"^-2, "s"^-1, ")")),
       title = "Automatic") +
  theme(axis.text = element_text(size = 18),
        axis.title = element_text(size = 18),
        plot.title = element_text(size = 18, hjust = 0.5)) + 
  scale_y_continuous(limits = c(-2, 4)) + 
  scale_x_continuous(limits = c(-2, 4)) 
auto_EC_annual

# combine plots into one figure

automanualEC <- ggarrange(
  auto_EC + rremove("ylab") + rremove("xlab"),
  manual_EC + rremove("ylab") + rremove("xlab"), 
  auto_EC_weekly + rremove("ylab") + rremove("xlab"),
  manual_EC_weekly + rremove("ylab") + rremove("xlab"), 
  auto_EC_monthly + rremove("ylab") + rremove("xlab"),
  manual_EC_monthly + rremove("ylab") + rremove("xlab"), 
  auto_EC_annual + rremove("ylab") + rremove("xlab"),
  manual_EC_annual + rremove("ylab") + rremove("xlab"), 
  ncol = 2, 
  nrow = 4,
  align="hv"
)

automanualEC <- annotate_figure(
  automanualEC, 
  left = textGrob(expression(paste("EC FCH"[4], " nmol m"^-2, " s"^-1)), rot = 90, vjust = 1, hjust = 0.3, gp = gpar(fontsize = 20)),
  bottom = textGrob(expression(paste("Chamber FCH"[4], " nmol m"^-2, " s"^-1)), gp = gpar(fontsize = 20))
)

automanualEC

#####################################################################################################################################

## Linear mixed models for estimating the slopes of plot-scale FCH4 and ecosystem-scale FCH4 relationships

## use dataset versions with cnhgu included

hh_noNA <- all_hh_aggr %>%
  filter(!is.na(EC_FCH4) & !is.na(CH_FCH4_MEDIAN))
hr_noNA <- all_hr_aggr %>%
  filter(!is.na(EC_FCH4_median) & !is.na(CH_FCH4_MEDIAN))
d_noNA <- all_d_aggr %>%
  filter(!is.na(EC_FCH4_MEDIAN) & !is.na(CH_FCH4_MEDIAN))
week_noNA <- all_week_aggr %>%
  filter(!is.na(EC_FCH4_MEDIAN) & !is.na(CH_FCH4_MEDIAN))
month_noNA_noCNHGU <- all_month_aggr_noCNHGU %>%
  filter(!is.na(EC_FCH4_MEDIAN) & !is.na(CH_FCH4_MEDIAN))
yr_noNA <- all_yr_aggr %>%
  filter(!is.na(EC_FCH4_MEDIAN) & !is.na(CH_FCH4_MEDIAN))

# Fit linear mixed models with inverse hyperbolic sine tranformation (to meet residual normality)
# Response-only transform
m_y <- lme(asinh(EC_FCH4_MEDIAN) ~ CH_FCH4_MEDIAN, random = ~ 1 | SITE, data = d_noNA, method = "REML")

# Both-sides transform
m_both <- lme(asinh(EC_FCH4_MEDIAN) ~ asinh(CH_FCH4_MEDIAN), random = ~ 1 | SITE, data = d_noNA, method = "REML")

# Compare (non-nested: just use AIC/BIC and diagnostics)
anova(m_y, m_both)

m0 <- lme(EC_FCH4_MEDIAN ~ CH_FCH4_MEDIAN, random = ~ 1 | SITE, data = d_noNA, method = "REML")


# Diagnostics
op <- par(mfrow = c(2,2))
# m_y
res_y <- residuals(m_y, type = "normalized"); fit_y <- fitted(m_y)
qqnorm(res_y, main = "QQ: asinh(EC) ~ CH"); qqline(res_y)
plot(fit_y, res_y, xlab="Fitted", ylab="Norm. residuals",
     main="Residuals vs Fitted (m_y)"); abline(h=0, lty=2)

# m_both
res_b <- residuals(m_both, type = "normalized"); fit_b <- fitted(m_both)
qqnorm(res_b, main = "QQ: asinh(EC) ~ asinh(CH)"); qqline(res_b)
plot(fit_b, res_b, xlab="Fitted", ylab="Norm. residuals",
     main="Residuals vs Fitted (m_both)"); abline(h=0, lty=2)
par(op)

summary(m_y)
summary(m_both)

# m_both
res_0 <- residuals(m0, type = "normalized"); fit_0 <- fitted(m0)
qqnorm(res_0, main = "QQ: EC ~ CH"); qqline(res_0)
plot(fit_0, res_0, xlab="Fitted", ylab="Norm. residuals",
     main="Residuals vs Fitted (m0)"); abline(h=0, lty=2)

anova(m0, m_y, m_both)

# response only transformed seems best

# get estimate of EC FCH4 at median chamber FCH4

# extract fixed effects
b  <- fixef(m_y)

# get the chamber FCH4 reference value (median)
x0 <- median(d_noNA$CH_FCH4_MEDIAN, na.rm = TRUE)

# predicted EC FCH4 at median chamber FCH4
eta0 <- b[1] + b[2]*x0

# predicted EC FCH4 at median chamber FCH4 +1 nmol m-2 s-1
eta1 <- b[1] + b[2]*(x0 + 1)

# back-transform to original scale
# --> the change in EC FCH4 (original units) for a +1 increase in chamber FCH₄, evaluated at the median
sinh(eta1) - sinh(eta0)   

summary(m_y)

### AMEs (to get average effect across all observations)

# Data vector of chamber flux for this aggregation:
CH <- d_noNA$CH_FCH4_MEDIAN

# Fixed effects and their covariance
b <- fixef(m_y)
V <- vcov(m_y)

# AME using fixed effects only (population-average)
ame <- mean( b[2] * cosh(b[1] + b[2]*CH), na.rm = TRUE )

# 95% CI by parametric draws of (beta0, beta1)
set.seed(1)
B <- mvrnorm(5000, mu = b, Sigma = V)
ame_draws <- apply(B, 1, function(bb) {
  mean( bb[2] * cosh(bb[1] + bb[2]*CH), na.rm = TRUE )
})
ci <- quantile(ame_draws, c(.025, .975))

ame
ci

#### weekly

# Fit linear mixed models with inverse hyperbolic sine tranformation (to meet residual normality)
# Response-only transform
m_y <- lme(asinh(EC_FCH4_MEDIAN) ~ CH_FCH4_MEDIAN, random = ~ 1 | SITE, data = week_noNA, method = "REML")

# Both-sides transform
m_both <- lme(asinh(EC_FCH4_MEDIAN) ~ asinh(CH_FCH4_MEDIAN), random = ~ 1 | SITE, data = week_noNA, method = "REML")

# Compare
anova(m_y, m_both)

m0 <- lme(EC_FCH4_MEDIAN ~ CH_FCH4_MEDIAN, random = ~ 1 | SITE, data = week_noNA, method = "REML")

# Diagnostics
op <- par(mfrow = c(2,2))
# m_y
res_y <- residuals(m_y, type = "normalized"); fit_y <- fitted(m_y)
qqnorm(res_y, main = "QQ: asinh(EC) ~ CH"); qqline(res_y)
plot(fit_y, res_y, xlab="Fitted", ylab="Norm. residuals",
     main="Residuals vs Fitted (m_y)"); abline(h=0, lty=2)

# m_both
res_b <- residuals(m_both, type = "normalized"); fit_b <- fitted(m_both)
qqnorm(res_b, main = "QQ: asinh(EC) ~ asinh(CH)"); qqline(res_b)
plot(fit_b, res_b, xlab="Fitted", ylab="Norm. residuals",
     main="Residuals vs Fitted (m_both)"); abline(h=0, lty=2)
par(op)

summary(m_y) 
summary(m_both)

# m_both
res_0 <- residuals(m0, type = "normalized"); fit_0 <- fitted(m0)
qqnorm(res_0, main = "QQ: EC ~ CH"); qqline(res_0)
plot(fit_0, res_0, xlab="Fitted", ylab="Norm. residuals",
     main="Residuals vs Fitted (m0)"); abline(h=0, lty=2)

anova(m0, m_y, m_both)

### residuals look good for m_y 

summary(m_y)

# get estimate of EC FCH4 at median chamber FCH4

# extract fixed effects
b  <- fixef(m_y)  

# get the chamber FCH4 reference value (median)
x0 <- median(week_noNA$CH_FCH4_MEDIAN, na.rm = TRUE)

# predicted EC FCH4 at median chamber FCH4
eta0 <- b[1] + b[2]*x0

# predicted EC FCH4 at median chamber FCH4 +1 nmol m-2 s-1
eta1 <- b[1] + b[2]*(x0 + 1)

# back-transform to original scale
# --> the change in EC FCH4 (original units) for a +1 increase in chamber FCH₄, evaluated at the median
sinh(eta1) - sinh(eta0)

### AMEs (to get average effect across all observations)

# Data vector of chamber flux for this aggregation:
CH <- week_noNA$CH_FCH4_MEDIAN

# Fixed effects and their covariance
b <- fixef(m_y)
V <- vcov(m_y)

# AME using fixed effects only (population-average)
ame <- mean( b[2] * cosh(b[1] + b[2]*CH), na.rm = TRUE )

# 95% CI by parametric draws of (beta0, beta1)
set.seed(1)
B <- mvrnorm(5000, mu = b, Sigma = V)
ame_draws <- apply(B, 1, function(bb) {
  mean( bb[2] * cosh(bb[1] + bb[2]*CH), na.rm = TRUE )
})
ci <- quantile(ame_draws, c(.025, .975))

ame
ci

#### monthly

# Fit linear mixed models with inverse hyperbolic sine tranformation (to meet residual normality)
# Response-only transform
m_y <- lme(asinh(EC_FCH4_MEDIAN) ~ CH_FCH4_MEDIAN, random = ~ 1 | SITE, data = month_noNA, method = "REML")

# Both-sides transform
m_both <- lme(asinh(EC_FCH4_MEDIAN) ~ asinh(CH_FCH4_MEDIAN), random = ~ 1 | SITE, data = week_noNA, method = "REML")

# Compare 
anova(m_y, m_both)

m0 <- lme(EC_FCH4_MEDIAN ~ CH_FCH4_MEDIAN, random = ~ 1 | SITE, data = week_noNA, method = "REML")

# Diagnostics
op <- par(mfrow = c(2,2))
# m_y
res_y <- residuals(m_y, type = "normalized"); fit_y <- fitted(m_y)
qqnorm(res_y, main = "QQ: asinh(EC) ~ CH"); qqline(res_y)
plot(fit_y, res_y, xlab="Fitted", ylab="Norm. residuals",
     main="Residuals vs Fitted (m_y)"); abline(h=0, lty=2)

# m_both
res_b <- residuals(m_both, type = "normalized"); fit_b <- fitted(m_both)
qqnorm(res_b, main = "QQ: asinh(EC) ~ asinh(CH)"); qqline(res_b)
plot(fit_b, res_b, xlab="Fitted", ylab="Norm. residuals",
     main="Residuals vs Fitted (m_both)"); abline(h=0, lty=2)
par(op)

summary(m_y) 
summary(m_both)

# m_both
res_0 <- residuals(m0, type = "normalized"); fit_0 <- fitted(m0)
qqnorm(res_0, main = "QQ: EC ~ CH"); qqline(res_0)
plot(fit_0, res_0, xlab="Fitted", ylab="Norm. residuals",
     main="Residuals vs Fitted (m0)"); abline(h=0, lty=2)

anova(m0, m_y, m_both)

### residuals look good for asinh

summary(m_y)

# extract fixed effects
b  <- fixef(m_y) 

# get the chamber FCH4 reference value (median)
x0 <- median(week_noNA$CH_FCH4_MEDIAN, na.rm = TRUE)

# predicted EC FCH4 at median chamber FCH4
eta0 <- b[1] + b[2]*x0

# predicted EC FCH4 at median chamber FCH4 +1 nmol m-2 s-1
eta1 <- b[1] + b[2]*(x0 + 1)

# back-transform to original scale
# --> the change in EC FCH4 (original units) for a +1 increase in chamber FCH₄, evaluated at the median
sinh(eta1) - sinh(eta0) 

### AMEs (to get average effect across all observations)

# Data vector of chamber flux for this aggregation:
CH <- month_noNA$CH_FCH4_MEDIAN

# Fixed effects and their covariance
b <- fixef(m_y)
V <- vcov(m_y)

# AME using fixed effects only (population-average)
ame <- mean( b[2] * cosh(b[1] + b[2]*CH), na.rm = TRUE )

# 95% CI by parametric draws of (beta0, beta1)
set.seed(1)
B <- mvrnorm(5000, mu = b, Sigma = V)
ame_draws <- apply(B, 1, function(bb) {
  mean( bb[2] * cosh(bb[1] + bb[2]*CH), na.rm = TRUE )
})
ci <- quantile(ame_draws, c(.025, .975))

ame
ci

### annual

# Fit linear mixed models with inverse hyperbolic sine tranformation (to meet residual normality)
# Response-only transform
m_y <- lme(asinh(EC_FCH4_MEDIAN) ~ CH_FCH4_MEDIAN, random = ~ 1 | SITE, data = yr_noNA, method = "REML")

# Both-sides transform
m_both <- lme(asinh(EC_FCH4_MEDIAN) ~ asinh(CH_FCH4_MEDIAN), random = ~ 1 | SITE, data = month_noNA, method = "REML")

# Compare 
anova(m_y, m_both)

m0 <- lme(EC_FCH4_MEDIAN ~ CH_FCH4_MEDIAN, random = ~ 1 | SITE, data = yr_noNA, method = "REML")

# Diagnostics
op <- par(mfrow = c(2,2))
# m_y
res_y <- residuals(m_y, type = "normalized"); fit_y <- fitted(m_y)
qqnorm(res_y, main = "QQ: asinh(EC) ~ CH"); qqline(res_y)
plot(fit_y, res_y, xlab="Fitted", ylab="Norm. residuals",
     main="Residuals vs Fitted (m_y)"); abline(h=0, lty=2)

# m_both
res_b <- residuals(m_both, type = "normalized"); fit_b <- fitted(m_both)
qqnorm(res_b, main = "QQ: asinh(EC) ~ asinh(CH)"); qqline(res_b)
plot(fit_b, res_b, xlab="Fitted", ylab="Norm. residuals",
     main="Residuals vs Fitted (m_both)"); abline(h=0, lty=2)
par(op)

summary(m_y) 
summary(m_both)

# m_both
res_0 <- residuals(m0, type = "normalized"); fit_0 <- fitted(m0)
qqnorm(res_0, main = "QQ: EC ~ CH"); qqline(res_0)
plot(fit_0, res_0, xlab="Fitted", ylab="Norm. residuals",
     main="Residuals vs Fitted (m0)"); abline(h=0, lty=2)

anova(m0, m_y, m_both)

# extract fixed effects
b  <- fixef(m_y)

# get the chamber FCH4 reference value (median)
x0 <- median(yr_noNA$CH_FCH4_MEDIAN, na.rm = TRUE)

# predicted EC FCH4 at median chamber FCH4
eta0 <- b[1] + b[2]*x0

# predicted EC FCH4 at median chamber FCH4 +1 nmol m-2 s-1
eta1 <- b[1] + b[2]*(x0 + 1)

# back-transform to original scale
# --> the change in EC FCH4 (original units) for a +1 increase in chamber FCH₄, evaluated at the median
sinh(eta1) - sinh(eta0)

summary(m_y)

### AMEs (to get average effect across all observations)

# Data vector of chamber flux for this aggregation:
CH <- yr_noNA$CH_FCH4_MEDIAN

# Fixed effects and their covariance
b <- fixef(m_y)
V <- vcov(m_y)

# AME using fixed effects only (population-average)
ame <- mean( b[2] * cosh(b[1] + b[2]*CH), na.rm = TRUE )

# 95% CI by parametric draws of (beta0, beta1)
set.seed(1)
B <- mvrnorm(5000, mu = b, Sigma = V)
ame_draws <- apply(B, 1, function(bb) {
  mean( bb[2] * cosh(bb[1] + bb[2]*CH), na.rm = TRUE )
})
ci <- quantile(ame_draws, c(.025, .975))

ame
ci

####################################################################################################################################
### testing the effect of DOMINANT_VEGETATION in the models with linear regressions (see Methods 2.4.2)

# monthly

# with veg
m_month_lm <- lm(
  DELTA_FCH4_abs_yjtrans_ctr ~ SITE_TS_TOP_ctr +
    SITE_WTL_ctr +
    USTAR_MEDIAN_ctr +
    U_WIND_F_MEAN_ctr +
    V_WIND_F_MEAN_ctr +
    PA_F_MEDIAN_ctr +
    NEE_F_MEAN_ctr +
    VPD_F_MEDIAN_ctr +
    DOMINANT_VEGETATION + 
    MONTH, 
  data = month_noNA_noCNHGU
)
summary(m_month_lm)
# R2 adj: 0.7535

# with site instead of veg
m_month_lm <- lm(
  DELTA_FCH4_abs_yjtrans_ctr ~ SITE_TS_TOP_ctr +
    SITE_WTL_ctr +
    USTAR_MEDIAN_ctr +
    U_WIND_F_MEAN_ctr +
    V_WIND_F_MEAN_ctr +
    PA_F_MEDIAN_ctr +
    NEE_F_MEAN_ctr +
    VPD_F_MEDIAN_ctr +
    SITE + 
    MONTH, 
  data = month_noNA_noCNHGU
)
summary(m_month_lm)
# R2 adj: 0.8712

# weekly

# with veg
m_week_lm <- lm(
  DELTA_FCH4_abs_yjtrans_ctr ~ SITE_TS_TOP_ctr +
    SITE_WTL_ctr +
    USTAR_MEDIAN_ctr +
    U_WIND_F_MEAN_ctr +
    V_WIND_F_MEAN_ctr +
    PA_F_MEDIAN_ctr +
    NEE_F_MEAN_ctr +
    VPD_F_MEDIAN_ctr +
    DOMINANT_VEGETATION + 
    MONTH,week_noNA_noCNHGU
)
summary(m_week_lm)
# R2 adj: 0.6422

# without veg, site instead
m_week_lm <- lm(
  DELTA_FCH4_abs_yjtrans ~ SITE_TS_TOP_ctr +
    SITE_WTL_ctr +
    USTAR_MEDIAN_ctr +
    U_WIND_F_MEAN_ctr +
    V_WIND_F_MEAN_ctr +
    PA_F_MEDIAN_ctr +
    NEE_F_MEAN_ctr +
    VPD_F_MEDIAN_ctr +
    SITE + 
    MONTH, 
  data = week_noNA_noCNHGU
)
summary(m_week_lm)
# R2 adj: 0.7604

# daily

# with veg, no site
m_day_lm <- lm(
  DELTA_FCH4_abs_yjtrans_ctr ~ SITE_TS_TOP_ctr +
    SITE_WTL_ctr +
    USTAR_MEDIAN_ctr +
    U_WIND_F_MEAN_ctr +
    V_WIND_F_MEAN_ctr +
    PA_F_MEDIAN_ctr +
    NEE_F_MEAN_ctr +
    VPD_F_MEDIAN_ctr +
    DOMINANT_VEGETATION + 
    MONTH,
  d_noNA_noCNHGU
)
summary(m_day_lm)
# R2 adj: 0.588

# without veg, with site
m_day_lm <- lm(
  DELTA_FCH4_abs_yjtrans_ctr ~ SITE_TS_TOP_ctr +
    SITE_WTL_ctr +
    USTAR_MEDIAN_ctr +
    U_WIND_F_MEAN_ctr +
    V_WIND_F_MEAN_ctr +
    PA_F_MEDIAN_ctr +
    NEE_F_MEAN_ctr +
    VPD_F_MEDIAN_ctr +
    SITE + 
    MONTH, 
  data = d_noNA_noCNHGU
)
summary(m_day_lm)
# R2 adj: 0.7866

# hourly

# with veg, no site
m_hr_lm <- lm(
  DELTA_FCH4_abs_yjtrans_ctr ~ SITE_TS_TOP_ctr +
    SITE_WTL_ctr +
    USTAR_MEDIAN_ctr +
    U_WIND_F_MEAN_ctr +
    V_WIND_F_MEAN_ctr +
    PA_F_MEDIAN_ctr +
    NEE_F_mean_ctr +
    VPD_F_MEDIAN_ctr +
    DOMINANT_VEGETATION + 
    HOUR +
    MONTH,
  hr_noNA_noCNHGU
)
summary(m_hr_lm)
# R2 adj: 0.4611

# without veg, with site
m_hr_lm <- lm(
  DELTA_FCH4_abs_yjtrans_ctr ~ SITE_TS_TOP_ctr +
    SITE_WTL_ctr +
    USTAR_MEDIAN_ctr +
    U_WIND_F_MEAN_ctr +
    V_WIND_F_MEAN_ctr +
    PA_F_MEDIAN_ctr +
    NEE_F_mean_ctr +
    VPD_F_MEDIAN_ctr +
    SITE + 
    HOUR +
    MONTH, 
  data = hr_noNA_noCNHGU
)
summary(m_hr_lm)
# R2 adj: 0.7051

# half-hourly

# with veg, no site
m_hh_lm <- lm(
  DELTA_FCH4_abs_yjtrans_ctr ~ SITE_TS_TOP_ctr +
    SITE_WTL_ctr +
    EC_USTAR_ctr +
    U_WIND_F_ctr +
    V_WIND_F_ctr +
    EC_PA_F_ctr +
    NEE_F_ctr +
    EC_VPD_F_ctr +
    DOMINANT_VEGETATION + 
    HOUR +
    MONTH,
  hh_noNA_noCNHGU
)
summary(m_hh_lm)
# R2 adj: 0.4605

# without veg, with site
m_hh_lm <- lm(
  DELTA_FCH4_abs_yjtrans_ctr ~ SITE_TS_TOP_ctr +
    SITE_WTL_ctr +
    EC_USTAR_ctr +
    U_WIND_F_ctr +
    V_WIND_F_ctr +
    EC_PA_F_ctr +
    NEE_F_ctr +
    EC_VPD_F_ctr +
    SITE + 
    HOUR +
    MONTH,
  hh_noNA_noCNHGU
)
summary(m_hh_lm)
# R2 adj: 0.7031




##########################################################################################################################################
## DAILY, WEEKLY, MONTHLY, ANNUAL SUMS FOR EC AND CHAMBER FLUXES

all_sites_nodupl <- read.csv("path/allsites_raw_no_duplicates.csv")

# need to first create a chamber-aggregated df to get to EC timestamp level
# FI-SI2, US-LA1 and US-LA2 have only DATE column

ch_aggr <- all_sites_nodupl %>%
  group_by(SITE, 
           time_var = case_when(
             SITE %in% c("FI-SI2", "US-LA1", "US-LA2") ~ as_date(DATE), 
             TRUE ~ as_datetime(TIMESTAMP_START)
           )) %>%
  summarise(
    ch_FCH4_mean   = mean(ch_FCH4_nmolCH4m2s1, na.rm = TRUE),
    ch_FCH4_median = median(ch_FCH4_nmolCH4m2s1, na.rm = TRUE),
    EC_FCH4_comb   = first(EC_FCH4_comb),
    EC_FCH4_F_ANNOPTLM = first(EC_FCH4_F_ANNOPTLM),
    .groups = "drop"
  ) %>%
  rename(time_ref = time_var)

# create the sums

### need to do this differently for each site because need to convert from nmol rate to mg CH4 m-2
### and this is done differently for n=7 sites than in the rest because of lack of time info for n=3 sites

# constants
# mw_CH4_g_per_mol <- 16.04        # molar mass of CH4 [g/mol]
# nmol_to_mol      <- 1e-9         # 1 nmol = 1e-9 mol
# g_to_mg          <- 1000         # 1 g = 1000 mg
# 
# sec_per_halfhour <- 1800         # 30 minutes
# sec_per_day      <- 86400        # 24 hours

# Conversion: (nmol m^-2 s^-1) * seconds -> mg CH4 m^-2
# mg per (nmol m^-2 s^-1) per half-hour window:
k_hh  <- 1e-9 * 16.04 * 1000 * 1800
# mg per (nmol m^-2 s^-1) per day:
k_day <- 1e-9 * 16.04 * 1000 * 86400


# Half-hour contributions in mg CH4 m^-2 (only sites with half-hour timestamps)
hh_v1 <- ch_aggr %>%
  filter(!(SITE %in% c("FI-SI2","US-LA1","US-LA2"))) %>%
  mutate(
    ch_mean_mgCH4m2   = ch_FCH4_mean   * k_hh,
    ch_median_mgCH4m2 = ch_FCH4_median * k_hh,
    ec_comb_mgCH4m2   = EC_FCH4_comb   * k_hh,
    ec_ann_mgCH4m2    = EC_FCH4_F_ANNOPTLM * k_hh,
    date = as_date(time_ref),
    iso_year = isoyear(time_ref),
    iso_week = isoweek(time_ref),
    year = year(time_ref),
    month = month(time_ref)
  )

# Daily / ISO-week / Monthly / Annual cumulative flux (mg CH4 m^-2)
sum_daily_v1 <- hh_v1 %>%
  group_by(SITE, date) %>%
  summarise(across(ends_with("mgCH4m2"), ~sum(.x, na.rm = TRUE)),
            n_halfhours = n(), .groups = "drop")

sum_isoweek_v1 <- hh_v1 %>%
  group_by(SITE, iso_year, iso_week) %>%
  summarise(across(ends_with("mgCH4m2"), ~sum(.x, na.rm = TRUE)),
            n_halfhours = n(), .groups = "drop")

sum_monthly_v1 <- hh_v1 %>%
  group_by(SITE, year, month) %>%
  summarise(across(ends_with("mgCH4m2"), ~sum(.x, na.rm = TRUE)),
            n_halfhours = n(), .groups = "drop")

sum_annual_v1 <- hh_v1 %>%
  group_by(SITE, year) %>%
  summarise(across(ends_with("mgCH4m2"), ~sum(.x, na.rm = TRUE)),
            n_halfhours = n(), .groups = "drop")


# FI-SI2, US-LA1, US-LA2
## (using EC data sets from Maatta_et_al_2026_EC_data_cleaning.R script)

FISI2EC_ECfilt <- read.csv("path/FISI2EC_ECfilt.csv")
USLA1EC_ECfilt <- read.csv("path/USLA1EC_ECfilt.csv")
USLA2EC_ECfilt <- read.csv("path/USLA2EC_ECfilt.csv")

sum_ec_daily_fisi2 <- FISI2EC_ECfilt %>%
  mutate(
    ec_mgCH4 = EC_FCH4_F_ANNOPTLM * k_hh,
    date = as.Date(EC_TIMESTAMP_START)
  ) %>%
  filter(is.finite(EC_FCH4_F_ANNOPTLM)) %>%
  group_by(SITE, date) %>%
  summarise(
    ec_daily_mgCH4 = sum(ec_mgCH4, na.rm = TRUE),
    n_windows      = n(),
    .groups = "drop"
  )

sum_ec_daily_usla1 <- USLA1EC_ECfilt %>%
  mutate(
    ec_mgCH4 = EC_FCH4_F_ANNOPTLM * k_hh,
    date = as.Date(EC_TIMESTAMP_START)
  ) %>%
  filter(is.finite(EC_FCH4_F_ANNOPTLM)) %>%
  group_by(SITE, date) %>%
  summarise(
    ec_daily_mgCH4 = sum(ec_mgCH4, na.rm = TRUE),
    n_windows      = n(),
    .groups = "drop"
  )


sum_ec_daily_usla2 <- USLA2EC_ECfilt %>%
  mutate(
    ec_mgCH4 = EC_FCH4_F_ANNOPTLM * k_hh,
    date = as.Date(EC_TIMESTAMP_START)
  ) %>%
  filter(is.finite(EC_FCH4_F_ANNOPTLM)) %>%
  group_by(SITE, date) %>%
  summarise(
    ec_daily_mgCH4 = sum(ec_mgCH4, na.rm = TRUE),
    n_windows      = n(),
    .groups = "drop"
  )

# Non-DATE-only sites: per half-hour -> daily totals
hh_nondate <- ch_aggr %>%
  filter(!(SITE %in% c("FI-SI2","US-LA1","US-LA2"))) %>%
  mutate(
    # contributions per half-hour window in mg CH4 m^-2
    ch_median_mgCH4m2 = ch_FCH4_median * k_hh,
    ch_mean_mgCH4m2   = ch_FCH4_mean   * k_hh,
    ec_comb_mgCH4m2   = EC_FCH4_comb   * k_hh,
    ec_F_ANNOPTLM_mgCH4m2    = EC_FCH4_F_ANNOPTLM * k_hh,
    date = as.Date(time_ref)
  ) %>%
  dplyr::select(SITE, date, time_ref, ch_median_mgCH4m2, ch_mean_mgCH4m2, ec_comb_mgCH4m2, ec_F_ANNOPTLM_mgCH4m2)

daily_nondate <- hh_nondate %>%
  group_by(SITE, date) %>%
  summarise(
    ch_median_mgCH4m2 = sum(ch_median_mgCH4m2, na.rm = TRUE),
    ch_mean_mgCH4m2   = sum(ch_mean_mgCH4m2,   na.rm = TRUE),
    ec_comb_mgCH4m2   = sum(ec_comb_mgCH4m2,   na.rm = TRUE),
    ec_F_ANNOPTLM_mgCH4m2    = sum(ec_F_ANNOPTLM_mgCH4m2,    na.rm = TRUE),
    .groups = "drop"
  )

# DATE-only chamber sites
# Chamber daily medians/means across collars (rates) -> daily totals via k_day
ch_daily_dateonly <- all_sites_nodupl %>%
  filter(SITE %in% c("FI-SI2","US-LA1","US-LA2")) %>%
  group_by(SITE, DATE = as.Date(DATE)) %>%
  summarise(
    ch_daily_median_rate = median(ch_FCH4_nmolCH4m2s1, na.rm = TRUE),
    ch_daily_mean_rate   = mean(ch_FCH4_nmolCH4m2s1, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    ch_median_mgCH4m2 = ch_daily_median_rate * k_day,  # median-as-constant over the full day
    ch_mean_mgCH4m2   = ch_daily_mean_rate   * k_day    # mean-as-constant over the full day
  ) %>%
  dplyr::select(SITE, date = DATE, ch_median_mgCH4m2, ch_mean_mgCH4m2)

# EC daily totals for those same dates

# combine fisi2, usla1, usla2
ec_daily_dateonly <- bind_rows(sum_ec_daily_fisi2, sum_ec_daily_usla1, sum_ec_daily_usla2)
ec_daily_dateonly <- as.data.frame(ec_daily_dateonly)

ec_daily_dateonly <- subset(ec_daily_dateonly, select = -n_windows)

# Join chamber + EC for DATE-only sites
daily_dateonly <- ch_daily_dateonly %>%
  full_join(ec_daily_dateonly, by = c("SITE","date"))


# Combine all sites (daily mg CH4 m^-2)
daily_nondate <- subset(daily_nondate, select = -ec_F_ANNOPTLM_mgCH4m2)

#rename the ec column
colnames(daily_nondate)[5] <- "ec_daily_mgCH4"

sum_daily_v2 <- bind_rows(daily_nondate, daily_dateonly)

sum_daily_v2 <- sum_daily_v2 %>% drop_na(ec_daily_mgCH4)

sum_daily_v2 <- as.data.frame(sum_daily_v2)

# Aggregate daily to ISO-week / Month / Year
sum_isoweek_v2 <- sum_daily_v2 %>%
  mutate(iso_year = isoyear(date), iso_week = isoweek(date)) %>%
  group_by(SITE, iso_year, iso_week) %>%
  summarise(
    ch_median_mgCH4m2 = sum(ch_median_mgCH4m2, na.rm = TRUE),
    ch_mean_mgCH4m2   = sum(ch_mean_mgCH4m2,   na.rm = TRUE),
    ec_comb_mgCH4m2   = sum(ec_daily_mgCH4,   na.rm = TRUE),
    .groups = "drop"
  )

sum_isoweek_v2 <- as.data.frame(sum_isoweek_v2)

sum_monthly_v2 <- sum_daily_v2 %>%
  mutate(year = year(date), month = month(date)) %>%
  group_by(SITE, year, month) %>%
  summarise(
    ch_median_mgCH4m2 = sum(ch_median_mgCH4m2, na.rm = TRUE),
    ch_mean_mgCH4m2   = sum(ch_mean_mgCH4m2,   na.rm = TRUE),
    ec_comb_mgCH4m2   = sum(ec_daily_mgCH4,   na.rm = TRUE),
    .groups = "drop"
  )

sum_monthly_v2 <- as.data.frame(sum_monthly_v2)

sum_annual_v2 <- sum_daily_v2 %>%
  mutate(year = year(date)) %>%
  group_by(SITE, year) %>%
  summarise(
    ch_median_mgCH4m2 = sum(ch_median_mgCH4m2, na.rm = TRUE),
    ch_mean_mgCH4m2   = sum(ch_mean_mgCH4m2,   na.rm = TRUE),
    ec_comb_mgCH4m2   = sum(ec_daily_mgCH4,   na.rm = TRUE),
    .groups = "drop"
  )

sum_annual_v2 <- as.data.frame(sum_annual_v2)

# calculate sum difference

sum_daily_v2$ECCH_diff_chmedian <- sum_daily_v2$ec_daily_mgCH4 - sum_daily_v2$ch_median_mgCH4m2
sum_isoweek_v2$ECCH_diff_chmedian <- sum_isoweek_v2$ec_comb_mgCH4m2 - sum_isoweek_v2$ch_median_mgCH4m2
sum_monthly_v2$ECCH_diff_chmedian <- sum_monthly_v2$ec_comb_mgCH4m2 - sum_monthly_v2$ch_median_mgCH4m2
sum_annual_v2$ECCH_diff_chmedian <- sum_annual_v2$ec_comb_mgCH4m2 - sum_annual_v2$ch_median_mgCH4m2

sum_daily_v2$ECCH_diff_chmean <- sum_daily_v2$ec_daily_mgCH4 - sum_daily_v2$ch_mean_mgCH4m2
sum_isoweek_v2$ECCH_diff_chmean <- sum_isoweek_v2$ec_comb_mgCH4m2 - sum_isoweek_v2$ch_mean_mgCH4m2
sum_monthly_v2$ECCH_diff_chmean <- sum_monthly_v2$ec_comb_mgCH4m2 - sum_monthly_v2$ch_mean_mgCH4m2
sum_annual_v2$ECCH_diff_chmean <- sum_annual_v2$ec_comb_mgCH4m2 - sum_annual_v2$ch_mean_mgCH4m2

## Cumulative per SITE-YEAR (resets each year)
sum_daily <- sum_daily_v2 %>%
  arrange(SITE, date) %>%
  mutate(year = year(date)) %>%
  group_by(SITE, year) %>%
  mutate(
    daily_sum_ec        = ec_daily_mgCH4,
    daily_sum_ch_median = ch_median_mgCH4m2,
    daily_sum_ch_mean   = ch_mean_mgCH4m2,
    cumulative_sum_ec          = cumsum(coalesce(ec_daily_mgCH4,    0)),
    cumulative_sum_ch_median   = cumsum(coalesce(ch_median_mgCH4m2, 0)),
    cumulative_sum_ch_mean     = cumsum(coalesce(ch_mean_mgCH4m2,   0))
  ) %>%
  ungroup()

# Long format with separate 'type' and 'method' (no collapsing) for ggplot
plot_df <- sum_daily %>%
  pivot_longer(
    cols = c(daily_sum_ec, daily_sum_ch_median, daily_sum_ch_mean,
             cumulative_sum_ec, cumulative_sum_ch_median, cumulative_sum_ch_mean),
    names_to = "series", values_to = "sum_mgCH4_m2"
  ) %>%
  mutate(
    type = ifelse(grepl("^daily_", series), "Daily", "Cumulative"),
    method = dplyr::recode(series,
                           "daily_sum_ec"             = "EC",
                           "daily_sum_ch_median"      = "Chamber median",
                           "daily_sum_ch_mean"        = "Chamber mean",
                           "cumulative_sum_ec"        = "EC",
                           "cumulative_sum_ch_median" = "Chamber median",
                           "cumulative_sum_ch_mean"   = "Chamber mean"
    ),
    method = factor(method, levels = c("EC","Chamber median","Chamber mean"))
  ) %>%
  distinct(SITE, date, type, method, .keep_all = TRUE)   # just in case


# plot only cumulative (Fig. B16)
plot_cum <- plot_df %>%
  filter(type == "Cumulative")

ggplot(plot_cum, aes(date, sum_mgCH4_m2, color = method, shape = method, fill = method)) +
  geom_point(size = 2.2, stroke = 0.9, alpha=0.7) +
  facet_wrap(~ SITE, scales = "free") +
  scale_x_date(
    breaks = date_breaks("1 month"), 
    labels = date_format("%Y-%m")
  ) +
  scale_shape_manual(values = c("EC"=17, "Chamber median"=19, "Chamber mean"=21)) +
  scale_color_manual(values = c("EC"="#6497bf","Chamber median"="#d8031c","Chamber mean"="black")) +
  scale_fill_manual(values  = c("EC"=NA, "Chamber median"="#d8031c", "Chamber mean"="white"), guide="none") +
  labs(x="Date", y="Cumulative CH4 (mg CH4 m-2)", color=NULL, shape=NULL) +
  theme_bw() +
  theme(panel.grid.minor = element_blank(), 
        legend.position = "bottom",
        strip.background = element_blank(),
        strip.text = element_text(size = 16),
        axis.text = element_text(size = 16),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 16),
        axis.title = element_text(size = 16),
        legend.text = element_text(size = 16))

### Table C5 ###
## without fisi2, usla1 and usla2

## DAILY

daily_stats <- sum_daily_v2 %>%
  filter(!(SITE %in% c("FI-SI2","US-LA1","US-LA2"))) %>%
  summarise_at(c("ECCH_diff_chmedian", "ECCH_diff_chmean"),
               list(mean = ~mean(.x, na.rm = TRUE),
                    sd = ~sd(.x, na.rm = TRUE),
                    median = ~median(.x, na.rm = TRUE),
                    IQR = ~IQR(.x, na.rm = TRUE),
                    CV = ~raster::cv(.x, na.rm = TRUE)))

daily_stats <- as.data.frame(daily_stats)

daily_stats

## WEEKLY 

weekly_stats <- sum_isoweek_v2 %>%
  filter(!(SITE %in% c("FI-SI2","US-LA1","US-LA2"))) %>%
  summarise_at(c("ECCH_diff_chmedian", "ECCH_diff_chmean"),
               list(mean = ~mean(.x, na.rm = TRUE),
                    sd = ~sd(.x, na.rm = TRUE),
                    median = ~median(.x, na.rm = TRUE),
                    IQR = ~IQR(.x, na.rm = TRUE),
                    CV = ~raster::cv(.x, na.rm = TRUE)))

weekly_stats <- as.data.frame(weekly_stats)

weekly_stats

## MONTHLY 

monthly_stats <- sum_monthly_v2 %>%
  filter(!(SITE %in% c("FI-SI2","US-LA1","US-LA2"))) %>%
  summarise_at(c("ECCH_diff_chmedian", "ECCH_diff_chmean"),
               list(mean = ~mean(.x, na.rm = TRUE),
                    sd = ~sd(.x, na.rm = TRUE),
                    median = ~median(.x, na.rm = TRUE),
                    IQR = ~IQR(.x, na.rm = TRUE),
                    CV = ~raster::cv(.x, na.rm = TRUE)))

monthly_stats <- as.data.frame(monthly_stats)

monthly_stats

## ANNUAL 

annual_stats <- sum_annual_v2 %>%
  filter(!(SITE %in% c("FI-SI2","US-LA1","US-LA2"))) %>%
  summarise_at(c("ECCH_diff_chmedian", "ECCH_diff_chmean"),
               list(mean = ~mean(.x, na.rm = TRUE),
                    sd = ~sd(.x, na.rm = TRUE),
                    median = ~median(.x, na.rm = TRUE),
                    IQR = ~IQR(.x, na.rm = TRUE),
                    CV = ~raster::cv(.x, na.rm = TRUE)))

annual_stats <- as.data.frame(annual_stats)

annual_stats

# Wilcoxon tests

# DAILY

daily_filtered <- sum_daily_v2 %>% filter(!(SITE %in% c("FI-SI2", "US-LA1", "US-LA2")))

# median-based (FI-SI2, US-LA1 and US-LA2 filtered out)
daily_ch <- subset(daily_filtered, select = c(date, ch_median_mgCH4m2))
daily_ch$method <- "chamber"
daily_ec <- subset(daily, select = c(date, ec_daily_mgCH4))
daily_ec$method <- "EC"
colnames(daily_ch)[2] <- "flux"
colnames(daily_ec)[2] <- "flux"
daily_comb <- rbind(daily_ch, daily_ec)
daily_comb <- daily_comb %>% drop_na(flux)

wilcox.test(flux ~ method, data = daily_comb) # p-value < 2.2e-16
daily_comb %>%
  filter(!is.na(flux) & method %in% c("EC", "chamber")) %>%
  group_by(method) %>%
  summarise(count = n())

# mean-based (FI-SI2, US-LA1 and US-LA2 filtered out)
daily_ch <- subset(daily_filtered, select = c(date, ch_mean_mgCH4m2))
daily_ch$method <- "chamber"
daily_ec <- subset(daily_filtered, select = c(date, ec_daily_mgCH4))
daily_ec$method <- "EC"
colnames(daily_ch)[2] <- "flux"
colnames(daily_ec)[2] <- "flux"
daily_comb <- rbind(daily_ch, daily_ec)
daily_comb <- daily_comb %>% drop_na(flux)

wilcox.test(flux ~ method, data = daily_comb) # p-value < 2.2e-16
daily_comb %>%
  filter(!is.na(flux) & method %in% c("EC", "chamber")) %>%
  group_by(method) %>%
  summarise(count = n())


# WEEKLY

weekly_filtered <- sum_isoweek_v2 %>% filter(!(SITE %in% c("FI-SI2", "US-LA1", "US-LA2")))

# median-based (FI-SI2, US-LA1 and US-LA2 filtered out)
weekly_ch <- subset(weekly_filtered, select = c(iso_year, iso_week, ch_median_mgCH4m2))
weekly_ch$method <- "chamber"
weekly_ec <- subset(weekly_filtered, select = c(iso_year, iso_week, ec_comb_mgCH4m2))
weekly_ec$method <- "EC"
colnames(weekly_ch)[3] <- "flux"
colnames(weekly_ec)[3] <- "flux"
weekly_comb <- rbind(weekly_ch, weekly_ec)
weekly_comb <- weekly_comb %>% drop_na(flux)

wilcox.test(flux ~ method, data = weekly_comb) # p-value < 2.2e-16
weekly_comb %>%
  filter(!is.na(flux) & method %in% c("EC", "chamber")) %>%
  group_by(method) %>%
  summarise(count = n())

# mean-based (FI-SI2, US-LA1 and US-LA2 filtered out)
weekly_ch <- subset(weekly, select = c(iso_year, iso_week, ch_mean_mgCH4m2))
weekly_ch$method <- "chamber"
weekly_ec <- subset(weekly, select = c(iso_year, iso_week, ec_comb_mgCH4m2))
weekly_ec$method <- "EC"
colnames(weekly_ch)[3] <- "flux"
colnames(weekly_ec)[3] <- "flux"
weekly_comb <- rbind(weekly_ch, weekly_ec)
weekly_comb <- weekly_comb %>% drop_na(flux)

wilcox.test(flux ~ method, data = weekly_comb) # p-value < 2.2e-16
weekly_comb %>%
  filter(!is.na(flux) & method %in% c("EC", "chamber")) %>%
  group_by(method) %>%
  summarise(count = n())


# MONTHLY

monthly_filtered <- sum_monthly_v2 %>% filter(!(SITE %in% c("FI-SI2", "US-LA1", "US-LA2")))

# median-based (FI-SI2, US-LA1 and US-LA2 filtered out)
m_ch <- subset(monthly_filtered, select = c(year, month, ch_median_mgCH4m2))
m_ch$method <- "chamber"
m_ec <- subset(monthly_filtered, select = c(year, month, ec_comb_mgCH4m2))
m_ec$method <- "EC"
colnames(m_ch)[3] <- "flux"
colnames(m_ec)[3] <- "flux"
m_comb <- rbind(m_ch, m_ec)
m_comb <- m_comb %>% drop_na(flux)

wilcox.test(flux ~ method, data = m_comb) # p-value < 2.2e-16
m_comb %>%
  filter(!is.na(flux) & method %in% c("EC", "chamber")) %>%
  group_by(method) %>%
  summarise(count = n())

# mean-based (FI-SI2, US-LA1 and US-LA2 filtered out)
m_ch <- subset(monthly_filtered, select = c(year, month, ch_mean_mgCH4m2))
m_ch$method <- "chamber"
m_ec <- subset(monthly_filtered, select = c(year, month, ec_comb_mgCH4m2))
m_ec$method <- "EC"
colnames(m_ch)[3] <- "flux"
colnames(m_ec)[3] <- "flux"
m_comb <- rbind(m_ch, m_ec)
m_comb <- m_comb %>% drop_na(flux)

wilcox.test(flux ~ method, data = m_comb) # p-value < 2.2e-16
weekly_comb %>%
  filter(!is.na(flux) & method %in% c("EC", "chamber")) %>%
  group_by(method) %>%
  summarise(count = n())

# ANNUAL

annual_filtered <- sum_annual_v2 %>% filter(!(SITE %in% c("FI-SI2", "US-LA1", "US-LA2")))

# # median-based (FI-SI2, US-LA1 and US-LA2 filtered out)
y_ch <- subset(annual_filtered, select = c(year, ch_median_mgCH4m2))
y_ch$method <- "chamber"
y_ec <- subset(annual_filtered, select = c(year, ec_comb_mgCH4m2))
y_ec$method <- "EC"
colnames(y_ch)[2] <- "flux"
colnames(y_ec)[2] <- "flux"
y_comb <- rbind(y_ch, y_ec)
y_comb <- y_comb %>% drop_na(flux)

wilcox.test(flux ~ method, data = y_comb) # p-value < 2.2e-16
y_comb %>%
  filter(!is.na(flux) & method %in% c("EC", "chamber")) %>%
  group_by(method) %>%
  summarise(count = n())

# # mean-based (FI-SI2, US-LA1 and US-LA2 filtered out)
y_ch <- subset(annual_filtered, select = c(year, ch_mean_mgCH4m2))
y_ch$method <- "chamber"
y_ec <- subset(annual_filtered, select = c(year, ec_comb_mgCH4m2))
y_ec$method <- "EC"
colnames(y_ch)[2] <- "flux"
colnames(y_ec)[2] <- "flux"
y_comb <- rbind(y_ch, y_ec)
y_comb <- y_comb %>% drop_na(flux)

wilcox.test(flux ~ method, data = y_comb) # p-value < 2.2e-16
y_comb %>%
  filter(!is.na(flux) & method %in% c("EC", "chamber")) %>%
  group_by(method) %>%
  summarise(count = n())



