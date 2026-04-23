#############################################################
### Code for: 
### Määttä et al. (2026) A cross-site comparison of ecosystem- and plot-scale methane fluxes across multiple sites
### Code created by Tiia Määttä, with parts written with GPT 4, 4o and 5.4

#############################################################
###--------------EC CH4 FLUX DATA CLEANING----------------###
#############################################################

# Remove all variables
rm(list=ls(all = TRUE))
# Clear the console window
cat("\014")
# Close all graphics windows
graphics.off()

# open libraries
library(dplyr)
library(lubridate)
library(ggplot2)
library(tidyr)

### NOTE: EC FCH4 datasets (except US-StJ) were obtained from FLUXNET-CH4 database (Delwiche et al. 2021, Knox et al., 2019) ###
### https://fluxnet.org/data/fluxnet-ch4-community-product/ ###

### CN-HGU ###

# read csv

CNHGU_EC <- read.csv("path/FLX_CN-Hgu_FLUXNET-CH4_HH_2015-2017_1-1.csv")

#create new column for Time2 with seconds
CNHGU_EC$Time2 <- paste(as.character(CNHGU_EC$TIMESTAMP_START),"00",sep="")

# convert to datetime

CNHGU_EC$TIMESTAMP_START <- as_datetime(CNHGU_EC$Time2)

#create new column for Time2 with seconds
CNHGU_EC$TimeEnd <- paste(as.character(CNHGU_EC$TIMESTAMP_END),"00",sep="")

# convert to datetime

CNHGU_EC$TIMESTAMP_END <- as_datetime(CNHGU_EC$TimeEnd)

# subset EC data to match with chamber data
CNHGU_EC_sub <- subset(CNHGU_EC, date(TIMESTAMP_START) >= "2015-08-01" & TIMESTAMP_END <= "2016-08-01 00:00:00")

# reset row names
rownames(CNHGU_EC_sub) <- NULL

# remove some columns that are not needed

CNHGU_EC_sub <- subset(CNHGU_EC_sub, select=-c(Time2, TimeEnd))

# replace -9999 with NA

CNHGU_EC_sub <- CNHGU_EC_sub %>%
  mutate(across(where(is.numeric), ~na_if(., -9999)))

# add EC_prefix to all columns
CNHGU_EC_sub <- CNHGU_EC_sub %>% rename_with( ~ paste0("EC_", .x))

# add Site

CNHGU_EC_sub$SITE <- "CN-HGU"

# move to front
CNHGU_EC_sub <- CNHGU_EC_sub %>%
  select(SITE, everything())

# save as .csv

write.csv(CNHGU_EC_sub, "path/CN_HGU_EC_subset.csv", row.names = FALSE)
CNHGU_EC_sub <- read.csv("path/CN_HGU_EC_subset.csv")

# calculate mean annual temperature and precipitation sum

# Summarize total annual precipitation
CNHGU_EC_sub$EC_TIMESTAMP_START <- as_datetime(CNHGU_EC_sub$EC_TIMESTAMP_START)

annual_precip <- CNHGU_EC_sub %>%
  group_by(year(EC_TIMESTAMP_START)) %>%
  summarize(Annual_Precip = sum(EC_P))

# Calculate the mean annual precipitation
mean(annual_precip$Annual_Precip)

# Calculate the annual mean temperature for each year
annual_mean_temp <- CNHGU_EC_sub %>%
  group_by(year(EC_TIMESTAMP_START)) %>%
  summarize(Annual_Mean_Temp = mean(EC_TA_F, na.rm = TRUE))

# Calculate the overall mean annual temperature across all years
mean(annual_mean_temp$Annual_Mean_Temp)

###### Filter out anomalous FCH4 data points with USTAR, empirical limits and precipitation filtering (see manuscript for details)

##### u*, -100 and dew point at monthly 99.5 percentile

# Dew point (Magnus; Alduchov & Eskridge constants)
dewpoint_c <- function(TC, RH){
  A <- 17.625; B <- 243.04
  RHc <- pmax(pmin(RH, 100), 1e-6)
  gamma <- log(RHc/100) + (A*TC)/(B + TC)
  (B * gamma) / (A - gamma)
}

adjacency_min <- 30  # exact 30-min adjacency

# get the thresholds
df0 <- CNHGU_EC_sub %>%
  arrange(EC_TIMESTAMP_START) %>%
  mutate(
    is_day     = EC_SW_IN_F > 10,
    is_night   = !is_day,
    Td         = dewpoint_c(EC_TA, EC_RH),
    near_sat   = !is.na(Td) & ((EC_TA - Td) <= 1),
    year_month = format(EC_TIMESTAMP_START, "%Y-%m")
  )

# remove extreme uptake (< -100) and apply nighttime u* filter
baseline <- df0 %>%
  filter(is.na(EC_FCH4) | EC_FCH4 >= -100) %>%
  mutate(ok_ustar = !(is_night & EC_USTAR < 0.1)) %>%
  filter(ok_ustar)

# monthly 99.5 percentile
p995 <- baseline %>%
  filter(!is.na(EC_FCH4)) %>%
  group_by(year_month) %>%
  summarise(p995 = quantile(EC_FCH4, 0.995, na.rm = TRUE), .groups = "drop")

# remove single extreme positive spikes under condensation conditions
baseline_tagged <- baseline %>%
  left_join(p995, by = "year_month") %>%
  mutate(
    extreme_pos = !is.na(p995) & !is.na(EC_FCH4) & EC_FCH4 > 0 & EC_FCH4 >= p995,
    
    prev_time = dplyr::lag(EC_TIMESTAMP_START),
    next_time = dplyr::lead(EC_TIMESTAMP_START),
    
    prev_adjacent = !is.na(prev_time) &
      as.numeric(difftime(EC_TIMESTAMP_START, prev_time, units = "mins")) == adjacency_min,
    next_adjacent = !is.na(next_time) &
      as.numeric(difftime(next_time, EC_TIMESTAMP_START, units = "mins")) == adjacency_min,
    
    prev_ext_close = dplyr::lag(extreme_pos, default = FALSE) & prev_adjacent,
    next_ext_close = dplyr::lead(extreme_pos, default = FALSE) & next_adjacent,
    
    singleton = extreme_pos & !(prev_ext_close | next_ext_close),
    suspect_condensation = singleton & is_night & near_sat
  )

CNHGU_EC_sub_clean <- baseline_tagged %>%
  filter(!suspect_condensation) %>%
  dplyr::select(
    -Td, -near_sat, -prev_time, -next_time,
    -prev_adjacent, -next_adjacent,
    -prev_ext_close, -next_ext_close,
    -singleton, -p995
  )

# remove known outlier (manual)
CNHGU_EC_sub_clean <- CNHGU_EC_sub_clean %>%
  filter(!(format(EC_TIMESTAMP_START, "%Y-%m-%d %H:%M:%S") == "2016-07-31 22:30:00"))

# save
write.csv(
  CNHGU_EC_sub_clean,
  "path/CN_HGU_EC_subset_28082025.csv",
  row.names = FALSE
)

### FI-Si2 ###

# read csv

FISi2_EC <- read.csv("path/FLX_FI-Si2_FLUXNET-CH4_HH_2012-2016_1-1.csv")

# convert timestamps to non-scientific format
FISi2_EC$TIMESTAMP_START <- format(FISi2_EC$TIMESTAMP_START, scientific = FALSE)
FISi2_EC$TIMESTAMP_END <- format(FISi2_EC$TIMESTAMP_END, scientific = FALSE)

# convert to datetime format

#create new column for TimeStart with seconds
FISi2_EC$TimeStart <- paste(as.character(FISi2_EC$TIMESTAMP_START),"00",sep="")

# convert to datetime

FISi2_EC$TIMESTAMP_START <- as_datetime(FISi2_EC$TimeStart)

#create new column for Time2 with seconds
FISi2_EC$TimeEnd <- paste(as.character(FISi2_EC$TIMESTAMP_END),"00",sep="")

# convert to datetime

FISi2_EC$TIMESTAMP_END <- as_datetime(FISi2_EC$TimeEnd)

# remove the extra columns

FISi2_EC <- subset(FISi2_EC, select=-c(TimeStart, TimeEnd))

# subset EC data to match with chamber data
FISi2_EC_sub <- subset(FISi2_EC, date(TIMESTAMP_START) >= "2012-06-26" & TIMESTAMP_END <= "2014-09-24 00:00:00")

# reset row names
rownames(FISi2_EC_sub) <- NULL

# replace -9999 with NA

FISi2_EC_sub <- replace(FISi2_EC_sub, FISi2_EC_sub == -9999, NA)

# add EC_prefix to all columns
FISi2_EC_sub <- FISi2_EC_sub %>% rename_with( ~ paste0("EC_", .x))

# add Site

FISi2_EC_sub$SITE <- "FI-SI2"

# move to front
FISi2_EC_sub <- FISi2_EC_sub %>%
  dplyr::select(SITE, everything())


# SET EC GAP-FILLED DATA AND CHAMBER DATA TO NA WHEN TIME GAP LONGER THAN 2 MONTHS

FISi2_EC_sub$EC_FCH4_F_ANNOPTLM_QC <- as.numeric(FISi2_EC_sub$EC_FCH4_F_ANNOPTLM_QC)

FISI2_new <- FISi2_EC_sub %>% mutate(across(c(EC_FCH4_F, EC_NEE_F, EC_H_F, EC_LE_F, EC_SW_IN_F,
                                              EC_SW_IN_F, EC_LW_IN_F, EC_NETRAD_F, EC_PPFD_IN_F, EC_VPD_F, EC_RH_F, EC_PA_F, EC_TA_F, EC_P_F, EC_WTD_F,
                                              EC_WS_F, EC_LE_F_ANNOPTLM, EC_NEE_F_ANNOPTLM, EC_FCH4_F_ANNOPTLM), ~ ifelse(EC_FCH4_F_ANNOPTLM_QC == 3, NA, .)))

# for calculating wind speed and direction averages, calculate u and v wind components:
# Calculate the u and v wind components
FISI2_new2$u.wind <- -FISI2_new2$EC_WS *sin(2*pi *FISI2_new2$EC_WD/360)
FISI2_new2$v.wind <- -FISI2_new2$EC_WS *cos(2*pi *FISI2_new2$EC_WD/360)
# gap-filled (no gap-filled WD)
FISI2_new2$u.wind.F <- -FISI2_new2$EC_WS_F *sin(2*pi *FISI2_new2$EC_WD/360)
FISI2_new2$v.wind.F <- -FISI2_new2$EC_WS_F *cos(2*pi *FISI2_new2$EC_WD/360)

# choose only the dates where chamber measurements available
ch_meas <- c(
  "2012-06-26", "2012-07-09", "2012-07-24", "2012-08-07",
  "2012-08-21", "2012-09-13", "2012-10-16",
  "2013-05-22", "2013-05-31", "2013-06-11", "2013-06-18",
  "2013-07-02", "2013-07-16", "2013-07-30", "2013-08-13",
  "2013-09-10",
  "2014-06-05", "2014-06-11", "2014-06-25", "2014-06-26",
  "2014-07-03", "2014-07-15", "2014-07-29",
  "2014-08-13", "2014-08-26", "2014-09-11"
)

# convert to Date
ch_meas_clean <- as.Date(ch_meas)

FISI2EC_ECfilt <- FISI2_new2[date(FISI2_new2$EC_TIMESTAMP_START) %in% ch_meas_clean,]

# check histograms for the predictor variables

FISI2_subset <- subset(FISI2EC_ECfilt, select = c(EC_NEE_F_ANNOPTLM, EC_LE_F_ANNOPTLM, EC_GPP_NT, EC_RECO_NT, EC_GPP_DT, EC_RECO_DT, 
                                                  EC_NEE_F, EC_H_F, EC_NETRAD_F, EC_PPFD_IN_F, EC_VPD_F, EC_RH_F, EC_PA_F, EC_TA_F,
                                                  EC_P_F, EC_WTD_F, EC_LE_F, EC_SW_IN_F, EC_LW_IN_F, EC_TS_1 ,EC_TS_2, EC_TS_3, EC_TS_4))


# transform to long format (2 columns)
df <- gather(FISI2_subset, key = "name", value = "value")

# plot histograms per name
ggplot(df) +
  geom_histogram(aes(value)) +
  facet_wrap(~name, ncol = 5, scales= "free")

# save
write.csv("path/FISI2EC_ECfilt.csv")

# calculate mean annual temperature and precipitation sum

# calculate precip sum for each year

# Summarize total annual precipitation
FISI2EC_ECfilt$EC_TIMESTAMP_START <- as_datetime(FISI2EC_ECfilt$EC_TIMESTAMP_START)

annual_precip <- FISI2EC_ECfilt %>%
  group_by(year(EC_TIMESTAMP_START)) %>%
  summarize(Annual_Precip = sum(EC_P, na.rm = T))

# Calculate the mean annual precipitation
mean(annual_precip$Annual_Precip)

# Calculate the annual mean temperature for each year
annual_mean_temp <- FISI2EC_ECfilt %>%
  group_by(year(EC_TIMESTAMP_START)) %>%
  summarize(Annual_Mean_Temp = mean(EC_TA_F, na.rm = TRUE))

# Calculate the overall mean annual temperature across all years
mean(annual_mean_temp$Annual_Mean_Temp)

# aggregate to weekly scale

FISI2_hh_week_aggr_ECfilt <- FISI2EC_ECfilt %>% 
  group_by(year(FISI2EC_ECfilt$EC_TIMESTAMP_START), week(FISI2EC_ECfilt$EC_TIMESTAMP_START)) %>%
  summarise_at(c("EC_FCH4", "EC_FCH4_F", "EC_FCH4_F_ANNOPTLM",
                 "EC_USTAR", "EC_GPP_NT", "EC_RECO_NT", "EC_GPP_DT", "EC_RECO_DT", "EC_NEE_F", "EC_H_F", "EC_LE_F", "EC_SW_IN_F", "EC_LW_IN_F", 
                 "EC_NETRAD_F", "EC_PPFD_IN_F", "EC_VPD_F", "EC_RH_F", 
                 "EC_PA_F", "EC_TA_F", "EC_P_F", "EC_WS_F", 
                 "EC_LE_F_ANNOPTLM", "EC_NEE_F_ANNOPTLM", "u.wind", "v.wind", "u.wind.F", "v.wind.F"), 
               list(mean = mean,
                    median = median), 
               na.rm = TRUE)

FISI2_hh_week_aggr_ECfilt <- as.data.frame(FISI2_hh_week_aggr_ECfilt)

# rename the wind component columns and date

FISI2_hh_week_aggr_ECfilt <- FISI2_hh_week_aggr_ECfilt %>% rename_with(~ c("Year", "Week_of_year","EC_WD_u.wind_mean", "EC_WD_v.wind_mean", "EC_WD_u.wind_mean_F", "EC_WD_v.wind_mean_F"), 
                                                                       all_of(c("year(FISI2EC_ECfilt$EC_TIMESTAMP_START)", "week(FISI2EC_ECfilt$EC_TIMESTAMP_START)","u.wind_mean", 
                                                                                "v.wind_mean", "u.wind.F_mean", "v.wind.F_mean")))

# calculate wind direction average
FISI2_hh_week_aggr_ECfilt$EC_WD_AVG <- (atan2(FISI2_hh_week_aggr_ECfilt$EC_WD_u.wind_mean, FISI2_hh_week_aggr_ECfilt$EC_WD_v.wind_mean) *360/2/pi) +180

# gap-filled:
FISI2_hh_week_aggr_ECfilt$EC_WD_AVG_F <- (atan2(FISI2_hh_week_aggr_ECfilt$EC_WD_u.wind_mean_F, FISI2_hh_week_aggr_ECfilt$EC_WD_v.wind_mean_F) *360/2/pi) +180

# aggregate to daily scale

FISI2_hh_d_aggr_ECfilt <- FISI2EC_ECfilt %>% 
  group_by(date(FISI2EC_ECfilt$EC_TIMESTAMP_START)) %>%
  summarise_at(c("EC_FCH4", "EC_FCH4_F", "EC_FCH4_F_ANNOPTLM",
                 "EC_USTAR", "EC_GPP_NT", "EC_RECO_NT", "EC_GPP_DT", "EC_RECO_DT", "EC_NEE_F", "EC_H_F", "EC_LE_F", "EC_SW_IN_F", "EC_LW_IN_F", 
                 "EC_NETRAD_F", "EC_PPFD_IN_F", "EC_VPD_F", "EC_RH_F", 
                 "EC_PA_F", "EC_TA_F", "EC_P_F", "EC_WS_F", 
                 "EC_LE_F_ANNOPTLM", "EC_NEE_F_ANNOPTLM", "u.wind", "v.wind", "u.wind.F", "v.wind.F"), 
               list(mean = mean,
                    median = median), 
               na.rm = TRUE)

FISI2_hh_d_aggr_ECfilt <- as.data.frame(FISI2_hh_d_aggr_ECfilt)

# rename columns

FISI2_hh_d_aggr_ECfilt <- FISI2_hh_d_aggr_ECfilt %>% rename_with(~ c("TIMESTAMP", "EC_WD_u.wind_mean", "EC_WD_v.wind_mean", "EC_WD_u.wind_mean_F", "EC_WD_v.wind_mean_F"), 
                                                                 all_of(c("date(FISI2EC_ECfilt$EC_TIMESTAMP_START)", "u.wind_mean", 
                                                                          "v.wind_mean", "u.wind.F_mean", "v.wind.F_mean")))

# calculate wind direction average
FISI2_hh_d_aggr_ECfilt$EC_WD_AVG <- (atan2(FISI2_hh_d_aggr_ECfilt$EC_WD_u.wind_mean, FISI2_hh_d_aggr_ECfilt$EC_WD_v.wind_mean) *360/2/pi) +180

# gap-filled:
FISI2_hh_d_aggr_ECfilt$EC_WD_AVG_F <- (atan2(FISI2_hh_d_aggr_ECfilt$EC_WD_u.wind_mean_F, FISI2_hh_d_aggr_ECfilt$EC_WD_v.wind_mean_F) *360/2/pi) +180

# aggregate to annual scale

FISI2_hh_yr_aggr_ECfilt <- FISI2EC_ECfilt %>% 
  group_by(year(FISI2EC_ECfilt$EC_TIMESTAMP_START)) %>%
  summarise_at(c("EC_FCH4", "EC_FCH4_F", "EC_FCH4_F_ANNOPTLM",
                 "EC_USTAR", "EC_GPP_NT", "EC_RECO_NT", "EC_GPP_DT", "EC_RECO_DT", "EC_NEE_F", "EC_H_F", "EC_LE_F", "EC_SW_IN_F", "EC_LW_IN_F", 
                 "EC_NETRAD_F", "EC_PPFD_IN_F", "EC_VPD_F", "EC_RH_F", 
                 "EC_PA_F", "EC_TA_F", "EC_P_F", "EC_WS_F", 
                 "EC_LE_F_ANNOPTLM", "EC_NEE_F_ANNOPTLM", "u.wind", "v.wind", "u.wind.F", "v.wind.F"), 
               list(mean = mean,
                    median = median), 
               na.rm = TRUE)

FISI2_hh_yr_aggr_ECfilt <- as.data.frame(FISI2_hh_yr_aggr_ECfilt)

FISI2_hh_yr_aggr_ECfilt <- FISI2_hh_yr_aggr_ECfilt %>% rename_with(~ c("EC_WD_u.wind_mean", "EC_WD_v.wind_mean", "EC_WD_u.wind_mean_F", "EC_WD_v.wind_mean_F"), 
                                                                   all_of(c("u.wind_mean", 
                                                                            "v.wind_mean", "u.wind.F_mean", "v.wind.F_mean")))
# calculate wind direction average
FISI2_hh_yr_aggr_ECfilt$EC_WD_AVG <- (atan2(FISI2_hh_yr_aggr_ECfilt$EC_WD_u.wind_mean, FISI2_hh_yr_aggr_ECfilt$EC_WD_v.wind_mean) *360/2/pi) +180

# gap-filled:
FISI2_hh_yr_aggr_ECfilt$EC_WD_AVG_F <- (atan2(FISI2_hh_yr_aggr_ECfilt$EC_WD_u.wind_mean_F, FISI2_hh_yr_aggr_ECfilt$EC_WD_v.wind_mean_F) *360/2/pi) +180


# aggregate to monthly scale

FISI2_hh_month_aggr_ECfilt <- FISI2EC_ECfilt %>% 
  group_by(year(FISI2EC_ECfilt$EC_TIMESTAMP_START), month(FISI2EC_ECfilt$EC_TIMESTAMP_START)) %>%
  summarise_at(c("EC_FCH4", "EC_FCH4_F", "EC_FCH4_F_ANNOPTLM",
                 "EC_USTAR", "EC_GPP_NT", "EC_RECO_NT", "EC_GPP_DT", "EC_RECO_DT", "EC_NEE_F", "EC_H_F", "EC_LE_F", "EC_SW_IN_F", "EC_LW_IN_F", 
                 "EC_NETRAD_F", "EC_PPFD_IN_F", "EC_VPD_F", "EC_RH_F", 
                 "EC_PA_F", "EC_TA_F", "EC_P_F", "EC_WS_F", 
                 "EC_LE_F_ANNOPTLM", "EC_NEE_F_ANNOPTLM", "u.wind", "v.wind", "u.wind.F", "v.wind.F"), 
               list(mean = mean,
                    median = median), 
               na.rm = TRUE)

FISI2_hh_month_aggr_ECfilt <- as.data.frame(FISI2_hh_month_aggr_ECfilt)

FISI2_hh_month_aggr_ECfilt <- FISI2_hh_month_aggr_ECfilt %>% rename_with(~ c("EC_WD_u.wind_mean", "EC_WD_v.wind_mean", "EC_WD_u.wind_mean_F", "EC_WD_v.wind_mean_F"), 
                                                                         all_of(c("u.wind_mean", 
                                                                                  "v.wind_mean", "u.wind.F_mean", "v.wind.F_mean")))

# calculate wind direction average
FISI2_hh_month_aggr_ECfilt$EC_WD_AVG <- (atan2(FISI2_hh_month_aggr_ECfilt$EC_WD_u.wind_mean, FISI2_hh_month_aggr_ECfilt$EC_WD_v.wind_mean) *360/2/pi) +180
 
# gap-filled:
FISI2_hh_month_aggr_ECfilt$EC_WD_AVG_F <- (atan2(FISI2_hh_month_aggr_ECfilt$EC_WD_u.wind_mean_F, FISI2_hh_month_aggr_ECfilt$EC_WD_v.wind_mean_F) *360/2/pi) +180


colnames(FISI2_hh_month_aggr_ECfilt )[2] <- "Month"
colnames(FISI2_hh_month_aggr_ECfilt )[1] <- "Year"

colnames(FISI2_hh_d_aggr_ECfilt)[1] <- "TIMESTAMP"
colnames(FISI2_hh_yr_aggr_ECfilt)[1] <- "Year"
colnames(FISI2_hh_month_aggr_ECfilt )[1] <- "Month"

# convert to datetime format
FISI2_hh_d_aggr_ECfilt$TIMESTAMP <- as_date(FISI2_hh_d_aggr_ECfilt$TIMESTAMP)

# add the missing extra columns

FISI2_hh_d_aggr_ECfilt$Year <- year(FISI2_hh_d_aggr_ECfilt$TIMESTAMP)
FISI2_hh_d_aggr_ECfilt$Month <- month(FISI2_hh_d_aggr_ECfilt$TIMESTAMP)
FISI2_hh_d_aggr_ECfilt$Day <- day(FISI2_hh_d_aggr_ECfilt$TIMESTAMP)
FISI2_hh_d_aggr_ECfilt$DOY <- yday(FISI2_hh_d_aggr_ECfilt$TIMESTAMP)

FISI2_hh_d_aggr_ECfilt$ch_method <- "manual"
FISI2_hh_d_aggr_ECfilt$SITE <- "FI-SI2"

FISI2_hh_yr_aggr_ECfilt$ch_method <- "manual"
FISI2_hh_yr_aggr_ECfilt$SITE <- "FI-SI2"

FISI2_hh_month_aggr_ECfilt$ch_method <- "manual"
FISI2_hh_month_aggr_ECfilt$SITE <- "FI-SI2"

FISI2_hh_week_aggr_ECfilt$ch_method <- "manual"
FISI2_hh_week_aggr_ECfilt$SITE <- "FI-SI2"


FISI2_hh_d_aggr_ECfilt <- FISI2_hh_d_aggr_ECfilt %>% relocate(SITE, .before = TIMESTAMP)
FISI2_hh_yr_aggr_ECfilt <- FISI2_hh_yr_aggr_ECfilt %>% relocate(SITE, .before = Year)
FISI2_hh_month_aggr_ECfilt <- FISI2_hh_month_aggr_ECfilt %>% relocate(SITE, .before = TIMESTAMP)
FISI2_hh_week_aggr_ECfilt <- FISI2_hh_week_aggr_ECfilt %>% relocate(SITE, .before = Year)


# save as .csv

write.csv(FISI2_hh_d_aggr_ECfilt, "path/FI_Si2_EC_subset_HH_TO_DD_AGGR_ECfilt.csv", row.names = FALSE)
write.csv(FISI2_hh_yr_aggr_ECfilt, "path/FI_Si2_EC_subset_HH_TO_YR_AGGR_ECfilt.csv", row.names = FALSE)
write.csv(FISI2_hh_month_aggr_ECfilt, "path/FI_Si2_EC_subset_HH_TO_MONTH_AGGR_ECfilt.csv", row.names = FALSE)
write.csv(FISI2_hh_week_aggr_ECfilt, "path/FI_Si2_EC_subset_HH_TO_WEEK_AGGR_ECfilt.csv", row.names = FALSE)

### SE-DEG ###

# read csv

SEDEG_EC <- read.csv("path/FLX_SE-Deg_FLUXNET-CH4_HH_2014-2018_1-1.csv")

SEDEG_EC$TIMESTAMP_END <- format(SEDEG_EC$TIMESTAMP_END, scientific = FALSE)

#create new column for Time2 with seconds
SEDEG_EC$TimeStart <- paste(as.character(SEDEG_EC$TIMESTAMP_START),"00",sep="")

# convert to datetime

SEDEG_EC$TIMESTAMP_START <- as_datetime(SEDEG_EC$TimeStart)

#create new column for Time2 with seconds
SEDEG_EC$TimeEnd <- paste(as.character(SEDEG_EC$TIMESTAMP_END),"00",sep="")

# convert to datetime

SEDEG_EC$TIMESTAMP_END <- as_datetime(SEDEG_EC$TimeEnd)

# subset EC data to match with chamber data
SEDEG_EC_sub <- subset(SEDEG_EC, date(TIMESTAMP_START) >= "2015-05-13" & TIMESTAMP_END < "2016-11-03 00:30:00")

# reset row names
rownames(SEDEG_EC_sub) <- NULL

# remove some columns that are not needed

SEDEG_EC_sub <- subset(SEDEG_EC_sub, select=-c(TimeStart, TimeEnd))

# replace -9999 with NA

SEDEG_EC_sub <- SEDEG_EC_sub %>% na_if(-9999)

# add EC_prefix to all columns
SEDEG_EC_sub <- SEDEG_EC_sub %>% rename_with( ~ paste0("EC_", .x))

# add Site

SEDEG_EC_sub$SITE <- "SE-DEG"

# move to front
SEDEG_EC_sub <- SEDEG_EC_sub %>%
  select(SITE, everything())

# save as .csv

write.csv(SEDEG_EC_sub, "path/SE_DEG_EC_subset.csv", row.names = FALSE)
SEDEG_EC_sub <- read.csv("path/SE_DEG_EC_subset.csv")

# calculate mean annual temperature and precipitation sum

# Summarize total annual precipitation
annual_precip <- SEDEG_EC_sub %>%
  group_by(year(EC_TIMESTAMP_START)) %>%
  summarize(Annual_Precip = sum(EC_P_F, na.rm = T))

# Calculate the mean annual precipitation
mean(annual_precip$Annual_Precip)

# Calculate the annual mean temperature for each year
annual_mean_temp <- SEDEG_EC_sub %>%
  group_by(year(EC_TIMESTAMP_START)) %>%
  summarize(Annual_Mean_Temp = mean(EC_TA_F, na.rm = TRUE))

# Calculate the overall mean annual temperature across all years
mean(annual_mean_temp$Annual_Mean_Temp)

### US-HO1 ###

# read csv

USHO1_EC <- read.csv("path/FLX_US-Ho1_FLUXNET-CH4_HH_2012-2018_1-1.csv")

USHO1_EC$TIMESTAMP_END <- format(USHO1_EC$TIMESTAMP_END, scientific = FALSE)

#create new column for Time2 with seconds
USHO1_EC$TimeStart <- paste(as.character(USHO1_EC$TIMESTAMP_START),"00",sep="")

# convert to datetime

USHO1_EC$TIMESTAMP_START <- as_datetime(USHO1_EC$TimeStart)

#create new column for Time2 with seconds
USHO1_EC$TimeEnd <- paste(as.character(USHO1_EC$TIMESTAMP_END),"00",sep="")

# convert to datetime

USHO1_EC$TIMESTAMP_END <- as_datetime(USHO1_EC$TimeEnd)

Sys.setenv(TZ = "UTC")

# subset EC data to match with chamber data
USHO1_EC_sub <- subset(USHO1_EC, date(TIMESTAMP_START) >= "2012-01-01" & TIMESTAMP_END <= "2017-01-01 00:00:00")

# delete last row
USHO1_EC_sub <- USHO1_EC_sub[-nrow(USHO1_EC_sub),]

# reset row names
rownames(USHO1_EC_sub) <- NULL

# remove some columns that are not needed

USHO1_EC_sub <- subset(USHO1_EC_sub, select=-c(TimeStart, TimeEnd))

# replace -9999 with NA

USHO1_EC_sub <- USHO1_EC_sub %>% na_if(-9999)

# add EC_prefix to all columns
USHO1_EC_sub <- USHO1_EC_sub %>% rename_with( ~ paste0("EC_", .x))

# add Site

USHO1_EC_sub$SITE <- "US-HO1"

# move to front
USHO1_EC_sub <- USHO1_EC_sub %>%
  select(SITE, everything())

# save as .csv

write.csv(USHO1_EC_sub, "path/US_HO1_EC_subset.csv", row.names = FALSE)

USHO1_EC_sub <- read.csv("path/US_HO1_EC_subset.csv")

# calculate mean annual temperature and precipitation sum

# calculate precip sum for each year

# Summarize total annual precipitation
annual_precip <- USHO1_EC_sub %>%
  group_by(year(EC_TIMESTAMP_START)) %>%
  summarize(Annual_Precip = sum(EC_P_F, na.rm = T))

# Calculate the mean annual precipitation
mean(annual_precip$Annual_Precip)

# Calculate the annual mean temperature for each year
annual_mean_temp <- USHO1_EC_sub %>%
  group_by(year(EC_TIMESTAMP_START)) %>%
  summarize(Annual_Mean_Temp = mean(EC_TA_F, na.rm = TRUE))

# Calculate the overall mean annual temperature across all years
mean(annual_mean_temp$Annual_Mean_Temp)


### US-LA1 ###

# read csv

USLA1_EC <- read.csv("path/FLX_US-LA1_FLUXNET-CH4_HH_2011-2012_1-1.csv")

# convert timestamps to non-scientific format
USLA1_EC$TIMESTAMP_START <- format(USLA1_EC$TIMESTAMP_START, scientific = FALSE)
USLA1_EC$TIMESTAMP_END <- format(USLA1_EC$TIMESTAMP_END, scientific = FALSE)

# convert to datetime format

#create new column for TimeStart with seconds
USLA1_EC$TimeStart <- paste(as.character(USLA1_EC$TIMESTAMP_START),"00",sep="")

# convert to datetime

USLA1_EC$TIMESTAMP_START <- as_datetime(USLA1_EC$TimeStart)

#create new column for Time2 with seconds
USLA1_EC$TimeEnd <- paste(as.character(USLA1_EC$TIMESTAMP_END),"00",sep="")

# convert to datetime

USLA1_EC$TIMESTAMP_END <- as_datetime(USLA1_EC$TimeEnd)

# remove the extra columns

USLA1_EC <- subset(USLA1_EC, select=-c(TimeStart, TimeEnd))

# convert to date format

USLA1_EC$TIMESTAMP_START <- as_datetime(USLA1_EC$TIMESTAMP_START)
USLA1_EC$TIMESTAMP_END <- as_datetime(USLA1_EC$TIMESTAMP_END)

# subset EC data to match with chamber data
USLA1_EC_sub <- subset(USLA1_EC, date(TIMESTAMP_START) >= "2012-03-19" & TIMESTAMP_END <= "2012-10-03 00:00:00")

# reset row names
rownames(USLA1_EC_sub) <- NULL

# replace -9999 with NA

USLA1_EC_sub <- replace(USLA1_EC_sub, USLA1_EC_sub == -9999, NA)


# add EC_prefix to all columns
USLA1_EC_sub <- USLA1_EC_sub %>% rename_with( ~ paste0("EC_", .x))

# add Site

USLA1_EC_sub$SITE <- "US-LA1"

# move to front
USLA1_EC_sub <- USLA1_EC_sub %>%
  dplyr::select(SITE, everything())



# SET EC GAP-FILLED DATA AND CHAMBER DATA TO NA WHEN TIME GAP LONGER THAN 2 MONTHS

USLA1_EC_sub$EC_FCH4_F_ANNOPTLM_QC <- as.numeric(USLA1_EC_sub$EC_FCH4_F_ANNOPTLM_QC)

USLA1_new <- USLA1_EC_sub %>% mutate(across(c(EC_FCH4_F, EC_NEE_F, EC_H_F, EC_LE_F, EC_SW_IN_F,
                                              EC_SW_IN_F, EC_LW_IN_F, EC_PPFD_IN_F, EC_VPD_F, EC_RH_F, EC_PA_F, EC_TA_F, EC_P_F,
                                              EC_WS_F, EC_TS_F, EC_LE_F_ANNOPTLM, EC_NEE_F_ANNOPTLM, EC_FCH4_F_ANNOPTLM), ~ ifelse(EC_FCH4_F_ANNOPTLM_QC == 3, NA, .)))

# choose only the dates where chamber measurements available

# define chamber measurement dates manually
ch_meas <- c(
  "2012-03-19",
  "2012-04-19",
  "2012-05-29",
  "2012-09-07",
  "2012-10-02"
)

# convert to Date
ch_meas_clean <- as.Date(ch_meas)

USLA1_new2$EC_TIMESTAMP_START <- as_datetime(USLA1_new2$EC_TIMESTAMP_START)

# create ECfilt df
USLA1EC_ECfilt <- USLA1_new2[date(USLA1_new2$EC_TIMESTAMP_START) %in% ch_meas,]

# calculate u and v wind components
USLA1EC_ECfilt$u.wind <- -USLA1EC_ECfilt$EC_WS *sin(2*pi *USLA1EC_ECfilt$EC_WD/360)
USLA1EC_ECfilt$v.wind <- -USLA1EC_ECfilt$EC_WS *cos(2*pi *USLA1EC_ECfilt$EC_WD/360)
# gap-filled (no gap-filled WD)
USLA1EC_ECfilt$u.wind.F <- -USLA1EC_ECfilt$EC_WS_F *sin(2*pi *USLA1EC_ECfilt$EC_WD/360)
USLA1EC_ECfilt$v.wind.F <- -USLA1EC_ECfilt$EC_WS_F *cos(2*pi *USLA1EC_ECfilt$EC_WD/360)

# save
write.csv("path/USLA1EC_ECfilt.csv")

# check histograms of the predictors
USLA1_subset <- subset(USLA1EC_ECfilt, select = c(EC_NEE_F_ANNOPTLM, EC_LE_F_ANNOPTLM, EC_GPP_NT, EC_RECO_NT, EC_GPP_DT, EC_RECO_DT, 
                                                  EC_NEE_F, EC_H_F, EC_PPFD_IN_F, EC_VPD_F, EC_RH_F, EC_PA_F, EC_TA_F,
                                                  EC_P_F, EC_VPD_F, EC_LE_F, EC_SW_IN_F, EC_LW_IN_F, EC_TS_F))


# transform to long format (2 columns)
df <- gather(USLA1_subset, key = "name", value = "value")

# plot histograms per name
ggplot(df) +
  geom_histogram(aes(value)) +
  facet_wrap(~name, ncol = 5, scales= "free")

# calculate mean annual temperature and precipitation sum

# calculate precip sum for each year

# Summarize total annual precipitation
annual_precip <- USLA1EC_ECfilt %>%
  group_by(year(EC_TIMESTAMP_START)) %>%
  summarize(Annual_Precip = sum(EC_P_F, na.rm = T))

# Calculate the mean annual precipitation
mean(annual_precip$Annual_Precip)

# Calculate the annual mean temperature for each year
annual_mean_temp <- USLA1EC_ECfilt %>%
  group_by(year(EC_TIMESTAMP_START)) %>%
  summarize(Annual_Mean_Temp = mean(EC_TA_F, na.rm = TRUE))

# Calculate the overall mean annual temperature across all years
mean(annual_mean_temp$Annual_Mean_Temp)

# aggregate to daily scale

USLA1_hh_d_aggr_ECfilt <- USLA1EC_ECfilt %>% 
  group_by(date(USLA1EC_ECfilt$EC_TIMESTAMP_START)) %>%
  summarise_at(c("EC_FCH4", "EC_FCH4_F", "EC_FCH4_F_ANNOPTLM",
                 "EC_USTAR", "EC_GPP_NT", "EC_RECO_NT", "EC_GPP_DT", "EC_RECO_DT", "EC_NEE_F", "EC_H_F", "EC_LE_F", "EC_SW_IN_F", "EC_LW_IN_F", 
                 "EC_PPFD_IN_F", "EC_VPD_F", "EC_RH_F", 
                 "EC_PA_F", "EC_TA_F", "EC_P_F", "EC_TS_F", "EC_WS_F", 
                 "EC_LE_F_ANNOPTLM", "EC_NEE_F_ANNOPTLM", "u.wind", "v.wind", "u.wind.F", "v.wind.F"), 
               list(mean = mean,
                    median = median), 
               na.rm = TRUE)

USLA1_hh_d_aggr_ECfilt <- as.data.frame(USLA1_hh_d_aggr_ECfilt)


# rename the wind component columns and date

USLA1_hh_d_aggr_ECfilt <- USLA1_hh_d_aggr_ECfilt %>% rename_with(~ c("TIMESTAMP", "EC_WD_u.wind_mean", "EC_WD_v.wind_mean", "EC_WD_u.wind_mean_F", "EC_WD_v.wind_mean_F"), 
                                                                 all_of(c("date(USLA1EC_ECfilt$EC_TIMESTAMP_START)", "u.wind_mean", 
                                                                          "v.wind_mean", "u.wind.F_mean", "v.wind.F_mean")))

# calculate wind direction average
USLA1_hh_d_aggr_ECfilt$EC_WD_AVG <- (atan2(USLA1_hh_d_aggr_ECfilt$EC_WD_u.wind_mean, USLA1_hh_d_aggr_ECfilt$EC_WD_v.wind_mean) *360/2/pi) +180

# gap-filled:
USLA1_hh_d_aggr_ECfilt$EC_WD_AVG_F <- (atan2(USLA1_hh_d_aggr_ECfilt$EC_WD_u.wind_mean_F, USLA1_hh_d_aggr_ECfilt$EC_WD_v.wind_mean_F) *360/2/pi) +180

# aggregate to annual scale

USLA1_hh_yr_aggr_ECfilt <- USLA1EC_ECfilt %>% 
  group_by(year(USLA1EC_ECfilt$EC_TIMESTAMP_START)) %>%
  summarise_at(c("EC_FCH4", "EC_FCH4_F", "EC_FCH4_F_ANNOPTLM",
                 "EC_USTAR", "EC_GPP_NT", "EC_RECO_NT", "EC_GPP_DT", "EC_RECO_DT", "EC_NEE_F", "EC_H_F", "EC_LE_F", "EC_SW_IN_F", "EC_LW_IN_F", 
                 "EC_PPFD_IN_F", "EC_VPD_F", "EC_RH_F", 
                 "EC_PA_F", "EC_TA_F", "EC_P_F", "EC_TS_F", "EC_WS_F", 
                 "EC_LE_F_ANNOPTLM", "EC_NEE_F_ANNOPTLM", "u.wind", "v.wind", "u.wind.F", "v.wind.F"), 
               list(mean = mean,
                    median = median), 
               na.rm = TRUE)

USLA1_hh_yr_aggr_ECfilt <- as.data.frame(USLA1_hh_yr_aggr_ECfilt)

# rename the wind component columns and date

USLA1_hh_yr_aggr_ECfilt <- USLA1_hh_yr_aggr_ECfilt %>% rename_with(~ c("Year", "EC_WD_u.wind_mean", "EC_WD_v.wind_mean", "EC_WD_u.wind_mean_F", "EC_WD_v.wind_mean_F"), 
                                                                   all_of(c("year(USLA1EC_ECfilt$EC_TIMESTAMP_START)", "u.wind_mean", 
                                                                            "v.wind_mean", "u.wind.F_mean", "v.wind.F_mean")))

# calculate wind direction average
USLA1_hh_yr_aggr_ECfilt$EC_WD_AVG <- (atan2(USLA1_hh_yr_aggr_ECfilt$EC_WD_u.wind_mean, USLA1_hh_yr_aggr_ECfilt$EC_WD_v.wind_mean) *360/2/pi) +180

# gap-filled:
USLA1_hh_yr_aggr_ECfilt$EC_WD_AVG_F <- (atan2(USLA1_hh_yr_aggr_ECfilt$EC_WD_u.wind_mean_F, USLA1_hh_yr_aggr_ECfilt$EC_WD_v.wind_mean_F) *360/2/pi) +180

USLA1_hh_yr_aggr_ECfilt <- subset(USLA1_hh_yr_aggr_ECfilt, select = -c(u.wind_median, v.wind_median, 
                                                                       u.wind.F_median, v.wind.F_median))

# aggregate to monthly scale

USLA1_hh_month_aggr_ECfilt <- USLA1EC_ECfilt %>% 
  group_by(year(USLA1EC_ECfilt$EC_TIMESTAMP_START), month(USLA1EC_ECfilt$EC_TIMESTAMP_START)) %>%
  summarise_at(c("EC_FCH4", "EC_FCH4_F", "EC_FCH4_F_ANNOPTLM",
                 "EC_USTAR", "EC_GPP_NT", "EC_RECO_NT", "EC_GPP_DT", "EC_RECO_DT", "EC_NEE_F", "EC_H_F", "EC_LE_F", "EC_SW_IN_F", "EC_LW_IN_F", 
                 "EC_PPFD_IN_F", "EC_VPD_F", "EC_RH_F", 
                 "EC_PA_F", "EC_TA_F", "EC_P_F", "EC_TS_F", "EC_WS_F", 
                 "EC_LE_F_ANNOPTLM", "EC_NEE_F_ANNOPTLM", "u.wind", "v.wind",
                 "u.wind.F", "v.wind.F"), 
               list(mean = mean,
                    median = median), 
               na.rm = TRUE)

USLA1_hh_month_aggr_ECfilt <- as.data.frame(USLA1_hh_month_aggr_ECfilt)

# rename the wind component columns and date

USLA1_hh_month_aggr_ECfilt <- USLA1_hh_month_aggr_ECfilt %>% rename_with(~ c("Year", "Month","EC_WD_u.wind_mean", "EC_WD_v.wind_mean", "EC_WD_u.wind_mean_F", "EC_WD_v.wind_mean_F"), 
                                                                         all_of(c("year(USLA1EC_ECfilt$EC_TIMESTAMP_START)", "month(USLA1EC_ECfilt$EC_TIMESTAMP_START)","u.wind_mean", 
                                                                                  "v.wind_mean", "u.wind.F_mean", "v.wind.F_mean")))

# calculate wind direction average
USLA1_hh_month_aggr_ECfilt$EC_WD_AVG <- (atan2(USLA1_hh_month_aggr_ECfilt$EC_WD_u.wind_mean, USLA1_hh_month_aggr_ECfilt$EC_WD_v.wind_mean) *360/2/pi) +180

# gap-filled:
USLA1_hh_month_aggr_ECfilt$EC_WD_AVG_F <- (atan2(USLA1_hh_month_aggr_ECfilt$EC_WD_u.wind_mean_F, USLA1_hh_month_aggr_ECfilt$EC_WD_v.wind_mean_F) *360/2/pi) +180

# weekly

USLA1_hh_week_aggr_ECfilt <- USLA1EC_ECfilt %>% 
  group_by(year(USLA1EC_ECfilt$EC_TIMESTAMP_START), week(USLA1EC_ECfilt$EC_TIMESTAMP_START)) %>%
  summarise_at(c("EC_FCH4", "EC_FCH4_F", "EC_FCH4_F_ANNOPTLM",
                 "EC_USTAR", "EC_GPP_NT", "EC_RECO_NT", "EC_GPP_DT", "EC_RECO_DT", "EC_NEE_F", "EC_H_F", "EC_LE_F", "EC_SW_IN_F", "EC_LW_IN_F", 
                 "EC_PPFD_IN_F", "EC_VPD_F", "EC_RH_F", 
                 "EC_PA_F", "EC_TA_F", "EC_P_F", "EC_TS_F", "EC_WS_F", 
                 "EC_LE_F_ANNOPTLM", "EC_NEE_F_ANNOPTLM", "u.wind", "v.wind",
                 "u.wind.F", "v.wind.F"), 
               list(mean = mean,
                    median = median), 
               na.rm = TRUE)

USLA1_hh_week_aggr_ECfilt <- as.data.frame(USLA1_hh_week_aggr_ECfilt)

# rename the wind component columns and date

USLA1_hh_week_aggr_ECfilt <- USLA1_hh_week_aggr_ECfilt %>% rename_with(~ c("Year", "Week_of_year","EC_WD_u.wind_mean", "EC_WD_v.wind_mean", "EC_WD_u.wind_mean_F", "EC_WD_v.wind_mean_F"), 
                                                                       all_of(c("year(USLA1EC_ECfilt$EC_TIMESTAMP_START)", "week(USLA1EC_ECfilt$EC_TIMESTAMP_START)","u.wind_mean", 
                                                                                "v.wind_mean", "u.wind.F_mean", "v.wind.F_mean")))

# calculate wind direction average
USLA1_hh_week_aggr_ECfilt$EC_WD_AVG <- (atan2(USLA1_hh_week_aggr_ECfilt$EC_WD_u.wind_mean, USLA1_hh_week_aggr_ECfilt$EC_WD_v.wind_mean) *360/2/pi) +180

# gap-filled:
USLA1_hh_week_aggr_ECfilt$EC_WD_AVG_F <- (atan2(USLA1_hh_week_aggr_ECfilt$EC_WD_u.wind_mean_F, USLA1_hh_week_aggr_ECfilt$EC_WD_v.wind_mean_F) *360/2/pi) +180

USLA1_hh_week_aggr_ECfilt <- subset(USLA1_hh_week_aggr_ECfilt, select = -c(u.wind_median, v.wind_median, 
                                                                           u.wind.F_median, v.wind.F_median))

# rename the date column

colnames(USLA1_hh_d_aggr_ECfilt)[1] <- "TIMESTAMP"

# convert to datetime format
USLA1_hh_d_aggr_ECfilt$TIMESTAMP <- as_date(USLA1_hh_d_aggr_ECfilt$TIMESTAMP)

# add the missing extra columns

USLA1_hh_d_aggr_ECfilt$Year <- year(USLA1_hh_d_aggr_ECfilt$TIMESTAMP)
USLA1_hh_d_aggr_ECfilt$Month <- month(USLA1_hh_d_aggr_ECfilt$TIMESTAMP)
USLA1_hh_d_aggr_ECfilt$Day <- day(USLA1_hh_d_aggr_ECfilt$TIMESTAMP)
USLA1_hh_d_aggr_ECfilt$DOY <- yday(USLA1_hh_d_aggr_ECfilt$TIMESTAMP)

USLA1_hh_d_aggr_ECfilt$ch_method <- "manual"
USLA1_hh_d_aggr_ECfilt$SITE <- "US-LA1"

USLA1_hh_yr_aggr_ECfilt$ch_method <- "manual"
USLA1_hh_yr_aggr_ECfilt$SITE <- "US-LA1"

USLA1_hh_month_aggr_ECfilt$ch_method <- "manual"
USLA1_hh_month_aggr_ECfilt$SITE <- "US-LA1"

USLA1_hh_week_aggr_ECfilt$ch_method <- "manual"
USLA1_hh_week_aggr_ECfilt$SITE <- "US-LA1"


USLA1_hh_d_aggr_ECfilt <- USLA1_hh_d_aggr_ECfilt %>% relocate(SITE, .before = TIMESTAMP)
USLA1_hh_yr_aggr_ECfilt <- USLA1_hh_yr_aggr_ECfilt %>% relocate(SITE, .before = Year)
USLA1_hh_month_aggr_ECfilt <- USLA1_hh_month_aggr_ECfilt %>% relocate(SITE, .before = Year)
USLA1_hh_week_aggr_ECfilt <- USLA1_hh_week_aggr_ECfilt %>% relocate(SITE, .before = Year)


# save as .csv

write.csv(USLA1_hh_d_aggr_ECfilt, "path/US_LA1_EC_subset_HH_TO_DD_AGGR_ECfilt.csv", row.names = FALSE)
write.csv(USLA1_hh_yr_aggr_ECfilt, "path/US_LA1_EC_subset_HH_TO_YR_AGGR_ECfilt.csv", row.names = FALSE)
write.csv(USLA1_hh_month_aggr_ECfilt, "path/US_LA1_EC_subset_HH_TO_MONTH_AGGR_ECfilt.csv", row.names = FALSE)
write.csv(USLA1_hh_week_aggr_ECfilt, "path/US_LA1_EC_subset_HH_TO_WEEK_AGGR_ECfilt.csv", row.names = FALSE)

### US-LA2 ###

# read csv

USLA2_EC <- read.csv("path/FLX_US-LA2_FLUXNET-CH4_HH_2011-2013_1-1.csv")

#create new column for TimeStart with seconds
USLA2_EC$TimeStart <- paste(as.character(USLA2_EC$TIMESTAMP_START),"00",sep="")

# convert to datetime

USLA2_EC$TIMESTAMP_START <- as_datetime(USLA2_EC$TimeStart)

#create new column for Time2 with seconds
USLA2_EC$TimeEnd <- paste(as.character(USLA2_EC$TIMESTAMP_END),"00",sep="")

# convert to datetime

USLA2_EC$TIMESTAMP_END <- as_datetime(USLA2_EC$TimeEnd)

# remove the extra columns

USLA2_EC <- subset(USLA2_EC, select=-c(TimeStart, TimeEnd))

# convert to date format

USLA2_EC$TIMESTAMP_START <- as_datetime(USLA2_EC$TIMESTAMP_START)
USLA2_EC$TIMESTAMP_END <- as_datetime(USLA2_EC$TIMESTAMP_END)

# subset EC data to match with chamber data
USLA2_EC_sub <- subset(USLA2_EC, date(TIMESTAMP_START) >= "2012-03-20" & TIMESTAMP_END <= "2013-10-30 00:00:00")

# reset row names
rownames(USLA2_EC_sub) <- NULL

# replace -9999 with NA

USLA2_EC_sub <- replace(USLA2_EC_sub, USLA2_EC_sub == -9999, NA)


# add EC_prefix to all columns
USLA2_EC_sub <- USLA2_EC_sub %>% rename_with( ~ paste0("EC_", .x))

# add Site

USLA2_EC_sub$SITE <- "US-LA2"

# move to front
USLA2_EC_sub <- USLA2_EC_sub %>%
  dplyr::select(SITE, everything())

# SET EC GAP-FILLED DATA AND CHAMBER DATA TO NA WHEN TIME GAP LONGER THAN 2 MONTHS 

USLA2_EC_sub$EC_FCH4_F_ANNOPTLM_QC <- as.numeric(USLA2_EC_sub$EC_FCH4_F_ANNOPTLM_QC)

USLA2_new <- USLA2_EC_sub %>% mutate(across(c(EC_FCH4_F, EC_NEE_F, EC_H_F, EC_LE_F, EC_SW_IN_F,
                                              EC_SW_IN_F, EC_LW_IN_F, EC_PPFD_IN_F, EC_VPD_F, EC_RH_F, EC_PA_F, EC_TA_F, EC_P_F,
                                              EC_WS_F, EC_TS_F, EC_LE_F_ANNOPTLM, EC_NEE_F_ANNOPTLM, EC_FCH4_F_ANNOPTLM), ~ ifelse(EC_FCH4_F_ANNOPTLM_QC == 3, NA, .)))

# choose only the dates where chamber measurements available
# define chamber measurement dates manually
ch_meas <- c(
  "2012-03-20",
  "2012-04-20",
  "2012-05-30",
  "2012-08-09",
  "2012-10-01",
  "2013-05-08",
  "2013-06-27",
  "2013-07-22",
  "2013-09-04",
  "2013-10-29"
)

# convert to Date
ch_meas <- as.Date(ch_meas)

# create ECfilt df
USLA2EC_ECfilt <- USLA2_new[date(USLA2_new2$EC_TIMESTAMP_START) %in% ch_meas,]

# # Calculate the u and v wind components
USLA2EC_ECfilt$u.wind <- -USLA2EC_ECfilt$EC_WS *sin(2*pi *USLA2EC_ECfilt$EC_WD/360)
USLA2EC_ECfilt$v.wind <- -USLA2EC_ECfilt$EC_WS *cos(2*pi *USLA2EC_ECfilt$EC_WD/360)
# gap-filled (no gap-filled WD)
USLA2EC_ECfilt$u.wind.F <- -USLA2EC_ECfilt$EC_WS_F *sin(2*pi *USLA2EC_ECfilt$EC_WD/360)
USLA2EC_ECfilt$v.wind.F <- -USLA2EC_ECfilt$EC_WS_F *cos(2*pi *USLA2EC_ECfilt$EC_WD/360)

# check predictor histograms
USLA2_subset <- subset(USLA2EC_ECfilt, select = c(EC_NEE_F_ANNOPTLM, EC_LE_F_ANNOPTLM, EC_GPP_NT, EC_RECO_NT, EC_GPP_DT, EC_RECO_DT, 
                                                  EC_NEE_F, EC_H_F, EC_PPFD_IN_F, EC_VPD_F, EC_RH_F, EC_PA_F, EC_TA_F,
                                                  EC_P_F, EC_VPD_F, EC_LE_F, EC_SW_IN_F, EC_LW_IN_F, EC_TS_F))

# transform to long format (2 columns)
df <- gather(USLA2_subset, key = "name", value = "value")

# plot histograms per name
ggplot(df) +
  geom_histogram(aes(value)) +
  facet_wrap(~name, ncol = 5, scales= "free")

# save
write.csv("path/USLA2EC_ECfilt.csv")

# calculate mean annual temperature and precipitation sum

# calculate precip sum for each year

# Summarize total annual precipitation
annual_precip <- USLA2EC_ECfilt %>%
  group_by(year(EC_TIMESTAMP_START)) %>%
  summarize(Annual_Precip = sum(EC_P_F, na.rm = T))

# Calculate the mean annual precipitation
mean(annual_precip$Annual_Precip)

# Calculate the annual mean temperature for each year
annual_mean_temp <- USLA2EC_ECfilt %>%
  group_by(year(EC_TIMESTAMP_START)) %>%
  summarize(Annual_Mean_Temp = mean(EC_TA_F, na.rm = TRUE))

# Calculate the overall mean annual temperature across all years
mean(annual_mean_temp$Annual_Mean_Temp)

# aggregate to daily scale

USLA2_hh_d_aggr_ECfilt <- USLA2EC_ECfilt %>% 
  group_by(date(USLA2EC_ECfilt$EC_TIMESTAMP_START)) %>%
  summarise_at(c("EC_FCH4", "EC_FCH4_F", "EC_FCH4_F_ANNOPTLM",
                 "EC_USTAR", "EC_GPP_NT", "EC_RECO_NT", "EC_GPP_DT", "EC_RECO_DT", "EC_NEE_F", "EC_H_F", "EC_LE_F", "EC_SW_IN_F", "EC_LW_IN_F", 
                 "EC_PPFD_IN_F", "EC_VPD_F", "EC_RH_F", 
                 "EC_PA_F", "EC_TA_F", "EC_P_F", "EC_TS_F", "EC_WS_F", 
                 "EC_LE_F_ANNOPTLM", "EC_NEE_F_ANNOPTLM", "u.wind", "v.wind", "u.wind.F", "v.wind.F"), 
               list(mean = mean,
                    median = median), 
               na.rm = TRUE)

USLA2_hh_d_aggr_ECfilt <- as.data.frame(USLA2_hh_d_aggr_ECfilt)

# rename the wind component columns and date

USLA2_hh_d_aggr_ECfilt <- USLA2_hh_d_aggr_ECfilt %>% rename_with(~ c("TIMESTAMP", "EC_WD_u.wind_mean", "EC_WD_v.wind_mean", "EC_WD_u.wind_mean_F", "EC_WD_v.wind_mean_F"), 
                                                                 all_of(c("date(USLA2EC_ECfilt$EC_TIMESTAMP_START)", "u.wind_mean", 
                                                                          "v.wind_mean", "u.wind.F_mean", "v.wind.F_mean")))
# calculate wind direction average
USLA2_hh_d_aggr_ECfilt$EC_WD_AVG <- (atan2(USLA2_hh_d_aggr_ECfilt$EC_WD_u.wind_mean, USLA2_hh_d_aggr_ECfilt$EC_WD_v.wind_mean) *360/2/pi) +180

# gap-filled:
USLA2_hh_d_aggr_ECfilt$EC_WD_AVG_F <- (atan2(USLA2_hh_d_aggr_ECfilt$EC_WD_u.wind_mean_F, USLA2_hh_d_aggr_ECfilt$EC_WD_v.wind_mean_F) *360/2/pi) +180

# aggregate to annual scale

USLA2_hh_yr_aggr_ECfilt <- USLA2EC_ECfilt %>% 
  group_by(year(USLA2EC_ECfilt$EC_TIMESTAMP_START)) %>%
  summarise_at(c("EC_FCH4", "EC_FCH4_F", "EC_FCH4_F_ANNOPTLM",
                 "EC_USTAR", "EC_GPP_NT", "EC_RECO_NT", "EC_GPP_DT", "EC_RECO_DT", "EC_NEE_F", "EC_H_F", "EC_LE_F", "EC_SW_IN_F", "EC_LW_IN_F", 
                 "EC_PPFD_IN_F", "EC_VPD_F", "EC_RH_F", 
                 "EC_PA_F", "EC_TA_F", "EC_P_F", "EC_TS_F", "EC_WS_F", 
                 "EC_LE_F_ANNOPTLM", "EC_NEE_F_ANNOPTLM", "u.wind", "v.wind", "u.wind.F", "v.wind.F"), 
               list(mean = mean,
                    median = median), 
               na.rm = TRUE)

USLA2_hh_yr_aggr_ECfilt <- as.data.frame(USLA2_hh_yr_aggr_ECfilt)

# rename the wind component columns and date

USLA2_hh_yr_aggr_ECfilt <- USLA2_hh_yr_aggr_ECfilt %>% rename_with(~ c("Year", "EC_WD_u.wind_mean", "EC_WD_v.wind_mean", "EC_WD_u.wind_mean_F", "EC_WD_v.wind_mean_F"), 
                                                                   all_of(c("year(USLA2EC_ECfilt$EC_TIMESTAMP_START)", "u.wind_mean", 
                                                                            "v.wind_mean", "u.wind.F_mean", "v.wind.F_mean")))
# calculate wind direction average
USLA2_hh_yr_aggr_ECfilt$EC_WD_AVG <- (atan2(USLA2_hh_yr_aggr_ECfilt$EC_WD_u.wind_mean, USLA2_hh_yr_aggr_ECfilt$EC_WD_v.wind_mean) *360/2/pi) +180

# gap-filled:
USLA2_hh_yr_aggr_ECfilt$EC_WD_AVG_F <- (atan2(USLA2_hh_yr_aggr_ECfilt$EC_WD_u.wind_mean_F, USLA2_hh_yr_aggr_ECfilt$EC_WD_v.wind_mean_F) *360/2/pi) +180

USLA2_hh_yr_aggr_ECfilt <- subset(USLA2_hh_yr_aggr_ECfilt, select = -c(u.wind_median, v.wind_median, 
                                                                       u.wind.F_median, v.wind.F_median))

# aggregate to monthly scale

USLA2_hh_month_aggr_ECfilt <- USLA2EC_ECfilt %>% 
  group_by(year(USLA2EC_ECfilt$EC_TIMESTAMP_START), month(USLA2EC_ECfilt$EC_TIMESTAMP_START)) %>%
  summarise_at(c("EC_FCH4", "EC_FCH4_F", "EC_FCH4_F_ANNOPTLM",
                 "EC_USTAR", "EC_GPP_NT", "EC_RECO_NT", "EC_GPP_DT", "EC_RECO_DT", "EC_NEE_F", "EC_H_F", "EC_LE_F", "EC_SW_IN_F", "EC_LW_IN_F", 
                 "EC_PPFD_IN_F", "EC_VPD_F", "EC_RH_F", 
                 "EC_PA_F", "EC_TA_F", "EC_P_F", "EC_TS_F", "EC_WS_F", 
                 "EC_LE_F_ANNOPTLM", "EC_NEE_F_ANNOPTLM", "u.wind", "v.wind", "u.wind.F", "v.wind.F"), 
               list(mean = mean,
                    median = median), 
               na.rm = TRUE)

USLA2_hh_month_aggr_ECfilt <- as.data.frame(USLA2_hh_month_aggr_ECfilt)

# rename the wind component columns and date

USLA2_hh_month_aggr_ECfilt <- USLA2_hh_month_aggr_ECfilt %>% rename_with(~ c("Year", "Month","EC_WD_u.wind_mean", "EC_WD_v.wind_mean", "EC_WD_u.wind_mean_F", "EC_WD_v.wind_mean_F"), 
                                                                         all_of(c("year(USLA2EC_ECfilt$EC_TIMESTAMP_START)", "month(USLA2EC_ECfilt$EC_TIMESTAMP_START)","u.wind_mean", 
                                                                                  "v.wind_mean", "u.wind.F_mean", "v.wind.F_mean")))
# calculate wind direction average
USLA2_hh_month_aggr_ECfilt$EC_WD_AVG <- (atan2(USLA2_hh_month_aggr_ECfilt$EC_WD_u.wind_mean, USLA2_hh_month_aggr_ECfilt$EC_WD_v.wind_mean) *360/2/pi) +180

# gap-filled:
USLA2_hh_month_aggr_ECfilt$EC_WD_AVG_F <- (atan2(USLA2_hh_month_aggr_ECfilt$EC_WD_u.wind_mean_F, USLA2_hh_month_aggr_ECfilt$EC_WD_v.wind_mean_F) *360/2/pi) +180

# aggregate to weekly scale

USLA2_hh_week_aggr_ECfilt <- USLA2EC_ECfilt %>% 
  group_by(year(USLA2EC_ECfilt$EC_TIMESTAMP_START), week(USLA2EC_ECfilt$EC_TIMESTAMP_START)) %>%
  summarise_at(c("EC_FCH4", "EC_FCH4_F", "EC_FCH4_F_ANNOPTLM",
                 "EC_USTAR", "EC_GPP_NT", "EC_RECO_NT", "EC_GPP_DT", "EC_RECO_DT", "EC_NEE_F", "EC_H_F", "EC_LE_F", "EC_SW_IN_F", "EC_LW_IN_F", 
                 "EC_PPFD_IN_F", "EC_VPD_F", "EC_RH_F", 
                 "EC_PA_F", "EC_TA_F", "EC_P_F", "EC_TS_F", "EC_WS_F", 
                 "EC_LE_F_ANNOPTLM", "EC_NEE_F_ANNOPTLM", "u.wind", "v.wind", "u.wind.F", "v.wind.F"), 
               list(mean = mean,
                    median = median), 
               na.rm = TRUE)

USLA2_hh_week_aggr_ECfilt <- as.data.frame(USLA2_hh_week_aggr_ECfilt)

# rename the wind component columns and date

USLA2_hh_week_aggr_ECfilt <- USLA2_hh_week_aggr_ECfilt %>% rename_with(~ c("Year", "Week_of_year","EC_WD_u.wind_mean", "EC_WD_v.wind_mean", "EC_WD_u.wind_mean_F", "EC_WD_v.wind_mean_F"), 
                                                                       all_of(c("year(USLA2EC_ECfilt$EC_TIMESTAMP_START)", "week(USLA2EC_ECfilt$EC_TIMESTAMP_START)","u.wind_mean", 
                                                                                "v.wind_mean", "u.wind.F_mean", "v.wind.F_mean")))
# calculate wind direction average
USLA2_hh_week_aggr_ECfilt$EC_WD_AVG <- (atan2(USLA2_hh_week_aggr_ECfilt$EC_WD_u.wind_mean, USLA2_hh_week_aggr_ECfilt$EC_WD_v.wind_mean) *360/2/pi) +180

# gap-filled:
USLA2_hh_week_aggr_ECfilt$EC_WD_AVG_F <- (atan2(USLA2_hh_week_aggr_ECfilt$EC_WD_u.wind_mean_F, USLA2_hh_week_aggr_ECfilt$EC_WD_v.wind_mean_F) *360/2/pi) +180

# rename the date column

colnames(USLA2_hh_d_aggr_ECfilt)[1] <- "TIMESTAMP"

# convert to datetime format
USLA2_hh_d_aggr_ECfilt$TIMESTAMP <- as_date(USLA2_hh_d_aggr_ECfilt$TIMESTAMP)

# add the missing extra columns

USLA2_hh_d_aggr_ECfilt$Year <- year(USLA2_hh_d_aggr_ECfilt$TIMESTAMP)
USLA2_hh_d_aggr_ECfilt$Month <- month(USLA2_hh_d_aggr_ECfilt$TIMESTAMP)
USLA2_hh_d_aggr_ECfilt$Day <- day(USLA2_hh_d_aggr_ECfilt$TIMESTAMP)
USLA2_hh_d_aggr_ECfilt$DOY <- yday(USLA2_hh_d_aggr_ECfilt$TIMESTAMP)

USLA2_hh_d_aggr_ECfilt$ch_method <- "manual"
USLA2_hh_d_aggr_ECfilt$SITE <- "US-LA2"

USLA2_hh_yr_aggr_ECfilt$ch_method <- "manual"
USLA2_hh_yr_aggr_ECfilt$SITE <- "US-LA2"

USLA2_hh_month_aggr_ECfilt$ch_method <- "manual"
USLA2_hh_month_aggr_ECfilt$SITE <- "US-LA2"

USLA2_hh_week_aggr_ECfilt$ch_method <- "manual"
USLA2_hh_week_aggr_ECfilt$SITE <- "US-LA2"


USLA2_hh_d_aggr_ECfilt <- USLA2_hh_d_aggr_ECfilt %>% relocate(SITE, .before = TIMESTAMP)
USLA2_hh_yr_aggr_ECfilt <- USLA2_hh_yr_aggr_ECfilt %>% relocate(SITE, .before = Year)
USLA2_hh_month_aggr_ECfilt <- USLA2_hh_month_aggr_ECfilt %>% relocate(SITE, .before = Year)
USLA2_hh_week_aggr_ECfilt <- USLA2_hh_week_aggr_ECfilt %>% relocate(SITE, .before = Year)


# save as .csv

write.csv(USLA2_hh_d_aggr_ECfilt, "path/US_LA2_EC_subset_HH_TO_DD_AGGR_ECfilt.csv", row.names = FALSE)
write.csv(USLA2_hh_yr_aggr_ECfilt, "path/US_LA2_EC_subset_HH_TO_YR_AGGR_ECfilt.csv", row.names = FALSE)
write.csv(USLA2_hh_month_aggr_ECfilt, "path/US_LA2_EC_subset_HH_TO_MONTH_AGGR_ECfilt.csv", row.names = FALSE)
write.csv(USLA2_hh_week_aggr_ECfilt, "path/US_LA2_EC_subset_HH_TO_WEEK_AGGR_ECfilt.csv", row.names = FALSE)

### US-LOS ###

# read csv

USLOS_EC <- read.csv("path/FLX_US-Los_FLUXNET-CH4_HH_2014-2018_1-1.csv")

USLOS_EC$TIMESTAMP_END <- format(USLOS_EC$TIMESTAMP_END, scientific = FALSE)

#create new column for Time2 with seconds
USLOS_EC$TimeStart <- paste(as.character(USLOS_EC$TIMESTAMP_START),"00",sep="")

# convert to datetime

USLOS_EC$TIMESTAMP_START <- as_datetime(USLOS_EC$TimeStart)

#create new column for Time2 with seconds
USLOS_EC$TimeEnd <- paste(as.character(USLOS_EC$TIMESTAMP_END),"00",sep="")

# convert to datetime

USLOS_EC$TIMESTAMP_END <- as_datetime(USLOS_EC$TimeEnd)

Sys.setenv(TZ = "UTC")
# subset EC data to match with chamber data
USLOS_EC_sub <- subset(USLOS_EC, date(TIMESTAMP_START) >= "2015-06-17" & date(TIMESTAMP_END) <= "2015-08-06")

# reset row names
rownames(USLOS_EC_sub) <- NULL

# remove some columns that are not needed

USLOS_EC_sub <- subset(USLOS_EC_sub, select=-c(TimeStart, TimeEnd))

# replace -9999 with NA

USLOS_EC_sub <- USLOS_EC_sub %>% na_if(-9999)

# add EC_prefix to all columns
USLOS_EC_sub <- USLOS_EC_sub %>% rename_with( ~ paste0("EC_", .x))

# add Site

USLOS_EC_sub$SITE <- "US-LOS"

# move to front
USLOS_EC_sub <- USLOS_EC_sub %>%
  select(SITE, everything())

# save as .csv

write.csv(USLOS_EC_sub, "path/US_LOS_EC_subset.csv", row.names = FALSE)
USLOS_EC_sub <- read.csv("path/US_LOS_EC_subset.csv")

# calculate mean annual temperature and precipitation sum

# calculate precip sum for each year

# Summarize total annual precipitation
annual_precip <- USLOS_EC_sub %>%
  group_by(year(EC_TIMESTAMP_START)) %>%
  summarize(Annual_Precip = sum(EC_P_F, na.rm = T))

# Calculate the mean annual precipitation
mean(annual_precip$Annual_Precip)

# Calculate the annual mean temperature for each year
annual_mean_temp <- USLOS_EC_sub %>%
  group_by(year(EC_TIMESTAMP_START)) %>%
  summarize(Annual_Mean_Temp = mean(EC_TA_F, na.rm = TRUE))

# Calculate the overall mean annual temperature across all years
mean(annual_mean_temp$Annual_Mean_Temp)


### US-OWC ###

# read csv

USOWC_EC <- read.csv("path/FLX_US-OWC_FLUXNET-CH4_HH_2015-2016_1-1.csv")

USOWC_EC$TIMESTAMP_START <- format(USOWC_EC$TIMESTAMP_START, scientific = FALSE)
USOWC_EC$TIMESTAMP_END <- format(USOWC_EC$TIMESTAMP_END, scientific = FALSE)

#create new column for Time2 with seconds
USOWC_EC$TimeStart <- paste(as.character(USOWC_EC$TIMESTAMP_START),"00",sep="")

# convert to datetime

USOWC_EC$TIMESTAMP_START <- as_datetime(USOWC_EC$TimeStart)

#create new column for Time2 with seconds
USOWC_EC$TimeEnd <- paste(as.character(USOWC_EC$TIMESTAMP_END),"00",sep="")

# convert to datetime

USOWC_EC$TIMESTAMP_END <- as_datetime(USOWC_EC$TimeEnd)

Sys.setenv(TZ = "UTC")
# subset EC data to match with chamber data
USOWC_EC_sub <- subset(USOWC_EC, TIMESTAMP_START >= "2015-06-17 06:00:00" & TIMESTAMP_END <= "2018-09-06 21:30:00")

# reset row names
rownames(USOWC_EC_sub) <- NULL

# remove some columns that are not needed

USOWC_EC_sub <- subset(USOWC_EC_sub, select=-c(TimeStart, TimeEnd))

# replace -9999 with NA

USOWC_EC_sub <- USOWC_EC_sub %>% na_if(-9999)

# add EC_prefix to all columns
USOWC_EC_sub <- USOWC_EC_sub %>% rename_with( ~ paste0("EC_", .x))

# add Site

USOWC_EC_sub$SITE <- "US-OWC"

# move to front
USOWC_EC_sub <- USOWC_EC_sub %>%
  select(SITE, everything())

# save as .csv

write.csv(USOWC_EC_sub, "path/US_OWC_EC_subset.csv", row.names = FALSE)
USOWC_EC_sub <- read.csv("path/US_OWC_EC_subset.csv")

# calculate mean annual temperature and precipitation sum

# calculate precip sum for each year

# Summarize total annual precipitation
annual_precip <- USOWC_EC_sub %>%
  group_by(year(EC_TIMESTAMP_START)) %>%
  summarize(Annual_Precip = sum(EC_P_F, na.rm = T))

# Calculate the mean annual precipitation
mean(annual_precip$Annual_Precip)

# Calculate the annual mean temperature for each year
annual_mean_temp <- USOWC_EC_sub %>%
  group_by(year(EC_TIMESTAMP_START)) %>%
  summarize(Annual_Mean_Temp = mean(EC_TA_F, na.rm = TRUE))

# Calculate the overall mean annual temperature across all years
mean(annual_mean_temp$Annual_Mean_Temp)


### US-UAF ###

# read csv

USUAF_EC <- read.csv("path/FLX_US-Uaf_FLUXNET-CH4_HH_2011-2018_1-1.csv")

USUAF_EC$TIMESTAMP_END <- format(USUAF_EC$TIMESTAMP_END, scientific = FALSE)

#create new column for Time2 with seconds
USUAF_EC$TimeStart <- paste(as.character(USUAF_EC$TIMESTAMP_START),"00",sep="")

# convert to datetime

USUAF_EC$TIMESTAMP_START <- as_datetime(USUAF_EC$TimeStart)

#create new column for Time2 with seconds
USUAF_EC$TimeEnd <- paste(as.character(USUAF_EC$TIMESTAMP_END),"00",sep="")

# convert to datetime

USUAF_EC$TIMESTAMP_END <- as_datetime(USUAF_EC$TimeEnd)

Sys.setenv(TZ = "UTC")
# subset EC data to match with chamber data
USUAF_EC_sub <- subset(USUAF_EC, date(TIMESTAMP_START) >= "2016-05-02" & TIMESTAMP_END < "2018-11-03 00:00:30")

# reset row names
rownames(USUAF_EC_sub) <- NULL

# remove some columns that are not needed

USUAF_EC_sub <- subset(USUAF_EC_sub, select=-c(TimeStart, TimeEnd))

# replace -9999 with NA

USUAF_EC_sub <- USUAF_EC_sub %>% na_if(-9999)

# add EC_prefix to all columns
USUAF_EC_sub <- USUAF_EC_sub %>% rename_with( ~ paste0("EC_", .x))

# add Site

USUAF_EC_sub$SITE <- "US-UAF"

# move to front
USUAF_EC_sub <- USUAF_EC_sub %>%
  select(SITE, everything())

# save as .csv

write.csv(USUAF_EC_sub, "path/US_UAF_EC_subset.csv", row.names = FALSE)

USUAF_EC_sub <- read.csv("path/US_UAF_EC_subset.csv")

# calculate mean annual temperature and precipitation sum

# calculate precip sum for each year

# Summarize total annual precipitation
annual_precip <- USUAF_EC_sub %>%
  group_by(year(EC_TIMESTAMP_START)) %>%
  summarize(Annual_Precip = sum(EC_P_F, na.rm = T))

# Calculate the mean annual precipitation
mean(annual_precip$Annual_Precip)

# Calculate the annual mean temperature for each year
annual_mean_temp <- USUAF_EC_sub %>%
  group_by(year(EC_TIMESTAMP_START)) %>%
  summarize(Annual_Mean_Temp = mean(EC_TA_F, na.rm = TRUE))

# Calculate the overall mean annual temperature across all years
mean(annual_mean_temp$Annual_Mean_Temp)


### US-STJ ###

# read csv

USSTJ_EC <- read.csv(
  "path/StJ_FULL_DATABASE_2020a.csv",
  colClasses = c(TIMESTAMP_START = "character",
                 TIMESTAMP_END   = "character")
)

# Parse as month-day-year hour:minute (US format)
USSTJ_EC$TIMESTAMP_START <- mdy_hm(USSTJ_EC$TIMESTAMP_START, tz = "UTC")
USSTJ_EC$TIMESTAMP_END   <- mdy_hm(USSTJ_EC$TIMESTAMP_END,   tz = "UTC")

# Parse DATE_TIME in US format m/d/yy H:M
USSTJ_EC$TIMESTAMP_END <- mdy_hm(USSTJ_EC$DATE_TIME, tz = "UTC")

# For 30-min EC data
USSTJ_EC$TIMESTAMP_START <- USSTJ_EC$TIMESTAMP_END - minutes(30)

# remove extra columns

USSTJ_EC <- subset(USSTJ_EC, select = c(TIMESTAMP_START, TIMESTAMP_END, DATE, DOY, ch4_flux, air_pressure, 
                                        air_temperature, VPD, wind_speed, wind_dir, ch4_scf, Level_YSI, Level_YSI_f, 
                                        Level_NOAA, Tsoil, Ustar, NEE_orig, NEE_f, CH4_orig, CH4_f, CH4_f_RF, FCH4_RF_filled,
                                        Precip_MET))


# Convert units
USSTJ_EC$air_pressure <- USSTJ_EC$air_pressure / 1000         # Pa → kPa
USSTJ_EC$Level_NOAA   <- USSTJ_EC$Level_NOAA   * 100           # m → cm
USSTJ_EC$Level_YSI    <- USSTJ_EC$Level_YSI    * 100           # m → cm
USSTJ_EC$Level_YSI_f  <- USSTJ_EC$Level_YSI_f  * 100           # m → cm


# cut based on chamber data period: 2.5 to 18.12 2020
Sys.setenv(TZ = "UTC")
# subset EC data to match chamber data
USSTJ_EC_sub <- subset(USSTJ_EC, date(TIMESTAMP_START) >= "2020-05-02" & TIMESTAMP_END <= "2020-12-19 00:00:00")

# reset row names
rownames(USSTJ_EC_sub) <- NULL

# rename columns

USSTJ_EC_sub <- USSTJ_EC_sub %>%
  rename(
    PA   = air_pressure,
    TA     = air_temperature,
    WS = wind_speed,
    WD = wind_dir, 
    TS = Tsoil,
    USTAR = Ustar
  )

# TA from Kelvin to Celsius

USSTJ_EC_sub$TA <- USSTJ_EC_sub$TA - 273.15

# add EC_prefix to all columns
USSTJ_EC_sub <- USSTJ_EC_sub %>% rename_with( ~ paste0("EC_", .x))

# add Site

USSTJ_EC_sub$SITE <- "US-STJ"

# move to front
USSTJ_EC_sub <- USSTJ_EC_sub %>%
  select(SITE, everything())

USSTJ_EC_sub <- subset(USSTJ_EC_sub, select = -c(EC_DATE, EC_DOY))

# save as .csv

write.csv(USSTJ_EC_sub, "path/US_STJ_EC_subset.csv", row.names = FALSE)

