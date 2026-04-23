#############################################################
### Code for: 
### Määttä et al. (2026) A cross-site comparison of ecosystem- and plot-scale methane fluxes across multiple sites
### Code created by Tiia Määttä, with parts written with GPT 4, 4o and 5.4

#############################################################
###----------EC AND CHAMBER DATA COMBINATION--------------###
###---------------TEMPORAL AGGREGATIONS-------------------###
#############################################################

# Remove all variables
rm(list=ls(all = TRUE))
# Clear the console window
cat("\014")
# Close all graphics windows
graphics.off()

# open libraries
library("ggplot2")
library("dplyr")
library("tidyr")
library("lubridate")
library("sqldf")

##########################################################################
### EC AND CHAMBER DATA COMBINATION AND TEMPORAL AGGREGATIONS PER SITE ###
##########################################################################

### CN-HGU ###

# read csv

# EC
CNHGU_EC_d <- read.csv("path/CN_HGU_EC_subset_28082025.csv")

# convert to datetime
CNHGU_EC_d$EC_TIMESTAMP_START <- as_datetime(CNHGU_EC_d$EC_TIMESTAMP_START, tz="UTC")
CNHGU_EC_d$EC_TIMESTAMP_END <- as_datetime(CNHGU_EC_d$EC_TIMESTAMP_END, tz="UTC")

# remove the filtering columns
CNHGU_EC_d <- subset(CNHGU_EC_d, select = -c(is_day, is_night, year_month, ok_ustar, extreme_pos, suspect_condensation))

# CHAMBER
CNHGU_ch_d <- read.csv("path/CN_HGU_flux_ST_SM_04042023.csv")

# convert to datetime
CNHGU_ch_d$ch_Datetime <- as_datetime(CNHGU_ch_d$ch_Datetime, tz="UTC")

# combine with SQL query
# set local environment time zone to UTC so the time stamps will match
Sys.setenv(TZ = "UTC")

CNHGU_all <- sqldf('SELECT ch_Datetime, EC_TIMESTAMP_START, EC_TIMESTAMP_END, ch_ID, ch_FCH4_nmolm2s1, EC_FCH4, EC_FCH4_F,
ch_method, ch_LC_class, ch_Collar_type,
EC_NEE, EC_H, EC_LE, EC_USTAR, EC_SW_IN, EC_SW_OUT,EC_LW_IN, EC_LW_OUT, EC_NETRAD, EC_PPFD_IN, EC_VPD, EC_RH, EC_TA, EC_P, EC_TS, EC_G, EC_SWC,
EC_GPP_NT, EC_RECO_NT, EC_GPP_DT, EC_RECO_DT, EC_WD, EC_WS, EC_NEE_F, EC_H_F, EC_LE_F, EC_SW_IN_F, EC_SW_OUT_F, EC_LW_IN_F, EC_LW_OUT_F, EC_NETRAD_F, EC_PPFD_IN_F, EC_VPD_F,
EC_RH_F, EC_PA_F, EC_TA_F, EC_P_F, EC_TS_F, EC_G_F, EC_SWC_F, EC_WS_F, EC_LE_F_ANNOPTLM, EC_NEE_F_ANNOPTLM, EC_FCH4_F_ANNOPTLM,
EC_FCH4_F_RANDUNC, EC_FCH4_F_ANNOPTLM_UNC, EC_FCH4_F_ANNOPTLM_QC
      FROM CNHGU_ch_d 
      LEFT JOIN CNHGU_EC_d ON ch_Datetime BETWEEN EC_TIMESTAMP_START and EC_TIMESTAMP_END')

# add Site column
CNHGU_all$SITE <- "CN-HGU"

# move to the front
CNHGU_all <- CNHGU_all %>%
  dplyr::select(SITE, everything())

# CN-HGU chamber data is already in nmol m2s1

# rename ch_FCH4 to include the unit
colnames(CNHGU_all)[6] <- "ch_FCH4_nmolCH4m2s1"

write.csv(CNHGU_all, "path/CN_HGU_combined.csv", row.names = FALSE)

### TEMPORAL AGGREGATIONS - CN-HGU

CNHGU$EC_TIMESTAMP_START <- as_datetime(CNHGU$EC_TIMESTAMP_START)
CNHGU$EC_TIMESTAMP_END <- as_datetime(CNHGU$EC_TIMESTAMP_END)
CNHGU$ch_Datetime <- as_datetime(CNHGU$ch_Datetime)

# subset to end at EC_TIMESTAMP_END  == 2016-07-31 00:00:00 because the chamber data ends at noon the next day
Sys.setenv(TZ = "UTC")
CNHGU$EC_TIMESTAMP_END <- as_datetime(CNHGU$EC_TIMESTAMP_END)

CNHGU <- subset(CNHGU, EC_TIMESTAMP_END <= "2016-07-31 00:00:00")

# SET EC GAP-FILLED DATA AND CHAMBER DATA TO NA WHEN TIME GAP LONGER THAN 2 MONTHS
# use the 2 months limit

CNHGU$EC_FCH4_F_ANNOPTLM_QC <- as.numeric(CNHGU$EC_FCH4_F_ANNOPTLM_QC)

CNHGU_new <- CNHGU %>% mutate(across(c(ch_FCH4_nmolCH4m2s1, EC_FCH4_F, EC_NEE_F, EC_H_F, EC_LE_F, EC_SW_IN_F, EC_SW_OUT_F,
                                       EC_LW_IN_F, EC_LW_OUT_F, EC_NETRAD_F,
                                       EC_PPFD_IN_F, EC_VPD_F, EC_RH_F, EC_PA_F, EC_TA_F, EC_P_F,
                                       EC_TS_F, EC_G_F, EC_SWC_F, EC_WS_F, EC_LE_F_ANNOPTLM, EC_NEE_F_ANNOPTLM, EC_FCH4_F_ANNOPTLM), ~ ifelse(EC_FCH4_F_ANNOPTLM_QC == 3, NA, .)))

# convert datetimes to datetime format
CNHGU_new$ch_Datetime <- as_datetime(CNHGU_new$ch_Datetime)
CNHGU_new$EC_TIMESTAMP_START <- as_datetime(CNHGU_new$EC_TIMESTAMP_START)
CNHGU_new$EC_TIMESTAMP_END <- as_datetime(CNHGU_new$EC_TIMESTAMP_END)

# set EC data to NA when chamber data is NA
CNHGU_ECfilt <- CNHGU_new %>% mutate(across(c(EC_FCH4, EC_FCH4_F, EC_FCH4_F_ANNOPTLM), 
                                                          ~ ifelse(ch_FCH4_nmolCH4m2s1 == "NA", NA, .)))

CNHGU_ECfilt$ch_FCH4_nmolCH4m2s1[is.na(CNHGU_ECfilt$EC_FCH4_F_ANNOPTLM)] <- NA

# drop rows where FCH4 is NA

CNHGU_ECfilt <- CNHGU_ECfilt %>% drop_na(ch_FCH4_nmolCH4m2s1)

CNHGU_ECfilt <- CNHGU_ECfilt %>% drop_na(EC_FCH4_F_ANNOPTLM)
CNHGU_ECfilt <- CNHGU_ECfilt %>% drop_na(ch_ID)

# HALF-HOURLY AGGREGATION (chamber medians and means per EC timestamp)
CNHGU_hh_aggr <- CNHGU_ECfilt %>%
  group_by(EC_TIMESTAMP_START) %>%
  summarise(
    ch_FCH4_nmolCH4m2s1_mean = mean(ch_FCH4_nmolCH4m2s1, na.rm = TRUE),
    ch_FCH4_nmolCH4m2s1_median = median(ch_FCH4_nmolCH4m2s1, na.rm = TRUE),
    EC_TIMESTAMP_END = first(EC_TIMESTAMP_END),
    across(starts_with("EC_"), first),
    .groups = "drop"
  )

CNHGU_hh_aggr <- as.data.frame(CNHGU_hh_aggr)

# move the new columns before EC_FCH4
CNHGU_hh_aggr <- CNHGU_hh_aggr %>% relocate(ch_FCH4_nmolCH4m2s1_mean, ch_FCH4_nmolCH4m2s1_median, .before = EC_FCH4)

write.csv(CNHGU_hh_aggr, "path/CN_HGU_combined_ch_datetimeaggr.csv", row.names = FALSE)

# HOURLY AGGREGATION

# create new date_hour column for grouping
CNHGU_ECfilt <- CNHGU_ECfilt %>%
  mutate(date_hour = format(EC_TIMESTAMP_START, "%Y-%m-%d %H"))

# set duplicated values to NA in EC columns
# create vector with EC column names
#rename EC_TIMESTAMP_START and END
colnames(CNHGU_ECfilt)[c(3,4)] <- c("TIMESTAMP_START", "TIMESTAMP_END")
cols <- CNHGU_ECfilt %>% dplyr::select(starts_with("EC_")) %>% colnames()

# convert duplicates to NA in EC columns
CNHGU_ECfilt[cols] <- sapply(CNHGU_ECfilt[cols], function(x) 
  ave(x, as_datetime(CNHGU_ECfilt$TIMESTAMP_START), FUN = function(x) replace(x, duplicated(x), NA)))

# convert character and logi to numeric
cols.num <- CNHGU_ECfilt %>% dplyr::select(starts_with("EC_")) %>% colnames()
CNHGU_ECfilt[cols.num] <- sapply(CNHGU_ECfilt[cols.num],as.numeric)

# Calculate the u and v wind components
CNHGU_ECfilt$u.wind <- -CNHGU_ECfilt$EC_WS *sin(2*pi *CNHGU_ECfilt$EC_WD/360)
CNHGU_ECfilt$v.wind <- -CNHGU_ECfilt$EC_WS *cos(2*pi *CNHGU_ECfilt$EC_WD/360)
# gap-filled (no gap-filled WD)
CNHGU_ECfilt$u.wind.F <- -CNHGU_ECfilt$EC_WS_F *sin(2*pi *CNHGU_ECfilt$EC_WD/360)
CNHGU_ECfilt$v.wind.F <- -CNHGU_ECfilt$EC_WS_F *cos(2*pi *CNHGU_ECfilt$EC_WD/360)

# aggregate
CNHGU_hr_aggr_ECfilt <- CNHGU_ECfilt %>% 
  group_by(date_hour) %>%
  summarise_at(c("ch_FCH4_nmolCH4m2s1", "EC_FCH4", "EC_FCH4_F", "EC_FCH4_F_ANNOPTLM",
                 "EC_USTAR", "EC_GPP_NT", "EC_RECO_NT", "EC_GPP_DT", "EC_RECO_DT", "EC_NEE_F", "EC_H_F", "EC_LE_F", "EC_SW_IN_F", "EC_SW_OUT_F", "EC_LW_IN_F", 
                 "EC_LW_OUT_F", "EC_NETRAD_F", "EC_PPFD_IN_F", "EC_VPD_F", "EC_RH_F", 
                 "EC_PA_F", "EC_TA_F", "EC_P_F", "EC_TS_F", "EC_G_F", "EC_SWC_F", "EC_WS_F", 
                 "EC_LE_F_ANNOPTLM", "EC_NEE_F_ANNOPTLM",
                 "u.wind", "v.wind", "u.wind.F", "v.wind.F"), 
               list(mean = mean, 
                    median = median), 
               na.rm = TRUE)

CNHGU_hr_aggr_ECfilt <- as.data.frame(CNHGU_hr_aggr_ECfilt)

# convert to datetime format
CNHGU_hr_aggr_ECfilt$date_hour <- as_datetime(ymd_h(CNHGU_hr_aggr_ECfilt$date_hour))

# calculate wind direction average
CNHGU_hr_aggr_ECfilt$EC_WD_AVG <- (atan2(CNHGU_hr_aggr_ECfilt$u.wind_mean, CNHGU_hr_aggr_ECfilt$v.wind_mean) *360/2/pi) +180

# gap-filled:
CNHGU_hr_aggr_ECfilt$EC_WD_AVG_F <- (atan2(CNHGU_hr_aggr_ECfilt$u.wind.F_mean, CNHGU_hr_aggr_ECfilt$v.wind.F_mean) *360/2/pi) +180

# add the missing extra columns

CNHGU_hr_aggr_ECfilt$Year <- year(CNHGU_hr_aggr_ECfilt$date_hour)
CNHGU_hr_aggr_ECfilt$Month <- month(CNHGU_hr_aggr_ECfilt$date_hour)
CNHGU_hr_aggr_ECfilt$Day <- day(CNHGU_hr_aggr_ECfilt$date_hour)
CNHGU_hr_aggr_ECfilt$Hour <- hour(CNHGU_hr_aggr_ECfilt$date_hour)
CNHGU_hr_aggr_ECfilt$DOY <- yday(CNHGU_hr_aggr_ECfilt$date_hour)

CNHGU_hr_aggr_ECfilt$ch_method <- "auto"
CNHGU_hr_aggr_ECfilt$SITE <- "CN-HGU"

CNHGU_hr_aggr_ECfilt <- CNHGU_hr_aggr_ECfilt %>% relocate(SITE, .before = date_hour)

# rename date_hour
colnames(CNHGU_hr_aggr_ECfilt)[2] <- "TIMESTAMP"

# add dom_veg
CNHGU_hr_aggr_ECfilt$DOM_VEG <- "AERENCHYMATOUS"

# add veg classes
CNHGU_hr_aggr_ECfilt$MOSS_BROWN <- 0
CNHGU_hr_aggr_ECfilt$MOSS_SPHAGNUM <- 0
CNHGU_hr_aggr_ECfilt$AERENCHYMATOUS <- 1
CNHGU_hr_aggr_ECfilt$ERI_SHRUB <- 0
CNHGU_hr_aggr_ECfilt$TREE <- 0

# add climate
CNHGU_hr_aggr_ECfilt$KOPPEN <- "Cwc"

# add ecosystem type
CNHGU_hr_aggr_ECfilt$SITE_CLASSIFICATION <- "upland"

# save as .csv

write.csv(CNHGU_hr_aggr_ECfilt, "path/CN_HGU_combined_filtered_ECfilt_with_outliers_ch_hraggr.csv", row.names = FALSE)

# DAILY AGGREGATION

# aggregate to daily scale

CNHGU_d_aggr_ECfilt <- CNHGU_ECfilt %>% 
  group_by(date(CNHGU_ECfilt$TIMESTAMP_START)) %>%
  summarise_at(c("ch_FCH4_nmolCH4m2s1", "EC_FCH4", "EC_FCH4_F", "EC_FCH4_F_ANNOPTLM",
                 "EC_USTAR", "EC_GPP_NT", "EC_RECO_NT", "EC_GPP_DT", "EC_RECO_DT", "EC_NEE_F", "EC_H_F", "EC_LE_F", "EC_SW_IN_F", "EC_SW_OUT_F", "EC_LW_IN_F", 
                 "EC_LW_OUT_F", "EC_NETRAD_F", "EC_PPFD_IN_F", "EC_VPD_F", "EC_RH_F", 
                 "EC_PA_F", "EC_TA_F", "EC_P_F", "EC_TS_F", "EC_G_F", "EC_SWC_F", "EC_WS_F", 
                 "EC_LE_F_ANNOPTLM", "EC_NEE_F_ANNOPTLM", "EC_FCH4_F_ANNOPTLM_UNC", "ECCH_diff",
                 "u.wind", "v.wind", "u.wind.F", "v.wind.F"), 
               list(mean = mean,
                    median = median), 
               na.rm = TRUE)

CNHGU_d_aggr_ECfilt <- as.data.frame(CNHGU_d_aggr_ECfilt)

# rename the date column

colnames(CNHGU_d_aggr_ECfilt)[1] <- "DATE"

# convert to datetime format
CNHGU_d_aggr_ECfilt$DATE <- as_date(CNHGU_d_aggr_ECfilt$DATE)

# calculate wind direction average
CNHGU_d_aggr_ECfilt$EC_WD_AVG <- (atan2(CNHGU_d_aggr_ECfilt$u.wind_mean, CNHGU_d_aggr_ECfilt$v.wind_mean) *360/2/pi) +180

# # gap-filled:
CNHGU_d_aggr_ECfilt$EC_WD_AVG_F <- (atan2(CNHGU_d_aggr_ECfilt$u.wind.F_mean, CNHGU_d_aggr_ECfilt$v.wind.F_mean) *360/2/pi) +180

# add the missing extra columns

CNHGU_d_aggr_ECfilt$Year <- year(CNHGU_d_aggr_ECfilt$DATE)
CNHGU_d_aggr_ECfilt$Month <- month(CNHGU_d_aggr_ECfilt$DATE)
CNHGU_d_aggr_ECfilt$Day <- day(CNHGU_d_aggr_ECfilt$DATE)
CNHGU_d_aggr_ECfilt$DOY <- yday(CNHGU_d_aggr_ECfilt$DATE)

CNHGU_d_aggr_ECfilt$ch_method <- "auto"
CNHGU_d_aggr_ECfilt$SITE <- "CN-HGU"

CNHGU_d_aggr_ECfilt <- CNHGU_d_aggr_ECfilt %>% relocate(SITE, .before = DATE)

# add dom_veg
CNHGU_d_aggr_ECfilt$DOM_VEG <- "AERENCHYMATOUS"

# add veg classes
CNHGU_d_aggr_ECfilt$MOSS_BROWN <- 0
CNHGU_d_aggr_ECfilt$MOSS_SPHAGNUM <- 0
CNHGU_d_aggr_ECfilt$AERENCHYMATOUS <- 1
CNHGU_d_aggr_ECfilt$ERI_SHRUB <- 0
CNHGU_d_aggr_ECfilt$TREE <- 0

# add climate
CNHGU_d_aggr_ECfilt$KOPPEN <- "Cwc"

# add ecosystem type
CNHGU_d_aggr_ECfilt$SITE_CLASSIFICATION <- "upland"

# save as .csv

write.csv(CNHGU_d_aggr_ECfilt, "path/CN_HGU_chec_daggr_ECfilt.csv", row.names = FALSE)

# WEEKLY AGGREGATION

CNHGU_week_aggr_ECfilt <- CNHGU_ECfilt %>% 
  group_by(year(CNHGU_ECfilt$TIMESTAMP_START), isoweek(CNHGU_ECfilt$TIMESTAMP_START)) %>%
  summarise_at(c("ch_FCH4_nmolCH4m2s1", "EC_FCH4", "EC_FCH4_F", "EC_FCH4_F_ANNOPTLM",
                 "EC_USTAR", "EC_GPP_NT", "EC_RECO_NT", "EC_GPP_DT", "EC_RECO_DT", "EC_NEE_F", "EC_H_F", "EC_LE_F", "EC_SW_IN_F", "EC_SW_OUT_F", "EC_LW_IN_F", 
                 "EC_LW_OUT_F", "EC_NETRAD_F", "EC_PPFD_IN_F", "EC_VPD_F", "EC_RH_F", 
                 "EC_PA_F", "EC_TA_F", "EC_P_F", "EC_TS_F", "EC_G_F", "EC_SWC_F", "EC_WS_F", 
                 "EC_LE_F_ANNOPTLM", "EC_NEE_F_ANNOPTLM",  
                 "u.wind", "v.wind", "u.wind.F", "v.wind.F"), 
               list(mean = mean,
                    median = median), 
               na.rm = TRUE)

CNHGU_week_aggr_ECfilt <- as.data.frame(CNHGU_week_aggr_ECfilt)

# rename the date column

colnames(CNHGU_week_aggr_ECfilt)[1] <- "Year"
colnames(CNHGU_week_aggr_ECfilt)[2] <- "Week_of_year"

# calculate wind direction average
CNHGU_week_aggr_ECfilt$EC_WD_AVG <- (atan2(CNHGU_week_aggr_ECfilt$u.wind_mean, CNHGU_week_aggr_ECfilt$v.wind_mean) *360/2/pi) +180

# gap-filled:
CNHGU_week_aggr_ECfilt$EC_WD_AVG_F <- (atan2(CNHGU_week_aggr_ECfilt$u.wind.F_mean, CNHGU_week_aggr_ECfilt$v.wind.F_mean) *360/2/pi) +180

# add the missing extra columns

CNHGU_week_aggr_ECfilt$ch_method <- "auto"
CNHGU_week_aggr_ECfilt$SITE <- "CN-HGU"

CNHGU_week_aggr_ECfilt <- CNHGU_week_aggr_ECfilt %>% relocate(SITE, .before = Year)

# add dom_veg
CNHGU_week_aggr_ECfilt$DOM_VEG <- "AERENCHYMATOUS"

# add veg classes
CNHGU_week_aggr_ECfilt$MOSS_BROWN <- 0
CNHGU_week_aggr_ECfilt$MOSS_SPHAGNUM <- 0
CNHGU_week_aggr_ECfilt$AERENCHYMATOUS <- 1
CNHGU_week_aggr_ECfilt$ERI_SHRUB <- 0
CNHGU_week_aggr_ECfilt$TREE <- 0

# add climate
CNHGU_week_aggr_ECfilt$KOPPEN <- "Cwc"

# add ecosystem type
CNHGU_week_aggr_ECfilt$SITE_CLASSIFICATION <- "upland"

# add month

# Extract the correct month directly using lubridate::make_date and ISO week (1st day of the week is Monday)
CNHGU_week_aggr_ECfilt$Month <- month(make_date(CNHGU_week_aggr_ECfilt$Year) + weeks(CNHGU_week_aggr_ECfilt$Week_of_year))

CNHGU_week_aggr_ECfilt <- CNHGU_week_aggr_ECfilt %>% relocate(Month, .after = Week_of_year)

# multiple weeks and months need to be fixed manually
CNHGU_week_aggr_ECfilt$Month[CNHGU_week_aggr_ECfilt$Week_of_year == 35 & CNHGU_week_aggr_ECfilt$Year == 2015] <- 8
CNHGU_week_aggr_ECfilt$Month[CNHGU_week_aggr_ECfilt$Week_of_year == 39 & CNHGU_week_aggr_ECfilt$Year == 2015] <- 9
CNHGU_week_aggr_ECfilt$Month[CNHGU_week_aggr_ECfilt$Week_of_year == 40 & CNHGU_week_aggr_ECfilt$Year == 2015] <- 9
CNHGU_week_aggr_ECfilt$Month[CNHGU_week_aggr_ECfilt$Week_of_year == 13 & CNHGU_week_aggr_ECfilt$Year == 2016] <- 3
CNHGU_week_aggr_ECfilt$Month[CNHGU_week_aggr_ECfilt$Week_of_year == 26 & CNHGU_week_aggr_ECfilt$Year == 2015] <- 6

# save as .csv

write.csv(CNHGU_week_aggr_ECfilt, "path/CN_HGU_chec_weekaggr_ECfilt.csv", row.names = FALSE)

# MONTHLY AGGREGATION

CNHGU_month_aggr_ECfilt <- CNHGU_ECfilt %>% 
  group_by(year(CNHGU_ECfilt$TIMESTAMP_START), month(CNHGU_ECfilt$TIMESTAMP_START)) %>%
  summarise_at(c("ch_FCH4_nmolCH4m2s1", "EC_FCH4", "EC_FCH4_F", "EC_FCH4_F_ANNOPTLM",
                 "EC_USTAR", "EC_GPP_NT", "EC_RECO_NT", "EC_GPP_DT", "EC_RECO_DT", "EC_NEE_F", "EC_H_F", "EC_LE_F", "EC_SW_IN_F", "EC_SW_OUT_F", "EC_LW_IN_F", 
                 "EC_LW_OUT_F", "EC_NETRAD_F", "EC_PPFD_IN_F", "EC_VPD_F", "EC_RH_F", 
                 "EC_PA_F", "EC_TA_F", "EC_P_F", "EC_TS_F", "EC_G_F", "EC_SWC_F", "EC_WS_F", 
                 "EC_LE_F_ANNOPTLM", "EC_NEE_F_ANNOPTLM", "EC_FCH4_F_ANNOPTLM_UNC", 
                 "u.wind", "v.wind", "u.wind.F", "v.wind.F"), 
               list(mean = mean,
                    median = median), 
               na.rm = TRUE)

CNHGU_month_aggr_ECfilt <- as.data.frame(CNHGU_month_aggr_ECfilt)

# rename the date column

colnames(CNHGU_month_aggr_ECfilt)[1] <- "TIMESTAMP"

# calculate wind direction average
CNHGU_month_aggr_ECfilt$EC_WD_AVG <- (atan2(CNHGU_month_aggr_ECfilt$u.wind_mean, CNHGU_month_aggr_ECfilt$v.wind_mean) *360/2/pi) +180

# gap-filled:
CNHGU_month_aggr_ECfilt$EC_WD_AVG_F <- (atan2(CNHGU_month_aggr_ECfilt$u.wind.F_mean, CNHGU_month_aggr_ECfilt$v.wind.F_mean) *360/2/pi) +180

# add the missing extra columns

CNHGU_month_aggr_ECfilt$ch_method <- "auto"
CNHGU_month_aggr_ECfilt$SITE <- "CN-HGU"

CNHGU_month_aggr_ECfilt <- CNHGU_month_aggr_ECfilt %>% relocate(SITE, .before = TIMESTAMP)

# add dom_veg
CNHGU_month_aggr_ECfilt$DOM_VEG <- "AERENCHYMATOUS"

# add veg classes
CNHGU_month_aggr_ECfilt$MOSS_BROWN <- 0
CNHGU_month_aggr_ECfilt$MOSS_SPHAGNUM <- 0
CNHGU_month_aggr_ECfilt$AERENCHYMATOUS <- 1
CNHGU_month_aggr_ECfilt$ERI_SHRUB <- 0
CNHGU_month_aggr_ECfilt$TREE <- 0

# add climate
CNHGU_month_aggr_ECfilt$KOPPEN <- "Cwc"

# add ecosystem type
CNHGU_month_aggr_ECfilt$SITE_CLASSIFICATION <- "upland"

colnames(CNHGU_month_aggr_ECfilt)[3] <- "Month"

# save as .csv

write.csv(CNHGU_month_aggr_ECfilt, "path/CN_HGU_chec_monthnaggr_ECfilt.csv", row.names = FALSE)

# ANNUAL AGGREGATION

# aggregate
CNHGU_yr_aggr_ECfilt <- CNHGU_ECfilt %>% 
  group_by(year(CNHGU_ECfilt$TIMESTAMP_START)) %>%
  summarise_at(c("ch_FCH4_nmolCH4m2s1", "EC_FCH4", "EC_FCH4_F", "EC_FCH4_F_ANNOPTLM",
                 "EC_USTAR", "EC_GPP_NT", "EC_RECO_NT", "EC_GPP_DT", "EC_RECO_DT", "EC_NEE_F", "EC_H_F", "EC_LE_F", "EC_SW_IN_F", "EC_SW_OUT_F", "EC_LW_IN_F", 
                 "EC_LW_OUT_F", "EC_NETRAD_F", "EC_PPFD_IN_F", "EC_VPD_F", "EC_RH_F", 
                 "EC_PA_F", "EC_TA_F", "EC_P_F", "EC_TS_F", "EC_G_F", "EC_SWC_F", "EC_WS_F", 
                 "EC_LE_F_ANNOPTLM", "EC_NEE_F_ANNOPTLM", "EC_FCH4_F_ANNOPTLM_UNC", 
                 "u.wind", "v.wind", "u.wind.F", "v.wind.F"), 
               list(mean = mean,
                    median = median), 
               na.rm = TRUE)

CNHGU_yr_aggr_ECfilt <- as.data.frame(CNHGU_yr_aggr_ECfilt)

# rename the date column

colnames(CNHGU_yr_aggr_ECfilt)[1] <- "TIMESTAMP"

# calculate wind direction average
CNHGU_yr_aggr_ECfilt$EC_WD_AVG <- (atan2(CNHGU_yr_aggr_ECfilt$u.wind_mean, CNHGU_yr_aggr_ECfilt$v.wind_mean) *360/2/pi) +180

# gap-filled:
CNHGU_yr_aggr_ECfilt$EC_WD_AVG_F <- (atan2(CNHGU_yr_aggr_ECfilt$u.wind.F_mean, CNHGU_yr_aggr_ECfilt$v.wind.F_mean) *360/2/pi) +180

# add the missing extra columns

CNHGU_yr_aggr_ECfilt$ch_method <- "auto"
CNHGU_yr_aggr_ECfilt$SITE <- "CN-HGU"

CNHGU_yr_aggr_ECfilt <- CNHGU_yr_aggr_ECfilt %>% relocate(SITE, .before = TIMESTAMP)

# add dom_veg
CNHGU_yr_aggr_ECfilt$DOM_VEG <- "AERENCHYMATOUS"

# add veg classes
CNHGU_yr_aggr_ECfilt$MOSS_BROWN <- 0
CNHGU_yr_aggr_ECfilt$MOSS_SPHAGNUM <- 0
CNHGU_yr_aggr_ECfilt$AERENCHYMATOUS <- 1
CNHGU_yr_aggr_ECfilt$ERI_SHRUB <- 0
CNHGU_yr_aggr_ECfilt$TREE <- 0

# add climate
CNHGU_yr_aggr_ECfilt$KOPPEN <- "Cwc"

# add ecosystem type
CNHGU_yr_aggr_ECfilt$SITE_CLASSIFICATION <- "upland"

# save as .csv

write.csv(CNHGU_yr_aggr_ECfilt, "path/CN_HGU_chec_yraggr_ECfilt.csv", row.names = FALSE)



### FI-Si2 ###

# read csv

# EC
FISI2_EC_d <- read.csv("path/FI_Si2_EC_subset_HH_TO_DD_AGGR_ECfilt.csv")

FISI2_EC_week <- read.csv("path/FI_Si2_EC_subset_HH_TO_WEEK_AGGR_ECfilt.csv")

FISI2_EC_month <- read.csv("path/FI_Si2_EC_subset_HH_TO_MONTH_AGGR_ECfilt.csv")

FISI2_EC_year <- read.csv("path/FI_Si2_EC_subset_HH_TO_YR_AGGR_ECfilt.csv")

# rename columns
FISI2_EC_month <- FISI2_EC_month %>% rename("u.wind_mean" = "EC_WD_u.wind_mean",
                                            "u.wind.F_mean" = "EC_WD_u.wind_mean_F",
                                            "v.wind_mean" = "EC_WD_v.wind_mean",
                                            "v.wind.F_mean" = "EC_WD_v.wind_mean_F"
)

FISI2_EC_month <- FISI2_EC_month %>% rename("v.wind.F_mean" = "v.wind_mean_F"
)

FISI2_EC_d <- FISI2_EC_d %>% rename("u.wind_mean" = "EC_WD_u.wind_mean",
                                "u.wind.F_mean" = "EC_WD_u.wind_mean_F",
                                "v.wind_mean" = "EC_WD_v.wind_mean",
                                "v.wind.F_mean" = "EC_WD_v.wind_mean_F"
)


# CHAMBER
FISI2_ch <- read.csv("path/FI_Si2.csv")

# convert to date format
FISI2_ch$Date <- as_date(FISI2_ch$Date)

# remove year, month, day, hour
FISI2_ch <- subset(FISI2_ch, select=-c(Time, Year, Month, Day, Hour))

# add ch_ prefix to chamber data
FISI2_ch <- FISI2_ch %>% rename_with( ~ paste0("ch_", .x))

colnames(FISI2_ch)[1] <- "SITE"

# CONVERT CHAMBER DATA TO nmol CH4 m-2 s-1 

# unit is mg m-2 day-1

# convert from day to second
FISI2_ch$ch_FCH4_mgCH4m2s1 <- FISI2_ch$ch_FCH4_mgCH4m2d1/86400

# convert from mg m2s1 to nmol m2s1
FISI2_ch$ch_FCH4_nmolCH4m2s1 <- (FISI2_ch$ch_FCH4_mgCH4m2s1/16.04)*1000000

# move next to mg
FISI2_ch <- FISI2_ch %>% relocate(c(ch_FCH4_nmolCH4m2s1), .after = ch_FCH4_mgCH4m2d1)

# remove the extra column
FISI2_ch <- subset(FISI2_ch, select=-c(ch_FCH4_mgCH4m2s1))

# rename the Date column
colnames(FISI2_ch)[2] <- "TIMESTAMP"

# combine chamber and EC
FISI2_ECall <- merge(FISI2_ch,FISI2_EC_d, by = "TIMESTAMP", all=TRUE, sort = FALSE) 

FISI2_ECall <- subset(FISI2_ECall, select = -c(SITE.x, SITE.y, ch_method.x, ch_method.y))

FISI2_ECall$ch_method <- "manual"
FISI2_ECall$SITE <- "FI-SI2"

FISI2_ECall <- FISI2_ECall %>% relocate(SITE, .before = TIMESTAMP)

# drop rows were ch is NA

FISI2_ECfilt <- FISI2_ECall %>% mutate(across(c(EC_FCH4_mean, EC_FCH4_median, EC_FCH4_sd, EC_FCH4_IQR, 
                                                EC_FCH4_F_mean, EC_FCH4_F_median, EC_FCH4_F_sd, EC_FCH4_F_IQR, 
                                                EC_FCH4_F_ANNOPTLM_mean, EC_FCH4_F_ANNOPTLM_median, EC_FCH4_F_ANNOPTLM_sd, EC_FCH4_F_ANNOPTLM_IQR), 
                                              ~ ifelse(ch_FCH4_nmolCH4m2s1 == "NA", NA, .)))

# save as .csv
write.csv(FISI2_ECfilt, "path/FI_SI2_combined_ECfilt.csv", row.names = FALSE)
FISI2_ECfilt <- read.csv("path/FI_SI2_combined_ECfilt.csv")

# remove the DOY etc columns
FISI2_ECfilt <- subset(FISI2_ECfilt, select = -c(DOY.x, DOY.y, Year.x, Year.y, Month.y, Month.x))
FISI2_ECfilt <- subset(FISI2_ECfilt, select = -c(Day.x, Day.y))

# DAILY AGGREGATION

# aggregate to daily level
FISI2_ch_d_aggr <- FISI2_ECfilt %>% 
  group_by(date(FISI2_ECfilt$TIMESTAMP)) %>%
  summarise_at(c("ch_FCH4_nmolCH4m2s1", "ch_WT", "ch_TS_2", "ch_TS_5", "ch_TS_10", "ch_TS_15", "ch_TS_30",
                 "ch_AerLAI", "ECCH_diff"), 
               list(mean = mean,
                    median = median), 
               na.rm = TRUE)

FISI2_ch_d_aggr <- as.data.frame(FISI2_ch_d_aggr)

# rename the date column

colnames(FISI2_ch_d_aggr)[1] <- "TIMESTAMP"

# convert to datetime format
FISI2_ch_d_aggr$TIMESTAMP <- as_date(FISI2_ch_d_aggr$TIMESTAMP)

# add the missing columns

FISI2_ch_d_aggr$ch_method <- "manual"
FISI2_ch_d_aggr$SITE <- "FI-SI2"

FISI2_ch_d_aggr <- FISI2_ch_d_aggr %>% relocate(SITE, .before = TIMESTAMP)

# remove the last row from 24.9.2014
FISI2_ch_d_aggr <- subset(FISI2_ch_d_aggr, FISI2_ch_d_aggr$TIMESTAMP != "2014-09-24")

# combine with the dd-aggregated EC data

FISI2_d_aggr <- merge(FISI2_ch_d_aggr, FISI2_EC_d, by=c("SITE", "TIMESTAMP", "ch_method"),all.x = TRUE, all.y=TRUE) 

# drop na

FISI2_d_aggr <- FISI2_d_aggr %>% drop_na(ch_FCH4_nmolCH4m2s1_mean)

# add dom_veg
FISI2_d_aggr$DOM_VEG <- "MOSS_SPHAGNUM"

# add veg classes
FISI2_d_aggr$MOSS_BROWN <- 0
FISI2_d_aggr$MOSS_SPHAGNUM <- 1
FISI2_d_aggr$AERENCHYMATOUS <- 1
FISI2_d_aggr$ERI_SHRUB <- 1
FISI2_d_aggr$TREE <- 1

# add climate
FISI2_d_aggr_ECfilt$KOPPEN <- "Dfc"

# add ecosystem type
FISI2_d_aggr_ECfilt$SITE_CLASSIFICATION <- "bog"

# save as. csv
write.csv(FISI2_d_aggr, "path/FI_SI2_chec_daggr.csv", row.names = FALSE)


# WEEKLY AGGREGATION

# aggregate to daily level

FISI2_ECfilt$TIMESTAMP <- as_date(FISI2_ECfilt$TIMESTAMP)

FISI2_ch_week_aggr <- FISI2_ECfilt %>% 
  group_by(year(FISI2_ECfilt$TIMESTAMP), isoweek(FISI2_ECfilt$TIMESTAMP)) %>%
  summarise_at(c("ch_FCH4_nmolCH4m2s1", "ch_WT", "ch_TS_5", "ch_TS_15", "ch_TS_30",
                 "ch_AerLAI"), 
               list(mean = mean,
                    median = median), 
               na.rm = TRUE)

FISI2_ch_week_aggr <- as.data.frame(FISI2_ch_week_aggr)

colnames(FISI2_ch_week_aggr)[1] <- "Year"
colnames(FISI2_ch_week_aggr)[2] <- "Week_of_year"

# drop NA

FISI2_ch_week_aggr <- FISI2_ch_week_aggr %>% drop_na(ch_FCH4_nmolCH4m2s1_median)

# add the missing columns

# add month

FISI2_ch_week_aggr$month <- month(as.Date(paste(FISI2_ch_week_aggr$Year, FISI2_ch_week_aggr$Week_of_year, 1, sep = "-"), format = "%Y-%U-%u"))

FISI2_ch_week_aggr$ch_method <- "manual"
FISI2_ch_week_aggr$SITE <- "FI-SI2"

FISI2_ch_week_aggr <- FISI2_ch_week_aggr %>% relocate(SITE, .before = Year)

# combine with the new week-aggregated EC data

FISI2_week_aggr <- merge(FISI2_ch_week_aggr, FISI2_EC_week, by=c("SITE", "Year", "Week_of_year", "ch_method"), all.x = TRUE, all.y=TRUE) 

# drop na
FISI2_week_aggr <- FISI2_week_aggr %>% drop_na(ch_FCH4_nmolCH4m2s1_mean)

# add dom_veg
FISI2_week_aggr$DOM_VEG <- "MOSS_SPHAGNUM"

# add veg classes
FISI2_week_aggr$MOSS_BROWN <- 0
FISI2_week_aggr$MOSS_SPHAGNUM <- 1
FISI2_week_aggr$AERENCHYMATOUS <- 1
FISI2_week_aggr$ERI_SHRUB <- 1
FISI2_week_aggr$TREE <- 1

# add climate
FISI2_week_aggr$KOPPEN <- "Dfc"

# add ecosystem type
FISI2_week_aggr$SITE_CLASSIFICATION <- "bog"

# rename month column

colnames(FISI2_week_aggr)[17] <- "Month"

FISI2_week_aggr <- FISI2_week_aggr %>% relocate(Month, .after = Week_of_year)

# fix wrong weeks and months
FISI2_week_aggr$Month[FISI2_week_aggr$Week_of_year == 22 & FISI2_week_aggr$Year == 2013] <- 5
FISI2_week_aggr$Month[FISI2_week_aggr$Week_of_year == 31 & FISI2_week_aggr$Year == 2013] <- 7

FISI2_week_aggr$Week_of_year[FISI2_week_aggr$Week_of_year == 28 & FISI2_week_aggr$Year == 2014] <- 29
FISI2_week_aggr$Week_of_year[FISI2_week_aggr$Week_of_year == 30 & FISI2_week_aggr$Year == 2014] <- 31
FISI2_week_aggr$Week_of_year[FISI2_week_aggr$Week_of_year == 34 & FISI2_week_aggr$Year == 2014] <- 35

# rename the u and v component columns
# rename the u and v component columns
FISI2_week_aggr <- FISI2_week_aggr %>% 
  rename(
    "u.wind_mean"   = "EC_WD_u.wind_mean",
    "u.wind.F_mean" = "EC_WD_u.wind_mean_F",
    "v.wind_mean"   = "EC_WD_v.wind_mean",
    "v.wind.F_mean" = "EC_WD_v.wind_mean_F"
  )

# save as .csv
write.csv(FISI2_week_aggr, "path/FI_SI2_chec_weekaggr_ECfilt.csv", row.names = FALSE)

# MONTHLY AGGREGATION

FISI2_ch_month_aggr_ECfilt <- FISI2_ECfilt %>% 
  group_by(year(FISI2_ECfilt$TIMESTAMP), month(FISI2_ECfilt$TIMESTAMP)) %>%
  summarise_at(c("ch_FCH4_nmolCH4m2s1", "ch_WT", "ch_TS_2", "ch_TS_5", "ch_TS_10", "ch_TS_15", "ch_TS_30",
                 "ch_AerLAI", "ECCH_diff"), 
               list(mean = mean,
                    median = median), 
               na.rm = TRUE)

FISI2_ch_month_aggr_ECfilt <- as.data.frame(FISI2_ch_month_aggr_ECfilt)

colnames(FISI2_ch_month_aggr_ECfilt)[1] <- "Year"
colnames(FISI2_ch_month_aggr_ECfilt)[2] <- "Month"

# add the missing columns

FISI2_ch_month_aggr_ECfilt$ch_method <- "manual"
FISI2_ch_month_aggr_ECfilt$SITE <- "FI-SI2"

FISI2_ch_month_aggr_ECfilt <- FISI2_ch_month_aggr_ECfilt %>% relocate(SITE, .before = Year)

colnames(FISI2_EC_month)[1] <- "Year"
colnames(FISI2_EC_month)[2] <- "Month"

# combine chamber and EC
FISI2_month_aggr_ECfilt <- merge(FISI2_ch_month_aggr_ECfilt, FISI2_EC_month, by=c("SITE", "Year", "Month" ,"ch_method"), all.x = TRUE, all.y=TRUE) 

# drop the NA rows
FISI2_month_aggr_ECfilt <- FISI2_month_aggr_ECfilt %>% drop_na(ch_FCH4_nmolCH4m2s1_median)

# add dom_veg
FISI2_month_aggr_ECfilt$DOM_VEG <- "MOSS_SPHAGNUM"

# add veg classes
FISI2_month_aggr_ECfilt$MOSS_BROWN <- 0
FISI2_month_aggr_ECfilt$MOSS_SPHAGNUM <- 1
FISI2_month_aggr_ECfilt$AERENCHYMATOUS <- 1
FISI2_month_aggr_ECfilt$ERI_SHRUB <- 1
FISI2_month_aggr_ECfilt$TREE <- 1

# add climate
FISI2_month_aggr_ECfilt$KOPPEN <- "Dfc"

# add ecosystem type
FISI2_month_aggr_ECfilt$SITE_CLASSIFICATION <- "bog"

# save as .csv
write.csv(FISI2_month_aggr_ECfilt, "path/FI_SI2_chec_monthaggr.csv", row.names = FALSE)


# ANNUAL AGGREGATION

FISI2_ch_yr_aggr <- FISI2_ECfilt %>% 
  group_by(year(FISI2_ECfilt$TIMESTAMP)) %>%
  summarise_at(c("ch_FCH4_nmolCH4m2s1", "ch_WT", "ch_TS_2", "ch_TS_5", "ch_TS_10", "ch_TS_15", "ch_TS_30",
                 "ch_AerLAI"), 
               list(mean = mean,
                    median = median), 
               na.rm = TRUE)

FISI2_ch_yr_aggr <- as.data.frame(FISI2_ch_yr_aggr)

colnames(FISI2_ch_yr_aggr)[1] <- "Year"

# add the missing columns

FISI2_ch_yr_aggr$ch_method <- "manual"
FISI2_ch_yr_aggr$SITE <- "FI-SI2"
FISI2_ch_yr_aggr <- FISI2_ch_yr_aggr %>% relocate(SITE, .before = Year)

FISI2_yr_aggr <- merge(FISI2_ch_yr_aggr, FISI2_EC_year, by=c("SITE", "Year", "ch_method"),all.x = TRUE, all.y=TRUE) 

# add dom_veg
FISI2_yr_aggr$DOM_VEG <- "MOSS_SPHAGNUM"

# add veg classes
FISI2_yr_aggr$MOSS_BROWN <- 0
FISI2_yr_aggr$MOSS_SPHAGNUM <- 1
FISI2_yr_aggr$AERENCHYMATOUS <- 1
FISI2_yr_aggr$ERI_SHRUB <- 1
FISI2_yr_aggr$TREE <- 1

# add climate
FISI2_yr_aggr$KOPPEN <- "Dfc"

# add ecosystem type
FISI2_yr_aggr$SITE_CLASSIFICATION <- "bog"

# save as .csv
write.csv(FISI2_yr_aggr_ECfilt, "path/FI_SI2_chec_yraggr.csv", row.names = FALSE)


### SE-DEG ###

# read csv

# EC
SEDEG_EC_d <- read.csv("path/SE_DEG_EC_subset.csv")

# convert to datetime
SEDEG_EC_d$EC_TIMESTAMP_START <- as_datetime(SEDEG_EC_d$EC_TIMESTAMP_START, tz = "UTC")
SEDEG_EC_d$EC_TIMESTAMP_END <- as_datetime(SEDEG_EC_d$EC_TIMESTAMP_END, tz = "UTC")

# CHAMBER
SEDEG_ch_d <- read.csv("path/SE_DEG.csv")

# convert to datetime
SEDEG_ch_d$ch_Datetime <- as_datetime(SEDEG_ch_d$ch_Datetime, tz = "UTC")

# combine the data frames using SQL:

# set local environment time zone to UTC so the time stamps will match
Sys.setenv(TZ = "UTC")

SEDEG_all <- sqldf('SELECT ch_Datetime, EC_TIMESTAMP_START, EC_TIMESTAMP_END, ch_ID, ch_FCH4_umolCH4m2s1, EC_FCH4, EC_FCH4_F,
ch_CH4_PAR_mean, ch_CH4_PAR_amb, ch_CH4_Ta_amb, ch_TS_2, ch_TS_5, ch_TS_10, ch_TS_15, ch_TS_30, ch_WTL, ch_PPT, ch_WS, ch_P, ch_VPD, ch_RH, ch_method, ch_LC_class,
EC_H, EC_LE, EC_USTAR, EC_SW_IN, EC_NETRAD, EC_PPFD_IN, EC_VPD, EC_RH, EC_TA, EC_PA, EC_P, EC_TS_1, EC_TS_2, EC_TS_3, EC_TS_4, EC_TS_5, EC_TS_6,
EC_WTD, EC_GPP_NT, EC_RECO_NT, EC_GPP_DT, EC_RECO_DT, EC_SWC_1, EC_WS, EC_WD, EC_NEE_F, EC_H_F, EC_LE_F, EC_SW_IN_F, EC_LW_IN_F, EC_NETRAD_F, EC_PPFD_IN_F, EC_VPD_F,
EC_RH_F, EC_PA_F, EC_TA_F, EC_WTD_F, EC_P_F, EC_WS_F, EC_LE_F_ANNOPTLM, EC_NEE_F_ANNOPTLM, EC_FCH4_F_ANNOPTLM,
EC_FCH4_F_RANDUNC, EC_FCH4_F_ANNOPTLM_UNC, EC_FCH4_F_ANNOPTLM_QC
      FROM SEDEG_ch_d 
      LEFT JOIN SEDEG_EC_d ON ch_Datetime BETWEEN EC_TIMESTAMP_START and EC_TIMESTAMP_END')

# add Site column
SEDEG_all$SITE <- "SE-DEG"

# move to the front
SEDEG_all <- SEDEG_all %>%
  select(SITE, everything())

# CONVERT CHAMBER DATA TO nmol CH4 m-2 s-1 

# in umol m2 s1 

# multiply umol by 1000 to get nmol m2 s1
SEDEG_all$ch_FCH4_nmolCH4m2s1 <- SEDEG_all$ch_FCH4_umolCH4m2s1*1000

# move next to umol
SEDEG_all <- SEDEG_all %>% relocate(ch_FCH4_nmolCH4m2s1, .after = ch_FCH4_umolCH4m2s1)

# remove ch_FCH4_umolCH4m2s1
SEDEG_all_new <- subset(SEDEG_all, select=-c(ch_FCH4_umolCH4m2s1))

# with converted units:
write.csv(SEDEG_all_new, "path/SE_DEG_combined.csv", row.names = FALSE)

# read in again
SEDEG <- read.csv("path/SE_DEG_combined.csv")

# SET EC GAP-FILLED DATA AND CHAMBER DATA TO NA WHEN TIME GAP LONGER THAN 2 MONTHS

SEDEG_new <- SEDEG %>% mutate(across(c(ch_FCH4_nmolCH4m2s1, 
                                              ch_CH4_PAR_mean, ch_CH4_PAR_amb, ch_CH4_Ta_amb, ch_TS_2, ch_TS_5, 
                                              ch_TS_10, ch_TS_15, ch_TS_30, ch_WTL, ch_PPT,
                                              ch_WS, ch_P, ch_VPD, ch_RH,
                                              EC_FCH4_F, EC_NEE_F, EC_H_F, EC_LE_F, EC_SW_IN_F, EC_LW_IN_F, EC_NETRAD_F, EC_PPFD_IN_F,
                                              EC_VPD_F, EC_RH_F, EC_PA_F, EC_TA_F, EC_P_F,EC_WTD_F,
                                              EC_WS_F, EC_LE_F_ANNOPTLM, EC_NEE_F_ANNOPTLM, EC_FCH4_F_ANNOPTLM, 
                                              EC_FCH4_F_RANDUNC, EC_FCH4_F_ANNOPTLM_UNC), ~ ifelse(EC_FCH4_F_ANNOPTLM_QC == 3, NA, .)))



# convert datetimes to datetime format
SEDEG_new$EC_TIMESTAMP_START <- as_datetime(SEDEG_new$EC_TIMESTAMP_START)
SEDEG_new$EC_TIMESTAMP_END <- as_datetime(SEDEG_new$EC_TIMESTAMP_END)

# move EC_FCH4_F_ANNOPTLM next to other flux columns
SEDEG_new <- SEDEG_new %>% relocate(EC_FCH4_F_ANNOPTLM, .after = EC_FCH4_F)

# subset to full days
Sys.setenv(TZ = "UTC")

SEDEG_new2 <- subset(SEDEG_new, EC_TIMESTAMP_START >= "2015-05-14 00:00:00" & EC_TIMESTAMP_END <= "2016-11-02 00:00:00")

SEDEG_ECfilt <- SEDEG_new2 %>% mutate(across(c(EC_FCH4, EC_FCH4_F, EC_FCH4_F_ANNOPTLM), 
                                                                          ~ ifelse(ch_FCH4_nmolCH4m2s1 == "NA", NA, .)))

# HALF-HOURLY AGGREGATION

SEDEG_hh_aggr <- SEDEG_ECfilt %>% 
  group_by(across(-c(ch_Datetime, ch_ID, ch_FCH4_nmolCH4m2s1,
                     ch_TS_2, ch_TS_10, ch_WTL))) %>%
  summarise_at(c("ch_FCH4_nmolCH4m2s1", "ch_TS_2", "ch_TS_10",
                 "ch_WTL", "ch_BM_AERENCHYMATOUS", "ch_BM_DWARFSHRUB", "ch_GAI_AERENCHYMATOUS", 
                 "ch_GAI_DWARFSHRUB", "ch_COVER_SPHAGNUM", "ECCH_diff"), 
               list(mean = mean,
                    median = median), 
               na.rm = TRUE)

SEDEG_hh_aggr <- as.data.frame(SEDEG_hh_aggr)

# remove rows where EC_FCH4_F_ANNOPTLM = NA

SEDEG_hh_aggr <- SEDEG_hh_aggr %>% drop_na(EC_FCH4_F_ANNOPTLM)

# move the new columns before EC_FCH4
SEDEG_hh_aggr <- SEDEG_hh_aggr %>% relocate(ch_FCH4_nmolCH4m2s1_mean, ch_FCH4_nmolCH4m2s1_sd, ch_FCH4_nmolCH4m2s1_median, ch_FCH4_nmolCH4m2s1_IQR, .before = EC_FCH4)

write.csv(SEDEG_hh_aggr, "path/SE_DEG_combined_ch_datetimeaggr.csv", row.names = FALSE)


# HOURLY AGGREGATION

# create new date_hour column for grouping
SEDEG_ECfilt$EC_TIMESTAMP_START <- as_datetime(SEDEG_ECfilt$EC_TIMESTAMP_START)

SEDEG_ECfilt <- SEDEG_ECfilt %>%
  mutate(date_hour = format(EC_TIMESTAMP_START, "%Y-%m-%d %H"))

# remove rows where ch_Datetime = NA

SEDEG_ECfilt <- SEDEG_ECfilt %>% drop_na(ch_Datetime)

# Calculate the u and v wind components
SEDEG_ECfilt$u.wind <- -SEDEG_ECfilt$EC_WS *sin(2*pi *SEDEG_ECfilt$EC_WD/360)
SEDEG_ECfilt$v.wind <- -SEDEG_ECfilt$EC_WS *cos(2*pi *SEDEG_ECfilt$EC_WD/360)
# # gap-filled (no gap-filled WD)
SEDEG_ECfilt$u.wind.F <- -SEDEG_ECfilt$EC_WS_F *sin(2*pi *SEDEG_ECfilt$EC_WD/360)
SEDEG_ECfilt$v.wind.F <- -SEDEG_ECfilt$EC_WS_F *cos(2*pi *SEDEG_ECfilt$EC_WD/360)

# set duplicated values to NA in EC columns
# create vector with EC column names
#rename EC_TIMESTAMP_START and END
colnames(SEDEG_ECfilt)[c(3,4)] <- c("TIMESTAMP_START", "TIMESTAMP_END")

cols <- SEDEG_ECfilt %>% dplyr::select(starts_with("EC_")) %>% colnames()

# convert duplicates to NA in EC columns
SEDEG_ECfilt[cols] <- sapply(SEDEG_ECfilt[cols], function(x) 
  ave(x, as_datetime(SEDEG_ECfilt$TIMESTAMP_START), FUN = function(x) replace(x, duplicated(x), NA)))

# convert character and logi to numeric
cols.num <- SEDEG_ECfilt %>% dplyr::select(starts_with(c("EC_", "ch_CH4_", "ch_TS_"))) %>% colnames()
SEDEG_ECfilt[cols.num] <- sapply(SEDEG_ECfilt[cols.num],as.numeric)
SEDEG_ECfilt$ch_P <- as.numeric(SEDEG_ECfilt$ch_P)

SEDEG_ECfilt <- SEDEG_ECfilt %>% drop_na(ch_FCH4_nmolCH4m2s1)

#aggregate

SEDEG_hr_aggr_ECfilt <- SEDEG_ECfilt %>% 
  group_by(date_hour) %>%
  summarise_at(c("ch_FCH4_nmolCH4m2s1", "ch_TS_2", "ch_TS_10",
                 "ch_WTL", "ch_CH4_PAR_amb", "ch_CH4_Ta_amb",
                 "ch_PPT", "ch_WS", "ch_P", "ch_VPD", "ch_RH",
                 "EC_GPP_NT", "EC_RECO_NT", "EC_GPP_DT", "EC_RECO_DT", "EC_SWC_1",
                 "EC_TS_1", "EC_TS_2",  "EC_TS_3",
                 "EC_TS_4", "EC_TS_5", "EC_TS_6",
                 "EC_FCH4", "EC_FCH4_F", "EC_FCH4_F_ANNOPTLM",
                 "EC_USTAR", "EC_NEE_F", "EC_H_F", "EC_LE_F", "EC_SW_IN_F", "EC_LW_IN_F", 
                 "EC_NETRAD_F", "EC_PPFD_IN_F", "EC_VPD_F", "EC_RH_F", 
                 "EC_PA_F", "EC_TA_F", "EC_WTD_F", "EC_P_F", "EC_WS_F", 
                 "EC_LE_F_ANNOPTLM", "EC_NEE_F_ANNOPTLM", "EC_FCH4_F_ANNOPTLM_UNC",
                 "u.wind", "v.wind", "u.wind.F", "v.wind.F"), 
               list(mean = mean,
                    median = median), 
               na.rm = TRUE)

SEDEG_hr_aggr_ECfilt <- as.data.frame(SEDEG_hr_aggr_ECfilt)

# convert to datetime format
SEDEG_hr_aggr_ECfilt$date_hour <- as_datetime(ymd_h(SEDEG_hr_aggr_ECfilt$date_hour))

# add the missing columns

SEDEG_hr_aggr_ECfilt$Year <- year(SEDEG_hr_aggr_ECfilt$date_hour)
SEDEG_hr_aggr_ECfilt$Month <- month(SEDEG_hr_aggr_ECfilt$date_hour)
SEDEG_hr_aggr_ECfilt$Day <- day(SEDEG_hr_aggr_ECfilt$date_hour)
SEDEG_hr_aggr_ECfilt$Hour <- hour(SEDEG_hr_aggr_ECfilt$date_hour)
SEDEG_hr_aggr_ECfilt$DOY <- yday(SEDEG_hr_aggr_ECfilt$date_hour)

SEDEG_hr_aggr_ECfilt$ch_method <- "auto"
SEDEG_hr_aggr_ECfilt$SITE <- "SE-DEG"

SEDEG_hr_aggr_ECfilt <- SEDEG_hr_aggr_ECfilt %>% relocate(SITE, .before = date_hour)
SEDEG_hr_aggr_ECfilt <- SEDEG_hr_aggr_ECfilt %>% relocate(ch_FCH4_nmolCH4m2s1_median, EC_FCH4_F_ANNOPTLM_mean, EC_FCH4_F_ANNOPTLM_median, .after = ch_FCH4_nmolCH4m2s1_mean)

# rename date_hour
colnames(SEDEG_hr_aggr_ECfilt)[2] <- "TIMESTAMP"

# add dom_veg
SEDEG_hr_aggr_ECfilt$DOM_VEG <- "MOSS_SPHAGNUM"

# add veg classes
SEDEG_hr_aggr_ECfilt$MOSS_BROWN <- 0
SEDEG_hr_aggr_ECfilt$MOSS_SPHAGNUM <- 1
SEDEG_hr_aggr_ECfilt$AERENCHYMATOUS <- 1
SEDEG_hr_aggr_ECfilt$ERI_SHRUB <- 1
SEDEG_hr_aggr_ECfilt$TREE <- 0

# add climate
SEDEG_hr_aggr_ECfilt$KOPPEN <- "Dfc"

# add ecosystem type
SEDEG_hr_aggr_ECfilt$SITE_CLASSIFICATION <- "fen"

# calculate wind direction average
SEDEG_hr_aggr_ECfilt$EC_WD_AVG <- (atan2(SEDEG_hr_aggr_ECfilt$u.wind_mean, SEDEG_hr_aggr_ECfilt$v.wind_mean) *360/2/pi) +180

# # gap-filled:
SEDEG_hr_aggr_ECfilt$EC_WD_AVG_F <- (atan2(SEDEG_hr_aggr_ECfilt$u.wind.F_mean, SEDEG_hr_aggr_ECfilt$v.wind.F_mean) *360/2/pi) +180

# save as .csv

write.csv(SEDEG_hr_aggr_ECfilt, "path/SE_DEG_combined_filtered_with_outliers_hraggr.csv", row.names = FALSE)



# DAILY AGGREGATION

SEDEG_d_aggr_ECfilt <- SEDEG_ECfilt %>% 
  group_by(date(SEDEG_ECfilt$TIMESTAMP_START)) %>%
  summarise_at(c("ch_FCH4_nmolCH4m2s1", "ch_TS_2", "ch_TS_10",
                 "ch_WTL", "ch_CH4_PAR_amb", "ch_CH4_Ta_amb",
                 "ch_PPT", "ch_WS", "ch_P", "ch_VPD", "ch_RH", "EC_SWC_1", "ECCH_diff",
                 "EC_TS_1", "EC_TS_2",  "EC_TS_3",
                 "EC_TS_4", "EC_TS_5", "EC_TS_6",
                 "EC_FCH4", "EC_FCH4_F", "EC_FCH4_F_ANNOPTLM", "EC_GPP_NT", "EC_RECO_NT", "EC_GPP_DT", "EC_RECO_DT",
                 "EC_USTAR", "EC_NEE_F", "EC_H_F", "EC_LE_F", "EC_SW_IN_F", "EC_LW_IN_F", 
                 "EC_NETRAD_F", "EC_PPFD_IN_F", "EC_VPD_F", "EC_RH_F", 
                 "EC_PA_F", "EC_TA_F", "EC_WTD_F", "EC_P_F", "EC_WS_F", 
                 "EC_LE_F_ANNOPTLM", "EC_NEE_F_ANNOPTLM", "EC_FCH4_F_ANNOPTLM_UNC",
                 "u.wind", "v.wind", "u.wind.F", "v.wind.F"), 
               list(mean = mean,
                    median = median), 
               na.rm = TRUE)

SEDEG_d_aggr_ECfilt <- as.data.frame(SEDEG_d_aggr_ECfilt)

# remove rows where CH4 ANN = NA

SEDEG_d_aggr_ECfilt <- SEDEG_d_aggr_ECfilt %>% drop_na(EC_FCH4_F_ANNOPTLM_mean)

# rename the date column

colnames(SEDEG_d_aggr_ECfilt)[1] <- "DATE"

# convert to date format
SEDEG_d_aggr_ECfilt$DATE <- as_date(SEDEG_d_aggr_ECfilt$DATE)

# calculate wind direction average
SEDEG_d_aggr_ECfilt$EC_WD_AVG <- (atan2(SEDEG_d_aggr_ECfilt$u.wind_mean, SEDEG_d_aggr_ECfilt$v.wind_mean) *360/2/pi) +180

# # gap-filled:
SEDEG_d_aggr_ECfilt$EC_WD_AVG_F <- (atan2(SEDEG_d_aggr_ECfilt$u.wind.F_mean, SEDEG_d_aggr_ECfilt$v.wind.F_mean) *360/2/pi) +180

# add the missing columns

SEDEG_d_aggr_ECfilt$Year <- year(SEDEG_d_aggr_ECfilt$DATE)
SEDEG_d_aggr_ECfilt$Month <- month(SEDEG_d_aggr_ECfilt$DATE)
SEDEG_d_aggr_ECfilt$Day <- day(SEDEG_d_aggr_ECfilt$DATE)
SEDEG_d_aggr_ECfilt$DOY <- yday(SEDEG_d_aggr_ECfilt$DATE)

SEDEG_d_aggr_ECfilt$ch_method <- "auto"
SEDEG_d_aggr_ECfilt$SITE <- "SE-DEG"

SEDEG_d_aggr_ECfilt <- SEDEG_d_aggr_ECfilt %>% relocate(SITE, .before = DATE)

# add dom_veg
SEDEG_d_aggr_ECfilt$DOM_VEG <- "MOSS_SPHAGNUM"

# add veg classes
SEDEG_d_aggr_ECfilt$MOSS_BROWN <- 0
SEDEG_d_aggr_ECfilt$MOSS_SPHAGNUM <- 1
SEDEG_d_aggr_ECfilt$AERENCHYMATOUS <- 1
SEDEG_d_aggr_ECfilt$ERI_SHRUB <- 1
SEDEG_d_aggr_ECfilt$TREE <- 0

# add climate
SEDEG_d_aggr_ECfilt$KOPPEN <- "Dfc"

# add ecosystem type
SEDEG_d_aggr_ECfilt$SITE_CLASSIFICATION <- "fen"

# save as .csv

write.csv(SEDEG_d_aggr_ECfilt, "path/SE_DEG_chec_daggr.csv", row.names = FALSE)

# WEEKLY AGGREGATION

SEDEG_week_aggr_ECfilt <- SEDEG_ECfilt %>% 
  group_by(year(SEDEG_ECfilt$TIMESTAMP_START), isoweek(SEDEG_ECfilt$TIMESTAMP_START)) %>%
  summarise_at(c("ch_FCH4_nmolCH4m2s1", "ch_TS_2", "ch_TS_10",
                 "ch_WTL", "ch_CH4_PAR_amb", "ch_CH4_Ta_amb",
                 "ch_PPT", "ch_WS", "ch_P", "ch_VPD", "ch_RH", "ECCH_diff" ,"EC_SWC_1",
                 "EC_TS_1", "EC_TS_2",  "EC_TS_3",
                 "EC_TS_4", "EC_TS_5", "EC_TS_6",
                 "EC_FCH4", "EC_FCH4_F", "EC_FCH4_F_ANNOPTLM", "EC_GPP_NT", "EC_RECO_NT", "EC_GPP_DT", "EC_RECO_DT",
                 "EC_USTAR", "EC_NEE_F", "EC_H_F", "EC_LE_F", "EC_SW_IN_F", "EC_LW_IN_F", 
                 "EC_NETRAD_F", "EC_PPFD_IN_F", "EC_VPD_F", "EC_RH_F", 
                 "EC_PA_F", "EC_TA_F", "EC_WTD_F", "EC_P_F", "EC_WS_F", 
                 "EC_LE_F_ANNOPTLM", "EC_NEE_F_ANNOPTLM", "EC_FCH4_F_ANNOPTLM_UNC",
                 "u.wind", "v.wind", "u.wind.F", "v.wind.F"), 
               list(mean = mean,
                    median = median), 
               na.rm = TRUE)

SEDEG_week_aggr_ECfilt <- as.data.frame(SEDEG_week_aggr_ECfilt)

# rename the date column

colnames(SEDEG_week_aggr_ECfilt)[1] <- "Year"
colnames(SEDEG_week_aggr_ECfilt)[2] <- "Week_of_year"

# add the missing extra columns

SEDEG_week_aggr_ECfilt$ch_method <- "auto"
SEDEG_week_aggr_ECfilt$SITE <- "SE-DEG"

SEDEG_week_aggr_ECfilt <- SEDEG_week_aggr_ECfilt %>% relocate(SITE, .before = Year)

# add dom_veg
SEDEG_week_aggr_ECfilt$DOM_VEG <- "MOSS_SPHAGNUM"

# add veg classes
SEDEG_week_aggr_ECfilt$MOSS_BROWN <- 0
SEDEG_week_aggr_ECfilt$MOSS_SPHAGNUM <- 1
SEDEG_week_aggr_ECfilt$AERENCHYMATOUS <- 1
SEDEG_week_aggr_ECfilt$ERI_SHRUB <- 1
SEDEG_week_aggr_ECfilt$TREE <- 0

# add climate
SEDEG_week_aggr_ECfilt$KOPPEN <- "Dfc"

# add ecosystem type
SEDEG_week_aggr_ECfilt$SITE_CLASSIFICATION <- "fen"

SEDEG_week_aggr_ECfilt$EC_SENSOR <- "closed-path"

# calculate wind direction average
SEDEG_week_aggr_ECfilt$EC_WD_AVG <- (atan2(SEDEG_week_aggr_ECfilt$u.wind_mean, SEDEG_week_aggr_ECfilt$v.wind_mean) *360/2/pi) +180

# # gap-filled:
SEDEG_week_aggr_ECfilt$EC_WD_AVG_F <- (atan2(SEDEG_week_aggr_ECfilt$u.wind.F_mean, SEDEG_week_aggr_ECfilt$v.wind.F_mean) *360/2/pi) +180

# Extract the correct month directly using lubridate::make_date and ISO week (1st day of the week is Monday)
SEDEG_week_aggr_ECfilt$Month <- month(make_date(SEDEG_week_aggr_ECfilt$Year) + weeks(SEDEG_week_aggr_ECfilt$Week_of_year))

SEDEG_week_aggr_ECfilt <- SEDEG_week_aggr_ECfilt %>% relocate(Month, .after = Week_of_year)

# multiple weeks and months are wrong
SEDEG_week_aggr_ECfilt$Month[SEDEG_week_aggr_ECfilt$Week_of_year == 22 & SEDEG_week_aggr_ECfilt$Year == 2015] <- 5
SEDEG_week_aggr_ECfilt$Month[SEDEG_week_aggr_ECfilt$Week_of_year == 26 & SEDEG_week_aggr_ECfilt$Year == 2015] <- 6
SEDEG_week_aggr_ECfilt$Month[SEDEG_week_aggr_ECfilt$Week_of_year == 31 & SEDEG_week_aggr_ECfilt$Year == 2015] <- 7
SEDEG_week_aggr_ECfilt$Month[SEDEG_week_aggr_ECfilt$Week_of_year == 35 & SEDEG_week_aggr_ECfilt$Year == 2015] <- 8
SEDEG_week_aggr_ECfilt$Month[SEDEG_week_aggr_ECfilt$Week_of_year == 39 & SEDEG_week_aggr_ECfilt$Year == 2015] <- 9
SEDEG_week_aggr_ECfilt$Month[SEDEG_week_aggr_ECfilt$Week_of_year == 44 & SEDEG_week_aggr_ECfilt$Year == 2015] <- 10
SEDEG_week_aggr_ECfilt$Month[SEDEG_week_aggr_ECfilt$Week_of_year == 48 & SEDEG_week_aggr_ECfilt$Year == 2015] <- 11

SEDEG_week_aggr_ECfilt$Month[SEDEG_week_aggr_ECfilt$Week_of_year == 26 & SEDEG_week_aggr_ECfilt$Year == 2016] <- 6
SEDEG_week_aggr_ECfilt$Month[SEDEG_week_aggr_ECfilt$Week_of_year == 44 & SEDEG_week_aggr_ECfilt$Year == 2016] <- 10
# now the weeks and months are correct

# save as .csv

write.csv(SEDEG_week_aggr_ECfilt, "path/SE_DEG_chec_weekaggr.csv", row.names = FALSE)

# MONTHLY AGGREGATION

SEDEG_month_aggr_ECfilt <- SEDEG_ECfilt %>% 
  group_by(year(SEDEG_ECfilt$TIMESTAMP_START), month(SEDEG_ECfilt$TIMESTAMP_START)) %>%
  summarise_at(c("ch_FCH4_nmolCH4m2s1", "ch_TS_2", "ch_TS_10",
                 "ch_WTL", "ch_CH4_PAR_amb", "ch_CH4_Ta_amb",
                 "ch_PPT", "ch_WS", "ch_P", "ch_VPD", "ch_RH", "ECCH_diff" ,"EC_SWC_1",
                 "EC_TS_1", "EC_TS_2",  "EC_TS_3",
                 "EC_TS_4", "EC_TS_5", "EC_TS_6",
                 "EC_FCH4", "EC_FCH4_F", "EC_FCH4_F_ANNOPTLM", "EC_GPP_NT", "EC_RECO_NT", "EC_GPP_DT", "EC_RECO_DT",
                 "EC_USTAR", "EC_NEE_F", "EC_H_F", "EC_LE_F", "EC_SW_IN_F", "EC_LW_IN_F", 
                 "EC_NETRAD_F", "EC_PPFD_IN_F", "EC_VPD_F", "EC_RH_F", 
                 "EC_PA_F", "EC_TA_F", "EC_WTD_F", "EC_P_F", "EC_WS_F", 
                 "EC_LE_F_ANNOPTLM", "EC_NEE_F_ANNOPTLM", "EC_FCH4_F_ANNOPTLM_UNC",
                 "u.wind", "v.wind", "u.wind.F", "v.wind.F"), 
               list(mean = mean,
                    median = median), 
               na.rm = TRUE)

SEDEG_month_aggr_ECfilt <- as.data.frame(SEDEG_month_aggr_ECfilt)

# rename the date column

colnames(SEDEG_month_aggr_ECfilt)[1] <- "Year"
colnames(SEDEG_month_aggr_ECfilt)[2] <- "Month"

# add the missing extra columns

SEDEG_month_aggr_ECfilt$ch_method <- "auto"
SEDEG_month_aggr_ECfilt$SITE <- "SE-DEG"

SEDEG_month_aggr_ECfilt <- SEDEG_month_aggr_ECfilt %>% relocate(SITE, .before = Year)

# add dom_veg
SEDEG_month_aggr_ECfilt$DOM_VEG <- "MOSS_SPHAGNUM"

# add veg classes
SEDEG_month_aggr_ECfilt$MOSS_BROWN <- 0
SEDEG_month_aggr_ECfilt$MOSS_SPHAGNUM <- 1
SEDEG_month_aggr_ECfilt$AERENCHYMATOUS <- 1
SEDEG_month_aggr_ECfilt$ERI_SHRUB <- 1
SEDEG_month_aggr_ECfilt$TREE <- 0

# add climate
SEDEG_month_aggr_ECfilt$KOPPEN <- "Dfc"

# add ecosystem type
SEDEG_month_aggr_ECfilt$SITE_CLASSIFICATION <- "fen"

# calculate wind direction average
SEDEG_month_aggr_ECfilt$EC_WD_AVG <- (atan2(SEDEG_month_aggr_ECfilt$u.wind_mean, SEDEG_month_aggr_ECfilt$v.wind_mean) *360/2/pi) +180

# # gap-filled:
SEDEG_month_aggr_ECfilt$EC_WD_AVG_F <- (atan2(SEDEG_month_aggr_ECfilt$u.wind.F_mean, SEDEG_month_aggr_ECfilt$v.wind.F_mean) *360/2/pi) +180

# save as .csv

write.csv(SEDEG_month_aggr_ECfilt, "path/SE_DEG_chec_monthaggr_ECfilt.csv", row.names = FALSE)

# ANNUAL AGGREGATION

SEDEG_yr_aggr_ECfilt <- SEDEG_ECfilt %>% 
  group_by(year(SEDEG_ECfilt$TIMESTAMP_START)) %>%
  summarise_at(c("ch_FCH4_nmolCH4m2s1", "ch_TS_2", "ch_TS_10",
                 "ch_WTL", "ch_CH4_PAR_amb", "ch_CH4_Ta_amb",
                 "ch_PPT", "ch_WS", "ch_P", "ch_VPD", "ch_RH",
                 "EC_SWC_1",
                 "EC_TS_1", "EC_TS_2",  "EC_TS_3",
                 "EC_TS_4", "EC_TS_5", "EC_TS_6",
                 "EC_FCH4", "EC_FCH4_F", "EC_FCH4_F_ANNOPTLM", "EC_GPP_NT", "EC_RECO_NT", "EC_GPP_DT", "EC_RECO_DT",
                 "EC_USTAR", "EC_NEE_F", "EC_H_F", "EC_LE_F", "EC_SW_IN_F", "EC_LW_IN_F", 
                 "EC_NETRAD_F", "EC_PPFD_IN_F", "EC_VPD_F", "EC_RH_F", 
                 "EC_PA_F", "EC_TA_F", "EC_WTD_F", "EC_P_F", "EC_WS_F", 
                 "EC_LE_F_ANNOPTLM", "EC_NEE_F_ANNOPTLM", "EC_FCH4_F_ANNOPTLM_UNC",
                 "u.wind", "v.wind", "u.wind.F", "v.wind.F"), 
               list(mean = mean,
                    median = median), 
               na.rm = TRUE)

SEDEG_yr_aggr_ECfilt <- as.data.frame(SEDEG_yr_aggr_ECfilt)

# rename the date column

colnames(SEDEG_yr_aggr_ECfilt)[1] <- "Year"

# add the missing extra columns

SEDEG_yr_aggr_ECfilt$ch_method <- "auto"
SEDEG_yr_aggr_ECfilt$SITE <- "SE-DEG"

SEDEG_yr_aggr_ECfilt <- SEDEG_yr_aggr_ECfilt %>% relocate(SITE, .before = Year)

# add dom_veg
SEDEG_yr_aggr_ECfilt$DOM_VEG <- "MOSS_SPHAGNUM"

# add veg classes
SEDEG_yr_aggr_ECfilt$MOSS_BROWN <- 0
SEDEG_yr_aggr_ECfilt$MOSS_SPHAGNUM <- 1
SEDEG_yr_aggr_ECfilt$AERENCHYMATOUS <- 1
SEDEG_yr_aggr_ECfilt$ERI_SHRUB <- 1
SEDEG_yr_aggr_ECfilt$TREE <- 0

# add climate
SEDEG_yr_aggr_ECfilt$KOPPEN <- "Dfc"

# add ecosystem type
SEDEG_yr_aggr_ECfilt$SITE_CLASSIFICATION <- "fen"

SEDEG_yr_aggr_ECfilt$EC_SENSOR <- "closed-path"

# calculate wind direction average
SEDEG_yr_aggr_ECfilt$EC_WD_AVG <- (atan2(SEDEG_yr_aggr_ECfilt$u.wind_mean, SEDEG_yr_aggr_ECfilt$v.wind_mean) *360/2/pi) +180

# # gap-filled:
SEDEG_yr_aggr_ECfilt$EC_WD_AVG_F <- (atan2(SEDEG_yr_aggr_ECfilt$u.wind.F_mean, SEDEG_yr_aggr_ECfilt$v.wind.F_mean) *360/2/pi) +180

# save as .csv

write.csv(SEDEG_yr_aggr_ECfilt, "path/SE_DEG_chec_yraggr_ECfilt.csv", row.names = FALSE)


### US-HO1 ###

# read csv

# EC
USHO1_EC_d <- read.csv("path/US_HO1_EC_subset.csv")

# convert to datetime
USHO1_EC_d$EC_TIMESTAMP_START <- as_datetime(USHO1_EC_d$EC_TIMESTAMP_START)
USHO1_EC_d$EC_TIMESTAMP_END <- as_datetime(USHO1_EC_d$EC_TIMESTAMP_END)

# CHAMBER
USHO1_ch_d <- read.csv("path/US_HO1.csv")

USHO1_ch_d$ch_Datetime <- as_datetime(USHO1_ch_d$ch_Datetime)

# combine using merge
USHO1_all <- merge(x = USHO1_ch_d, y = USHO1_EC_d, by.x = c("ch_Datetime","SITE"), by.y = c("EC_TIMESTAMP_START","SITE"),
                   all.y = TRUE)

# move EC_TIMESTAMP_END and SITE
USHO1_all <- USHO1_all %>% relocate(EC_TIMESTAMP_END, .before = SITE)
USHO1_all <- USHO1_all %>% relocate(SITE, .before = ch_Datetime)

# rename ch_Datetime
colnames(USHO1_all)[2] <- "ch_EC_TIMESTAMP_START"
colnames(USHO1_all)[3] <- "ch_EC_TIMESTAMP_END"

USHO1_all$ch_EC_TIMESTAMP_START <- as_datetime(USHO1_all$ch_EC_TIMESTAMP_START, tz="UTC")
USHO1_all$EC_TIMESTAMP_END <- as_datetime(USHO1_all$ch_EC_TIMESTAMP_END, tz="UTC")

# subset to start from 16.4.2012 13.30.00 because before that chamber measurements NA 
# and end 7.11.2016 23.30.00 because after that chamber measurements NA

Sys.setenv(TZ = "UTC")

USHO1_all_sub <- subset(USHO1_all, ch_EC_TIMESTAMP_START >= "2012-04-16 13:30:00" & ch_EC_TIMESTAMP_END <= "2016-11-08 00:00:00")

# remove the extra column 
USHO1_all_sub <- select(USHO1_all_sub, -58)

# CONVERT CHAMBER TO nmol CH4 m-2 s-1
# IN UMOL CH4 M-2 S-1

# add new column for nmol chamber unit
USHO1_all_sub$ch_FCH4_nmolCH4m2s1 <- USHO1_all_sub$ch_FCH4_umolCH4m2s1*1000

# move next to umol
USHO1_all_sub <- USHO1_all_sub %>% relocate(ch_FCH4_nmolCH4m2s1, .after = ch_FCH4_umolCH4m2s1)

write.csv(USHO1_all_sub, "path/US_HO1_combined.csv", row.names = FALSE)

# read in the file again
USHO1 <- read.csv("path/US_HO1_combined.csv")

# make start from 2012-04-17 because 16th is not a full day
Sys.setenv(TZ = "UTC")
USHO1$ch_EC_TIMESTAMP_START <- as_datetime(USHO1$ch_EC_TIMESTAMP_START)
USHO1 <- subset(USHO1, ch_EC_TIMESTAMP_START >= "2012-04-17 00:00:00")

# SET EC GAP-FILLED DATA AND CHAMBER DATA TO NA WHEN TIME GAP LONGER THAN 2 MONTHS
USHO1_new <- USHO1 %>% mutate(across(c(ch_FCH4_umolCH4m2s1, ch_FCH4_nmolCH4m2s1, ch_SM_X, ch_TS_X, 
                                       EC_FCH4_F, EC_NEE_F, EC_H_F, EC_LE_F, EC_SW_IN_F, EC_SW_OUT_F, EC_LW_IN_F, EC_LW_OUT_F,
                                       EC_NETRAD_F, EC_PPFD_IN_F,
                                       EC_VPD_F, EC_RH_F, EC_PA_F, EC_TA_F, EC_P_F,EC_WTD_F,
                                       EC_WS_F, EC_LE_F_ANNOPTLM, EC_NEE_F_ANNOPTLM, EC_FCH4_F_ANNOPTLM, 
                                       EC_FCH4_F_RANDUNC, EC_FCH4_F_ANNOPTLM_UNC), ~ ifelse(EC_FCH4_F_ANNOPTLM_QC == 3, NA, .)))

# convert datetimes to datetime format
USHO1_new$ch_EC_TIMESTAMP_START <- as_datetime(USHO1_new$ch_EC_TIMESTAMP_START)
USHO1_new$ch_EC_TIMESTAMP_END <- as_datetime(USHO1_new$ch_EC_TIMESTAMP_END)

# move EC_FCH4_F_ANNOPTLM next to other flux columns
USHO1_new <- USHO1_new %>% relocate(c(EC_FCH4, EC_FCH4_F, EC_FCH4_F_ANNOPTLM), .after = ch_FCH4_nmolCH4m2s1)

USHO1_with_outliers_ECall <- USHO1_new

USHO1_ECfilt <- USHO1_with_outliers_ECall %>% mutate(across(c(ch_FCH4_nmolCH4m2s1, EC_FCH4, EC_FCH4_F, EC_FCH4_F_ANNOPTLM), 
                                                                          ~ ifelse(ch_FCH4_umolCH4m2s1 == "NA", NA, .)))

write.csv(USHO1_ECfilt, "path/US_HO1_combined_ECfilt.csv", row.names = FALSE)

# HALF-HOURLY AGGREGATION

# note: ch_FCH4_umolCH4m2s1 is actually in nmol CH4 m-2 s-1, and that is why that is being used here and in the rest of the US-HO1 code
USHO1_ECfilt <- USHO1_ECfilt %>% drop_na(ch_FCH4_umolCH4m2s1)

# ch_TS_X appears to be character, change to numeric
USHO1_ECfilt$ch_TS_X <- as.numeric(USHO1_ECfilt$ch_TS_X)

USHO1_hh_aggr <- USHO1_ECfilt %>% 
  group_by(across(-c(ch_FCH4_umolCH4m2s1, ch_FCH4_nmolCH4m2s1, ch_ID,
                     ch_SM_X, ch_TS_X, ch_TS_2, ch_TS_5, ch_TS_10, ch_TS_15, ch_TS_30, ch_WTL,
                     ch_LC_class, ch_OutlierFlag, ECCH_diff))) %>%
  summarise_at(c("ch_FCH4_umolCH4m2s1", "ch_SM_X", "ch_TS_X", "ECCH_diff"), 
               list(mean = mean,
                    median = median), 
               na.rm = TRUE)

USHO1_hh_aggr <- as.data.frame(USHO1_hh_aggr)

# remove rows where EC_FCH4_F_ANNOPTLM = NA
USHO1_hh_aggr <- USHO1_hh_aggr %>% drop_na(EC_FCH4_F_ANNOPTLM)

# move the new columns before EC_FCH4
USHO1_hh_aggr <- USHO1_hh_aggr %>% relocate(ch_FCH4_umolCH4m2s1_mean, ch_FCH4_umolCH4m2s1_sd, ch_FCH4_umolCH4m2s1_median, ch_FCH4_umolCH4m2s1_IQR, .before = EC_FCH4)

write.csv(USHO1_hh_aggr, "path/US_HO1_combined_ch_datetimeaggr.csv", row.names = FALSE)



# HOURLY AGGREGATION

# create new date_hour column for grouping
USHO1_ECfilt$ch_EC_TIMESTAMP_START <- as_datetime(USHO1_ECfilt$ch_EC_TIMESTAMP_START)
USHO1_ECfilt$ch_EC_TIMESTAMP_END <- as_datetime(USHO1_ECfilt$ch_EC_TIMESTAMP_END)

USHO1_ECfilt <- USHO1_ECfilt %>%
  mutate(date_hour = format(ch_EC_TIMESTAMP_START, "%Y-%m-%d %H"))


# set duplicated values to NA in EC columns
# create vector with EC column names
#rename EC_TIMESTAMP_START and END
colnames(USHO1_ECfilt)[c(2,3)] <- c("TIMESTAMP_START", "TIMESTAMP_END")

cols <- USHO1_ECfilt %>% dplyr::select(starts_with("EC_")) %>% colnames()

USHO1_ECfilt$ch_WTL <- as.numeric(USHO1_ECfilt$ch_WTL)

USHO1_ECfilt <- USHO1_ECfilt %>% drop_na(ch_FCH4_umolCH4m2s1)

# convert duplicates to NA in EC columns
USHO1_ECfilt[cols] <- sapply(USHO1_ECfilt[cols], function(x) 
  ave(x, as_datetime(USHO1_ECfilt$TIMESTAMP_START), FUN = function(x) replace(x, duplicated(x), NA)))

# convert character and logi to numeric
cols.num <- USHO1_ECfilt %>% dplyr::select(starts_with(c("EC_", "ch_TS_"))) %>% dplyr::select(-EC_FCH4_MEASTYPE) %>% colnames()
USHO1_ECfilt[cols.num] <- sapply(USHO1_ECfilt[cols.num],as.numeric)

# Calculate the u and v wind components
USHO1_ECfilt$u.wind <- -USHO1_ECfilt$EC_WS *sin(2*pi *USHO1_ECfilt$EC_WD/360)
USHO1_ECfilt$v.wind <- -USHO1_ECfilt$EC_WS *cos(2*pi *USHO1_ECfilt$EC_WD/360)
# gap-filled (no gap-filled WD)
USHO1_ECfilt$u.wind.F <- -USHO1_ECfilt$EC_WS_F *sin(2*pi *USHO1_ECfilt$EC_WD/360)
USHO1_ECfilt$v.wind.F <- -USHO1_ECfilt$EC_WS_F *cos(2*pi *USHO1_ECfilt$EC_WD/360)

# aggregate
USHO1_hr_aggr_ECfilt <- USHO1_ECfilt %>% 
  group_by(date_hour) %>%
  summarise_at(c("ch_FCH4_umolCH4m2s1", "ch_SM_X", "ch_TS_X",
                 "EC_FCH4", "EC_FCH4_F", "EC_FCH4_F_ANNOPTLM",
                 "EC_USTAR", "EC_GPP_NT", "EC_RECO_NT", "EC_GPP_DT", "EC_RECO_DT", "EC_NEE_F", "EC_H_F", "EC_LE_F", "EC_SW_IN_F", "EC_SW_OUT_F", "EC_LW_IN_F", 
                 "EC_LW_OUT_F", "EC_NETRAD_F", "EC_PPFD_IN_F", "EC_VPD_F", "EC_RH_F", 
                 "EC_PA_F", "EC_TA_F", "EC_WTD_F", "EC_P_F", "EC_WS_F", 
                 "EC_LE_F_ANNOPTLM", "EC_NEE_F_ANNOPTLM", "EC_FCH4_F_ANNOPTLM_UNC",
                 "u.wind", "v.wind", "u.wind.F", "v.wind.F"), 
               list(mean = mean,
                    median = median), 
               na.rm = TRUE)

USHO1_hr_aggr_ECfilt <- as.data.frame(USHO1_hr_aggr_ECfilt)

# convert to datetime format
USHO1_hr_aggr_ECfilt$date_hour <- as_datetime(ymd_h(USHO1_hr_aggr_ECfilt$date_hour))

# add the missing columns

USHO1_hr_aggr_ECfilt$Year <- year(USHO1_hr_aggr_ECfilt$date_hour)
USHO1_hr_aggr_ECfilt$Month <- month(USHO1_hr_aggr_ECfilt$date_hour)
USHO1_hr_aggr_ECfilt$Day <- day(USHO1_hr_aggr_ECfilt$date_hour)
USHO1_hr_aggr_ECfilt$Hour <- hour(USHO1_hr_aggr_ECfilt$date_hour)
USHO1_hr_aggr_ECfilt$DOY <- yday(USHO1_hr_aggr_ECfilt$date_hour)

USHO1_hr_aggr_ECfilt$ch_method <- "auto"
USHO1_hr_aggr_ECfilt$SITE <- "US-HO1"

USHO1_hr_aggr_ECfilt <- USHO1_hr_aggr_ECfilt %>% relocate(SITE, .before = date_hour)

# rename date_hour
colnames(USHO1_hr_aggr_ECfilt)[2] <- "TIMESTAMP"

USHO1_hr_aggr_ECfilt <- USHO1_hr_aggr_ECfilt %>% drop_na(EC_FCH4_F_ANNOPTLM_mean)

# add dom_veg
USHO1_hr_aggr_ECfilt$DOM_VEG <- "TREE"

# add veg classes
USHO1_hr_aggr_ECfilt$MOSS_BROWN <- 0
USHO1_hr_aggr_ECfilt$MOSS_SPHAGNUM <- 0
USHO1_hr_aggr_ECfilt$AERENCHYMATOUS <- 0
USHO1_hr_aggr_ECfilt$ERI_SHRUB <- 0
USHO1_hr_aggr_ECfilt$TREE <- 1

# add climate
USHO1_hr_aggr_ECfilt$KOPPEN <- "Dfb"

# add ecosystem type
USHO1_hr_aggr_ECfilt$SITE_CLASSIFICATION <- "upland"

# calculate wind direction average
USHO1_hr_aggr_ECfilt$EC_WD_AVG <- (atan2(USHO1_hr_aggr_ECfilt$u.wind_mean, USHO1_hr_aggr_ECfilt$v.wind_mean) *360/2/pi) +180

# # gap-filled:
USHO1_hr_aggr_ECfilt$EC_WD_AVG_F <- (atan2(USHO1_hr_aggr_ECfilt$u.wind.F_mean, USHO1_hr_aggr_ECfilt$v.wind.F_mean) *360/2/pi) +180

# save as .csv

write.csv(USHO1_hr_aggr_ECfilt, "path/US_HO1_combined_filtered_with_outliers_ch_hraggr.csv", row.names = FALSE)


# DAILY AGGREGATION

USHO1_d_aggr_ECfilt <- USHO1_ECfilt %>% 
  group_by(date(USHO1_ECfilt$TIMESTAMP_START)) %>%
  summarise_at(c("ch_FCH4_umolCH4m2s1", "ch_SM_X", "ch_TS_X",
                 "EC_FCH4", "EC_FCH4_F", "EC_FCH4_F_ANNOPTLM",
                 "EC_USTAR", "EC_GPP_NT", "EC_RECO_NT", "EC_GPP_DT", "EC_RECO_DT", "EC_NEE_F", "EC_H_F", "EC_LE_F", "EC_SW_IN_F", "EC_SW_OUT_F", "EC_LW_IN_F", 
                 "EC_LW_OUT_F", "EC_NETRAD_F", "EC_PPFD_IN_F", "EC_VPD_F", "EC_RH_F", 
                 "EC_PA_F", "EC_TA_F", "EC_WTD_F", "EC_P_F", "EC_WS_F", 
                 "EC_LE_F_ANNOPTLM", "EC_NEE_F_ANNOPTLM", "EC_FCH4_F_ANNOPTLM_UNC",
                 "u.wind", "v.wind", "u.wind.F", "v.wind.F"), 
               list(mean = mean,
                    median = median), 
               na.rm = TRUE)

USHO1_d_aggr_ECfilt <- as.data.frame(USHO1_d_aggr_ECfilt)

# rename the date column

colnames(USHO1_d_aggr_ECfilt)[1] <- "DATE"

# convert to date format
USHO1_d_aggr_ECfilt$DATE <- as_date(USHO1_d_aggr_ECfilt$DATE)

# calculate wind direction average
USHO1_d_aggr_ECfilt$EC_WD_AVG <- (atan2(USHO1_d_aggr_ECfilt$u.wind_mean, USHO1_d_aggr_ECfilt$v.wind_mean) *360/2/pi) +180

# gap-filled:
USHO1_d_aggr_ECfilt$EC_WD_AVG_F <- (atan2(USHO1_d_aggr_ECfilt$u.wind.F_mean, USHO1_d_aggr_ECfilt$v.wind.F_mean) *360/2/pi) +180

# add the missing columns

USHO1_d_aggr_ECfilt$Year <- year(USHO1_d_aggr_ECfilt$DATE)
USHO1_d_aggr_ECfilt$Month <- month(USHO1_d_aggr_ECfilt$DATE)
USHO1_d_aggr_ECfilt$Day <- day(USHO1_d_aggr_ECfilt$DATE)
USHO1_d_aggr_ECfilt$DOY <- yday(USHO1_d_aggr_ECfilt$DATE)

USHO1_d_aggr_ECfilt$ch_method <- "auto"
USHO1_d_aggr_ECfilt$SITE <- "US-HO1"

USHO1_d_aggr_ECfilt <- USHO1_d_aggr_ECfilt %>% relocate(SITE, .before = DATE)

# remove rows with FCH4_ANNOTPLTM = NA

USHO1_d_aggr_ECfilt <- USHO1_d_aggr_ECfilt %>% drop_na(EC_FCH4_F_ANNOPTLM_mean) 

# add dom_veg
USHO1_d_aggr_ECfilt$DOM_VEG <- "TREE"

# add veg classes
USHO1_d_aggr_ECfilt$MOSS_BROWN <- 0
USHO1_d_aggr_ECfilt$MOSS_SPHAGNUM <- 0
USHO1_d_aggr_ECfilt$AERENCHYMATOUS <- 0
USHO1_d_aggr_ECfilt$ERI_SHRUB <- 0
USHO1_d_aggr_ECfilt$TREE <- 1

# add climate
USHO1_d_aggr_ECfilt$KOPPEN <- "Dfb"

# add ecosystem type
USHO1_d_aggr_ECfilt$SITE_CLASSIFICATION <- "upland"

# save as .csv
write.csv(USHO1_d_aggr_ECfilt, "path/US_HO1_chec_daggr_ECfilt.csv", row.names = FALSE)


# WEEKLY AGGREGATION

USHO1_week_aggr_ECfilt <- USHO1_ECfilt %>% 
  group_by(year(USHO1_ECfilt$TIMESTAMP_START), isoweek(USHO1_ECfilt$TIMESTAMP_START)) %>%
  summarise_at(c("ch_FCH4_umolCH4m2s1", "ch_SM_X", "ch_TS_X",
                 "EC_FCH4", "EC_FCH4_F", "EC_FCH4_F_ANNOPTLM",
                 "EC_USTAR", "EC_GPP_NT", "EC_RECO_NT", "EC_GPP_DT", "EC_RECO_DT", "EC_NEE_F", "EC_H_F", "EC_LE_F", "EC_SW_IN_F", "EC_SW_OUT_F", "EC_LW_IN_F", 
                 "EC_LW_OUT_F", "EC_NETRAD_F", "EC_PPFD_IN_F", "EC_VPD_F", "EC_RH_F", 
                 "EC_PA_F", "EC_TA_F", "EC_WTD_F", "EC_P_F", "EC_WS_F", 
                 "EC_LE_F_ANNOPTLM", "EC_NEE_F_ANNOPTLM", "EC_FCH4_F_ANNOPTLM_UNC",
                 "u.wind", "v.wind", "u.wind.F", "v.wind.F"), 
               list(mean = mean,
                    median = median), 
               na.rm = TRUE)

USHO1_week_aggr_ECfilt <- as.data.frame(USHO1_week_aggr_ECfilt)

# rename the date column

colnames(USHO1_week_aggr_ECfilt)[1] <- "Year"
colnames(USHO1_week_aggr_ECfilt)[2] <- "Week_of_year"

# add the missing extra columns

USHO1_week_aggr_ECfilt$ch_method <- "auto"
USHO1_week_aggr_ECfilt$SITE <- "US-HO1"

USHO1_week_aggr_ECfilt <- USHO1_week_aggr_ECfilt %>% relocate(SITE, .before = Year)

# add dom_veg
USHO1_week_aggr_ECfilt$DOM_VEG <- "TREE"

# add veg classes
USHO1_week_aggr_ECfilt$MOSS_BROWN <- 0
USHO1_week_aggr_ECfilt$MOSS_SPHAGNUM <- 0
USHO1_week_aggr_ECfilt$AERENCHYMATOUS <- 0
USHO1_week_aggr_ECfilt$ERI_SHRUB <- 0
USHO1_week_aggr_ECfilt$TREE <- 1

# add climate
USHO1_week_aggr_ECfilt$KOPPEN <- "Dfb"

# add ecosystem type
USHO1_week_aggr_ECfilt$SITE_CLASSIFICATION <- "upland"

# calculate wind direction average
USHO1_week_aggr_ECfilt$EC_WD_AVG <- (atan2(USHO1_week_aggr_ECfilt$u.wind_mean, USHO1_week_aggr_ECfilt$v.wind_mean) *360/2/pi) +180

# # gap-filled:
USHO1_week_aggr_ECfilt$EC_WD_AVG_F <- (atan2(USHO1_week_aggr_ECfilt$u.wind.F_mean, USHO1_week_aggr_ECfilt$v.wind.F_mean) *360/2/pi) +180

# add month

# Extract the correct month directly using lubridate::make_date and ISO week (1st day of the week is Monday)
USHO1_week_aggr_ECfilt$Month <- month(make_date(USHO1_week_aggr_ECfilt$Year) + weeks(USHO1_week_aggr_ECfilt$Week_of_year))
USHO1_week_aggr_ECfilt <- USHO1_week_aggr_ECfilt %>% relocate(Month, .after = Week_of_year)

# multiple weeks and months are wrong
USHO1_week_aggr_ECfilt$Month[USHO1_week_aggr_ECfilt$Week_of_year == 26 & USHO1_week_aggr_ECfilt$Year == 2012] <- 7
USHO1_week_aggr_ECfilt$Month[USHO1_week_aggr_ECfilt$Week_of_year == 35 & USHO1_week_aggr_ECfilt$Year == 2012] <- 9

USHO1_week_aggr_ECfilt$Month[USHO1_week_aggr_ECfilt$Week_of_year == 26 & USHO1_week_aggr_ECfilt$Year == 2013] <- 6
USHO1_week_aggr_ECfilt$Month[USHO1_week_aggr_ECfilt$Week_of_year == 35 & USHO1_week_aggr_ECfilt$Year == 2013] <- 8
USHO1_week_aggr_ECfilt$Month[USHO1_week_aggr_ECfilt$Week_of_year == 39 & USHO1_week_aggr_ECfilt$Year == 2013] <- 9
USHO1_week_aggr_ECfilt$Month[USHO1_week_aggr_ECfilt$Week_of_year == 44 & USHO1_week_aggr_ECfilt$Year == 2013] <- 10

USHO1_week_aggr_ECfilt$Month[USHO1_week_aggr_ECfilt$Week_of_year == 26 & USHO1_week_aggr_ECfilt$Year == 2014] <- 6
USHO1_week_aggr_ECfilt$Month[USHO1_week_aggr_ECfilt$Week_of_year == 31 & USHO1_week_aggr_ECfilt$Year == 2014] <- 7
USHO1_week_aggr_ECfilt$Month[USHO1_week_aggr_ECfilt$Week_of_year == 35 & USHO1_week_aggr_ECfilt$Year == 2014] <- 8
USHO1_week_aggr_ECfilt$Month[USHO1_week_aggr_ECfilt$Week_of_year == 39 & USHO1_week_aggr_ECfilt$Year == 2014] <- 9
USHO1_week_aggr_ECfilt$Month[USHO1_week_aggr_ECfilt$Week_of_year == 44 & USHO1_week_aggr_ECfilt$Year == 2014] <- 10

USHO1_week_aggr_ECfilt$Month[USHO1_week_aggr_ECfilt$Week_of_year == 26 & USHO1_week_aggr_ECfilt$Year == 2015] <- 6
USHO1_week_aggr_ECfilt$Month[USHO1_week_aggr_ECfilt$Week_of_year == 31 & USHO1_week_aggr_ECfilt$Year == 2015] <- 7
USHO1_week_aggr_ECfilt$Month[USHO1_week_aggr_ECfilt$Week_of_year == 35 & USHO1_week_aggr_ECfilt$Year == 2015] <- 8
USHO1_week_aggr_ECfilt$Month[USHO1_week_aggr_ECfilt$Week_of_year == 39 & USHO1_week_aggr_ECfilt$Year == 2016] <- 9

USHO1_week_aggr_ECfilt$Month[USHO1_week_aggr_ECfilt$Week_of_year == 26 & USHO1_week_aggr_ECfilt$Year == 2016] <- 6

# save as .csv

write.csv(USHO1_week_aggr_ECfilt, "path/US_HO1_chec_weekaggr_ECfilt.csv", row.names = FALSE)

# MONTHLY AGGREGATION

USHO1_month_aggr_ECfilt <- USHO1_ECfilt %>% 
  group_by(year(USHO1_ECfilt$TIMESTAMP_START), month(USHO1_ECfilt$TIMESTAMP_START)) %>%
  summarise_at(c("ch_FCH4_umolCH4m2s1", "ch_SM_X", "ch_TS_X",
                 "EC_FCH4", "EC_FCH4_F", "EC_FCH4_F_ANNOPTLM",
                 "EC_USTAR", "EC_GPP_NT", "EC_RECO_NT", "EC_GPP_DT", "EC_RECO_DT", "EC_NEE_F", "EC_H_F", "EC_LE_F", "EC_SW_IN_F", "EC_SW_OUT_F", "EC_LW_IN_F", 
                 "EC_LW_OUT_F", "EC_NETRAD_F", "EC_PPFD_IN_F", "EC_VPD_F", "EC_RH_F", 
                 "EC_PA_F", "EC_TA_F", "EC_WTD_F", "EC_P_F", "EC_WS_F", 
                 "EC_LE_F_ANNOPTLM", "EC_NEE_F_ANNOPTLM", "EC_FCH4_F_ANNOPTLM_UNC",
                 "u.wind", "v.wind", "u.wind.F", "v.wind.F"), 
               list(mean = mean,
                    median = median), 
               na.rm = TRUE)

USHO1_month_aggr_ECfilt <- as.data.frame(USHO1_month_aggr_ECfilt)

# rename the date column

colnames(USHO1_month_aggr_ECfilt)[1] <- "Year"
colnames(USHO1_month_aggr_ECfilt)[2] <- "Month"

# add the missing extra columns

USHO1_month_aggr_ECfilt$ch_method <- "auto"
USHO1_month_aggr_ECfilt$SITE <- "US-HO1"

USHO1_month_aggr_ECfilt <- USHO1_month_aggr_ECfilt %>% relocate(SITE, .before = Year)

# add dom_veg
USHO1_month_aggr_ECfilt$DOM_VEG <- "TREE"

# add veg classes
USHO1_month_aggr_ECfilt$MOSS_BROWN <- 0
USHO1_month_aggr_ECfilt$MOSS_SPHAGNUM <- 0
USHO1_month_aggr_ECfilt$AERENCHYMATOUS <- 0
USHO1_month_aggr_ECfilt$ERI_SHRUB <- 0
USHO1_month_aggr_ECfilt$TREE <- 1

# add climate
USHO1_month_aggr_ECfilt$KOPPEN <- "Dfb"

# add ecosystem type
USHO1_month_aggr_ECfilt$SITE_CLASSIFICATION <- "upland"

# calculate wind direction average
USHO1_month_aggr_ECfilt$EC_WD_AVG <- (atan2(USHO1_month_aggr_ECfilt$u.wind_mean, USHO1_month_aggr_ECfilt$v.wind_mean) *360/2/pi) +180

# # gap-filled:
USHO1_month_aggr_ECfilt$EC_WD_AVG_F <- (atan2(USHO1_month_aggr_ECfilt$u.wind.F_mean, USHO1_month_aggr_ECfilt$v.wind.F_mean) *360/2/pi) +180

# save as .csv

write.csv(USHO1_month_aggr_ECfilt, "path/US_HO1_chec_monthaggr_ECfilt.csv", row.names = FALSE)


# ANNUAL AGGREGATION

USHO1_yr_aggr_ECfilt <- USHO1_ECfilt %>% 
  group_by(year(USHO1_ECfilt$TIMESTAMP_START)) %>%
  summarise_at(c("ch_FCH4_umolCH4m2s1", "ch_SM_X", "ch_TS_X",  
                 "EC_FCH4", "EC_FCH4_F", "EC_FCH4_F_ANNOPTLM",
                 "EC_USTAR", "EC_GPP_NT", "EC_RECO_NT", "EC_GPP_DT", "EC_RECO_DT", "EC_NEE_F", "EC_H_F", "EC_LE_F", "EC_SW_IN_F", "EC_SW_OUT_F", "EC_LW_IN_F", 
                 "EC_LW_OUT_F", "EC_NETRAD_F", "EC_PPFD_IN_F", "EC_VPD_F", "EC_RH_F", 
                 "EC_PA_F", "EC_TA_F", "EC_WTD_F", "EC_P_F", "EC_WS_F", 
                 "EC_LE_F_ANNOPTLM", "EC_NEE_F_ANNOPTLM", "EC_FCH4_F_ANNOPTLM_UNC",
                 "u.wind", "v.wind", "u.wind.F", "v.wind.F"), 
               list(mean = mean,
                    median = median), 
               na.rm = TRUE)

USHO1_yr_aggr_ECfilt <- as.data.frame(USHO1_yr_aggr_ECfilt)

# rename the date column

colnames(USHO1_yr_aggr_ECfilt)[1] <- "Year"

# add the missing extra columns

USHO1_yr_aggr_ECfilt$ch_method <- "auto"
USHO1_yr_aggr_ECfilt$SITE <- "US-HO1"

USHO1_yr_aggr_ECfilt <- USHO1_yr_aggr_ECfilt %>% relocate(SITE, .before = Year)

# add dom_veg
USHO1_yr_aggr_ECfilt$DOM_VEG <- "TREE"

# add veg classes
USHO1_yr_aggr_ECfilt$MOSS_BROWN <- 0
USHO1_yr_aggr_ECfilt$MOSS_SPHAGNUM <- 0
USHO1_yr_aggr_ECfilt$AERENCHYMATOUS <- 0
USHO1_yr_aggr_ECfilt$ERI_SHRUB <- 0
USHO1_yr_aggr_ECfilt$TREE <- 1

# add climate
USHO1_yr_aggr_ECfilt$KOPPEN <- "Dfb"

# add ecosystem type
USHO1_yr_aggr_ECfilt$SITE_CLASSIFICATION <- "upland"

# calculate wind direction average
USHO1_yr_aggr_ECfilt$EC_WD_AVG <- (atan2(USHO1_yr_aggr_ECfilt$u.wind_mean, USHO1_yr_aggr_ECfilt$v.wind_mean) *360/2/pi) +180

# # gap-filled:
USHO1_yr_aggr_ECfilt$EC_WD_AVG_F <- (atan2(USHO1_yr_aggr_ECfilt$u.wind.F_mean, USHO1_yr_aggr_ECfilt$v.wind.F_mean) *360/2/pi) +180

# save as .csv

write.csv(USHO1_yr_aggr_ECfilt, "path/US_HO1_chec_yraggr.csv", row.names = FALSE)

### US-LA1 ###

# read csv

# EC

USLA1_EC_d <- read.csv("path/US_LA1_EC_subset_HH_TO_DD_AGGR_ECfilt.csv")

USLA1_EC_year <- read.csv("path/US_LA1_EC_subset_HH_TO_YR_AGGR_ECfilt.csv")

USLA1_EC_week <- read.csv("path/US_LA1_EC_subset_HH_TO_WEEK_AGGR_ECfilt.csv")

USLA1_EC_month <- read.csv("path/US_LA1_EC_subset_HH_TO_MONTH_AGGR_ECfilt.csv")

USLA1_EC_month <- USLA1_EC_month %>% rename("u.wind_mean" = "EC_WD_u.wind_mean",
                                            "u.wind.F_mean" = "EC_WD_u.wind_mean_F",
                                            "v.wind_mean" = "EC_WD_v.wind_mean",
                                            "v.wind.F_mean" = "EC_WD_v.wind_mean_F"
)

USLA1_EC_month <- USLA1_EC_month %>% rename("v.wind.F_mean" = "v.wind_mean_F"
)


USLA1_EC <- USLA1_EC %>% rename("u.wind_mean" = "EC_WD_u.wind_mean",
                                "u.wind.F_mean" = "EC_WD_u.wind_mean_F",
                                "v.wind_mean" = "EC_WD_v.wind_mean",
                                "v.wind.F_mean" = "EC_WD_v.wind_mean_F"
)


# CHAMBER
USLA1_ch <- read.csv("path/US_LA1.csv")

# convert to date format
USLA1_ch$ch_Date <- as_date(USLA1_ch$ch_Date)

# rename the Date column
colnames(USLA1_ch)[1] <- "TIMESTAMP"
colnames(USLA1_ch)[2] <- "SITE"

# CONVERT CHAMBER TO nmol CH4 m-2 s-1

# unit is ug m-2 hr-1

# first convert to ug m-2 s-1 by dividing by 3600
USLA1_ch$ch_FCH4_ugCH4m2s1 <- USLA1_ch$ch_FCH4_ugCH4m2hr/3600

# convert to nmol (divide ug m2s1 by CH4 molecular weight and multiply by 1000 to convert from umol to nmol):
USLA1_ch$ch_FCH4_nmolCH4m2s1 <- (USLA1_ch$ch_FCH4_ugCH4m2s1/16.04)*1000

# remove unnecessary columns and move the new one
USLA1_ch <- subset(USLA1_ch, select=-c(ch_FCH4_ugCH4m2s1))
USLA1_ch <- USLA1_ch %>% relocate(ch_FCH4_nmolCH4m2s1, .after = ch_FCH4_ugCH4m2hr)

# combine chamber and EC to create the unaggregated dataset

USLA1_ECall <- merge(USLA1_ch, USLA1_EC_d, by = c("TIMESTAMP", "Year", "Month", "DOY", "Day"), all=TRUE, sort = FALSE) 

USLA1_ECall$ch_method <- "manual"
USLA1_ECall$SITE <- "US-LA1"

USLA1_ECall <- USLA1_ECall %>% relocate(SITE, .before = TIMESTAMP)

# drop rows were ch is NA

USLA1_ECfilt <- USLA1_ECall %>% mutate(across(c(EC_FCH4_mean, EC_FCH4_median, 
                                                EC_FCH4_F_mean, EC_FCH4_F_median, 
                                                EC_FCH4_F_ANNOPTLM_mean, EC_FCH4_F_ANNOPTLM_median), 
                                              ~ ifelse(ch_FCH4_nmolCH4m2s1 == "NA", NA, .)))
# drop na from ch FCH4
USLA1_ECfilt <- USLA1_ECfilt %>% drop_na(ch_FCH4_nmolCH4m2s1)

# save as .csv
write.csv(USLA1_ECfilt, "path/US_LA1_combined_ECfilt.csv", row.names = FALSE)

# DAILY AGGREGATION

USLA1_ch_d_aggr <- USLA1_ECfilt %>% 
  group_by(date(USLA1_ECfilt$TIMESTAMP)) %>%
  summarise_at(c("ch_FCH4_nmolCH4m2s1", "ch_WTL", "ch_TS_2", "ch_TS_5", "ch_TS_10", "ch_TS_15", "ch_TS_30", "ECCH_diff"), 
               list(mean = mean,
                    median = median), 
               na.rm = TRUE)

USLA1_ch_d_aggr <- as.data.frame(USLA1_ch_d_aggr)

# rename the date column

colnames(USLA1_ch_d_aggr)[1] <- "TIMESTAMP"

# convert to datetime format
USLA1_ch_d_aggr$TIMESTAMP <- as_date(USLA1_ch_d_aggr$TIMESTAMP)

# add the missing columns

USLA1_ch_d_aggr$ch_method <- "manual"
USLA1_ch_d_aggr$SITE <- "US-LA1"

USLA1_ch_d_aggr <- USLA1_ch_d_aggr %>% relocate(SITE, .before = TIMESTAMP)

# combine with the new dd-aggregated EC data

USLA1_d_aggr_ECall <- merge(USLA1_ch_d_aggr, USLA1_EC_d, by=c("SITE", "TIMESTAMP", "ch_method"),all.x = TRUE, all.y=TRUE) 

# drop na

USLA1_d_aggr_ECfilt <- USLA1_d_aggr_ECall %>% drop_na(ch_FCH4_nmolCH4m2s1_mean)

# add dom_veg
USLA1_d_aggr_ECfilt$DOM_VEG <- "AERENCHYMATOUS"

# add veg classes
USLA1_d_aggr_ECfilt$MOSS_BROWN <- 0
USLA1_d_aggr_ECfilt$MOSS_SPHAGNUM <- 0
USLA1_d_aggr_ECfilt$AERENCHYMATOUS <- 1
USLA1_d_aggr_ECfilt$ERI_SHRUB <- 0
USLA1_d_aggr_ECfilt$TREE <- 0

# add climate
USLA1_d_aggr_ECfilt$KOPPEN <- "Cfa"

# add ecosystem type
USLA1_d_aggr_ECfilt$SITE_CLASSIFICATION <- "salt marsh"

# add EC sensor info
USLA1_d_aggr_ECfilt$EC_SENSOR <- "open-path"

# save as. csv
write.csv(USLA1_d_aggr_ECfilt, "path/US_LA1_chec_daggr.csv", row.names = FALSE)


# WEEKLY AGGREGATION

USLA1_ch_week_aggr_ECfilt <- USLA1_ECfilt %>% 
  group_by(year(USLA1_ECfilt$TIMESTAMP), isoweek(USLA1_ECfilt$TIMESTAMP)) %>%
  summarise_at(c("ch_FCH4_nmolCH4m2s1", "ch_WTL", "ch_TS_2", "ch_TS_5", "ch_TS_10", "ch_TS_15", "ch_TS_30"), 
               list(mean = mean,
                    median = median), 
               na.rm = TRUE)

USLA1_ch_week_aggr_ECfilt <- as.data.frame(USLA1_ch_week_aggr_ECfilt)

# rename the date column

colnames(USLA1_ch_week_aggr_ECfilt)[1] <- "Year"
colnames(USLA1_ch_week_aggr_ECfilt)[2] <- "Week_of_year"

# add the missing extra columns

USLA1_ch_week_aggr_ECfilt$ch_method <- "manual"
USLA1_ch_week_aggr_ECfilt$SITE <- "US-LA1"

USLA1_ch_week_aggr_ECfilt <- USLA1_ch_week_aggr_ECfilt %>% relocate(SITE, .before = Year)

# combine chamber and EC data
USLA1_week_aggr_ECfilt <- merge(USLA1_ch_week_aggr_ECfilt,USLA1_EC_week,by=c("SITE", "Year", "Week_of_year","ch_method"),all.x = TRUE, all.y=TRUE) 

# add dom_veg
USLA1_week_aggr_ECfilt$DOM_VEG <- "AERENCHYMATOUS"

# add veg classes
USLA1_week_aggr_ECfilt$MOSS_BROWN <- 0
USLA1_week_aggr_ECfilt$MOSS_SPHAGNUM <- 0
USLA1_week_aggr_ECfilt$AERENCHYMATOUS <- 1
USLA1_week_aggr_ECfilt$ERI_SHRUB <- 0
USLA1_week_aggr_ECfilt$TREE <- 0

# add climate
USLA1_week_aggr_ECfilt$KOPPEN <- "Cfa"

# add ecosystem type
USLA1_week_aggr_ECfilt$SITE_CLASSIFICATION <- "salt marsh"

# drop NA
USLA1_week_aggr_ECfilt <- USLA1_week_aggr_ECfilt %>% drop_na(ch_FCH4_nmolCH4m2s1_median)

# move SITE 
USLA1_week_aggr_ECfilt <- USLA1_week_aggr_ECfilt %>% relocate(SITE, .before = Year)

# add month
USLA1_week_aggr_ECfilt$Month <- month(as.Date(paste(USLA1_week_aggr_ECfilt$Year, USLA1_week_aggr_ECfilt$Week_of_year, 1, sep = "-"), format = "%Y-%U-%u"))

# rename the u and v component columns
USLA1_week_aggr_ECfilt <- USLA1_week_aggr_ECfilt %>% rename("u.wind_mean" = "EC_WD_u.wind_mean",
                                                            "u.wind.F_mean" = "EC_WD_u.wind_mean_F",
                                                            "v.wind_mean" = "EC_WD_v.wind_mean",
                                                            "v.wind.F_mean" = "EC_WD_v.wind_mean_F"
)

# save as .csv
write.csv(USLA1_week_aggr_ECfilt, "path/US_LA1_chec_weekaggr.csv", row.names = FALSE)


# MONTHLY AGGREGATION

USLA1_ch_month_aggr_ECfilt <- USLA1_ECfilt %>% 
  group_by(year(USLA1_ECfilt$TIMESTAMP), month(USLA1_ECfilt$TIMESTAMP)) %>%
  summarise_at(c("ch_FCH4_nmolCH4m2s1", "ch_WTL", "ch_TS_2", "ch_TS_5", "ch_TS_10", "ch_TS_15", "ch_TS_30"), 
               list(mean = mean,
                    median = median), 
               na.rm = TRUE)

USLA1_ch_month_aggr_ECfilt <- as.data.frame(USLA1_ch_month_aggr_ECfilt)

# rename the date column

colnames(USLA1_ch_month_aggr_ECfilt)[1] <- "Year"
colnames(USLA1_ch_month_aggr_ECfilt)[2] <- "Month"

# add the missing extra columns

USLA1_ch_month_aggr_ECfilt$ch_method <- "manual"
USLA1_ch_month_aggr_ECfilt$SITE <- "US-LA1"

USLA1_ch_month_aggr_ECfilt <- USLA1_ch_month_aggr_ECfilt %>% relocate(SITE, .before = Year)

# combine chamber and EC data
USLA1_month_aggr_ECfilt <- merge(USLA1_ch_month_aggr_ECfilt,USLA1_EC_month,by=c("SITE", "Year", "Month","ch_method"),all.x = TRUE, all.y=TRUE) 

# add dom_veg
USLA1_month_aggr_ECfilt$DOM_VEG <- "AERENCHYMATOUS"

# add veg classes
USLA1_month_aggr_ECfilt$MOSS_BROWN <- 0
USLA1_month_aggr_ECfilt$MOSS_SPHAGNUM <- 0
USLA1_month_aggr_ECfilt$AERENCHYMATOUS <- 1
USLA1_month_aggr_ECfilt$ERI_SHRUB <- 0
USLA1_month_aggr_ECfilt$TREE <- 0

# add climate
USLA1_month_aggr_ECfilt$KOPPEN <- "Cfa"

# add ecosystem type
USLA1_month_aggr_ECfilt$SITE_CLASSIFICATION <- "salt marsh"

# save as .csv
write.csv(USLA1_month_aggr_ECfilt, "path/US_LA1_chec_monthaggr.csv", row.names = FALSE)

# ANNUAL AGGREGATION

USLA1_ch_yr_aggr_ECfilt <- USLA1_ECfilt %>% 
  group_by(year(USLA1_ECfilt$TIMESTAMP)) %>%
  summarise_at(c("ch_FCH4_nmolCH4m2s1", "ch_WTL", "ch_TS_2", "ch_TS_5", "ch_TS_10", "ch_TS_15", "ch_TS_30"), 
               list(mean = mean,
                    median = median), 
               na.rm = TRUE)

USLA1_ch_yr_aggr_ECfilt <- as.data.frame(USLA1_ch_yr_aggr_ECfilt)

colnames(USLA1_ch_yr_aggr_ECfilt)[1] <- "Year"

# add the missing columns

USLA1_ch_yr_aggr_ECfilt$ch_method <- "manual"
USLA1_ch_yr_aggr_ECfilt$SITE <- "US-LA1"

USLA1_ch_yr_aggr_ECfilt <- USLA1_ch_yr_aggr_ECfilt %>% relocate(SITE, .before = Year)

# combine chamber and EC data
USLA1_yr_aggr_ECfilt <- merge(USLA1_ch_yr_aggr_ECfilt, USLA1_EC_year ,by=c("SITE", "Year", "ch_method"), all.x = TRUE, all.y=TRUE) 

# add dom_veg
USLA1_yr_aggr_ECfilt$DOM_VEG <- "AERENCHYMATOUS"

# add veg classes
USLA1_yr_aggr_ECfilt$MOSS_BROWN <- 0
USLA1_yr_aggr_ECfilt$MOSS_SPHAGNUM <- 0
USLA1_yr_aggr_ECfilt$AERENCHYMATOUS <- 1
USLA1_yr_aggr_ECfilt$ERI_SHRUB <- 0
USLA1_yr_aggr_ECfilt$TREE <- 0

# add climate
USLA1_yr_aggr_ECfilt$KOPPEN <- "Cfa"

# add ecosystem type
USLA1_yr_aggr_ECfilt$SITE_CLASSIFICATION <- "salt marsh"

# save as .csv

write.csv(USLA1_yr_aggr_ECfilt, "path/US_LA1_chec_yraggr_ECfilt.csv", row.names = FALSE)


### US-LA2 ###

# read csv

# EC

USLA2_EC_d <- read.csv("path/US_LA2_EC_subset_HH_TO_DD_AGGR_ECfilt.csv")

USLA2_EC_week <- read.csv("path/US_LA2_EC_subset_HH_TO_WEEK_AGGR_ECfilt.csv")

USLA2_EC_month <- read.csv("path/US_LA2_EC_subset_HH_TO_MONTH_AGGR_ECfilt.csv")

USLA2_EC_year <- read.csv("path/US_LA2_EC_subset_HH_TO_YR_AGGR_ECfilt.csv")

USLA2_EC_d <- USLA2_EC_d %>% rename("u.wind_mean" = "EC_WD_u.wind_mean",
                                "u.wind.F_mean" = "EC_WD_u.wind_mean_F",
                                "v.wind_mean" = "EC_WD_v.wind_mean",
                                "v.wind.F_mean" = "EC_WD_v.wind_mean_F"
)

USLA2_EC_month <- USLA2_EC_month %>% rename("u.wind_mean" = "EC_WD_u.wind_mean",
                                            "u.wind.F_mean" = "EC_WD_u.wind_mean_F",
                                            "v.wind_mean" = "EC_WD_v.wind_mean",
                                            "v.wind.F_mean" = "EC_WD_v.wind_mean_F"
)

# CHAMBER
USLA2_ch <- read.csv("path/US_LA2.csv")

# convert to date format
USLA2_ch$ch_Date <- as_date(USLA2_ch$ch_Date)

# CONVERT CHAMBER TO nmol CH4 m-2 s-1
# unit is ug m-2 hr-1

# first convert to ug m-2 s-1 by dividing by 3600
USLA2_ch$ch_FCH4_ugCH4m2s1 <- USLA2_ch$ch_FCH4_ugCH4m2hr/3600

# convert to nanomoles (divide ug m2s1 by CH4 molecular weight and multiply by 1000 to convert from umol to nmol):
USLA2_ch$ch_FCH4_nmolCH4m2s1 <- (USLA2_ch$ch_FCH4_ugCH4m2s1/16.04)*1000

# remove unnecessary columns and move the new one
USLA2_ch <- subset(USLA2_ch, select=-c(ch_FCH4_ugCH4m2s1))
USLA2_ch <- USLA2_ch %>% relocate(ch_FCH4_nmolCH4m2s1, .after = ch_FCH4_ugCH4m2hr)

# combine to create an unaggregated dataset
USLA2_ECall <- merge(USLA2_ch, USLA2_EC_d, by = c("SITE","TIMESTAMP"), all=TRUE, sort = FALSE) 

USLA2_ECall$ch_method <- "manual"
USLA2_ECall$SITE <- "US-LA2"

USLA2_ECall <- USLA2_ECall %>% relocate(SITE, .before = TIMESTAMP)

# drop rows were ch is NA

USLA2_ECfilt <- USLA2_ECall %>% mutate(across(c(EC_FCH4_mean, EC_FCH4_median, EC_FCH4_sd, EC_FCH4_IQR, 
                                                EC_FCH4_F_mean, EC_FCH4_F_median, EC_FCH4_F_sd, EC_FCH4_F_IQR, 
                                                EC_FCH4_F_ANNOPTLM_mean, EC_FCH4_F_ANNOPTLM_median, EC_FCH4_F_ANNOPTLM_sd, EC_FCH4_F_ANNOPTLM_IQR), 
                                              ~ ifelse(ch_FCH4_nmolCH4m2s1 == "NA", NA, .)))
# save as .csv
write.csv(USLA2_ECfilt, "path/US_LA2_combined_ECfilt.csv", row.names = FALSE)

USLA2_ECfilt <- read.csv("path/US_LA2_combined_ECfilt.csv")

USLA2_ECfilt <- USLA2_ECfilt %>% drop_na(ch_FCH4_nmolCH4m2s1)

# DAILY AGGREGATION

USLA2_ch_d_aggr <- USLA2_ECfilt %>% 
  group_by(date(USLA2_ECfilt$TIMESTAMP)) %>%
  summarise_at(c("ch_FCH4_nmolCH4m2s1", "ch_WTL", "ch_TS_2", "ch_TS_5", "ch_TS_10", "ch_TS_15", "ch_TS_30"), 
               list(mean = mean,
                    median = median),
               na.rm = TRUE)

USLA2_ch_d_aggr <- as.data.frame(USLA2_ch_d_aggr)

# rename the date column

colnames(USLA2_ch_d_aggr)[1] <- "TIMESTAMP"

# convert to datetime format
USLA2_ch_d_aggr$TIMESTAMP <- as_date(USLA2_ch_d_aggr$TIMESTAMP)

# add the missing columns

USLA2_ch_d_aggr$ch_method <- "manual"
USLA2_ch_d_aggr$SITE <- "US-LA2"

USLA2_ch_d_aggr <- USLA2_ch_d_aggr %>% relocate(SITE, .before = TIMESTAMP)

# combine with the new dd-aggregated EC data

USLA2_d_aggr_ECall <- merge(USLA2_ch_d_aggr, USLA2_EC_d, by=c("SITE", "TIMESTAMP", "ch_method"), all.x = TRUE, all.y=TRUE) 

# one row with NAs in EC column --> remove

USLA2_d_aggr_ECall <- USLA2_d_aggr_ECall %>% drop_na(EC_FCH4_F_ANNOPTLM_median)

# drop na

USLA2_d_aggr_ECfilt <- USLA2_d_aggr_ECall %>% drop_na(ch_FCH4_nmolCH4m2s1_mean)

# add dom_veg
USLA2_d_aggr_ECfilt$DOM_VEG <- "AERENCHYMATOUS"

# add veg classes
USLA2_d_aggr_ECfilt$MOSS_BROWN <- 0
USLA2_d_aggr_ECfilt$MOSS_SPHAGNUM <- 0
USLA2_d_aggr_ECfilt$AERENCHYMATOUS <- 1
USLA2_d_aggr_ECfilt$ERI_SHRUB <- 0
USLA2_d_aggr_ECfilt$TREE <- 0

# add climate
USLA2_d_aggr_ECfilt$KOPPEN <- "Cfa"

# add ecosystem type
USLA2_d_aggr_ECfilt$SITE_CLASSIFICATION <- "marsh"

write.csv(USLA2_d_aggr_ECfilt, "path/US_LA2_chec_daggr_ECfilt.csv", row.names = FALSE)

# WEEKLY AGGREGATION

USLA2_ch_week_aggr_ECfilt <- USLA2_ECfilt %>% 
  group_by(year(USLA2_ECfilt$TIMESTAMP), isoweek(USLA2_ECfilt$TIMESTAMP)) %>%
  summarise_at(c("ch_FCH4_nmolCH4m2s1", "ch_WTL", "ch_TS_2", "ch_TS_5", "ch_TS_10", "ch_TS_15", "ch_TS_30"), 
               list(mean = mean,
                    median = median), 
               na.rm = TRUE)

USLA2_ch_week_aggr_ECfilt <- as.data.frame(USLA2_ch_week_aggr_ECfilt)


# rename the date column

colnames(USLA2_ch_week_aggr_ECfilt)[1] <- "Year"
colnames(USLA2_ch_week_aggr_ECfilt)[2] <- "Week_of_year"

# add the missing extra columns

USLA2_ch_week_aggr_ECfilt$ch_method <- "manual"
USLA2_ch_week_aggr_ECfilt$SITE <- "US-LA2"

USLA2_ch_week_aggr_ECfilt <- USLA2_ch_week_aggr_ECfilt %>% relocate(SITE, .before = Year)

# combine chamber and EC data
USLA2_week_aggr_ECfilt <- merge(USLA2_ch_week_aggr_ECfilt,USLA2_EC_week,by=c("SITE", "Year", "Week_of_year","ch_method"),all.x = TRUE, all.y=TRUE) 

# add dom_veg
USLA2_week_aggr_ECfilt$DOM_VEG <- "AERENCHYMATOUS"

# add veg classes
USLA2_week_aggr_ECfilt$MOSS_BROWN <- 0
USLA2_week_aggr_ECfilt$MOSS_SPHAGNUM <- 0
USLA2_week_aggr_ECfilt$AERENCHYMATOUS <- 1
USLA2_week_aggr_ECfilt$ERI_SHRUB <- 0
USLA2_week_aggr_ECfilt$TREE <- 0

# add climate
USLA2_week_aggr_ECfilt$KOPPEN <- "Cfa"

# add ecosystem type
USLA2_week_aggr_ECfilt$SITE_CLASSIFICATION <- "marsh"

# add month
 
USLA2_week_aggr_ECfilt$Month <- month(as.Date(paste(USLA2_week_aggr_ECfilt$Year, USLA2_week_aggr_ECfilt$Week_of_year, 1, sep = "-"), format = "%Y-%U-%u"))

# fix wrong months

USLA2_week_aggr_ECfilt$Month[USLA2_week_aggr_ECfilt$Month == 11] <- 10
USLA2_week_aggr_ECfilt$Month[USLA2_week_aggr_ECfilt$Week_of_year == 26] <- 6

# rename the u and v component columns
USLA2_week_aggr_ECfilt <- USLA2_week_aggr_ECfilt %>% rename("u.wind_mean" = "EC_WD_u.wind_mean",
                                                            "u.wind.F_mean" = "EC_WD_u.wind_mean_F",
                                                            "v.wind_mean" = "EC_WD_v.wind_mean",
                                                            "v.wind.F_mean" = "EC_WD_v.wind_mean_F"
)

# save as .csv
write.csv(USLA2_week_aggr_ECfilt, "path/US_LA2_chec_weekaggr_ECfilt.csv", row.names = FALSE)

# MONTHLY AGGREGATION

USLA2_ch_month_aggr_ECfilt <- USLA2_ECfilt %>% 
  group_by(year(USLA2_ECfilt$TIMESTAMP), month(USLA2_ECfilt$TIMESTAMP)) %>%
  summarise_at(c("ch_FCH4_nmolCH4m2s1", "ch_WTL", "ch_TS_2", "ch_TS_5", "ch_TS_10", "ch_TS_15", "ch_TS_30"), 
               list(mean = mean,
                    median = median), 
               na.rm = TRUE)

USLA2_ch_month_aggr_ECfilt <- as.data.frame(USLA2_ch_month_aggr_ECfilt)

# rename the date column
colnames(USLA2_ch_month_aggr_ECfilt)[1] <- "Year"
colnames(USLA2_ch_month_aggr_ECfilt)[2] <- "Month"

# add the missing extra columns

USLA2_ch_month_aggr_ECfilt$ch_method <- "manual"
USLA2_ch_month_aggr_ECfilt$SITE <- "US-LA2"

USLA2_ch_month_aggr_ECfilt <- USLA2_ch_month_aggr_ECfilt %>% relocate(SITE, .before = Year)

# combine chamber and EC data
USLA2_month_aggr_ECfilt <- merge(USLA2_ch_month_aggr_ECfilt,USLA2_EC_month,by=c("SITE", "Year", "Month","ch_method"),all.x = TRUE, all.y=TRUE) 

# add dom_veg
USLA2_month_aggr_ECfilt$DOM_VEG <- "AERENCHYMATOUS"

# add veg classes
USLA2_month_aggr_ECfilt$MOSS_BROWN <- 0
USLA2_month_aggr_ECfilt$MOSS_SPHAGNUM <- 0
USLA2_month_aggr_ECfilt$AERENCHYMATOUS <- 1
USLA2_month_aggr_ECfilt$ERI_SHRUB <- 0
USLA2_month_aggr_ECfilt$TREE <- 0

# add climate
USLA2_month_aggr_ECfilt$KOPPEN <- "Cfa"

# add ecosystem type
USLA2_month_aggr_ECfilt$SITE_CLASSIFICATION <- "marsh"

# save as .csv
write.csv(USLA2_month_aggr_ECfilt, "path/US_LA2_chec_monthaggr_ECfilt.csv", row.names = FALSE)

# ANNUAL AGGREGATION

USLA2_ch_yr_aggr_ECfilt <- USLA2_ECfilt %>% 
  group_by(year(USLA2_ECfilt$TIMESTAMP)) %>%
  summarise_at(c("ch_FCH4_nmolCH4m2s1", "ch_WTL", "ch_TS_2", "ch_TS_5", "ch_TS_10", "ch_TS_15", "ch_TS_30"), 
               list(mean = mean,
                    median = median), 
               na.rm = TRUE)

USLA2_ch_yr_aggr_ECfilt <- as.data.frame(USLA2_ch_yr_aggr_ECfilt)

colnames(USLA2_ch_yr_aggr_ECfilt)[1] <- "Year"

# add the missing columns

USLA2_ch_yr_aggr_ECfilt$ch_method <- "manual"
USLA2_ch_yr_aggr_ECfilt$SITE <- "US-LA2"

USLA2_ch_yr_aggr_ECfilt <- USLA2_ch_yr_aggr_ECfilt %>% relocate(SITE, .before = Year)

# combine chamber and EC data
USLA2_yr_aggr_ECfilt <- merge(USLA2_ch_yr_aggr_ECfilt,USLA2_EC_year,by=c("SITE", "Year", "ch_method"),all.x = TRUE, all.y=TRUE) 

# save as .csv

write.csv(USLA2_yr_aggr_ECfilt, "path/US_LA2_chec_yraggr.csv", row.names = FALSE)


### US-LOS ###

# read csv

# EC
USLOS_EC_d <- read.csv("path/US_LOS_EC_subset.csv")

# convert to datetime
USLOS_EC_d$EC_TIMESTAMP_START <- as_datetime(USLOS_EC_d$EC_TIMESTAMP_START)
USLOS_EC_d$EC_TIMESTAMP_END <- as_datetime(USLOS_EC_d$EC_TIMESTAMP_END)

# CHAMBER
USLOS_ch_d <- read.csv("path/US_LOS.csv")

USLOS_ch_d$ch_Datetime <- as_datetime(USLOS_ch_d$ch_Datetime)

# combine the data frames using the SQL:

# set local environment time zone to UTC so the time stamps will match
Sys.setenv(TZ = "UTC")

USLOS_all <- sqldf('SELECT ch_Datetime, EC_TIMESTAMP_START, EC_TIMESTAMP_END, ch_ID, ch_FCH4_nmolm2s1, EC_FCH4, EC_FCH4_F, 
ch_AirTemp_C, ch_AtmCO2_ppm, ch_AtmCH4_LGR, 
ch_TS_X, ch_TS_2, ch_TS_5, ch_TS_10, ch_TS_15, ch_TS_30, ch_SM_X, ch_WTL, ch_method, ch_LC_class, ch_Collar_type,
EC_NEE, EC_H, EC_LE, EC_USTAR, EC_SW_IN, EC_SW_OUT, EC_LW_IN, EC_LW_OUT, EC_NETRAD, EC_PPFD_IN, EC_VPD, EC_TA, EC_P, EC_TS_1, EC_TS_2, EC_TS_3, EC_TS_4, EC_TS_5,
EC_G, EC_WTD, EC_WD, EC_GPP_NT, EC_RECO_NT, EC_GPP_DT, EC_RECO_DT, EC_WS, EC_NEE_F, EC_H_F, EC_LE_F, EC_SW_IN_F, EC_SW_OUT_F, EC_LW_IN_F, EC_LW_OUT_F, EC_NETRAD_F, EC_PPFD_IN_F, EC_VPD_F,
EC_PA_F, EC_TA_F, EC_WTD_F, EC_P_F, EC_G_F, EC_WS_F, EC_LE_F_ANNOPTLM, EC_NEE_F_ANNOPTLM, EC_FCH4_F_ANNOPTLM,
EC_FCH4_F_RANDUNC, EC_FCH4_F_ANNOPTLM_UNC, EC_FCH4_F_ANNOPTLM_QC
      FROM USLOS_EC_d 
      LEFT JOIN USLOS_ch_d ON ch_Datetime BETWEEN EC_TIMESTAMP_START and EC_TIMESTAMP_END')

# add Site column
USLOS_all$SITE <- "US-LOS"

# move to the front
USLOS_all <- USLOS_all %>%
  select(SITE, everything())

# chamber FCH4 is already in nmol m-2 s-1

write.csv(USLOS_all, "path/US_LOS_combined.csv", row.names = FALSE)

# read file again
USLOS <- read.csv("path/US_LOS_combined.csv")

# remove h_AirTemp_C and ch_AtmCO2_ppm and ch_AtmCH4_LGR

USLOS <- subset(USLOS, select=-c(ch_AirTemp_C, ch_AtmCO2_ppm, ch_AtmCH4_LGR))

# SET EC GAP-FILLED DATA AND CHAMBER DATA TO NA WHEN TIME GAP LONGER THAN 2 MONTHS

USLOS_new <- USLOS %>% mutate(across(c(ch_FCH4_nmolm2s1, ch_SM_X, ch_TS_X, ch_TS_2, ch_TS_5, ch_TS_10, ch_TS_15, ch_TS_30,ch_WTL,
                                       EC_FCH4_F, EC_NEE_F, EC_H_F, EC_LE_F, EC_SW_IN_F, EC_SW_OUT_F,EC_LW_IN_F, EC_LW_OUT_F, EC_NETRAD_F, EC_PPFD_IN_F,
                                       EC_VPD_F,EC_PA_F, EC_TA_F, EC_P_F,EC_WTD_F,EC_G_F,
                                       EC_WS_F, EC_LE_F_ANNOPTLM, EC_NEE_F_ANNOPTLM, EC_FCH4_F_ANNOPTLM, 
                                       EC_FCH4_F_RANDUNC, EC_FCH4_F_ANNOPTLM_UNC), ~ ifelse(EC_FCH4_F_ANNOPTLM_QC == 3, NA, .)))

# convert datetimes to datetime format
USLOS_new$EC_TIMESTAMP_START <- as_datetime(USLOS_new$EC_TIMESTAMP_START)
USLOS_new$EC_TIMESTAMP_END <- as_datetime(USLOS_new$EC_TIMESTAMP_END)
USLOS_new$ch_Datetime <- as_datetime(USLOS_new$ch_Datetime)

# move EC_FCH4_F_ANNOPTLM next to other flux columns
USLOS_new <- USLOS_new %>% relocate(EC_FCH4_F_ANNOPTLM, .after = EC_FCH4_F)

# set EC data to NA when chamber data is NA
USLOS_ECfilt <- USLOS_new %>% mutate(across(c(EC_FCH4, EC_FCH4_F, EC_FCH4_F_ANNOPTLM), 
                                           ~ ifelse(ch_FCH4_nmolm2s1 == "NA", NA, .)))

colnames(USLOS_ECfilt)[7] <- "ch_FCH4_nmolCH4m2s1"

write.csv(USLOS_ECfilt, "path/US_LOS_combined_ECfilt.csv", row.names = FALSE)


# chambers aggregated to EC timestamp

USLOS_ch_aggr <- USLOS_ECfilt %>% 
  group_by(across(-c(ch_Datetime, ch_FCH4_nmolCH4m2s1, ch_ID,
                     ch_TS_2, ch_TS_X, ch_TS_5, ch_TS_10, ch_TS_15, ch_TS_30, ch_SM_X, ch_WTL,
                     ch_Collar_type, ch_LC_class))) %>%
  summarise_at(c("ch_FCH4_nmolCH4m2s1", "ch_TS_X", "ch_SM_X"), 
               list(mean = mean,
                    median = median), 
               na.rm = TRUE)

USLOS_ch_aggr <- as.data.frame(USLOS_ch_aggr)

# remove rows where EC_FCH4_F_ANNOPTLM = NA

USLOS_ch_aggr <- USLOS_ch_aggr %>% drop_na(EC_FCH4_F_ANNOPTLM)

# move the new columns before EC_FCH4
USLOS_ch_aggr <- USLOS_ch_aggr %>% relocate(ch_FCH4_nmolCH4m2s1_mean, ch_FCH4_nmolCH4m2s1_median, .before = EC_FCH4)

write.csv(USLOS_ch_aggr, "path/US_LOS_combined_ch_datetimeaggr.csv", row.names = FALSE)

# DAILY AGGREGATION

# convert date to datetime
USLOS_ECfilt$EC_TIMESTAMP_START <- as_datetime(USLOS_ECfilt$EC_TIMESTAMP_START)

# set duplicated values to NA in EC columns
# create vector with EC column names
#rename EC_TIMESTAMP_START and END
colnames(USLOS_ECfilt)[c(3,4)] <- c("TIMESTAMP_START", "TIMESTAMP_END")

cols <- USLOS_ECfilt %>% dplyr::select(starts_with("EC_")) %>% colnames()

# convert duplicates to NA in EC columns
USLOS_ECfilt[cols] <- sapply(USLOS_ECfilt[cols], function(x) 
  ave(x, as_datetime(USLOS_ECfilt$TIMESTAMP_START), FUN = function(x) replace(x, duplicated(x), NA)))

# convert character and logi to numeric
cols.num <- USLOS_ECfilt %>% dplyr::select(starts_with(c("EC_", "ch_TS_"))) %>% colnames()
USLOS_ECfilt[cols.num] <- sapply(USLOS_ECfilt[cols.num],as.numeric)
USLOS_ECfilt$ch_WTL <- as.numeric(USLOS_ECfilt$ch_WTL)
USLOS_ECfilt <- USLOS_ECfilt %>% drop_na(ch_FCH4_nmolCH4m2s1)

# Calculate the u and v wind components
USLOS_ECfilt$u.wind <- -USLOS_ECfilt$EC_WS *sin(2*pi *USLOS_ECfilt$EC_WD/360)
USLOS_ECfilt$v.wind <- -USLOS_ECfilt$EC_WS *cos(2*pi *USLOS_ECfilt$EC_WD/360)
# # gap-filled (no gap-filled WD)
USLOS_ECfilt$u.wind.F <- -USLOS_ECfilt$EC_WS_F *sin(2*pi *USLOS_ECfilt$EC_WD/360)
USLOS_ECfilt$v.wind.F <- -USLOS_ECfilt$EC_WS_F *cos(2*pi *USLOS_ECfilt$EC_WD/360)

# aggregate
USLOS_d_aggr_ECfilt <- USLOS_ECfilt %>% 
  group_by(date(USLOS_ECfilt$TIMESTAMP_START)) %>%
  summarise_at(c("ch_FCH4_nmolCH4m2s1", "ch_SM_X", "ch_TS_X", 
                 "EC_FCH4", "EC_FCH4_F", "EC_FCH4_F_ANNOPTLM",
                 "EC_USTAR", "EC_GPP_NT", "EC_RECO_NT", "EC_GPP_DT", "EC_RECO_DT", "EC_NEE_F", "EC_H_F", "EC_LE_F", "EC_SW_IN_F", "EC_SW_OUT_F", "EC_LW_IN_F", 
                 "EC_LW_OUT_F", "EC_NETRAD_F", "EC_PPFD_IN_F", "EC_VPD_F",
                 "EC_PA_F", "EC_TA_F", "EC_WTD_F", "EC_P_F", "EC_WS_F", 
                 "EC_LE_F_ANNOPTLM", "EC_NEE_F_ANNOPTLM", "EC_TS_1",
                 "u.wind", "v.wind", "u.wind.F", "v.wind.F"), 
               list(mean = mean,
                    median = median), 
               na.rm = TRUE)

USLOS_d_aggr_ECfilt <- as.data.frame(USLOS_d_aggr_ECfilt)

# rename the date column

colnames(USLOS_d_aggr_ECfilt)[1] <- "DATE"

# convert to date format
USLOS_d_aggr_ECfilt$DATE <- as_date(USLOS_d_aggr_ECfilt$DATE)

# calculate wind direction average
USLOS_d_aggr_ECfilt$EC_WD_AVG <- (atan2(USLOS_d_aggr_ECfilt$u.wind_mean, USLOS_d_aggr_ECfilt$v.wind_mean) *360/2/pi) +180

# gap-filled:
USLOS_d_aggr_ECfilt$EC_WD_AVG_F <- (atan2(USLOS_d_aggr_ECfilt$u.wind.F_mean, USLOS_d_aggr_ECfilt$v.wind.F_mean) *360/2/pi) +180

# add the missing columns

USLOS_d_aggr_ECfilt$Year <- year(USLOS_d_aggr_ECfilt$DATE)
USLOS_d_aggr_ECfilt$Month <- month(USLOS_d_aggr_ECfilt$DATE)
USLOS_d_aggr_ECfilt$Day <- day(USLOS_d_aggr_ECfilt$DATE)
USLOS_d_aggr_ECfilt$DOY <- yday(USLOS_d_aggr_ECfilt$DATE)

USLOS_d_aggr_ECfilt$ch_method <- "manual"
USLOS_d_aggr_ECfilt$SITE <- "US-LOS"

USLOS_d_aggr_ECfilt <- USLOS_d_aggr_ECfilt %>% relocate(SITE, .before = DATE)

# remove rows with FCH4_ANNOTPLTM = NA

USLOS_d_aggr_ECfilt <- USLOS_d_aggr_ECfilt %>% drop_na(EC_FCH4_F_ANNOPTLM_mean) 

# add dom_veg
USLOS_d_aggr_ECfilt$DOM_VEG <- "ERI_SHRUB"

# add veg classes
USLOS_d_aggr_ECfilt$MOSS_BROWN <- 0
USLOS_d_aggr_ECfilt$MOSS_SPHAGNUM <- 0
USLOS_d_aggr_ECfilt$AERENCHYMATOUS <- 1
USLOS_d_aggr_ECfilt$ERI_SHRUB <- 1
USLOS_d_aggr_ECfilt$TREE <- 1

# add climate
USLOS_d_aggr_ECfilt$KOPPEN <- "Dfb"

# add ecosystem type
USLOS_d_aggr_ECfilt$SITE_CLASSIFICATION <- "fen"

# save as .csv

write.csv(USLOS_d_aggr_ECfilt, "path/US_LOS_chec_daggr.csv", row.names = FALSE)


# WEEKLY AGGREGATION

USLOS_week_aggr_ECfilt <- USLOS_ECfilt %>% 
  group_by(year(USLOS_ECfilt$TIMESTAMP_START), week(USLOS_ECfilt$TIMESTAMP_START)) %>%
  summarise_at(c("ch_FCH4_nmolCH4m2s1", "ch_SM_X", "ch_TS_X", 
                 "EC_FCH4", "EC_FCH4_F", "EC_FCH4_F_ANNOPTLM",
                 "EC_USTAR", "EC_GPP_NT", "EC_RECO_NT", "EC_GPP_DT", "EC_RECO_DT", "EC_NEE_F", "EC_H_F", "EC_LE_F", "EC_SW_IN_F", "EC_SW_OUT_F", "EC_LW_IN_F", 
                 "EC_LW_OUT_F", "EC_NETRAD_F", "EC_PPFD_IN_F", "EC_VPD_F",
                 "EC_PA_F", "EC_TA_F", "EC_WTD_F", "EC_P_F", "EC_WS_F", 
                 "EC_LE_F_ANNOPTLM", "EC_NEE_F_ANNOPTLM",
                 "u.wind", "v.wind", "u.wind.F", "v.wind.F"), 
               list(mean = mean,
                    median = median), 
               na.rm = TRUE)

USLOS_week_aggr_ECfilt <- as.data.frame(USLOS_week_aggr_ECfilt)

# rename the date column

colnames(USLOS_week_aggr_ECfilt)[1] <- "Year"
colnames(USLOS_week_aggr_ECfilt)[2] <- "Week_of_year"

# add the missing extra columns

USLOS_week_aggr_ECfilt$ch_method <- "manual"
USLOS_week_aggr_ECfilt$SITE <- "US-LOS"

USLOS_week_aggr_ECfilt <- USLOS_week_aggr_ECfilt %>% relocate(SITE, .before = Year)

# add dom_veg
USLOS_week_aggr_ECfilt$DOM_VEG <- "ERI_SHRUB"

# add veg classes
USLOS_week_aggr_ECfilt$MOSS_BROWN <- 0
USLOS_week_aggr_ECfilt$MOSS_SPHAGNUM <- 0
USLOS_week_aggr_ECfilt$AERENCHYMATOUS <- 1
USLOS_week_aggr_ECfilt$ERI_SHRUB <- 1
USLOS_week_aggr_ECfilt$TREE <- 1

# add climate
USLOS_week_aggr_ECfilt$KOPPEN <- "Dfb"

# add ecosystem type
USLOS_week_aggr_ECfilt$SITE_CLASSIFICATION <- "fen"

# calculate wind direction average
USLOS_week_aggr_ECfilt$EC_WD_AVG <- (atan2(USLOS_week_aggr_ECfilt$u.wind_mean, USLOS_week_aggr_ECfilt$v.wind_mean) *360/2/pi) +180

# # gap-filled:
USLOS_week_aggr_ECfilt$EC_WD_AVG_F <- (atan2(USLOS_week_aggr_ECfilt$u.wind.F_mean, USLOS_week_aggr_ECfilt$v.wind.F_mean) *360/2/pi) +180

# add month
USLOS_week_aggr_ECfilt$Month <- month(as.Date(paste(USLOS_week_aggr_ECfilt$Year, USLOS_week_aggr_ECfilt$Week_of_year, 1, sep = "-"), format = "%Y-%U-%u"))

USLOS_week_aggr_ECfilt <- USLOS_week_aggr_ECfilt %>% relocate(Month, .after = Week_of_year)

# week 26 is month 7, not 6
USLOS_week_aggr_ECfilt$Month[USLOS_week_aggr_ECfilt$Week_of_year == 26] <- 7

# save as .csv

write.csv(USLOS_week_aggr_ECfilt, "path/US_LOS_chec_weekaggr.csv", row.names = FALSE)


# MONTHLY AGGREGATION

USLOS_month_aggr_ECfilt <- USLOS_ECfilt %>% 
  group_by(year(USLOS_ECfilt$TIMESTAMP_START), month(USLOS_ECfilt$TIMESTAMP_START)) %>%
  summarise_at(c("ch_FCH4_nmolCH4m2s1", "ch_SM_X", "ch_TS_X", 
                 "EC_FCH4", "EC_FCH4_F", "EC_FCH4_F_ANNOPTLM",
                 "EC_USTAR", "EC_GPP_NT", "EC_RECO_NT", "EC_GPP_DT", "EC_RECO_DT", "EC_NEE_F", "EC_H_F", "EC_LE_F", "EC_SW_IN_F", "EC_SW_OUT_F", "EC_LW_IN_F", 
                 "EC_LW_OUT_F", "EC_NETRAD_F", "EC_PPFD_IN_F", "EC_VPD_F",
                 "EC_PA_F", "EC_TA_F", "EC_WTD_F", "EC_P_F", "EC_WS_F", 
                 "EC_LE_F_ANNOPTLM", "EC_NEE_F_ANNOPTLM",
                 "u.wind", "v.wind", "u.wind.F", "v.wind.F"), 
               list(mean = mean,
                    median = median), 
               na.rm = TRUE)

USLOS_month_aggr_ECfilt <- as.data.frame(USLOS_month_aggr_ECfilt)

# rename the date column

colnames(USLOS_month_aggr_ECfilt)[1] <- "Year"
colnames(USLOS_month_aggr_ECfilt)[2] <- "Month"

# add the missing extra columns

USLOS_month_aggr_ECfilt$ch_method <- "manual"
USLOS_month_aggr_ECfilt$SITE <- "US-LOS"

USLOS_month_aggr_ECfilt <- USLOS_month_aggr_ECfilt %>% relocate(SITE, .before = Year)

# add dom_veg
USLOS_month_aggr_ECfilt$DOM_VEG <- "ERI_SHRUB"

# add veg classes
USLOS_month_aggr_ECfilt$MOSS_BROWN <- 0
USLOS_month_aggr_ECfilt$MOSS_SPHAGNUM <- 0
USLOS_month_aggr_ECfilt$AERENCHYMATOUS <- 1
USLOS_month_aggr_ECfilt$ERI_SHRUB <- 1
USLOS_month_aggr_ECfilt$TREE <- 1

# add climate
USLOS_month_aggr_ECfilt$KOPPEN <- "Dfb"

# add ecosystem type
USLOS_month_aggr_ECfilt$SITE_CLASSIFICATION <- "fen"

# calculate wind direction average
USLOS_month_aggr_ECfilt$EC_WD_AVG <- (atan2(USLOS_month_aggr_ECfilt$u.wind_mean, USLOS_month_aggr_ECfilt$v.wind_mean) *360/2/pi) +180

# # gap-filled:
USLOS_month_aggr_ECfilt$EC_WD_AVG_F <- (atan2(USLOS_month_aggr_ECfilt$u.wind.F_mean, USLOS_month_aggr_ECfilt$v.wind.F_mean) *360/2/pi) +180

# save as .csv

write.csv(USLOS_month_aggr_ECfilt, "path/US_LOS_chec_monthaggr.csv", row.names = FALSE)

# ANNUAL AGGREGATION

USLOS_yr_aggr_ECfilt <- USLOS_ECfilt %>% 
  group_by(year(USLOS_ECfilt$TIMESTAMP_START)) %>%
  summarise_at(c("ch_FCH4_nmolCH4m2s1", "ch_SM_X", "ch_TS_X", 
                 "EC_FCH4", "EC_FCH4_F", "EC_FCH4_F_ANNOPTLM",
                 "EC_USTAR", "EC_GPP_NT", "EC_RECO_NT", "EC_GPP_DT", "EC_RECO_DT", "EC_NEE_F", "EC_H_F", "EC_LE_F", "EC_SW_IN_F", "EC_SW_OUT_F", "EC_LW_IN_F", 
                 "EC_LW_OUT_F", "EC_NETRAD_F", "EC_PPFD_IN_F", "EC_VPD_F",
                 "EC_PA_F", "EC_TA_F", "EC_WTD_F", "EC_P_F", "EC_WS_F", 
                 "EC_LE_F_ANNOPTLM", "EC_NEE_F_ANNOPTLM",
                 "u.wind", "v.wind", "u.wind.F", "v.wind.F"), 
               list(mean = mean,
                    median = median), 
               na.rm = TRUE)

USLOS_yr_aggr_ECfilt <- as.data.frame(USLOS_yr_aggr_ECfilt)

# rename the date column

colnames(USLOS_yr_aggr_ECfilt)[1] <- "Year"

# add the missing extra columns

USLOS_yr_aggr_ECfilt$ch_method <- "manual"
USLOS_yr_aggr_ECfilt$SITE <- "US-LOS"

USLOS_yr_aggr_ECfilt <- USLOS_yr_aggr_ECfilt %>% relocate(SITE, .before = Year)

# add dom_veg
USLOS_yr_aggr_ECfilt$DOM_VEG <- "ERI_SHRUB"

# add veg classes
USLOS_yr_aggr_ECfilt$MOSS_BROWN <- 0
USLOS_yr_aggr_ECfilt$MOSS_SPHAGNUM <- 0
USLOS_yr_aggr_ECfilt$AERENCHYMATOUS <- 1
USLOS_yr_aggr_ECfilt$ERI_SHRUB <- 1
USLOS_yr_aggr_ECfilt$TREE <- 1

# add climate
USLOS_yr_aggr_ECfilt$KOPPEN <- "Dfb"

# add ecosystem type
USLOS_yr_aggr_ECfilt$SITE_CLASSIFICATION <- "fen"

# calculate wind direction average
USLOS_yr_aggr_ECfilt$EC_WD_AVG <- (atan2(USLOS_yr_aggr_ECfilt$u.wind_mean, USLOS_yr_aggr_ECfilt$v.wind_mean) *360/2/pi) +180

# # gap-filled:
USLOS_yr_aggr_ECfilt$EC_WD_AVG_F <- (atan2(USLOS_yr_aggr_ECfilt$u.wind.F_mean, USLOS_yr_aggr_ECfilt$v.wind.F_mean) *360/2/pi) +180

# save as .csv

write.csv(USLOS_yr_aggr_ECfilt, "path/US_LOS_chec_yraggr.csv", row.names = FALSE)


### US-OWC ###

# read csv

# EC
USOWC_EC_d <- read.csv("path/US_OWC_EC_subset.csv")

# convert to datetime
USOWC_EC_d$EC_TIMESTAMP_START <- as_datetime(USOWC_EC_d$EC_TIMESTAMP_START)
USOWC_EC_d$EC_TIMESTAMP_END <- as_datetime(USOWC_EC_d$EC_TIMESTAMP_END)

# CHAMBER
USOWC_ch_d <- read.csv("path/US_OWC.csv")

USOWC_ch_d$ch_Datetime <- as_datetime(USOWC_ch_d$ch_Datetime)

# subset to EC scale --> until 2016-10-30 23:30:00
Sys.setenv(TZ = "UTC")
USOWC_ch_sub <- subset(USOWC_ch_d, ch_Datetime <= "2016-10-30 23:30:00")

# combine the data frames using SQL

# set local environment time zone to UTC so the time stamps will match
Sys.setenv(TZ = "UTC")

USOWC_all <- sqldf('SELECT ch_Datetime, EC_TIMESTAMP_START, EC_TIMESTAMP_END, ch_LOC_LATITUDE, ch_LOC_LONGITUDE, ch_ID, 
ch_FCH4_nmolCH4m2s1, EC_FCH4, EC_FCH4_F,
ch_TS_2, ch_TS_5, ch_TS_10, ch_TS_15, ch_TS_30, ch_WTL, ch_method, ch_LC_class, ch_Collar_type,
EC_NEE, EC_H, EC_LE, EC_USTAR, EC_SW_IN, EC_SW_OUT, EC_LW_IN, EC_LW_OUT, EC_NETRAD, EC_PPFD_IN, EC_VPD, EC_RH, EC_PA, EC_TA, EC_P, EC_TS_1, EC_TS_2,
EC_WTD, EC_GPP_NT, EC_RECO_NT, EC_GPP_DT, EC_RECO_DT, EC_WD, EC_WS, EC_NEE_F, EC_H_F, EC_LE_F, EC_SW_IN_F, EC_SW_OUT_F, EC_VPD_F, EC_RH_F,
EC_LW_IN_F, EC_LW_OUT_F, EC_NETRAD_F, EC_PPFD_IN_F,
EC_PA_F, EC_TA_F, EC_WTD_F, EC_P_F, EC_WS_F, EC_LE_F_ANNOPTLM, EC_FCH4_F_ANNOPTLM, EC_NEE_F_ANNOPTLM,
EC_FCH4_F_RANDUNC, EC_FCH4_F_ANNOPTLM_UNC, EC_FCH4_F_ANNOPTLM_QC
      FROM USOWC_EC_d 
      LEFT JOIN USOWC_ch_d ON ch_Datetime BETWEEN EC_TIMESTAMP_START and EC_TIMESTAMP_END')

# add Site column
USOWC_all$SITE <- "US-OWC"

# move to the front
USOWC_all <- USOWC_all %>%
  select(SITE, everything())

cols <- c("ch_LOC_LATITUDE", "ch_LOC_LONGITUDE" ,"ch_ID", "ch_FCH4_nmolCH4m2s1", "ch_TS_2", "ch_TS_5", "ch_TS_10", "ch_TS_15", "ch_TS_30", "ch_WTL", "ch_method", "ch_LC_class", "ch_Collar_type")

USOWC_all[cols] <- sapply(USOWC_all[cols], function(x) 
  ave(x, as_datetime(USOWC_all$ch_Datetime), FUN = function(x) replace(x, duplicated(x, fromLast=T), NA)))

# set the first duplicate to NA in the timestamp columns
cols2 <- "ch_Datetime"

USOWC_all[cols2] <- sapply(USOWC_all[cols2], function(x) 
  ave(x, USOWC_all$ch_Datetime, FUN = function(x) replace(x, duplicated(x, fromLast=T), NA)))

# since the datetimes are numeric, convert them back to date time
USOWC_all$ch_Datetime <- as.POSIXct(USOWC_all$ch_Datetime, origin="1970-01-01")

# subset to chamber ending date 2016-10-18 11:28:00
USOWC_all_sub <- subset(USOWC_all, date(EC_TIMESTAMP_END) <= "2016-10-18")

# chamber FCH4 data is already in nmol CH4 m-2 s-1

USOWC_all <- USOWC_all %>% drop_na(ch_FCH4_nmolCH4m2s1)

write.csv(USOWC_all_sub, "path/US_OWC_combined.csv", row.names = FALSE)

# read the file again
USOWC <- read.csv("path/US_OWC_combined.csv")

# SET EC GAP-FILLED DATA AND CHAMBER DATA TO NA WHEN TIME GAP LONGER THAN 2 MONTHS

USOWC_new <- USOWC %>% mutate(across(c(ch_FCH4_nmolCH4m2s1, ch_TS_2, ch_TS_5, ch_TS_10, ch_TS_15, ch_TS_30,ch_WTL,
                                       EC_FCH4_F, EC_NEE_F, EC_H_F, EC_LE_F, EC_SW_IN_F, EC_SW_OUT_F,
                                       EC_LW_IN_F, EC_LW_OUT_F, EC_NETRAD_F, EC_PPFD_IN_F,
                                       EC_VPD_F,EC_RH_F, EC_PA_F, EC_TA_F, EC_P_F,EC_WTD_F,
                                       EC_GPP_NT, EC_GPP_DT,
                                       EC_WS_F, EC_LE_F_ANNOPTLM, EC_FCH4_F_ANNOPTLM, EC_NEE_F_ANNOPTLM,
                                       EC_FCH4_F_RANDUNC, EC_FCH4_F_ANNOPTLM_UNC), ~ ifelse(EC_FCH4_F_ANNOPTLM_QC == 3, NA, .)))


# convert datetimes to datetime format
USOWC_new$EC_TIMESTAMP_START <- as_datetime(USOWC_new$EC_TIMESTAMP_START)
USOWC_new$EC_TIMESTAMP_END <- as_datetime(USOWC_new$EC_TIMESTAMP_END)
USOWC_new$ch_Datetime <- as_datetime(USOWC_new$ch_Datetime)

# move EC_FCH4_F_ANNOPTLM next to other flux columns
USOWC_new <- USOWC_new %>% relocate(EC_FCH4_F_ANNOPTLM, .after = EC_FCH4_F)

# there is one row where the ch_Datetime is NA --> bug in the combination
# set it manually

USOWC_new["18837", "ch_Datetime"] <- as_datetime("2016-07-12 17:00:00")
USOWC_new["18837", "ch_LOC_LONGITUDE"] <- -82.51120
USOWC_new["18837", "ch_ID"] <- 1
USOWC_new["18837", "ch_method"] <- "manual"
USOWC_new["18837", "ch_LC_class"] <- "Open Water"
USOWC_new["18837", "ch_Collar_type"] <- "floating"

# set EC data to NA when chamber data is NA
USOWC_ECfilt <- USOWC_new %>% mutate(across(c(EC_FCH4, EC_FCH4_F, EC_FCH4_F_ANNOPTLM), 
                                                           ~ ifelse(ch_FCH4_nmolCH4m2s1 == "NA", NA, .)))

USOWC_ECfilt <- USOWC_ECfilt %>% drop_na(ch_FCH4_nmolCH4m2s1)

# aggregate chamber measurements
USOWC_ch_aggr <- USOWC_ECfilt %>% 
  group_by(across(-c(ch_Datetime, ch_FCH4_nmolCH4m2s1, ch_LOC_LATITUDE, ch_LOC_LONGITUDE, ch_ID,
                     ch_TS_2, ch_TS_5, ch_TS_10, ch_TS_15, ch_TS_30, ch_WTL,
                     ch_Collar_type, ch_LC_class))) %>%
  summarise_at(c("ch_FCH4_nmolCH4m2s1"),
               list(mean = mean,
                    median = median), 
               na.rm = TRUE)

USOWC_ch_aggr <- as.data.frame(USOWC_ch_aggr)

# remove rows where EC_FCH4_F_ANNOPTLM = NA

USOWC_ch_aggr <- USOWC_ch_aggr %>% drop_na(EC_FCH4_F_ANNOPTLM)

# move the new columns before EC_FCH4
USOWC_ch_aggr <- USOWC_ch_aggr %>% relocate(ch_FCH4_nmolCH4m2s1_mean, ch_FCH4_nmolCH4m2s1_sd, ch_FCH4_nmolCH4m2s1_median, ch_FCH4_nmolCH4m2s1_IQR, .before = EC_FCH4)

# save as .csv
write.csv(USOWC_ch_aggr, "path/US_OWC_combined_ch_datetimeaggr.csv", row.names = FALSE)

# DAILY AGGREGATION

USOWC_ECfilt$EC_TIMESTAMP_START <- as_datetime(USOWC_ECfilt$EC_TIMESTAMP_START)

# set duplicated values to NA in EC columns because they affect the stats
# create vector with EC column names
#rename EC_TIMESTAMP_START and END
colnames(USOWC_ECfilt)[c(3,4)] <- c("TIMESTAMP_START", "TIMESTAMP_END")

cols <- USOWC_ECfilt %>% dplyr::select(starts_with("EC_")) %>% colnames()

# convert duplicates to NA in EC columns
USOWC_ECfilt[cols] <- sapply(USOWC_ECfilt[cols], function(x) 
  ave(x, as_datetime(USOWC_ECfilt$TIMESTAMP_START), FUN = function(x) replace(x, duplicated(x), NA)))

# convert character and logi to numeric
cols.num <- USOWC_ECfilt %>% dplyr::select(starts_with(c("EC_", "ch_TS_"))) %>% colnames()
USOWC_ECfilt[cols.num] <- sapply(USOWC_ECfilt[cols.num],as.numeric)
USOWC_ECfilt$ch_WTL <- as.numeric(USOWC_ECfilt$ch_WTL)

USOWC_ECfilt <- USOWC_ECfilt %>% drop_na(ch_FCH4_nmolCH4m2s1)

# Calculate the u and v wind components
USOWC_ECfilt$u.wind <- -USOWC_ECfilt$EC_WS *sin(2*pi *USOWC_ECfilt$EC_WD/360)
USOWC_ECfilt$v.wind <- -USOWC_ECfilt$EC_WS *cos(2*pi *USOWC_ECfilt$EC_WD/360)

# gap-filled (no gap-filled WD)
USOWC_ECfilt$u.wind.F <- -USOWC_ECfilt$EC_WS_F *sin(2*pi *USOWC_ECfilt$EC_WD/360)
USOWC_ECfilt$v.wind.F <- -USOWC_ECfilt$EC_WS_F *cos(2*pi *USOWC_ECfilt$EC_WD/360)

# aggregate
USOWC_d_aggr_ECfilt <- USOWC_ECfilt %>% 
  group_by(date(USOWC_ECfilt$TIMESTAMP_START)) %>%
  summarise_at(c("ch_FCH4_nmolCH4m2s1",
                 "EC_FCH4", "EC_FCH4_F", "EC_FCH4_F_ANNOPTLM", "EC_GPP_NT", "EC_GPP_DT", "EC_RECO_NT", "EC_RECO_DT",
                 "EC_USTAR", "EC_NEE_F", "EC_H_F", "EC_LE_F", "EC_SW_IN_F", "EC_SW_OUT_F", "EC_LW_IN_F", 
                 "EC_LW_OUT_F", "EC_NETRAD_F", "EC_PPFD_IN_F", "EC_VPD_F", "EC_RH_F",
                 "EC_PA_F", "EC_TA_F", "EC_WTD_F", "EC_P_F", "EC_WS_F", 
                 "EC_LE_F_ANNOPTLM", "EC_TS_1", "EC_NEE_F_ANNOPTLM",
                 "u.wind", "v.wind", "u.wind.F", "v.wind.F"), 
               list(mean = mean,
                    median = median), 
               na.rm = TRUE)

USOWC_d_aggr_ECfilt <- as.data.frame(USOWC_d_aggr_ECfilt)

# remove rows where CH4 ANN = NA

USOWC_d_aggr_ECfilt <- USOWC_d_aggr_ECfilt %>% drop_na(EC_FCH4_F_ANNOPTLM_mean)

# rename the date column

colnames(USOWC_d_aggr_ECfilt)[1] <- "DATE"

# convert to date format
USOWC_d_aggr_ECfilt$DATE <- as_date(USOWC_d_aggr_ECfilt$DATE)

# calculate wind direction average
USOWC_d_aggr_ECfilt$EC_WD_AVG <- (atan2(USOWC_d_aggr_ECfilt$u.wind_mean, USOWC_d_aggr_ECfilt$v.wind_mean) *360/2/pi) +180

# gap-filled:
USOWC_d_aggr_ECfilt$EC_WD_AVG_F <- (atan2(USOWC_d_aggr_ECfilt$u.wind.F_mean, USOWC_d_aggr_ECfilt$v.wind.F_mean) *360/2/pi) +180

# add the missing columns
USOWC_d_aggr_ECfilt$Year <- year(USOWC_d_aggr_ECfilt$DATE)
USOWC_d_aggr_ECfilt$Month <- month(USOWC_d_aggr_ECfilt$DATE)
USOWC_d_aggr_ECfilt$Day <- day(USOWC_d_aggr_ECfilt$DATE)
USOWC_d_aggr_ECfilt$DOY <- yday(USOWC_d_aggr_ECfilt$DATE)

USOWC_d_aggr_ECfilt$ch_method <- "manual"
USOWC_d_aggr_ECfilt$SITE <- "US-OWC"

USOWC_d_aggr_ECfilt <- USOWC_d_aggr_ECfilt %>% relocate(SITE, .before = DATE)

# add dom_veg
USOWC_d_aggr_ECfilt$DOM_VEG <- "AERENCHYMATOUS"

# add veg classes
USOWC_d_aggr_ECfilt$MOSS_BROWN <- 0
USOWC_d_aggr_ECfilt$MOSS_SPHAGNUM <- 0
USOWC_d_aggr_ECfilt$AERENCHYMATOUS <- 1
USOWC_d_aggr_ECfilt$ERI_SHRUB <- 0
USOWC_d_aggr_ECfilt$TREE <- 0

# add climate
USOWC_d_aggr_ECfilt$KOPPEN <- "Dfb"

# add ecosystem type
USOWC_d_aggr_ECfilt$SITE_CLASSIFICATION <- "marsh"

# save as .csv
write.csv(USOWC_d_aggr_ECfilt, "path/US_OWC_chec_daggr_ECfilt.csv", row.names = FALSE)

# WEEKLY AGGREGATION

USOWC_week_aggr_ECfilt <- USOWC_ECfilt %>% 
  group_by(year(USOWC_ECfilt$TIMESTAMP_START), isoweek(USOWC_ECfilt$TIMESTAMP_START)) %>%
  summarise_at(c("ch_FCH4_nmolCH4m2s1",
                 "EC_FCH4", "EC_FCH4_F", "EC_FCH4_F_ANNOPTLM", "EC_GPP_NT", "EC_GPP_DT", "EC_RECO_NT", "EC_RECO_DT",
                 "EC_USTAR", "EC_NEE_F", "EC_H_F", "EC_LE_F", "EC_SW_IN_F", "EC_SW_OUT_F", "EC_LW_IN_F", 
                 "EC_LW_OUT_F", "EC_NETRAD_F", "EC_PPFD_IN_F", "EC_VPD_F", "EC_RH_F",
                 "EC_PA_F", "EC_TA_F", "EC_WTD_F", "EC_P_F", "EC_WS_F", 
                 "EC_LE_F_ANNOPTLM", "EC_NEE_F_ANNOPTLM", "EC_TS_1",
                 "u.wind", "v.wind", "u.wind.F", "v.wind.F"), 
               list(mean = mean,
                    median = median), 
               na.rm = TRUE)

USOWC_week_aggr_ECfilt <- as.data.frame(USOWC_week_aggr_ECfilt)

# rename the date column

colnames(USOWC_week_aggr_ECfilt)[1] <- "Year"
colnames(USOWC_week_aggr_ECfilt)[2] <- "Week_of_year"

# add the missing extra columns

USOWC_week_aggr_ECfilt$ch_method <- "manual"
USOWC_week_aggr_ECfilt$SITE <- "US-OWC"

USOWC_week_aggr_ECfilt <- USOWC_week_aggr_ECfilt %>% relocate(SITE, .before = Year)

# add dom_veg
USOWC_week_aggr_ECfilt$DOM_VEG <- "AERENCHYMATOUS"

# add veg classes
USOWC_week_aggr_ECfilt$MOSS_BROWN <- 0
USOWC_week_aggr_ECfilt$MOSS_SPHAGNUM <- 0
USOWC_week_aggr_ECfilt$AERENCHYMATOUS <- 1
USOWC_week_aggr_ECfilt$ERI_SHRUB <- 0
USOWC_week_aggr_ECfilt$TREE <- 0

# add climate
USOWC_week_aggr_ECfilt$KOPPEN <- "Dfb"

# add ecosystem type
USOWC_week_aggr_ECfilt$SITE_CLASSIFICATION <- "marsh"

USOWC_week_aggr_ECfilt$EC_SENSOR <- "open-path"

# calculate wind direction average
USOWC_week_aggr_ECfilt$EC_WD_AVG <- (atan2(USOWC_week_aggr_ECfilt$u.wind_mean, USOWC_week_aggr_ECfilt$v.wind_mean) *360/2/pi) +180

# gap-filled:
USOWC_week_aggr_ECfilt$EC_WD_AVG_F <- (atan2(USOWC_week_aggr_ECfilt$u.wind.F_mean, USOWC_week_aggr_ECfilt$v.wind.F_mean) *360/2/pi) +180

# add month
USOWC_week_aggr_ECfilt$Month <- month(as.Date(paste(USOWC_week_aggr_ECfilt$Year, USOWC_week_aggr_ECfilt$Week_of_year, 1, sep = "-"), format = "%Y-%U-%u"))

USOWC_week_aggr_ECfilt <- USOWC_week_aggr_ECfilt %>% relocate(Month, .after = Week_of_year)
# weeks and months are correct

USOWC_week_aggr_ECfilt$ECCH_diff <- USOWC_week_aggr_ECfilt$EC_FCH4_F_ANNOPTLM_median - USOWC_week_aggr_ECfilt$ch_FCH4_nmolCH4m2s1_median

# save as .csv

write.csv(USOWC_week_aggr_ECfilt, "path/US_OWC_chec_weekaggr.csv", row.names = FALSE)


# MONTHLY AGGREGATION

USOWC_month_aggr_ECfilt <- USOWC_ECfilt %>% 
  group_by(year(USOWC_ECfilt$TIMESTAMP_START), month(USOWC_ECfilt$TIMESTAMP_START)) %>%
  summarise_at(c("ch_FCH4_nmolCH4m2s1",
                 "EC_FCH4", "EC_FCH4_F", "EC_FCH4_F_ANNOPTLM", "EC_GPP_NT", "EC_GPP_DT", "EC_RECO_NT", "EC_RECO_DT",
                 "EC_USTAR", "EC_NEE_F", "EC_H_F", "EC_LE_F", "EC_SW_IN_F", "EC_SW_OUT_F", "EC_LW_IN_F", 
                 "EC_LW_OUT_F", "EC_NETRAD_F", "EC_PPFD_IN_F", "EC_VPD_F", "EC_RH_F",
                 "EC_PA_F", "EC_TA_F", "EC_WTD_F", "EC_P_F", "EC_WS_F", 
                 "EC_LE_F_ANNOPTLM", "EC_NEE_F_ANNOPTLM", "EC_TS_1",
                 "u.wind", "v.wind", "u.wind.F", "v.wind.F"), 
               list(mean = mean,
                    median = median), 
               na.rm = TRUE)

USOWC_month_aggr_ECfilt <- as.data.frame(USOWC_month_aggr_ECfilt)

# rename the date column

colnames(USOWC_month_aggr_ECfilt)[1] <- "Year"
colnames(USOWC_month_aggr_ECfilt)[2] <- "Month"

# add the missing extra columns

USOWC_month_aggr_ECfilt$ch_method <- "manual"
USOWC_month_aggr_ECfilt$SITE <- "US-OWC"

USOWC_month_aggr_ECfilt <- USOWC_month_aggr_ECfilt %>% relocate(SITE, .before = Year)

# add dom_veg
USOWC_month_aggr_ECfilt$DOM_VEG <- "AERENCHYMATOUS"

# add veg classes
USOWC_month_aggr_ECfilt$MOSS_BROWN <- 0
USOWC_month_aggr_ECfilt$MOSS_SPHAGNUM <- 0
USOWC_month_aggr_ECfilt$AERENCHYMATOUS <- 1
USOWC_month_aggr_ECfilt$ERI_SHRUB <- 0
USOWC_month_aggr_ECfilt$TREE <- 0

# add climate
USOWC_month_aggr_ECfilt$KOPPEN <- "Dfb"

# add ecosystem type
USOWC_month_aggr_ECfilt$SITE_CLASSIFICATION <- "marsh"

# calculate wind direction average
USOWC_month_aggr_ECfilt$EC_WD_AVG <- (atan2(USOWC_month_aggr_ECfilt$u.wind_mean, USOWC_month_aggr_ECfilt$v.wind_mean) *360/2/pi) +180

# gap-filled:
USOWC_month_aggr_ECfilt$EC_WD_AVG_F <- (atan2(USOWC_month_aggr_ECfilt$u.wind.F_mean, USOWC_month_aggr_ECfilt$v.wind.F_mean) *360/2/pi) +180

# save as .csv

write.csv(USOWC_month_aggr_ECfilt, "path/US_OWC_chec_monthaggr.csv", row.names = FALSE)

# ANNUAL AGGREGATION

USOWC_yr_aggr_ECfilt <- USOWC_ECfilt %>% 
  group_by(year(USOWC_ECfilt$TIMESTAMP_START)) %>%
  summarise_at(c("ch_FCH4_nmolCH4m2s1", 
                 "EC_FCH4", "EC_FCH4_F", "EC_FCH4_F_ANNOPTLM", "EC_GPP_NT", "EC_GPP_DT", "EC_RECO_NT", "EC_RECO_DT",
                 "EC_USTAR", "EC_NEE_F", "EC_H_F", "EC_LE_F", "EC_SW_IN_F", "EC_SW_OUT_F", "EC_LW_IN_F", 
                 "EC_LW_OUT_F", "EC_NETRAD_F", "EC_PPFD_IN_F", "EC_VPD_F", "EC_RH_F",
                 "EC_PA_F", "EC_TA_F", "EC_WTD_F", "EC_P_F", "EC_WS_F", 
                 "EC_LE_F_ANNOPTLM", "EC_TS_1", "EC_NEE_F_ANNOPTLM",
                 "u.wind", "v.wind", "u.wind.F", "v.wind.F"), 
               list(mean = mean,
                    median = median), 
               na.rm = TRUE)

USOWC_yr_aggr_ECfilt <- as.data.frame(USOWC_yr_aggr_ECfilt)

# rename the date column

colnames(USOWC_yr_aggr_ECfilt)[1] <- "Year"

# add the missing extra columns

USOWC_yr_aggr_ECfilt$ch_method <- "manual"
USOWC_yr_aggr_ECfilt$SITE <- "US-OWC"

USOWC_yr_aggr_ECfilt <- USOWC_yr_aggr_ECfilt %>% relocate(SITE, .before = Year)

# add dom_veg
USOWC_yr_aggr_ECfilt$DOM_VEG <- "AERENCHYMATOUS"

# add veg classes
USOWC_yr_aggr_ECfilt$MOSS_BROWN <- 0
USOWC_yr_aggr_ECfilt$MOSS_SPHAGNUM <- 0
USOWC_yr_aggr_ECfilt$AERENCHYMATOUS <- 1
USOWC_yr_aggr_ECfilt$ERI_SHRUB <- 0
USOWC_yr_aggr_ECfilt$TREE <- 0

# add climate
USOWC_yr_aggr_ECfilt$KOPPEN <- "Dfb"

# add ecosystem type
USOWC_yr_aggr_ECfilt$SITE_CLASSIFICATION <- "marsh"

# calculate wind direction average
USOWC_yr_aggr_ECfilt$EC_WD_AVG <- (atan2(USOWC_yr_aggr_ECfilt$u.wind_mean, USOWC_yr_aggr_ECfilt$v.wind_mean) *360/2/pi) +180

# gap-filled:
USOWC_yr_aggr_ECfilt$EC_WD_AVG_F <- (atan2(USOWC_yr_aggr_ECfilt$u.wind.F_mean, USOWC_yr_aggr_ECfilt$v.wind.F_mean) *360/2/pi) +180

# save as .csv

write.csv(USOWC_yr_aggr_ECfilt, "path/US_OWC_chec_yraggr.csv", row.names = FALSE)


### US-UAF ###

# read csv

# EC
USUAF_EC_d <- read.csv("path/US_UAF_EC_subset.csv")

# convert to datetime
USUAF_EC_d$EC_TIMESTAMP_START <- as_datetime(USUAF_EC_d$EC_TIMESTAMP_START, tz = "UTC")
USUAF_EC_d$EC_TIMESTAMP_END <- as_datetime(USUAF_EC_d$EC_TIMESTAMP_END, tz = "UTC")

# CHAMBER
USUAF_ch_d <- read.csv("path/US_UAF.csv")

# convert to datetime
USUAF_ch_d$ch_Datetime <- as_datetime(USUAF_ch_d$ch_Datetime, tz = "UTC")

# combine using merge
Sys.setenv(TZ = "UTC")
USUAF_all <- merge(x = USUAF_ch_d, y = USUAF_EC_d, by.x = "ch_Datetime", by.y =  "EC_TIMESTAMP_START",
                   all.y = TRUE)

# move EC_TIMESTAMP_END and SITE.x
USUAF_all <- USUAF_all %>% relocate(EC_TIMESTAMP_END, .before = SITE.x)
USUAF_all <- USUAF_all %>% relocate(SITE.x, .before = ch_Datetime)

# rename SITE.x
colnames(USUAF_all)[1] <- "SITE"

# rename ch_Datetime to ch_EC_TIMESTAMP_START
colnames(USUAF_all)[2] <- "ch_EC_TIMESTAMP_START"

# remove SITE.y
USUAF_all <- select(USUAF_all, -15)

# chamber FCH4 data is already in nmol m-2 s-1 

write.csv(USUAF_all, "path/US_UAF_combined.csv", row.names = FALSE)

# read the file again
USUAF <- read.csv("path/US_UAF_combined.csv")

# SET EC GAP-FILLED DATA AND CHAMBER DATA TO NA WHEN TIME GAP LONGER THAN 2 MONTHS

USUAF_new <- USUAF %>% mutate(across(c(ch_FCH4_nmolCH4m2s1, ch_TS_2, 
                                       ch_TS_10, ch_TS_20, ch_TS_30, ch_TS_40, ch_VWC,
                                       EC_FCH4_F, EC_NEE_F, EC_H_F, EC_LE_F, EC_SW_IN_F, EC_SW_OUT_F,
                                       EC_LW_IN_F, EC_LW_OUT_F, EC_NETRAD_F, EC_PPFD_IN_F, EC_PPFD_OUT_F,
                                       EC_VPD_F, EC_RH_F, EC_PA_F, EC_TA_F, EC_P_F,EC_WTD_F,EC_G_F,
                                       EC_WS_F, EC_LE_F_ANNOPTLM, EC_NEE_F_ANNOPTLM, EC_FCH4_F_ANNOPTLM, 
                                       EC_FCH4_F_RANDUNC, EC_FCH4_F_ANNOPTLM_UNC), ~ ifelse(EC_FCH4_F_ANNOPTLM_QC == 3, NA, .)))

# convert datetimes to datetime format
USUAF_new$ch_EC_TIMESTAMP_START <- as_datetime(USUAF_new$ch_EC_TIMESTAMP_START)
USUAF_new$ch_EC_TIMESTAMP_END <- as_datetime(USUAF_new$EC_TIMESTAMP_END)

# move EC_FCH4_F_ANNOPTLM next to other flux columns
USUAF_new <- USUAF_new %>% relocate(c(EC_FCH4, EC_FCH4_F, EC_FCH4_F_ANNOPTLM), .after = ch_FCH4_nmolCH4m2s1)

# remove the accidental tiemstamp at the end
USUAF_new2 <- subset(USUAF_new, select = -c(ch_EC_TIMESTAMP_END))

# subset to full days
Sys.setenv(TZ = "UTC")
USUAF_new2$EC_TIMESTAMP_END <- as_datetime(USUAF_new2$EC_TIMESTAMP_END)
USUAF_new2$ch_EC_TIMESTAMP_START <- as_datetime(USUAF_new2$ch_EC_TIMESTAMP_START)

USUAF_new3 <- subset(USUAF_new2, EC_TIMESTAMP_END <= "2018-11-01 00:00:00")

# set EC data to NA when chamber data is NA
USUAF_ECfilt <- USUAF_new3 %>% mutate(across(c(EC_FCH4, EC_FCH4_F, EC_FCH4_F_ANNOPTLM), 
                                                                          ~ ifelse(ch_FCH4_nmolCH4m2s1 == "NA", NA, .)))
# HALF-HOURLY AGGREGATION

USUAF_ECfilt$ch_VWC <- as.numeric(USUAF_ECfilt$ch_VWC)

USUAF_hh_aggr <- USUAF_ECfilt %>% 
  group_by(across(-c(ch_FCH4_nmolCH4m2s1, ch_ID, ch_veg_cover,
                     ch_TS_2, ch_TS_10, ch_TS_20, ch_TS_30, ch_TS_40, ch_VWC,
                     ch_Collar_type))) %>%
  summarise_at(c("ch_FCH4_nmolCH4m2s1", "ch_TS_2", "ch_TS_10", "ch_TS_20", "ch_TS_30", "ch_TS_40",
                 "ch_VWC", "ECCH_diff"), 
               list(mean = mean,
                    median = median), 
               na.rm = TRUE)

USUAF_hh_aggr <- as.data.frame(USUAF_hh_aggr)

# remove rows where EC_FCH4_F_ANNOPTLM = NA

USUAF_hh_aggr <- USUAF_hh_aggr %>% drop_na(EC_FCH4_F_ANNOPTLM)

# move the new columns before EC_FCH4
USUAF_hh_aggr <- USUAF_hh_aggr %>% relocate(ch_FCH4_nmolCH4m2s1_mean, ch_FCH4_nmolCH4m2s1_median, .before = EC_FCH4)

# save as .csv
write.csv(USUAF_hh_aggr, "path/US_UAF_combined_ch_datetimeaggr.csv", row.names = FALSE)


# HOURLY AGGREGATION

# create new date_hour column for grouping
USUAF_ECfilt$ch_EC_TIMESTAMP_START <- as_datetime(USUAF_ECfilt$ch_EC_TIMESTAMP_START)

USUAF_ECfilt <- USUAF_ECfilt %>%
  mutate(date_hour = format(ch_EC_TIMESTAMP_START, "%Y-%m-%d %H"))

# set duplicated values to NA in EC columns
# create vector with EC column names
#rename EC_TIMESTAMP_START and END
colnames(USUAF_ECfilt)[c(2,3)] <- c("TIMESTAMP_START", "TIMESTAMP_END")

cols <- USUAF_ECfilt %>% dplyr::select(starts_with("EC_")) %>% colnames()

# convert duplicates to NA in EC columns
USUAF_ECfilt[cols] <- sapply(USUAF_ECfilt[cols], function(x) 
  ave(x, as_datetime(USUAF_ECfilt$TIMESTAMP_START), FUN = function(x) replace(x, duplicated(x), NA)))

# Calculate the u and v wind components
USUAF_ECfilt$u.wind <- -USUAF_ECfilt$EC_WS *sin(2*pi *USUAF_ECfilt$EC_WD/360)
USUAF_ECfilt$v.wind <- -USUAF_ECfilt$EC_WS *cos(2*pi *USUAF_ECfilt$EC_WD/360)
# # gap-filled (no gap-filled WD)
USUAF_ECfilt$u.wind.F <- -USUAF_ECfilt$EC_WS_F *sin(2*pi *USUAF_ECfilt$EC_WD/360)
USUAF_ECfilt$v.wind.F <- -USUAF_ECfilt$EC_WS_F *cos(2*pi *USUAF_ECfilt$EC_WD/360)

# convert character and logi to numeric
cols.num <- USUAF_ECfilt %>% dplyr::select(starts_with("EC_")) %>% dplyr::select(-EC_FCH4_MEASTYPE) %>% colnames()
USUAF_ECfilt[cols.num] <- sapply(USUAF_ECfilt[cols.num],as.numeric)
USUAF_ECfilt$ch_VWC <- as.numeric(USUAF_ECfilt$ch_VWC)

# aggregate
USUAF_hr_aggr_ECfilt <- USUAF_ECfilt %>% 
  group_by(date_hour) %>%
  summarise_at(c("ch_FCH4_nmolCH4m2s1", "ch_TS_2", "ch_TS_10", "ch_TS_20", "ch_TS_30",
                 "ch_TS_40", "ch_VWC",
                 "EC_GPP_NT", "EC_RECO_NT", "EC_GPP_DT", "EC_RECO_DT",
                 "EC_TS_1", "EC_TS_2",  "EC_TS_3",
                 "EC_TS_4", "EC_TS_5", "EC_TS_6", "EC_TS_7",  "EC_TS_8", "EC_TS_9",
                 "EC_SWC_1", "EC_SWC_2", "EC_SWC_3",
                 "EC_FCH4", "EC_FCH4_F", "EC_FCH4_F_ANNOPTLM",
                 "EC_USTAR", "EC_NEE_F", "EC_H_F", "EC_LE_F", "EC_SW_IN_F", "EC_SW_OUT_F", "EC_LW_IN_F", 
                 "EC_LW_OUT_F", "EC_NETRAD_F", "EC_PPFD_IN_F", "EC_PPFD_OUT_F", "EC_VPD_F", "EC_RH_F", 
                 "EC_PA_F", "EC_TA_F", "EC_G_F", "EC_WTD_F", "EC_P_F", "EC_WS_F", 
                 "EC_LE_F_ANNOPTLM", "EC_NEE_F_ANNOPTLM", "EC_FCH4_F_ANNOPTLM_UNC",
                 "u.wind", "v.wind", "u.wind.F", "v.wind.F"), 
               list(mean = mean,
                    median = median), 
               na.rm = TRUE)

USUAF_hr_aggr_ECfilt <- as.data.frame(USUAF_hr_aggr_ECfilt)

# convert to datetime format
USUAF_hr_aggr_ECfilt$date_hour <- as_datetime(ymd_h(USUAF_hr_aggr_ECfilt$date_hour))

# add the missing columns

USUAF_hr_aggr_ECfilt$Year <- year(USUAF_hr_aggr_ECfilt$date_hour)
USUAF_hr_aggr_ECfilt$Month <- month(USUAF_hr_aggr_ECfilt$date_hour)
USUAF_hr_aggr_ECfilt$Day <- day(USUAF_hr_aggr_ECfilt$date_hour)
USUAF_hr_aggr_ECfilt$Hour <- hour(USUAF_hr_aggr_ECfilt$date_hour)
USUAF_hr_aggr_ECfilt$DOY <- yday(USUAF_hr_aggr_ECfilt$date_hour)

USUAF_hr_aggr_ECfilt$ch_method <- "auto"
USUAF_hr_aggr_ECfilt$SITE <- "US-UAF"

USUAF_hr_aggr_ECfilt <- USUAF_hr_aggr_ECfilt %>% relocate(SITE, .before = date_hour)

# rename date_hour
colnames(USUAF_hr_aggr_ECfilt)[2] <- "TIMESTAMP"

# add dom_veg
USUAF_hr_aggr_ECfilt$DOM_VEG <- "MOSS_SPHAGNUM"

# add veg classes
USUAF_hr_aggr_ECfilt$MOSS_BROWN <- 1
USUAF_hr_aggr_ECfilt$MOSS_SPHAGNUM <- 1
USUAF_hr_aggr_ECfilt$AERENCHYMATOUS <- 1
USUAF_hr_aggr_ECfilt$ERI_SHRUB <- 1
USUAF_hr_aggr_ECfilt$TREE <- 1

# add climate
USUAF_hr_aggr_ECfilt$KOPPEN <- "Dwc"

# add ecosystem type
USUAF_hr_aggr_ECfilt$SITE_CLASSIFICATION <- "bog"

# calculate wind direction average
USUAF_hr_aggr_ECfilt$EC_WD_AVG <- (atan2(USUAF_hr_aggr_ECfilt$u.wind_mean, USUAF_hr_aggr_ECfilt$v.wind_mean) *360/2/pi) +180

# gap-filled:
USUAF_hr_aggr_ECfilt$EC_WD_AVG_F <- (atan2(USUAF_hr_aggr_ECfilt$u.wind.F_mean, USUAF_hr_aggr_ECfilt$v.wind.F_mean) *360/2/pi) +180

# save as .csv

write.csv(USUAF_hr_aggr_ECfilt, "path/US_UAF_combined_filtered_with_outliers_hraggr_ECfilt.csv", row.names = FALSE)

# DAILY AGGREGATION

USUAF_d_aggr_ECfilt <- USUAF_ECfilt %>% 
  group_by(date(USUAF_ECfilt$TIMESTAMP_START)) %>%
  summarise_at(c("ch_FCH4_nmolCH4m2s1", "ch_TS_2", "ch_TS_10", "ch_TS_20", "ch_TS_30",
                 "ch_TS_40", "ch_VWC", 
                 "EC_GPP_NT", "EC_RECO_NT", "EC_GPP_DT", "EC_RECO_DT",
                 "EC_TS_1", "EC_TS_2",  "EC_TS_3",
                 "EC_TS_4", "EC_TS_5", "EC_TS_6", "EC_TS_7",  "EC_TS_8", "EC_TS_9",
                 "EC_SWC_1", "EC_SWC_2", "EC_SWC_3",
                 "EC_FCH4", "EC_FCH4_F", "EC_FCH4_F_ANNOPTLM",
                 "EC_USTAR", "EC_NEE_F", "EC_H_F", "EC_LE_F", "EC_SW_IN_F", "EC_SW_OUT_F", "EC_LW_IN_F", 
                 "EC_LW_OUT_F", "EC_NETRAD_F", "EC_PPFD_IN_F", "EC_PPFD_OUT_F", "EC_VPD_F", "EC_RH_F", 
                 "EC_PA_F", "EC_TA_F", "EC_G_F", "EC_WTD_F", "EC_P_F", "EC_WS_F", 
                 "EC_LE_F_ANNOPTLM", "EC_NEE_F_ANNOPTLM", "EC_FCH4_F_ANNOPTLM_UNC",
                 "u.wind", "v.wind", "u.wind.F", "v.wind.F"), 
               list(mean = mean,
                    median = median), 
               na.rm = TRUE)

USUAF_d_aggr_ECfilt <- as.data.frame(USUAF_d_aggr_ECfilt)

# remove rows where CH4 ANN = NA

USUAF_d_aggr_ECfilt <- USUAF_d_aggr_ECfilt %>% drop_na(EC_FCH4_F_ANNOPTLM_mean)

# rename the date column

colnames(USUAF_d_aggr_ECfilt)[1] <- "DATE"

# convert to date format
USUAF_d_aggr_ECfilt$DATE <- as_date(USUAF_d_aggr_ECfilt$DATE)

# calculate wind direction average
USUAF_d_aggr_ECfilt$EC_WD_AVG <- (atan2(USUAF_d_aggr_ECfilt$u.wind_mean, USUAF_d_aggr_ECfilt$v.wind_mean) *360/2/pi) +180

# gap-filled:
USUAF_d_aggr_ECfilt$EC_WD_AVG_F <- (atan2(USUAF_d_aggr_ECfilt$u.wind.F_mean, USUAF_d_aggr_ECfilt$v.wind.F_mean) *360/2/pi) +180

# add the missing columns

USUAF_d_aggr_ECfilt$Year <- year(USUAF_d_aggr_ECfilt$DATE)
USUAF_d_aggr_ECfilt$Month <- month(USUAF_d_aggr_ECfilt$DATE)
USUAF_d_aggr_ECfilt$Day <- day(USUAF_d_aggr_ECfilt$DATE)
USUAF_d_aggr_ECfilt$DOY <- yday(USUAF_d_aggr_ECfilt$DATE)

USUAF_d_aggr_ECfilt$ch_method <- "auto"
USUAF_d_aggr_ECfilt$SITE <- "US-UAF"

USUAF_d_aggr_ECfilt <- USUAF_d_aggr_ECfilt %>% relocate(SITE, .before = DATE)

# add dom_veg
USUAF_d_aggr_ECfilt$DOM_VEG <- "MOSS_SPHAGNUM"

# add veg classes
USUAF_d_aggr_ECfilt$MOSS_BROWN <- 1
USUAF_d_aggr_ECfilt$MOSS_SPHAGNUM <- 1
USUAF_d_aggr_ECfilt$AERENCHYMATOUS <- 1
USUAF_d_aggr_ECfilt$ERI_SHRUB <- 1
USUAF_d_aggr_ECfilt$TREE <- 1

# add climate
USUAF_d_aggr_ECfilt$KOPPEN <- "Dwc"

# add ecosystem type
USUAF_d_aggr_ECfilt$SITE_CLASSIFICATION <- "bog"

# save as .csv

write.csv(USUAF_d_aggr_ECfilt, "path/US_UAF_chec_daggr.csv", row.names = FALSE)


# WEEKLY AGGREGATION

USUAF_week_aggr_ECfilt <- USUAF_ECfilt %>% 
  group_by(year(USUAF_ECfilt$TIMESTAMP_START), isoweek(USUAF_ECfilt$TIMESTAMP_START)) %>%
  summarise_at(c("ch_FCH4_nmolCH4m2s1", "ch_TS_2", "ch_TS_10", "ch_TS_20", "ch_TS_30",
                 "ch_TS_40", "ch_VWC", 
                 "EC_GPP_NT", "EC_RECO_NT", "EC_GPP_DT", "EC_RECO_DT",
                 "EC_TS_1", "EC_TS_2",  "EC_TS_3",
                 "EC_TS_4", "EC_TS_5", "EC_TS_6", "EC_TS_7",  "EC_TS_8", "EC_TS_9",
                 "EC_SWC_1", "EC_SWC_2", "EC_SWC_3",
                 "EC_FCH4", "EC_FCH4_F", "EC_FCH4_F_ANNOPTLM",
                 "EC_USTAR", "EC_NEE_F", "EC_H_F", "EC_LE_F", "EC_SW_IN_F", "EC_SW_OUT_F", "EC_LW_IN_F", 
                 "EC_LW_OUT_F", "EC_NETRAD_F", "EC_PPFD_IN_F", "EC_PPFD_OUT_F", "EC_VPD_F", "EC_RH_F", 
                 "EC_PA_F", "EC_TA_F", "EC_G_F", "EC_WTD_F", "EC_P_F", "EC_WS_F", 
                 "EC_LE_F_ANNOPTLM", "EC_NEE_F_ANNOPTLM", "EC_FCH4_F_ANNOPTLM_UNC",
                 "u.wind", "v.wind", "u.wind.F", "v.wind.F"), 
               list(mean = mean,
                    median = median), 
               na.rm = TRUE)

USUAF_week_aggr_ECfilt <- as.data.frame(USUAF_week_aggr_ECfilt)

# rename the date column

colnames(USUAF_week_aggr_ECfilt)[1] <- "Year"
colnames(USUAF_week_aggr_ECfilt)[2] <- "Week_of_year"

# add the missing extra columns

USUAF_week_aggr_ECfilt$ch_method <- "auto"
USUAF_week_aggr_ECfilt$SITE <- "US-UAF"

USUAF_week_aggr_ECfilt <- USUAF_week_aggr_ECfilt %>% relocate(SITE, .before = Year)

# add dom_veg
USUAF_week_aggr_ECfilt$DOM_VEG <- "MOSS_SPHAGNUM"

# add veg classes
USUAF_week_aggr_ECfilt$MOSS_BROWN <- 1
USUAF_week_aggr_ECfilt$MOSS_SPHAGNUM <- 1
USUAF_week_aggr_ECfilt$AERENCHYMATOUS <- 1
USUAF_week_aggr_ECfilt$ERI_SHRUB <- 1
USUAF_week_aggr_ECfilt$TREE <- 1

# add climate
USUAF_week_aggr_ECfilt$KOPPEN <- "Dwc"

# add ecosystem type
USUAF_week_aggr_ECfilt$SITE_CLASSIFICATION <- "bog"

# calculate wind direction average
USUAF_week_aggr_ECfilt$EC_WD_AVG <- (atan2(USUAF_week_aggr_ECfilt$u.wind_mean, USUAF_week_aggr_ECfilt$v.wind_mean) *360/2/pi) +180

# gap-filled:
USUAF_week_aggr_ECfilt$EC_WD_AVG_F <- (atan2(USUAF_week_aggr_ECfilt$u.wind.F_mean, USUAF_week_aggr_ECfilt$v.wind.F_mean) *360/2/pi) +180

# add month

# Extract the correct month directly using lubridate::make_date and ISO week (1st day of the week is Monday)
USUAF_week_aggr_ECfilt$Month <- month(make_date(USUAF_week_aggr_ECfilt$Year) + weeks(USUAF_week_aggr_ECfilt$Week_of_year))
USUAF_week_aggr_ECfilt <- USUAF_week_aggr_ECfilt %>% relocate(Month, .after = Week_of_year)
# weeks and months are correct

# save as .csv

write.csv(USUAF_week_aggr_ECfilt, "path/US_UAF_chec_weekaggr.csv", row.names = FALSE)


# MONHTLY AGGREGATION

USUAF_month_aggr_ECfilt <- USUAF_ECfilt %>% 
  group_by(year(USUAF_ECfilt$TIMESTAMP_START), month(USUAF_ECfilt$TIMESTAMP_START)) %>%
  summarise_at(c("ch_FCH4_nmolCH4m2s1", "ch_TS_2", "ch_TS_10", "ch_TS_20", "ch_TS_30",
                 "ch_TS_40", "ch_VWC", 
                 "EC_GPP_NT", "EC_RECO_NT", "EC_GPP_DT", "EC_RECO_DT",
                 "EC_TS_1", "EC_TS_2",  "EC_TS_3",
                 "EC_TS_4", "EC_TS_5", "EC_TS_6", "EC_TS_7",  "EC_TS_8", "EC_TS_9",
                 "EC_SWC_1", "EC_SWC_2", "EC_SWC_3",
                 "EC_FCH4", "EC_FCH4_F", "EC_FCH4_F_ANNOPTLM",
                 "EC_USTAR", "EC_NEE_F", "EC_H_F", "EC_LE_F", "EC_SW_IN_F", "EC_SW_OUT_F", "EC_LW_IN_F", 
                 "EC_LW_OUT_F", "EC_NETRAD_F", "EC_PPFD_IN_F", "EC_PPFD_OUT_F", "EC_VPD_F", "EC_RH_F", 
                 "EC_PA_F", "EC_TA_F", "EC_G_F", "EC_WTD_F", "EC_P_F", "EC_WS_F", 
                 "EC_LE_F_ANNOPTLM", "EC_NEE_F_ANNOPTLM", "EC_FCH4_F_ANNOPTLM_UNC",
                 "u.wind", "v.wind", "u.wind.F", "v.wind.F"), 
               list(mean = mean,
                    median = median), 
               na.rm = TRUE)

USUAF_month_aggr_ECfilt <- as.data.frame(USUAF_month_aggr_ECfilt)

# rename the date column

colnames(USUAF_month_aggr_ECfilt)[1] <- "Year"
colnames(USUAF_month_aggr_ECfilt)[2] <- "Month"

# add the missing extra columns

USUAF_month_aggr_ECfilt$ch_method <- "auto"
USUAF_month_aggr_ECfilt$SITE <- "US-UAF"

USUAF_month_aggr_ECfilt <- USUAF_month_aggr_ECfilt %>% relocate(SITE, .before = Year)

# add dom_veg
USUAF_month_aggr_ECfilt$DOM_VEG <- "MOSS_SPHAGNUM"

# add veg classes
USUAF_month_aggr_ECfilt$MOSS_BROWN <- 1
USUAF_month_aggr_ECfilt$MOSS_SPHAGNUM <- 1
USUAF_month_aggr_ECfilt$AERENCHYMATOUS <- 1
USUAF_month_aggr_ECfilt$ERI_SHRUB <- 1
USUAF_month_aggr_ECfilt$TREE <- 1

# add climate
USUAF_month_aggr_ECfilt$KOPPEN <- "Dwc"

# add ecosystem type
USUAF_month_aggr_ECfilt$SITE_CLASSIFICATION <- "bog"

# calculate wind direction average
USUAF_month_aggr_ECfilt$EC_WD_AVG <- (atan2(USUAF_month_aggr_ECfilt$u.wind_mean, USUAF_month_aggr_ECfilt$v.wind_mean) *360/2/pi) +180

# gap-filled:
USUAF_month_aggr_ECfilt$EC_WD_AVG_F <- (atan2(USUAF_month_aggr_ECfilt$u.wind.F_mean, USUAF_month_aggr_ECfilt$v.wind.F_mean) *360/2/pi) +180

# save as .csv

write.csv(USUAF_month_aggr_ECfilt, "path/US_UAF_chec_monthaggr_ECfilt.csv", row.names = FALSE)


# ANNUAL AGGREGATION

USUAF_yr_aggr_ECfilt <- USUAF_ECfilt %>% 
  group_by(year(USUAF_ECfilt$TIMESTAMP_START)) %>%
  summarise_at(c("ch_FCH4_nmolCH4m2s1", "ch_TS_2", "ch_TS_10", "ch_TS_20", "ch_TS_30",
                 "ch_TS_40", "ch_VWC", 
                 "EC_GPP_NT", "EC_RECO_NT", "EC_GPP_DT", "EC_RECO_DT",
                 "EC_TS_1", "EC_TS_2",  "EC_TS_3",
                 "EC_TS_4", "EC_TS_5", "EC_TS_6", "EC_TS_7",  "EC_TS_8", "EC_TS_9",
                 "EC_SWC_1", "EC_SWC_2", "EC_SWC_3",
                 "EC_FCH4", "EC_FCH4_F", "EC_FCH4_F_ANNOPTLM",
                 "EC_USTAR", "EC_NEE_F", "EC_H_F", "EC_LE_F", "EC_SW_IN_F", "EC_SW_OUT_F", "EC_LW_IN_F", 
                 "EC_LW_OUT_F", "EC_NETRAD_F", "EC_PPFD_IN_F", "EC_PPFD_OUT_F", "EC_VPD_F", "EC_RH_F", 
                 "EC_PA_F", "EC_TA_F", "EC_G_F", "EC_WTD_F", "EC_P_F", "EC_WS_F", 
                 "EC_LE_F_ANNOPTLM", "EC_NEE_F_ANNOPTLM", "EC_FCH4_F_ANNOPTLM_UNC",
                 "u.wind", "v.wind", "u.wind.F", "v.wind.F"), 
               list(mean = mean,
                    median = median), 
               na.rm = TRUE)

USUAF_yr_aggr_ECfilt <- as.data.frame(USUAF_yr_aggr_ECfilt)

# rename the date column

colnames(USUAF_yr_aggr_ECfilt)[1] <- "Year"

# add the missing extra columns

USUAF_yr_aggr_ECfilt$ch_method <- "auto"
USUAF_yr_aggr_ECfilt$SITE <- "US-UAF"

USUAF_yr_aggr_ECfilt <- USUAF_yr_aggr_ECfilt %>% relocate(SITE, .before = Year)

# add dom_veg
USUAF_yr_aggr_ECfilt$DOM_VEG <- "MOSS_SPHAGNUM"

# add veg classes
USUAF_yr_aggr_ECfilt$MOSS_BROWN <- 1
USUAF_yr_aggr_ECfilt$MOSS_SPHAGNUM <- 1
USUAF_yr_aggr_ECfilt$AERENCHYMATOUS <- 1
USUAF_yr_aggr_ECfilt$ERI_SHRUB <- 1
USUAF_yr_aggr_ECfilt$TREE <- 1

# add climate
USUAF_yr_aggr_ECfilt$KOPPEN <- "Dwc"

# add ecosystem type
USUAF_yr_aggr_ECfilt$SITE_CLASSIFICATION <- "bog"

# calculate wind direction average
USUAF_yr_aggr_ECfilt$EC_WD_AVG <- (atan2(USUAF_yr_aggr_ECfilt$u.wind_mean, USUAF_yr_aggr_ECfilt$v.wind_mean) *360/2/pi) +180

# gap-filled:
USUAF_yr_aggr_ECfilt$EC_WD_AVG_F <- (atan2(USUAF_yr_aggr_ECfilt$u.wind.F_mean, USUAF_yr_aggr_ECfilt$v.wind.F_mean) *360/2/pi) +180

# save as .csv

write.csv(USUAF_yr_aggr_ECfilt, "path/US_UAF_chec_yraggr.csv", row.names = FALSE)



### US-STJ ###

# read csv

# EC
USSTJ_EC_d <- read.csv("path/US_STJ_EC_subset.csv")

# convert to datetime
USSTJ_EC_d$EC_TIMESTAMP_START <- as_datetime(USSTJ_EC_d$EC_TIMESTAMP_START)
USSTJ_EC_d$EC_TIMESTAMP_END <- as_datetime(USSTJ_EC_d$EC_TIMESTAMP_END)

# CHAMBER
USSTJ_ch_d <- read.csv("path/US-STJ.csv")

colnames(USSTJ_ch_d)[2] <- "ch_Datetime"

USSTJ_ch_d <- subset(USSTJ_ch_d, select = -c(Date, Time))

USSTJ_ch_d$ch_Datetime <- as_datetime(USSTJ_ch_d$ch_Datetime)

# combine the data frames using SQL:

# set local environment time zone to UTC so the time stamps will match
Sys.setenv(TZ = "UTC")

USSTJ_all <- sqldf('SELECT ch_Datetime, EC_TIMESTAMP_START, EC_TIMESTAMP_END, ch_ID, ch_FCH4_nmolCH4m2s1, EC_ch4_flux, 
EC_PA, EC_TA, EC_VPD, EC_WS, EC_WD, EC_ch4_scf, EC_Level_YSI, 
EC_Level_YSI_f, EC_Level_NOAA, EC_TS, EC_USTAR, EC_NEE_orig, EC_NEE_f, EC_CH4_orig, EC_CH4_f, EC_CH4_f_RF, EC_FCH4_RF_filled,
ch_method
      FROM USSTJ_EC_d 
      LEFT JOIN USSTJ_ch_d ON ch_Datetime BETWEEN EC_TIMESTAMP_START and EC_TIMESTAMP_END')

# remove NA rows
USSTJ_all <- USSTJ_all %>% drop_na(ch_Datetime)

# add Site column
USSTJ_all$SITE <- "US-STJ"

# move to the front
USSTJ_all <- USSTJ_all %>%
  select(SITE, everything())

# remove some extra columns
USSTJ_all <- subset(USSTJ_all, select = -c(EC_ch4_flux, EC_Level_YSI, EC_Level_NOAA, EC_ch4_scf))

# rename columns
USSTJ_all <- USSTJ_all %>%
  rename(
    EC_FCH4_nmolCH4m2s1   = EC_CH4_orig,
    EC_FCH4_F_MDS     = EC_CH4_f,
    EC_FCH4_F_RF = EC_CH4_f_RF,
    EC_FCH4_F_RF_filled = EC_FCH4_RF_filled, 
    EC_WTD_F = EC_Level_YSI_f,
    EC_NEE_F_MDS = EC_NEE_f,
    EC_NEE = EC_NEE_orig
  )

USSTJ_all <- subset(USSTJ_all, select = -EC_ch4_scf)

write.csv(USSTJ_all, "path/US_STJ_combined.csv", row.names = FALSE)

# read in the file again
USSTJ <- read.csv("path/US_STJ_combined.csv")

# remove one column
USSTJ <- subset(USSTJ, select = -EC_FCH4_F_RF)

colnames(USSTJ)[18] <- "EC_FCH4_F_RF"

USSTJ <- USSTJ %>% drop_na(EC_FCH4_nmolCH4m2s1)

# DAILY AGGREGATION (also corresponds to chambers aggregated to EC timestamp)

# convert date to datetime

USSTJ$EC_TIMESTAMP_START <- as_datetime(USSTJ$EC_TIMESTAMP_START)

# set duplicated values to NA in EC columns
# create vector with EC column names
#rename EC_TIMESTAMP_START and END
colnames(USSTJ)[c(3,4)] <- c("TIMESTAMP_START", "TIMESTAMP_END")

cols <- USSTJ %>% dplyr::select(starts_with("EC_")) %>% colnames()

# convert duplicates to NA in EC columns
USSTJ[cols] <- sapply(USSTJ[cols], function(x) 
  ave(x, as_datetime(USSTJ$TIMESTAMP_START), FUN = function(x) replace(x, duplicated(x), NA)))

# convert character and logi to numeric
cols.num <- USSTJ %>% dplyr::select(starts_with(c("EC_", "ch_FCH4"))) %>% colnames()
USSTJ[cols.num] <- sapply(USSTJ[cols.num],as.numeric)

# Calculate the u and v wind components
USSTJ$u.wind <- -USSTJ$EC_WS *sin(2*pi *USSTJ$EC_WD/360)
USSTJ$v.wind <- -USSTJ$EC_WS *cos(2*pi *USSTJ$EC_WD/360)

USSTJ <- subset(USSTJ, select = -c(EC_FCH4_F_RF, EC_FCH4_F_RF_filled))

# aggregate
USSTJ_d_aggr <- USSTJ %>% 
  group_by(date(USSTJ$TIMESTAMP_START)) %>%
  summarise_at(c("ch_FCH4_nmolCH4m2s1",
                 "EC_FCH4_nmolCH4m2s1",
                 "EC_USTAR", "EC_NEE", "EC_NEE_F_MDS",
                 "EC_VPD",
                 "EC_PA", "EC_TA", "EC_TS" , "EC_WTD_F", "EC_WS",
                 "u.wind", "v.wind"), 
               list(mean = mean,
                    median = median), 
               na.rm = TRUE)

USSTJ_d_aggr <- as.data.frame(USSTJ_d_aggr)


# rename the date column

colnames(USSTJ_d_aggr)[1] <- "DATE"

# convert to date format
USSTJ_d_aggr$DATE <- as_date(USSTJ_d_aggr$DATE)

# calculate wind direction average
USSTJ_d_aggr$EC_WD_AVG <- (atan2(USSTJ_d_aggr$u.wind_mean, USSTJ_d_aggr$v.wind_mean) *360/2/pi) +180

# add the missing columns

USSTJ_d_aggr$Year <- year(USSTJ_d_aggr$DATE)
USSTJ_d_aggr$Month <- month(USSTJ_d_aggr$DATE)
USSTJ_d_aggr$Day <- day(USSTJ_d_aggr$DATE)
USSTJ_d_aggr$DOY <- yday(USSTJ_d_aggr$DATE)

USSTJ_d_aggr$ch_method <- "manual"
USSTJ_d_aggr$SITE <- "US-STJ"

USSTJ_d_aggr <- USSTJ_d_aggr %>% relocate(SITE, .before = DATE)

# remove rows with FCH4_ANNOTPLTM = NA

USSTJ_d_aggr <- USSTJ_d_aggr %>% drop_na(EC_FCH4_nmolCH4m2s1_mean) 

# add dom_veg
USSTJ_d_aggr$DOM_VEG <- "AERENCHYMATOUS"

# add veg classes
USSTJ_d_aggr$MOSS_BROWN <- 0
USSTJ_d_aggr$MOSS_SPHAGNUM <- 0
USSTJ_d_aggr$AERENCHYMATOUS <- 1
USSTJ_d_aggr$ERI_SHRUB <- 0
USSTJ_d_aggr$TREE <- 0

# add climate
USSTJ_d_aggr$KOPPEN <- "Cfa"

# add ecosystem type
USSTJ_d_aggr$SITE_CLASSIFICATION <- "salt marsh"

# save as .csv

write.csv(USSTJ_d_aggr, "path/US_STJ_chec_daggr.csv", row.names = FALSE)



# WEEKLY AGGREGATION

USSTJ_week_aggr <- USSTJ %>% 
  group_by(year(USSTJ$TIMESTAMP_START), isoweek(USSTJ$TIMESTAMP_START)) %>%
  summarise_at(c("ch_FCH4_nmolCH4m2s1",
                 "EC_FCH4_nmolCH4m2s1",
                 "EC_USTAR", "EC_NEE", "EC_NEE_F_MDS",
                 "EC_VPD",
                 "EC_PA", "EC_TA", "EC_TS","EC_WTD_F", "EC_WS",
                 "u.wind", "v.wind"), 
               list(mean = mean,
                    median = median), 
               na.rm = TRUE)

USSTJ_week_aggr <- as.data.frame(USSTJ_week_aggr)

# rename the date column

colnames(USSTJ_week_aggr)[1] <- "Year"
colnames(USSTJ_week_aggr)[2] <- "Week_of_year"

# add the missing extra columns

USSTJ_week_aggr$ch_method <- "manual"
USSTJ_week_aggr$SITE <- "US-STJ"

USSTJ_week_aggr <- USSTJ_week_aggr %>% relocate(SITE, .before = Year)

# add dom_veg
USSTJ_week_aggr$DOM_VEG <- "AERENCHYMATOUS"

# add veg classes
USSTJ_week_aggr$MOSS_BROWN <- 0
USSTJ_week_aggr$MOSS_SPHAGNUM <- 0
USSTJ_week_aggr$AERENCHYMATOUS <- 1
USSTJ_week_aggr$ERI_SHRUB <- 0
USSTJ_week_aggr$TREE <- 0

# add climate
USSTJ_week_aggr$KOPPEN <- "Cfa"

# add ecosystem type
USSTJ_week_aggr$SITE_CLASSIFICATION <- "salt marsh"

# calculate wind direction average
USSTJ_week_aggr$EC_WD_AVG <- (atan2(USSTJ_week_aggr$u.wind_mean, USSTJ_week_aggr$v.wind_mean) *360/2/pi) +180


# add month
# 
USSTJ_week_aggr$Month <- month(as.Date(paste(USSTJ_week_aggr$Year, USSTJ_week_aggr$Week_of_year, 1, sep = "-"), format = "%Y-%U-%u"))

USSTJ_week_aggr <- USSTJ_week_aggr %>% relocate(Month, .after = Week_of_year)

# fix some months
USSTJ_week_aggr$Month[USSTJ_week_aggr$Week_of_year == 18 & USSTJ_week_aggr$Year == 2020] <- 4  # was 5, should be Apr
USSTJ_week_aggr$Month[USSTJ_week_aggr$Week_of_year == 22 & USSTJ_week_aggr$Year == 2020] <- 5  # was 6, should be May
USSTJ_week_aggr$Month[USSTJ_week_aggr$Week_of_year == 44 & USSTJ_week_aggr$Year == 2020] <- 10 # was 11, should be Oct

# save as .csv

write.csv(USSTJ_week_aggr, "path/US_STJ_chec_weekaggr.csv", row.names = FALSE)


# MONTHLY AGGREGATION

USSTJ_month_aggr <- USSTJ %>% 
  group_by(year(USSTJ$TIMESTAMP_START), month(USSTJ$TIMESTAMP_START)) %>%
  summarise_at(c("ch_FCH4_nmolCH4m2s1", 
                 "EC_FCH4_nmolCH4m2s1",
                 "EC_USTAR", "EC_NEE", "EC_NEE_F_MDS",
                 "EC_VPD",
                 "EC_PA", "EC_TA", "EC_TS", "EC_WTD_F", "EC_WS",
                 "u.wind", "v.wind"), 
               list(mean = mean,
                    median = median), 
               na.rm = TRUE)

USSTJ_month_aggr <- as.data.frame(USSTJ_month_aggr)

# rename the date column

colnames(USSTJ_month_aggr)[1] <- "Year"
colnames(USSTJ_month_aggr)[2] <- "Month"

# add the missing extra columns

USSTJ_month_aggr$ch_method <- "manual"
USSTJ_month_aggr$SITE <- "US-STJ"

USSTJ_month_aggr <- USSTJ_month_aggr %>% relocate(SITE, .before = Year)

# add dom_veg
USSTJ_month_aggr$DOM_VEG <- "AERENCHYMATOUS"

# add veg classes
USSTJ_month_aggr$MOSS_BROWN <- 0
USSTJ_month_aggr$MOSS_SPHAGNUM <- 0
USSTJ_month_aggr$AERENCHYMATOUS <- 1
USSTJ_month_aggr$ERI_SHRUB <- 0
USSTJ_month_aggr$TREE <- 0

# add climate
USSTJ_month_aggr$KOPPEN <- "Cfa"

# add ecosystem type
USSTJ_month_aggr$SITE_CLASSIFICATION <- "salt marsh"

# calculate wind direction average
USSTJ_month_aggr$EC_WD_AVG <- (atan2(USSTJ_month_aggr$u.wind_mean, USSTJ_month_aggr$v.wind_mean) *360/2/pi) +180

# save as .csv

write.csv(USSTJ_month_aggr, "path/US_STJ_chec_monthaggr.csv", row.names = FALSE)


# ANNUAL AGGREGATION

USSTJ_yr_aggr <- USSTJ %>% 
  group_by(year(USSTJ$TIMESTAMP_START)) %>%
  summarise_at(c("ch_FCH4_nmolCH4m2s1", 
                 "EC_FCH4_nmolCH4m2s1",
                 "EC_USTAR", "EC_NEE", "EC_NEE_F_MDS",
                 "EC_VPD",
                 "EC_PA", "EC_TA", "EC_TS","EC_WTD_F", "EC_WS",
                 "u.wind", "v.wind"), 
               list(mean = mean,
                    median = median), 
               na.rm = TRUE)

USSTJ_yr_aggr <- as.data.frame(USSTJ_yr_aggr)

# rename the date column

colnames(USSTJ_yr_aggr)[1] <- "Year"

# add the missing extra columns

USSTJ_yr_aggr$ch_method <- "manual"
USSTJ_yr_aggr$SITE <- "US-STJ"

USSTJ_yr_aggr <- USSTJ_yr_aggr %>% relocate(SITE, .before = Year)

# add dom_veg
USSTJ_yr_aggr$DOM_VEG <- "AERENCHYMATOUS"

# add veg classes
USSTJ_yr_aggr$MOSS_BROWN <- 0
USSTJ_yr_aggr$MOSS_SPHAGNUM <- 0
USSTJ_yr_aggr$AERENCHYMATOUS <- 1
USSTJ_yr_aggr$ERI_SHRUB <- 0
USSTJ_yr_aggr$TREE <- 0

# add climate
USSTJ_yr_aggr$KOPPEN <- "Cfa"

# add ecosystem type
USSTJ_yr_aggr$SITE_CLASSIFICATION <- "salt marsh"

# calculate wind direction average
USSTJ_yr_aggr$EC_WD_AVG <- (atan2(USSTJ_yr_aggr$u.wind_mean, USSTJ_yr_aggr$v.wind_mean) *360/2/pi) +180

# save as .csv

write.csv(USSTJ_yr_aggr, "path/US_STJ_chec_yraggr.csv", row.names = FALSE)


########################################################################################
### COMBINE SITE-SPECIFIC DATASETS TO MULTI-SITE TEMPORAL AGGREGATION DATASETS (N=6) ###
########################################################################################


###### ALL SITES COMBINED DATASET ######

# This script is organized from raw / half-hourly data preparation to coarser aggregations
# and finally the publication-ready datasets. Updated site datasets added later
# (US-STJ, SE-DEG, CN-HGU) are already included in the published input files;
# if you re-run the workflow from raw site files, they can be added analogously
# to the other site-specific sections above.


#### RAW DATA SET, NOTHING AGGREGATED ####
#### ALL INDIVIDUAL CHAMBER MEASUREMENTS INCLUDED ####
#### NOTE: these datasets have not been published because not all raw chamber data sets have been published by data providers ####
#### please contact chamber data providers for these data ####

# read in datasets

CNHGU_ECfilt <- read.csv("path/CN_HGU_combined_ECfilt.csv")

FISI2_ECfilt <- read.csv("path/FI_Si2_combined_ECfilt.csv")

SEDEG_ECfilt <- read.csv("path/SE_DEG_combined_ECfilt.csv")

USHO1_ECfilt <- read.csv("path/US_HO1_combined_ECfilt.csv")

USLA1_ECfilt <- read.csv("path/US_LA1_combined_ECfilt.csv")

USLA2_ECfilt <- read.csv("path/US_LA2_combined_ECfilt.csv")

USLOS_ECfilt <- read.csv("path/US_LOS_combined_ECfilt.csv")

USOWC_ECfilt <- read.csv("path/US_OWC_combined_ECfilt.csv")

USUAF_ECfilt <- read.csv("path/US_UAF_combined_ECfilt.csv")

USSTJ_ECfilt <- read.csv("path/US_STJ_combined.csv")

CNHGU_ECfilt$ECCH_diff <- CNHGU_ECfilt$EC_FCH4_F_ANNOPTLM - CNHGU_ECfilt$ch_FCH4_nmolCH4m2s1

# harmonize the datasets so that column names are the same in all of them

# to make combining easier, rename TIMESTAMP_START and TIMESTAMP_END columns to just TIMESTAMP_START and TIMESTAMP_END

colnames(FISI2_ECfilt)[2] <- "DATE"
colnames(FISI2_ECfilt)[6] <- "ch_WTL"

SEDEG_ECfilt <- SEDEG_ECfilt %>% rename_with(~ c("TIMESTAMP_START", "TIMESTAMP_END"), all_of(c("EC_TIMESTAMP_START", "EC_TIMESTAMP_END")))
CNHGU_ECfilt <- CNHGU_ECfilt %>% rename_with(~ c("TIMESTAMP_START", "TIMESTAMP_END"), all_of(c("EC_TIMESTAMP_START", "EC_TIMESTAMP_END")))

colnames(USLA1_ECfilt)[1] <- "DATE"
colnames(USLA1_ECfilt)[2] <- "SITE"

USLA1_ECfilt <- USLA1_ECfilt %>% relocate(SITE, .before = DATE)

colnames(USLA2_ECfilt)[2] <- "DATE"

USLOS_ECfilt <- USLOS_ECfilt %>% rename_with(~ c("TIMESTAMP_START", "TIMESTAMP_END"), all_of(c("EC_TIMESTAMP_START", "EC_TIMESTAMP_END")))

USHO1_ECfilt <- USHO1_ECfilt %>% rename_with(~ c("TIMESTAMP_START", "TIMESTAMP_END"), all_of(c("ch_EC_TIMESTAMP_START", "ch_EC_TIMESTAMP_END")))

USOWC_ECfilt <- USOWC_ECfilt %>% rename_with(~ c("TIMESTAMP_START", "TIMESTAMP_END"), all_of(c("EC_TIMESTAMP_START", "EC_TIMESTAMP_END")))

USUAF_ECfilt <- USUAF_ECfilt %>% rename_with(~ c("TIMESTAMP_START", "TIMESTAMP_END"), all_of(c("ch_EC_TIMESTAMP_START", "EC_TIMESTAMP_END")))

USSTJ_ECfilt <- USSTJ_ECfilt %>% rename_with(~ c("TIMESTAMP_START", "TIMESTAMP_END"), all_of(c("EC_TIMESTAMP_START", "EC_TIMESTAMP_END")))

# change USHO1 ch umol column name to nmol

# delete the current nmol

USHO1_ECfilt <- subset(USHO1_ECfilt, select = -c(ch_FCH4_nmolCH4m2s1))

USHO1_ECfilt <- USHO1_ECfilt %>% rename_with(~ "ch_FCH4_nmolCH4m2s1", all_of("ch_FCH4_umolCH4m2s1"))

# add ch_ to USUAF ch_method
USUAF_ECfilt <- USUAF_ECfilt %>% rename_with(~ "ch_method", all_of("method"))

# remove ch_Outlierflag
USHO1_ECfilt <- subset(USHO1_ECfilt, select = -c(ch_OutlierFlag))
USLOS_ECfilt <- subset(USLOS_ECfilt, select = -c(ch_OutlierFlag))
USOWC_ECfilt <- subset(USOWC_ECfilt, select = -c(ch_OutlierFlag))
USUAF_ECfilt <- subset(USUAF_ECfilt, select = -c(ch_OutlierFlag))
USLA1_ECfilt <- subset(USLA1_ECfilt, select = -c(ch_OutlierFlag))
FISI2_ECfilt <- subset(FISI2_ECfilt, select = -c(ch_OutlierFlag))
USLA2_ECfilt <- subset(USLA2_ECfilt, select = -c(ch_OutlierFlag))

# drop na from USHO1
USHO1_ECfilt <- USHO1_ECfilt %>% drop_na("ch_FCH4_nmolCH4m2s1")
SEDEG_ECfilt <- SEDEG_ECfilt %>% drop_na("ch_ID")

CNHGU_ECfilt <- CNHGU_ECfilt %>% drop_na("ch_ID")
USLOS_ECfilt <- USLOS_ECfilt %>% drop_na("ch_ID")
USOWC_ECfilt <- USOWC_ECfilt %>% drop_na("ch_FCH4_nmolCH4m2s1")

# remove plot_group from FISI2

FISI2_ECfilt <- subset(FISI2_ECfilt, select = -ch_Plot_group)

# add DATE to all dfs
CNHGU_ECfilt$TIMESTAMP_START <- as_datetime(CNHGU_ECfilt$TIMESTAMP_START)
CNHGU_ECfilt$DATE <- date(CNHGU_ECfilt$TIMESTAMP_START)

USHO1_ECfilt$TIMESTAMP_START <- as_datetime(USHO1_ECfilt$TIMESTAMP_START)
USHO1_ECfilt$DATE <- date(USHO1_ECfilt$TIMESTAMP_START)

USLOS_ECfilt$TIMESTAMP_START <- as_datetime(USLOS_ECfilt$TIMESTAMP_START)
USLOS_ECfilt$DATE <- date(USLOS_ECfilt$TIMESTAMP_START)

SEDEG_ECfilt$TIMESTAMP_START <- as_datetime(SEDEG_ECfilt$TIMESTAMP_START)
SEDEG_ECfilt$DATE <- date(SEDEG_ECfilt$TIMESTAMP_START)

USOWC_ECfilt$TIMESTAMP_START <- as_datetime(USOWC_ECfilt$TIMESTAMP_START)
USOWC_ECfilt$DATE <- date(USOWC_ECfilt$TIMESTAMP_START)

USUAF_ECfilt$TIMESTAMP_START <- as_datetime(USUAF_ECfilt$TIMESTAMP_START)
USUAF_ECfilt$DATE <- date(USUAF_ECfilt$TIMESTAMP_START)

USSTJ_ECfilt$TIMESTAMP_START <- as_datetime(USSTJ_ECfilt$TIMESTAMP_START)
USSTJ_ECfilt$DATE <- date(USSTJ_ECfilt$TIMESTAMP_START)

# drop na from EC FCH4 nmol CH4 m2s1

USSTJ_ECfilt <- USSTJ_ECfilt %>% drop_na(EC_FCH4_nmolCH4m2s1)

# rename US-STJ columns
USSTJ_ECfilt <- USSTJ_ECfilt %>%
  rename(
    EC_NEE_F   = EC_NEE_F_MDS,
    EC_FCH4 = EC_FCH4_nmolCH4m2s1
  )

# remove TA and EC_NEE (gap-filled NEE better)

USSTJ_ECfilt <- subset(USSTJ_ECfilt, select = -c(EC_TA,  EC_NEE))

USSTJ_ECfilt <- USSTJ_ECfilt %>%
  rename(
    EC_WS_F   = EC_WS, 
    EC_VPD_F = EC_VPD
  )

### combine FCH4 values from USSTJ to ANNOPTLM values from other sites
# --> create new column for this

USSTJ_ECfilt$EC_FCH4_comb <- USSTJ_ECfilt$EC_FCH4

CNHGU_ECfilt$EC_FCH4_comb <- CNHGU_ECfilt$EC_FCH4_F_ANNOPTLM
FISI2_ECfilt$EC_FCH4_comb <- FISI2_ECfilt$EC_FCH4_F_ANNOPTLM
SEDEG_ECfilt$EC_FCH4_comb <- SEDEG_ECfilt$EC_FCH4_F_ANNOPTLM
USHO1_ECfilt$EC_FCH4_comb <- USHO1_ECfilt$EC_FCH4_F_ANNOPTLM
USLA1_ECfilt$EC_FCH4_comb <- USLA1_ECfilt$EC_FCH4_F_ANNOPTLM
USLA2_ECfilt$EC_FCH4_comb <- USLA2_ECfilt$EC_FCH4_F_ANNOPTLM
USLOS_ECfilt$EC_FCH4_comb <- USLOS_ECfilt$EC_FCH4_F_ANNOPTLM
USOWC_ECfilt$EC_FCH4_comb <- USOWC_ECfilt$EC_FCH4_F_ANNOPTLM
USUAF_ECfilt$EC_FCH4_comb <- USUAF_ECfilt$EC_FCH4_F_ANNOPTLM

USSTJ_ECfilt <- USSTJ_ECfilt %>%
  rename(
    EC_PA_F   = EC_PA
  )

USSTJ_ECfilt$ECCH_diff <- USSTJ_ECfilt$EC_FCH4_comb - USSTJ_ECfilt$ch_FCH4_nmolCH4m2s1

USSTJ_ECfilt <- subset(USSTJ_ECfilt, select = -c(EC_FCH4_F_RF, EC_FCH4_F_RF_filled, 
                                                                             EC_FCH4_F_MDS))

# add other temporal columns

CNHGU_ECfilt$Year <- year(CNHGU_ECfilt$TIMESTAMP_START)
CNHGU_ECfilt$DOY <- yday(CNHGU_ECfilt$TIMESTAMP_START)
CNHGU_ECfilt$Month <- month(CNHGU_ECfilt$TIMESTAMP_START)
CNHGU_ECfilt$Day <- day(CNHGU_ECfilt$TIMESTAMP_START)
CNHGU_ECfilt$Hour <- hour(CNHGU_ECfilt$TIMESTAMP_START)

FISI2_ECfilt$Hour <- NA
USLA1_ECfilt$Hour <- NA
USLA2_ECfilt$Hour <- NA

FISI2_ECfilt$Year <- year(FISI2_ECfilt$DATE)
FISI2_ECfilt$DOY <- yday(FISI2_ECfilt$DATE)
FISI2_ECfilt$Month <- month(FISI2_ECfilt$DATE)
FISI2_ECfilt$Day <- day(FISI2_ECfilt$DATE)

USLA1_ECfilt$Year <- year(USLA1_ECfilt$DATE)
USLA1_ECfilt$DOY <- yday(USLA1_ECfilt$DATE)
USLA1_ECfilt$Month <- month(USLA1_ECfilt$DATE)
USLA1_ECfilt$Day <- day(USLA1_ECfilt$DATE)

USLA2_ECfilt$Year <- year(USLA2_ECfilt$DATE)
USLA2_ECfilt$DOY <- yday(USLA2_ECfilt$DATE)
USLA2_ECfilt$Month <- month(USLA2_ECfilt$DATE)
USLA2_ECfilt$Day <- day(USLA2_ECfilt$DATE)

SEDEG_ECfilt$Year <- year(SEDEG_ECfilt$TIMESTAMP_START)
SEDEG_ECfilt$DOY <- yday(SEDEG_ECfilt$TIMESTAMP_START)
SEDEG_ECfilt$Month <- month(SEDEG_ECfilt$TIMESTAMP_START)
SEDEG_ECfilt$Day <- day(SEDEG_ECfilt$TIMESTAMP_START)
SEDEG_ECfilt$Hour <- hour(SEDEG_ECfilt$TIMESTAMP_START)

USHO1_ECfilt$Year <- year(USHO1_ECfilt$TIMESTAMP_START)
USHO1_ECfilt$DOY <- yday(USHO1_ECfilt$TIMESTAMP_START)
USHO1_ECfilt$Month <- month(USHO1_ECfilt$TIMESTAMP_START)
USHO1_ECfilt$Day <- day(USHO1_ECfilt$TIMESTAMP_START)
USHO1_ECfilt$Hour <- hour(USHO1_ECfilt$TIMESTAMP_START)

USOWC_ECfilt$Year <- year(USOWC_ECfilt$TIMESTAMP_START)
USOWC_ECfilt$DOY <- yday(USOWC_ECfilt$TIMESTAMP_START)
USOWC_ECfilt$Month <- month(USOWC_ECfilt$TIMESTAMP_START)
USOWC_ECfilt$Day <- day(USOWC_ECfilt$TIMESTAMP_START)
USOWC_ECfilt$Hour <- hour(USOWC_ECfilt$TIMESTAMP_START)

USUAF_ECfilt$Year <- year(USUAF_ECfilt$TIMESTAMP_START)
USUAF_ECfilt$DOY <- yday(USUAF_ECfilt$TIMESTAMP_START)
USUAF_ECfilt$Month <- month(USUAF_ECfilt$TIMESTAMP_START)
USUAF_ECfilt$Day <- day(USUAF_ECfilt$TIMESTAMP_START)
USUAF_ECfilt$Hour <- hour(USUAF_ECfilt$TIMESTAMP_START)

USLOS_ECfilt$Year <- year(USLOS_ECfilt$TIMESTAMP_START)
USLOS_ECfilt$DOY <- yday(USLOS_ECfilt$TIMESTAMP_START)
USLOS_ECfilt$Month <- month(USLOS_ECfilt$TIMESTAMP_START)
USLOS_ECfilt$Day <- day(USLOS_ECfilt$TIMESTAMP_START)
USLOS_ECfilt$Hour <- hour(USLOS_ECfilt$TIMESTAMP_START)

USSTJ_ECfilt$Year <- year(USSTJ_ECfilt$TIMESTAMP_START)
USSTJ_ECfilt$DOY <- yday(USSTJ_ECfilt$TIMESTAMP_START)
USSTJ_ECfilt$Month <- month(USSTJ_ECfilt$TIMESTAMP_START)
USSTJ_ECfilt$Day <- day(USSTJ_ECfilt$TIMESTAMP_START)
USSTJ_ECfilt$Hour <- hour(USSTJ_ECfilt$TIMESTAMP_START)

# rename the WD columns for FISI, USLA
FISI2_ECfilt <- FISI2_ECfilt %>% rename_with(~ c("EC_WD", "EC_WD_F"), all_of(c("EC_WD_AVG", "EC_WD_AVG_F")))
USLA1_ECfilt <- USLA1_ECfilt %>% rename_with(~ c("EC_WD", "EC_WD_F"), all_of(c("EC_WD_AVG", "EC_WD_AVG_F")))
USLA2_ECfilt <- USLA2_ECfilt %>% rename_with(~ c("EC_WD", "EC_WD_F"), all_of(c("EC_WD_AVG", "EC_WD_AVG_F")))

# remove the original flux values from dfs that have them
USLA1_ECfilt <- subset(USLA1_ECfilt, select = -c(ch_FCH4_ugCH4m2hr))
USLA2_ECfilt <- subset(USLA2_ECfilt, select = -c(ch_FCH4_ugCH4m2hr))
FISI2_ECfilt <- subset(FISI2_ECfilt, select = -c(ch_FCH4_mgCH4m2d1))

# remove the repeated EC data at ch_ columns in SEDEG
SEDEG_ECfilt <- subset(SEDEG_ECfilt, select = -c(ch_CH4_Ta_amb, ch_WS, ch_PPT, ch_P, ch_VPD, ch_RH))
SEDEG_ECfilt <- subset(SEDEG_ECfilt, select = -c(ch_CH4_PAR_mean, ch_CH4_PAR_amb))

# rename SEDEG PAR columns to make it clearer

# remove _UNC from all dfs

CNHGU_ECfilt <- subset(CNHGU_ECfilt, select = -c(EC_FCH4_F_ANNOPTLM_UNC, EC_FCH4_F_RANDUNC))
SEDEG_ECfilt <- subset(SEDEG_ECfilt, select = -c(EC_FCH4_F_ANNOPTLM_UNC, EC_FCH4_F_RANDUNC))
USHO1_ECfilt <- subset(USHO1_ECfilt, select = -c(EC_FCH4_F_ANNOPTLM_UNC, EC_FCH4_F_RANDUNC))
USLOS_ECfilt <- subset(USLOS_ECfilt, select = -c(EC_FCH4_F_ANNOPTLM_UNC, EC_FCH4_F_RANDUNC))
USOWC_ECfilt <- subset(USOWC_ECfilt, select = -c(EC_FCH4_F_ANNOPTLM_UNC, EC_FCH4_F_RANDUNC))
USUAF_ECfilt <- subset(USUAF_ECfilt, select = -c(EC_FCH4_F_ANNOPTLM_UNC, EC_FCH4_F_RANDUNC))

# FI-SI2, US-LA1, US-LA2:
# use medians, so delete means

USLA1_ECfilt <- USLA1_ECfilt %>% dplyr::select(-contains(c('_mean')))
USLA2_ECfilt <- USLA2_ECfilt %>% dplyr::select(-contains(c('_mean')))
FISI2_ECfilt <- FISI2_ECfilt %>% dplyr::select(-contains(c('_mean')))

# remove the _median from column names
names(FISI2_ECfilt) <- sub("_median", "", names(FISI2_ECfilt))
names(USLA1_ECfilt) <- sub("_median", "", names(USLA1_ECfilt))
names(USLA2_ECfilt) <- sub("_median", "", names(USLA2_ECfilt))

# add site info

# add dom_veg
CNHGU_ECfilt$DOM_VEG <- "AERENCHYMATOUS"

# add veg classes
CNHGU_ECfilt$MOSS_BROWN <- 0
CNHGU_ECfilt$MOSS_SPHAGNUM <- 0
CNHGU_ECfilt$AERENCHYMATOUS <- 1
CNHGU_ECfilt$ERI_SHRUB <- 0
CNHGU_ECfilt$TREE <- 0

# add climate
CNHGU_ECfilt$KOPPEN <- "Cwc"

# add ecosystem type
CNHGU_ECfilt$SITE_CLASSIFICATION <- "upland"

# add dom_veg
FISI2_ECfilt$DOM_VEG <- "MOSS_SPHAGNUM"

# add veg classes
FISI2_ECfilt$MOSS_BROWN <- 0
FISI2_ECfilt$MOSS_SPHAGNUM <- 1
FISI2_ECfilt$AERENCHYMATOUS <- 1
FISI2_ECfilt$ERI_SHRUB <- 1
FISI2_ECfilt$TREE <- 1

# add climate
FISI2_ECfilt$KOPPEN <- "Dfc"

# add ecosystem type
FISI2_ECfilt$SITE_CLASSIFICATION <- "bog"

# add dom_veg
USHO1_ECfilt$DOM_VEG <- "TREE"

# add veg classes
USHO1_ECfilt$MOSS_BROWN <- 0
USHO1_ECfilt$MOSS_SPHAGNUM <- 0
USHO1_ECfilt$AERENCHYMATOUS <- 0
USHO1_ECfilt$ERI_SHRUB <- 0
USHO1_ECfilt$TREE <- 1

# add climate
USHO1_ECfilt$KOPPEN <- "Dfb"

# add ecosystem type
USHO1_ECfilt$SITE_CLASSIFICATION <- "upland"

# add dom_veg
USLA1_ECfilt$DOM_VEG <- "AERENCHYMATOUS"

# add veg classes
USLA1_ECfilt$MOSS_BROWN <- 0
USLA1_ECfilt$MOSS_SPHAGNUM <- 0
USLA1_ECfilt$AERENCHYMATOUS <- 1
USLA1_ECfilt$ERI_SHRUB <- 0
USLA1_ECfilt$TREE <- 0

# add climate
USLA1_ECfilt$KOPPEN <- "Cfa"

# add ecosystem type
USLA1_ECfilt$SITE_CLASSIFICATION <- "salt marsh"

# # add dom_veg
USLA2_ECfilt$DOM_VEG <- "AERENCHYMATOUS"

# add veg classes
USLA2_ECfilt$MOSS_BROWN <- 0
USLA2_ECfilt$MOSS_SPHAGNUM <- 0
USLA2_ECfilt$AERENCHYMATOUS <- 1
USLA2_ECfilt$ERI_SHRUB <- 0
USLA2_ECfilt$TREE <- 0

# add climate
USLA2_ECfilt$KOPPEN <- "Cfa"

# add ecosystem type
USLA2_ECfilt$SITE_CLASSIFICATION <- "marsh"

# add dom_veg
USLOS_ECfilt$DOM_VEG <- "ERI_SHRUB"

# add veg classes
USLOS_ECfilt$MOSS_BROWN <- 0
USLOS_ECfilt$MOSS_SPHAGNUM <- 0
USLOS_ECfilt$AERENCHYMATOUS <- 1
USLOS_ECfilt$ERI_SHRUB <- 1
USLOS_ECfilt$TREE <- 1

# add climate
USLOS_ECfilt$KOPPEN <- "Dfb"

# add ecosystem type
USLOS_ECfilt$SITE_CLASSIFICATION <- "fen"

# add dom_veg
USOWC_ECfilt$DOM_VEG <- "AERENCHYMATOUS"

# add veg classes
USOWC_ECfilt$MOSS_BROWN <- 0
USOWC_ECfilt$MOSS_SPHAGNUM <- 0
USOWC_ECfilt$AERENCHYMATOUS <- 1
USOWC_ECfilt$ERI_SHRUB <- 0
USOWC_ECfilt$TREE <- 0

# add climate
USOWC_ECfilt$KOPPEN <- "Dfb"

# add ecosystem type
USOWC_ECfilt$SITE_CLASSIFICATION <- "marsh"

# add dom_veg
USUAF_ECfilt$DOM_VEG <- "MOSS_SPHAGNUM"

# add veg classes
USUAF_ECfilt$MOSS_BROWN <- 1
USUAF_ECfilt$MOSS_SPHAGNUM <- 1
USUAF_ECfilt$AERENCHYMATOUS <- 1
USUAF_ECfilt$ERI_SHRUB <- 1
USUAF_ECfilt$TREE <- 1

# add climate
USUAF_ECfilt$KOPPEN <- "Dwc"

# add ecosystem type
USUAF_ECfilt$SITE_CLASSIFICATION <- "bog"

# add dom_veg
SEDEG_ECfilt$DOM_VEG <- "MOSS_SPHAGNUM"

# add veg classes
SEDEG_ECfilt$MOSS_BROWN <- 0
SEDEG_ECfilt$MOSS_SPHAGNUM <- 1
SEDEG_ECfilt$AERENCHYMATOUS <- 1
SEDEG_ECfilt$ERI_SHRUB <- 1
SEDEG_ECfilt$TREE <- 0

# add climate
SEDEG_ECfilt$KOPPEN <- "Dfc"

# add ecosystem type
SEDEG_ECfilt$SITE_CLASSIFICATION <- "fen"

# add dom_veg
USSTJ_ECfilt$DOM_VEG <- "AERENCHYMATOUS"

# add veg classes
USSTJ_ECfilt$MOSS_BROWN <- 0
USSTJ_ECfilt$MOSS_SPHAGNUM <- 0
USSTJ_ECfilt$AERENCHYMATOUS <- 1
USSTJ_ECfilt$ERI_SHRUB <- 0
USSTJ_ECfilt$TREE <- 0

# add climate
USSTJ_ECfilt$KOPPEN <- "Cfa"

# add ecosystem type
USSTJ_ECfilt$SITE_CLASSIFICATION <- "salt marsh"

# add missing columns to dataframes
setdiff(names(SEDEG_ECfilt), names(CNHGU_ECfilt))
setdiff(names(USHO1_ECfilt), names(CNHGU_ECfilt))
setdiff(names(USLOS_ECfilt), names(CNHGU_ECfilt))
setdiff(names(USOWC_ECfilt), names(CNHGU_ECfilt))
setdiff(names(USUAF_ECfilt), names(CNHGU_ECfilt))

setdiff(names(USLA1_ECfilt), names(CNHGU_ECfilt))
setdiff(names(USLA2_ECfilt), names(CNHGU_ECfilt))
setdiff(names(FISI2_ECfilt), names(CNHGU_ECfilt))

USLA2_ECfilt <- subset(USLA2_ECfilt, select= -ch_AirTEddy)

# Define the common columns used to harmonize site-specific dataframes before row-binding

cols_1 <- c("ch_TS_2", "ch_TS_5","ch_TS_10", "ch_TS_15","ch_TS_20", "ch_TS_30", "ch_TS_40", 
            "EC_SWC_1", "EC_PA", "EC_TS_1", "EC_TS_2", 
            "EC_TS_3", "EC_TS_4", "EC_TS_5", "EC_TS_6", "EC_TS_7", "EC_TS_8", "EC_TS_9", 
            "EC_WTD", "EC_WTD_F", "EC_PPFD_OUT",
            "EC_SWC_2", "EC_SWC_3", "EC_PPFD_OUT_F", "ch_SM_X", 
            "ch_TS_X", 
            "ch_AerLAI", "u.wind.F", "v.wind.F")
CNHGU_ECfilt[cols_1] <- NA

CNHGU_ECfilt$ch_ID <- as.character(CNHGU_ECfilt$ch_ID)
CNHGU_ECfilt$DATE <- as_datetime(CNHGU_ECfilt$DATE)

cols <- CNHGU_ECfilt %>% dplyr::select(starts_with(c("EC_", "ch_TS", "ch_WTL", "ch_PAR", "ch_BM", "ch_GAI", "ch_FCH4", "ch_AerLAI", "ch_VWC",
                                                                   "ch_SM"))) %>% colnames()

CNHGU_ECfilt[cols] <- sapply(CNHGU_ECfilt[cols],as.numeric)

# now CNHGU is a df that includes all columns --> compare the rest to this

setdiff(names(CNHGU_ECfilt), names(SEDEG_ECfilt))

cols_2 <- setdiff(names(CNHGU_ECfilt), names(SEDEG_ECfilt))

SEDEG_ECfilt[cols_2] <- NA

SEDEG_ECfilt$ch_ID <- as.character(SEDEG_ECfilt$ch_ID)
SEDEG_ECfilt$DATE <- as_datetime(SEDEG_ECfilt$DATE)

SEDEG_ECfilt[cols] <- sapply(SEDEG_ECfilt[cols],as.numeric)

#---

setdiff(names(CNHGU_ECfilt), names(USHO1_ECfilt))

cols_3 <- setdiff(names(CNHGU_ECfilt), names(USHO1_ECfilt))

USHO1_ECfilt[cols_3] <- NA

USHO1_ECfilt$ch_ID <- as.character(USHO1_ECfilt$ch_ID)
USHO1_ECfilt$DATE <- as_datetime(USHO1_ECfilt$DATE)

USHO1_ECfilt[cols] <- sapply(USHO1_ECfilt[cols],as.numeric)

#---

setdiff(names(CNHGU_ECfilt), names(USLOS_ECfilt))

cols_4 <- setdiff(names(CNHGU_ECfilt), names(USLOS_ECfilt))

USLOS_ECfilt[cols_4] <- NA

USLOS_ECfilt$ch_ID <- as.character(USLOS_ECfilt$ch_ID)
USLOS_ECfilt$DATE <- as_datetime(USLOS_ECfilt$DATE)
USLOS_ECfilt[cols] <- sapply(USLOS_ECfilt[cols],as.numeric)

#---

setdiff(names(CNHGU_ECfilt), names(USOWC_ECfilt))

cols_5 <- setdiff(names(CNHGU_ECfilt), names(USOWC_ECfilt))

USOWC_ECfilt[cols_5] <- NA

USOWC_ECfilt$ch_ID <- as.character(USOWC_ECfilt$ch_ID)
USOWC_ECfilt$DATE <- as_datetime(USOWC_ECfilt$DATE)

USOWC_ECfilt[cols] <- sapply(USOWC_ECfilt[cols],as.numeric)

#---

setdiff(names(CNHGU_ECfilt), names(USUAF_ECfilt))

cols_6 <- setdiff(names(CNHGU_ECfilt), names(USUAF_ECfilt))

USUAF_ECfilt[cols_6] <- NA

USUAF_ECfilt$ch_ID <- as.character(USUAF_ECfilt$ch_ID)

USUAF_ECfilt$DATE<- as_datetime(USUAF_ECfilt$DATE)

USUAF_ECfilt[cols] <- sapply(USUAF_ECfilt[cols],as.numeric)

#---

setdiff(names(CNHGU_ECfilt), names(FISI2_ECfilt))

cols_7 <- setdiff(names(CNHGU_ECfilt), names(FISI2_ECfilt))

FISI2_ECfilt[cols_7] <- NA

FISI2_ECfilt$ch_ID <- as.character(FISI2_ECfilt$ch_ID)
FISI2_ECfilt$DATE <- as_datetime(FISI2_ECfilt$DATE)
cols <- FISI2_ECfilt %>% dplyr::select(starts_with(c("EC_", "ch_TS", "ch_WTL", "ch_PAR", "ch_BM", "ch_GAI", "ch_FCH4", "ch_AerLAI", "ch_VWC",
                                                                   "ch_SM"))) %>% colnames()
FISI2_ECfilt[cols] <- sapply(FISI2_ECfilt[cols],as.numeric)

#---

setdiff(names(CNHGU_ECfilt), names(USLA1_ECfilt))

cols_8 <- setdiff(names(CNHGU_ECfilt), names(USLA1_ECfilt))

USLA1_ECfilt[cols_8] <- NA

USLA1_ECfilt$ch_ID <- as.character(USLA1_ECfilt$ch_ID)
USLA1_ECfilt$DATE <- as_datetime(USLA1_ECfilt$DATE)
cols <- USLA1_ECfilt %>% dplyr::select(starts_with(c("EC_", "ch_TS", "ch_WTL", "ch_PAR", "ch_BM", "ch_GAI", "ch_FCH4", "ch_AerLAI", "ch_VWC",
                                                                   "ch_SM"))) %>% colnames()
USLA1_ECfilt[cols] <- sapply(USLA1_ECfilt[cols],as.numeric)

sapply(USLA1_ECfilt, class)

#---

setdiff(names(CNHGU_ECfilt), names(USLA2_ECfilt))

cols_9 <- setdiff(names(CNHGU_ECfilt), names(USLA2_ECfilt))

USLA2_ECfilt[cols_9] <- NA

USLA2_ECfilt$ch_ID <- as.character(USLA2_ECfilt$ch_ID)
USLA2_ECfilt$DATE <- as_datetime(USLA2_ECfilt$DATE)
cols <- USLA1_ECfilt %>% dplyr::select(starts_with(c("EC_", "ch_TS", "ch_WTL", "ch_PAR", "ch_BM", "ch_GAI", "ch_FCH4", "ch_AerLAI", "ch_VWC",
                                                                   "ch_SM"))) %>% colnames()
USLA2_ECfilt[cols] <- sapply(USLA2_ECfilt[cols],as.numeric)

sapply(USLA2_ECfilt, class)

#---

setdiff(names(CNHGU_ECfilt), names(USSTJ_ECfilt))

cols_9 <- setdiff(names(CNHGU_ECfilt), names(USSTJ_ECfilt))

USSTJ_ECfilt[cols_9] <- NA

USSTJ_ECfilt$ch_ID <- as.character(USSTJ_ECfilt$ch_ID)
USSTJ_ECfilt$DATE <- as_datetime(USSTJ_ECfilt$DATE)
cols <- USLA1_ECfilt %>% dplyr::select(starts_with(c("EC_", "ch_TS", "ch_WTL", "ch_PAR", "ch_BM", "ch_GAI", "ch_FCH4", "ch_AerLAI", "ch_VWC",
                                                                   "ch_SM"))) %>% colnames()
USSTJ_ECfilt[cols] <- sapply(USSTJ_ECfilt[cols],as.numeric)

#--------------------------------------------------------------

# combine

# only sites with automated chambers
auto_sites_raw <- bind_rows(CNHGU_ECfilt, SEDEG_ECfilt, USHO1_ECfilt,
                            USUAF_ECfilt)

# all sites:
all_sites_preds <- bind_rows(CNHGU_ECfilt, FISI2_ECfilt, SEDEG_ECfilt, USHO1_ECfilt,
                             USLA1_ECfilt, 
                             USLA2_ECfilt, 
                             USLOS_ECfilt, USOWC_ECfilt, USUAF_ECfilt,
                             USSTJ_ECfilt)

# remove rows where ch_FCH4 is NA

all_sites_preds <- all_sites_preds %>% drop_na(ch_FCH4_nmolCH4m2s1)
auto_sites_raw <- auto_sites_raw %>% drop_na(ch_FCH4_nmolCH4m2s1)

# add u and v components

# gap-filled (no gap-filled WD)
auto_sites_raw$u.wind.F <- -auto_sites_raw$EC_WS_F *sin(2*pi *auto_sites_raw$EC_WD/360)
auto_sites_raw$v.wind.F <- -auto_sites_raw$EC_WS_F *cos(2*pi *auto_sites_raw$EC_WD/360)

auto_sites_raw$u.wind <- -auto_sites_raw$EC_WS *sin(2*pi *auto_sites_raw$EC_WD/360)
auto_sites_raw$v.wind <- -auto_sites_raw$EC_WS *cos(2*pi *auto_sites_raw$EC_WD/360)

# save as .csv

write.csv(all_sites_preds, "path/allsites_raw_with_duplicates.csv", row.names = FALSE)
write.csv(auto_sites_raw, "path/auto_sites_raw_28082025.csv", row.names = FALSE)

# make a version without duplicates

# convert duplicates to NA in EC columns
# because fisi2, usla1 and usla2 are in date scale but the rest are at timestamp scale, 
# need to remove duplicates based on their different timescales

# Define the study sites that require DATE scale
specific_sites <- c("FI-SI2", "US-LA1", "US-LA2")

all_sites_preds[cols] <- sapply(all_sites_preds[cols], function(x) {
  ave(x, 
      ifelse(all_sites_preds$SITE %in% specific_sites, 
             as_datetime(all_sites_preds$DATE), 
             as_datetime(all_sites_preds$TIMESTAMP_START)),
      FUN = function(x) replace(x, duplicated(x), NA))
})

# save as .csv

write.csv(all_sites_preds, "path/allsites_raw_no_duplicates.csv", row.names = FALSE)

# with duplicates
all_sites_preds <- read.csv("path/allsites_raw_with_duplicates.csv")

# create soil temp and soil moist columns

# Count the number of non-NA values in the ch_TS columns based on SITE
result <- all_sites_preds %>%
  group_by(SITE) %>%
  summarize(across(starts_with("ch_TS_"), ~ sum(!is.na(.)), .names = "non_na_{col}"))

result <- as.data.frame(result)

# Create the new columns based on the "SITE" group
all_sites_preds <- all_sites_preds %>%
  mutate(
    ch_TS_TOP = case_when(
      SITE == 'FI-SI2' ~ ch_TS_5,
      SITE == 'CN-HGU' ~ NA,
      SITE == 'SE-DEG' ~ ch_TS_2,
      SITE == 'US-HO1' ~ ch_TS_X,
      SITE == 'US-LA1' ~ ch_TS_10,
      SITE == 'US-LA2' ~ ch_TS_10,
      SITE == 'US-LOS' ~ ch_TS_X,
      SITE == 'US-OWC' ~ NA,
      SITE == 'US-UAF' ~ ch_TS_2,
      TRUE ~ NA_real_
    ),
    EC_TS_TOP = case_when(
      SITE == 'FI-SI2' ~ EC_TS_1,
      SITE == 'CN-HGU' ~ EC_TS_F,
      SITE == 'SE-DEG' ~ EC_TS_1,
      SITE == 'US-HO1' ~ EC_TS_1,
      SITE == 'US-LA1' ~ EC_TS_F,
      SITE == 'US-LA2' ~ EC_TS_F,
      SITE == 'US-LOS' ~ EC_TS_2,
      SITE == 'US-OWC' ~ EC_TS_1,
      SITE == 'US-UAF' ~ EC_TS_1,
      TRUE ~ NA_real_
    )
  )

auto_sites_raw <- auto_sites_raw %>%
  mutate(
    ch_TS_TOP = case_when(
      SITE == 'CN-HGU' ~ NA,
      SITE == 'SE-DEG' ~ ch_TS_2,
      SITE == 'US-HO1' ~ ch_TS_X,
      SITE == 'US-UAF' ~ ch_TS_2
    ),
    EC_TS_TOP = case_when(
      SITE == 'CN-HGU' ~ EC_TS_F,
      SITE == 'SE-DEG' ~ EC_TS_1,
      SITE == 'US-HO1' ~ EC_TS_1,
      SITE == 'US-UAF' ~ EC_TS_1
    )
  )

# wtd

# create SITE_WTD

# EC_WTD is in m so convert EC_WTD to cm

auto_sites_raw$EC_WTD_F <- auto_sites_raw$EC_WTD_F * 100

# Create the new column SITE_WTD
auto_sites_raw <- auto_sites_raw %>%
  mutate(SITE_WTD = case_when(
    SITE %in% c("US-HO1", "US-UAF") ~ EC_WTD_F,
    SITE == "CN-HGU" ~ NA,
    TRUE ~ rowMeans(cbind(ch_WTL, EC_WTD_F), na.rm = TRUE)
  ))

# do the same for EC_SWC

# Count the number of non-NA values in the EC_SWC columns based on SITE
result <- all_sites_preds %>%
  group_by(SITE) %>%
  summarize(across(starts_with("EC_SWC_"), ~ sum(!is.na(.)), .names = "non_na_{col}"))

result <- as.data.frame(result)

# check ch_SM_X and ch_VWC

# Count the number of non-NA values in the EC_SWC columns based on SITE
result <- all_sites_preds %>%
  group_by(SITE) %>%
  summarize(across(c("ch_SM_X", "ch_VWC"), ~ sum(!is.na(.)), .names = "non_na_{col}"))

result <- as.data.frame(result)

# Create the new columns based on the "SITE" group
all_sites_preds <- all_sites_preds %>%
  mutate(
    ch_SWC_TOP = case_when(
      SITE == 'FI-SI2' ~ NA,
      SITE == 'CN-HGU' ~ EC_SWC_F,
      SITE == 'SE-DEG' ~ EC_SWC_1,
      SITE == 'US-HO1' ~ NA,
      SITE == 'US-LA1' ~ NA,
      SITE == 'US-LA2' ~ NA,
      SITE == 'US-LOS' ~ NA,
      SITE == 'US-OWC' ~ NA,
      SITE == 'US-UAF' ~ EC_SWC_1,
      TRUE ~ NA_real_
    ),
    EC_SWC_TOP = case_when(
      SITE == 'FI-SI2' ~ NA,
      SITE == 'CN-HGU' ~ NA,
      SITE == 'SE-DEG' ~ NA,
      SITE == 'US-HO1' ~ ch_SM_X,
      SITE == 'US-LA1' ~ NA,
      SITE == 'US-LA2' ~ NA,
      SITE == 'US-LOS' ~ ch_SM_X,
      SITE == 'US-OWC' ~ NA,
      SITE == 'US-UAF' ~ ch_VWC,
      TRUE ~ NA_real_
    )
  )

# Create the new columns based on the "SITE" group
auto_sites_raw <- auto_sites_raw %>%
  mutate(
    EC_SWC_TOP = case_when(
      SITE == 'CN-HGU' ~ EC_SWC_F,
      SITE == 'SE-DEG' ~ EC_SWC_1,
      SITE == 'US-HO1' ~ NA,
      SITE == 'US-UAF' ~ EC_SWC_1,
      TRUE ~ NA_real_
    ),
    ch_SWC_TOP = case_when(
      SITE == 'CN-HGU' ~ NA,
      SITE == 'SE-DEG' ~ NA,
      SITE == 'US-HO1' ~ ch_SM_X,
      SITE == 'US-UAF' ~ ch_VWC
    )
  )

auto_sites_raw <- auto_sites_raw %>%
  mutate(SITE_SWC = case_when(
    SITE %in% c("US-HO1") ~ ch_SWC_TOP,
    SITE %in% c("SE-DEG", "CN-HGU") ~ EC_SWC_TOP,
    TRUE ~ rowMeans(cbind(ch_SWC_TOP, EC_SWC_TOP), na.rm = TRUE)
  ))


# remove the extra columns

all_sites_preds <- subset(all_sites_preds, select = -c(EC_TS_1, EC_TS_2, EC_TS_3, EC_TS_4, EC_TS_5,
                                                       EC_TS_6, EC_TS_7, EC_TS_8, EC_TS_9,
                                                       ch_SM_X, ch_TS_X, ch_TS_2, ch_TS_5, ch_TS_10, 
                                                       ch_TS_15, ch_TS_20, ch_TS_30, ch_TS_40,
                                                       ch_VWC, EC_SWC_1, EC_SWC_2, EC_SWC_3))

all_sites_preds <- subset(all_sites_preds, select = -EC_SWC_F)

auto_sites_raw_new <- subset(auto_sites_raw, select = -c(EC_TS_1, EC_TS_2, EC_TS_3, EC_TS_4, EC_TS_5,
                                                         EC_TS_6, EC_TS_7, EC_TS_8, EC_TS_9,
                                                         ch_SM_X, ch_TS_X, ch_TS_2, ch_TS_5, ch_TS_10, 
                                                         ch_TS_15, ch_TS_20, ch_TS_30, ch_TS_40,
                                                         ch_VWC, EC_SWC_1, EC_SWC_2, EC_SWC_3))

#### HALF-HOURLY AGGREGATION ####

# CREATE A DF WITH CHAMBERS AGGREGATED TO EC TIMESTAMP LEVEL (for automated chamber sites and all sites separately)

# Define the specific study sites that require DATE scale
specific_sites <- c("FI-SI2", "US-LA1", "US-LA2")  

# Columns to aggregate
cols_to_aggregate <- c("ch_FCH4_nmolCH4m2s1", "ch_TS_TOP", "ch_SWC_TOP", "ch_WTL")

# Columns to exclude
cols_to_exclude <- c("ch_Datetime", "ch_ID", "ch_LC_class", "ch_Collar_type", "ch_AerLAI")

# Aggregation function
aggregate_site <- function(data, group_cols) {
  data %>%
    group_by(across(all_of(group_cols))) %>%
    summarize(across(all_of(cols_to_aggregate), 
                     list(mean = ~ mean(.x, na.rm = TRUE), median = ~ median(.x, na.rm = TRUE)),
                     .names = "{col}_{fn}"),
              across(-all_of(cols_to_aggregate), first),
              .groups = "drop")
}

# Aggregate by SITE and DATE for specific sites (first i only included DATE which removed a lot of data!)
date_aggregated <- all_sites_preds %>%
  filter(SITE %in% specific_sites) %>%
  mutate(DATE = as_date(ch_Datetime)) %>%
  aggregate_site(c("SITE", "DATE"))

# Aggregate by SITE and TIMESTAMP_START for other sites
timestamp_aggregated <- all_sites_preds %>%
  filter(!SITE %in% specific_sites) %>%
  mutate(TIMESTAMP_START = as_datetime(TIMESTAMP_START)) %>%
  aggregate_site(c("SITE", "TIMESTAMP_START"))

# AUTOMATED CHAMBER SITES:
# Aggregate by SITE and TIMESTAMP_START for other sites
timestamp_aggregated_auto <- auto_sites_raw_new %>%
  mutate(TIMESTAMP_START = as_datetime(TIMESTAMP_START)) %>%
  aggregate_site(c("SITE", "TIMESTAMP_START"))

timestamp_aggregated_auto <- as.data.frame(timestamp_aggregated_auto)

timestamp_aggregated_auto$TIMESTAMP_START <- as_datetime(timestamp_aggregated_auto$TIMESTAMP_START)
timestamp_aggregated_auto$TIMESTAMP_END <- as_datetime(timestamp_aggregated_auto$TIMESTAMP_END)
date_aggregated$DATE <- as_date(date_aggregated$DATE)
timestamp_aggregated$DATE <- as_date(timestamp_aggregated$DATE)
timestamp_aggregated_auto$DATE <- as_date(timestamp_aggregated_auto$DATE)
date_aggregated$TIMESTAMP_START <- as_datetime(date_aggregated$TIMESTAMP_START)
date_aggregated$TIMESTAMP_END <- as_datetime(date_aggregated$TIMESTAMP_END)
timestamp_aggregated$TIMESTAMP_START <- as_datetime(timestamp_aggregated$TIMESTAMP_START)
timestamp_aggregated$TIMESTAMP_END <- as_datetime(timestamp_aggregated$TIMESTAMP_END)

# Combine the two aggregated data frames
all_sites_preds_chaggr <- bind_rows(date_aggregated, timestamp_aggregated)

# Exclude the specified columns
all_sites_preds_chaggr <- all_sites_preds_chaggr %>%
  dplyr::select(-all_of(cols_to_exclude))

# AUTOMATED CHAMBER SITES:
# Exclude the specified columns
timestamp_aggregated_auto <- timestamp_aggregated_auto %>%
  dplyr::select(-all_of(cols_to_exclude))

# Ensure ch_WTL_mean and ch_WTL_median for site "CN-HGU" are NA
timestamp_aggregated_auto <- timestamp_aggregated_auto %>%
  mutate(ch_WTL_mean = ifelse(SITE == "CN-HGU", NA, ch_WTL_mean),
         ch_WTL_median = ifelse(SITE == "CN-HGU", NA, ch_WTL_median))

# Ordering the dataset
timestamp_aggregated_auto <-timestamp_aggregated_auto %>%
  arrange(SITE, as_datetime(TIMESTAMP_START))

timestamp_aggregated_auto <- as.data.frame(timestamp_aggregated_auto)

# relocate columns
timestamp_aggregated_auto <- timestamp_aggregated_auto %>% relocate(c(SITE, DATE, TIMESTAMP_START, TIMESTAMP_END, season, Year, DOY, Month, Day, Hour), .before = ch_FCH4_nmolCH4m2s1_mean)

# create a SITE_TS_TOP

timestamp_aggregated_auto <- timestamp_aggregated_auto %>%
  mutate(SITE_TS_TOP = case_when(
    SITE == "CN-HGU" ~ EC_TS_TOP,
    TRUE ~ rowMeans(cbind(ch_TS_TOP_median, EC_TS_TOP), na.rm = TRUE)
  ))

# Create the new column SITE_WTD

timestamp_aggregated_auto <- timestamp_aggregated_auto %>%
  mutate(SITE_WTD = case_when(
    SITE %in% c("US-HO1", "US-UAF") ~ EC_WTD_F,
    SITE == "CN-HGU" ~ NA,
    TRUE ~ rowMeans(cbind(ch_WTL_median, EC_WTD_F), na.rm = TRUE)
  ))

# Create SITE_SWC
timestamp_aggregated_auto <- timestamp_aggregated_auto %>%
  mutate(SITE_SWC = case_when(
    SITE %in% c("US-HO1") ~ ch_SWC_TOP_median,
    SITE %in% c("SE-DEG", "CN-HGU") ~ EC_SWC_TOP,
    TRUE ~ rowMeans(cbind(ch_SWC_TOP_median, EC_SWC_TOP), na.rm = TRUE)
  ))

# Remove TS with "_median" and SWC with mean

all_sites_preds_chaggr <- subset(all_sites_preds_chaggr, select = -c(ch_TS_TOP_median, ch_SWC_TOP_mean))

# remove ch_WTL_median

all_sites_preds_chaggr <- subset(all_sites_preds_chaggr, select = -ch_WTL_median)

# save

# all sites
write.csv(all_sites_preds_chaggr, "path/allsites_preds_chaggr.csv", row.names = FALSE)

# calculate delta FCH4
timestamp_aggregated_auto$deltaFCH4 <- timestamp_aggregated_auto$EC_FCH4_F_ANNOPTLM - timestamp_aggregated_auto$ch_FCH4_nmolCH4m2s1_median

# HALF-HOURLY AGGREGATION:
write.csv(timestamp_aggregated_auto, "path/autosites_preds_chaggr_NEW_VERSION_28082025.csv", row.names = FALSE)

#########################

###### HOURLY AGGREGATION ######

# read in datasets

CNHGU_hr_aggr_ECfilt <- read.csv("path/CN_HGU_combined_filtered_ECfilt_with_outliers_ch_hraggr.csv")

SEDEG_hr_aggr_ECfilt <- read.csv("path/SE_DEG_combined_filtered_with_outliers_hraggr.csv")

USHO1_hr_aggr_ECfilt <- read.csv("path/US_HO1_combined_filtered_with_outliers_ch_hraggr.csv")

USUAF_hr_aggr_ECfilt <- read.csv("path/US_UAF_combined_filtered_with_outliers_hraggr_ECfilt.csv")

# harmonize datasets

USHO1_hr_aggr_ECfilt <- USHO1_hr_aggr_ECfilt %>% rename("ch_FCH4_nmolCH4m2s1_mean" = "ch_FCH4_umolCH4m2s1_mean",
                                                        "ch_FCH4_nmolCH4m2s1_median" = "ch_FCH4_umolCH4m2s1_median"
)

SEDEG_hr_aggr_ECfilt <- subset(SEDEG_hr_aggr_ECfilt, select=-c(ch_PPT_mean, ch_WS_mean, ch_P_mean, ch_VPD_mean, ch_RH_mean,
                                                               ch_PPT_median, ch_WS_median, ch_P_median, ch_VPD_median, ch_RH_median))

# remove UNC columns

CNHGU_hr_aggr_ECfilt <- subset(CNHGU_hr_aggr_ECfilt, select = -c(EC_FCH4_F_ANNOPTLM_UNC_median, EC_FCH4_F_ANNOPTLM_UNC_mean))
SEDEG_hr_aggr_ECfilt <- subset(SEDEG_hr_aggr_ECfilt, select = -c(EC_FCH4_F_ANNOPTLM_UNC_median, EC_FCH4_F_ANNOPTLM_UNC_mean))
USHO1_hr_aggr_ECfilt <- subset(USHO1_hr_aggr_ECfilt, select = -c(EC_FCH4_F_ANNOPTLM_UNC_median, EC_FCH4_F_ANNOPTLM_UNC_mean))
USUAF_hr_aggr_ECfilt <- subset(USUAF_hr_aggr_ECfilt, select = -c(EC_FCH4_F_ANNOPTLM_UNC_median, EC_FCH4_F_ANNOPTLM_UNC_mean))

# add missing columns to dataframes
setdiff(names(SEDEG_hr_aggr_ECfilt), names(CNHGU_hr_aggr_ECfilt))
setdiff(names(USHO1_hr_aggr_ECfilt), names(CNHGU_hr_aggr_ECfilt))
setdiff(names(USUAF_hr_aggr_ECfilt), names(CNHGU_hr_aggr_ECfilt))
setdiff(names(USOWC_hh_aggr_ECfilt), names(CNHGU_hh_aggr_ECfilt))
setdiff(names(USUAF_hh_aggr_ECfilt), names(CNHGU_hh_aggr_ECfilt))

cols_1 <- c("ch_TS_2_mean", "ch_TS_10_mean", "ch_TS_20_mean", "ch_TS_30_mean", "ch_TS_40_mean", 
            "ch_CH4_Ta_amb_mean", "ch_TS_2_median",
            "EC_SWC_1_mean", "EC_PA_mean", "EC_TS_1_mean", "EC_TS_2_mean", 
            "EC_TS_3_mean", "EC_TS_4_mean", "EC_TS_5_mean", "EC_TS_6_mean", "EC_TS_7_mean", "EC_TS_8_mean", "EC_TS_9_mean", 
            "EC_WTD_F_mean", "EC_PPFD_OUT_mean",
            "EC_SWC_2_mean", "EC_SWC_3_mean", "EC_PPFD_OUT_F_mean", 
            "EC_SWC_1_median", "EC_PA_median", "EC_TS_1_median", "EC_TS_2_median", 
            "EC_TS_3_median", "EC_TS_4_median", "EC_TS_5_median", "EC_TS_6_median", "EC_TS_7_median", "EC_TS_8_median", "EC_TS_9_median", 
            "EC_WTD_F_median", "EC_PPFD_OUT_median",
            "EC_SWC_2_median", "EC_SWC_3_median", "EC_PPFD_OUT_F_median",
            "ch_SM_X_mean", 
            "ch_TS_X_mean", "ch_SM_X_median", "ch_TS_X_median",
            "ch_VWC_mean",
            "ch_TS_20_median", "ch_TS_30_median", "ch_TS_40_median", "ch_VWC_median"
)
CNHGU_hr_aggr_ECfilt[cols_1] <- NA

# now CNHGU is a df that includes all columns --> compare the rest to this

setdiff(names(CNHGU_hr_aggr_ECfilt), names(SEDEG_hr_aggr_ECfilt))

cols_2 <- setdiff(names(CNHGU_hr_aggr_ECfilt), names(SEDEG_hr_aggr_ECfilt))

SEDEG_hr_aggr_ECfilt[cols_2] <- NA

#---

setdiff(names(CNHGU_hr_aggr_ECfilt), names(USHO1_hr_aggr_ECfilt))

cols_3 <- setdiff(names(CNHGU_hr_aggr_ECfilt), names(USHO1_hr_aggr_ECfilt))

USHO1_hr_aggr_ECfilt[cols_3] <- NA

#---

setdiff(names(CNHGU_hr_aggr_ECfilt), names(USUAF_hr_aggr_ECfilt))

cols_4 <- setdiff(names(CNHGU_hr_aggr_ECfilt), names(USUAF_hr_aggr_ECfilt))

USUAF_hr_aggr_ECfilt[cols_4] <- NA

# combine!

all_sites_hourly_ECfilt <- bind_rows(CNHGU_hr_aggr_ECfilt, SEDEG_hr_aggr_ECfilt, USHO1_hr_aggr_ECfilt,
                                     USUAF_hr_aggr_ECfilt)

# add Year, Month, Day, Hour and DOY columns

all_sites_hourly_ECfilt$TIMESTAMP <- as_datetime(all_sites_hourly_ECfilt$TIMESTAMP)

all_sites_hourly_ECfilt$Year <- year(all_sites_hourly_ECfilt$TIMESTAMP)

all_sites_hourly_ECfilt$Month <- month(all_sites_hourly_ECfilt$TIMESTAMP)

all_sites_hourly_ECfilt$Day <- day(all_sites_hourly_ECfilt$TIMESTAMP)

all_sites_hourly_ECfilt$Hour <- hour(all_sites_hourly_ECfilt$TIMESTAMP)

all_sites_hourly_ECfilt$DOY <- yday(all_sites_hourly_ECfilt$TIMESTAMP)

# combine TS info into one column

# Create the new columns based on the "SITE" group
all_sites_hourly_ECfilt <- all_sites_hourly_ECfilt %>%
  mutate(
    ch_TS_TOP = case_when(
      SITE == 'CN-HGU' ~ NA,
      SITE == 'SE-DEG' ~ ch_TS_2_mean,
      SITE == 'US-HO1' ~ ch_TS_X_mean,
      SITE == 'US-UAF' ~ ch_TS_2_mean,
      TRUE ~ NA_real_
    ),
    EC_TS_TOP = case_when(
      SITE == 'CN-HGU' ~ EC_TS_F_mean,
      SITE == 'SE-DEG' ~ EC_TS_1_mean,
      SITE == 'US-HO1' ~ EC_TS_1_mean,
      SITE == 'US-UAF' ~ EC_TS_1_mean,
      TRUE ~ NA_real_
    )
  )

# create a SITE_TS_TOP

all_sites_hourly_ECfilt <- all_sites_hourly_ECfilt %>%
  mutate(SITE_TS_TOP = case_when(
    SITE == "CN-HGU" ~ EC_TS_TOP,
    TRUE ~ rowMeans(cbind(ch_TS_TOP, EC_TS_TOP), na.rm = TRUE)
  ))

# do the same for EC_SWC

# Create the new columns based on the "SITE" group
all_sites_hourly_ECfilt <- all_sites_hourly_ECfilt %>%
  mutate(
    EC_SWC_TOP = case_when(
      SITE == 'CN-HGU' ~ EC_SWC_F_mean,
      SITE == 'SE-DEG' ~ EC_SWC_1_mean,
      SITE == 'US-HO1' ~ NA,
      SITE == 'US-UAF' ~ EC_SWC_1_mean,
      TRUE ~ NA_real_
    ),
    ch_SWC_TOP = case_when(
      SITE == 'CN-HGU' ~ NA,
      SITE == 'SE-DEG' ~ NA,
      SITE == 'US-HO1' ~ ch_SM_X_mean,
      SITE == 'US-UAF' ~ ch_VWC_mean,
      TRUE ~ NA_real_
    )
  )

# create SITE_SWC

all_sites_hourly_ECfilt <- all_sites_hourly_ECfilt %>%
  mutate(SITE_SWC = case_when(
    SITE %in% c("US-HO1") ~ ch_SWC_TOP,
    SITE %in% c("SE-DEG", "CN-HGU") ~ EC_SWC_TOP,
    TRUE ~ rowMeans(cbind(ch_SWC_TOP, EC_SWC_TOP), na.rm = TRUE)
  ))

# create SITE_WTD

# EC_WTD is in m so convert EC_WTD to cm

all_sites_hourly_ECfilt$EC_WTD_F <- all_sites_hourly_ECfilt$EC_WTD_F_mean * 100

# Create the new column SITE_WTD
all_sites_hourly_ECfilt <- all_sites_hourly_ECfilt %>%
  mutate(SITE_WTD = case_when(
    SITE %in% c("US-HO1", "US-UAF") ~ EC_WTD_F,
    SITE == "CN-HGU" ~ NA,
    TRUE ~ rowMeans(cbind(ch_WTL_mean, EC_WTD_F), na.rm = TRUE)
  ))

# remove the extra columns

all_sites_hourly_ECfilt <- subset(all_sites_hourly_ECfilt, select = -c(EC_TS_1_mean, EC_TS_2_mean, EC_TS_3_mean, EC_TS_4_mean, EC_TS_5_mean,
                                                                       EC_TS_6_mean, EC_TS_7_mean, EC_TS_8_mean, EC_TS_9_mean,
                                                                       ch_SM_X_mean, ch_TS_X_mean, ch_TS_2_mean, ch_TS_10_mean, 
                                                                       ch_TS_20_mean, ch_TS_30_mean, ch_TS_40_mean,
                                                                       ch_VWC_mean, EC_SWC_1_mean, EC_SWC_2_mean, EC_SWC_3_mean))

all_sites_hourly_ECfilt <- subset(all_sites_hourly_ECfilt, select = -c(EC_TS_1_median, EC_TS_2_median, EC_TS_3_median, EC_TS_4_median, EC_TS_5_median,
                                                                       EC_TS_6_median, EC_TS_7_median, EC_TS_8_median, EC_TS_9_median,
                                                                       ch_SM_X_median, ch_TS_X_median, ch_TS_2_median, ch_TS_10_median, 
                                                                       ch_TS_20_median, ch_TS_30_median, ch_TS_40_median,
                                                                       ch_VWC_median, EC_SWC_1_median, EC_SWC_2_median, EC_SWC_3_median))

# calculate delta FCH4
all_sites_hourly_ECfilt$deltaFCH4 <- all_sites_hourly_ECfilt$EC_FCH4_F_ANNOPTLM_median - all_sites_hourly_ECfilt$ch_FCH4_nmolCH4m2s1_median

# save as .csv

write.csv(all_sites_hourly_ECfilt, "path/autosites_hraggr_ECfilt.csv", row.names = FALSE)

#########################

##### DAILY AGGREGATION #####

# read in datasets

CNHGU_d_aggr_ECfilt <- read.csv("path/CN_HGU_chec_daggr_ECfilt.csv")

SEDEG_d_aggr_ECfilt <- read.csv("path/SE_DEG_chec_daggr_ECfilt.csv")

USHO1_d_aggr_ECfilt <- read.csv("path/US_HO1_chec_daggr_ECfilt.csv")

USUAF_d_aggr_ECfilt <- read.csv("path/US_UAF_chec_daggr_ECfilt.csv")

FISI2_d_aggr_ECfilt <- read.csv("path/FI_SI2_chec_daggr_ECfilt.csv")

USLA1_d_aggr_ECfilt <- read.csv("path/US_LA1_chec_daggr_ECfilt.csv")

USLA2_d_aggr_ECfilt <- read.csv("path/US_LA2_chec_daggr_ECfilt.csv")

USLOS_d_aggr_ECfilt <- read.csv("path/US_LOS_chec_daggr_ECfilt.csv")

USOWC_d_aggr_ECfilt <- read.csv("path/US_OWC_chec_daggr_ECfilt.csv")

# harmonize datasets

USHO1_d_aggr_ECfilt <- USHO1_d_aggr_ECfilt %>% dplyr::rename("ch_FCH4_nmolCH4m2s1_mean" = "ch_FCH4_umolCH4m2s1_mean",
                                                             
                                                             "ch_FCH4_nmolCH4m2s1_median" = "ch_FCH4_umolCH4m2s1_median"
                                                             
)

SEDEG_d_aggr_ECfilt <- subset(SEDEG_d_aggr_ECfilt, select=-c(ch_PPT_mean, ch_WS_mean, ch_P_mean, ch_VPD_mean, ch_RH_mean,
                                                             ch_PPT_median, ch_WS_median, ch_P_median, ch_VPD_median, ch_RH_median,
                                                             ch_CH4_Ta_amb_median, ch_CH4_Ta_amb_mean, ch_CH4_PAR_amb_median, ch_CH4_PAR_amb_mean))

# rename the date columns

colnames(SEDEG_d_aggr_ECfilt)[2] <- "TIMESTAMP"
colnames(CNHGU_d_aggr_ECfilt)[2] <- "TIMESTAMP"
colnames(USHO1_d_aggr_ECfilt)[2] <- "TIMESTAMP"
colnames(USLOS_d_aggr_ECfilt)[2] <- "TIMESTAMP"
colnames(USOWC_d_aggr_ECfilt)[2] <- "TIMESTAMP"
colnames(USUAF_d_aggr_ECfilt)[2] <- "TIMESTAMP"

FISI2_d_aggr_ECfilt <- FISI2_d_aggr_ECfilt %>% dplyr::rename("ch_WTL_mean" = "ch_WT_mean",
                                                             "ch_WTL_median" = "ch_WT_median"
)

# remove ch_TS_5 and 15

FISI2_d_aggr_ECfilt <- subset(FISI2_d_aggr_ECfilt, select=-c(ch_TS_15_mean,  ch_TS_15_median))
USLA1_d_aggr_ECfilt <- subset(USLA1_d_aggr_ECfilt, select=-c(ch_TS_5_mean,  ch_TS_5_median,  ch_TS_15_mean,  ch_TS_15_median))
USLA2_d_aggr_ECfilt <- subset(USLA2_d_aggr_ECfilt, select=-c(ch_TS_5_mean,  ch_TS_5_median,  ch_TS_15_mean,  ch_TS_15_median))

cols_1 <- c("ch_TS_2_mean", "ch_TS_10_mean", "ch_TS_20_mean", "ch_TS_30_mean", "ch_TS_40_mean", 
            "ch_TS_2_median",
            "ch_TS_10_median",  "ch_WTL_median",
            
            
            "EC_SWC_1_mean", "EC_PA_mean", "EC_TS_1_mean", "EC_TS_2_mean", 
            "EC_TS_3_mean", "EC_TS_4_mean", "EC_TS_5_mean", "EC_TS_6_mean", "EC_TS_7_mean", "EC_TS_8_mean", "EC_TS_9_mean", 
            "EC_WTD_F_mean", "EC_PPFD_OUT_mean",
            "EC_SWC_2_mean", "EC_SWC_3_mean", "EC_PPFD_OUT_F_mean", 
            "EC_SWC_1_median", "EC_PA_median", "EC_TS_1_median", "EC_TS_2_median", 
            "EC_TS_3_median", "EC_TS_4_median", "EC_TS_5_median", "EC_TS_6_median", "EC_TS_7_median", "EC_TS_8_median", "EC_TS_9_median", 
            "EC_WTD_F_median", "EC_PPFD_OUT_median",
            "EC_SWC_2_median", "EC_SWC_3_median", "EC_PPFD_OUT_F_median",
            
            "ch_SM_X_mean", 
            "ch_TS_X_mean",   "ch_SM_X_median", "ch_TS_X_median", 
            "ch_VWC_mean",
            "ch_TS_20_median", "ch_TS_30_median", "ch_TS_40_median", "ch_VWC_median")

CNHGU_d_aggr_ECfilt[cols_1] <- NA

# now CNHGU is a df that includes all columns --> compare the rest to this

setdiff(names(CNHGU_d_aggr_ECfilt), names(SEDEG_d_aggr_ECfilt))

cols_2 <- setdiff(names(CNHGU_d_aggr_ECfilt), names(SEDEG_d_aggr_ECfilt))

SEDEG_d_aggr_ECfilt[cols_2] <- NA

#---

setdiff(names(CNHGU_d_aggr_ECfilt), names(USHO1_d_aggr_ECfilt))

cols_3 <- setdiff(names(CNHGU_d_aggr_ECfilt), names(USHO1_d_aggr_ECfilt))

USHO1_d_aggr_ECfilt[cols_3] <- NA

#---
setdiff(names(CNHGU_d_aggr_ECfilt), names(USUAF_d_aggr_ECfilt))

cols_4 <- setdiff(names(CNHGU_d_aggr_ECfilt), names(USUAF_d_aggr_ECfilt))

USUAF_d_aggr_ECfilt[cols_4] <- NA

#---
setdiff(names(CNHGU_d_aggr_ECfilt), names(FISI2_d_aggr_ECfilt))

cols_5 <- setdiff(names(CNHGU_d_aggr_ECfilt), names(FISI2_d_aggr_ECfilt))

FISI2_d_aggr_ECfilt[cols_5] <- NA

#---
setdiff(names(CNHGU_d_aggr_ECfilt), names(USLA1_d_aggr_ECfilt))

cols_6 <- setdiff(names(CNHGU_d_aggr_ECfilt), names(USLA1_d_aggr_ECfilt))

USLA1_d_aggr_ECfilt[cols_6] <- NA

#---
setdiff(names(CNHGU_d_aggr_ECfilt), names(USLA2_d_aggr_ECfilt))

cols_7 <- setdiff(names(CNHGU_d_aggr_ECfilt), names(USLA2_d_aggr_ECfilt))

USLA2_d_aggr_ECfilt[cols_7] <- NA

#---
setdiff(names(CNHGU_d_aggr_ECfilt), names(USLOS_d_aggr_ECfilt))

cols_8 <- setdiff(names(CNHGU_d_aggr_ECfilt), names(USLOS_d_aggr_ECfilt))

USLOS_d_aggr_ECfilt[cols_8] <- NA

#---

setdiff(names(CNHGU_d_aggr_ECfilt), names(USOWC_d_aggr_ECfilt))

cols_9 <- setdiff(names(CNHGU_d_aggr_ECfilt), names(USOWC_d_aggr_ECfilt))

USOWC_d_aggr_ECfilt[cols_9] <- NA

# combine

all_sites_daily_ECfilt <- bind_rows(CNHGU_d_aggr_ECfilt, FISI2_d_aggr_ECfilt, SEDEG_d_aggr_ECfilt, USHO1_d_aggr_ECfilt,
                                    USLA1_d_aggr_ECfilt, USLA2_d_aggr_ECfilt, USLOS_d_aggr_ECfilt, USOWC_d_aggr_ECfilt, USUAF_d_aggr_ECfilt)

# add Year, Month, Day, Hour and DOY columns

all_sites_daily_ECfilt$TIMESTAMP <- as_datetime(all_sites_daily_ECfilt$TIMESTAMP)

all_sites_daily_ECfilt$Year <- year(all_sites_daily_ECfilt$TIMESTAMP)

all_sites_daily_ECfilt$Month <- month(all_sites_daily_ECfilt$TIMESTAMP)

all_sites_daily_ECfilt$Day <- day(all_sites_daily_ECfilt$TIMESTAMP)

all_sites_daily_ECfilt$DOY <- yday(all_sites_daily_ECfilt$TIMESTAMP)

# add ch_method

all_sites_daily_ECfilt <- all_sites_daily_ECfilt %>% mutate(ch_method =
                                                              case_when(all_sites_daily_ECfilt$SITE == "CN-HGU" ~ "auto", 
                                                                        all_sites_daily_ECfilt$SITE == "FI-SI2" ~ "manual",
                                                                        all_sites_daily_ECfilt$SITE == "SE-DEG" ~ "auto",
                                                                        all_sites_daily_ECfilt$SITE == "US-HO1" ~ "auto",
                                                                        all_sites_daily_ECfilt$SITE == "US-LA1" ~ "manual",
                                                                        all_sites_daily_ECfilt$SITE == "US-LA2" ~ "manual",
                                                                        all_sites_daily_ECfilt$SITE == "US-LOS" ~ "manual",
                                                                        all_sites_daily_ECfilt$SITE == "US-OWC" ~ "manual",
                                                                        all_sites_daily_ECfilt$SITE == "US-UAF" ~ "auto")
)

#relocate 

all_sites_daily_ECfilt <- all_sites_daily_ECfilt %>% relocate(ch_method, .after = TIMESTAMP)

# remove the ANN UNC columns
all_sites_daily_ECfilt <- subset(all_sites_daily_ECfilt, select= -c(EC_FCH4_F_ANNOPTLM_UNC_mean,  EC_FCH4_F_ANNOPTLM_UNC_median))

# Create the new columns based on the "SITE" group
all_sites_daily_ECfilt <- all_sites_daily_ECfilt %>%
  mutate(
    ch_TS_TOP = case_when(
      SITE == 'FI-SI2' ~ ch_TS_5_mean,
      SITE == 'CN-HGU' ~ NA,
      SITE == 'SE-DEG' ~ ch_TS_2_mean,
      SITE == 'US-HO1' ~ ch_TS_X_mean,
      SITE == 'US-LA1' ~ ch_TS_10_mean,
      SITE == 'US-LA2' ~ ch_TS_10_mean,
      SITE == 'US-LOS' ~ ch_TS_X_mean,
      SITE == 'US-OWC' ~ NA,
      SITE == 'US-UAF' ~ ch_TS_2_mean,
      TRUE ~ NA_real_
    ),
    EC_TS_TOP = case_when(
      SITE == 'FI-SI2' ~ EC_TS_1_mean,
      SITE == 'CN-HGU' ~ EC_TS_F_mean,
      SITE == 'SE-DEG' ~ EC_TS_1_mean,
      SITE == 'US-HO1' ~ EC_TS_1_mean,
      SITE == 'US-LA1' ~ EC_TS_F_mean,
      SITE == 'US-LA2' ~ EC_TS_F_mean,
      SITE == 'US-LOS' ~ EC_TS_1_mean,
      SITE == 'US-OWC' ~ EC_TS_1_mean,
      SITE == 'US-UAF' ~ EC_TS_1_mean,
      TRUE ~ NA_real_
    )
  )

# do the same for EC_SWC

# Create the new columns based on the "SITE" group
all_sites_daily_ECfilt <- all_sites_daily_ECfilt %>%
  mutate(
    EC_SWC_TOP = case_when(
      SITE == 'FI-SI2' ~ NA,
      SITE == 'CN-HGU' ~ EC_SWC_F_mean,
      SITE == 'SE-DEG' ~ EC_SWC_1_mean,
      SITE == 'US-HO1' ~ NA,
      SITE == 'US-LA1' ~ NA,
      SITE == 'US-LA2' ~ NA,
      SITE == 'US-LOS' ~ NA,
      SITE == 'US-OWC' ~ NA,
      SITE == 'US-UAF' ~ EC_SWC_1_mean,
      TRUE ~ NA_real_
    ),
    ch_SWC_TOP = case_when(
      SITE == 'FI-SI2' ~ NA,
      SITE == 'CN-HGU' ~ NA,
      SITE == 'SE-DEG' ~ NA,
      SITE == 'US-HO1' ~ ch_SM_X_mean,
      SITE == 'US-LA1' ~ NA,
      SITE == 'US-LA2' ~ NA,
      SITE == 'US-LOS' ~ ch_SM_X_mean,
      SITE == 'US-OWC' ~ NA,
      SITE == 'US-UAF' ~ ch_VWC_mean,
      TRUE ~ NA_real_
    )
  )

# remove the extra columns

all_sites_daily_ECfilt <- subset(all_sites_daily_ECfilt, select = -c(EC_TS_1_mean, EC_TS_2_mean, EC_TS_3_mean, EC_TS_4_mean, EC_TS_5_mean,
                                                                     EC_TS_6_mean, EC_TS_7_mean, EC_TS_8_mean, EC_TS_9_mean,
                                                                     ch_SM_X_mean, ch_TS_X_mean, ch_TS_2_mean, ch_TS_5_mean, ch_TS_10_mean, 
                                                                     ch_TS_20_mean, ch_TS_30_mean, ch_TS_40_mean,
                                                                     ch_VWC_mean, EC_SWC_1_mean, EC_SWC_2_mean, EC_SWC_3_mean))
all_sites_daily_ECfilt <- subset(all_sites_daily_ECfilt, select = -c(EC_TS_1_median, EC_TS_2_median, EC_TS_3_median, EC_TS_4_median, EC_TS_5_median,
                                                                     EC_TS_6_median, EC_TS_7_median, EC_TS_8_median, EC_TS_9_median,
                                                                     ch_SM_X_median, ch_TS_X_median, ch_TS_2_median, ch_TS_5_median, ch_TS_10_median, 
                                                                     ch_TS_20_median, ch_TS_30_median, ch_TS_40_median,
                                                                     ch_VWC_median, EC_SWC_1_median, EC_SWC_2_median, EC_SWC_3_median))

all_sites_daily_ECfilt <- subset(all_sites_daily_ECfilt, select = -c(ch_TS_15_mean, ch_WT_median, ch_TS_15_median))

all_sites_daily_ECfilt <- subset(all_sites_daily_ECfilt, select = -EC_SWC_F_mean)

# order based on alphabetical order
all_sites_daily_ECfilt <- all_sites_daily_ECfilt %>%
  arrange(SITE, TIMESTAMP)

###### UPDATED DATA: ADD US-STJ
# This section was run with updated US-STJ data (already included in the published datasets in aggregated format)
# To rebuild from source files, add the updated US-STJ data analogously to the site-specific sections above.

USSTJ_daggr <- read.csv("path/US_STJ_chec_daggr.csv")

all_sites_preds_daggr <- all_sites_daily_ECfilt

# change some column names

USSTJ_daggr <- USSTJ_daggr %>%
  rename(
    EC_NEE_F_median   = EC_NEE_F_MDS_median,
    EC_NEE_F_mean     = EC_NEE_F_MDS_mean,
    EC_FCH4_median = EC_FCH4_nmolCH4m2s1_median,
    EC_FCH4_mean = EC_FCH4_nmolCH4m2s1_mean, 
    TIMESTAMP = DATE
  )

# remove TA and EC_NEE (gap-filled NEE better)

USSTJ_daggr <- subset(USSTJ_daggr, select = -c(EC_TA_mean, EC_TA_median, EC_NEE_mean, EC_NEE_median))

USSTJ_daggr <- USSTJ_daggr %>%
  rename(
    EC_WS_F_mean   = EC_WS_mean, # not really gap-filled but needs to be in the same column as in rest of the sites
    EC_WS_F_median     = EC_WS_median, # --"--
    EC_VPD_F_mean = EC_VPD_mean, #--"--
    EC_VPD_F_median = EC_VPD_median # --"--
  )

all_sites_preds_daggr %>%
  filter(!is.na(EC_WD_u.wind_mean))

# need to move values from EC_WD_u.wind etc columns to the right columns in all sites preds
# Pairs of (source column, target column)
col_pairs <- list(
  c("EC_WD_u.wind_mean",      "u.wind_mean"),
  c("EC_WD_u.wind_median",    "u.wind_median"),
  c("EC_WD_u.wind_mean_F",    "u.wind.F_mean"),
  c("EC_WD_u.wind_median_F",  "u.wind.F_median"),
  c("EC_WD_v.wind_mean",      "v.wind_mean"),
  c("EC_WD_v.wind_median",    "v.wind_median"),
  c("EC_WD_v.wind_mean_F",    "v.wind.F_mean"),
  c("EC_WD_v.wind_median_F",  "v.wind.F_median")
)

# Fill missing target values from source values when available
for (pair in col_pairs) {
  src <- pair[1]
  tgt <- pair[2]
  
  idx <- is.na(all_sites_preds_daggr[[tgt]]) & !is.na(all_sites_preds_daggr[[src]])
  all_sites_preds_daggr[[tgt]][idx] <- all_sites_preds_daggr[[src]][idx]
}

# Remove EC_WD_... columns except EC_WD_AVG and EC_WD_AVG_F
cols_to_remove <- grep("^EC_WD_", names(all_sites_preds_daggr), value = TRUE)
cols_to_keep   <- c("EC_WD_AVG", "EC_WD_AVG_F")
cols_to_remove <- setdiff(cols_to_remove, cols_to_keep)

all_sites_preds_daggr <- all_sites_preds_daggr[ , !(names(all_sites_preds_daggr) %in% cols_to_remove)]

USSTJ_daggr <- USSTJ_daggr %>%
  rename(
    v.wind.F_mean   = v.wind_mean, # not really gap-filled but needs to be in the same column as in rest of the sites
    u.wind.F_mean     = u.wind_mean
  )

USSTJ_daggr <- USSTJ_daggr %>%
  rename(
    v.wind.F_median   = v.wind_median, # not really gap-filled but needs to be in the same column as in rest of the sites
    u.wind.F_median     = u.wind_median
  )

### need to combine FCH4 values from USSTJ to ANNOPTLM values from other sites
# --> create new column for this

USSTJ_daggr$EC_FCH4_comb_median <- USSTJ_daggr$EC_FCH4_median
all_sites_preds_daggr$EC_FCH4_comb_median <- all_sites_preds_daggr$EC_FCH4_F_ANNOPTLM_median

USSTJ_daggr$EC_FCH4_comb_mean <- USSTJ_daggr$EC_FCH4_mean
all_sites_preds_daggr$EC_FCH4_comb_mean <- all_sites_preds_daggr$EC_FCH4_F_ANNOPTLM_mean

# EC_TS_TOP = EC_TS

USSTJ_daggr$EC_TS_TOP <- USSTJ_daggr$EC_TS_mean

USSTJ_daggr <- USSTJ_daggr %>%
  rename(
    EC_PA_F_mean   = EC_PA_mean, # not really gap-filled but needs to be in the same column as in rest of the sites
    EC_PA_F_median     = EC_PA_median
  )

setdiff(names(all_sites_preds_daggr), names(USSTJ_daggr))
setdiff(names(USSTJ_daggr), names(all_sites_preds_daggr))

cols_2 <- setdiff(names(all_sites_preds_daggr), names(USSTJ_daggr))

USSTJ_daggr[cols_2] <- NA

# combine
all_sites_daily_2 <- bind_rows(all_sites_preds_daggr, USSTJ_daggr)

##### update: fixing some US-Owc missing data (already fixed in the published dataset)

setdiff(names(all_sites_daily_2), names(USOWC_d_aggr_ECfilt))

USOWC_d_aggr_ECfilt$EC_TS_TOP <- USOWC_d_aggr_ECfilt$EC_TS_1_mean
USOWC_d_aggr_ECfilt$EC_SWC_TOP <- NA
USOWC_d_aggr_ECfilt$ch_SWC_TOP <- NA

cols_2 <- setdiff(names(all_sites_daily_2), names(USOWC_d_aggr_ECfilt))

USOWC_d_aggr_ECfilt[cols_2] <- NA

USOWC_d_aggr_ECfilt$EC_FCH4_comb_mean <- USOWC_d_aggr_ECfilt$EC_FCH4_F_ANNOPTLM_mean
USOWC_d_aggr_ECfilt$EC_FCH4_comb_median <- USOWC_d_aggr_ECfilt$EC_FCH4_F_ANNOPTLM_median

# 1. Remove US-OWC from the main df
all_sites_daily_3 <- subset(all_sites_daily_2, SITE != "US-OWC")

# 2. Ensure column names match between dataframes
USOWC_d_aggr_ECfilt <- USOWC_d_aggr_ECfilt[ , names(all_sites_daily_3)]

# 3. Add the rows from USOWC_d_aggr_ECfilt
all_sites_daily_3 <- rbind(all_sites_daily_3, USOWC_d_aggr_ECfilt)

##### update: US-LA2 fix (already fixed in the published dataset)

setdiff(names(all_sites_daily_3), names(USLA2_d_aggr_ECfilt))

cols_2 <- setdiff(names(all_sites_daily_5), names(USLA2_d_aggr_ECfilt))

USLA2_d_aggr_ECfilt[cols_2] <- NA

USLA2_d_aggr_ECfilt$EC_FCH4_comb_mean <- USLA2_d_aggr_ECfilt$EC_FCH4_F_ANNOPTLM_mean
USLA2_d_aggr_ECfilt$EC_FCH4_comb_median <- USLA2_d_aggr_ECfilt$EC_FCH4_F_ANNOPTLM_median

USLA2_d_aggr_ECfilt$ch_TS_TOP <- USLA2_d_aggr_ECfilt$ch_TS_10_mean
USLA2_d_aggr_ECfilt$EC_TS_TOP <- USLA2_d_aggr_ECfilt$EC_TS_F_mean

# 1. Remove US-OWC from the main df
all_sites_daily_3 <- subset(all_sites_daily_3, SITE != "US-LA2")

# 2. Ensure column names match between dataframes
USLA2_d_aggr_ECfilt <- USLA2_d_aggr_ECfilt[ , names(all_sites_daily_3)]

# 3. Add the rows from USLA2_d_aggr_ECfilt
all_sites_daily_5 <- rbind(all_sites_daily_3, USLA2_d_aggr_ECfilt)

all_sites_daily_5 <- read.csv("path/allsites_daily_aggr_11082025.csv")

##### UPDATED DATA: add fixed CN-HGU
# This section was run with updated CN-HGU data already included in the published datasets.
# To rebuild from source files, add the updated CN-HGU data analogously to the site-specific sections above.

setdiff(names(all_sites_daily_5), names(CNHGU_d_aggr_ECfilt))

cols_2 <- setdiff(names(all_sites_daily_5), names(CNHGU_d_aggr_ECfilt))

CNHGU_d_aggr_ECfilt[cols_2] <- NA

CNHGU_d_aggr_ECfilt$EC_FCH4_comb_mean <- CNHGU_d_aggr_ECfilt$EC_FCH4_F_ANNOPTLM_mean
CNHGU_d_aggr_ECfilt$EC_FCH4_comb_median <- CNHGU_d_aggr_ECfilt$EC_FCH4_F_ANNOPTLM_median

CNHGU_d_aggr_ECfilt$EC_TS_TOP <- CNHGU_d_aggr_ECfilt$EC_TS_F_mean
CNHGU_d_aggr_ECfilt$SITE_TS_TOP <- CNHGU_d_aggr_ECfilt$EC_TS_TOP

# 1. Remove CN-HGU from the main df
all_sites_daily_6 <- subset(all_sites_daily_5, SITE != "CN-HGU")

# 2. Ensure column names match between dataframes
CNHGU_d_aggr_ECfilt <- CNHGU_d_aggr_ECfilt[ , names(all_sites_daily_6)]

# 3. Add the rows from CNHGU_d_aggr_ECfilt
all_sites_daily_6 <- rbind(all_sites_daily_6, CNHGU_d_aggr_ECfilt)

# create a combination NEE F column 

all_sites_daily_6 <- all_sites_daily_6 %>%
  mutate(
    EC_NEE_F_comb_mean   = if_else(SITE == "US-STJ", EC_NEE_F_mean,   EC_NEE_F_ANNOPTLM_mean),
    EC_NEE_F_comb_median = if_else(SITE == "US-STJ", EC_NEE_F_median, EC_NEE_F_ANNOPTLM_median)
  )

### update: fix for SE-DEG data (already fixed in the published dataset)

setdiff(names(all_sites_daily_5), names(SEDEG_d_aggr_ECfilt))

cols_2 <- setdiff(names(all_sites_daily_4), names(SEDEG_d_aggr_ECfilt))

SEDEG_d_aggr_ECfilt[cols_2] <- NA

SEDEG_d_aggr_ECfilt$EC_FCH4_comb_mean <- SEDEG_d_aggr_ECfilt$EC_FCH4_F_ANNOPTLM_mean
SEDEG_d_aggr_ECfilt$EC_FCH4_comb_median <- SEDEG_d_aggr_ECfilt$EC_FCH4_F_ANNOPTLM_median

SEDEG_d_aggr_ECfilt$ch_TS_TOP <- SEDEG_d_aggr_ECfilt$ch_TS_2_mean
SEDEG_d_aggr_ECfilt$EC_TS_TOP <- SEDEG_d_aggr_ECfilt$EC_TS_1_mean
SEDEG_d_aggr_ECfilt$SITE_TS_TOP <- mean(c(SEDEG_d_aggr_ECfilt$EC_TS_TOP, SEDEG_d_aggr_ECfilt$ch_TS_TOP), na.rm=T)

SEDEG_d_aggr_ECfilt$EC_WTD_F <- SEDEG_d_aggr_ECfilt$EC_WTD_F_mean * 100
SEDEG_d_aggr_ECfilt$SITE_WTD <- mean(c(SEDEG_d_aggr_ECfilt$EC_WTD_F, SEDEG_d_aggr_ECfilt$ch_WTL_mean))

# 1. Remove SE-DEG from the main df
all_sites_daily_5 <- subset(all_sites_daily_5, SITE != "SE-DEG")

# 2. Ensure column names match between dataframes
SEDEG_d_aggr_ECfilt <- SEDEG_d_aggr_ECfilt[ , names(all_sites_daily_5)]

# 3. Add the rows from SEDEG_d_aggr_ECfilt
all_sites_daily_5 <- rbind(all_sites_daily_5, SEDEG_d_aggr_ECfilt)

# create a combination NEE F column 

all_sites_daily_5 <- all_sites_daily_5 %>%
  mutate(
    EC_NEE_F_comb_mean   = if_else(SITE == "US-STJ", EC_NEE_F_mean,   EC_NEE_F_ANNOPTLM_mean),
    EC_NEE_F_comb_median = if_else(SITE == "US-STJ", EC_NEE_F_median, EC_NEE_F_ANNOPTLM_median)
  )

# add EC_TS_TOP to US-STJ

all_sites_daily_5 <- all_sites_daily_5 %>%
  mutate(EC_TS_TOP = if_else(SITE == "US-STJ", EC_TS_mean, EC_TS_TOP))

# create a SITE_TS_TOP

all_sites_daily_5 <- all_sites_daily_5 %>%
  mutate(SITE_TS_TOP = case_when(
    SITE %in% c("FI-SI2", "US-HO1") ~ ch_TS_TOP,
    SITE %in% c("US-OWC", "CN-HGU", "US-STJ") ~ EC_TS_TOP,
    TRUE ~ rowMeans(cbind(ch_TS_TOP, EC_TS_TOP), na.rm = TRUE)
  ))

# create SITE_WTD

# EC_WTD is in m so convert EC_WTD to cm

all_sites_daily_5$EC_WTD_F <- ifelse(
  all_sites_daily_5$SITE == "US-STJ",
  all_sites_daily_5$EC_WTD_F_mean,
  all_sites_daily_5$EC_WTD_F_mean * 100
)

# Create the new column SITE_WTD
all_sites_daily_5 <- all_sites_daily_5 %>%
  mutate(SITE_WTD = case_when(
    SITE %in% c("FI-SI2", "US-LA1", "US-LA2") ~ ch_WTL_mean,
    SITE %in% c("US-HO1", "US-OWC", "US-UAF", "US-STJ") ~ EC_WTD_F,
    SITE == "CN-HGU" ~ NA,
    TRUE ~ rowMeans(cbind(ch_WTL_mean, EC_WTD_F), na.rm = TRUE)
  ))

# save as .csv
write.csv(all_sites_daily_6, "path/allsites_daily_aggr_28082025.csv", row.names = FALSE)

##### WEEKLY DATA SET #####

# read in datasets

CNHGU_week_aggr_ECfilt <- read.csv("path/CN_HGU_chec_weekaggr_ECfilt.csv")

SEDEG_week_aggr_ECfilt <- read.csv("path/SE_DEG_chec_weekaggr_ECfilt.csv")

USHO1_week_aggr_ECfilt <- read.csv("path/US_HO1_chec_weekaggr_ECfilt.csv")

USUAF_week_aggr_ECfilt <- read.csv("path/US_UAF_chec_weekaggr_ECfilt.csv")

FISI2_week_aggr_ECfilt <- read.csv("path/FI_SI2_chec_weekaggr_ECfilt.csv")

USLA1_week_aggr_ECfilt <- read.csv("path/US_LA1_chec_weekaggr_ECfilt.csv")

USLA2_week_aggr_ECfilt <- read.csv("path/US_LA2_chec_weekaggr_ECfilt.csv")

USLOS_week_aggr_ECfilt <- read.csv("path/US_LOS_chec_weekaggr_ECfilt.csv")

USOWC_week_aggr_ECfilt <- read.csv("path/US_OWC_chec_weekaggr_ECfilt.csv")

# harmonize datasets

USHO1_week_aggr_ECfilt <- USHO1_week_aggr_ECfilt %>% rename("ch_FCH4_nmolCH4m2s1_mean" = "ch_FCH4_umolCH4m2s1_mean",
                                                            "ch_FCH4_nmolCH4m2s1_median" = "ch_FCH4_umolCH4m2s1_median"
)

# remove the ch measurements of the same EC neasurennents

SEDEG_week_aggr_ECfilt <- subset(SEDEG_week_aggr_ECfilt, select=-c(ch_PPT_mean, ch_WS_mean, ch_P_mean, ch_VPD_mean, ch_RH_mean,
                                                                   ch_PPT_median, ch_WS_median, ch_P_median, ch_VPD_median, ch_RH_median))

# add missing columns to dataframes

# rename ch_WT to ch_WTL

FISI2_week_aggr_ECfilt <- FISI2_week_aggr_ECfilt %>% rename("ch_WTL_mean" = "ch_WT_mean",
                                                            "ch_WTL_median" = "ch_WT_median"
)

# remove ch_TS_5 and 15

FISI2_week_aggr_ECfilt <- subset(FISI2_week_aggr_ECfilt, select=-c(ch_TS_15_mean,  ch_TS_15_median))
USLA1_week_aggr_ECfilt <- subset(USLA1_week_aggr_ECfilt, select=-c(ch_TS_5_mean, ch_TS_5_median,  ch_TS_15_mean,  ch_TS_15_median))
USLA2_week_aggr_ECfilt <- subset(USLA2_week_aggr_ECfilt, select=-c(ch_TS_5_mean, ch_TS_5_median, ch_TS_15_mean, ch_TS_15_median))

cols_1 <- c("ch_TS_2_mean", "ch_TS_5_mean", "ch_TS_5_median", "ch_TS_10_mean", "ch_TS_20_mean", "ch_TS_30_mean", "ch_TS_40_mean", 
            "ch_CH4_Ta_amb_mean", 
            "ch_TS_2_median",
            "EC_SWC_1_mean", "EC_PA_mean", "EC_TS_1_mean", "EC_TS_2_mean", 
            "EC_TS_3_mean", "EC_TS_4_mean", "EC_TS_5_mean", "EC_TS_6_mean", "EC_TS_7_mean", "EC_TS_8_mean", "EC_TS_9_mean", 
            "EC_WTD_F_mean", "EC_PPFD_OUT_mean",
            "EC_SWC_2_mean", "EC_SWC_3_mean", "EC_PPFD_OUT_F_mean", 
            "EC_SWC_1_median", "EC_PA_median", "EC_TS_1_median", "EC_TS_2_median", 
            "EC_TS_3_median", "EC_TS_4_median", "EC_TS_5_median", "EC_TS_6_median", "EC_TS_7_median", "EC_TS_8_median", "EC_TS_9_median", 
            "EC_WTD_F_median", "EC_PPFD_OUT_median",
            "EC_SWC_2_median", "EC_SWC_3_median", "EC_PPFD_OUT_F_median", 
            "ch_SM_X_mean", 
            "ch_TS_X_mean", "ch_SM_X_median", "ch_TS_X_median", 
            "ch_VWC_mean",
            "ch_TS_20_median", "ch_TS_30_median", "ch_TS_40_median", "ch_VWC_median", 
            "ch_AerLAI_mean", "ch_AerLAI_median", 
            "salinity_mean",  "salinity_median")

CNHGU_week_aggr_ECfilt[cols_1] <- NA

# now CNHGU is a df that includes all columns --> compare the rest to this

setdiff(names(CNHGU_week_aggr_ECfilt), names(SEDEG_week_aggr_ECfilt))

cols_2 <- setdiff(names(CNHGU_week_aggr_ECfilt), names(SEDEG_week_aggr_ECfilt))

SEDEG_week_aggr_ECfilt[cols_2] <- NA

#---

setdiff(names(CNHGU_week_aggr_ECfilt), names(USHO1_week_aggr_ECfilt))

cols_3 <- setdiff(names(CNHGU_week_aggr_ECfilt), names(USHO1_week_aggr_ECfilt))

USHO1_week_aggr_ECfilt[cols_3] <- NA

#---

setdiff(names(CNHGU_week_aggr_ECfilt), names(USUAF_week_aggr_ECfilt))

cols_4 <- setdiff(names(CNHGU_week_aggr_ECfilt), names(USUAF_week_aggr_ECfilt))

USUAF_week_aggr_ECfilt[cols_4] <- NA

#---

setdiff(names(CNHGU_week_aggr_ECfilt), names(FISI2_week_aggr_ECfilt))

cols_5 <- setdiff(names(CNHGU_week_aggr_ECfilt), names(FISI2_week_aggr_ECfilt))

FISI2_week_aggr_ECfilt[cols_5] <- NA

#---

setdiff(names(CNHGU_week_aggr_ECfilt), names(USLA1_week_aggr_ECfilt))

cols_6 <- setdiff(names(CNHGU_week_aggr_ECfilt), names(USLA1_week_aggr_ECfilt))

USLA1_week_aggr_ECfilt[cols_6] <- NA

#---

setdiff(names(CNHGU_week_aggr_ECfilt), names(USLA2_week_aggr_ECfilt))

cols_7 <- setdiff(names(CNHGU_week_aggr_ECfilt), names(USLA2_week_aggr_ECfilt))

USLA2_week_aggr_ECfilt[cols_7] <- NA

#---

setdiff(names(CNHGU_week_aggr_ECfilt), names(USLOS_week_aggr_ECfilt))

cols_8 <- setdiff(names(CNHGU_week_aggr_ECfilt), names(USLOS_week_aggr_ECfilt))

USLOS_week_aggr_ECfilt[cols_8] <- NA

#---

setdiff(names(CNHGU_week_aggr_ECfilt), names(USOWC_week_aggr_ECfilt))

cols_9 <- setdiff(names(CNHGU_week_aggr_ECfilt), names(USOWC_week_aggr_ECfilt))

USOWC_week_aggr_ECfilt[cols_9] <- NA

# combine!

all_sites_week_ECfilt <- bind_rows(CNHGU_week_aggr_ECfilt, FISI2_week_aggr_ECfilt, SEDEG_week_aggr_ECfilt, USHO1_week_aggr_ECfilt,
                                   USLA1_week_aggr_ECfilt, USLA2_week_aggr_ECfilt, USLOS_week_aggr_ECfilt, USOWC_week_aggr_ECfilt, USUAF_week_aggr_ECfilt)

# remove UNC

all_sites_week_ECfilt <- subset(all_sites_week_ECfilt, select = -c(EC_FCH4_F_ANNOPTLM_UNC_mean, EC_FCH4_F_ANNOPTLM_UNC_median,
                                                                   ECCH_diff_mean, ECCH_diff_median))

# combine TS info into one column

# Create the new columns based on the "SITE" group
all_sites_week_ECfilt <- all_sites_week_ECfilt %>%
  mutate(
    ch_TS_TOP = case_when(
      SITE == 'FI-SI2' ~ ch_TS_5_mean,
      SITE == 'CN-HGU' ~ NA,
      SITE == 'SE-DEG' ~ ch_TS_2_mean,
      SITE == 'US-HO1' ~ ch_TS_X_mean,
      SITE == 'US-LA1' ~ ch_TS_10_mean,
      SITE == 'US-LA2' ~ ch_TS_10_mean,
      SITE == 'US-LOS' ~ ch_TS_X_mean,
      SITE == 'US-OWC' ~ NA,
      SITE == 'US-UAF' ~ ch_TS_2_mean,
      TRUE ~ NA_real_
    ),
    EC_TS_TOP = case_when(
      SITE == 'FI-SI2' ~ EC_TS_1_mean,
      SITE == 'CN-HGU' ~ EC_TS_F_mean,
      SITE == 'SE-DEG' ~ EC_TS_1_mean,
      SITE == 'US-HO1' ~ EC_TS_1_mean,
      SITE == 'US-LA1' ~ EC_TS_F_mean,
      SITE == 'US-LA2' ~ EC_TS_F_mean,
      SITE == 'US-LOS' ~ EC_TS_1_mean,
      SITE == 'US-OWC' ~ EC_TS_1_mean,
      SITE == 'US-UAF' ~ EC_TS_1_mean,
      TRUE ~ NA_real_
    )
  )

# create a SITE_TS_TOP

all_sites_week_ECfilt <- all_sites_week_ECfilt %>%
  mutate(SITE_TS_TOP = case_when(
    SITE %in% c("FI-SI2", "US-HO1") ~ ch_TS_TOP,
    SITE %in% c("US-OWC", "CN-HGU") ~ EC_TS_TOP,
    TRUE ~ rowMeans(cbind(ch_TS_TOP, EC_TS_TOP), na.rm = TRUE)
  ))

# do the same for EC_SWC 

# Create the new columns based on the "SITE" group
all_sites_week_ECfilt <- all_sites_week_ECfilt %>%
  mutate(
    EC_SWC_TOP = case_when(
      SITE == 'FI-SI2' ~ NA,
      SITE == 'CN-HGU' ~ EC_SWC_F_mean,
      SITE == 'SE-DEG' ~ EC_SWC_1_mean,
      SITE == 'US-HO1' ~ NA,
      SITE == 'US-LA1' ~ NA,
      SITE == 'US-LA2' ~ NA,
      SITE == 'US-LOS' ~ NA,
      SITE == 'US-OWC' ~ NA,
      SITE == 'US-UAF' ~ EC_SWC_1_mean,
      TRUE ~ NA_real_
    ),
    ch_SWC_TOP = case_when(
      SITE == 'FI-SI2' ~ NA,
      SITE == 'CN-HGU' ~ NA,
      SITE == 'SE-DEG' ~ NA,
      SITE == 'US-HO1' ~ ch_SM_X_mean,
      SITE == 'US-LA1' ~ NA,
      SITE == 'US-LA2' ~ NA,
      SITE == 'US-LOS' ~ ch_SM_X_mean,
      SITE == 'US-OWC' ~ NA,
      SITE == 'US-UAF' ~ ch_VWC_mean,
      TRUE ~ NA_real_
    )
  )

all_sites_preds_dlcaggr <- all_sites_preds_dlcaggr %>%
  mutate(
    EC_SWC_TOP = case_when(
      SITE == 'FI-SI2' ~ NA,
      SITE == 'US-HO1' ~ NA,
      SITE == 'US-LOS' ~ NA,
      SITE == 'US-OWC' ~ NA,
      SITE == 'US-UAF' ~ EC_SWC_1_mean,
      TRUE ~ NA_real_
    ),
    ch_SWC_TOP = case_when(
      SITE == 'FI-SI2' ~ NA,
      SITE == 'US-HO1' ~ ch_SM_X_mean,
      SITE == 'US-LOS' ~ ch_SM_X_mean,
      SITE == 'US-OWC' ~ NA,
      SITE == 'US-UAF' ~ ch_VWC_mean,
      TRUE ~ NA_real_
    )
  )

# create SITE_SWC

all_sites_week_ECfilt <- all_sites_week_ECfilt %>%
  mutate(SITE_SWC = case_when(
    SITE %in% c("US-HO1") ~ ch_SWC_TOP,
    TRUE ~ rowMeans(cbind(ch_SWC_TOP, EC_SWC_TOP), na.rm = TRUE)
  ))

all_sites_preds_daggr <- all_sites_preds_daggr %>%
  mutate(SITE_SWC = case_when(
    SITE %in% c("US-HO1") ~ ch_SWC_TOP,
    TRUE ~ rowMeans(cbind(ch_SWC_TOP, EC_SWC_TOP), na.rm = TRUE)
  ))

all_sites_preds_dlcaggr <- all_sites_preds_dlcaggr %>%
  mutate(SITE_SWC = case_when(
    SITE %in% c("US-HO1") ~ ch_SWC_TOP,
    TRUE ~ rowMeans(cbind(ch_SWC_TOP, EC_SWC_TOP), na.rm = TRUE)
  ))

# create SITE_WTD

# EC_WTD is in m so convert EC_WTD to cm

all_sites_week_ECfilt$EC_WTD_F <- all_sites_week_ECfilt$EC_WTD_F_mean * 100

# Create the new column SITE_WTD
all_sites_week_ECfilt <- all_sites_week_ECfilt %>%
  mutate(SITE_WTD = case_when(
    SITE %in% c("FI-SI2", "US-LA1", "US-LA2") ~ ch_WTL_mean,
    SITE %in% c("US-HO1", "US-OWC", "US-UAF") ~ EC_WTD_F,
    SITE == "CN-HGU" ~ NA,
    TRUE ~ rowMeans(cbind(ch_WTL_mean, EC_WTD_F), na.rm = TRUE)
  ))

# remove the extra columns

all_sites_week_ECfilt <- subset(all_sites_week_ECfilt, select = -c(EC_TS_1_mean, EC_TS_2_mean, EC_TS_3_mean, EC_TS_4_mean, EC_TS_5_mean,
                                                                   EC_TS_6_mean, EC_TS_7_mean, EC_TS_8_mean, EC_TS_9_mean,
                                                                   ch_SM_X_mean, ch_TS_X_mean, ch_TS_2_mean, ch_TS_10_mean, 
                                                                   ch_TS_20_mean, ch_TS_30_mean, ch_TS_40_mean,
                                                                   ch_VWC_mean, EC_SWC_1_mean, EC_SWC_2_mean, EC_SWC_3_mean))

all_sites_week_ECfilt <- subset(all_sites_week_ECfilt, select = -c(EC_TS_1_median, EC_TS_2_median, EC_TS_3_median, EC_TS_4_median, EC_TS_5_median,
                                                                   EC_TS_6_median, EC_TS_7_median, EC_TS_8_median, EC_TS_9_median,
                                                                   ch_SM_X_median, ch_TS_X_median, ch_TS_2_median, ch_TS_10_median, 
                                                                   ch_TS_20_median, ch_TS_30_median, ch_TS_40_median,
                                                                   ch_VWC_median, EC_SWC_1_median, EC_SWC_2_median, EC_SWC_3_median))

# calculate ECCH_diff
all_sites_week_ECfilt$ECCH_diff <- all_sites_week_ECfilt$EC_FCH4_F_ANNOPTLM_median - all_sites_week_ECfilt$ch_FCH4_nmolCH4m2s1_median

# remove NA

all_sites_week_ECfilt <- all_sites_week_ECfilt %>% drop_na(ECCH_diff)
median(all_sites_week_ECfilt$ECCH_diff, na.rm = T)

###### UPDATED DATA: ADD US-STJ
# This section was run with updated US-STJ data already included in the published datasets.
# To rebuild from source files, add the updated US-STJ data analogously to the site-specific sections above.

USSTJ_week <- read.csv("path/US_STJ_chec_weekaggr.csv")

all_week_aggr <- all_sites_week_ECfilt

# change some column names

USSTJ_week <- USSTJ_week %>%
  rename(
    EC_NEE_F_median   = EC_NEE_F_MDS_median,
    EC_NEE_F_mean     = EC_NEE_F_MDS_mean,
    EC_FCH4_median = EC_FCH4_nmolCH4m2s1_median,
    EC_FCH4_mean = EC_FCH4_nmolCH4m2s1_mean
  )

# remove TA and EC_NEE (gap-filled NEE better)

USSTJ_week <- subset(USSTJ_week, select = -c(EC_TA_mean, EC_TA_median, EC_NEE_mean, EC_NEE_median))

USSTJ_week <- USSTJ_week %>%
  rename(
    EC_WS_F_mean   = EC_WS_mean, # not really gap-filled but needs to be in the same column as in rest of the sites
    EC_WS_F_median     = EC_WS_median, # --"--
    EC_VPD_F_mean = EC_VPD_mean, #--"--
    EC_VPD_F_median = EC_VPD_median # --"--
  )

USSTJ_week <- USSTJ_week %>%
  rename(
    v.wind.F_mean   = v.wind_mean, # not really gap-filled but needs to be in the same column as in rest of the sites
    u.wind.F_mean     = u.wind_mean
  )

USSTJ_week <- USSTJ_week %>%
  rename(
    v.wind.F_median   = v.wind_median, # not really gap-filled but needs to be in the same column as in rest of the sites
    u.wind.F_median     = u.wind_median
  )

### need to combine FCH4 values from USSTJ to ANNOPTLM values from other sites
# --> create new column for this

USSTJ_week$EC_FCH4_comb_median <- USSTJ_week$EC_FCH4_median
all_week_aggr$EC_FCH4_comb_median <- all_week_aggr$EC_FCH4_F_ANNOPTLM_median

USSTJ_week$EC_FCH4_comb_mean <- USSTJ_week$EC_FCH4_mean
all_week_aggr$EC_FCH4_comb_mean <- all_week_aggr$EC_FCH4_F_ANNOPTLM_mean

# EC_TS_TOP = EC_TS

USSTJ_week$EC_TS_TOP <- USSTJ_week$EC_TS_mean

USSTJ_week <- USSTJ_week %>%
  rename(
    EC_PA_F_mean   = EC_PA_mean, # not really gap-filled but needs to be in the same column as in rest of the sites
    EC_PA_F_median     = EC_PA_median
  )

USSTJ_week <- subset(USSTJ_week, select = -c(ECCH_diff_mean, ECCH_diff_median, EC_TS_median))

USSTJ_week$ECCH_diff <- USSTJ_week$EC_FCH4_comb_median - USSTJ_week$ch_FCH4_nmolCH4m2s1_median

USSTJ_week$SITE_TS_TOP <- USSTJ_week$EC_TS_TOP

USSTJ_week <- USSTJ_week %>%
  rename(
    EC_WTD_F   = EC_WTD_F_mean
  )

USSTJ_week$SITE_WTD <- USSTJ_week$EC_WTD_F

setdiff(names(all_week_aggr), names(USSTJ_week))
setdiff(names(USSTJ_week), names(all_week_aggr))

cols_2 <- setdiff(names(all_week_aggr), names(USSTJ_week))

USSTJ_week[cols_2] <- NA

# combine
all_sites_week_2 <- bind_rows(all_week_aggr, USSTJ_week)

#### update 12.08.2025: fixing missing US-Owc data (already fixed in the published dataset)

USOWC_week_aggr_ECfilt$EC_TS_TOP <- USOWC_week_aggr_ECfilt$EC_TS_1_mean

USOWC_week_aggr_ECfilt$EC_FCH4_comb_median <- USOWC_week_aggr_ECfilt$EC_FCH4_F_ANNOPTLM_median
USOWC_week_aggr_ECfilt$EC_FCH4_comb_mean <- USOWC_week_aggr_ECfilt$EC_FCH4_F_ANNOPTLM_mean

USOWC_week_aggr_ECfilt$EC_WTD_F <- USOWC_week_aggr_ECfilt$EC_WTD_F_mean * 100

USOWC_week_aggr_ECfilt$SITE_WTD <- USOWC_week_aggr_ECfilt$EC_WTD_F

USOWC_week_aggr_ECfilt$SITE_TS_TOP <- USOWC_week_aggr_ECfilt$EC_TS_TOP

USOWC_week_aggr_ECfilt$ECCH_diff <- USOWC_week_aggr_ECfilt$EC_FCH4_comb_median - USOWC_week_aggr_ECfilt$ch_FCH4_nmolCH4m2s1_median

setdiff(names(all_sites_week_2), names(USOWC_week_aggr_ECfilt))

cols_9 <- setdiff(names(all_sites_week_2), names(USOWC_week_aggr_ECfilt))

USOWC_week_aggr_ECfilt[cols_9] <- NA

# 1. Remove US-OWC from the main df
all_sites_week_3 <- subset(all_sites_week_2, SITE != "US-OWC")

# 2. Ensure column names match between dataframes
USOWC_week_aggr_ECfilt <- USOWC_week_aggr_ECfilt[ , names(all_sites_week_3)]

# 3. Add the rows from USOWC_d_aggr_ECfilt
all_sites_week_3 <- rbind(all_sites_week_3, USOWC_week_aggr_ECfilt)

# create a combination NEE F column 

all_sites_week_4 <- all_sites_week_4 %>%
  mutate(
    EC_NEE_F_comb_mean   = if_else(SITE == "US-STJ", EC_NEE_F_mean,   EC_NEE_F_ANNOPTLM_mean),
    EC_NEE_F_comb_median = if_else(SITE == "US-STJ", EC_NEE_F_median, EC_NEE_F_ANNOPTLM_median)
  )

##### UPDATED DATA: add SE-DEG
# This section was run with updated SE-DEG data already included in the published datasets.
# To rebuild from source files, add the updated SE-DEG data analogously to the site-specific sections above.

setdiff(names(all_sites_week_3), names(SEDEG_week_aggr_ECfilt))

cols_2 <- setdiff(names(all_sites_week_3), names(SEDEG_week_aggr_ECfilt))

SEDEG_week_aggr_ECfilt[cols_2] <- NA

SEDEG_week_aggr_ECfilt$EC_FCH4_comb_mean <- SEDEG_week_aggr_ECfilt$EC_FCH4_F_ANNOPTLM_mean
SEDEG_week_aggr_ECfilt$EC_FCH4_comb_median <- SEDEG_week_aggr_ECfilt$EC_FCH4_F_ANNOPTLM_median

SEDEG_week_aggr_ECfilt$ch_TS_TOP <- SEDEG_week_aggr_ECfilt$ch_TS_2_mean
SEDEG_week_aggr_ECfilt$EC_TS_TOP <- SEDEG_week_aggr_ECfilt$EC_TS_1_mean
SEDEG_week_aggr_ECfilt$SITE_TS_TOP <- mean(c(SEDEG_week_aggr_ECfilt$EC_TS_TOP, SEDEG_week_aggr_ECfilt$ch_TS_TOP), na.rm=T)

SEDEG_week_aggr_ECfilt$EC_WTD_F <- SEDEG_week_aggr_ECfilt$EC_WTD_F_mean * 100
SEDEG_week_aggr_ECfilt$SITE_WTD <- mean(c(SEDEG_week_aggr_ECfilt$EC_WTD_F, SEDEG_week_aggr_ECfilt$ch_WTL_mean))

# 1. Remove SE-DEG from the main df
all_sites_week_4 <- subset(all_sites_week_3, SITE != "SE-DEG")

# 2. Ensure column names match between dataframes
SEDEG_week_aggr_ECfilt <- SEDEG_week_aggr_ECfilt[ , names(all_sites_week_4)]

# 3. Add the rows from SEDEG_d_aggr_ECfilt
all_sites_week_4 <- rbind(all_sites_week_4, SEDEG_week_aggr_ECfilt)

write.csv(all_sites_week_4, "path/allsites_ECfilt_weekaggr_12082025.csv", row.names = FALSE)

### UPDATED DATA: add fixed CN-HGU data
# This section was run later with updated CN-HGU data already included in the published datasets.

all_sites_week_4 <- read.csv("path/allsites_ECfilt_weekaggr_12082025.csv")

##### UPDATED DATA: add fixed CN-HGU
# This section was run with updated CN-HGU data already included in the published datasets.
# To rebuild from source files, add the updated CN-HGU data analogously to the site-specific sections above.

setdiff(names(all_sites_week_4), names(CNHGU_week_aggr_ECfilt))

cols_2 <- setdiff(names(all_sites_week_4), names(CNHGU_week_aggr_ECfilt))

CNHGU_week_aggr_ECfilt[cols_2] <- NA

CNHGU_week_aggr_ECfilt$EC_FCH4_comb_mean <- CNHGU_week_aggr_ECfilt$EC_FCH4_F_ANNOPTLM_mean
CNHGU_week_aggr_ECfilt$EC_FCH4_comb_median <- CNHGU_week_aggr_ECfilt$EC_FCH4_F_ANNOPTLM_median

CNHGU_week_aggr_ECfilt$EC_TS_TOP <- CNHGU_week_aggr_ECfilt$EC_TS_F_mean
CNHGU_week_aggr_ECfilt$SITE_TS_TOP <- CNHGU_week_aggr_ECfilt$EC_TS_TOP

CNHGU_week_aggr_ECfilt$EC_NEE_F_comb_mean <- CNHGU_week_aggr_ECfilt$EC_NEE_F_ANNOPTLM_mean
CNHGU_week_aggr_ECfilt$EC_NEE_F_comb_median <- CNHGU_week_aggr_ECfilt$EC_NEE_F_ANNOPTLM_median

# 1. Remove CN-HGU from the main df
all_sites_week_5 <- subset(all_sites_week_4, SITE != "CN-HGU")

# 2. Ensure column names match between dataframes
CNHGU_week_aggr_ECfilt <- CNHGU_week_aggr_ECfilt[ , names(all_sites_week_5)]

# 3. Add the rows from CNHGU_d_aggr_ECfilt
all_sites_week_5 <- rbind(all_sites_week_5, CNHGU_week_aggr_ECfilt)

write.csv(all_sites_week_5, "path/allsites_ECfilt_weekaggr_28082025.csv", row.names = FALSE)

##### MONTHLY DATA SETS #####

# read in datasets

CNHGU_month_aggr_ECfilt <- read.csv("path/CN_HGU_chec_monthnaggr_ECfilt.csv")

SEDEG_month_aggr_ECfilt <- read.csv("path/SE_DEG_chec_monthaggr_ECfilt.csv")

USHO1_month_aggr_ECfilt <- read.csv("path/US_HO1_chec_monthaggr_ECfilt.csv")

USUAF_month_aggr_ECfilt <- read.csv("path/US_UAF_chec_monthaggr_ECfilt.csv")

FISI2_month_aggr_ECfilt <- read.csv("path/FI_SI2_chec_monthaggr_ECfilt.csv")

USLA1_month_aggr_ECfilt <- read.csv("path/US_LA1_chec_monthaggr_ECfilt.csv")

USLA2_month_aggr_ECfilt <- read.csv("path/US_LA2_chec_monthaggr_ECfilt.csv")

USLOS_month_aggr_ECfilt <- read.csv("path/US_LOS_chec_monthaggr_ECfilt.csv")

USOWC_month_aggr_ECfilt <- read.csv("path/US_OWC_chec_monthaggr_ECfilt.csv")

# harmonize datasets

USHO1_month_aggr_ECfilt <- USHO1_month_aggr_ECfilt %>% rename("ch_FCH4_nmolCH4m2s1_mean" = "ch_FCH4_umolCH4m2s1_mean",
                                                              "ch_FCH4_nmolCH4m2s1_median" = "ch_FCH4_umolCH4m2s1_median"
)

# remove the ch measurements of the same EC neasurennents

SEDEG_month_aggr_ECfilt <- subset(SEDEG_month_aggr_ECfilt, select=-c(ch_PPT_mean, ch_WS_mean, ch_P_mean, ch_VPD_mean, ch_RH_mean,
                                                                     ch_PPT_median, ch_WS_median, ch_P_median, ch_VPD_median, ch_RH_median))

# add missing columns to dataframes

# rename ch_WT to ch_WTL

FISI2_month_aggr_ECfilt <- FISI2_month_aggr_ECfilt %>% rename("ch_WTL_mean" = "ch_WT_mean",
                                                              "ch_WTL_median" = "ch_WT_median"
)

# remove ch_TS_5 and 15

FISI2_month_aggr_ECfilt <- subset(FISI2_month_aggr_ECfilt, select=-c(ch_TS_15_mean,  ch_TS_15_median))
USLA1_month_aggr_ECfilt <- subset(USLA1_month_aggr_ECfilt, select=-c(ch_TS_5_mean, ch_TS_5_median,  ch_TS_15_mean,  ch_TS_15_median))
USLA2_month_aggr_ECfilt <- subset(USLA2_month_aggr_ECfilt, select=-c(ch_TS_5_mean, ch_TS_5_median, ch_TS_15_mean, ch_TS_15_median))

cols_1 <- c("ch_TS_2_mean", "ch_TS_5_mean", "ch_TS_5_median", "ch_TS_10_mean", "ch_TS_20_mean", "ch_TS_30_mean", "ch_TS_40_mean", 
            "ch_CH4_Ta_amb_mean", 
            "ch_TS_2_median",
            "EC_SWC_1_mean", "EC_PA_mean", "EC_TS_1_mean", "EC_TS_2_mean", 
            "EC_TS_3_mean", "EC_TS_4_mean", "EC_TS_5_mean", "EC_TS_6_mean", "EC_TS_7_mean", "EC_TS_8_mean", "EC_TS_9_mean", 
            "EC_WTD_F_mean", "EC_PPFD_OUT_mean",
            "EC_SWC_2_mean", "EC_SWC_3_mean", "EC_PPFD_OUT_F_mean", 
            "EC_SWC_1_median", "EC_PA_median", "EC_TS_1_median", "EC_TS_2_median", 
            "EC_TS_3_median", "EC_TS_4_median", "EC_TS_5_median", "EC_TS_6_median", "EC_TS_7_median", "EC_TS_8_median", "EC_TS_9_median", 
            "EC_WTD_F_median", "EC_PPFD_OUT_median",
            "EC_SWC_2_median", "EC_SWC_3_median", "EC_PPFD_OUT_F_median", 
            "ch_SM_X_mean", 
            "ch_TS_X_mean", "ch_SM_X_median", "ch_TS_X_median", 
            "ch_VWC_mean",
            "ch_TS_20_median", "ch_TS_30_median", "ch_TS_40_median", "ch_VWC_median", 
            "ch_AerLAI_mean", "ch_AerLAI_median", 
            "salinity_mean",  "salinity_median")

CNHGU_month_aggr_ECfilt[cols_1] <- NA

# now CNHGU is a df that includes all columns --> compare the rest to this

setdiff(names(CNHGU_month_aggr_ECfilt), names(SEDEG_month_aggr_ECfilt))

cols_2 <- setdiff(names(CNHGU_month_aggr_ECfilt), names(SEDEG_month_aggr_ECfilt))

SEDEG_month_aggr_ECfilt[cols_2] <- NA

#---

setdiff(names(CNHGU_month_aggr_ECfilt), names(USHO1_month_aggr_ECfilt))

cols_3 <- setdiff(names(CNHGU_month_aggr_ECfilt), names(USHO1_month_aggr_ECfilt))

USHO1_month_aggr_ECfilt[cols_3] <- NA

#---

setdiff(names(CNHGU_month_aggr_ECfilt), names(USUAF_month_aggr_ECfilt))

cols_4 <- setdiff(names(CNHGU_month_aggr_ECfilt), names(USUAF_month_aggr_ECfilt))

USUAF_month_aggr_ECfilt[cols_4] <- NA

#---

setdiff(names(CNHGU_month_aggr_ECfilt), names(FISI2_month_aggr_ECfilt))

cols_5 <- setdiff(names(CNHGU_month_aggr_ECfilt), names(FISI2_month_aggr_ECfilt))

FISI2_month_aggr_ECfilt[cols_5] <- NA

#---

setdiff(names(CNHGU_month_aggr_ECfilt), names(USLA1_month_aggr_ECfilt))

cols_6 <- setdiff(names(CNHGU_month_aggr_ECfilt), names(USLA1_month_aggr_ECfilt))

USLA1_month_aggr_ECfilt[cols_6] <- NA

#---

setdiff(names(CNHGU_month_aggr_ECfilt), names(USLA2_month_aggr_ECfilt))

cols_7 <- setdiff(names(CNHGU_month_aggr_ECfilt), names(USLA2_month_aggr_ECfilt))

USLA2_month_aggr_ECfilt[cols_7] <- NA

#---

setdiff(names(CNHGU_month_aggr_ECfilt), names(USLOS_month_aggr_ECfilt))

cols_8 <- setdiff(names(CNHGU_month_aggr_ECfilt), names(USLOS_month_aggr_ECfilt))

USLOS_month_aggr_ECfilt[cols_8] <- NA

#---

setdiff(names(CNHGU_month_aggr_ECfilt), names(USOWC_month_aggr_ECfilt))

cols_9 <- setdiff(names(CNHGU_month_aggr_ECfilt), names(USOWC_month_aggr_ECfilt))

USOWC_month_aggr_ECfilt[cols_9] <- NA

# combine!

all_sites_month_ECfilt <- bind_rows(CNHGU_month_aggr_ECfilt, FISI2_month_aggr_ECfilt, SEDEG_month_aggr_ECfilt, USHO1_month_aggr_ECfilt,
                                    USLA1_month_aggr_ECfilt, USLA2_month_aggr_ECfilt, USLOS_month_aggr_ECfilt, USOWC_month_aggr_ECfilt, USUAF_month_aggr_ECfilt)

# CN-HGU accidentally has timestamp instead of year etc

all_sites_month_ECfilt <- all_sites_month_ECfilt %>%
  # Move values from TIMESTAMP to Year, and month.CNHGU_ECfilt.TIMESTAMP_START. to Month when SITE == "CN-HGU"
  mutate(Year = ifelse(SITE == "CN-HGU", TIMESTAMP, Year),
         Month = ifelse(SITE == "CN-HGU", `month.CNHGU_ECfilt.TIMESTAMP_START.`, Month)) %>%
  # Remove the TIMESTAMP and month.CNHGU_ECfilt.TIMESTAMP_START. columns
  dplyr::select(-TIMESTAMP, -`month.CNHGU_ECfilt.TIMESTAMP_START.`)

all_sites_month_ECfilt <- all_sites_month_ECfilt %>% relocate(c(Year, Month), .after = SITE)

# combine TS info into one column

# Create the new columns based on the "SITE" group
all_sites_month_ECfilt <- all_sites_month_ECfilt %>%
  mutate(
    ch_TS_TOP = case_when(
      SITE == 'FI-SI2' ~ ch_TS_5_mean,
      SITE == 'CN-HGU' ~ NA,
      SITE == 'SE-DEG' ~ ch_TS_2_mean,
      SITE == 'US-HO1' ~ ch_TS_X_mean,
      SITE == 'US-LA1' ~ ch_TS_10_mean,
      SITE == 'US-LA2' ~ ch_TS_10_mean,
      SITE == 'US-LOS' ~ ch_TS_X_mean,
      SITE == 'US-OWC' ~ NA,
      SITE == 'US-UAF' ~ ch_TS_2_mean,
      TRUE ~ NA_real_
    ),
    EC_TS_TOP = case_when(
      SITE == 'FI-SI2' ~ EC_TS_1_mean,
      SITE == 'CN-HGU' ~ EC_TS_F_mean,
      SITE == 'SE-DEG' ~ EC_TS_1_mean,
      SITE == 'US-HO1' ~ EC_TS_1_mean,
      SITE == 'US-LA1' ~ EC_TS_F_mean,
      SITE == 'US-LA2' ~ EC_TS_F_mean,
      SITE == 'US-LOS' ~ EC_TS_1_mean,
      SITE == 'US-OWC' ~ EC_TS_1_mean,
      SITE == 'US-UAF' ~ EC_TS_1_mean,
      TRUE ~ NA_real_
    )
  )

# create a SITE_TS_TOP

all_sites_month_ECfilt <- all_sites_month_ECfilt %>%
  mutate(SITE_TS_TOP = case_when(
    SITE %in% c("FI-SI2", "US-HO1") ~ ch_TS_TOP,
    SITE %in% c("US-OWC", "CN-HGU") ~ EC_TS_TOP,
    TRUE ~ rowMeans(cbind(ch_TS_TOP, EC_TS_TOP), na.rm = TRUE)
  ))

# do the same for EC_SWC

# Create the new columns based on the "SITE" group
all_sites_month_ECfilt <- all_sites_month_ECfilt %>%
  mutate(
    EC_SWC_TOP = case_when(
      SITE == 'FI-SI2' ~ NA,
      SITE == 'CN-HGU' ~ EC_SWC_F_mean,
      SITE == 'SE-DEG' ~ EC_SWC_1_mean,
      SITE == 'US-HO1' ~ NA,
      SITE == 'US-LA1' ~ NA,
      SITE == 'US-LA2' ~ NA,
      SITE == 'US-LOS' ~ NA,
      SITE == 'US-OWC' ~ NA,
      SITE == 'US-UAF' ~ EC_SWC_1_mean,
      TRUE ~ NA_real_
    ),
    ch_SWC_TOP = case_when(
      SITE == 'FI-SI2' ~ NA,
      SITE == 'CN-HGU' ~ NA,
      SITE == 'SE-DEG' ~ NA,
      SITE == 'US-HO1' ~ ch_SM_X_mean,
      SITE == 'US-LA1' ~ NA,
      SITE == 'US-LA2' ~ NA,
      SITE == 'US-LOS' ~ ch_SM_X_mean,
      SITE == 'US-OWC' ~ NA,
      SITE == 'US-UAF' ~ ch_VWC_mean,
      TRUE ~ NA_real_
    )
  )

# create SITE_SWC

all_sites_month_ECfilt <- all_sites_month_ECfilt %>%
  mutate(SITE_SWC = case_when(
    SITE %in% c("US-HO1") ~ ch_SWC_TOP,
    TRUE ~ rowMeans(cbind(ch_SWC_TOP, EC_SWC_TOP), na.rm = TRUE)
  ))

# create SITE_WTD

# EC_WTD is in m so convert EC_WTD to cm

all_sites_month_ECfilt$EC_WTD_F <- all_sites_month_ECfilt$EC_WTD_F_mean * 100

# Create the new column SITE_WTD
all_sites_month_ECfilt <- all_sites_month_ECfilt %>%
  mutate(SITE_WTD = case_when(
    SITE %in% c("FI-SI2", "US-LA1", "US-LA2") ~ ch_WTL_mean,
    SITE %in% c("US-HO1", "US-OWC", "US-UAF") ~ EC_WTD_F,
    SITE == "CN-HGU" ~ NA,
    TRUE ~ rowMeans(cbind(ch_WTL_mean, EC_WTD_F), na.rm = TRUE)
  ))

# remove the extra columns

all_sites_month_ECfilt <- subset(all_sites_month_ECfilt, select = -c(EC_TS_1_mean, EC_TS_2_mean, EC_TS_3_mean, EC_TS_4_mean, EC_TS_5_mean,
                                                                     EC_TS_6_mean, EC_TS_7_mean, EC_TS_8_mean, EC_TS_9_mean,
                                                                     ch_SM_X_mean, ch_TS_X_mean, ch_TS_2_mean, ch_TS_10_mean, 
                                                                     ch_TS_20_mean, ch_TS_30_mean, ch_TS_40_mean,
                                                                     ch_VWC_mean, EC_SWC_1_mean, EC_SWC_2_mean, EC_SWC_3_mean))

all_sites_month_ECfilt <- subset(all_sites_month_ECfilt, select = -c(EC_TS_1_median, EC_TS_2_median, EC_TS_3_median, EC_TS_4_median, EC_TS_5_median,
                                                                     EC_TS_6_median, EC_TS_7_median, EC_TS_8_median, EC_TS_9_median,
                                                                     ch_SM_X_median, ch_TS_X_median, ch_TS_2_median, ch_TS_10_median, 
                                                                     ch_TS_20_median, ch_TS_30_median, ch_TS_40_median,
                                                                     ch_VWC_median, EC_SWC_1_median, EC_SWC_2_median, EC_SWC_3_median))

# calculate ECCH_diff
all_sites_month_ECfilt$ECCH_diff <- all_sites_month_ECfilt$EC_FCH4_F_ANNOPTLM_median - all_sites_month_ECfilt$ch_FCH4_nmolCH4m2s1_median

all_sites_month_ECfilt <- subset(all_sites_month_ECfilt, select = -c(ECCH_diff_mean, ECCH_diff_median))

###### UPDATED DATA: ADD US-STJ
# This section was run with updated US-STJ data already included in the published datasets.
# To rebuild from source files, add the updated US-STJ data analogously to the site-specific sections above.

USSTJ_month <- read.csv("path/US_STJ_chec_monthaggr.csv")

all_month_aggr <- all_sites_month_ECfilt

# change some column names

USSTJ_month <- USSTJ_month %>%
  rename(
    EC_NEE_F_median   = EC_NEE_F_MDS_median,
    EC_NEE_F_mean     = EC_NEE_F_MDS_mean,
    EC_FCH4_median = EC_FCH4_nmolCH4m2s1_median,
    EC_FCH4_mean = EC_FCH4_nmolCH4m2s1_mean
  )

# remove TA and EC_NEE (gap-filled NEE better)

USSTJ_month <- subset(USSTJ_month, select = -c(EC_TA_mean, EC_TA_median, EC_NEE_mean, EC_NEE_median))

USSTJ_month <- USSTJ_month %>%
  rename(
    EC_WS_F_mean   = EC_WS_mean, # not really gap-filled but needs to be in the same column as in rest of the sites
    EC_WS_F_median     = EC_WS_median, # --"--
    EC_VPD_F_mean = EC_VPD_mean, #--"--
    EC_VPD_F_median = EC_VPD_median # --"--
  )

# move values from EC_WD_u.wind etc columns to the right columns in all sites preds
# Pairs of (source column, target column)
col_pairs <- list(
  c("EC_WD_u.wind_mean",      "u.wind_mean"),
  c("EC_WD_u.wind_median",    "u.wind_median"),
  c("EC_WD_u.wind_mean_F",    "u.wind.F_mean"),
  c("EC_WD_u.wind_median_F",  "u.wind.F_median"),
  c("EC_WD_v.wind_mean",      "v.wind_mean"),
  c("EC_WD_v.wind_median",    "v.wind_median"),
  c("EC_WD_v.wind_mean_F",    "v.wind.F_mean"),
  c("EC_WD_v.wind_median_F",  "v.wind.F_median")
)

# Fill missing target values from source values when available
for (pair in col_pairs) {
  src <- pair[1]
  tgt <- pair[2]
  
  idx <- is.na(all_month_aggr[[tgt]]) & !is.na(all_month_aggr[[src]])
  all_month_aggr[[tgt]][idx] <- all_month_aggr[[src]][idx]
}

# Remove EC_WD_... columns except EC_WD_AVG and EC_WD_AVG_F
cols_to_remove <- grep("^EC_WD_", names(all_month_aggr), value = TRUE)
cols_to_keep   <- c("EC_WD_AVG", "EC_WD_AVG_F")
cols_to_remove <- setdiff(cols_to_remove, cols_to_keep)

all_month_aggr <- all_month_aggr[ , !(names(all_month_aggr) %in% cols_to_remove)]

USSTJ_month <- USSTJ_month %>%
  rename(
    v.wind.F_mean   = v.wind_mean, # not really gap-filled but needs to be in the same column as in rest of the sites
    u.wind.F_mean     = u.wind_mean
  )

USSTJ_month <- USSTJ_month %>%
  rename(
    v.wind.F_median   = v.wind_median, # not really gap-filled but needs to be in the same column as in rest of the sites
    u.wind.F_median     = u.wind_median
  )

### need to combine FCH4 values from USSTJ to ANNOPTLM values from other sites
# --> create new column for this

USSTJ_month$EC_FCH4_comb_median <- USSTJ_month$EC_FCH4_median
all_month_aggr$EC_FCH4_comb_median <- all_month_aggr$EC_FCH4_F_ANNOPTLM_median

USSTJ_month$EC_FCH4_comb_mean <- USSTJ_month$EC_FCH4_mean
all_month_aggr$EC_FCH4_comb_mean <- all_month_aggr$EC_FCH4_F_ANNOPTLM_mean

# EC_TS_TOP = EC_TS

USSTJ_month$EC_TS_TOP <- USSTJ_month$EC_TS_mean

USSTJ_month <- USSTJ_month %>%
  rename(
    EC_PA_F_mean   = EC_PA_mean, # not really gap-filled but needs to be in the same column as in rest of the sites
    EC_PA_F_median     = EC_PA_median
  )

USSTJ_month <- subset(USSTJ_month, select = -c(ECCH_diff_mean, ECCH_diff_median, EC_TS_median))

USSTJ_month$ECCH_diff <- USSTJ_month$EC_FCH4_comb_median - USSTJ_month$ch_FCH4_nmolCH4m2s1_median

USSTJ_month$SITE_TS_TOP <- USSTJ_month$EC_TS_TOP

USSTJ_month <- USSTJ_month %>%
  rename(
    EC_WTD_F   = EC_WTD_F_mean
  )

USSTJ_month$SITE_WTD <- USSTJ_month$EC_WTD_F

setdiff(names(all_month_aggr), names(USSTJ_month))
setdiff(names(USSTJ_month), names(all_month_aggr))

cols_2 <- setdiff(names(all_month_aggr), names(USSTJ_month))

USSTJ_month[cols_2] <- NA

# combine
all_sites_month_2 <- bind_rows(all_month_aggr, USSTJ_month)

#### update 11.08.2025: fixing missing US-Owc data (already fixed in the published dataset)

USOWC_month_aggr_ECfilt$EC_TS_TOP <- USOWC_month_aggr_ECfilt$EC_TS_1_mean

USOWC_month_aggr_ECfilt$EC_FCH4_comb_median <- USOWC_month_aggr_ECfilt$EC_FCH4_F_ANNOPTLM_median
USOWC_month_aggr_ECfilt$EC_FCH4_comb_mean <- USOWC_month_aggr_ECfilt$EC_FCH4_F_ANNOPTLM_mean

USOWC_month_aggr_ECfilt$EC_WTD_F <- USOWC_month_aggr_ECfilt$EC_WTD_F_mean * 100

USOWC_month_aggr_ECfilt$SITE_WTD <- USOWC_month_aggr_ECfilt$EC_WTD_F

USOWC_month_aggr_ECfilt$SITE_TS_TOP <- USOWC_month_aggr_ECfilt$EC_TS_TOP

USOWC_month_aggr_ECfilt$ECCH_diff <- USOWC_month_aggr_ECfilt$EC_FCH4_comb_median - USOWC_month_aggr_ECfilt$ch_FCH4_nmolCH4m2s1_median

setdiff(names(all_sites_month_2), names(USOWC_month_aggr_ECfilt))

cols_9 <- setdiff(names(all_sites_month_2), names(USOWC_month_aggr_ECfilt))

USOWC_month_aggr_ECfilt[cols_9] <- NA

# 1. Remove US-OWC from the main df
all_sites_month_3 <- subset(all_sites_month_2, SITE != "US-OWC")

# 2. Ensure column names match between dataframes
USOWC_month_aggr_ECfilt <- USOWC_month_aggr_ECfilt[ , names(all_sites_month_3)]

# 3. Add the rows from USOWC_d_aggr_ECfilt
all_sites_month_3 <- rbind(all_sites_month_3, USOWC_month_aggr_ECfilt)

##### UPDATED DATA: add SE-DEG
# This section was run with updated SE-DEG data already included in the published datasets.
# To rebuild from source files, add the updated SE-DEG data analogously to the site-specific sections above.

setdiff(names(all_sites_month_3), names(SEDEG_month_aggr_ECfilt))

cols_2 <- setdiff(names(all_sites_month_3), names(SEDEG_month_aggr_ECfilt))

SEDEG_month_aggr_ECfilt[cols_2] <- NA

SEDEG_month_aggr_ECfilt$EC_FCH4_comb_mean <- SEDEG_month_aggr_ECfilt$EC_FCH4_F_ANNOPTLM_mean
SEDEG_month_aggr_ECfilt$EC_FCH4_comb_median <- SEDEG_month_aggr_ECfilt$EC_FCH4_F_ANNOPTLM_median

SEDEG_month_aggr_ECfilt$ch_TS_TOP <- SEDEG_month_aggr_ECfilt$ch_TS_2_mean
SEDEG_month_aggr_ECfilt$EC_TS_TOP <- SEDEG_month_aggr_ECfilt$EC_TS_1_mean
SEDEG_month_aggr_ECfilt$SITE_TS_TOP <- mean(c(SEDEG_month_aggr_ECfilt$EC_TS_TOP, SEDEG_month_aggr_ECfilt$ch_TS_TOP), na.rm=T)

SEDEG_month_aggr_ECfilt$EC_WTD_F <- SEDEG_month_aggr_ECfilt$EC_WTD_F_mean * 100
SEDEG_month_aggr_ECfilt$SITE_WTD <- mean(c(SEDEG_month_aggr_ECfilt$EC_WTD_F, SEDEG_month_aggr_ECfilt$ch_WTL_mean))

# 1. Remove US-OWC from the main df
all_sites_month_4 <- subset(all_sites_month_3, SITE != "SE-DEG")

# 2. Ensure column names match between dataframes
SEDEG_month_aggr_ECfilt <- SEDEG_month_aggr_ECfilt[ , names(all_sites_month_4)]

# 3. Add the rows from SEDEG_d_aggr_ECfilt
all_sites_month_4 <- rbind(all_sites_month_4, SEDEG_month_aggr_ECfilt)

# create a combination NEE F column 

all_sites_month_4 <- all_sites_month_4 %>%
  mutate(
    EC_NEE_F_comb_mean   = if_else(SITE == "US-STJ", EC_NEE_F_mean,   EC_NEE_F_ANNOPTLM_mean),
    EC_NEE_F_comb_median = if_else(SITE == "US-STJ", EC_NEE_F_median, EC_NEE_F_ANNOPTLM_median)
  )

write.csv(all_sites_month_4, "path/allsites_ECfilt_monthaggr_12082025.csv", row.names = FALSE)

##### UPDATED DATA: add fixed CN-HGU
# This section was run with updated CN-HGU data already included in the published datasets.
# To rebuild from source files, add the updated CN-HGU data analogously to the site-specific sections above.

all_sites_month_4 <- read.csv("path/allsites_ECfilt_monthaggr_12082025.csv")

setdiff(names(all_sites_month_4), names(CNHGU_month_aggr_ECfilt))

colnames(CNHGU_month_aggr_ECfilt)[2] <- "Year"

cols_2 <- setdiff(names(all_sites_month_4), names(CNHGU_month_aggr_ECfilt))

CNHGU_month_aggr_ECfilt[cols_2] <- NA

CNHGU_month_aggr_ECfilt$EC_FCH4_comb_mean <- CNHGU_month_aggr_ECfilt$EC_FCH4_F_ANNOPTLM_mean
CNHGU_month_aggr_ECfilt$EC_FCH4_comb_median <- CNHGU_month_aggr_ECfilt$EC_FCH4_F_ANNOPTLM_median

CNHGU_month_aggr_ECfilt$EC_TS_TOP <- CNHGU_month_aggr_ECfilt$EC_TS_F_mean
CNHGU_month_aggr_ECfilt$SITE_TS_TOP <- CNHGU_month_aggr_ECfilt$EC_TS_TOP

CNHGU_month_aggr_ECfilt$EC_NEE_F_comb_mean <- CNHGU_month_aggr_ECfilt$EC_NEE_F_ANNOPTLM_mean
CNHGU_month_aggr_ECfilt$EC_NEE_F_comb_median <- CNHGU_month_aggr_ECfilt$EC_NEE_F_ANNOPTLM_median

# 1. Remove CN-HGU from the main df
all_sites_month_5 <- subset(all_sites_month_4, SITE != "CN-HGU")

# 2. Ensure column names match between dataframes
CNHGU_month_aggr_ECfilt <- CNHGU_month_aggr_ECfilt[ , names(all_sites_month_5)]

# 3. Add the rows from CNHGU_d_aggr_ECfilt
all_sites_month_5 <- rbind(all_sites_month_5, CNHGU_month_aggr_ECfilt)

write.csv(all_sites_month_5, "path/allsites_ECfilt_monthaggr_28082025.csv", row.names = FALSE)

# save as .csv

write.csv(all_sites_month_ECfilt, "path/allsites_ECfilt_monthaggr.csv", row.names = FALSE)

##### ANNUAL DATA SETS #####

# read in datasets

CNHGU_yr_aggr_ECfilt <- read.csv("path/CN_HGU_chec_yraggr_ECfilt.csv")

SEDEG_yr_aggr_ECfilt <- read.csv("path/SE_DEG_chec_yraggr_ECfilt.csv")

USHO1_yr_aggr_ECfilt <- read.csv("path/US_HO1_chec_yraggr_ECfilt.csv")

USUAF_yr_aggr_ECfilt <- read.csv("path/US_UAF_chec_yraggr_ECfilt.csv")

FISI2_yr_aggr_ECfilt <- read.csv("path/FI_SI2_chec_yraggr_ECfilt.csv")

USLA1_yr_aggr_ECfilt <- read.csv("path/US_LA1_chec_yraggr_ECfilt.csv")

USLA2_yr_aggr_ECfilt <- read.csv("path/US_LA2_chec_yraggr_ECfilt.csv")

USLOS_yr_aggr_ECfilt <- read.csv("path/US_LOS_chec_yraggr_ECfilt.csv")

USOWC_yr_aggr_ECfilt <- read.csv("path/US_OWC_chec_yraggr_ECfilt.csv")

# harmonize datasets

USHO1_yr_aggr_ECfilt <- USHO1_yr_aggr_ECfilt %>% rename("ch_FCH4_nmolCH4m2s1_mean" = "ch_FCH4_umolCH4m2s1_mean",
                                                        "ch_FCH4_nmolCH4m2s1_median" = "ch_FCH4_umolCH4m2s1_median"
)

# remove the ch measurements of the same EC neasurennents

SEDEG_yr_aggr_ECfilt <- subset(SEDEG_yr_aggr_ECfilt, select=-c(ch_PPT_mean, ch_WS_mean, ch_P_mean, ch_VPD_mean, ch_RH_mean,
                                                               ch_PPT_median, ch_WS_median, ch_P_median, ch_VPD_median, ch_RH_median))

# add missing columns to dataframes

# rename ch_WT to ch_WTL

FISI2_yr_aggr_ECfilt <- FISI2_yr_aggr_ECfilt %>% rename("ch_WTL_mean" = "ch_WT_mean",
                                                        "ch_WTL_median" = "ch_WT_median"
)

# remove ch_TS_5 and 15

FISI2_yr_aggr_ECfilt <- subset(FISI2_yr_aggr_ECfilt, select=-c(ch_TS_15_mean,  ch_TS_15_median))
USLA1_yr_aggr_ECfilt <- subset(USLA1_yr_aggr_ECfilt, select=-c(ch_TS_5_mean, ch_TS_5_median,  ch_TS_15_mean,  ch_TS_15_median))
USLA2_yr_aggr_ECfilt <- subset(USLA2_yr_aggr_ECfilt, select=-c(ch_TS_5_mean, ch_TS_5_median, ch_TS_15_mean, ch_TS_15_median))

# define missing columns for CN-HGU
cols_1 <- c("ch_TS_2_mean", "ch_TS_5_mean", "ch_TS_5_median", "ch_TS_10_mean", "ch_TS_20_mean", "ch_TS_30_mean", "ch_TS_40_mean", 
            "ch_CH4_Ta_amb_mean", 
            "ch_TS_2_median",
            "EC_SWC_1_mean", "EC_PA_mean", "EC_TS_1_mean", "EC_TS_2_mean", 
            "EC_TS_3_mean", "EC_TS_4_mean", "EC_TS_5_mean", "EC_TS_6_mean", "EC_TS_7_mean", "EC_TS_8_mean", "EC_TS_9_mean", 
            "EC_WTD_F_mean", "EC_PPFD_OUT_mean",
            "EC_SWC_2_mean", "EC_SWC_3_mean", "EC_PPFD_OUT_F_mean", 
            "EC_SWC_1_median", "EC_PA_median", "EC_TS_1_median", "EC_TS_2_median", 
            "EC_TS_3_median", "EC_TS_4_median", "EC_TS_5_median", "EC_TS_6_median", "EC_TS_7_median", "EC_TS_8_median", "EC_TS_9_median", 
            "EC_WTD_F_median", "EC_PPFD_OUT_median",
            "EC_SWC_2_median", "EC_SWC_3_median", "EC_PPFD_OUT_F_median", 
            "ch_SM_X_mean", 
            "ch_TS_X_mean", "ch_SM_X_median", "ch_TS_X_median", 
            "ch_VWC_mean",
            "ch_TS_20_median", "ch_TS_30_median", "ch_TS_40_median", "ch_VWC_median", 
            "ch_AerLAI_mean", "ch_AerLAI_median", 
            "salinity_mean",  "salinity_median")

CNHGU_yr_aggr_ECfilt[cols_1] <- NA

# now CNHGU is a df that includes all columns --> compare the rest to this

setdiff(names(CNHGU_yr_aggr_ECfilt), names(SEDEG_yr_aggr_ECfilt))

cols_2 <- setdiff(names(CNHGU_yr_aggr_ECfilt), names(SEDEG_yr_aggr_ECfilt))

SEDEG_yr_aggr_ECfilt[cols_2] <- NA

#---

setdiff(names(CNHGU_yr_aggr_ECfilt), names(USHO1_yr_aggr_ECfilt))

cols_3 <- setdiff(names(CNHGU_yr_aggr_ECfilt), names(USHO1_yr_aggr_ECfilt))

USHO1_yr_aggr_ECfilt[cols_3] <- NA

#---

setdiff(names(CNHGU_yr_aggr_ECfilt), names(USUAF_yr_aggr_ECfilt))

cols_4 <- setdiff(names(CNHGU_yr_aggr_ECfilt), names(USUAF_yr_aggr_ECfilt))

USUAF_yr_aggr_ECfilt[cols_4] <- NA

#---

setdiff(names(CNHGU_yr_aggr_ECfilt), names(FISI2_yr_aggr_ECfilt))

cols_5 <- setdiff(names(CNHGU_yr_aggr_ECfilt), names(FISI2_yr_aggr_ECfilt))

FISI2_yr_aggr_ECfilt[cols_5] <- NA

#---

setdiff(names(CNHGU_yr_aggr_ECfilt), names(USLA1_yr_aggr_ECfilt))

cols_6 <- setdiff(names(CNHGU_yr_aggr_ECfilt), names(USLA1_yr_aggr_ECfilt))

USLA1_yr_aggr_ECfilt[cols_6] <- NA

#---

setdiff(names(CNHGU_yr_aggr_ECfilt), names(USLA2_yr_aggr_ECfilt))

cols_7 <- setdiff(names(CNHGU_yr_aggr_ECfilt), names(USLA2_yr_aggr_ECfilt))

USLA2_yr_aggr_ECfilt[cols_7] <- NA

#---

setdiff(names(CNHGU_yr_aggr_ECfilt), names(USLOS_yr_aggr_ECfilt))

cols_8 <- setdiff(names(CNHGU_yr_aggr_ECfilt), names(USLOS_yr_aggr_ECfilt))

USLOS_yr_aggr_ECfilt[cols_8] <- NA

#---

setdiff(names(CNHGU_yr_aggr_ECfilt), names(USOWC_yr_aggr_ECfilt))

cols_9 <- setdiff(names(CNHGU_yr_aggr_ECfilt), names(USOWC_yr_aggr_ECfilt))

USOWC_yr_aggr_ECfilt[cols_9] <- NA

# combine

all_sites_annual_ECfilt <- bind_rows(CNHGU_yr_aggr_ECfilt, FISI2_yr_aggr_ECfilt, SEDEG_yr_aggr_ECfilt, USHO1_yr_aggr_ECfilt,
                                     USLA1_yr_aggr_ECfilt, USLA2_yr_aggr_ECfilt, USLOS_yr_aggr_ECfilt, USOWC_yr_aggr_ECfilt, USUAF_yr_aggr_ECfilt)

# CN-HGU accidentally has timestamp instead of year

all_sites_annual_ECfilt$Year[1] <- 2015
all_sites_annual_ECfilt$Year[2] <- 2016

all_sites_annual_ECfilt <- subset(all_sites_annual_ECfilt, select = -TIMESTAMP)

all_sites_annual_ECfilt <- all_sites_annual_ECfilt %>% relocate(Year, .after = SITE)

# combine TS info into one column

# Create the new columns based on the "SITE" group
all_sites_annual_ECfilt <- all_sites_annual_ECfilt %>%
  mutate(
    ch_TS_TOP = case_when(
      SITE == 'FI-SI2' ~ ch_TS_5_mean,
      SITE == 'CN-HGU' ~ NA,
      SITE == 'SE-DEG' ~ ch_TS_2_mean,
      SITE == 'US-HO1' ~ ch_TS_X_mean,
      SITE == 'US-LA1' ~ ch_TS_10_mean,
      SITE == 'US-LA2' ~ ch_TS_10_mean,
      SITE == 'US-LOS' ~ ch_TS_X_mean,
      SITE == 'US-OWC' ~ NA,
      SITE == 'US-UAF' ~ ch_TS_2_mean,
      TRUE ~ NA_real_
    ),
    EC_TS_TOP = case_when(
      SITE == 'FI-SI2' ~ EC_TS_1_mean,
      SITE == 'CN-HGU' ~ EC_TS_F_mean,
      SITE == 'SE-DEG' ~ EC_TS_1_mean,
      SITE == 'US-HO1' ~ EC_TS_1_mean,
      SITE == 'US-LA1' ~ EC_TS_F_mean,
      SITE == 'US-LA2' ~ EC_TS_F_mean,
      SITE == 'US-LOS' ~ EC_TS_1_mean,
      SITE == 'US-OWC' ~ EC_TS_1_mean,
      SITE == 'US-UAF' ~ EC_TS_1_mean,
      TRUE ~ NA_real_
    )
  )

# create a SITE_TS_TOP

all_sites_annual_ECfilt <- all_sites_annual_ECfilt %>%
  mutate(SITE_TS_TOP = case_when(
    SITE %in% c("FI-SI2", "US-HO1") ~ ch_TS_TOP,
    SITE %in% c("US-OWC", "CN-HGU") ~ EC_TS_TOP,
    TRUE ~ rowMeans(cbind(ch_TS_TOP, EC_TS_TOP), na.rm = TRUE)
  ))

# do the same for EC_SWC

# Create the new columns based on the "SITE" group
all_sites_annual_ECfilt <- all_sites_annual_ECfilt %>%
  mutate(
    EC_SWC_TOP = case_when(
      SITE == 'FI-SI2' ~ NA,
      SITE == 'CN-HGU' ~ EC_SWC_F_mean,
      SITE == 'SE-DEG' ~ EC_SWC_1_mean,
      SITE == 'US-HO1' ~ NA,
      SITE == 'US-LA1' ~ NA,
      SITE == 'US-LA2' ~ NA,
      SITE == 'US-LOS' ~ NA,
      SITE == 'US-OWC' ~ NA,
      SITE == 'US-UAF' ~ EC_SWC_1_mean,
      TRUE ~ NA_real_
    ),
    ch_SWC_TOP = case_when(
      SITE == 'FI-SI2' ~ NA,
      SITE == 'CN-HGU' ~ NA,
      SITE == 'SE-DEG' ~ NA,
      SITE == 'US-HO1' ~ ch_SM_X_mean,
      SITE == 'US-LA1' ~ NA,
      SITE == 'US-LA2' ~ NA,
      SITE == 'US-LOS' ~ ch_SM_X_mean,
      SITE == 'US-OWC' ~ NA,
      SITE == 'US-UAF' ~ ch_VWC_mean,
      TRUE ~ NA_real_
    )
  )

# create SITE_SWC

all_sites_annual_ECfilt <- all_sites_annual_ECfilt %>%
  mutate(SITE_SWC = case_when(
    SITE %in% c("US-HO1") ~ ch_SWC_TOP,
    TRUE ~ rowMeans(cbind(ch_SWC_TOP, EC_SWC_TOP), na.rm = TRUE)
  ))

# create SITE_WTD

# EC_WTD is in m so convert EC_WTD to cm

all_sites_annual_ECfilt$EC_WTD_F <- all_sites_annual_ECfilt$EC_WTD_F_mean * 100

# Create the new column SITE_WTD
all_sites_annual_ECfilt <- all_sites_annual_ECfilt %>%
  mutate(SITE_WTD = case_when(
    SITE %in% c("FI-SI2", "US-LA1", "US-LA2") ~ ch_WTL_mean,
    SITE %in% c("US-HO1", "US-OWC", "US-UAF") ~ EC_WTD_F,
    SITE == "CN-HGU" ~ NA,
    TRUE ~ rowMeans(cbind(ch_WTL_mean, EC_WTD_F), na.rm = TRUE)
  ))

# remove the extra columns

all_sites_annual_ECfilt <- subset(all_sites_annual_ECfilt, select = -c(EC_TS_1_mean, EC_TS_2_mean, EC_TS_3_mean, EC_TS_4_mean, EC_TS_5_mean,
                                                                       EC_TS_6_mean, EC_TS_7_mean, EC_TS_8_mean, EC_TS_9_mean,
                                                                       ch_SM_X_mean, ch_TS_X_mean, ch_TS_2_mean, ch_TS_10_mean, 
                                                                       ch_TS_20_mean, ch_TS_30_mean, ch_TS_40_mean,
                                                                       ch_VWC_mean, EC_SWC_1_mean, EC_SWC_2_mean, EC_SWC_3_mean))

all_sites_annual_ECfilt <- subset(all_sites_annual_ECfilt, select = -c(EC_TS_1_median, EC_TS_2_median, EC_TS_3_median, EC_TS_4_median, EC_TS_5_median,
                                                                       EC_TS_6_median, EC_TS_7_median, EC_TS_8_median, EC_TS_9_median,
                                                                       ch_SM_X_median, ch_TS_X_median, ch_TS_2_median, ch_TS_10_median, 
                                                                       ch_TS_20_median, ch_TS_30_median, ch_TS_40_median,
                                                                       ch_VWC_median, EC_SWC_1_median, EC_SWC_2_median, EC_SWC_3_median))

# save as .csv

write.csv(all_sites_annual_ECfilt, "path/allsites_ECfilt_yraggr.csv", row.names = FALSE)

###### UPDATED DATA: ADD US-STJ
# This section was run with updated US-STJ data already included in the published datasets.
# To rebuild from source files, add the updated US-STJ data analogously to the site-specific sections above.

USSTJ_yr <- read.csv("path/US_STJ_chec_yraggr.csv")

all_yr_aggr <- all_sites_annual_ECfilt

# change some column names

USSTJ_yr <- USSTJ_yr %>%
  rename(
    EC_NEE_F_median   = EC_NEE_F_MDS_median,
    EC_NEE_F_mean     = EC_NEE_F_MDS_mean,
    EC_FCH4_median = EC_FCH4_nmolCH4m2s1_median,
    EC_FCH4_mean = EC_FCH4_nmolCH4m2s1_mean
  )

# remove TA and EC_NEE (gap-filled NEE better)

USSTJ_yr <- subset(USSTJ_yr, select = -c(EC_TA_mean, EC_TA_median, EC_NEE_mean, EC_NEE_median))

USSTJ_yr <- USSTJ_yr %>%
  rename(
    EC_WS_F_mean   = EC_WS_mean, # not really gap-filled but needs to be in the same column as in rest of the sites
    EC_WS_F_median     = EC_WS_median, # --"--
    EC_VPD_F_mean = EC_VPD_mean, #--"--
    EC_VPD_F_median = EC_VPD_median # --"--
  )

# need to move values from EC_WD_u.wind etc columns to the right columns in all sites preds
# Pairs of (source column, target column)
col_pairs <- list(
  c("EC_WD_u.wind_mean",      "u.wind_mean"),
  c("EC_WD_u.wind_median",    "u.wind_median"),
  c("EC_WD_u.wind_mean_F",    "u.wind.F_mean"),
  c("EC_WD_u.wind_median_F",  "u.wind.F_median"),
  c("EC_WD_v.wind_mean",      "v.wind_mean"),
  c("EC_WD_v.wind_median",    "v.wind_median"),
  c("EC_WD_v.wind_mean_F",    "v.wind.F_mean"),
  c("EC_WD_v.wind_median_F",  "v.wind.F_median")
)

# Fill missing target values from source values when available
for (pair in col_pairs) {
  src <- pair[1]
  tgt <- pair[2]
  
  idx <- is.na(all_yr_aggr[[tgt]]) & !is.na(all_yr_aggr[[src]])
  all_yr_aggr[[tgt]][idx] <- all_yr_aggr[[src]][idx]
}

# Remove EC_WD_... columns except EC_WD_AVG and EC_WD_AVG_F
cols_to_remove <- grep("^EC_WD_", names(all_yr_aggr), value = TRUE)
cols_to_keep   <- c("EC_WD_AVG", "EC_WD_AVG_F")
cols_to_remove <- setdiff(cols_to_remove, cols_to_keep)

all_yr_aggr <- all_yr_aggr[ , !(names(all_yr_aggr) %in% cols_to_remove)]

USSTJ_yr <- USSTJ_yr %>%
  rename(
    v.wind.F_mean   = v.wind_mean, # not really gap-filled but needs to be in the same column as in rest of the sites
    u.wind.F_mean     = u.wind_mean
  )

USSTJ_yr <- USSTJ_yr %>%
  rename(
    v.wind.F_median   = v.wind_median, # not really gap-filled but needs to be in the same column as in rest of the sites
    u.wind.F_median     = u.wind_median
  )

USSTJ_yr$EC_height_m <- 3.5

### need to combine FCH4 values from USSTJ to ANNOPTLM values from other sites
# --> create new column for this

USSTJ_yr$EC_FCH4_comb_median <- USSTJ_yr$EC_FCH4_median
all_yr_aggr$EC_FCH4_comb_median <- all_yr_aggr$EC_FCH4_F_ANNOPTLM_median

USSTJ_yr$EC_FCH4_comb_mean <- USSTJ_yr$EC_FCH4_mean
all_yr_aggr$EC_FCH4_comb_mean <- all_yr_aggr$EC_FCH4_F_ANNOPTLM_mean

# EC_TS_TOP = EC_TS

USSTJ_yr$EC_TS_TOP <- USSTJ_yr$EC_TS_mean

USSTJ_yr <- USSTJ_yr %>%
  rename(
    EC_PA_F_mean   = EC_PA_mean, # not really gap-filled but needs to be in the same column as in rest of the sites
    EC_PA_F_median     = EC_PA_median
  )

USSTJ_yr <- subset(USSTJ_yr, select = -c(ECCH_diff_mean, ECCH_diff_median, EC_height_m, EC_TS_median))

USSTJ_yr$ECCH_diff <- USSTJ_yr$EC_FCH4_comb_median - USSTJ_yr$ch_FCH4_nmolCH4m2s1_median

USSTJ_yr$SITE_TS_TOP <- USSTJ_yr$EC_TS_TOP

USSTJ_yr <- USSTJ_yr %>%
  rename(
    EC_WTD_F   = EC_WTD_F_mean
  )

USSTJ_yr$SITE_WTD <- USSTJ_yr$EC_WTD_F

setdiff(names(all_yr_aggr), names(USSTJ_yr))
setdiff(names(USSTJ_yr), names(all_yr_aggr))

cols_2 <- setdiff(names(all_yr_aggr), names(USSTJ_yr))

USSTJ_yr[cols_2] <- NA

# combine
all_sites_yr_2 <- bind_rows(all_yr_aggr, USSTJ_yr)

#### update 11.08.2025: fixing missing US-Owc data (already fixed in the published dataset)

USOWC_yr_aggr_ECfilt$EC_TS_TOP <- USOWC_yr_aggr_ECfilt$EC_TS_1_mean

USOWC_yr_aggr_ECfilt$EC_FCH4_comb_median <- USOWC_yr_aggr_ECfilt$EC_FCH4_F_ANNOPTLM_median
USOWC_yr_aggr_ECfilt$EC_FCH4_comb_mean <- USOWC_yr_aggr_ECfilt$EC_FCH4_F_ANNOPTLM_mean

USOWC_yr_aggr_ECfilt$EC_WTD_F <- USOWC_yr_aggr_ECfilt$EC_WTD_F_mean * 100

USOWC_yr_aggr_ECfilt$SITE_WTD <- USOWC_yr_aggr_ECfilt$EC_WTD_F

USOWC_yr_aggr_ECfilt$SITE_TS_TOP <- USOWC_yr_aggr_ECfilt$EC_TS_TOP

USOWC_yr_aggr_ECfilt$ECCH_diff <- USOWC_yr_aggr_ECfilt$EC_FCH4_comb_median - USOWC_yr_aggr_ECfilt$ch_FCH4_nmolCH4m2s1_median

setdiff(names(all_sites_yr_2), names(USOWC_yr_aggr_ECfilt))

cols_9 <- setdiff(names(all_sites_yr_2), names(USOWC_yr_aggr_ECfilt))

USOWC_yr_aggr_ECfilt[cols_9] <- NA

# 1. Remove US-OWC from the main df
all_sites_yr_3 <- subset(all_sites_yr_2, SITE != "US-OWC")

# 2. Ensure column names match between dataframes
USOWC_yr_aggr_ECfilt <- USOWC_yr_aggr_ECfilt[ , names(all_sites_yr_3)]

# 3. Add the rows from USOWC_d_aggr_ECfilt
all_sites_yr_3 <- rbind(all_sites_yr_3, USOWC_yr_aggr_ECfilt)

# create a combination NEE F column 

all_sites_yr_4 <- all_sites_yr_4 %>%
  mutate(
    EC_NEE_F_comb_mean   = if_else(SITE == "US-STJ", EC_NEE_F_mean,   EC_NEE_F_ANNOPTLM_mean),
    EC_NEE_F_comb_median = if_else(SITE == "US-STJ", EC_NEE_F_median, EC_NEE_F_ANNOPTLM_median)
  )

write.csv(all_sites_yr_4, "path/allsites_ECfilt_yraggr_11082025.csv", row.names = FALSE)

##### UPDATED DATA: add SE-DEG
# This section was run with updated SE-DEG data already included in the published datasets.
# To rebuild from source files, add the updated SE-DEG data analogously to the site-specific sections above.

setdiff(names(all_sites_yr_3), names(SEDEG_yr_aggr_ECfilt))

cols_2 <- setdiff(names(all_sites_yr_3), names(SEDEG_yr_aggr_ECfilt))

SEDEG_yr_aggr_ECfilt[cols_2] <- NA

SEDEG_yr_aggr_ECfilt$EC_FCH4_comb_mean <- SEDEG_yr_aggr_ECfilt$EC_FCH4_F_ANNOPTLM_mean
SEDEG_yr_aggr_ECfilt$EC_FCH4_comb_median <- SEDEG_yr_aggr_ECfilt$EC_FCH4_F_ANNOPTLM_median

SEDEG_yr_aggr_ECfilt$ch_TS_TOP <- SEDEG_yr_aggr_ECfilt$ch_TS_2_mean
SEDEG_yr_aggr_ECfilt$EC_TS_TOP <- SEDEG_yr_aggr_ECfilt$EC_TS_1_mean
SEDEG_yr_aggr_ECfilt$SITE_TS_TOP <- mean(c(SEDEG_yr_aggr_ECfilt$EC_TS_TOP, SEDEG_yr_aggr_ECfilt$ch_TS_TOP), na.rm=T)

SEDEG_yr_aggr_ECfilt$EC_WTD_F <- SEDEG_yr_aggr_ECfilt$EC_WTD_F_mean * 100
SEDEG_yr_aggr_ECfilt$SITE_WTD <- mean(c(SEDEG_yr_aggr_ECfilt$EC_WTD_F, SEDEG_yr_aggr_ECfilt$ch_WTL_mean))

# 1. Remove SE-DEG from the main df
all_sites_yr_4 <- subset(all_sites_yr_3, SITE != "SE-DEG")

# 2. Ensure column names match between dataframes
SEDEG_yr_aggr_ECfilt <- SEDEG_yr_aggr_ECfilt[ , names(all_sites_yr_4)]

# 3. Add the rows from SEDEG_d_aggr_ECfilt
all_sites_yr_4 <- rbind(all_sites_yr_4, SEDEG_yr_aggr_ECfilt)

##### UPDATED DATA: add fixed CN-HGU
# This section was run later with updated CN-HGU data already included in the published datasets.
# To rebuild from source files, add the updated CN-HGU data analogously to the site-specific sections above.

all_sites_yr_4 <- read.csv("path/allsites_ECfilt_yraggr_11082025.csv")

setdiff(names(all_sites_yr_4), names(CNHGU_yr_aggr_ECfilt))

colnames(CNHGU_yr_aggr_ECfilt)[2] <- "Year"

cols_2 <- setdiff(names(all_sites_yr_4), names(CNHGU_yr_aggr_ECfilt))

CNHGU_yr_aggr_ECfilt[cols_2] <- NA

CNHGU_yr_aggr_ECfilt$EC_FCH4_comb_mean <- CNHGU_yr_aggr_ECfilt$EC_FCH4_F_ANNOPTLM_mean
CNHGU_yr_aggr_ECfilt$EC_FCH4_comb_median <- CNHGU_yr_aggr_ECfilt$EC_FCH4_F_ANNOPTLM_median

CNHGU_yr_aggr_ECfilt$EC_TS_TOP <- CNHGU_yr_aggr_ECfilt$EC_TS_F_mean
CNHGU_yr_aggr_ECfilt$SITE_TS_TOP <- CNHGU_yr_aggr_ECfilt$EC_TS_TOP

# 1. Remove CN-HGU from the main df
all_sites_yr_5 <- subset(all_sites_yr_4, SITE != "CN-HGU")

# 2. Ensure column names match between dataframes
CNHGU_yr_aggr_ECfilt <- CNHGU_yr_aggr_ECfilt[ , names(all_sites_yr_5)]

# 3. Add the rows from CNHGU_d_aggr_ECfilt
all_sites_yr_5 <- rbind(all_sites_yr_5, CNHGU_yr_aggr_ECfilt)

all_sites_yr_5 <- all_sites_yr_5 %>%
  mutate(
    EC_NEE_F_comb_mean   = if_else(SITE == "US-STJ", EC_NEE_F_mean,   EC_NEE_F_ANNOPTLM_mean),
    EC_NEE_F_comb_median = if_else(SITE == "US-STJ", EC_NEE_F_median, EC_NEE_F_ANNOPTLM_median)
  )

write.csv(all_sites_yr_5, "path/allsites_ECfilt_yraggr_28082025.csv", row.names = FALSE)

#########################################################################################################################################
### PREPARE DATASETS FOR PUBLICATION ###
### the resulting datasets were used in Maatta_et_al_2026_Statistical_analyses.R ###

## HALF-HOURLY AGGREGATION ##
auto_sites_preds_chaggr <- read.csv("path/autosites_preds_chaggr_NEW_VERSION_28082025.csv")

# remove extra columns not used in the analyses

auto_sites_preds_chaggr <- subset(auto_sites_preds_chaggr, select = -c(DOY, ch_TS_TOP_mean, ch_TS_TOP_median, 
                                                                       ch_SWC_TOP_mean, ch_SWC_TOP_median, 
                                                                       EC_FCH4_F, EC_NEE, EC_H, EC_LE, EC_SW_OUT, 
                                                                       EC_SW_IN, EC_LW_IN, EC_LW_OUT, 
                                                                       EC_NETRAD, EC_PPFD_IN, EC_VPD, EC_RH, 
                                                                       EC_TA, EC_TS, EC_G, EC_SWC, EC_GPP_NT, 
                                                                       EC_RECO_NT, EC_GPP_DT, EC_RECO_DT, 
                                                                       EC_NEE_F, EC_H_F, EC_LE_F, 
                                                                       EC_SW_IN_F, EC_SW_OUT_F, 
                                                                       EC_LW_IN_F, EC_LW_OUT_F, 
                                                                       EC_NETRAD_F, EC_PPFD_IN_F, 
                                                                       EC_RH_F, 
                                                                       EC_TA_F, EC_TS_F, EC_G_F, EC_SWC_F, 
                                                                       EC_LE_F_ANNOPTLM, EC_WTD, EC_WTD_F, EC_PPFD_OUT, EC_PPFD_OUT_F, 
                                                                       EC_height_m, u.wind, v.wind, EC_TS_TOP, 
                                                                       EC_SWC_TOP, ECCH_diff_abs))

auto_sites_preds_chaggr <- subset(auto_sites_preds_chaggr, select = -c(EC_SENSOR, SITE_SWC, EC_FCH4, EC_FCH4_F_ANNOPTLM,
                                                                       ch_WTL_mean, ch_WTL_median,
                                                                       KOPPEN, MOSS_BROWN, MOSS_SPHAGNUM, AERENCHYMATOUS,
                                                                       ERI_SHRUB, TREE))

auto_sites_preds_chaggr <- subset(auto_sites_preds_chaggr, select = -season)
auto_sites_preds_chaggr <- subset(auto_sites_preds_chaggr, select = -EC_FCH4_F_ANNOPTLM_QC)

# reorder columns
auto_sites_preds_chaggr <- auto_sites_preds_chaggr %>%
  dplyr::select(
    SITE,
    starts_with("TIMESTAMP"),
    DATE,
    Year, Month, Day, Hour,
    starts_with("ch_FCH4"),
    EC_FCH4,
    deltaFCH4,
    SITE_TS_TOP, SITE_WTL,
    EC_USTAR, EC_NEE_F_ANNOPTLM,
    EC_PA, EC_PA_F,
    EC_VPD_F,
    EC_WD, EC_WD_F,
    EC_WS, EC_WS_F,
    u_wind_F, v_wind_F,
    DOM_VEG, ch_method,
    SITE_CLASSIFICATION,
    EC_P, EC_P_F
  )

names(auto_sites_preds_chaggr) <- toupper(names(auto_sites_preds_chaggr))

auto_sites_preds_chaggr <- auto_sites_preds_chaggr %>% dplyr::rename("CH_FCH4_MEAN" = "CH_FCH4_NMOLCH4M2S1_MEAN",
                                                                     "CH_FCH4_MEDIAN" = "CH_FCH4_NMOLCH4M2S1_MEDIAN",
                                                                     "DELTA_FCH4" = "DELTAFCH4",
                                                                     "NEE_F" = "EC_NEE_F_ANNOPTLM",
                                                                     "VEG" = "DOM_VEG"
                                                                     
)

auto_sites_preds_chaggr <- auto_sites_preds_chaggr %>% dplyr::rename("DOMINANT_VEGETATION" = "VEG"
                                                                     
)

auto_sites_preds_chaggr <- auto_sites_preds_chaggr %>%
  rename_with(
    ~ gsub("^EC_", "", .x),       # remove "EC_" from the start
    .cols = !contains("EC_FCH4")  # but skip columns that contain "EC_FCH4"
  )

auto_sites_preds_chaggr <- auto_sites_preds_chaggr %>%
  arrange(SITE)

auto_sites_preds_chaggr <- auto_sites_preds_chaggr %>%
  mutate(CH_METHOD = if_else(CH_METHOD == "auto", "automated", CH_METHOD))

auto_sites_preds_chaggr <- auto_sites_preds_chaggr %>%
  relocate(EC_FCH4, .before = CH_FCH4_MEAN)

## save

write.csv(auto_sites_preds_chaggr, "path/halfhourly_aggregation.csv", row.names = FALSE)

## HOURLY AGGREGATION ##
all_hr_aggr <- read.csv("path/autosites_hraggr_ECfilt.csv")

# remove extra columns not used in the analyses

all_hr_aggr <- subset(all_hr_aggr, select = c(SITE, TIMESTAMP, Year, Month, Day, Hour, ch_FCH4_nmolCH4m2s1_mean,
                                              ch_FCH4_nmolCH4m2s1_median, EC_FCH4_F_ANNOPTLM_mean, 
                                              EC_FCH4_F_ANNOPTLM_median, deltaFCH4, 
                                              SITE_TS_TOP, SITE_WTD, EC_NEE_F_ANNOPTLM_mean, EC_USTAR_median, 
                                              EC_VPD_F_median, EC_PA_median, EC_PA_F_median, u.wind.F_mean, v.wind.F_mean, DOM_VEG,
                                              ch_method, SITE_CLASSIFICATION
                                              
))

# rename columns

all_hr_aggr <- all_hr_aggr %>% dplyr::rename("SITE_WTL" = "SITE_WTD",
                                             "EC_FCH4_mean" = "EC_FCH4_F_ANNOPTLM_mean",
                                             "EC_FCH4_median" = "EC_FCH4_F_ANNOPTLM_median",
                                             "u_wind_F_mean" = "u.wind.F_mean",
                                             "v_wind_F_mean" = "v.wind.F_mean",
                                             "DOMINANT_VEGETATION" = "DOM_VEG",
                                             "ch_FCH4_mean" = "ch_FCH4_nmolCH4m2s1_mean",
                                             "ch_FCH4_median" = "ch_FCH4_nmolCH4m2s1_median"
                                             
)

all_hr_aggr <- all_hr_aggr %>%
  rename_with(
    ~ gsub("^EC_", "", .x),       # remove "EC_" from the start
    .cols = !contains("EC_FCH4")  # but skip columns that contain "EC_FCH4"
  )

all_hr_aggr <- all_hr_aggr %>% dplyr::rename("NEE_F_mean" = "NEE_F_ANNOPTLM_mean"
                                             
)

names(all_hr_aggr) <- toupper(names(all_hr_aggr))

all_hr_aggr <- all_hr_aggr %>% dplyr::rename("DELTA_FCH4" = "DELTAFCH4"
                                             
)

# alphabetical order according to site
all_hr_aggr <- all_hr_aggr %>%
  arrange(SITE)

all_hr_aggr <- all_hr_aggr %>%
  mutate(CH_METHOD = if_else(CH_METHOD == "auto", "automated", CH_METHOD))

all_hr_aggr <- all_hr_aggr %>%
  relocate(EC_FCH4_MEAN, EC_FCH4_MEDIAN, .before = CH_FCH4_MEAN)

all_hr_aggr <- all_hr_aggr %>% drop_na(CH_FCH4_MEDIAN)

all_hr_aggr$DATE <- date(all_hr_aggr$TIMESTAMP)

all_hr_aggr <- all_hr_aggr %>%
  relocate(DATE, .before = YEAR)

## save

write.csv(all_hr_aggr, "path/hourly_aggregation.csv", row.names = FALSE)

## DAILY AGGREGATION ##
all_d_aggr <- read.csv("path/allsites_daily_aggr_28082025.csv")

all_d_aggr$DELTA_FCH4 <- all_d_aggr$EC_FCH4_comb_median - all_d_aggr$ch_FCH4_nmolCH4m2s1_median

# remove extra columns not used in the analyses

all_d_aggr <- subset(all_d_aggr, select = c(SITE, TIMESTAMP, Year, Month, Day, ch_FCH4_nmolCH4m2s1_mean,
                                            ch_FCH4_nmolCH4m2s1_median, EC_FCH4_comb_mean, 
                                            EC_FCH4_comb_median, DELTA_FCH4, 
                                            SITE_TS_TOP, SITE_WTD, EC_NEE_F_comb_mean, EC_USTAR_median, 
                                            EC_VPD_F_median, EC_PA_median, EC_PA_F_median, u.wind.F_mean, v.wind.F_mean, DOM_VEG,
                                            ch_method, SITE_CLASSIFICATION
                                            
))

# rename columns

all_d_aggr <- all_d_aggr %>% dplyr::rename("SITE_WTL" = "SITE_WTD",
                                           "EC_FCH4_mean" = "EC_FCH4_comb_mean",
                                           "EC_FCH4_median" = "EC_FCH4_comb_median",
                                           "u_wind_F_mean" = "u.wind.F_mean",
                                           "v_wind_F_mean" = "v.wind.F_mean",
                                           "DOMINANT_VEGETATION" = "DOM_VEG",
                                           "ch_FCH4_mean" = "ch_FCH4_nmolCH4m2s1_mean",
                                           "ch_FCH4_median" = "ch_FCH4_nmolCH4m2s1_median"
                                           
)

all_d_aggr <- all_d_aggr %>%
  rename_with(
    ~ gsub("^EC_", "", .x),       # remove "EC_" from the start
    .cols = !contains("EC_FCH4")  # but skip columns that contain "EC_FCH4"
  )

all_d_aggr <- all_d_aggr %>% dplyr::rename("NEE_F_mean" = "NEE_F_comb_mean"
                                           
)

names(all_d_aggr) <- toupper(names(all_d_aggr))

# alphabetical order according to site
all_d_aggr <- all_d_aggr %>%
  arrange(SITE)

all_d_aggr <- all_d_aggr %>%
  mutate(CH_METHOD = if_else(CH_METHOD == "auto", "automated", CH_METHOD))

all_d_aggr <- all_d_aggr %>%
  relocate(EC_FCH4_MEAN, EC_FCH4_MEDIAN, .before = CH_FCH4_MEAN)

## save

write.csv(all_d_aggr, "path/daily_aggregation.csv", row.names = FALSE)

## WEEKLY AGGREGATION ##
all_w_aggr <- read.csv("path/allsites_ECfilt_weekaggr_28082025.csv")

all_w_aggr$DELTA_FCH4 <- all_w_aggr$EC_FCH4_comb_median - all_w_aggr$ch_FCH4_nmolCH4m2s1_median

# remove extra columns not used in the analyses

all_w_aggr <- subset(all_w_aggr, select = c(SITE, Year, Week_of_year, Month, ch_FCH4_nmolCH4m2s1_mean,
                                            ch_FCH4_nmolCH4m2s1_median, EC_FCH4_comb_mean, 
                                            EC_FCH4_comb_median, DELTA_FCH4, 
                                            SITE_TS_TOP, SITE_WTD, EC_NEE_F_comb_mean, EC_USTAR_median, 
                                            EC_VPD_F_median, EC_PA_median, EC_PA_F_median, u.wind.F_mean, v.wind.F_mean, DOM_VEG,
                                            ch_method, SITE_CLASSIFICATION
                                            
))

# rename columns

all_w_aggr <- all_w_aggr %>% dplyr::rename("SITE_WTL" = "SITE_WTD",
                                           "EC_FCH4_mean" = "EC_FCH4_comb_mean",
                                           "EC_FCH4_median" = "EC_FCH4_comb_median",
                                           "u_wind_F_mean" = "u.wind.F_mean",
                                           "v_wind_F_mean" = "v.wind.F_mean",
                                           "DOMINANT_VEGETATION" = "DOM_VEG",
                                           "ch_FCH4_mean" = "ch_FCH4_nmolCH4m2s1_mean",
                                           "ch_FCH4_median" = "ch_FCH4_nmolCH4m2s1_median",
                                           "EC_NEE_F" = "EC_NEE_F_comb_mean"
                                           
)

all_w_aggr <- all_w_aggr %>%
  rename_with(
    ~ gsub("^EC_", "", .x),       # remove "EC_" from the start
    .cols = !contains("EC_FCH4")  # but skip columns that contain "EC_FCH4"
  )

names(all_w_aggr) <- toupper(names(all_w_aggr))

# alphabetical order according to site
all_w_aggr <- all_w_aggr %>%
  arrange(SITE)

all_w_aggr <- all_w_aggr %>%
  mutate(CH_METHOD = if_else(CH_METHOD == "auto", "automated", CH_METHOD))

all_w_aggr <- all_w_aggr %>%
  relocate(EC_FCH4_MEAN, EC_FCH4_MEDIAN, .before = CH_FCH4_MEAN)

all_w_aggr <- all_w_aggr %>% dplyr::rename("NEE_F_MEAN" = "NEE_F"
                                           
)

## save

write.csv(all_w_aggr, "path/weekly_aggregation.csv", row.names = FALSE)

## MONTHLY AGGREGATION ##
all_m_aggr <- read.csv("path/allsites_ECfilt_monthaggr_28082025.csv")

all_m_aggr$DELTA_FCH4 <- all_m_aggr$EC_FCH4_comb_median - all_m_aggr$ch_FCH4_nmolCH4m2s1_median

# remove extra columns not used in the analyses

all_m_aggr <- subset(all_m_aggr, select = c(SITE, Year, Month, ch_FCH4_nmolCH4m2s1_mean,
                                            ch_FCH4_nmolCH4m2s1_median, EC_FCH4_comb_mean, 
                                            EC_FCH4_comb_median, DELTA_FCH4, 
                                            SITE_TS_TOP, SITE_WTD, EC_NEE_F_comb_mean, EC_USTAR_median, 
                                            EC_VPD_F_median, EC_PA_median, EC_PA_F_median, u.wind.F_mean, v.wind.F_mean, DOM_VEG,
                                            ch_method, SITE_CLASSIFICATION
                                            
))

# rename columns

all_m_aggr <- all_m_aggr %>% dplyr::rename("SITE_WTL" = "SITE_WTD",
                                           "EC_FCH4_mean" = "EC_FCH4_comb_mean",
                                           "EC_FCH4_median" = "EC_FCH4_comb_median",
                                           "u_wind_F_mean" = "u.wind.F_mean",
                                           "v_wind_F_mean" = "v.wind.F_mean",
                                           "DOMINANT_VEGETATION" = "DOM_VEG",
                                           "ch_FCH4_mean" = "ch_FCH4_nmolCH4m2s1_mean",
                                           "ch_FCH4_median" = "ch_FCH4_nmolCH4m2s1_median",
                                           "EC_NEE_F_mean" = "EC_NEE_F_comb_mean"
                                           
)

all_m_aggr <- all_m_aggr %>%
  rename_with(
    ~ gsub("^EC_", "", .x),       # remove "EC_" from the start
    .cols = !contains("EC_FCH4")  # but skip columns that contain "EC_FCH4"
  )

names(all_m_aggr) <- toupper(names(all_m_aggr))

# alphabetical order according to site
all_m_aggr <- all_m_aggr %>%
  arrange(SITE)

all_m_aggr <- all_m_aggr %>%
  mutate(CH_METHOD = if_else(CH_METHOD == "auto", "automated", CH_METHOD))

all_m_aggr <- all_m_aggr %>%
  relocate(EC_FCH4_MEAN, EC_FCH4_MEDIAN, .before = CH_FCH4_MEAN)

all_m_aggr <- all_m_aggr %>% drop_na(CH_FCH4_MEDIAN)

## save

write.csv(all_m_aggr, "path/monthly_aggregation.csv", row.names = FALSE)

## ANNUAL AGGREGATION ##
all_y_aggr <- read.csv("path/allsites_ECfilt_yraggr_28082025.csv")

all_y_aggr$DELTA_FCH4 <- all_y_aggr$EC_FCH4_comb_median - all_y_aggr$ch_FCH4_nmolCH4m2s1_median

# remove extra columns not used in the analyses

all_y_aggr <- subset(all_y_aggr, select = c(SITE, Year, ch_FCH4_nmolCH4m2s1_mean,
                                            ch_FCH4_nmolCH4m2s1_median, EC_FCH4_comb_mean, 
                                            EC_FCH4_comb_median, DELTA_FCH4, 
                                            SITE_TS_TOP, SITE_WTD, EC_NEE_F_comb_mean, EC_USTAR_median, 
                                            EC_VPD_F_median, EC_PA_median, EC_PA_F_median, u.wind.F_mean, v.wind.F_mean, DOM_VEG,
                                            ch_method, SITE_CLASSIFICATION
                                            
))

# rename columns

all_y_aggr <- all_y_aggr %>% dplyr::rename("SITE_WTL" = "SITE_WTD",
                                           "EC_FCH4_mean" = "EC_FCH4_comb_mean",
                                           "EC_FCH4_median" = "EC_FCH4_comb_median",
                                           "u_wind_F_mean" = "u.wind.F_mean",
                                           "v_wind_F_mean" = "v.wind.F_mean",
                                           "DOMINANT_VEGETATION" = "DOM_VEG",
                                           "ch_FCH4_mean" = "ch_FCH4_nmolCH4m2s1_mean",
                                           "ch_FCH4_median" = "ch_FCH4_nmolCH4m2s1_median",
                                           "EC_NEE_F_mean" = "EC_NEE_F_comb_mean"
                                           
)

all_y_aggr <- all_y_aggr %>%
  rename_with(
    ~ gsub("^EC_", "", .x),       # remove "EC_" from the start
    .cols = !contains("EC_FCH4")  # but skip columns that contain "EC_FCH4"
  )

names(all_y_aggr) <- toupper(names(all_y_aggr))

# alphabetical order according to site
all_y_aggr <- all_y_aggr %>%
  arrange(SITE)

all_y_aggr <- all_y_aggr %>%
  mutate(CH_METHOD = if_else(CH_METHOD == "auto", "automated", CH_METHOD))

all_y_aggr <- all_y_aggr %>%
  relocate(EC_FCH4_MEAN, EC_FCH4_MEDIAN, .before = CH_FCH4_MEAN)

## fix US-LA2

all_y_aggr <- all_y_aggr %>%
  mutate(
    DOMINANT_VEGETATION = if_else(SITE == "US-LA2" & is.na(DOMINANT_VEGETATION),
                                  "AERENCHYMATOUS", DOMINANT_VEGETATION),
    SITE_CLASSIFICATION = if_else(SITE == "US-LA2" & is.na(SITE_CLASSIFICATION),
                                  "marsh", SITE_CLASSIFICATION)
  )

## save

write.csv(all_y_aggr, "path/annual_aggregation.csv", row.names = FALSE)








