#############################################################
### Code for: 
### Määttä et al. (2026) A cross-site comparison of ecosystem- and plot-scale methane fluxes across multiple sites
### Code created by Tiia Määttä, with parts written with GPT 4, 4o and 5.4

#############################################################
###-----------CHAMBER CH4 FLUX DATA CLEANING--------------###
#############################################################

# Remove all variables
rm(list=ls(all = TRUE))
# Clear the console window
cat("\014")
# Close all graphics windows
graphics.off()

# open libraries
library(dplyr)
library(readxl)
library(rio)
library(lubridate)
library(stringr)
library(tidyr)
library(sqldf)

#### CN-HGU

# CHAMBER 1

# read csv

CNHGU_d_C1_F <- read.csv("path")

head(CNHGU_d_C1_F)
tail(CNHGU_d_C1_F)

# rename some columns
colnames(CNHGU_d_C1_F)[1] <- "ch_ID"
colnames(CNHGU_d_C1_F)[2] <- "Date"
colnames(CNHGU_d_C1_F)[3] <- "Time"
colnames(CNHGU_d_C1_F)[4] <- "ch_FCH4_nmolm2s1"

# delete column X
CNHGU_d_C1_F <- subset(CNHGU_d_C1_F, select=-c(X))
head(CNHGU_d_C1_F)

# CHAMBER 2

# read csv

CNHGU_d_C2_F <- read.csv("path")

head(CNHGU_d_C2_F)
tail(CNHGU_d_C2_F)

# rename some columns
colnames(CNHGU_d_C2_F)[1] <- "ch_ID"
colnames(CNHGU_d_C2_F)[2] <- "Date"
colnames(CNHGU_d_C2_F)[3] <- "Time"
colnames(CNHGU_d_C2_F)[4] <- "ch_FCH4_nmolm2s1"

# CHAMBER 3

# read csv

CNHGU_d_C3_F <- read.csv("path")

head(CNHGU_d_C3_F)
tail(CNHGU_d_C3_F)

# rename some columns
colnames(CNHGU_d_C3_F)[1] <- "ch_ID"
colnames(CNHGU_d_C3_F)[2] <- "Date"
colnames(CNHGU_d_C3_F)[3] <- "Time"
colnames(CNHGU_d_C3_F)[4] <- "ch_FCH4_nmolm2s1"

# convert datetimes in these dfs
# replace / with - in Date
CNHGU_d_C1_F$Date <- gsub("/", "-", CNHGU_d_C1_F$Date)
head(CNHGU_d_C1_F)
tail(CNHGU_d_C1_F)
CNHGU_d_C1_F$Date <- mdy(CNHGU_d_C1_F$Date)

head(CNHGU_d_C2_F)
CNHGU_d_C2_F$Date <- gsub("/", "-", CNHGU_d_C2_F$Date)
head(CNHGU_d_C2_F)
tail(CNHGU_d_C2_F)
CNHGU_d_C2_F$Date <- mdy(CNHGU_d_C2_F$Date)

head(CNHGU_d_C3_F)
CNHGU_d_C3_F$Date <- gsub("/", "-", CNHGU_d_C3_F$Date)
head(CNHGU_d_C3_F)
tail(CNHGU_d_C3_F)
CNHGU_d_C3_F$Date <- mdy(CNHGU_d_C3_F$Date)

CNHGU_d_C1_F$Time <- hms::as_hms(CNHGU_d_C1_F$Time)
CNHGU_d_C2_F$Time <- hms::as_hms(CNHGU_d_C2_F$Time)
CNHGU_d_C3_F$Time <- hms::as_hms(CNHGU_d_C3_F$Time)

# combine Date and Time to Datetime
# specifying the format
date_format <- "%Y-%m-%d %H:%M:%S"
# combine Date and Time to Datetime
CNHGU_d_C1_F$Datetime <- as.POSIXct(paste(CNHGU_d_C1_F$Date, CNHGU_d_C1_F$Time), format=date_format)
# delete columns Date and Time
CNHGU_d_C1_F <- subset(CNHGU_d_C1_F, select=-c(Date, Time))
head(CNHGU_d_C1_F)
# move to front
CNHGU_d_C1_F<- CNHGU_d_C1_F %>%
  select(Datetime, everything())

# combine Date and Time to Datetime
CNHGU_d_C2_F$Datetime <- as.POSIXct(paste(CNHGU_d_C2_F$Date, CNHGU_d_C2_F$Time), format=date_format)
# delete columns Date and Time
CNHGU_d_C2_F <- subset(CNHGU_d_C2_F, select=-c(Date, Time))
head(CNHGU_d_C2_F)
# move to front
CNHGU_d_C2_F<- CNHGU_d_C2_F %>%
  select(Datetime, everything())

# combine Date and Time to Datetime
CNHGU_d_C3_F$Datetime <- as.POSIXct(paste(CNHGU_d_C3_F$Date, CNHGU_d_C3_F$Time), format=date_format)
# delete columns Date and Time
CNHGU_d_C3_F <- subset(CNHGU_d_C3_F, select=-c(Date, Time))
head(CNHGU_d_C3_F)
# move to front
CNHGU_d_C3_F<- CNHGU_d_C3_F %>%
  select(Datetime, everything())


# combine the dfs
CNHGU_d_F_C1C2 <- rbind(CNHGU_d_C1_F, CNHGU_d_C2_F)
head(CNHGU_d_F_C1C2)
tail(CNHGU_d_F_C1C2)

CNHGU_d_F_all_comb <- rbind(CNHGU_d_F_C1C2, CNHGU_d_C3_F)
head(CNHGU_d_F_all_comb)
tail(CNHGU_d_F_all_comb)

# remove C from chamber IDs
CNHGU_d_F_all_comb$ch_ID <- gsub('C', '', CNHGU_d_F_all_comb$ch_ID)

# add ch_ to Datetime column name
colnames(CNHGU_d_F_all_comb)[1] <- "ch_Datetime"


# soil moisture (note: not included in analyses for the manuscript)

# read csv

CNHGU_d_SM <- read.csv("path")

head(CNHGU_d_SM)
tail(CNHGU_d_SM)

# rename some columns
colnames(CNHGU_d_SM)[1] <- "Datetime"
colnames(CNHGU_d_SM)[2] <- "C1"
colnames(CNHGU_d_SM)[3] <- "C2"
colnames(CNHGU_d_SM)[4] <- "C3"

# delete column Unit
CNHGU_d_SM <- subset(CNHGU_d_SM, select=-c(Unit))
head(CNHGU_d_SM)
tail(CNHGU_d_SM)

# replace / with - in Datetime
CNHGU_d_SM$Datetime <- gsub("/", "-", CNHGU_d_SM$Datetime)

# separate Datetime into Date and Time
CNHGU_d_SM <- CNHGU_d_SM %>% separate(Datetime, c("Date", "Time"), " ")

# convert to Date format
CNHGU_d_SM$Date <- as_date(CNHGU_d_SM$Date)
class(CNHGU_d_SM$Date)

# add :00 at the end of Time
CNHGU_d_SM$Time <- paste(as.character(CNHGU_d_SM$Time),":00",sep="")
# convert to hms format
# has to be first converted to the 2H 0M 0S format 
CNHGU_d_SM$Time <- hms(CNHGU_d_SM$Time) # weird format
CNHGU_d_SM$Time <- hms::hms(CNHGU_d_SM$Time) # this converts it to the correct format

# specifying the format
date_format <- "%Y-%m-%d %H:%M:%S"
# combine Date and Time to Datetime
CNHGU_d_SM$Datetime <- as.POSIXct(paste(CNHGU_d_SM$Date, CNHGU_d_SM$Time), format=date_format)

# delete columns Date and Time
CNHGU_d_SM <- subset(CNHGU_d_SM, select=-c(Date, Time))
head(CNHGU_d_SM)

# move to front
CNHGU_d_SM<- CNHGU_d_SM %>%
  select(Datetime, everything())

# make C1, C2 and C3 into separate dfs
CNHGU_d_SM_C1 <- CNHGU_d_SM[, c("Datetime", "C1")]
head(CNHGU_d_SM_C1)
tail(CNHGU_d_SM_C1)
CNHGU_d_SM_C2 <- CNHGU_d_SM[, c("Datetime", "C2")]
head(CNHGU_d_SM_C2)
tail(CNHGU_d_SM_C2)
CNHGU_d_SM_C3 <- CNHGU_d_SM[, c("Datetime", "C3")]
head(CNHGU_d_SM_C3)
tail(CNHGU_d_SM_C3)

# add ch_ID column in all of them
CNHGU_d_SM_C1$ch_ID <- "1"
CNHGU_d_SM_C2$ch_ID <- "2"
CNHGU_d_SM_C3$ch_ID <- "3"

# rename the chamber columns
colnames(CNHGU_d_SM_C1)[2] <- "ch_SM"
colnames(CNHGU_d_SM_C2)[2] <- "ch_SM"
colnames(CNHGU_d_SM_C3)[2] <- "ch_SM"

# merge dfs
CNHGU_d_SM_comb_1 <- rbind(CNHGU_d_SM_C1, CNHGU_d_SM_C2)
head(CNHGU_d_SM_comb_1)
tail(CNHGU_d_SM_comb_1)

CNHGU_d_SM_comb_2 <- rbind(CNHGU_d_SM_comb_1, CNHGU_d_SM_C3)

# flux data ends at 11:23:50 --> subset SM data to end at 2016-07-31 11:00:00
CNHGU_SM_sub <- subset(CNHGU_d_SM_comb_2, Datetime <= "2016-07-31 11:00:00")
head(CNHGU_SM_sub)
tail(CNHGU_SM_sub)

# save csv
write.csv(CNHGU_SM_sub, "path", row.names = FALSE)

# soil temperature

# read csv

CNHGU_d_ST <- read.csv("path")

head(CNHGU_d_ST)
tail(CNHGU_d_ST)

# rename some columns
colnames(CNHGU_d_ST)[1] <- "Datetime"
colnames(CNHGU_d_ST)[2] <- "C1"
colnames(CNHGU_d_ST)[3] <- "C2"
colnames(CNHGU_d_ST)[4] <- "C3"

# delete column Unit
CNHGU_d_ST <- subset(CNHGU_d_ST, select=-c(Unit))
head(CNHGU_d_ST)
tail(CNHGU_d_ST)

# replace / with - in Datetime
CNHGU_d_ST$Datetime <- gsub("/", "-", CNHGU_d_ST$Datetime)

# separate Datetime into Date and Time
CNHGU_d_ST <- CNHGU_d_ST %>% separate(Datetime, c("Date", "Time"), " ")

# convert to Date format
CNHGU_d_ST$Date <- as_date(CNHGU_d_ST$Date)
class(CNHGU_d_ST$Date)

# add :00 at the end of Time
CNHGU_d_ST$Time <- paste(as.character(CNHGU_d_ST$Time),":00",sep="")
# convert to hms format
# has to be first converted to the 2H 0M 0S format
CNHGU_d_ST$Time <- hms(CNHGU_d_ST$Time) # weird format
CNHGU_d_ST$Time <- hms::hms(CNHGU_d_ST$Time) # this converts it to the correct format

# specifying the format
date_format <- "%Y-%m-%d %H:%M:%S"
# combine Date and Time to Datetime
CNHGU_d_ST$Datetime <- as.POSIXct(paste(CNHGU_d_ST$Date, CNHGU_d_ST$Time), format=date_format)

# delete columns Date and Time
CNHGU_d_ST <- subset(CNHGU_d_ST, select=-c(Date, Time))
head(CNHGU_d_ST)

# move to front
CNHGU_d_ST<- CNHGU_d_ST %>%
  select(Datetime, everything())

# make C1, C2 and C3 into separate dfs
CNHGU_d_ST_C1 <- CNHGU_d_ST[, c("Datetime", "C1")]
head(CNHGU_d_ST_C1)
tail(CNHGU_d_ST_C1)
CNHGU_d_ST_C2 <- CNHGU_d_ST[, c("Datetime", "C2")]
head(CNHGU_d_ST_C2)
tail(CNHGU_d_ST_C2)
CNHGU_d_ST_C3 <- CNHGU_d_ST[, c("Datetime", "C3")]
head(CNHGU_d_ST_C3)
tail(CNHGU_d_ST_C3)

# add ch_ID column in all of them
CNHGU_d_ST_C1$ch_ID <- "1"
CNHGU_d_ST_C2$ch_ID <- "2"
CNHGU_d_ST_C3$ch_ID <- "3"

# rename the chamber columns
colnames(CNHGU_d_ST_C1)[2] <- "ch_ST"
colnames(CNHGU_d_ST_C2)[2] <- "ch_ST"
colnames(CNHGU_d_ST_C3)[2] <- "ch_ST"

# merge dfs
CNHGU_d_ST_comb_1 <- rbind(CNHGU_d_ST_C1, CNHGU_d_ST_C2)
head(CNHGU_d_ST_comb_1)
tail(CNHGU_d_ST_comb_1)

CNHGU_d_ST_comb_2 <- rbind(CNHGU_d_ST_comb_1, CNHGU_d_ST_C3)

# flux data ends at 11:23:50 --> subset SM data to end at 2016-07-31 11:00:00
CNHGU_ST_sub <- subset(CNHGU_d_ST_comb_2, Datetime <= "2016-07-31 11:00:00")
head(CNHGU_ST_sub)
tail(CNHGU_ST_sub)

# save csv
write.csv(CNHGU_ST_sub, "path", row.names = FALSE)

# merge SM and ST

CNHGU_env_vars <- merge(CNHGU_SM_sub, CNHGU_ST_sub, by=c("Datetime", "ch_ID"))
head(CNHGU_env_vars)
tail(CNHGU_env_vars)

# save as csv
write.csv(CNHGU_env_vars, "path", row.names = FALSE)

head(CNHGU_d_F_all_comb)

# add prefix to env vars ID
colnames(CNHGU_env_vars)[2] <- "env_ch_ID"

# combine the dfs
# ch dataset has different timestamps than env vars
# create a column with + 30 mins to env vars df, called Datetime_end

CNHGU_env_vars$Datetime_end <- CNHGU_env_vars$Datetime + minutes(30)
head(CNHGU_env_vars)
tail(CNHGU_env_vars)

# SQL query
# set local environment time zone to UTC so the time stamps will match
#Sys.setenv(TZ = "UTC")

# rename some columns for SQL
colnames(CNHGU_d_SM_C1)[3] <- "env_ch_ID"
colnames(CNHGU_d_C1_F)[1] <- "ch_Datetime"
colnames(CNHGU_d_ST_C1)[3] <- "env_ch_ID"

# combine ST and SM

CNHGU_env_vars_C1 <- merge(CNHGU_d_SM_C1, CNHGU_d_ST_C1, by=c("Datetime", "env_ch_ID"))
head(CNHGU_env_vars_C1)
tail(CNHGU_env_vars_C1)

# combine ch and env

CNHGU_ch_env_C1 <- sqldf('SELECT ch_Datetime, Datetime, Datetime_end, ch_ID, env_ch_ID, ch_FCH4_nmolm2s1, ch_SM, ch_ST
      FROM CNHGU_d_C1_F
      LEFT JOIN CNHGU_env_vars_C1 ON ch_Datetime BETWEEN Datetime and Datetime_end')

head(CNHGU_ch_env_C1)
tail(CNHGU_ch_env_C1)


# Chamber 2

# rename some columns for SQL
colnames(CNHGU_d_SM_C2)[3] <- "env_ch_ID"
colnames(CNHGU_d_C2_F)[1] <- "ch_Datetime"

colnames(CNHGU_d_ST_C2)[3] <- "env_ch_ID"

# combine ST and SM

CNHGU_env_vars_C2 <- merge(CNHGU_d_SM_C2, CNHGU_d_ST_C2, by=c("Datetime", "env_ch_ID"))
head(CNHGU_env_vars_C2)
tail(CNHGU_env_vars_C2)

# SQL query
# set local environment time zone to UTC so the time stamps will match
#Sys.setenv(TZ = "UTC")

CNHGU_ch_env_C2 <- sqldf('SELECT ch_Datetime, Datetime, Datetime_end, ch_ID, env_ch_ID, ch_FCH4_nmolm2s1, ch_SM, ch_ST
      FROM CNHGU_d_C2_F
      LEFT JOIN CNHGU_env_vars_C2 ON ch_Datetime BETWEEN Datetime and Datetime_end')

head(CNHGU_ch_env_C2)
tail(CNHGU_ch_env_C2)


# chamber 3

# rename some columns for SQL
colnames(CNHGU_d_SM_C3)[3] <- "env_ch_ID"
colnames(CNHGU_d_C3_F)[1] <- "ch_Datetime"

colnames(CNHGU_d_ST_C3)[3] <- "env_ch_ID"

# flux data ends at 11:23:50 --> subset SM data to end at 2016-07-31 11:00:00
CNHGU_d_SM_C3 <- subset(CNHGU_d_SM_C3, Datetime <= "2016-07-31 11:00:00")
tail(CNHGU_d_SM_C3)

CNHGU_d_ST_C3 <- subset(CNHGU_d_ST_C3, Datetime <= "2016-07-31 11:00:00")
tail(CNHGU_d_ST_C3)

# combine ST and SM

CNHGU_env_vars_C3 <- merge(CNHGU_d_SM_C3, CNHGU_d_ST_C3, by=c("Datetime", "env_ch_ID"))
head(CNHGU_env_vars_C3)
tail(CNHGU_env_vars_C3)

# combine the dfs
# ch dataset has different timestamps than env vars
# create a column with + 30 mins to env vars df, called Datetime_end

CNHGU_env_vars_C3$Datetime_end <- CNHGU_env_vars_C3$Datetime + minutes(30)
head(CNHGU_env_vars_C3)
tail(CNHGU_env_vars_C3)


# SQL query
# set local environment time zone to UTC so the time stamps will match
#Sys.setenv(TZ = "UTC")

CNHGU_ch_env_C3 <- sqldf('SELECT ch_Datetime, Datetime, Datetime_end, ch_ID, env_ch_ID, ch_FCH4_nmolm2s1, ch_SM, ch_ST
      FROM CNHGU_d_C3_F
      LEFT JOIN CNHGU_env_vars_C3 ON ch_Datetime BETWEEN Datetime and Datetime_end')

head(CNHGU_ch_env_C3)
tail(CNHGU_ch_env_C3)


# merge dfs
CNHGU_d_env_ch_comb_1 <- rbind(CNHGU_ch_env_C1, CNHGU_ch_env_C2)
head(CNHGU_d_env_ch_comb_1)
tail(CNHGU_d_env_ch_comb_1)

CNHGU_d_env_ch_comb_all <- rbind(CNHGU_d_env_ch_comb_1, CNHGU_ch_env_C3)

# clean column names

colnames(CNHGU_d_env_ch_comb_all_sub)[5] <- "ch_ST"
colnames(CNHGU_d_env_ch_comb_all_sub)[2] <- "ch_ID"


CNHGU_d_env_ch_comb_all_sub <- subset(CNHGU_d_env_ch_comb_all, select=-c(ch_ID))

CNHGU_d_env_ch_comb_all_sub <- subset(CNHGU_d_env_ch_comb_all_sub, select=-c(Datetime, Datetime_end))

# add ch_method etc

CNHGU_d_env_ch_comb_all_sub$ch_method <- "auto"

# save as .csv
write.csv(CNHGU_d_env_ch_comb_all_sub, "path", row.names = FALSE)

#### FI-Si2

# import FI-Si2

FI_Si2_tibble <- read_excel("path")

# convert to dataframe

FI_Si2_df <- as.data.frame(FI_Si2_tibble)

# add column for method

FI_Si2_df$method <- "manual"

# rename original chamber ID column
colnames(FI_Si2_df)[3] <- "ch_ID"

# chamge DOY to year-month-day
# add new column for date

FI_Si2_df$Date <- NA

FI_Si2_df <- FI_Si2_df %>% mutate(Date =
                                    case_when(Year == 2012 ~ as.Date(FI_Si2_df$DOY, origin = "2012-01-01"), 
                                              Year == 2013 ~ as.Date(FI_Si2_df$DOY, origin = "2013-01-01"),
                                              Year == 2014 ~ as.Date(FI_Si2_df$DOY, origin = "2014-01-01"))
)

# add Month and Day columns

FI_Si2_df$Month <- month(FI_Si2_df$Date)

FI_Si2_df$Day <- day(FI_Si2_df$Date)

# move to the front

FI_Si2_df <- FI_Si2_df %>% relocate(Date, Month, Day, .before = ch_ID)

# move year before month
FI_Si2_df <- FI_Si2_df %>% relocate(Year, .before = Month)

# delete DOY

FI_Si2_df = subset(FI_Si2_df, select = -c(DOY) )

# rename original chamber ID column
colnames(FI_Si2_df)[8] <- "ch_FCH4"

# delete Tcha

FI_Si2_df = subset(FI_Si2_df, select = -c(Tcha) )

# add Time and Hour columns --> NA

FI_Si2_df$Time <- NA
FI_Si2_df$Hour <- NA

#move Time after Date and Hour after Day
FI_Si2_df <- FI_Si2_df %>% relocate(Time, .after = Date)
FI_Si2_df <- FI_Si2_df %>% relocate(Hour, .after = Day)

# save as .csv

write.csv(FI_Si2_df, "FI_Si2.csv", row.names=FALSE)

# read csv

FISi2_d <- read.csv("path/FI_Si2.csv")

# add LC_class column

FISi2_d <- FISi2_d %>% mutate(LC_class = case_when(
  startsWith(ch_ID, "HHU") ~ "HHU", # high hummock
  startsWith(ch_ID, "HU") ~ "HU", # hummock
  startsWith(ch_ID, "HL") ~ "HL", # high lawn
  startsWith(ch_ID, "L") ~ "L", # lawn
  startsWith(ch_ID, "HO") ~ "HO", # hollow
  startsWith(ch_ID, "MB") ~ "BP",) # bare peat
)

# remove Surface_type because this information is not needed

FISi2_d <- select(FISi2_d, -8)

# add site

FISi2_d$site <- "FI-Si2"

# move site name to the front

FISi2_d<- FISi2_d %>%
  select(site, everything())

# change the names of soil moist and temp 

colnames(FISi2_d)[which(names(FISi2_d) == "T5")] <- "TS_5"
colnames(FISi2_d)[which(names(FISi2_d) == "T15")] <- "TS_15"
colnames(FISi2_d)[which(names(FISi2_d) == "T30")] <- "TS_30"

# add columns for each available TS across the data sets for each site

FISi2_d$TS_2 <- NA
FISi2_d$TS_10 <- NA

# move TS_2 columns after WT

FISi2_d <- FISi2_d %>% relocate(TS_2, .before = TS_5)

FISi2_d <- FISi2_d %>% relocate(TS_10, .before = TS_15)

# 2012 was leap year --> DOY 191 = 9.7.2012, not 10.7.2012
# subtract 1 day from all 2012 dates

FISi2_d$Date <- as_date(FISi2_d$Date)
FISi2_d$Date[year(FISi2_d$Date) == 2012] <- FISi2_d$Date[year(FISi2_d$Date) == 2012] - 1

FISi2_d <- arrange(FISi2_d, Date)
## save as .csv

write.csv(FISi2_d, "path", row.names = FALSE)

#### SE-DEG ####

# includes multiple Excel sheets for different ch4 plots

# for reading multiple sheets at once:
# specifying the path name
path <- "path/file_name_2015.xlsx"

# reading data from all sheets
data <- import_list(path)

# print data
print (data)

class(data) # list

# extract the dfs from the data list and convert them to dfs

SE_DEG_ch1_15 <- as.data.frame(data["ch1"])

# remove the "ch1." from the column names

names(SE_DEG_ch1_15) <- sub('^ch1.', '', names(SE_DEG_ch1_15))

# rename "chamber." to ch_ID

names(SE_DEG_ch1_15)[names(SE_DEG_ch1_15) == "chamber."] <- "ch_ID"


SE_DEG_ch2_15 <- as.data.frame(data["ch2"])

# remove the "ch2." from the column names

names(SE_DEG_ch2_15) <- sub('^ch2.', '', names(SE_DEG_ch2_15))

# rename "chamber." to ch_ID

names(SE_DEG_ch2_15)[names(SE_DEG_ch2_15) == "chamber."] <- "ch_ID"


SE_DEG_ch4_15 <- as.data.frame(data["ch4"])

# remove the "ch2." from the column names

names(SE_DEG_ch4_15) <- sub('^ch4.', '', names(SE_DEG_ch4_15))

# rename "chamber." to ch_ID

names(SE_DEG_ch4_15)[names(SE_DEG_ch4_15) == "chamber."] <- "ch_ID"



SE_DEG_ch5_15 <- as.data.frame(data["ch5"])

# remove the "ch5." from the column names

names(SE_DEG_ch5_15) <- sub('^ch5.', '', names(SE_DEG_ch5_15))

# rename "chamber." to ch_ID

names(SE_DEG_ch5_15)[names(SE_DEG_ch5_15) == "chamber."] <- "ch_ID"


SE_DEG_ch7_15 <- as.data.frame(data["ch7"])

# remove the "ch2." from the column names

names(SE_DEG_ch7_15) <- sub('^ch7.', '', names(SE_DEG_ch7_15))

# rename "chamber." to ch_ID

names(SE_DEG_ch7_15)[names(SE_DEG_ch7_15) == "chamber."] <- "ch_ID"




SE_DEG_ch8_15 <- as.data.frame(data["ch8"])

# remove the "ch8." from the column names

names(SE_DEG_ch8_15) <- sub('^ch8.', '', names(SE_DEG_ch8_15))

# rename "chamber." to ch_ID

names(SE_DEG_ch8_15)[names(SE_DEG_ch8_15) == "chamber."] <- "ch_ID"


SE_DEG_ch10_15 <- as.data.frame(data["ch10"])

# remove the "ch2." from the column names

names(SE_DEG_ch10_15) <- sub('^ch10.', '', names(SE_DEG_ch10_15))

# rename "chamber." to ch_ID

names(SE_DEG_ch10_15)[names(SE_DEG_ch10_15) == "chamber."] <- "ch_ID"




SE_DEG_ch11_15 <- as.data.frame(data["ch11"])

# remove the "ch11." from the column names

names(SE_DEG_ch11_15) <- sub('^ch11.', '', names(SE_DEG_ch11_15))

# rename "chamber." to ch_ID

names(SE_DEG_ch11_15)[names(SE_DEG_ch11_15) == "chamber."] <- "ch_ID"

# combine the ch-specific dfs into one df

SE_DEG_15 <- bind_rows(SE_DEG_ch1_15, SE_DEG_ch2_15, SE_DEG_ch4_15, 
                       SE_DEG_ch5_15, SE_DEG_ch7_15, SE_DEG_ch8_15, 
                       SE_DEG_ch10_15, SE_DEG_ch11_15)

# ---------------------------------
# then do the same for the 2016 data

# specifying the path name
path2 <- "path/file_name_2016.xlsx"

# reading data from all sheets
data2 <- import_list(path2)

# print data
print(data2)

class(data2) # list

# extract the dfs from the data list and convert them to dfs

SE_DEG_ch1_16 <- as.data.frame(data2["ch1"])

# remove the "ch1." from the column names

names(SE_DEG_ch1_16) <- sub('^ch1.', '', names(SE_DEG_ch1_16))

# rename "chamber." to ch_ID

names(SE_DEG_ch1_16)[names(SE_DEG_ch1_16) == "chamber."] <- "ch_ID"

# check the class of ts_2

typeof(SE_DEG_ch1_16$Ts_2) # double
#change to character

SE_DEG_ch1_16$Ts_2 <- as.character(SE_DEG_ch1_16$Ts_2) #now character

# check the class of ts_10

typeof(SE_DEG_ch1_16$Ts_10) # double
#change to character

SE_DEG_ch1_16$Ts_10 <- as.character(SE_DEG_ch1_16$Ts_10) #now character



SE_DEG_ch2_16 <- as.data.frame(data2["ch2"])

# remove the "ch2." from the column names

names(SE_DEG_ch2_16) <- sub('^ch2.', '', names(SE_DEG_ch2_16))

# rename "chamber." to ch_ID

names(SE_DEG_ch2_16)[names(SE_DEG_ch2_16) == "chamber."] <- "ch_ID"

# check the class of ts_2

typeof(SE_DEG_ch2_16$Ts_2) # double
#change to character

SE_DEG_ch2_16$Ts_2 <- as.character(SE_DEG_ch2_16$Ts_2) #now character

# check the class of ts_10

typeof(SE_DEG_ch2_16$Ts_10) # double
#change to character

SE_DEG_ch2_16$Ts_10 <- as.character(SE_DEG_ch2_16$Ts_10) #now character


SE_DEG_ch4_16 <- as.data.frame(data2["ch4"])

# remove the "ch4." from the column names

names(SE_DEG_ch4_16) <- sub('^ch4.', '', names(SE_DEG_ch4_16))

# rename "chamber." to ch_ID

names(SE_DEG_ch4_16)[names(SE_DEG_ch4_16) == "chamber."] <- "ch_ID"

# check the class of ts_2

typeof(SE_DEG_ch4_16$Ts_2) # double
#change to character

SE_DEG_ch4_16$Ts_2 <- as.character(SE_DEG_ch4_16$Ts_2) #now character

# check the class of ts_10

typeof(SE_DEG_ch4_16$Ts_10) # double
#change to character

SE_DEG_ch4_16$Ts_10 <- as.character(SE_DEG_ch4_16$Ts_10) #now character


SE_DEG_ch5_16 <- as.data.frame(data2["ch5"])


# remove the "ch5." from the column names

names(SE_DEG_ch5_16) <- sub('^ch5.', '', names(SE_DEG_ch5_16))

# rename "chamber." to ch_ID

names(SE_DEG_ch5_16)[names(SE_DEG_ch5_16) == "chamber."] <- "ch_ID"

#check the type of Ts_2

typeof(SE_DEG_ch5_16$Ts_2) # "double"
#change to character

SE_DEG_ch5_16$Ts_2 <- as.character(SE_DEG_ch5_16$Ts_2) #now character

# check the class of ts_10

typeof(SE_DEG_ch5_16$Ts_10) # double
#change to character

SE_DEG_ch5_16$Ts_10 <- as.character(SE_DEG_ch5_16$Ts_10) #now character


SE_DEG_ch7_16 <- as.data.frame(data2["ch7"])

# remove the "ch2." from the column names

names(SE_DEG_ch7_16) <- sub('^ch7.', '', names(SE_DEG_ch7_16))

# rename "chamber." to ch_ID

names(SE_DEG_ch7_16)[names(SE_DEG_ch7_16) == "chamber."] <- "ch_ID"

# check the class of ts_2

typeof(SE_DEG_ch7_16$Ts_2) # double
#change to character

SE_DEG_ch7_16$Ts_2 <- as.character(SE_DEG_ch7_16$Ts_2) #now character

# check the class of ts_10

typeof(SE_DEG_ch7_16$Ts_10) # double
#change to character

SE_DEG_ch7_16$Ts_10 <- as.character(SE_DEG_ch7_16$Ts_10) #now character



SE_DEG_ch8_16 <- as.data.frame(data2["ch8"])

# remove the "ch8." from the column names

names(SE_DEG_ch8_16) <- sub('^ch8.', '', names(SE_DEG_ch8_16))

# rename "chamber." to ch_ID

names(SE_DEG_ch8_16)[names(SE_DEG_ch8_16) == "chamber."] <- "ch_ID"

#check the type of Ts_2

typeof(SE_DEG_ch8_16$Ts_2) # "double"
#change to character

SE_DEG_ch8_16$Ts_2 <- as.character(SE_DEG_ch8_16$Ts_2) #now character

# check the class of ts_10

typeof(SE_DEG_ch8_16$Ts_10) # double
#change to character

SE_DEG_ch8_16$Ts_10 <- as.character(SE_DEG_ch8_16$Ts_10) #now character



SE_DEG_ch10_16 <- as.data.frame(data2["ch10"])

# remove the "ch2." from the column names

names(SE_DEG_ch10_16) <- sub('^ch10.', '', names(SE_DEG_ch10_16))

# rename "chamber." to ch_ID

names(SE_DEG_ch10_16)[names(SE_DEG_ch10_16) == "chamber."] <- "ch_ID"

# check the class of ts_2

typeof(SE_DEG_ch10_16$Ts_2) # double
#change to character

SE_DEG_ch10_16$Ts_2 <- as.character(SE_DEG_ch10_16$Ts_2) #now character

# check the class of ts_10

typeof(SE_DEG_ch10_16$Ts_10) # double
#change to character

SE_DEG_ch10_16$Ts_10 <- as.character(SE_DEG_ch10_16$Ts_10) #now character




SE_DEG_ch11_16 <- as.data.frame(data2["ch11"])

# remove the "ch11." from the column names

names(SE_DEG_ch11_16) <- sub('^ch11.', '', names(SE_DEG_ch11_16))

# rename "chamber." to ch_ID

names(SE_DEG_ch11_16)[names(SE_DEG_ch11_16) == "chamber."] <- "ch_ID"

#check the type of Ts_2

typeof(SE_DEG_ch11_16$Ts_2) # "double"
#change to character

SE_DEG_ch11_16$Ts_2 <- as.character(SE_DEG_ch11_16$Ts_2) #now character
# check the class of ts_10

typeof(SE_DEG_ch11_16$Ts_10) # double
#change to character

SE_DEG_ch11_16$Ts_10 <- as.character(SE_DEG_ch11_16$Ts_10) #now character

# combine the ch-specific dfs into one df

SE_DEG_16 <- bind_rows(SE_DEG_ch1_16, SE_DEG_ch2_16, SE_DEG_ch4_16, 
                       SE_DEG_ch5_16, SE_DEG_ch7_16,
                       SE_DEG_ch8_16, SE_DEG_ch10_16,
                       SE_DEG_ch11_16)

# combine SE_DEG_15 and SE_DEG_16

SE_DEG_df <- rbind(SE_DEG_15, SE_DEG_16)

# add column for method (manual or auto)

SE_DEG_df$method <- "auto"

# save as .csv

write.csv(SE_DEG_df, "SE_DEG.csv", row.names = FALSE)

#### import new SE_DEG.csv

SE_DEG_df <- read.csv("path/SE_DEG.csv")

# add Date and Time

SE_DEG_df$Date <- as.Date(with(SE_DEG_df,paste(year,month,day,sep="-")),"%Y-%m-%d")

SE_DEG_df$Time <- paste(SE_DEG_df$hour,":",SE_DEG_df$min,":00",sep="")

# set Time to hms format
SE_DEG_df$Time <- hms::as_hms(SE_DEG_df$Time)


# move to the front

SE_DEG_df <- SE_DEG_df %>% relocate(Date, Time, .before = year)

# rename "CH4.flux" to FCH4

names(SE_DEG_df)[names(SE_DEG_df) == "CH4.flux"] <- "ch_FCH4"

# remove min and sec

SE_DEG_df = subset(SE_DEG_df, select = -c(min, sec) )

# rename hour to Hour

# rename "CH4.flux" to FCH4

names(SE_DEG_df)[names(SE_DEG_df) == "hour"] <- "Hour"
names(SE_DEG_df)[names(SE_DEG_df) == "year"] <- "Year"
names(SE_DEG_df)[names(SE_DEG_df) == "month"] <- "Month"
names(SE_DEG_df)[names(SE_DEG_df) == "day"] <- "Day"

# save this version as SE_DEG_stats.csv

write.csv(SE_DEG_df, "SE_DEG_stats.csv", row.names = FALSE)

# for the other version, remove the stats

SE_DEG_df = subset(SE_DEG_df, select = -c(CH4.R2, CH4.RMSE, CH4_PAR_std, CH4_Ta_mean, CH4_Ta_std) )

# save as .csv
write.csv(SE_DEG_df, "path/SE_DEG.csv", row.names=FALSE)

# read again
SEDEG_d <- read.csv("path/SE_DEG.csv")

# add land cover class

SEDEG_d$LC_class <- NA

# add site

SEDEG_d$site <- "SE-DEG"

# move site name to the front

SEDEG_d <- SEDEG_d %>%
  select(site, everything())

# change the names of soil moist and temp

colnames(SEDEG_d)[which(names(SEDEG_d) == "Ts_2")] <- "TS_2"
colnames(SEDEG_d)[which(names(SEDEG_d) == "Ts_10")] <- "TS_10"

# add columns for each available TS across the data sets for each site

SEDEG_d$TS_5 <- NA
SEDEG_d$TS_15 <- NA
SEDEG_d$TS_30 <- NA

# move TS columns to the right spot

SEDEG_d <- SEDEG_d %>% relocate(TS_5, .before = TS_10)
SEDEG_d <- SEDEG_d %>% relocate(TS_15, TS_30, .after = TS_10)

# combine Date and Time to datetime

# specifying the format
date_format <- "%Y-%m-%d %H:%M:%S"

# combining date and time into single object
SEDEG_d$Datetime <- as.POSIXct(paste(SEDEG_d$Date, SEDEG_d$Time), format=date_format)

# remove some columns that are not needed

SEDEG_d <- subset(SEDEG_d, select=-c(Date, Time, Year, Month, Day, Hour))

# convert to datetime
SEDEG_d$Datetime <- as_datetime(SEDEG_d$Datetime, tz = "UTC")

# move to the right spot
SEDEG_d <- SEDEG_d %>% relocate(Datetime, .before = ch_ID)

# add ch_ prefix to all columns

SEDEG_d <- SEDEG_d %>% rename_with( ~ paste0("ch_", .x))

# rename SITE, ch_ch_FCH4

colnames(SEDEG_d)[3] <- "ch_ID"
colnames(SEDEG_d)[4] <- "ch_FCH4_nmol_or_umolCH4m2s1"
colnames(SEDEG_d)[1] <- "SITE"

## save as .csv

write.csv(SEDEG_d, "path", row.names = FALSE)

#### US-HO1 ####

# import US-HO1

US_HO1_tibble <- read.csv("path")

# convert to dataframe

US_HO1_df <- as.data.frame(US_HO1_tibble)

# replace -9999 values with NA
US_HO1_df <- US_HO1_df %>% na_if(-9999)

# add column for method (manual or auto)

US_HO1_df$method <- "auto"
US_HO1_df

# drop the time_start and time_end columns

US_HO1_df = subset(US_HO1_df, select = -c(CMB_TIME_START_2_1_1, CMB_TIME_END_2_1_1,
                                          CMB_TIME_START_4_1_1, CMB_TIME_END_4_1_1,
                                          CMB_TIME_START_6_1_1, CMB_TIME_END_6_1_1,
                                          CMB_TIME_START_7_1_1, CMB_TIME_END_7_1_1,
                                          CMB_TIME_START_8_1_1, CMB_TIME_END_8_1_1,
                                          CMB_TIME_START_9_1_1, CMB_TIME_END_9_1_1,
                                          CMB_TIME_START_10_1_1, CMB_TIME_END_10_1_1,
                                          CMB_TIME_START_11_1_1, CMB_TIME_END_11_1_1,
                                          CMB_TIME_START_12_1_1, CMB_TIME_END_12_1_1,
                                          CMB_TIME_START_13_1_1, CMB_TIME_END_13_1_1,
                                          CMB_TIME_START_14_1_1, CMB_TIME_END_14_1_1,
                                          CMB_TIME_START_15_1_1, CMB_TIME_END_15_1_1,
                                          CMB_TIME_START_16_1_1, CMB_TIME_END_16_1_1,
                                          CMB_TIME_START_17_1_1, CMB_TIME_END_17_1_1,
                                          CMB_TIME_START_18_1_1, CMB_TIME_END_18_1_1,
                                          CMB_TIME_START_19_1_1, CMB_TIME_END_19_1_1,
                                          CMB_TIME_START_20_1_1, CMB_TIME_END_20_1_1,
                                          CMB_TIME_START_21_1_1, CMB_TIME_END_21_1_1,
                                          CMB_TIME_START_22_1_1, CMB_TIME_END_22_1_1,
                                          CMB_TIME_START_23_1_1, CMB_TIME_END_23_1_1) )

# create new column for chamber IDs, flux and soil temp

US_HO1_df$ch_ID <- NA
US_HO1_df$FCH4 <- NA
US_HO1_df$TS <- NA

#convert the timestamp to non-scientific format

US_HO1_df$Timestamp_Start <- format(US_HO1_df$Timestamp_Start, scientific = FALSE)
US_HO1_df$Timestamp_End <- format(US_HO1_df$Timestamp_End, scientific = FALSE)

# assigning new names to the columns of the data frame
colnames(US_HO1_df) <- c('Timestamp_Start','Timestamp_End','2', '4', '6', '7', '8', '9', '10', 
                         '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22',
                         '23','SWC_2', 'SWC_4', 'SWC_6', 'SWC_20', 'SWC_21', 'SWC_22', 'TS_2', 'TS_4',
                         'TS_6', 'TS_20', 'TS_21', 'TS_22', 'method', 'ch_ID', 'FCH4', 'TS')

# choose only the flux data into flux df
US_HO1_flux <- select(US_HO1_df, "2","4","6","7","8","9","10","11","12","13","14","15",
                      "16","17","18","19","20","22","23")

# choose only SWC data
US_HO1_SWC <- select(US_HO1_df, "SWC_2", "SWC_4", "SWC_6", "SWC_20", "SWC_21", "SWC_22")

# choose only TS data
US_HO1_TS <- select(US_HO1_df, "TS_2", "TS_4", "TS_6", "TS_20", "TS_21", "TS_22")

# choose only the timestamps
US_HO1_time <- select(US_HO1_df, "Timestamp_Start", "Timestamp_End")

###### choose by chamber

US_HO1_flux_2 <- select(US_HO1_df, "Timestamp_Start", "Timestamp_End","2")

# add ch_ID column

US_HO1_flux_2$ch_ID <- 2

# rename original chamber ID column
colnames(US_HO1_flux_2)[3] <- "FCH4"


# 4

US_HO1_flux_4 <- select(US_HO1_df,"Timestamp_Start", "Timestamp_End","4")

# add ch_ID column

US_HO1_flux_4$ch_ID <- 4

# rename original chamber ID column
colnames(US_HO1_flux_4)[3] <- "FCH4"

# 6

US_HO1_flux_6 <- select(US_HO1_df, "Timestamp_Start", "Timestamp_End","6")

# add ch_ID column

US_HO1_flux_6$ch_ID <- 6

# rename original chamber ID column
colnames(US_HO1_flux_6)[3] <- "FCH4"

# 7

US_HO1_flux_7 <- select(US_HO1_df, "Timestamp_Start", "Timestamp_End","7")

# add ch_ID column

US_HO1_flux_7$ch_ID <- 7

# rename original chamber ID column
colnames(US_HO1_flux_7)[3] <- "FCH4"

# 8

US_HO1_flux_8 <- select(US_HO1_df, "Timestamp_Start", "Timestamp_End","8")

# add ch_ID column

US_HO1_flux_8$ch_ID <- 8

# rename original chamber ID column
colnames(US_HO1_flux_8)[3] <- "FCH4"

# 9

US_HO1_flux_9 <- select(US_HO1_df, "Timestamp_Start", "Timestamp_End","9")

# add ch_ID column

US_HO1_flux_9$ch_ID <- 9

# rename original chamber ID column
colnames(US_HO1_flux_9)[3] <- "FCH4"

# 10

US_HO1_flux_10 <- select(US_HO1_df, "Timestamp_Start", "Timestamp_End","10")

# add ch_ID column

US_HO1_flux_10$ch_ID <- 10

# rename original chamber ID column
colnames(US_HO1_flux_10)[3] <- "FCH4"

# 11

US_HO1_flux_11 <- select(US_HO1_df, "Timestamp_Start", "Timestamp_End","11")

# add ch_ID column

US_HO1_flux_11$ch_ID <- 11

# rename original chamber ID column
colnames(US_HO1_flux_11)[3] <- "FCH4"

# 12

US_HO1_flux_12 <- select(US_HO1_df, "Timestamp_Start", "Timestamp_End","12")

# add ch_ID column

US_HO1_flux_12$ch_ID <- 12

# rename original chamber ID column
colnames(US_HO1_flux_12)[3] <- "FCH4"

# 13

US_HO1_flux_13 <- select(US_HO1_df, "Timestamp_Start", "Timestamp_End","13")

# add ch_ID column

US_HO1_flux_13$ch_ID <- 13

# rename original chamber ID column
colnames(US_HO1_flux_13)[3] <- "FCH4"

# 14

US_HO1_flux_14 <- select(US_HO1_df, "Timestamp_Start", "Timestamp_End","14")

# add ch_ID column

US_HO1_flux_14$ch_ID <- 14

# rename original chamber ID column
colnames(US_HO1_flux_14)[3] <- "FCH4"

# 15

US_HO1_flux_15 <- select(US_HO1_df, "Timestamp_Start", "Timestamp_End","15")

# add ch_ID column

US_HO1_flux_15$ch_ID <- 15

# rename original chamber ID column
colnames(US_HO1_flux_15)[3] <- "FCH4"


# 16

US_HO1_flux_16 <- select(US_HO1_df, "Timestamp_Start", "Timestamp_End","16")

# add ch_ID column

US_HO1_flux_16$ch_ID <- 16

# rename original chamber ID column
colnames(US_HO1_flux_16)[3] <- "FCH4"

# 17

US_HO1_flux_17 <- select(US_HO1_df, "Timestamp_Start", "Timestamp_End","17")

# add ch_ID column

US_HO1_flux_17$ch_ID <- 17

# rename original chamber ID column
colnames(US_HO1_flux_17)[3] <- "FCH4"

# 18

US_HO1_flux_18 <- select(US_HO1_df, "Timestamp_Start", "Timestamp_End","18")

# add ch_ID column

US_HO1_flux_18$ch_ID <- 18

# rename original chamber ID column
colnames(US_HO1_flux_18)[3] <- "FCH4"

# 19

US_HO1_flux_19 <- select(US_HO1_df, "Timestamp_Start", "Timestamp_End","19")

# add ch_ID column

US_HO1_flux_19$ch_ID <- 19

# rename original chamber ID column
colnames(US_HO1_flux_19)[3] <- "FCH4"

# 20

US_HO1_flux_20 <- select(US_HO1_df, "Timestamp_Start", "Timestamp_End","20")

# add ch_ID column

US_HO1_flux_20$ch_ID <- 20

# rename original chamber ID column
colnames(US_HO1_flux_20)[3] <- "FCH4"

# 21

US_HO1_flux_21 <- select(US_HO1_df, "Timestamp_Start", "Timestamp_End","21")

# add ch_ID column

US_HO1_flux_21$ch_ID <- 21

# rename original chamber ID column
colnames(US_HO1_flux_21)[3] <- "FCH4"

# 22

US_HO1_flux_22 <- select(US_HO1_df, "Timestamp_Start", "Timestamp_End","22")

# add ch_ID column

US_HO1_flux_22$ch_ID <- 22

# rename original chamber ID column
colnames(US_HO1_flux_22)[3] <- "FCH4"

# 23

US_HO1_flux_23 <- select(US_HO1_df, "Timestamp_Start", "Timestamp_End","23")

# add ch_ID column

US_HO1_flux_23$ch_ID <- 23

# rename original chamber ID column
colnames(US_HO1_flux_23)[3] <- "FCH4"

# now combine the flux dfs row after row

# combine the ch-specific dfs into one df

US_HO1_flux_comb <- bind_rows(US_HO1_flux_2, US_HO1_flux_4, US_HO1_flux_6, US_HO1_flux_7, US_HO1_flux_8, US_HO1_flux_9, US_HO1_flux_10,
                              US_HO1_flux_11, US_HO1_flux_12, US_HO1_flux_13, US_HO1_flux_14, US_HO1_flux_15, US_HO1_flux_16,
                              US_HO1_flux_17, US_HO1_flux_18, US_HO1_flux_19, US_HO1_flux_20, US_HO1_flux_21,
                              US_HO1_flux_22, US_HO1_flux_23)

###### SWC (soil water content) data ######

# 2
# choose only SWC data
US_HO1_SWC_2 <- select(US_HO1_df, "Timestamp_Start", "Timestamp_End", "SWC_2")

# add ch_ID column

US_HO1_SWC_2$ch_ID <- 2

# rename original chamber ID column
colnames(US_HO1_SWC_2)[3] <- "SWC"

# 4

# choose only SWC data
US_HO1_SWC_4 <- select(US_HO1_df, "Timestamp_Start", "Timestamp_End", "SWC_4")

# add ch_ID column

US_HO1_SWC_4$ch_ID <- 4

# rename original chamber ID column
colnames(US_HO1_SWC_4)[3] <- "SWC"

# 6

# choose only SWC data
US_HO1_SWC_6 <- select(US_HO1_df, "Timestamp_Start", "Timestamp_End", "SWC_6")

# add ch_ID column

US_HO1_SWC_6$ch_ID <- 6

# rename original chamber ID column
colnames(US_HO1_SWC_6)[3] <- "SWC"

# 20

# choose only SWC data
US_HO1_SWC_20 <- select(US_HO1_df, "Timestamp_Start", "Timestamp_End", "SWC_20")

# add ch_ID column

US_HO1_SWC_20$ch_ID <- 20

# rename original chamber ID column
colnames(US_HO1_SWC_20)[3] <- "SWC"

# 21

# choose only SWC data
US_HO1_SWC_21 <- select(US_HO1_df, "Timestamp_Start", "Timestamp_End", "SWC_21")

# add ch_ID column

US_HO1_SWC_21$ch_ID <- 21

# rename original chamber ID column
colnames(US_HO1_SWC_21)[3] <- "SWC"


# 22

# choose only SWC data
US_HO1_SWC_22 <- select(US_HO1_df, "Timestamp_Start", "Timestamp_End", "SWC_22")

# add ch_ID column

US_HO1_SWC_22$ch_ID <- 22

# rename original chamber ID column
colnames(US_HO1_SWC_22)[3] <- "SWC"

# now combine the SWC dfs into one

# combine the ch-specific dfs into one df

US_HO1_SWC_comb <- bind_rows(US_HO1_SWC_2, US_HO1_SWC_4, US_HO1_SWC_6, US_HO1_SWC_20, US_HO1_SWC_21,
                             US_HO1_SWC_22)

######### now the same for TS (soil temperature) ##########

# 2
# choose only SWC data
US_HO1_TS_2 <- select(US_HO1_df, "Timestamp_Start", "Timestamp_End", "TS_2")

# add ch_ID column

US_HO1_TS_2$ch_ID <- 2

# rename original chamber ID column
colnames(US_HO1_TS_2)[3] <- "TS"

# 4

# choose only SWC data
US_HO1_TS_4 <- select(US_HO1_df, "Timestamp_Start", "Timestamp_End", "TS_4")

# add ch_ID column

US_HO1_TS_4$ch_ID <- 4

# rename original chamber ID column
colnames(US_HO1_TS_4)[3] <- "TS"

# 6

# choose only SWC data
US_HO1_TS_6 <- select(US_HO1_df, "Timestamp_Start", "Timestamp_End", "TS_6")

# add ch_ID column

US_HO1_TS_6$ch_ID <- 6

# rename original chamber ID column
colnames(US_HO1_TS_6)[3] <- "TS"

# 20

# choose only SWC data
US_HO1_TS_20 <- select(US_HO1_df, "Timestamp_Start", "Timestamp_End", "TS_20")

# add ch_ID column

US_HO1_TS_20$ch_ID <- 20

# rename original chamber ID column
colnames(US_HO1_TS_20)[3] <- "TS"

# 21

# choose only SWC data
US_HO1_TS_21 <- select(US_HO1_df, "Timestamp_Start", "Timestamp_End", "TS_21")

# add ch_ID column

US_HO1_TS_21$ch_ID <- 21

# rename original chamber ID column
colnames(US_HO1_TS_21)[3] <- "TS"


# 22

# choose only SWC data
US_HO1_TS_22 <- select(US_HO1_df, "Timestamp_Start", "Timestamp_End", "TS_22")

# add ch_ID column

US_HO1_TS_22$ch_ID <- 22

# rename original chamber ID column
colnames(US_HO1_TS_22)[3] <- "TS"

# now combine the TS dfs into one

# combine the ch-specific dfs into one df

US_HO1_TS_comb <- bind_rows(US_HO1_TS_2, US_HO1_TS_4, US_HO1_TS_6, US_HO1_TS_20, US_HO1_TS_21,
                            US_HO1_TS_22)

###############

# combine fluxes and SWC

US_HO1_flux_SWC <- merge(US_HO1_flux_comb, US_HO1_SWC_comb, by=c("Timestamp_Start", "Timestamp_End", "ch_ID"), all=TRUE)

# add TS df

US_HO1_all <- merge(US_HO1_flux_SWC, US_HO1_TS_comb, by=c("Timestamp_Start", "Timestamp_End", "ch_ID"), all=TRUE)

# add method and site columns

US_HO1_all$method <- "auto"

US_HO1_all$site <- "US-HO1"

US_HO1 <- US_HO1_all

class(US_HO1$Timestamp_Start)
class(US_HO1$Timestamp_End)

#convert the timestamp to non-scientific format

US_HO1$Timestamp_start <- format(US_HO1$Timestamp_start, scientific = FALSE)

# rename FCH4

colnames(US_HO1)[4] <- "ch_FCH4"

# remove Timestamp end

US_HO1 = subset(US_HO1, select = -c(Timestamp_End) )

# rename timestamp

colnames(US_HO1)[1] <- "Datetime"

class(US_HO1$Datetime) #character

US_HO1$Time <- str_sub(US_HO1$Datetime, - 4, - 1)

US_HO1$Hour <- str_sub(US_HO1$Datetime, - 4, - 3)

US_HO1$Min <- str_sub(US_HO1$Datetime, - 2, - 1)

#create new column for Time2 with seconds

US_HO1$Time2 <- paste(as.character(US_HO1$Datetime),"00",sep="")

US_HO1$Time2 <- paste(US_HO1$Hour,US_HO1$Min,"00",sep=":")

#add new column for Sec
US_HO1$Sec <- "00"

#add new column for Year
US_HO1$Year <- str_sub(US_HO1$Datetime, 1, 4)

#add new column for Month
US_HO1$Month <- str_sub(US_HO1$Datetime, 5, 6)

#add new column for Day
US_HO1$Day <- str_sub(US_HO1$Datetime, 7, 8)

# remove Time

US_HO1 = subset(US_HO1, select = -c(Time) )

# remove Min

US_HO1 = subset(US_HO1, select = -c(Min) )

# remove Sec

US_HO1 = subset(US_HO1, select = -c(Sec) )

# rename Time2

colnames(US_HO1)[10] <- "Time"

# convert to hms

US_HO1$Time <- hms::as_hms(US_HO1$Time)

# cut Datetime into date

US_HO1$Date <- str_sub(US_HO1$Datetime, - 12, - 5)
class(US_HO1$Date)
US_HO1$Date <- as.numeric(US_HO1$Date)

# convert to date

US_HO1$Date <- ymd(US_HO1$Date)

# add seconds to datetime

US_HO1$Datetime <- paste(US_HO1$Datetime,"00",sep="")

# convert Datetime to datetime

US_HO1$Datetime <- as_datetime(US_HO1$Time2)

# remove Time2
US_HO1 = subset(US_HO1, select = -c(Time2) )

# rename original datetime
colnames(US_HO1)[1] <- "Datetime_old"

# rename Datetime

colnames(US_HO1)[16] <- "Datetime"

# move columns to front
US_HO1 <- US_HO1 %>% relocate(Datetime, Date, Time, Year, Month, Day, Hour, .before = Datetime_old)

# save as .csv

write.csv(US_HO1, "path/US_HO1.csv", row.names = FALSE)

# read the df again

USHO1_d <- read.csv("path/US_HO1.csv")

# move site name to the front

USHO1_d <- USHO1_d %>%
  select(site, everything())

# remove datetime column as it is not needed

USHO1_d<- select(USHO1_d, -2)

# remove Datetime_old

USHO1_d<- select(USHO1_d, -8)

# add land cover class

USHO1_d <- USHO1_d %>% mutate(LC_class = case_when(
  ch_ID == 2 ~ "upland-control",
  ch_ID == 4 ~ "upland-control",
  ch_ID == 6 ~ "upland-control", 
  ch_ID == 7 ~ "transitional", 
  ch_ID == 8 ~ "transitional",
  ch_ID == 9 ~ "transitional",
  ch_ID == 10 ~ "wetland",
  ch_ID == 11 ~ "wetland",
  ch_ID == 12 ~ "wetland",
  ch_ID == 13 ~ "upland-control",
  ch_ID == 14 ~ "upland-control",
  ch_ID == 15 ~ "wetland",
  ch_ID == 16 ~ "wetland",
  ch_ID == 17 ~ "transitional",
  ch_ID == 18 ~ "transitional",
  ch_ID == 19 ~ "transitional",
  ch_ID == 20 ~ "upland-trenched",
  ch_ID == 21 ~ "upland-trenched",
  ch_ID == 22 ~ "upland-trenched",
  ch_ID == 23 ~ "upland-control",)
)

# change the names of soil moist and temp

colnames(USHO1_d)[which(names(USHO1_d) == "TS")] <- "TS_X" # depth not known
colnames(USHO1_d)[which(names(USHO1_d) == "SWC")] <- "SM_X" # depth not known

# add columns for each available TS across the data sets for each site

USHO1_d$TS_2 <- NA
USHO1_d$TS_5 <- NA
USHO1_d$TS_10 <- NA
USHO1_d$TS_15 <- NA
USHO1_d$TS_30 <- NA

# move TS columns to the right spot

USHO1_d <- USHO1_d %>% relocate(TS_2, TS_5, TS_10, TS_15, TS_30, .after = TS_X)

# add WTL column

USHO1_d$WTL <- NA
USHO1_d <- USHO1_d %>% relocate(WTL, .after = TS_30)

# specifying the format
format <- "%Y-%m-%d %H:%M:%S"

USHO1_d$Date <- as.character(USHO1_d$Date)
USHO1_d$Time <- as.character(USHO1_d$Time)

# combining date and time into single object
USHO1_d$Datetime <- as_datetime(paste(USHO1_d$Date, USHO1_d$Time), format=format)

# move datetime to front

USHO1_d<- USHO1_d %>% relocate(Datetime, .after = site)

# add ch_ prefix to all columns

USHO1_d <- USHO1_d %>% rename_with( ~ paste0("ch_", .x))

# rename SITE, ch_ch_FCH4

colnames(USHO1_d)[3] <- "ch_ID"
colnames(USHO1_d)[4] <- "ch_FCH4_umolCH4m2s1"
colnames(USHO1_d)[1] <- "SITE"

# remove unnecessary columns
USHO1_d <- subset(USHO1_d, select=-c(ch_Date, ch_Time, ch_Year, ch_Month, ch_Day, ch_Hour))

## save as .csv

write.csv(USHO1_d, "path", row.names = FALSE)

#### US-LA1 ####

# multiple sheets for different ch4 plots

# for reading multiple sheets at once:

# specifying the path name
path <- "path"

# reading data from all sheets
data <- import_list(path)

# print data
print (data)

class(data) # list

# extract the dfs from the data list and convert them to dfs

US_LA1 <- as.data.frame(data["US-LA1"])

# remove the "US.LA1." from the column names

names(US_LA1) <- sub('^US.LA1.', '', names(US_LA1))

# rename "Chamber" to ch_ID

names(US_LA1)[names(US_LA1) == "Chamber"] <- "ch_ID"

# remove unnecessary columns

US_LA1 <- subset(US_LA1, select = -c(Site, Stderr...10, CO2.Flux.rate...ug.CO2.C.m2.hr., CO2.Flux.rate...mg.CO2.C.m2.hr., 
                                     Stderr...15, CO2.C.mg.day...16, CO2.C.mg.day...17, N2O.Flux.rate...ug.N20.N.m2.hr.,
                                     N2O.Flux.rate...microg.N2O.N.m2.hr., Stderr...20, N2O.N.microg.day...21, N2O.N.microg.day...22) )
US_LA1 <- subset(US_LA1, select = -c(Chamber.H2O.depth..cm., Bdwlk) ) # bdwlk=boardwalk

US_LA1 <- subset(US_LA1, select = -c(CH4.C.mg.day...12) )

US_LA1 <- subset(US_LA1, select = -c(B1.SoilT, B1.AirT) ) # remove columns in fahrenheit (B1)

US_LA1 <- subset(US_LA1, select = -c(B2.SoilT, B2.AirT) ) 

US_LA1 <- subset(US_LA1, select = -c(ch_FCH4_HH_mg) ) 
US_LA1 <- subset(US_LA1, select = -c(ch_FCH4_DD_mg) ) 

# rename columns

names(US_LA1)[names(US_LA1) == "Chamber.Soil.Temp...10cm...C."] <- "TS_10"
names(US_LA1)[names(US_LA1) == "CH4.Flux.rate...ug.CH4.C.m2.hr."] <- "ch_FCH4_HH_ug"
names(US_LA1)[names(US_LA1) == "CH4.Flux.rate...mg.CH4.C.m2.hr."] <- "ch_FCH4_HH_mg"
names(US_LA1)[names(US_LA1) == "CH4.C.mg.day...11"] <- "ch_FCH4_DD_mg"
names(US_LA1)[names(US_LA1) == "Water.Level"] <- "WTL"
names(US_LA1)[names(US_LA1) == "ch_FCH4_HH_ug"] <- "ch_FCH4"

# convert Date to date format

US_LA1$Date <- ymd(US_LA1$Date)

# add year and time (NA)

US_LA1$Year <- year(US_LA1$Date)

US_LA1$Time <- NA

US_LA1$Hour <- NA

US_LA1$Day <- day(US_LA1$Date)

US_LA1 <- subset(US_LA1, select = -c(Min) ) 

# move columns to front
US_LA1 <- US_LA1 %>% relocate(Time, Year, Month, Day, Hour, .after = Date)

US_LA1$method <- "manual"

# save as .csv

write.csv(US_LA1, "US_LA1.csv", row.names = FALSE)

# read the df again
USLA1_d <- read.csv("path/US_LA1.csv")

# add site

USLA1_d$site <- "US-LA1"

# move site name to the front

USLA1_d <- USLA1_d %>%
  select(site, everything())

# add land cover class

USLA1_d$LC_class <- NA

# add columns for each available TS across the data sets for each site

USLA1_d$TS_2 <- NA
USLA1_d$TS_5 <- NA
USLA1_d$TS_15 <- NA
USLA1_d$TS_30 <- NA

# move TS columns to the right spot

USLA1_d <- USLA1_d %>% relocate(TS_2, TS_5, .before = TS_10)
USLA1_d <- USLA1_d %>% relocate(TS_15, TS_30, .after = TS_10)

# add ch_ prefix to all columns

USLA1_d <- USLA1_d %>% rename_with( ~ paste0("ch_", .x))

# rename SITE, ch_ch_FCH4

colnames(USLA1_d)[8] <- "ch_ID"
colnames(USLA1_d)[14] <- "ch_FCH4_ugCH4m2hr"
colnames(USLA1_d)[1] <- "SITE"

# remove unnecessary columns
USLA1_d <- subset(USLA1_d, select=-c(ch_Time, ch_Year, ch_Month, ch_Day, ch_Hour))

# move the flux column to right spot
USLA1_d <- USLA1_d %>% relocate(ch_FCH4_ugCH4m2hr, .before = ch_TS_2)

## save as .csv

write.csv(USLA1_d, "path", row.names = FALSE)

#### US-LA2 ####

# includes multiple sheets for different ch4 plots

# for reading multiple sheets at once:
# specifying the path name
path <- "path"

# reading data from all sheets
data <- import_list(path)

# print data
print (data)

class(data) # list

# extract the dfs from the data list and convert them to dfs

US_LA2 <- as.data.frame(data["US-LA2"])

# remove the "US.LA2." from the column names

names(US_LA2) <- sub('^US.LA2.', '', names(US_LA2))

# rename "Chamber" to ch_ID

names(US_LA2)[names(US_LA2) == "Chamber"] <- "ch_ID"

# remove unnecessary columns

US_LA2 <- subset(US_LA2, select = -c(Site, Stderr...10, CO2.Flux.rate...ug.CO2.C.m2.hr., CO2.Flux.rate...mg.CO2.C.m2.hr., 
                                     Stderr...15, CO2.C.mg.day...16, CO2.C.mg.day...17, N2O.Flux.rate...ug.N20.N.m2.hr.,
                                     N2O.Flux.rate...microg.N2O.N.m2.hr., Stderr...20, N2O.N.microg.day...21, N2O.N.microg.day...22,
                                     Chamber.H2O.depth..cm., Bdwlk, CH4.C.mg.day...12, B1.soilT, B1.AirT, B2.soilT, B2.AirT) )

US_LA2 <- subset(US_LA2, select = -c(ch_FCH4_HH_mg) ) 
US_LA2 <- subset(US_LA2, select = -c(ch_FCH4_DD_mg) )

# rename columns

names(US_LA2)[names(US_LA2) == "Chamber.Soil.Temp...10cm...C."] <- "TS_10"
names(US_LA2)[names(US_LA2) == "CH4.Flux.rate...ug.CH4.C.m2.hr."] <- "ch_FCH4_HH_ug"
names(US_LA2)[names(US_LA2) == "CH4.Flux.rate...mg.CH4.C.m2.hr."] <- "ch_FCH4_HH_mg"
names(US_LA2)[names(US_LA2) == "CH4.C.mg.day...11"] <- "ch_FCH4_DD_mg"
names(US_LA2)[names(US_LA2) == "Water.Level"] <- "WTL"
names(US_LA2)[names(US_LA2) == "ch_FCH4_HH_ug"] <- "ch_FCH4" 

# convert Date to date format

US_LA2$Date <- ymd(US_LA2$Date)

# add year and time (NA)

US_LA2$Year <- year(US_LA2$Date)

US_LA2$Time <- NA

US_LA2$Hour <- NA

US_LA2$Day <- day(US_LA2$Date)

# move columns to front
US_LA2 <- US_LA2 %>% relocate(Time, Year, Month, Day, Hour, .after = Date)

# add method

US_LA2$method <- "manual"

# save as .csv

write.csv(US_LA2, "US_LA2.csv", row.names = FALSE)

# read in the df again

USLA2_d <- read.csv("path/US_LA2.csv")

# add site

USLA2_d$site <- "US-LA2"

# move site name to the front

USLA2_d <- USLA2_d %>%
  select(site, everything())

# add land cover class

USLA2_d$LC_class <- NA

# add columns for each available TS across the data sets for each site

USLA2_d$TS_2 <- NA
USLA2_d$TS_5 <- NA
USLA2_d$TS_15 <- NA
USLA2_d$TS_30 <- NA

# move TS columns to the right spot

USLA2_d <- USLA2_d %>% relocate(TS_2, TS_5, .before = TS_10)
USLA2_d <- USLA2_d %>% relocate(TS_15, TS_30, .after = TS_10)

# add collar_type

USLA2_d$Collar_type <- NA

# add ch_ prefix to all columns

USLA2_d <- USLA2_d %>% rename_with( ~ paste0("ch_", .x))

# rename SITE, ch_ch_FCH4

colnames(USLA2_d)[8] <- "ch_ID"
colnames(USLA2_d)[14] <- "ch_FCH4_ugCH4m2hr"
colnames(USLA2_d)[1] <- "SITE"

# remove unnecessary columns
USLA2_d <- subset(USLA2_d, select=-c(ch_Time, ch_Year, ch_Month, ch_Day, ch_Hour))

# move the flux column to right spot
USLA2_d <- USLA2_d %>% relocate(ch_FCH4_ugCH4m2hr, .before = ch_TS_2)


## save as .csv

write.csv(USLA2_d, "path", row.names = FALSE)

#### US-LOS ####

US_LOS_df <- read.csv("path")

# remove co2 columns

US_LOS_df = subset(US_LOS_df, select = -c(flux_co2_umol_m2_s) )

# remove chamberht

US_LOS_df = subset(US_LOS_df, select = -c(chamberht) )

# rename some columns

names(US_LOS_df)[names(US_LOS_df) == "ï..Collar_no"] <- "ch_ID"
names(US_LOS_df)[names(US_LOS_df) == "flux_ch4_nmol_m2_s"] <- "ch_FCH4"

# move ch_FCH4 and ch_ID

US_LOS_df <- US_LOS_df %>% relocate(ch_ID, .before = AirTemp_C)
US_LOS_df <- US_LOS_df %>% relocate(ch_FCH4, .after = ch_ID)

# add Date column

US_LOS_df$Date <- as.Date(with(US_LOS_df,paste(Year,Month,Day,sep="-")),"%Y-%m-%d")

# create new column for HHMMSS (Time)
US_LOS_df$Min <- (US_LOS_df$Hour_LT - trunc(US_LOS_df$Hour_LT))*60

US_LOS_df$Sec <- (US_LOS_df$Min - trunc(US_LOS_df$Min))*60

# remove the decimals in Min
US_LOS_df$Min <- sub("\\..*", "", as.character(US_LOS_df$Min))

# remove the decimals in Hour_LT
US_LOS_df$Hour_LT <- sub("\\..*", "", as.character(US_LOS_df$Hour_LT))

# rename Hour_LT to Hour
names(US_LOS_df)[names(US_LOS_df) == "Hour_LT"] <- "Hour"

#convert Sec to character
class(US_LOS_df$Sec)

US_LOS_df$Sec <- as.numeric(US_LOS_df$Sec)

#round the Sec column
US_LOS_df$Sec <- round(US_LOS_df$Sec)

US_LOS_df$Time <- paste(US_LOS_df$Hour,":",US_LOS_df$Min,":",US_LOS_df$Sec,sep="")

# convert Time to hms

US_LOS_df$Time <- hms::as_hms(US_LOS_df$Time)

# move columns to front

US_LOS_df <- US_LOS_df %>% relocate(Date, Time, .before = Year)

#remove min and sec

US_LOS_df = subset(US_LOS_df, select = -c(Min, Sec) )

# add column for method (manual or auto)

US_LOS_df$method <- "manual"

# replace -9999 values with NA
US_LOS_df <- US_LOS_df %>% na_if(-9999.0)

# save as .csv

write.csv(US_LOS_df, "US_LOS.csv", row.names = FALSE)

# read df again

USLOS_d <- read.csv("path/US_LOS.csv")

# add site

USLOS_d$site <- "US-LOS"

# move site name to the front

USLOS_d <- USLOS_d %>%
  select(site, everything())

# change the names of soil moist and temp

colnames(USLOS_d)[which(names(USLOS_d) == "soil_temp")] <- "TS_X" 
colnames(USLOS_d)[which(names(USLOS_d) == "soil_moist")] <- "SM_X" 

# add columns for each available TS across the data sets for each site

USLOS_d$TS_2 <- NA
USLOS_d$TS_5 <- NA
USLOS_d$TS_10 <- NA
USLOS_d$TS_15 <- NA
USLOS_d$TS_30 <- NA

# move TS columns to the right spot

USLOS_d <- USLOS_d %>% relocate(TS_2, TS_5, TS_10, TS_15, TS_30, .after = TS_X)

# add WTL column

USLOS_d$WTL <- NA
USLOS_d <- USLOS_d %>% relocate(WTL, .after = SM_X)

# add land cover class

USLOS_d <- USLOS_d %>% mutate(LC_class =
                                case_when(ch_ID == 2 ~ "Willow", 
                                          ch_ID == 3 ~ "Alder",
                                          ch_ID == 4 ~ "Alder",
                                          ch_ID == 5 ~ "Alder",
                                          ch_ID == 6 ~ "Alder", 
                                          ch_ID == 7 ~ "Willow",
                                          ch_ID == 8 ~ "Willow",
                                          ch_ID == 9 ~ "Mixed",
                                          ch_ID == 10 ~ "Mixed",
                                          ch_ID == 11 ~ "Mixed",
                                          ch_ID == 12 ~ "Mixed",
                                          ch_ID == 13 ~ "Mixed",
                                          ch_ID == 14 ~ "Mixed",
                                          ch_ID == 15 ~ "Mixed")
)

# add ch_ prefix to all columns

USLOS_d <- USLOS_d %>% rename_with( ~ paste0("ch_", .x))

# rename SITE, ch_ch_FCH4

colnames(USLOS_d)[8] <- "ch_ID"
colnames(USLOS_d)[9] <- "ch_FCH4_nmolm2s1"
colnames(USLOS_d)[1] <- "SITE"

# combine date and time to datetime

# specifying the format
format <- "%Y-%m-%d %H:%M:%S"

USLOS_d$ch_Date <- as.character(USLOS_d$ch_Date)
USLOS_d$ch_Time <- as.character(USLOS_d$ch_Time)

# combining date and time into single object
USLOS_d$ch_Datetime <- as_datetime(paste(USLOS_d$ch_Date, USLOS_d$ch_Time), format=format)

# move datetime to front

USLOS_d<- USLOS_d %>% relocate(ch_Datetime, .after = SITE)

# remove unnecessary columns
USLOS_d <- subset(USLOS_d, select=-c(ch_Date, ch_Time, ch_Year, ch_Month, ch_Day, ch_Hour))

# arrange the df according to datetime
USLOS_d <- arrange(USLOS_d, ch_Datetime)

## save as .csv

write.csv(USLOS_d, "path", row.names = FALSE)

#### US-OWC ####

# US-OWC

US_OWC_tibble <- read.csv("path")

# convert to dataframe

US_OWC <- as.data.frame(US_OWC_tibble)

# replace -9999 values with NA
US_OWC <- US_OWC %>% na_if(-9999)

# add ch_ID column indicating the "chamber ID" which didn't exist here as they were not necessarily fixed.
# 1=open water, 2=Typha, 3=Nelumbo, 4=mud flat

US_OWC <- US_OWC %>% mutate(ch_ID =
                                    case_when(CMB_SPP == "Mud Flat" ~ 4, 
                                              CMB_SPP == "Open Water" ~ 1,
                                              CMB_SPP == "Typha" ~ 2,
                                              CMB_SPP == "Nelumbo" ~ 3)
)

# rename LOC_latitude column

colnames(US_OWC)[1] <- "LOC_LATITUDE"

# convert Datetime to datetime

US_OWC$Datetime <- ymd_hm(US_OWC$Datetime)

# new column for year, date, month, hour, min, sec

US_OWC$Year <- year(US_OWC$Datetime)

US_OWC$Month <- month(US_OWC$Datetime)

US_OWC$Date <- date(US_OWC$Datetime)

US_OWC$Hour <- hour(US_OWC$Datetime)

US_OWC$Min <- minute(US_OWC$Datetime)

US_OWC$Sec <- second(US_OWC$Datetime)

# combine hour min sec into Time column

US_OWC$Time <- paste(US_OWC$Hour,US_OWC$Min,"00",sep=":")

# convert to hms

US_OWC$Time <- hms::as_hms(US_OWC$Time)

# remove min and sec

US_OWC = subset(US_OWC, select = -c(Min, Sec) )

# rename FCH4
colnames(US_OWC)[8] <- "ch_FCH4"

# remove CMB_TA

US_OWC = subset(US_OWC, select = -c(CMB_TA, CMB_AREA_CM2) )
US_OWC = subset(US_OWC, select = -c(CMB_VEGTYPE) )
US_OWC = subset(US_OWC, select = -c(CMB_SPP) )

# move columns to front
US_OWC <- US_OWC %>% relocate(Datetime, Date, Time, Year, Month, Hour, .before = LOC_LATITUDE)

US_OWC <- US_OWC %>% relocate(ch_ID, .after = LOC_LONGITUDE)

# save as .csv

write.csv(US_OWC, "US_OWC.csv", row.names=FALSE)

# read the df again

USOWC_d <- read.csv("path/US_OWC.csv")

# add site

USOWC_d$site <- "US-OWC"

# move site name to the front

USOWC_d <- USOWC_d %>%
  select(SITE, everything())

# add columns for each available TS across the data sets for each site

USOWC_d$TS_2 <- NA
USOWC_d$TS_5 <- NA
USOWC_d$TS_10 <- NA
USOWC_d$TS_15 <- NA
USOWC_d$TS_30 <- NA

# move TS columns to the right spot

USOWC_d <- USOWC_d %>% relocate(TS_2, TS_5, TS_10, TS_15, TS_30, .after = ch_FCH4)

# add WTL column

USOWC_d$WTL <- NA
USOWC_d <- USOWC_d %>% relocate(WTL, .after = TS_30)

# add land cover class

USOWC_d <- USOWC_d %>% mutate(LC_class =
                                case_when(ch_ID == 4 ~ "Mud Flat", 
                                          ch_ID == 1 ~ "Open Water",
                                          ch_ID == 2 ~ "Typha",
                                          ch_ID == 3 ~ "Nelumbo")
)

# add ch_ prefix to all columns

USOWC_d <- USOWC_d %>% rename_with( ~ paste0("ch_", .x))

# rename SITE, ch_ch_FCH4

colnames(USOWC_d)[9] <- "ch_ID"
colnames(USOWC_d)[10] <- "ch_FCH4_nmolCH4m2s1"
USOWC_d$SITE <- "US-OWC"

# remove some columns that are  not needed
USOWC_d <- subset(USOWC_d, select=-c(ch_Date, ch_Time, ch_Year, ch_Month, ch_Hour))

# arrange according to ch_Datetime
USOWC_d <- arrange(USOWC_d, ch_Datetime)

## save as .csv

write.csv(USOWC_d, "path", row.names = FALSE)

######### US-StJ ###########

USSTJ_df <- read.csv("path")

# remove extra columns 

USSTJ_df <- subset(USSTJ_df, select = -c(Tide, CO2_flux_dark))

colnames(USSTJ_df)[2] <- "ch_ID"

# change Date from US to international

USSTJ_df$Date <- as.Date(USSTJ_df$Date, format = "%m/%d/%Y")

# combine Date and Time

USSTJ_df$datetime <- as.POSIXct(
  paste(USSTJ_df$Date, USSTJ_df$Time),
  format = "%Y-%m-%d %H:%M:%S"
)

USSTJ_df <- USSTJ_df %>% relocate(datetime, .before = Date)
USSTJ_df <- USSTJ_df %>% relocate(Time, .after = Date)

USSTJ_df$SITE <- "US-STJ"
USSTJ_df <- USSTJ_df %>% relocate(SITE, .before = datetime)

USSTJ_df$ch_method <- "manual"

colnames(USSTJ_df)[5] <- "ch_ID"

colnames(USSTJ_df)[6] <- "ch_FCH4_nmolCH4m2s1"

write.csv(USSTJ_df, "path", row.names = FALSE)


#### US-UAF ####
# (separate dfs for each year)

USUAF_2016 <- read.csv("path")

# replace #N/A with NA
USUAF_2016 <- USUAF_2016 %>% na_if("#N/A")

# remove Fco2 columns
USUAF_2016 = subset(USUAF_2016, select = -c(Fco2_1_ws1, Fco2_2_dl1, Fco2_3_wc1) )

#create new column for Time2 with seconds
USUAF_2016$TimeStart2 <- paste(as.character(USUAF_2016$StartTime),":00",sep="")

# convert to datetime

USUAF_2016$TimeStart2 <- dmy_hms(as.character(USUAF_2016$TimeStart2))

# remove Timestamp

USUAF_2016 = subset(USUAF_2016, select = -c(StartTime) )

# move to the front

USUAF_2016 <- USUAF_2016 %>% relocate(TimeStart2, .before = Fch4_1_ws1)

# rename
colnames(USUAF_2016)[1] <- "ch_Datetime"

###### choose by chamber

# 1

USUAF_2016_flux_1ws1 <- select(USUAF_2016, "ch_Datetime", "Fch4_1_ws1")

# add ch_ID column

USUAF_2016_flux_1ws1$ch_ID <- 1

# add veg column
USUAF_2016_flux_1ws1$veg_cover <- "wet_sphagnum"

# rename veg column
colnames(USUAF_2016_flux_1ws1)[4] <- "ch_veg_cover"

# rename original chamber ID column
colnames(USUAF_2016_flux_1ws1)[2] <- "FCH4"

# 2

USUAF_2016_flux_2dl1 <- select(USUAF_2016, "ch_Datetime", "Fch4_2_dl1")

# add ch_ID column

USUAF_2016_flux_2dl1$ch_ID <- 2

# add veg column
USUAF_2016_flux_2dl1$veg_cover <- "dry_lichen"

# rename veg column
colnames(USUAF_2016_flux_2dl1)[4] <- "ch_veg_cover"

# rename original chamber ID column
colnames(USUAF_2016_flux_2dl1)[2] <- "FCH4"

# 3

USUAF_2016_flux_3wc1 <- select(USUAF_2016, "ch_Datetime", "Fch4_3_wc1")

# add ch_ID column

USUAF_2016_flux_3wc1$ch_ID <- 3

# add veg column
USUAF_2016_flux_3wc1$ch_veg_cover <- "wet_carex"

# rename original chamber ID column
colnames(USUAF_2016_flux_3wc1)[2] <- "FCH4"

# combine the flux dfs

USUAF_2016_flux_all <- bind_rows(USUAF_2016_flux_1ws1, USUAF_2016_flux_2dl1, USUAF_2016_flux_3wc1)

## VWC (volumetric water content; soil moisture)

# 1

USUAF_2016_vwc_1ws1 <- select(USUAF_2016, "ch_Datetime", "VWC1_ws1")

# add ch_ID column

USUAF_2016_vwc_1ws1$ch_ID <- 1

# add veg column
USUAF_2016_vwc_1ws1$ch_veg_cover <- "wet_sphagnum"

# rename original VWC column
colnames(USUAF_2016_vwc_1ws1)[2] <- "ch_VWC"

# 2

USUAF_2016_vwc_2dl1 <- select(USUAF_2016, "ch_Datetime", "VWC2_dl1")

# add ch_ID column

USUAF_2016_vwc_2dl1$ch_ID <- 2

# add veg column
USUAF_2016_vwc_2dl1$ch_veg_cover <- "dry_lichen"

# rename original VWC column
colnames(USUAF_2016_vwc_2dl1)[2] <- "ch_VWC"

# 3

USUAF_2016_vwc_3wc1 <- select(USUAF_2016, "ch_Datetime", "VWC3_wc1")

# add ch_ID column

USUAF_2016_vwc_3wc1$ch_ID <- 3

# add veg column
USUAF_2016_vwc_3wc1$ch_veg_cover <- "wet_carex"

# rename original VWC column
colnames(USUAF_2016_vwc_3wc1)[2] <- "ch_VWC"

# combine the dfs into one

USUAF_2016_VWC_all <- bind_rows(USUAF_2016_vwc_1ws1, USUAF_2016_vwc_2dl1, USUAF_2016_vwc_3wc1)

## soil temperature

# first coalesce Tsoil_10cm_ws1 and Tsoil_10cm_ws1.1
USUAF_2016$Tsoil_10cm_ws1_vol2 <- coalesce(USUAF_2016$Tsoil_10cm_ws1, USUAF_2016$Tsoil_10cm_ws1.1)

# remove the old Tsoil 10 cm columns
USUAF_2016 = subset(USUAF_2016, select = -c(Tsoil_10cm_ws1, Tsoil_10cm_ws1.1) )

# rename 
colnames(USUAF_2016)[16] <- "Tsoil_10cm_ws1"

# 1

USUAF_2016_tsoil_1ws1_10 <- select(USUAF_2016, "ch_Datetime", "Tsoil_10cm_ws1")

# add ch_ID column

USUAF_2016_tsoil_1ws1_10$ch_ID <- 1

# add veg column
USUAF_2016_tsoil_1ws1_10$ch_veg_cover <- "wet_sphagnum"

# rename original VWC column
colnames(USUAF_2016_tsoil_1ws1_10)[2] <- "ch_TS_10"

# 1, 2 cm

USUAF_2016_tsoil_1ws1_2 <- select(USUAF_2016, "ch_Datetime", "Tsoil_2cm_ws1")

# add ch_ID column

USUAF_2016_tsoil_1ws1_2$ch_ID <- 1

# add veg column
USUAF_2016_tsoil_1ws1_2$ch_veg_cover <- "wet_sphagnum"

# rename original VWC column
colnames(USUAF_2016_tsoil_1ws1_2)[2] <- "ch_TS_2"

# 2, 2 cm

USUAF_2016_tsoil_2dl1_2 <- select(USUAF_2016, "ch_Datetime", "Tsoil2_2cm_dl1")

# add ch_ID column

USUAF_2016_tsoil_2dl1_2$ch_ID <- 2

# add veg column
USUAF_2016_tsoil_2dl1_2$ch_veg_cover <- "dry_lichen"

# rename original VWC column
colnames(USUAF_2016_tsoil_2dl1_2)[2] <- "ch_TS_2"

# 3, 2 cm

USUAF_2016_tsoil_3wc1_2 <- select(USUAF_2016, "ch_Datetime", "Tsoil3_2cm_wc1")

# add ch_ID column

USUAF_2016_tsoil_3wc1_2$ch_ID <- 3

# add veg column
USUAF_2016_tsoil_3wc1_2$ch_veg_cover <- "wet_carex"

# rename original tsoil column
colnames(USUAF_2016_tsoil_3wc1_2)[2] <- "ch_TS_2"

# 1, 20 cm

USUAF_2016_tsoil_1ws1_20 <- select(USUAF_2016, "ch_Datetime", "Tsoil1_20cm_ws1")

# add ch_ID column

USUAF_2016_tsoil_1ws1_20$ch_ID <- 1

# add veg column
USUAF_2016_tsoil_1ws1_20$ch_veg_cover <- "wet_sphagnum"

# rename original VWC column
colnames(USUAF_2016_tsoil_1ws1_20)[2] <- "ch_TS_20"

# 1, 30 cm

USUAF_2016_tsoil_1ws1_30 <- select(USUAF_2016, "ch_Datetime", "Tsoil1_30cm_ws1")

# add ch_ID column

USUAF_2016_tsoil_1ws1_30$ch_ID <- 1

# add veg column
USUAF_2016_tsoil_1ws1_30$ch_veg_cover <- "wet_sphagnum"

# rename original VWC column
colnames(USUAF_2016_tsoil_1ws1_30)[2] <- "ch_TS_30"

# 1, 40 cm

USUAF_2016_tsoil_1ws1_40 <- select(USUAF_2016, "ch_Datetime", "Tsoil1_40cm_ws1")

# add ch_ID column

USUAF_2016_tsoil_1ws1_40$ch_ID <- 1

# add veg column
USUAF_2016_tsoil_1ws1_40$ch_veg_cover <- "wet_sphagnum"

# rename original VWC column
colnames(USUAF_2016_tsoil_1ws1_40)[2] <- "ch_TS_40"

# combine the chamber-specific measurements into one
USUAF_tsoil_1ws1_all_1 <- merge(x = USUAF_2016_tsoil_1ws1_2, y = USUAF_2016_tsoil_1ws1_10, by = c("ch_Datetime", "ch_ID", "ch_veg_cover"))

USUAF_tsoil_1ws1_all_2 <- merge(x = USUAF_tsoil_1ws1_all_1, y = USUAF_2016_tsoil_1ws1_20, by = c("ch_Datetime", "ch_ID", "ch_veg_cover"))

USUAF_tsoil_1ws1_all_3 <- merge(x = USUAF_tsoil_1ws1_all_2, y = USUAF_2016_tsoil_1ws1_30, by = c("ch_Datetime", "ch_ID", "ch_veg_cover"))

USUAF_tsoil_1ws1_all_4 <- merge(x = USUAF_tsoil_1ws1_all_3, y = USUAF_2016_tsoil_1ws1_40, by = c("ch_Datetime", "ch_ID", "ch_veg_cover"))

# merge tsoil with vwc
USUAF_2016_tsoil_vwc <- merge(x = USUAF_2016_tsoil_all, y = USUAF_2016_VWC_all, by = c("ch_Datetime", "ch_ID", "ch_veg_cover"))

# combine all dfs into one

USUAF_2016_all <- merge(x = USUAF_2016_flux_all, y = USUAF_2016_tsoil_vwc, by = c("ch_Datetime", "ch_ID", "ch_veg_cover"))

### 2017 ###

# import 2017
USUAF_2017 <- read.csv("path")

# replace #N/A with NA
USUAF_2017 <- USUAF_2017 %>% na_if("#N/A")

# remove Fco2 columns
USUAF_2017 = subset(USUAF_2017, select = -c(Fco2_1_ws1, Fco2_4_wc2, Fco2_3_wc1) )

#create new column for Time2 with seconds
USUAF_2017$TimeStart2 <- paste(as.character(USUAF_2017$StartTime),":00",sep="")

# convert to datetime

USUAF_2017$TimeStart2 <- dmy_hms(as.character(USUAF_2017$TimeStart2))

# remove Timestamp

USUAF_2017 = subset(USUAF_2017, select = -c(StartTime) )

# move to the front

USUAF_2017 <- USUAF_2017 %>% relocate(TimeStart2, .before = Fch4_1_ws1)

# rename
colnames(USUAF_2017)[1] <- "ch_Datetime"

# fluxes
###### choose by chamber

# 1

USUAF_2017_flux_1ws1 <- select(USUAF_2017, "ch_Datetime", "Fch4_1_ws1")

# add ch_ID column

USUAF_2017_flux_1ws1$ch_ID <- 1

# add veg column
USUAF_2017_flux_1ws1$veg_cover <- "wet_sphagnum"

# rename veg column
colnames(USUAF_2017_flux_1ws1)[4] <- "ch_veg_cover"

# rename original chamber ID column
colnames(USUAF_2017_flux_1ws1)[2] <- "FCH4"

# 4

USUAF_2017_flux_4wc2 <- select(USUAF_2017, "ch_Datetime", "Fch4_4_wc2")

# add ch_ID column

USUAF_2017_flux_4wc2$ch_ID <- 4

# add veg column
USUAF_2017_flux_4wc2$veg_cover <- "wet_carex"

# rename veg column
colnames(USUAF_2017_flux_4wc2)[4] <- "ch_veg_cover"

# rename original chamber ID column
colnames(USUAF_2017_flux_4wc2)[2] <- "FCH4"

# 3

USUAF_2017_flux_3wc1 <- select(USUAF_2017, "ch_Datetime", "Fch4_3_wc1")

# add ch_ID column

USUAF_2017_flux_3wc1$ch_ID <- 3

# add veg column
USUAF_2017_flux_3wc1$ch_veg_cover <- "wet_carex"

# rename original chamber ID column
colnames(USUAF_2017_flux_3wc1)[2] <- "FCH4"

# combine the flux dfs

USUAF_2017_flux_all <- bind_rows(USUAF_2017_flux_1ws1, USUAF_2017_flux_4wc2, USUAF_2017_flux_3wc1)

# add to the comb flux df
USUAF_2017_flux_all_2 <- bind_rows(USUAF_2017_flux_all, USUAF_2017_flux_2dl1)

## VWC

# 1

USUAF_2017_vwc_1ws1 <- select(USUAF_2017, "ch_Datetime", "VWC1_ws1")

# add ch_ID column

USUAF_2017_vwc_1ws1$ch_ID <- 1

# add veg column
USUAF_2017_vwc_1ws1$ch_veg_cover <- "wet_sphagnum"

# rename original VWC column
colnames(USUAF_2017_vwc_1ws1)[2] <- "ch_VWC"

# 2

USUAF_2017_vwc_2dl1 <- select(USUAF_2017, "ch_Datetime", "VWC2_dl1")

# add ch_ID column

USUAF_2017_vwc_2dl1$ch_ID <- 2

# add veg column
USUAF_2017_vwc_2dl1$ch_veg_cover <- "dry_lichen"

# rename original VWC column
colnames(USUAF_2017_vwc_2dl1)[2] <- "ch_VWC"

# 3

USUAF_2017_vwc_3wc1 <- select(USUAF_2017, "ch_Datetime", "VWC3_wc1")

# add ch_ID column

USUAF_2017_vwc_3wc1$ch_ID <- 3

# add veg column
USUAF_2017_vwc_3wc1$ch_veg_cover <- "wet_carex"

# rename original VWC column
colnames(USUAF_2017_vwc_3wc1)[2] <- "ch_VWC"

# combine the dfs into one

USUAF_2017_VWC_all <- bind_rows(USUAF_2017_vwc_1ws1, USUAF_2017_vwc_2dl1, USUAF_2017_vwc_3wc1)

## Tsoil

# 1, 2 cm

USUAF_2017_tsoil_1ws1_2 <- select(USUAF_2017, "ch_Datetime", "Tsoil_2cm_ws1")

# add ch_ID column

USUAF_2017_tsoil_1ws1_2$ch_ID <- 1

# add veg column
USUAF_2017_tsoil_1ws1_2$ch_veg_cover <- "wet_sphagnum"

# rename original VWC column
colnames(USUAF_2017_tsoil_1ws1_2)[2] <- "ch_TS_2"

# 2, 2 cm

USUAF_2017_tsoil_2dl1_2 <- select(USUAF_2017, "ch_Datetime", "Tsoil2_2cm_dl1")

# add ch_ID column

USUAF_2017_tsoil_2dl1_2$ch_ID <- 2

# add veg column
USUAF_2017_tsoil_2dl1_2$ch_veg_cover <- "dry_lichen"

# rename original VWC column
colnames(USUAF_2017_tsoil_2dl1_2)[2] <- "ch_TS_2"

# 3, 2 cm

USUAF_2017_tsoil_3wc1_2 <- select(USUAF_2017, "ch_Datetime", "Tsoil3_2cm_wc1")

# add ch_ID column

USUAF_2017_tsoil_3wc1_2$ch_ID <- 3

# add veg column
USUAF_2017_tsoil_3wc1_2$ch_veg_cover <- "wet_carex"

# rename original VWC column
colnames(USUAF_2017_tsoil_3wc1_2)[2] <- "ch_TS_2"

# 1, 10 cm

USUAF_2017_tsoil_1ws1_10 <- select(USUAF_2017, "ch_Datetime", "Tsoil_10cm_ws1")

# add ch_ID column

USUAF_2017_tsoil_1ws1_10$ch_ID <- 1

# add veg column
USUAF_2017_tsoil_1ws1_10$ch_veg_cover <- "wet_sphagnum"

# rename original tsoil column
colnames(USUAF_2017_tsoil_1ws1_10)[2] <- "ch_TS_10"

# 1, 20 cm

USUAF_2017_tsoil_1ws1_20 <- select(USUAF_2017, "ch_Datetime", "Tsoil1_20cm_ws1")

# add ch_ID column

USUAF_2017_tsoil_1ws1_20$ch_ID <- 1

# add veg column
USUAF_2017_tsoil_1ws1_20$ch_veg_cover <- "wet_sphagnum"

# rename original VWC column
colnames(USUAF_2017_tsoil_1ws1_20)[2] <- "ch_TS_20"

# 1, 30 cm

USUAF_2017_tsoil_1ws1_30 <- select(USUAF_2017, "ch_Datetime", "Tsoil1_30cm_ws1")

# add ch_ID column

USUAF_2017_tsoil_1ws1_30$ch_ID <- 1

# add veg column
USUAF_2017_tsoil_1ws1_30$ch_veg_cover <- "wet_sphagnum"

# rename original VWC column
colnames(USUAF_2017_tsoil_1ws1_30)[2] <- "ch_TS_30"

# 1, 40 cm

USUAF_2017_tsoil_1ws1_40 <- select(USUAF_2017, "ch_Datetime", "Tsoil1_40cm_ws1")

# add ch_ID column

USUAF_2017_tsoil_1ws1_40$ch_ID <- 1

# add veg column
USUAF_2017_tsoil_1ws1_40$ch_veg_cover <- "wet_sphagnum"

# rename original VWC column
colnames(USUAF_2017_tsoil_1ws1_40)[2] <- "ch_TS_40"

# combine the chamber-specific measurements into one
USUAF_2017_tsoil_1ws1_all_1 <- merge(x = USUAF_2017_tsoil_1ws1_2, y = USUAF_2017_tsoil_1ws1_10, by = c("ch_Datetime", "ch_ID", "ch_veg_cover"))

USUAF_2017_tsoil_1ws1_all_2 <- merge(x = USUAF_2017_tsoil_1ws1_all_1, y = USUAF_2017_tsoil_1ws1_20, by = c("ch_Datetime", "ch_ID", "ch_veg_cover"))

USUAF_2017_tsoil_1ws1_all_3 <- merge(x = USUAF_2017_tsoil_1ws1_all_2, y = USUAF_2017_tsoil_1ws1_30, by = c("ch_Datetime", "ch_ID", "ch_veg_cover"))

USUAF_2017_tsoil_1ws1_all_4 <- merge(x = USUAF_2017_tsoil_1ws1_all_3, y = USUAF_2017_tsoil_1ws1_40, by = c("ch_Datetime", "ch_ID", "ch_veg_cover"))

# combine all temps
USUAF_2017_tsoil_all <- bind_rows(USUAF_2017_tsoil_1ws1_all_4, USUAF_2017_tsoil_2dl1_2, USUAF_2017_tsoil_3wc1_2)

# merge tsoil with vwc
USUAF_2017_tsoil_vwc <- merge(x = USUAF_2017_tsoil_all, y = USUAF_2017_VWC_all, by = c("ch_Datetime", "ch_ID", "ch_veg_cover"))

# combine all dfs into one

USUAF_2017_all <- merge(x = USUAF_2017_tsoil_vwc, y = USUAF_2017_flux_all, by = c("ch_Datetime", "ch_ID", "ch_veg_cover"), all=TRUE)

# 2018

USUAF_2018 <- read.csv("path")

# replace #N/A with NA
USUAF_2018 <- USUAF_2018 %>% na_if("#N/A")

# remove Fco2 columns
USUAF_2018 = subset(USUAF_2018, select = -c(Fco2_1_ws1, Fco2_4_wc2, Fco2_3_wc1, Fco2_5_dc1) )

#create new column for Time2 with seconds
USUAF_2018$TimeStart2 <- paste(as.character(USUAF_2018$StartTime),":00",sep="")

# convert to datetime

USUAF_2018$TimeStart2 <- dmy_hms(as.character(USUAF_2018$TimeStart2))

# remove Timestamp

USUAF_2018 = subset(USUAF_2018, select = -c(StartTime) )

# move to the front

USUAF_2018 <- USUAF_2018 %>% relocate(TimeStart2, .before = Fch4_1_ws1)

# rename
colnames(USUAF_2018)[1] <- "ch_Datetime"

# fluxes
###### choose by chamber

# 1

USUAF_2018_flux_1ws1 <- select(USUAF_2018, "ch_Datetime", "Fch4_1_ws1")

# add ch_ID column

USUAF_2018_flux_1ws1$ch_ID <- 1

# add veg column
USUAF_2018_flux_1ws1$veg_cover <- "wet_sphagnum"

# rename veg column
colnames(USUAF_2018_flux_1ws1)[4] <- "ch_veg_cover"

# rename original chamber ID column
colnames(USUAF_2018_flux_1ws1)[2] <- "FCH4"

# 3

USUAF_2018_flux_3wc1 <- select(USUAF_2018, "ch_Datetime", "Fch4_3_wc1")

# add ch_ID column

USUAF_2018_flux_3wc1$ch_ID <- 3

# add veg column
USUAF_2018_flux_3wc1$veg_cover <- "wet_carex"

# rename veg column
colnames(USUAF_2018_flux_3wc1)[4] <- "ch_veg_cover"

# rename original chamber ID column
colnames(USUAF_2018_flux_3wc1)[2] <- "FCH4"

# 3

USUAF_2018_flux_4wc2 <- select(USUAF_2018, "ch_Datetime", "Fch4_4_wc2")

# add ch_ID column

USUAF_2018_flux_4wc2$ch_ID <- 4

# add veg column
USUAF_2018_flux_4wc2$ch_veg_cover <- "wet_carex"

# rename original chamber ID column
colnames(USUAF_2018_flux_4wc2)[2] <- "FCH4"

# 5

USUAF_2018_flux_5dc1 <- select(USUAF_2018, "ch_Datetime", "Fch4_5_dc1")

# add ch_ID column

USUAF_2018_flux_5dc1$ch_ID <- 5

# add veg column
USUAF_2018_flux_5dc1$ch_veg_cover <- "dry_carex"

# rename original chamber ID column
colnames(USUAF_2018_flux_5dc1)[2] <- "FCH4"

# combine the flux dfs

USUAF_2018_flux_all <- bind_rows(USUAF_2018_flux_1ws1, USUAF_2018_flux_4wc2, USUAF_2018_flux_3wc1, USUAF_2018_flux_5dc1)

## VWC

# 1

USUAF_2018_vwc_1ws1 <- select(USUAF_2018, "ch_Datetime", "VWC1_ws1")

# add ch_ID column

USUAF_2018_vwc_1ws1$ch_ID <- 1

# add veg column
USUAF_2018_vwc_1ws1$ch_veg_cover <- "wet_sphagnum"

# rename original VWC column
colnames(USUAF_2018_vwc_1ws1)[2] <- "ch_VWC"

# 2

USUAF_2018_vwc_2dl1 <- select(USUAF_2018, "ch_Datetime", "VWC2_dl1")

# add ch_ID column

USUAF_2018_vwc_2dl1$ch_ID <- 2

# add veg column
USUAF_2018_vwc_2dl1$ch_veg_cover <- "dry_lichen"

# rename original VWC column
colnames(USUAF_2018_vwc_2dl1)[2] <- "ch_VWC"

# 3

USUAF_2018_vwc_3wc1 <- select(USUAF_2018, "ch_Datetime", "VWC3_wc1")

# add ch_ID column

USUAF_2018_vwc_3wc1$ch_ID <- 3

# add veg column
USUAF_2018_vwc_3wc1$ch_veg_cover <- "wet_carex"

# rename original VWC column
colnames(USUAF_2018_vwc_3wc1)[2] <- "ch_VWC"

# combine the dfs into one

USUAF_2018_VWC_all <- bind_rows(USUAF_2018_vwc_1ws1, USUAF_2018_vwc_2dl1, USUAF_2018_vwc_3wc1)

## Tsoil

# 1, 2 cm

USUAF_2018_tsoil_1ws1_2 <- select(USUAF_2018, "ch_Datetime", "Tsoil_2cm_ws1")

# add ch_ID column

USUAF_2018_tsoil_1ws1_2$ch_ID <- 1

# add veg column
USUAF_2018_tsoil_1ws1_2$ch_veg_cover <- "wet_sphagnum"

# rename original VWC column
colnames(USUAF_2018_tsoil_1ws1_2)[2] <- "ch_TS_2"

# 2, 2 cm

USUAF_2018_tsoil_2dl1_2 <- select(USUAF_2018, "ch_Datetime", "Tsoil2_2cm_dl1")

# add ch_ID column

USUAF_2018_tsoil_2dl1_2$ch_ID <- 2

# add veg column
USUAF_2018_tsoil_2dl1_2$ch_veg_cover <- "dry_lichen"

# rename original VWC column
colnames(USUAF_2018_tsoil_2dl1_2)[2] <- "ch_TS_2"

# 3, 2 cm

USUAF_2018_tsoil_3wc1_2 <- select(USUAF_2018, "ch_Datetime", "Tsoil3_2cm_wc1")

# add ch_ID column

USUAF_2018_tsoil_3wc1_2$ch_ID <- 3

# add veg column
USUAF_2018_tsoil_3wc1_2$ch_veg_cover <- "wet_carex"

# rename original VWC column
colnames(USUAF_2018_tsoil_3wc1_2)[2] <- "ch_TS_2"

# 1, 10 cm

USUAF_2018_tsoil_1ws1_10 <- select(USUAF_2018, "ch_Datetime", "Tsoil_10cm_ws1")

# add ch_ID column

USUAF_2018_tsoil_1ws1_10$ch_ID <- 1

# add veg column
USUAF_2018_tsoil_1ws1_10$ch_veg_cover <- "wet_sphagnum"

# rename original tsoil column
colnames(USUAF_2018_tsoil_1ws1_10)[2] <- "ch_TS_10"

# 1, 20 cm

USUAF_2018_tsoil_1ws1_20 <- select(USUAF_2018, "ch_Datetime", "Tsoil1_20cm_ws1")

# add ch_ID column

USUAF_2018_tsoil_1ws1_20$ch_ID <- 1

# add veg column
USUAF_2018_tsoil_1ws1_20$ch_veg_cover <- "wet_sphagnum"

# rename original VWC column
colnames(USUAF_2018_tsoil_1ws1_20)[2] <- "ch_TS_20"

# 1, 30 cm

USUAF_2018_tsoil_1ws1_30 <- select(USUAF_2018, "ch_Datetime", "Tsoil1_30cm_ws1")

# add ch_ID column

USUAF_2018_tsoil_1ws1_30$ch_ID <- 1

# add veg column
USUAF_2018_tsoil_1ws1_30$ch_veg_cover <- "wet_sphagnum"

# rename original VWC column
colnames(USUAF_2018_tsoil_1ws1_30)[2] <- "ch_TS_30"

# 1, 40 cm

USUAF_2018_tsoil_1ws1_40 <- select(USUAF_2018, "ch_Datetime", "Tsoil1_40cm_ws1")

# add ch_ID column

USUAF_2018_tsoil_1ws1_40$ch_ID <- 1

# add veg column
USUAF_2018_tsoil_1ws1_40$ch_veg_cover <- "wet_sphagnum"

# rename original VWC column
colnames(USUAF_2018_tsoil_1ws1_40)[2] <- "ch_TS_40"

# combine the chamber-specific measurements into one
USUAF_2018_tsoil_1ws1_all_1 <- merge(x = USUAF_2018_tsoil_1ws1_2, y = USUAF_2018_tsoil_1ws1_10, by = c("ch_Datetime", "ch_ID", "ch_veg_cover"))

USUAF_2018_tsoil_1ws1_all_2 <- merge(x = USUAF_2018_tsoil_1ws1_all_1, y = USUAF_2018_tsoil_1ws1_20, by = c("ch_Datetime", "ch_ID", "ch_veg_cover"))

USUAF_2018_tsoil_1ws1_all_3 <- merge(x = USUAF_2018_tsoil_1ws1_all_2, y = USUAF_2018_tsoil_1ws1_30, by = c("ch_Datetime", "ch_ID", "ch_veg_cover"))

USUAF_2018_tsoil_1ws1_all_4 <- merge(x = USUAF_2018_tsoil_1ws1_all_3, y = USUAF_2018_tsoil_1ws1_40, by = c("ch_Datetime", "ch_ID", "ch_veg_cover"))

# combine all temps
USUAF_2018_tsoil_all <- bind_rows(USUAF_2018_tsoil_1ws1_all_4, USUAF_2018_tsoil_2dl1_2, USUAF_2018_tsoil_3wc1_2)

# merge tsoil with vwc
USUAF_2018_tsoil_vwc <- merge(x = USUAF_2018_tsoil_all, y = USUAF_2018_VWC_all, by = c("ch_Datetime", "ch_ID", "ch_veg_cover"))

# combine all dfs into one

USUAF_2018_all <- merge(x = USUAF_2018_tsoil_vwc, y = USUAF_2018_flux_all, by = c("ch_Datetime", "ch_ID", "ch_veg_cover"), all=TRUE)

# combine the year-specific dfs into one df

USUAF_all <- bind_rows(USUAF_2016_all, USUAF_2017_all, USUAF_2018_all)

write.csv(USUAF_all, "path/US_UAF_chamber.csv", row.names = FALSE)

# read csv

USUAF_d <- read.csv("path/US_UAF_chamber.csv")

# add site

USUAF_d$SITE <- "US-UAF"

# move site name to the front

USUAF_d <- USUAF_d %>%
  select(SITE, everything())

# change column names
colnames(USUAF_d)[5] <- "ch_FCH4_nmolCH4m2s1"
colnames(USUAF_d)[12] <- "ch_Collar_type"

# subset to start from 02.05.2016 for full days and end 2018-11-01
Sys.setenv(TZ = "UTC")
USUAF_d$ch_Datetime <- as_datetime(USUAF_d$ch_Datetime)
USUAF_d_sub <- subset(USUAF_d, date(ch_Datetime) >= "2016-05-02" & date(ch_Datetime) <= "2018-11-02")

# add method
USUAF_d_sub$method <- "auto"

## save as .csv

write.csv(USUAF_d_sub, "path", row.names = FALSE)
