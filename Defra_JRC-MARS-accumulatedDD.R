# R code to automatically calculate degree days for the JRC-MARS gridded climate data

# ================================
# Anastasia Korycinska, Defra Risk and Horizon Scanning Team
# Animal and Plant Health Directorate, Defra, UK
# ================================

# SET THE THRESHOLD TEMPERATURE FOR DEVELOPMENT (oC) IN "threshold" [code line 30]
# SET THE ACCUMULATED DAY DEGREE THRESHOLD IN "accumulated_threshold" [code line 31]
# Set the input file folder by changing the input file path [code line 35]
# Set the output file folder  by changing the output file path [code line 205]
#     just remember to use / as a separator and enclose the path in ""

# Required fields: GRID_NO is a location identifier
#                  DAY is the date in YYYYMMDD format (no separators)
#                  TEMPERATURE_MAX and TEMPERATURE_MIN are self-explanatory

# Accumulated day degrees by grid square for each year in the data set are generated as output
# Min, mean and maximum accumulated day degrees for each grid square across all input years also created
# The number and percentage of years out of the total years analysed which have an actual accumulated threshold greater than
# or equal to the accumulated threshold also included in the output.
# Maximum gap & max. run (max no. of consecutive years unsuitable/suitable) also calculated courtesy of Matt Upson

# The output csv file will be saved with the filename "output_for_threshold_temp-XX_accumulatedDD-YY.csv" and can be imported into ArcGIS for mapping
# Using the JRC map in ArcGIS, the "GRID_NO" field in the output can be linked to the "Grid_code" field (i.e. lat/long are not required)

# remove all pre-existing variables in the workspace
rm(list=ls())

threshold = 10
accumulated_threshold = 500

# load files: put the desired files for analysis into one folder and include the file path for that folder below.
# Multiple files can be read at once.
#allfiles <- list.files(path = "set input file path", full.names = TRUE)
#temperature <- do.call(rbind, lapply(allfiles, read.csv, header = TRUE, sep = ";"))
library(tidyverse)
t1 <- read_delim("../pest-risk/input/jrc-gridded-agrometeo/efsa20180420.csv",delim = ";")

temperature <-t1 %>%
  bind_rows(read_delim("../pest-risk/input/jrc-gridded-agrometeo/efsa20180613.csv",delim = ";")) %>%
  rename(GRID_NO=IDGRID,
         TEMPERATURE_MIN=TMIN,
         TEMPERATURE_MAX=TMAX) %>%
  select(GRID_NO,TEMPERATURE_MIN,TEMPERATURE_MAX,DAY) %>% 
  distinct(GRID_NO,DAY,.keep_all = T)
rm(t1)


# add mean temp to dataset
temperature$mean <- (temperature$TEMPERATURE_MAX + temperature$TEMPERATURE_MIN)/2

# extract the text date string for year and add this to the dataset
temperature$year <- substr(temperature$DAY, 1, 4)

# filter years
#temperature <- temperature %>% 
#  filter(as.numeric(year)<2015) %>% 
#  filter(as.numeric(year)>1999)

# select days where max temperature equal to or below threshold
below <- temperature[temperature$TEMPERATURE_MAX <= threshold, ]
# accumulated degree days for days max temperature below threshold (= zero)
dd_below <- cbind(below, acc_dd = (below$TEMPERATURE_MAX - below$TEMPERATURE_MAX))

# select days where min temperature is equal to or above threshold
above <- temperature[temperature$TEMPERATURE_MIN >= threshold, ]
# accumulated degree days for days min temp above threshold (= average max + min, minus the threshold)
dd_above <- cbind(above,acc_dd = ((above$TEMPERATURE_MAX + above$TEMPERATURE_MIN)*0.5)-threshold)

# select days where the mean temp is greater than or equal to threshold
mean_above <- temperature[temperature$mean >= threshold & 
                            temperature$TEMPERATURE_MIN < threshold & temperature$TEMPERATURE_MAX > threshold, ]
# accummulated degree days where mean is over threshold (= ((max-threshold)/2) - ((threshold-min)/4)
dd_mean_above <- cbind(mean_above, acc_dd = (((mean_above$TEMPERATURE_MAX-threshold)/2)
                                             - ((threshold-mean_above$TEMPERATURE_MIN)/4)))

# select days where the mean is less than threshold
mean_below <- temperature[temperature$mean < threshold &
                            temperature$TEMPERATURE_MIN < threshold & temperature$TEMPERATURE_MAX > threshold,]
# acccumulated degree days where mean less than threshold (= (max-threshold)/4)
dd_mean_below <- cbind(mean_below, acc_dd = ((mean_below$TEMPERATURE_MAX - threshold)/4))
remove(temperature)

# accumulated degree days, all temperatures recombined
dd_temperature <- rbind(dd_below, dd_above, dd_mean_above, dd_mean_below)

# accumulated degree days, by year and by location
final <- aggregate(acc_dd ~ GRID_NO + year, data = dd_temperature, sum)

# some tidying up
# reshaping the dataset
final_tidy <- reshape(final, v.names = "acc_dd", idvar = "GRID_NO", timevar = "year", direction = "wide" )

# sort by grid code
final_sort <- final_tidy[order(final_tidy$GRID_NO),]
final_data <- final_sort[, 2:ncol(final_sort)]

# Replace the default column names with year-only labels
col_labels <- unique(final[,2])
names (final_sort) <- c("GRID_NO", col_labels)

# some data analysis
# Minimum accumulated DD, all years
final_sort$min <- apply(final_data, 1, min, na.rm=TRUE)

# Average all years (mean)
final_sort$mean <- apply(final_data, 1, mean, na.rm=TRUE)

# Maximum accumulated DD, all years
final_sort$max <- apply(final_data, 1, max, na.rm=TRUE)

# Total number of years analysed
final_sort$no_years_analysed <- length(col_labels)

# Number of years over accumulated threshold
function_years_over <- function (x) {
  return (length(which(x >= accumulated_threshold)))
}
final_sort$count_years_over <- apply(final_data, 1, function_years_over)

# Percentage of years over accumulated threshold
final_sort$percent_years_over <- (final_sort$count_years_over/final_sort$no_years_analysed)*100

# Maximum gap analysis - what is the greatest number of years threshold not reached?
# identify which years over and under threshold
function_over <- function(x) {
  return (x >= accumulated_threshold)
}

# apply true/false to whether over or under threshold
gap <- apply(final_sort[,2:(ncol(final_sort)-6)], 2, function_over)

# ==================
# This section of code courtesy of Matt Upson
# ==================

# Define the max_rl function ----

max_rl <- function(x, val) {
  
  # x is the target vector, val is the value for which we want to know the longest
  # run length.
  
  # Test that we are passing a vector a list to the function
  
  if (!is.vector(x) && !is.list(x)) stop("'x' must be a vector of an atomic type")
  
  # Caclulate length of vector
  
  n <- length(x)
  
  # Offset the vector with itself minus the last value
  
  y <- x[-1L] != x[-n]
  
  # Calculate index
  
  i <- c(which(y | is.na(y)), n)
  
  # Calculate length and values as in rle
  
  lengths <- diff(c(0L, i))
  values <-  x[i]
  
  # Determine the max which equals val. Need to first check that val is present
  # in values, if not return 0.
  
  if (val %in% values) {
    
    max_val <- max(lengths[values == val])
    
  } else {
    
    max_val <- 0
    
  }
  
  # This step appears to be necessary to avoid issues when calling map_int later
  
  max_val <- as.integer(max_val)
  
  return(max_val)
  
}

# Now lets wrap this up into a function that can handle dataframes.

max_rl_df <- function(x, val) {
  
  # First check that x is a dataframe or a matrix (note this encompasses data_frames/
  # tibbles).
  
  if (!is.data.frame(x) && !is.matrix(x)) stop("'x' must be a data.frame or matrix")
  
  # Takes a df an val as an argument to pass to max_rl
  
  x <- apply(x, MARGIN = 1, FUN = function(x) max_rl(x, val))
  
  return(x)
  
}

# ============

# Now bind the dataframe to the max length of the TRUEs. 

final_sort <- cbind(
  final_sort,
  max_run = max_rl_df(gap, TRUE),
  max_gap = max_rl_df(gap, FALSE)
)

# Mean annual number of generations possible
final_sort$mean_annual_generation <- final_sort$mean/accumulated_threshold

# export the file: default location as specified, including the two threshold values in the filename
filename <- paste("set output file path",
  threshold, "_accumulatedDD-", accumulated_threshold, "-20y.csv", sep = "")
write.csv(final_sort, file = filename)

