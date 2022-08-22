### functions KG
library("zoo")  # for na.locf function

# Extract a specific time period from a timeline
# Parameters: a single-hen timeline
# Output: a filtered timeline between start and end, with added records at start and end exactly
extractInterval <- function(henData, start, end) {
  stopifnot(start < end)
  prevRecord <- henData[findInterval(start, henData[,Time])]
  lastRecord <- henData[findInterval(end, henData[,Time])]
  inside <- henData[between(Time, start, end)]
  if (nrow(inside)) {
    if (nrow(prevRecord)) {
      # We have a record before start time - take it, but truncate to start time
      prevRecord[,Time := start]
    } else {
      # We don't have a record - duplicate first row with start time and NA zone
      prevRecord <- first(inside)
      prevRecord[,Time := start]
      prevRecord[,Zone := NA]
    }
    
    if (lastRecord$Time == end) {
      return(rbind(prevRecord, inside, fill=TRUE))
    } else {
      lastRecord[,Time := end]
      return(rbind(prevRecord, inside, lastRecord, fill=TRUE))
    }
  } else {
    return(inside)
  }
}


# take a timeline and fill all seconds 
# Parameters: a single hen timeline
# Output: a list of full sequence days for one hen
fillSeqHenData <- function(data) {
  if (nrow(data)) {
    days <- as.character(unique(data[, Date]))
    splitData <- list()
    for (day in days) {
      seq = data.table(Time = seq(ymd_hms(paste(day, "00:00:00")), 
                                  ymd_hms(paste(day, "23:59:59")), by = "sec"))
      #extract data for that day 
      relDat = data[Date == day,]
      #merge timelines together
      relDat = relDat[seq, on = "Time"]
      relDat = na.locf(relDat, na.rm = FALSE)
      splitData <- append(splitData, list(relDat))
    }
  } else { # Data is empty - return a list with one empty dataframe
    splitData <- list(data) 
  }
  return(splitData)
}

#take a list of full sequence days and extract most common zone for time of day
#Parameters: a list of full sequence days for one hen
#Output: a data.table with the most common zone per time for this hen
averageDailyZone <- function(data){
    nDays <- length(data)
    hen = unique(data[[ndays]][,HenID])
    vectorTime = data.table(Time = as_hms(seq(ymd_hms("2000-01-01 00:00:00"), ymd_hms("2000-01-01 23:59:59"), by = "sec")))
    for(ndays in 1:nDays){
      day = unlist(data[[ndays]][1,2])
      zones = data[[ndays]][,Zone]
      vectorTime = cbind(vectorTime, zones)
      colnames(vectorTime)[ndays+1] = paste0("zones_",ndays)
    }
    #extract most common zone for corresponding time
     averageZone = melt(vectorTime, measure = patterns('zones'))[, .(Zone = Mode(value)), by = Time]
     averageZone[, HenID := hen]
     averageZone = rbind(sparseSeqHenData(averageZone),averageZone[.N,.(Time,Zone,HenID)])

  return(averageZone)
  
}

#returns most common element in vector 
Mode <- function(x) {
      ux <- unique(x)
     ux[which.max(tabulate(match(x, ux)))]
}

#delete unnecessary entries again
#Parameters: a full sequence of entries 
sparseSeqHenData <- function(data){
  data[, Dupl := duplicated(rleid(Zone))]
  sparseData = data[Dupl == F,c(1,2,3)]
  return(sparseData)
}

### Data-wrangling functions

# Filter a timeline to a single hen
# Parameters: a timeline and a hen ID
# Output: a single-hen timeline
filterHen <- function(data, henID) {
  data[HenID == henID]
}

# Convert a timeline into an list of single-hen timelines
# Parameter: a timeline
# Output: a list of single-hen timelines
splitHenData <- function(data) {
  if (nrow(data)) {
    hens <- sort(unique(data[, HenID]))
    splitData <- list()
    for (henID in hens) {
      splitData <- append(splitData, list(filterHen(data, henID)))
    }
  } else { # Data is empty - return a list with one empty dataframe
    splitData <- list(data) 
  }
  return(splitData)
}

# Merge a list of timelines into a single timeline
# Parameter: a list of timelines
# Output: a single timeline
mergeHenData <- function(henDataList) {
  combined <- rbindlist(henDataList)
  setkeyv(combined, c("Time", "Hen"))
  setindex(combined, Hen)
  return(combined)
}




### Duration-related functions

# Add a Duration column to a timeline
# Calculated for each row as the time difference to next row; 0 for the last row
# Parameter: a single-hen timeline
# Output: a timeline with durations
addDurations <- function(henData) {
  newData <- copy(henData)
  # Calculate duration of entry as difference to next entry
  newData[, Duration := (shift(Time, type="lead") - Time)]
  # Set duration for last entry as 0
  newData[nrow(newData), Duration := as.difftime(0, units="secs")]
  return(data.table(newData))
}

# Return a table of zone / cumulative duration for the zone in a timeline
# Parameters: timeline data, and OPTIONAL list of zones (to display 0 in missing ones)
sumDurations <- function(data, zones) {
  # Add durations
  dataD <- applyPerHen(addDurations, data)
  # Sum durations by zone
  sums <- dataD[, list(TotalDuration=sum(Duration)), by=c("Zone")]
  if (!missing(zones)) {
    zonesDT <- data.table(Zone = zones)
    # Replace zone list with "zones" argument
    sums <- sums[zonesDT, on="Zone"]
    # Fill in NA for missing zones
    sums[is.na(TotalDuration), TotalDuration := as.difftime(0, units="secs")]
  }
  return(sums)
}
