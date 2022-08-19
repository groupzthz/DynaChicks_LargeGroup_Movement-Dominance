### functions KG

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
#TODO: make this function work
fillSeqHenData <- function(data) {
  if (nrow(data)) {
    days <- sort(unique(data[, Date]))
    splitData <- list()
    for (day in days) {
      seq = seq(ymd_hms(paste(day, "00:00:00")), ymd_hms(paste(day, "23:59:59")), by = "sec")
      #continue here with merging 
      splitData <- append(splitData, list(filterHen(data, henID)))
    }
  } else { # Data is empty - return a list with one empty dataframe
    splitData <- list(data) 
  }
  return(splitData)
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
