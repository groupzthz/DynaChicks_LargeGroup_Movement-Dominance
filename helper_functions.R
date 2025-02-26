### helper functions

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
fillSeqHenData <- function(data, start = "00:00:00", end = "23:59:59", interval = "sec") {
  if (nrow(data)) {
    days <- as.character(unique(data[, Date]))
    splitData <- list()
    for (day in days) {
      seq = data.table(Time = seq(ymd_hms(paste(day, start)), 
                                  ymd_hms(paste(day, end)), by = interval))
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


#insert a row at a certain index in the data 
# Parameter: the data.table, the index of the row to insert, and the row to be inserted
# Output: the full data.table

insertRow <- function(henData, index, newRow){
  partA = henData[1:index-1,]
  partB = henData[index:dim(henData)[1],]
  combined = rbind(partA,newRow, partB, fill = T)
}



### Vertical travel distance 

#Return the distance traveled between two zones
#Parameters: vector containing two zones to compare
defineDistance = function(zones) {
  
  zone1 = zones[1]
  zone2 = zones[2]
  
  if((zone1 == zone2)|
     is.na(zone1) | is.na(zone2)){
    distance = 0
  } 
  else if((zone1 == "Tier_4" & zone2 == "Ramp_Nestbox")|
     (zone2 == "Tier_4" & zone1 == "Ramp_Nestbox")|
     (zone1 == "Tier_2" & zone2 == "Ramp_Nestbox")|
     (zone2 == "Tier_2" & zone1 == "Ramp_Nestbox")|
     (zone1 == "Tier_2" & zone2 == "Litter")|
     (zone1 == "Litter" & zone2 == "Tier_2")){
    distance = 1
  }
  else if((zone1 == "Tier_4" & zone2 == "Tier_2")|
          (zone2 == "Tier_4" & zone1 == "Tier_2")|
          (zone1 == "Litter" & zone2 == "Ramp_Nestbox")|
          (zone2 == "Litter" & zone1 == "Ramp_Nestbox")
          ){
    distance = 2
  }
  else if((zone1 == "Tier_4" & zone2 == "Litter")|
          (zone2 == "Tier_4" & zone1 == "Litter")){
    distance = 3
  }
  else {
    stop("Error: the zones contained an unknown zone")
  }
  return(distance)
}


## Feeder activity

#Feed zone presence during feed chain run
#Parameter: tracking data.table
#Output: daily presence within feed zones during feed chain runs by individual
infeedZone = function(trackingData) {
  
  # is hen in any of the two feeder zones during the run and if so in which?
  #run with 1min delay in different pens but noise projected throughout whole barn 
  #all relevant runs
  freshFeed = as_hms(c("01:59:00","03:59:00","05:59:00", "07:59:00", "9:59:00", "12:59:00", "14:59:00", "16:14:00"))
  freshFeed = rbindlist(lapply(freshFeed, function(x){data.table(Time = as_hms(seq(x, x+300)))}))
  
  #split data by hen
  splitHen = splitHenData(trackingData)
  days = unique(trackingData$Date)
  #create empty frame for data
  dailyFeed = data.table()
  #create list with feeding events
  fullFeed = rbindlist(lapply(days, function(x){data.table(Time = ymd_hms(paste(x, as.character(freshFeed$Time))))}))
  #create identifiers for each feeder run time
  fullFeed[, Run := rep(1:8, each = 301, times = length(days))]
  fullFeed[, Obs := rep(1:(length(days)*8), each = 301)]
  #delete all feeder runs before the 21.11. at 2 am
  fullFeed = fullFeed[!(hour(Time) < 3 & day(Time) <21 & month(Time) == 11),]

  for (i in 1:length(splitHen)){
    #create full time series with all seconds of target hen
    cat("Working on hen", i,"\n")
    fullHen = mergeHenData(fillSeqHenData(splitHen[[i]]))
    # extract only the seconds during feed runs
    joinFeed = fullHen[fullFeed, on = "Time"]
    joinFeed[, Date := as.IDate(Time)]
    #calculate the duration in the feed zone during feed runs
    DurFeedZone = joinFeed[, .(Duration2 = sum(Zone == "Tier_2"),
                               Duration4 = sum(Zone == "Tier_4")), by = .(Date, Run, Obs)]
    name = unique(fullHen$Hen)
    dailyFeed = rbind(dailyFeed, DurFeedZone[, .(Total2 = sum(Duration2),
                                      Total4 = sum(Duration4),
                                      Runs = length(Run),
                                      Hen = name), by= Date][,NotFeedZone := Runs*301-(Total2+Total4), by = Date])
    
  }
  return(dailyFeed)

}


# Reorder correlation matrix
reorder_cormat <- function(cormat){
  # Use correlation between variables as distance
  dd <- as.dist((1-cormat)/2)
  hc <- hclust(dd)
  cormat <-cormat[hc$order, hc$order]
}

# Get upper triangle of a correlation matrix
get_upper_tri <- function(cormat){
  cormat[lower.tri(cormat)]<- NA
  return(cormat)
}


