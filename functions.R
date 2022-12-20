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

# Filter a timeline to a single date
# Parameters: a timeline and a date
# Output: a single-date timeline
filterDate <- function(data, date) {
  data[Date == date]
}

# Convert a timeline into an list of single-date timelines
# Parameter: a timeline
# Output: a list of single-date timelines
splitDateData <- function(data) {
  if (nrow(data)) {
    dates <- sort(unique(data[, Date]))
    splitData <- list()
    for (date in dates) {
      splitData <- append(splitData, list(filterDate(data, date)))
    }
  } else { # Data is empty - return a list with one empty dataframe
    splitData <- list(data) 
  }
  return(splitData)
}



#insert a row at a certain index in the data 
# Parameter: the data.table, the index of the row to insert, and the row to be inserted
# Output: the full data.table

insertRow <- function(henData, index, newRow){
  partA = henData[1:index-1,]
  partB = henData[index:dim(henData)[1],]
  combined = rbind(partA,newRow, partB, fill = T)
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

### Feed reactivity

# calculates how much hens react to the feeder run, do they spend more time in feeding zones
# during runs than outside runs
# input: trackingData
# output: daily estimate per hen of feed reactivity

feedReactivity = function(trackingData) {
  
  # take out first and last run because randomisation plus/minus one hour not possible
  #run with 1min delay in different pens but noise projected throughout whole barn 
  #TODO: Take out if hen is in wintergarden during run? 
  #all relevant runs
  freshFeed = as_hms(c("03:59:00","05:59:00", "07:59:00", "9:59:00", "12:59:00", "14:59:00"))
  freshFeed = rbindlist(lapply(freshFeed, function(x){data.table(Time = as_hms(seq(x, x+300)))}))
  
  
  splitHen = splitHenData(trackingData)
  days = unique(trackingData$Date)
  dailyFeed = data.table(Date = days)
  fullFeed = rbindlist(lapply(days, function(x){data.table(Time = ymd_hms(paste(x, as.character(freshFeed$Time))))}))
  fullFeed[, Run := rep(1:6, each = 301, times = length(days))]
  fullFeed[, Obs := rep(1:(length(days)*6), each = 301)]
  
  # randomisation standardised across hens
  # to randomise 50 times around the time point
  # minimum plus/minus 20 min -> 1200 s 
  # maximum plus/minus 60 min -> 3600s,
  randomSelect = matrix(nrow = 50, ncol = max(fullFeed$Obs))
  randomSelect = apply(randomSelect, 2, function(x){x = sample(c(-3600:-1200, 1200:3600), 50)})
  
  for (i in 1:length(splitHen)){
    #create full time series with all seconds of target hen
    cat("Working on hen", i,"\n")
    fullHen = mergeHenData(fillSeqHenData(splitHen[[i]]))
    # extract only the seconds during feed runs
    joinFeed = fullHen[fullFeed, on = "Time"]
    joinFeed[, Date := as.IDate(Time)]
    #delete those feed runs at 4 where the light did not go on at 2
    joinFeed = joinFeed[!(Date < as.IDate("2019-11-22") & Run == 1),]
    #calculate the duration in the feed zone during feed runs
    DurFeedZone = joinFeed[, .(Duration = sum(Zone == "Tier_4" | Zone == "Tier_2")), by = .(Date, Run, Obs)]
    
    #calculate duration in feed zone outside of feed runs 50 times randomised around the time point of the run
    for (j in 1:50){
      #plus/minus to the time by Obs
      fullFeed[, randTime := Time + randomSelect[j,Obs], by = Obs] 
      # extract only the seconds during feed runs
      joinFeedRand = fullHen[fullFeed, on = .(Time == randTime)]
      joinFeedRand[, Date := as.IDate(Time)]
      #delete those feed runs at 4 where the light did not go on at 2
      joinFeedRand = joinFeedRand[!(Date < as.IDate("2019-11-22") & Run == 1),]
      #calculate the duration in the feed zone outside feed runs
      duration = joinFeedRand[, .(Duration = sum(Zone == "Tier_4" | Zone == "Tier_2")), by = .(Date, Run, Obs)][, Duration]
      name = paste0("Eval_", j)
      DurFeedZone[, (name) := (Duration-duration)/(Duration+duration)]
    }
    DurFeedZone[, Total := rowSums(.SD, na.rm = T)/50 , .SDcols = grep("Eval", colnames(DurFeedZone))]
    name = unique(fullHen$Hen)
    dailyFeed[, (name) := DurFeedZone[, mean(Total, na.rm = T), by = Date][,V1]]
    
  }
  return(dailyFeed)
  
}



### Sequence similarity

#calculates similarity of movement sequences between individuals 
#input: trackingdata
# output: data.table containing similarity for each day between each hen pair
similarityBetween = function(data, zone = F, interval = "day"){
  cat("Calculating similarity of movement between individuals on days \n")
  splitDate = splitDateData(data)
  splitHen = vector(mode='list', length= length(unique(data$HenID)))
  allPairs = vector(mode='list', length= length(splitDate))
  comblist <- combn(unique(data$Hen),2,FUN=list)
  cat("To do for", length(unique(data$Hen))*(length(unique(data$Hen))-1)/2, "pairs of individuals and", length(splitDate), "days. \n")
  for (i in 1:length(splitDate)){
    splitHen[[i]] = splitHenData(splitDate[[i]])
    if(interval == "day"){
      start = splitHen[[i]][[1]][LightIndic == T,as.character(as.ITime(Time))][1]
      end = splitHen[[i]][[1]][LightIndic == T,as.character(as.ITime(Time))][2]
    }
    else if(interval == "morning"){
      start = splitHen[[i]][[2]][LightIndic == T,as.character(as.ITime(Time))][1]
      end = "09:00:00"
    }
    for (j in 1: length(splitHen[[i]])){
      splitHen[[i]][[j]] = fillSeqHenData(splitHen[[i]][[j]], start, end)[[1]]
      splitHen[[i]][[j]][, dupl := duplicated(Time)]
      splitHen[[i]][[j]] = splitHen[[i]][[j]][dupl ==F,] 
    }
    fullDate = mergeHenData(splitHen[[i]])
    fullDate[, TimePure := as.ITime(Time)]
    date = fullDate[,unique(Date)]
    names(allPairs)[i] = as.character(date)
    #save all possible pair combinations
    allPairs[[i]] = data.table(Pair1 = combn(unique(data$Hen),2)[1,], 
                               Pair2 = combn(unique(data$Hen),2)[2,])
    #if a certain zone should be isolated
    if(zone != F){
      fullDate[Zone == zone, ZoneIso := zone]
      fullDate[Zone != zone, ZoneIso := "Other"]
      fullDate = dcast(fullDate, formula = TimePure ~ Hen, value.var = "ZoneIso")
      #if no zone needs to be isolated
      } else{
        fullDate = dcast(fullDate, formula = TimePure ~ Hen, value.var = "Zone")
      }
    if (dim(fullDate)[2]-1 != length(unique(data$Hen))){
        hens = unique(data$Hen)[!(unique(data$Hen) %in% colnames(fullDate)[-1])]
        for(hen in hens){
          fullDate[, (hen) := NA]
        }
    }
    #calculate the similarity by dividing the amount of zones shared at the same time 
    #by the total amount of time
    allPairs[[i]]  = cbind(allPairs[[i]], t(fullDate[,lapply(comblist,
                                         function(x) sum(get(x[1]) == get(x[2]))/length(get(x[1])))]))
    allPairs[[i]] = rbind(allPairs[[i]], data.table(Pair1 = unique(data$Hen),
                                                    Pair2 = unique(data$Hen)), fill = T)
    allPairs[[i]] = as.matrix(dcast(allPairs[[i]], Pair1 ~ Pair2, value.var = "V1", drop = F)[, -1])
    rownames(allPairs[[i]]) = colnames(allPairs[[i]]) 
    cat("Done for", i, "days. \n")
    }
  
  return(allPairs)
}

#calculates similarity of movement sequences within individuals
# input: trackingdata
# output: data.table with similarity for consecutive measure days for each hen
#TODO: so far cannot compare days that are non-consecutive
#idea: maybe fill all days until they include all times starting at 2? 

similarityWithin = function(data, interval = "day"){
  cat("Calculating similarity of movement within individuals between days \n")
  splitHen = splitHenData(data)
  fullDate = vector(mode='list', length= length(unique(data$Date)))
  allDates = data.table(Date = unique(data$Date))
  cat("To do for", length(splitHen), "individuals and", dim(allDates)[1], "days. \n")
  for (i in 1:length(splitHen)){
    
    fullHen = mergeHenData(fillSeqHenData(splitHen[[i]]))
    fullHen[, dupl := duplicated(Time)]
    fullHen = fullHen[dupl == F,]
    if(interval == "day"){
      fullHen = fullHen[Light == T,]
    }
    
    entry = fullHen[ , ZoneComp := Zone == Zone[match(Time - 24*3600, Time)]][
      , (sum(ZoneComp == T, na.rm = TRUE)/(sum(ZoneComp == T, na.rm = TRUE)+ sum(ZoneComp == F, na.rm = TRUE))), by = .(Date = as.Date(Time))]
    names(entry)[2] = unique(fullHen$Hen)
    allDates = entry[allDates, on = "Date"]
  cat("Done for", i, "individuals. \n")
  }
  return(allDates)
}


#create randomised transition sequences

#1. only shuffle within day within individual
#2. only shuffle within time span of pophole open/closed
#3. control for unlikely transitions
#4. change also timepoint of transition?
#5. what if bird stays on the same zone 