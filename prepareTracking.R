source("functions.R")

#Function that takes tracking data input and outputs relevant data
# input: trackingData = data table containing tracking data
#        hens         = vector containing relevant animals 

prepareTrackingdata <- function(trackingData, hens){

  #add numerical hen identifier to tracking data
  trackingData[, HenID := unlist(regmatches(Hen, gregexpr('\\(?[0-9,.]+', Hen)))]
  trackingData[, HenID := as.numeric(HenID)]
  #extract only relevant hens
  trackingData = trackingData[HenID %in% hens,]
  
  trackingData[, Time := ymd_hms(Time)]
  trackingData[, Date:= as_date(Time)]
  trackingData[, Pen := unlist(regmatches(PackID, gregexpr('\\(?[0-9,.]+', PackID)))]
  
  trackingData = trackingData[!(Date < as.IDate("2019-11-08")),]
  
  ###### Data preparations ################
  
  #data cleaning step
  #remove Wintergarten entries outside time that Wintergarten can be open
  #delete all Wintergarten entries at night (not possible, tracking system errors)
  trackingData = trackingData[!(hour(Time) > 16 & Zone == "Wintergarten"), ]
  #delete all Wintergarten entries before 10 in the morning (not possible, tracking system error)
  trackingData = trackingData[!(hour(Time) < 10 & Zone == "Wintergarten"), ]
  
  
  #include light-dark cycle
  trackingData[, Light := TRUE]
  trackingData[hour(Time) > 16, Light := FALSE]
  trackingData[month(Time)== 11 & day(Time) < 15 & hour(Time) < 4, Light := FALSE]
  trackingData[month(Time)== 11 & day(Time) < 22 & day(Time) > 14 & hour(Time) < 3, Light := FALSE]
  trackingData[(month(Time)!= 11| (month(Time)== 11 & day(Time) > 21)) & hour(Time) < 2, Light := FALSE]
  #time shift on 29.03. (in tracking system) but done on the 30th for hens (delete days)
  
  
  # mark true transitions
  trackingData[, TrueTransition := TRUE]
  # mark transitions as not indicating day start and end
  trackingData[, DayIndic := F]
  # mark transitions as not indicating lights on and off
  trackingData[, LightIndic := F]
  
  #splitting data into Hens
  splitHen = splitHenData(trackingData)
  
  #add beginning and end of each day by looping through every hen and every day and adding the stamps
  hen_list <- vector(mode='list', length=length(unique(trackingData$hens)))
  days = as.character(unique(trackingData$Date))
  
  cat("Inserting day and light indicators \n")
  cat("To do for", length(hens), "individuals and", length(days), "days. \n")
  for (i in 1: length(hens)) {
    hen_list[[i]] <- vector(mode='list', length=length(days))
    j = 1
    for (day in days){
      start = ymd_hms(paste(day, "00:00:00"))
      end = ymd_hms(paste(day, "23:59:59"))
      relDat = extractInterval(splitHen[[i]], start, end)
      #to keep track if an entry already exists at lights on
      start = F
      
      #query to figure out the time the lights went on
      if (nrow(relDat)){
        relDat[c(1,.N), `:=`( Light = F, TrueTransition = F, DayIndic = T, LightIndic = F)]
        if( month(day) == 11 & day(day) < 15){
          if(any(hour(relDat$Time) == 4 & minute(relDat$Time) == 0 & second(relDat$Time) == 0)){
            start = T
          } else{
            dayStart = "04:00:00"
          }
          
        } else if (month(day)== 11 & day(day) < 22 & day(day) >14){
          if(any(hour(relDat$Time) == 3 & minute(relDat$Time) == 0 & second(relDat$Time) == 0)){
            start = T
          } else{
            dayStart = "03:00:00"
          }
        } else{
          if(any(hour(relDat$Time) == 2 & minute(relDat$Time) == 0 & second(relDat$Time) == 0)){
            start = T
          } else{
            dayStart = "02:00:00"
          }
        }
        #inserting the lights on stamp
        if(start == F){
          #use first indication of lights as index
          indexStart = which(relDat$Light == T)[1]
          startEntry = relDat[indexStart-1,]
          startEntry[, Time := ymd_hms(paste(day, dayStart))]
          startEntry[, Light := T]
          startEntry[, TrueTransition := F]
          startEntry[, DayIndic := F]
          relDat = insertRow(relDat, indexStart, startEntry)
        } 
        #inserting the lights off stamp
        if(!any(hour(relDat$Time) == 17 & minute(relDat$Time) == 0 & second(relDat$Time) == 0)){
          dayEnd = "17:00:00"
          #use last entry before lights out time
          indexEnd = which(hour(relDat$Time) > 17)[1]
          endEntry = relDat[indexEnd-1,]
          endEntry[, Time := ymd_hms(paste(day, dayEnd))]
          endEntry[, Light := F]
          endEntry[, TrueTransition := F]
          endEntry[, DayIndic := F]
          relDat = insertRow(relDat, indexEnd, endEntry)
        }
        
        relDat[which(Light == T)[1], LightIndic := T]
        relDat[which(hour(Time) == 17)[1], LightIndic := T]
        
        hen_list[[i]][[j]] = relDat
        j = j+1
       }
    }
    hen_list[[i]] = mergeHenData(hen_list[[i]])
    cat("Done for", i, "hens. \n")
  } 
  trackingData = mergeHenData(hen_list)
  
  #remove unnecessary data from workspace
  rm(splitHen)
  
  #make sure date is correct (now with inserting of end and start of day)
  trackingData[, Date:= as_date(Time)]

  
  ###### data checks #############################
  
  allDays = seq(ymd("2019-11-09"), to = ymd("2020-06-28"), by = "day")
  
  #select only relevant days
  trackingData = trackingData[Date %in% allDays, ]
  #there is data for every day?
  length(unique(trackingData$Date)) == length(allDays)
  #which days are missing?
  allDays[which(!(allDays %in% unique(trackingData$Date)))]
  #missing days:
  #"2019-12-05" "2019-12-22" "2019-12-27" "2019-12-28" "2019-12-29" "2019-12-30" 
  #"2019-12-31" "2020-01-01" "2020-01-02" "2020-01-03" "2020-01-04" "2020-05-19" "2020-05-24"
  
  #days to exclude
  excl = c(ymd("2019-12-09"), ymd("2020-01-06"), ymd("2020-02-10"),
           ymd("2020-03-16"), ymd("2020-05-04"), ymd("2020-06-02"), #health assessments
           ymd("2019-12-04"), ymd("2019-12-06"), ymd("2019-12-21"),
           ymd("2019-12-26"), ymd("2020-01-12"), ymd("2020-01-20"), ymd("2020-06-20"), #missing data
           seq(ymd("2020-01-21"), to = ymd("2020-02-03"), by = "day"), #different configuration
           seq(ymd("2020-05-17"), to = ymd("2020-05-26"), by = "day"), # experimental light change
           ymd("2020-03-29"), ymd("2020-03-30"), #time shift winter -> summer
           ymd("2020-01-07"),ymd("2020-02-25"),ymd("2020-03-11") #other disturbances
  )
  
  trackingData = trackingData[!(Date %in% excl),]

  length(unique(trackingData$Date)) #-> now we have 180 days in total of data
  
  #do all individuals have data?
  length(unique(trackingData$HenID)) == 36
  #all individuals have data every day?
  entriesPerDay = trackingData[, length(unique(Date)), by = HenID][V1 != 180,][, miss := 180-V1]
  
  #which days are missing for those hens?
  days = unique(trackingData$Date)
  
  trackingData[HenID == 1, days[which(!(days %in% Date))]] #39 days miss, no data between the 22.11. - 26.12. + 05.01. + 09.- 20.01.
  trackingData[HenID == 17, days[which(!(days %in% Date))]] # 2 days miss, no data on the 23.12. & 05.01.
  trackingData[HenID == 33, days[which(!(days %in% Date))]] # 1 day miss, no data on the 17.03.
  trackingData[HenID == 35, days[which(!(days %in% Date))]] # 1 day miss, no data on the 23.12.
  trackingData[HenID == 39, days[which(!(days %in% Date))]] # 1 day miss, no data on the 05.01.
  trackingData[HenID == 65,days[which(!(days %in% Date))]] #7 day miss, no data between 2.12. - 08.12. + 23.12. - 26.12.
  trackingData[HenID == 84,days[which(!(days %in% Date))]] #2 days miss, no data between 07.12. and 08.12.
  
  
  #exclude identified problem days
  trackingData = trackingData[!(Date == as.Date("2020-01-05")|
                                Date == as.Date("2019-12-23")),]
  entriesPerDay = trackingData[, length(unique(Date)), by = HenID][V1 != 178,][, miss := 178-V1]
  length(unique(trackingData$Date)) #-> now we have 178 days in total of data
  #TODO: what to do with hens with missing days? (especially HenID 1 (37 miss) & 65 (6 miss))
  
  ###### last data additions and checks ######
  
  #catching Wintergarten errors of inserting time stamps
  trackingData[(LightIndic == T | DayIndic == T)& Zone == "Wintergarten", Zone := "Litter"]
  
  
  #add duration
  trackingData[, Duration := (shift(Time, type="lead") - Time), by = HenID]
  #set all last day entries to 1 sec so that missing days don't count into the duration
  trackingData[hour(Time)== 23 & minute(Time)== 59 & second(Time) == 59, Duration := 1]
  
  #indicate which night times belong together
  trackingData[Light == T, NightCycle := Date]
  trackingData[hour(Time) > 17|hour(Time) == 17 , NightCycle := Date]
  trackingData[hour(Time) < 4 & Light == F, NightCycle := Date-1]
  #delete duplicated entries
  trackingData = trackingData[Duration > 0,]
  
  return(trackingData)
  
}
