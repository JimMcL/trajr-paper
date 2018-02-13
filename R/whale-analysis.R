# Whale trajectory analysis.
#
# This file defines a number of functions which read in and analyse multiple
# whale trajectories from CSV files. See generate-figures-report.R for how the
# functions may be used.

library(rgdal)
library(sp)
library(ggmap)
library(trajr)


# Some basic utility functions

# Converts radians to degrees
rad2deg <- function(rad) rad * 180 / pi

# Converts metres/sec to knots (i.e. 1852 m / hour)
msToKnts <- function(ms) 60 * 60 * ms / 1852

# Converts metres/sec to km/h
msToKmh <- function(ms) 60 * 60 * ms / 1000

# Returns the month that the trajectory was recorded
trjMonth <- function(trj) as.numeric(strftime(trj$date[1], "%m"))

# Returns TRUE if the trajectory was part of the northern migration (i.e. May - Aug)
trjNorthernMigration <- function(trj) trjMonth(trj) < 9

# Returns the duration of the trajectory.
# Implemented as TrajDuration in trajr v1.0.1
trjDuration <- function(trj, startIdx = 1, endIdx = nrow(trj)) {
  trj$displacementTime[endIdx]
}

# Returns the overall velocity vector of a trajectory, units are metres / sec.
# Value, v, is returned as a complex number, Re(v) is x-component, Im(v) is y-component.
# Implemented as TrajMeanVelocity in trajr v1.0.1
trjVelocity <- function(trj, startIdx = 1, endIdx = nrow(trj)) {
  d <- (trj[endIdx, c("x", "y")] - trj[startIdx, c("x", "y")]) / trjDuration(trj, startIdx, endIdx)
  complex(re = d[1], im = d[2])
}

# Returns mean speed along the trajectory.
trjSpeed <- function(trj) {
  TrajLength(trj) / trjDuration(trj)
}

# Returns the colour to plot a trajectory
NORTH_COL <- "black"
SOUTH_COL <- "red"
colForTrj <- function(trj) ifelse(trjNorthernMigration(trj), NORTH_COL, SOUTH_COL)

# Returns the line style to plot a trajectory
NORTH_LTY <- 2
SOUTH_LTY <- 1
ltyForTrj <- function(trj) ifelse(trjNorthernMigration(trj), NORTH_LTY, SOUTH_LTY)

# Basic long/lat projection
wgs84.proj4 <- "+proj=longlat +ellps=WGS84 +datum=WGS84"

# Albers equal area projection is a square coordinate system, 
# with these parameters, suitable for use in Australia
albers.proj4 <- "+proj=aea +lat_1=-18 +lat_2=-36 +lat_0=0 +lon_0=134 +x_0=0 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs";



# Whale trajectories are in CSV format with a row for each point and the
# following columns: Pod, Latitude, Longitude, Time, Behaviour,
# Closest.boat..m., Boats...100, Boats...300, Boats...300.1. Since latitude and
# longitude are not square coordinates, they must be transformed before use.
# This is a function that reads in the CSV file and converts everything to units
# that can be used by `trajr`
ReadWhaleCSV <- function(filename, ...) {
  # Read in the file
  # Handle comments, even though they are not part of CSV
  coords <- read.csv(filename, comment.char = "#", stringsAsFactors = FALSE)
  
  # Parse time column to a POSIXlt object
  coords$date <- strptime(coords$Time, "%Y-%m-%d %H:%M:%S")
  # Also convert date/time to seconds for easy comparison and duration calculations
  coords$Time <- as.numeric(coords$date)
  
  # --- Project to a square coordinate system, in this case, Albers equal area
  # 1. Create a SpatialPointsDataFrame with just the longitude and latitude
  pts <- SpatialPoints(coords[, c("Longitude", "Latitude")], 
                       proj4string=CRS("+proj=longlat"))
  # 2. Project to Albers equal-area
  pts <- sp::spTransform(pts, CRS(albers.proj4))
  # Stick x & y columns onto the original coordinates
  coords$x <- coordinates(pts)[, 1]
  coords$y <- coordinates(pts)[, 2]
  # ---
  
  # data cleaning - since observations should be monotonically increasing in
  # time, get rid of any entry which is not later than the previous entry
  coords <- coords[diff(c(0, coords$Time)) > 0, ]
  
  # Each CSV file may contain more than 1 pod. Return a  list containing 1 track
  # for each pod
  r <- lapply(unique(coords$Pod), function(pod) {
    coords[coords$Pod == pod, ]
  })
  
  # Discard tracks with less than 30 coordinates
  r <- r[sapply(r, nrow) > 30]
  
  r
}


# Function to find, read in and convert all of the whale trajectory files
ReadWhaleTrajectories <- function() {
  # Find all files called "trajectories*.csv"
  files <- list.files("..", "trajectories-.*.csv", recursive = TRUE, full.names = TRUE)
  # Read them in, convert them to trajectories, then smooth them 
  trjs <- TrajsBuild(files, spatialUnits = "m", 
                     csvReadFn = ReadWhaleCSV, 
                     csvStruct = list(x = "x", y = "y", time = "Time"),
                     smoothP = 3, smoothN = 21)
  
  # Discard trajectories which were observed for less than 15 minutes
  trjs <- trjs[(sapply(trjs, function(trj) trj[nrow(trj),]$displacementTime) / 60) >= 15]
  
  trjs
}


# Calculate indices for each trajectory
CalculateWhaleIndices <- function(trjs) {

  # Define a function to calculate whatever values we need for a single trajectory
  statsFn <- function(trj) {
    duration <- trjDuration(trj)
    meanSpeed <- TrajLength(trj) / duration
    
    # A trajectory consists of possibly many observations during each up-time,
    # and no observations during a downtime. Assume a dive is any down time
    # longer than 120 seconds, following Gulesserian (2011)
    diveCutoff <- 120
    stepTimes <- diff(trj$displacementTime)
    downTimes <- stepTimes[stepTimes > diveCutoff]

    # Items of interest are:
    list(
      # Northern or southern migration
      isNorthern = trjNorthernMigration(trj),
      # Were there ever more than 3 boats within 300 m
      Too.many.boats = max(trj$Boats...300 + trj$Boats...100) > 3,
      # Number of breaches / sec
      Breach.freq = length(grep(".*breach.*", trj$Behaviour, ignore.case = TRUE)) / duration,
      # Number of tail swipes / sec
      Tail.swipe.freq = length(grep(".*swipe.*", trj$Behaviour, ignore.case = TRUE)) / duration,
      # Mean speed (km/h) along entire trajectory
      Mean.speed = msToKmh(meanSpeed),
      # Mean downtime
      Mean.downtime = ifelse(length(downTimes) == 0, 0, mean(downTimes)),
      # Straightness
      Straightness = TrajStraightness(trj),
      # Sinuosity
      Sinuosity = TrajSinuosity2(trj)
    )
  }
  
  # Calculate movement indices for each trajectory 
  TrajsMergeStats(trjs, statsFn)
}

# Plots all trajectories on a map. Trajectories are coloured based on the migration direction.
PlotWhaleTrajectories <- function(trjs) {
  # Draw a map of the area. Hardwire zoom level because auto doesn't do a good job
  xlim <- range(sapply(trjs, function(trj) range(trj$Longitude)))
  ylim <- range(sapply(trjs, function(trj) range(trj$Latitude)))
  m <- get_map(location = c(mean(xlim), mean(ylim)), zoom = 11, messaging = FALSE)
  m <- ggmap(m, messaging = FALSE, legend = "none") +
    xlab("Longitude") + ylab("Latitude")
  # Add the trajectories to the map
  for(trj in trjs) {
    d <- trj
    d$Direction <- ifelse(trjNorthernMigration(trj), "Northern migration", "Southern migration")
    m <- m + geom_path(data = d, aes(x = Longitude, y = Latitude, color = Direction))
  }
  m <- m + scale_color_manual(values=c(NORTH_COL, SOUTH_COL))
  plot(m)
}

# Prints a summary of mean vector of trajectories for northern and southern migrations
SummariseWhalesByDirection <- function(trjs, stats) {
  # Calculate average northern and southern velocities in m / sec
  isNorthern <- sapply(trjs, trjNorthernMigration)
  meanNorthernVelocity <- mean(sapply(trjs[isNorthern], trjVelocity))
  meanSouthernVelocity <- mean(sapply(trjs[!isNorthern], trjVelocity))
  
  nn <- sum(isNorthern)
  ns <- sum(!isNorthern)
  ntm <- sum(stats[isNorthern, "Too.many.boats"])
  stm <- sum(stats[!isNorthern, "Too.many.boats"])
  cat(sprintf("Total %d trajectories, %d north, %d south\n", length(trjs), nn, ns))
  cat(sprintf("Number of trajectories approached by > 3 boats within 300m: %d north (%g%%), %d south (%g%%)\n",
    ntm, 100 * ntm / nn, stm, 100 * stm / ns))
  cat(sprintf("Mean velocities:\n    Northern migration %g km/h, direction ~%g degrees\n    Southern migration %g km/h, direction ~%g degrees\n", 
              signif(msToKmh(Mod(meanNorthernVelocity)), 2),
              signif(90 - rad2deg(Arg(meanNorthernVelocity)), 2),
              signif(msToKmh(Mod(meanSouthernVelocity)), 2),
              signif(90 - rad2deg(Arg(meanSouthernVelocity)), 2)))
  
  
  meanNorthernSpeed <- mean(sapply(trjs[isNorthern], trjSpeed))
  meanSouthernSpeed <- mean(sapply(trjs[!isNorthern], trjSpeed))
  cat(sprintf("Mean swimming speeds: north %g km/h, south %g km/h", 
              round(msToKmh(meanNorthernSpeed), 1), round(msToKmh(meanSouthernSpeed), 1)))
}

# Performs a t-test comparing term to every behavioural variable, looking for
# statistical significant comparisons after adjustment for multiple comparisons.
# Returns a list of statistically significant t-tests.
CompareWhaleTrajectories <- function(stats, term, alpha = 0.05, correction = "holm") {
  
  prettify <- function(str) gsub("freq", "freq.", gsub("\\.", " ", str))
  
  behaviouralVars <- names(stats)[3:ncol(stats)]
  t <- lapply(behaviouralVars, function(response) t.test(reformulate(term, response), stats))
  # Adjust p values for multiple comparisons (reduce type I errors)
  p <- sapply(t, function(tt) tt$p.value)
  p <- p.adjust(p, correction)
  names(p) <- prettify(behaviouralVars)
  
  # Save t values
  t.st <- sapply(t, function(tt) tt$statistic)
  names(t.st) <- prettify(behaviouralVars)
  
  diffs <- t[p < alpha]

  # Return adjusted p value and list of significant test results
  list(p = round(p, 2), t = round(t.st, 2), diffs = diffs)
}

ReportAllWhaleParams <- function(trjs) {
  params <- CalculateWhaleIndices(trjs)
  # Add a column for date-AM/PM
  params$Date <- sapply(trjs, function(trj) strftime(trj$date[1], "%Y-%m-%d %P"))
  names(params) <- gsub("\\.", " ", names(params))
  print(params)
}

# Print a report of various whale statistics
ReportWhaleStats <- function(trjs) {
  # summarise t test result
  .ttPrint <- function(tt) {
    group <- gsub("\\.", " ", sub(".* by ", "", tt$data.name))
    cat(sprintf("%s: %s\n    t = %g, df = %g, unadjusted p-value = %g\n    Means: %s %g, not %s %g\n", 
                gsub(" .*", "", tt$data.name), tt$method,
                round(tt$statistic, 1), round(tt$parameter, 1), round(tt$p.value, 3),
                group, signif(tt$estimate[2], 1), group, signif(tt$estimate[1], 1)))
  }
  .ttlPrint <- function(ttl) sapply(ttl, .ttPrint)
  
  whaleIndices <- CalculateWhaleIndices(trjs)

  cat(sprintf("-----------------------------------------------------------------\nSummary of trajectories:\n"))
  SummariseWhalesByDirection(trjs, whaleIndices)
  
  alpha <- 0.05
  correction <- "holm"
  #correction <- "none"
  cat(sprintf("\n-----------------------------------------------------------------\n"))
  cat(sprintf("Significant behavioural differences between northern and southern migrations.\n"))
  cat(sprintf("%s correction for multiple comparisons.\n", correction))
  .ttlPrint(CompareWhaleTrajectories(whaleIndices, "isNorthern", alpha = alpha, correction = correction)$diff)
  
  north <- CompareWhaleTrajectories(whaleIndices[whaleIndices$isNorthern,], "Too.many.boats", alpha = alpha, correction = correction)
  south <- CompareWhaleTrajectories(whaleIndices[!whaleIndices$isNorthern,], "Too.many.boats", alpha = alpha, correction = correction)
  
  cat(sprintf("\n=================================================================\n"))
  cat(sprintf("Significant behavioural differences between 4 or more boats approaching within 300 m\n"))
  cat(sprintf("%s correction for multiple comparisons.\n", correction))
  
  # Write out a table of adjusted p values which can be copied and pasted into a Word document
  cat("Adjusted p values\n\n")
  pTable <- cbind(Direction = c("Northern", "Southern"), rbind(north$p, south$p))
  pTable <- ifelse(pTable > alpha, pTable, paste0(pTable, "*"))
  write.table(pTable, sep = "\t", quote = FALSE, row.names = FALSE)
  cat("\n")
  
  cat(sprintf("-----------------------------------------------------------------\nNorthern migration:\n"))
  .ttlPrint(north$diffs)
  cat(sprintf("-----------------------------------------------------------------\nSouthern migration:\n"))
  .ttlPrint(south$diffs)
  
  invisible()
}
