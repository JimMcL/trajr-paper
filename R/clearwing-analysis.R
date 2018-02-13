# Clearwing moth trajectory analysis


.colFromCategory <- function(category) {
  # Bees are 2 (which defaults to red), wasps are 1 (black)
  ifelse(grepl("bee", category), 2, 1)
}

# Returns a pch (point symbol) to use to plot a trajectory, based on category
.pchFromCategory <- function(category) {
  isMimic <- grepl("mimic", category)
  isWasp <- grepl("wasp", category)

  ifelse(isMimic & !isWasp, 15,                       # Bee mimic
         ifelse(!isMimic & !isWasp, 17,               # Bee
                ifelse(isMimic & isWasp, 22,          # Wasp mimic
                       ifelse(!isMimic & isWasp, 24,  # Wasp
                              8))))                   # Bug ;)
}

# Any slower than this speed is considered to be hovering
MAX_HOVER_SPEED <- 0.1

# Returns the longest single hovering time for a trajectory, or 0 if there is no
# hovering
longest_hover_time <- function(trj) {
  intervals <- TrajSpeedIntervals(trj, slowerThan = MAX_HOVER_SPEED)
  # Return time of longest interval, or 0 if there are no intervals
  max(c(0, intervals$duration))
}

# A function to calculate every index that we are interested in for a single trajectory.
clearwing_traj_indices <- function(trj, fractal_step_range){
  derivs <- TrajDerivatives(trj)
  #plot(derivs$speed, type = 'l', col = 'red')
  #plot(derivs$acceleration, type = 'l')
  max_hover <- longest_hover_time(trj)
  mean_speed <- mean(derivs$speed)
  min_speed <- min(derivs$speed)
  max_speed <- max(derivs$speed)
  sd_speed <- sd(derivs$speed)
  sinuosity <- TrajSinuosity(trj)
  straightness <- TrajStraightness(trj)
  fractal_dimension <- TrajFractalDimension(trj, fractal_step_range)
  resampled <- TrajRediscretize(trj, .001)
  corr <- TrajDirectionAutocorrelations(resampled, round(nrow(resampled) / 4))
  #plot(corr, type='l')
  first_min <- TrajDAFindFirstMinimum(corr)
  Emax <- TrajEmax(trj)
  directional_change_mean <- mean(TrajDirectionalChange(trj))
  directional_change_sd <- sd(TrajDirectionalChange(trj))
  list(longest_hover_time = max_hover, mean_speed = mean_speed, min_speed = min_speed, max_speed = max_speed, sd_speed = sd_speed, sinuosity = sinuosity, straightness = straightness, fractal_dimension = fractal_dimension, first_min_deltaS = first_min[1], first_min_C = first_min[2], Emax = Emax, directional_change_mean = directional_change_mean, directional_change_sd = directional_change_sd)
}

# Reads in clearwing moth trajectores as specified in the file
# "data/clearwing-moths/files_for_analyses.csv", then scales and smooths them.
# Calculates various indices, then returns a data frame with a row for each
# trajectory and various indices for each column.
# 
# Note that this function takes minutes to run because the fractal dimension
# calculation is slow.
ReadClearwingTrajectories <- function() {
  # Read in the CSV file which identifies trajectory files together with indices describing the trajectories
  file_list <- read.csv("data/clearwing-moths/files_for_analyses.csv", sep = ";", stringsAsFactors = FALSE)
  # Omit blank lines
  file_list <- na.omit(file_list)
  # Remove 1 wasp trajectory which is atypical because the wasp landed 
  file_list <- file_list[!grepl("Pyrophleps_C0073xypts.csv", file_list$csv.file.name), ]
  # Read in paths from CSV files, convert to Trajectory objects, scale and smooth them
  trjs <- TrajsBuild(file_list$csv.file.name, file_list$fps, file_list$scale, "m", rootDir = "data/clearwing-moths", smoothP = 3, smoothN = 101)

  # Calculate various indices
  # Use the same scale range for all fractal dimension calculations
  step_lengths <- TrajsStepLengths(trjs)
  mean_step_length <- mean(step_lengths)
  fractal_step_range <- TrajLogSequence(mean_step_length / 2, mean_step_length * 5, 10)
  indices <- TrajsMergeStats(trjs, clearwing_traj_indices, fractal_step_range)
  
  # Return a data frame containing both the file list (for category, species
  # etc) and the calculated indices
  cbind(file_list, indices)
}


# Plots a cluster analysis of clearwing trajectories.
#
PlotClearwingClusters <- function(indices) {
  # First columns are from the file list and contain species, category etc.
  # There are 13 calculated indices (defined in clearwing_traj_indices),
  # and they are the last columns
  nc <- ncol(indices)
  behaviour <- indices[, (nc-13+1):nc]

  # For both PCA and cluster analysis, first remove NAs  
  pca_able <- TrajsStatsReplaceNAs(behaviour, "first_min_deltaS", flagColumn = "no_first_min")
  pca_able <- TrajsStatsReplaceNAs(pca_able, "first_min_C")
  
  # Group into 2 clusters based on scaled indices (with NAs removed)
  set.seed(1)
  centres <- 2
  km <- kmeans(scale(pca_able), centres)
  #table(data.frame(indices$category, km$cluster))
  
  # PCA so we can plot in 2 dimensions
  par(mar = c(4, 4, 0, 0) + 0.1)
  pca <- prcomp(pca_able, scale. = TRUE)
  bg <- "grey"
  plot(pca$x[, c(1, 2)], col = .colFromCategory(indices$category), pch = .pchFromCategory(indices$category), bg = bg)
  # We assume that there is 1 cluster containing wasps and wasp mimics, and
  # another containing bees and bee mimics. Highlight any trajectories which
  # aren't in their expected cluster
  mismatched <- which(km$cluster != .colFromCategory(indices$category))
  points(pca$x[mismatched, c(1, 2)], pch = 1, cex = 3)
  legend("topleft", c("Bee", "Bee mimic", "Wasp", "Wasp mimic"), inset = c(0.01, 0.01),
         col = .colFromCategory(c("bee", "bee-mimic", "wasp", "wasp-mimic")),
         pt.bg = bg,
         pch = .pchFromCategory(c("bee", "bee-mimic", "wasp", "wasp-mimic"))
  )
}