# Clearwing moth trajectory analysis


.colFromCategory <- function(category) {
  ifelse(grepl("bee", category), 2, 1)
}

.pchFromCategory <- function(category, hollow = FALSE) {
  if (hollow)
    ifelse(grepl("mimic", category), 0, 2)
  else
    ifelse(grepl("mimic", category), 15, 17)
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

# A function to calculate every parameter that we are interested in for a single trajectory.
clearwing_traj_parameters <- function(trj, fractal_step_range){
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
# Calculates various parameters, then returns a data frame with a row for each
# trajectory and various parameters for each column.
# 
# Note that this function takes minutes to run because the fractal dimension
# calculation is slow.
ReadClearwingTrajectories <- function() {
  # Read in the CSV file which identifies trajectory files together with parameters describing the trajectories
  file_list <- read.csv("data/clearwing-moths/files_for_analyses.csv", sep = ";", stringsAsFactors = FALSE)
  # Omit blank lines
  file_list <- na.omit(file_list)
  # Remove 1 wasp trajectory which is atypical because the wasp landed 
  file_list <- file_list[!grepl("Pyrophleps_C0073xypts.csv", file_list$csv.file.name), ]
  # Read in paths from CSV (actually tab-separated) files, convert to Trajectory objects, scale and smooth them
  trjs <- TrajsBuild(file_list$csv.file.name, file_list$fps, file_list$scale, "m", rootDir = "data/clearwing-moths", smoothP = 3, smoothN = 101)

  # Calculate various parameters
  # Use the same scale range for all fractal dimension calculations
  step_lengths <- TrajsStepLengths(trjs)
  mean_step_length <- mean(step_lengths)
  fractal_step_range <- TrajLogSequence(mean_step_length / 2, mean_step_length * 5, 10)
  params <- TrajsMergeStats(trjs, clearwing_traj_parameters, fractal_step_range)
  
  # Return a data frame containing both the file list (for category, species
  # etc) and the calculated parameters
  cbind(file_list, params)
}


# Plots a cluster analysis of clearwing trajectories.
#
PlotClearwingClusters <- function(params) {
  # First columns are from the file list and contain species, category etc.
  # There are 13 calculated parameters (defined in clearwing_traj_parameters),
  # and they are the last columns
  nc <- ncol(params)
  behaviour <- params[, (nc-13+1):nc]

  # For both PCA and cluster analysis, first remove NAs  
  pca_able <- TrajsStatsReplaceNAs(behaviour, "first_min_deltaS", flagColumn = "no_first_min")
  pca_able <- TrajsStatsReplaceNAs(pca_able, "first_min_C")
  
  # Group into 2 clusters based on scaled parameters (with NAs removed)
  set.seed(1)
  centres <- 2
  km <- kmeans(scale(pca_able), centres)
  #table(data.frame(file_list$category, km$cluster))
  
  # PCA so we can plot in 2 dimensions
  par(mar = c(4, 4, 0, 0) + 0.1)
  pca <- prcomp(pca_able, scale. = TRUE)
  plot(pca$x[, c(1, 2)], col = .colFromCategory(file_list$category), pch = .pchFromCategory(file_list$category))
  # We assume that there is 1 cluster containing wasps and wasp mimics, and
  # another containing bees and bee mimics. Highlight any trajectories which
  # aren't in their expected cluster
  mismatched <- km$cluster != .colFromCategory(file_list$category)
  points(pca$x[mismatched, c(1, 2)], col = km$cluster[mismatched], pch = .pchFromCategory(file_list$category[mismatched], hollow = TRUE), cex = 3)
  legend("topleft", c("Bee", "Bee mimic", "Wasp", "Wasp mimic"), inset = c(0.01, 0.01),
         col = .colFromCategory(c("bee", "bee-mimic", "wasp", "wasp-mimic")),
         pch = .pchFromCategory(c("bee", "bee-mimic", "wasp", "wasp-mimic"))
  )
}