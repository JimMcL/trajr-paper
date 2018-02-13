# Script to generate document figures

library(trajr)
source("R/whale-analysis.R")
source("R/clearwing-analysis.R")

# ============= General Functions ===============

# "Private"
FigureFile <- function(number, extension) file.path(".", "figures", sprintf("Figure_%d.%s", number, extension))

ReportFile <- function(name) file.path(".", "reports", name)


.JplotToDevice <- function(filename, plotFn, onlyIfDoesntExist, openDeviceFn, closeDevFn = dev.off) {
  if (!onlyIfDoesntExist || !file.exists(filename)) {
    openDeviceFn()
    tryCatch({
      if (is.function(plotFn))
        plotFn()
      else
        invisible(eval(plotFn))
    }, finally = {
      closeDevFn()
    })
  }
}

# Default resolution seems to be 72, so to increase image pixel size without decreasing text size, 
# line width etc, increase resolution accordingly.
# Doesn't work with ggplot since it must be evaluated to work.
# Try using ggsave or print(plotFn())
PlotToPng <- function(filename, plotFn, width=800, aspectRatio = 3 / 2, type = 'cairo', onlyIfDoesntExist = F, ...) {
  .JplotToDevice(filename, plotFn, onlyIfDoesntExist, function () {
    # type = 'cairo' seems to produce _much_ nicer graphics with antialiasing
    png(filename, width = width, height = width / aspectRatio, type = type, ...)
  })
}

# Beware: width and height are in inches, not pixels!
PlotToPDF <- function(filename, plotFn, width=8, aspectRatio = 3 / 2, onlyIfDoesntExist = F) {
  .JplotToDevice(filename, plotFn, onlyIfDoesntExist, function () {
    pdf(filename, width = width, height = width / aspectRatio, paper = "special")
  })
}

PlotToEPS <- function(filename, plotFn, width=8, aspectRatio = 3 / 2, onlyIfDoesntExist = F, ...) {
  .JplotToDevice(filename, plotFn, onlyIfDoesntExist, function () {
    setEPS()
    postscript(filename, width = width, height = width / aspectRatio, paper = "special", ...)
  })
}

ReportToFile <- function(filename, expr) {
  oldOptions <- options()
  on.exit(options(oldOptions))
  options(width = 10000)
  
  filename <- ReportFile(filename)
  .JplotToDevice(filename, expr, onlyIfDoesntExist = FALSE, 
                 openDeviceFn = function () { sink(filename) },
                 closeDevFn = sink)
}


doLabel <- function(label) {
  if (!missing(label))
    mtext(label, line = -1.6, adj = 0.92)
}

plotTrajectory <- function(filename, scale, label){
  ## Daerlac nigricans
  trj <- TrajFromCoords(read.csv(filename, stringsAsFactors = FALSE), 
                        xCol = "x", yCol = "y", timeCol = "Time")
  trj <- TrajScale(trj, scale = scale, units = "m")
  
  plot(trj)
  doLabel(label)
}

# ============= Figure-specific functions ===============

sampleTrajectories <- function() {
  par(mfrow = c(1, 3))
  par(mar = c(4, 4, 0, 0) + 0.1)
  # Polyrhachis sp.
  #plotTrajectory("data/ant-mimics/3543.csv", scale = .220 / 720)
  # Crematogaster sp.
  plotTrajectory("data/ant-mimics/3548.csv", scale = .220 / 720, "a)")
  ## Daerlac nigricans
  plotTrajectory("data/ant-mimics/3530.csv", scale = .220 / 720, "b)")
  # Zodariidae sp.
  plotTrajectory("data/ant-mimics/3527.csv", scale = .220 / 720, "c)")
}

scaleVaryingFractalIndices <- function() {
  .pf <- function(trj, min = 1, max = 100, adjustD = TRUE, ...) {
    v <- TrajFractalDimensionValues(trj, TrajLogSequence(min, max, 20), adjustD = adjustD)
    plot(v, pch = 16, log = "xy", ...)
    y <- log(v[,2])
    x <- log(v[,1])
    abline(lm(v[, 2] ~ v[, 1]), untf = TRUE)
  }
  par(mfrow = c(3, 1))
  crw <- TrajGenerate(10)
  .pf(crw, max = 5)
  fw <- TrajGenerate(linearErrorDist = function(n) runif(n, .1, 100))
  .pf(fw)
  bw <- TrajGenerate(angularErrorDist = function(n) runif(n, -pi, pi), linearErrorDist = function(n) runif(n, .1, 100))
  .pf(bw)
}

randomTrajs <- function() {
  set.seed(41)
  
  par(mfrow = c(1, 3), mar = c(4, 4, 1, 0) + .1)
  trj <- TrajGenerate(n = 100)
  plot(trj)
  doLabel("a)")
  
  trj <- TrajGenerate(n = 50, random = FALSE)
  plot(trj)
  doLabel("b)")
  
  set.seed(2)
  trj <- TrajGenerate(n = 50, linearErrorDist = stats::rcauchy, angularErrorDist = function(n) runif(n, -pi, pi))
  plot(trj)
  doLabel("c)")
}

# __________ End functions ______________


#### Generate figures and reports ####

ar <- 2
PlotToEPS(FigureFile(1, "eps"), sampleTrajectories, width = 6, aspectRatio = ar, bg = "white")
# Word can't embed EPS anymore (security problem), so create PNG for temporarily adding to Word
PlotToPng(FigureFile(1, "png"), sampleTrajectories, width = 1800, aspectRatio = ar, res = 300)

ar = 2.5
PlotToEPS(FigureFile(2, "eps"), randomTrajs, width = 6, aspectRatio = ar, bg = "white")
PlotToPng(FigureFile(2, "png"), randomTrajs, width = 1800, aspectRatio = ar, res = 300)

# Whale analysis example
ar = 3 / 2
whaleTrjs <- ReadWhaleTrajectories()
PlotToEPS(FigureFile(3, "eps"), { PlotWhaleTrajectories(whaleTrjs) }, width = 6, aspectRatio = ar, bg = "white")
PlotToPng(FigureFile(3, "png"), { PlotWhaleTrajectories(whaleTrjs) }, width = 1800, aspectRatio = ar, res = 300)

ReportToFile("whale-all-indices.txt", { ReportAllWhaleParams(whaleTrjs) })
ReportToFile("whale-report.txt", { ReportWhaleStats(whaleTrjs) })

# Clearwing moth example
ar = 3 / 2
clearwingParams <- ReadClearwingTrajectories()    # VERY slow to run
PlotToEPS(FigureFile(4, "eps"), { PlotClearwingClusters(clearwingParams) }, width = 6, aspectRatio = ar, bg = "white")
PlotToPng(FigureFile(4, "png"), { PlotClearwingClusters(clearwingParams) }, width = 1800, aspectRatio = ar, res = 300)
