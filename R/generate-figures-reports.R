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
  trj <- TrajFromCoords(read.csv(filename, stringsAsFactors = FALSE), 
                        xCol = "x", yCol = "y", timeCol = "Time")
  trj <- TrajScale(trj, scale = scale, units = "m")
  
  plot(trj, lwd = 2)
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
  # Daerlac nigricans
  plotTrajectory("data/ant-mimics/3530.csv", scale = .220 / 720, "b)")
  # Zodariidae sp.
  plotTrajectory("data/ant-mimics/3527.csv", scale = .220 / 720, "c)")
}


randomTrajs <- function() {
  set.seed(41)
  
  lwd = 2
  par(mfrow = c(1, 3), mar = c(4, 4, 1, 0) + .1)
  trj <- TrajGenerate(n = 100)
  plot(trj, lwd = lwd)
  doLabel("a)")
  
  trj <- TrajGenerate(n = 50, random = FALSE)
  plot(trj, lwd = lwd)
  doLabel("b)")
  
  set.seed(2)
  trj <- TrajGenerate(n = 50, linearErrorDist = stats::rcauchy, angularErrorDist = function(n) runif(n, -pi, pi))
  plot(trj, lwd = lwd)
  doLabel("c)")
}

# __________ End functions ______________


#### Generate figures and reports ####

ar <- 2
PlotToEPS(FigureFile(1, "eps"), sampleTrajectories, width = 18, aspectRatio = ar, pointsize = 36, bg = "white")
# Word can't embed EPS anymore (security problem), so create PNG for temporarily adding to Word
PlotToPng(FigureFile(1, "png"), sampleTrajectories, width = 1800, aspectRatio = ar, res = 300)

ar <- 2.5
PlotToEPS(FigureFile(2, "eps"), randomTrajs, width = 18, aspectRatio = ar, pointsize = 36, bg = "white")
PlotToPng(FigureFile(2, "png"), randomTrajs, width = 1800, aspectRatio = ar, res = 300)

# Whale analysis example
ar <- 3 / 2
whaleTrjs <- ReadWhaleTrajectories()
PlotToEPS(FigureFile(3, "eps"), { PlotWhaleTrajectories(whaleTrjs, 24, 1.1) }, width = 18, aspectRatio = ar, bg = "white")
PlotToPng(FigureFile(3, "png"), { PlotWhaleTrajectories(whaleTrjs) }, width = 1800, aspectRatio = ar, res = 300)

ReportToFile("whale-all-indices.txt", { ReportAllWhaleParams(whaleTrjs) })
ReportToFile("whale-all-indices.csv", { ReportAllWhaleParams(whaleTrjs, asCSV = TRUE) })
ReportToFile("whale-report.txt", { ReportWhaleStats(whaleTrjs) })

# Clearwing moth example
ar <- 3 / 2
clearwingParams <- ReadClearwingTrajectories()    # VERY slow to run
PlotToEPS(FigureFile(4, "eps"), { PlotClearwingPCA(clearwingParams, lwd = 2) }, width = 18, aspectRatio = ar, pointsize = 32, bg = "white")
PlotToPng(FigureFile(4, "png"), { PlotClearwingPCA(clearwingParams) }, width = 1800, aspectRatio = ar, res = 300)
