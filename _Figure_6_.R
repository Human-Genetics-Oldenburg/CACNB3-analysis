#
# R script for analyzing Ca2+ imaging data from ATP stimulation, only beta 3 subunit
#
# please cite:
# Jüschke C et al. CACNB3 defects are associated with infantile idiopathic nystagmus.
# Brain communications. 2025
#
# Figure 6, panels A-D
#
#
library(stringr)
library(pracma)


today <- gsub("-", "", Sys.Date())
today <- "fig6"

# folder with the experimental data
folder <- "data/"

# subfolder/file which lists the ROIs to exclude from analysis
# based on visual inspection of 340 nm and 380 nm channels
EXCLUDE <- "filter/cavbeta.csv"


# constant values defining the measurements
MAX_TIME = 300           # duration of measurement 2*300 = 600 sek = 10 min
MAX_ROIs = 280           # max. number of ROIs to consider

UNSTIM = 1:90            # range for determining resting level
STIMUL = 115             # time of ATP stimulation (exact time 120 +/- 1)

# constant values defined for analysis
Z = 6

# color definition
colmock <- "#A9A9A96F"
colwt   <- "#FF00006F"
colmut  <- "#0000FF6F"

#
# function to read the data set into an array
#
read_data_array <- function(folder, subfolder, ratio = TRUE) {
  # files containing the Calcium imaging data
  files <- list.files(paste0(folder, subfolder),
                      pattern = "*.csv",
                      recursive = F)[c(-(1:6))] #,-(28:30))] # select/exclude specific dates
  
  # assign type of experiment
  exp_type <- rep("unknown", length(files))
  exp_type[grep("MOCK", files)] <- "MOCK"
  exp_type[grep("WT", files)] <- "WT"
  exp_type[grep("MUT", files)] <- "MUT"
  
  # assign date of experiment
  exp_date <- substring(files, 1, 8)
  
  meta_data <- data.frame(
    file      = files,
    exp_type  = exp_type,
    exp_date  = exp_date,
    roi_num   = rep(NA, length(files)),
    stim_time = rep(STIMUL, length(files)),
    stringsAsFactors = F
  )
  
  # prepare data array
  data_array <- array(
    data = NA,
    dim = c(MAX_TIME, MAX_ROIs, length(unique(exp_type)), length(unique(exp_date))),
    dimnames = list(
      time = 2 * c(0:(MAX_TIME - 1)),
      roi = paste0("ROI.", str_pad(c(1:MAX_ROIs), 4, pad = "0")),
      condition = unique(exp_type),
      experiment = unique(exp_date)
    )
  )
  
  # read experimental data into data_array
  for (fileno in 1:length(files)) {
    file <- files[fileno]
    
    ca_flow <- read.csv(
      paste0(folder, subfolder, file),
      header = T,
      row.names = 1,
      colClasses = "numeric",
      skip = 1,
      skipNul = T,
      sep = ",",
      dec = ".",
      nrows = MAX_TIME
    )
    if (ratio) {
      back <- 0
      colnames(ca_flow) <- paste0("ROI.", str_pad(c(1:ncol(ca_flow)), 4, pad = "0"))
    }
    else {
      colnames(ca_flow) <- c("background", paste0("ROI.", str_pad(c(1:(ncol(ca_flow) - 1)), 4, pad = "0")))
      back <- ca_flow$background
      ca_flow <- ca_flow[, -1]
    }
    
    meta_data[fileno, "roi_num"] <- ncol(ca_flow)
    mat = as.matrix(ca_flow)
    for (roi in (dimnames(mat)[[2]])[1:min(ncol(ca_flow), MAX_ROIs)]) {
      signal <- mat[1:min(nrow(ca_flow):MAX_TIME), roi] - back
      if (sum(abs(signal)) == 0)
        # if only zero, then set to NA
        signal <- rep(NA, min(nrow(ca_flow):MAX_TIME))
      data_array[1:min(nrow(ca_flow):MAX_TIME), roi, meta_data[fileno, "exp_type"], meta_data[fileno, "exp_date"]] = signal
    }
  }
  list(dat = data_array, met = meta_data)
}

# read all data
subfolder <- "Data for Analysis_cavbeta_340_380/"
data_array_ratio <- read_data_array(folder, subfolder)


# read list of ROIs that are manually excluded
desel <- read.csv(
  paste0(folder, EXCLUDE, ""),
  header = T,
  sep = ",",
  dec = ".",
  colClasses = "character"
)
desel$ROI <- paste0("ROI.", str_pad(desel$ROI, 4, pad = "0"))


########
########
########

meta_data <- data_array_ratio$met
data_array <- data_array_ratio$dat

#
# the command "which(data_array==0, arr.ind=T)" can be used to find the zeros in the array
# (which should not be in, but can be ignored)
# ROIs with only zeros are set to NA
#

sink(paste0(today, "_output_beta.txt"))

cat("Data Summary\n")
meta_data

Ns <- length(unique(meta_data$exp_date))
cat("\nNumber of experimental days: ", Ns, "\n")
cat("\nNumber of ROIs (in total): ", sum(meta_data$roi_num), "\n")
cat("\nNumber of ROIs (split by condition): \n")
tapply(meta_data$roi_num, meta_data$exp_type, "sum")[c("MOCK", "WT", "MUT")]


cat("\n==========\n")


###
#
# determine data characteristics of all ROIs (prior to selection)
#
cat("\n\nData analysis using all valid ROIs:\n")

# calculate statistics for pre-stimulation
unstimulated_mean <- apply(data_array[UNSTIM, , , ], c(2, 3, 4), mean, na.rm = T)
unstimulated_std  <- apply(data_array[UNSTIM, , , ], c(2, 3, 4), sd, na.rm = T)


roi_valid <- !is.na(unstimulated_mean)
cat("Number of valid ROIs: ", sum(roi_valid, na.rm = T), "\n")

cat("\nNumber of valid ROIs (split by condition): \n")
apply(roi_valid, c(2), sum, na.rm = T)



# substract pre-stimulation mean
data_unst_substr <- sweep(data_array, STATS = unstimulated_mean, MARGIN = c(2, 3, 4))



###
#
# for all ROIs (uses data after substraction of pre-stimulation mean)
#

# determine max peak height
peak_ampl <- apply(data_unst_substr[STIMUL:MAX_TIME, , , ], c(2, 3, 4), max, na.rm = T)
peak_ampl[peak_ampl <= 0] <- NA

# determine max peak time point
whichmax <- function(x, na.rm = T) {
  res <- unname(which.max(x))
  if (length(res) > 0)
    res
  else
    NA
}
peak_maxt <- apply(data_unst_substr[1:MAX_TIME, , , ], c(2, 3, 4), whichmax, na.rm = T)

# determine peak start time point (works with and without "unstimulated" substraction)
whichstart <- function(x, z = Z, na.rm = T) {
  x_mean <- mean(x[UNSTIM], trim = 0, na.rm)
  x_sd   <- sd(x[UNSTIM], na.rm)
  res <- min(which(x > (z * x_sd + x_mean)), MAX_TIME)
  if ((length(res) > 0) && (res < MAX_TIME))
    res
  else
    NA
}
peak_start <- apply(data_unst_substr[1:MAX_TIME, , , ], c(2, 3, 4), whichstart, na.rm = T)


# determine peak end time point (works with and without "unstimulated" substraction)
# with mult. peaks it gives the end of the last peak
whichend <- function(x, z = Z, na.rm = T) {
  x_mean <- mean(x[UNSTIM], trim = 0, na.rm)
  x_sd   <- sd(x[UNSTIM], na.rm)
  res <- max(which(x > (z * x_sd + x_mean)), 0)
  if ((length(res) > 0) && (res > 0))
    res
  else
    NA
}
peak_end <- apply(data_unst_substr[1:MAX_TIME, , , ], c(2, 3, 4), whichend, na.rm = T)

#
#
#
###


#####
#
# selection of appropriate ROIs based on data characteristics
#

set_NA <- function(x, y) {
  tmp <- x * y
  tmp[tmp == 0] <- NA
  return(tmp)
}

# automatic selection if peak is present:
# - peak must be higher than mean (=0, after pre-stim. substr.) plus 6x SD pre-stimulation "noise"
# - peak must be at least 20% above background level
# - ...

peak_valid <- (peak_ampl > (unstimulated_std * Z)) &
  (peak_ampl > (unstimulated_mean * 0.2)) &
  (unstimulated_std < 0.05) &
  (unstimulated_mean > 0.1) &
  (peak_start > STIMUL) & (peak_start < (STIMUL + 120)) &
  (peak_maxt > STIMUL) & (peak_maxt < STIMUL + 120)
#
#
#
#####



#####
#
# de-selection of in-appropriate ROIs based on visual inspection of 340 nm and 380 nm channels
#

for (i in 1:nrow(desel))
  if (desel[i, ]$date %in% dimnames(data_array)$experiment)
    peak_valid[desel[i, ]$ROI, desel[i, ]$condition, desel[i, ]$date] <- FALSE

#
#
#
#####



cat("\n==========\n==========\n")

cat("\n\nData analysis using ROIs with valid peaks:\n")

cat("Number of ROIs with valid peaks: ", sum(peak_valid, na.rm = T), "\n")

cat("\nNumber of ROIs with valid peaks (split by condition): \n")
(roisn <- apply(peak_valid, c(2), sum, na.rm = T))

cat("\n==========\n==========\n")

# ROIs that do not show a valid peak are set to NA
data_clean <- sweep(data_array,
                    STATS = peak_valid,
                    FUN = set_NA,
                    MARGIN = c(2, 3, 4))

#
#
#
#####


#####
#
#
#

unstimulated_mean <- apply(data_clean[UNSTIM, , , ], c(2, 3, 4), mean, na.rm = T)
unstimulated_std <- apply(data_clean[UNSTIM, , , ], c(2, 3, 4), sd, na.rm = T)

data_clean_substr <- sweep(data_clean, STATS = unstimulated_mean, MARGIN = c(2, 3, 4))


#
# resting Calcium
#
cat("\nResting Calcium\n")
cat("\nResting Calcium levels mean (by condition):\n")
(val <- apply(unstimulated_mean, c(2), mean, na.rm = T))
# SEM
err <- apply(unstimulated_mean, c(2), sd, na.rm = T) / sqrt(apply(peak_valid, c(2), sum, na.rm = T))

# Statistics
pairwise.wilcox.test(
  c(unstimulated_mean[, "MOCK", ], unstimulated_mean[, "WT", ], unstimulated_mean[, "MUT", ]),
  rep(c("MOCK", "WT", "MUT"), each = length(unstimulated_mean[, "MOCK", ]))
)

#
# Figure 6B
#
cairo_pdf(
  file = paste0(today, "_Ca_ATP_resting.pdf"),
  width = 3,
  height = 4
)
par(
  mar = c(5, 3.5, 2, 1),
  mgp = c(2.4, 0.8, 0),
  las = 1,
  family = "sans"
)
bar <- barplot(
  val,
  legend = F,
  las = 1,
  ylim = c(0.23, 0.262),
  ylab = "F340/F380",
  col = c(colmock, colwt, colmut),
  main = expression("resting Ca"^{2 * "+"})
)
arrows(
  bar,
  val + err,
  bar,
  val - err,
  lwd = 1.5,
  angle = 90,
  code = 3,
  length = 0.05
)

p <- 0.26
lines(c(bar[1, 1], bar[3, 1]), c(p, p), lwd = 2)
text((bar[1, 1] + bar[3, 1]) / 2, p + 0.001, "**", cex = 1.2)
p <- 0.255
lines(c(bar[1, 1] + 0.05, bar[2, 1] - 0.05), c(p, p), lwd = 2)
text((bar[1, 1] + bar[2, 1]) / 2, p + 0.001, "***", cex = 1.2)
lines(c(bar[2, 1] + 0.05, bar[3, 1] - 0.05), c(p, p), lwd = 2)
text((bar[2, 1] + bar[3, 1]) / 2, p + 0.001, "n.s.", cex = 0.8)

dev.off()



#
# determine peak height
#
cat("\nCalcium peak amplitude\n")
peak_ampl <- apply(data_clean_substr[STIMUL:MAX_TIME, , , ], c(2, 3, 4), max, na.rm = T)
peak_ampl[peak_ampl <= 0] <- NA

cat("\nCalcium peak amplitude mean (by condition):\n")
(val <- apply(peak_ampl, c(2), mean, na.rm = T))
# SEM
err <- apply(peak_ampl, c(2), sd, na.rm = T) / sqrt(apply(peak_valid, c(2), sum, na.rm = T))

# Statistics
pairwise.wilcox.test(c(peak_ampl[, "MOCK", ], peak_ampl[, "WT", ], peak_ampl[, "MUT", ]), 
                     rep(c("MOCK", "WT", "MUT"), each = length(peak_ampl[, "MOCK", ])))


#
# Figure 6C
#
cairo_pdf(
  file = paste0(today, "_Ca_ATP_amplitude.pdf"),
  width = 3,
  height = 4
)
par(
  mar = c(5, 3.5, 2, 1),
  mgp = c(2.4, 0.8, 0),
  las = 1,
  family = "sans"
)
bar <- barplot(
  val,
  legend = F,
  las = 1,
  ylim = c(0.3, 0.49),
  ylab = "ΔF340/F380",
  col = c(colmock, colwt, colmut),
  main = expression("peak amplitude")
)
arrows(
  bar,
  val + err,
  bar,
  val - err,
  lwd = 1.5,
  angle = 90,
  code = 3,
  length = 0.05
)

p <- 0.47
lines(c(bar[1, 1], bar[3, 1]), c(p, p), lwd = 2)
text((bar[1, 1] + bar[3, 1]) / 2, p + 0.01, "n.s.", cex = 0.8)
p <- 0.45
lines(c(bar[1, 1] + 0.05, bar[2, 1] - 0.05), c(p, p), lwd = 2)
text((bar[1, 1] + bar[2, 1]) / 2, p + 0.01, "**", cex = 1.2)
lines(c(bar[2, 1] + 0.05, bar[3, 1] - 0.05), c(p, p), lwd = 2)
text((bar[2, 1] + bar[3, 1]) / 2, p + 0.01, "***", cex = 1.2)

dev.off()


# determine peak volume
cat("\nCalcium peak volume\n")

peak_vol <- apply(data_clean_substr[STIMUL:MAX_TIME, , , ], c(2, 3, 4), sum, na.rm = T) * 2
peak_vol[peak_vol <= 0] <- NA

cat("\nCalcium peak volume mean (by condition):\n")
(val <- apply(peak_vol, c(2), mean, na.rm = T))
# SEM
err <- apply(peak_vol, c(2), sd, na.rm = T) / sqrt(apply(peak_valid, c(2), sum, na.rm = T))

# Statistics
pairwise.wilcox.test(c(peak_vol[, "MOCK", ], peak_vol[, "WT", ], peak_vol[, "MUT", ]), 
                     rep(c("MOCK", "WT", "MUT"), each = length(peak_vol[, "MOCK", ])))

#
# Figure 6D
#
cairo_pdf(
  file = paste0(today, "_Ca_ATP_volume.pdf"),
  width = 3,
  height = 4
)
par(
  mar = c(5, 3.5, 2, 1),
  mgp = c(2.4, 0.8, 0),
  las = 1,
  family = "sans"
)
bar <- barplot(
  val,
  legend = F,
  las = 1,
  ylim = c(7, 13),
  ylab = "ΔF340/F380*s",
  col = c(colmock, colwt, colmut),
  main = expression("area under the curve")
)
arrows(
  bar,
  val + err,
  bar,
  val - err,
  lwd = 1.5,
  angle = 90,
  code = 3,
  length = 0.05
)

p <- 12
lines(c(bar[1, 1], bar[3, 1]), c(p, p), lwd = 2)
text((bar[1, 1] + bar[3, 1]) / 2, p + 0.5, "n.s.", cex = 0.8)
p <- 11
lines(c(bar[1, 1] + 0.1, bar[2, 1] - 0.1), c(p, p), lwd = 2)
text((bar[1, 1] + bar[2, 1]) / 2, p + 0.5, "n.s.", cex = 0.8)
lines(c(bar[2, 1] + 0.1, bar[3, 1] - 0.1), c(p, p), lwd = 2)
text((bar[2, 1] + bar[3, 1]) / 2, p + 0.5, "n.s.", cex = 0.8)

dev.off()



# Align data with respect to peak starting time
data_aligned <- array(
  data = NA,
  dim = c(
    MAX_TIME,
    MAX_ROIs,
    length(dimnames(data_array)$condition),
    length(dimnames(data_array)$experiment)
  ),
  dimnames = list(
    time = 2 * c(0:(MAX_TIME - 1)),
    roi = paste0("ROI.", str_pad(c(1:MAX_ROIs), 4, pad = "0")),
    condition = dimnames(data_array)$condition,
    experiment = dimnames(data_array)$experiment
  )
)
data_dim <- dimnames(data_aligned)
for (roi in data_dim$roi)
  for (cond in data_dim$condition)
    for (exp in data_dim$experiment)
      if (!is.na(peak_start[roi, cond, exp]) &&
          peak_start[roi, cond, exp] > STIMUL) {
        ran <- max(1, peak_start[roi, cond, exp] - STIMUL):min(peak_start[roi, cond, exp] + 100, MAX_TIME)
        #        cat(range(ran),roi,cond,exp,"\n")
        val <- data_array[ran, roi, cond, exp]
        data_aligned[1:length(ran), roi, cond, exp] <- val
      }

#
# Figure 6A
#
cairo_pdf(
  file = paste0(today, "_Ca_ATP_curves_aln.pdf"),
  width = 5,
  height = 4
)
STATIME = 50
ENDTIME = 200
par(
  mar = c(5, 3.5, 2, 1),
  mgp = c(2.4, 0.8, 0),
  las = 1,
  family = "sans"
)
rmean <- rowMeans(data_aligned[, , "MUT", ], na.rm = T)
rsdev <- apply(data_aligned[, , "MUT", ], 1, sd, na.rm = T)
plot(
  STATIME:ENDTIME * 2,
  rmean[STATIME:ENDTIME],
  xlab = "time (s)",
  ylab = "F340/F380",
  main = "Summary (aligned)",
  las = 1,
  ylim = c(0.2, 0.7),
  type = "l",
  col = colmut,
  lwd = 3
)
rmean <- rowMeans(data_aligned[, , "WT", ], na.rm = T)
rsdev <- apply(data_aligned[, , "WT", ], 1, sd, na.rm = T)
points(
  STATIME:ENDTIME * 2,
  rmean[STATIME:ENDTIME],
  pch = ".",
  type = "l",
  col = colwt,
  lwd = 3
)
rmean <- rowMeans(data_aligned[, , "MOCK", ], na.rm = T)
rsdev <- apply(data_aligned[, , "MOCK", ], 1, sd, na.rm = T)
points(
  STATIME:ENDTIME * 2,
  rmean[STATIME:ENDTIME],
  pch = ".",
  type = "l",
  col = colmock,
  lwd = 3
)

legend(
  "topleft",
  paste0(c("mock (", "Cavβ3 wt (", "Cavβ3 mut ("), roisn[c("MOCK", "WT", "MUT")], rep(paste0("/", Ns, ")"), 3)),
  fill = c(colmock, colwt, colmut),
  bty = "n"
)

arrows(
  STIMUL * 2,
  0.65,
  ENDTIME * 2,
  0.65,
  lwd = 1.5,
  angle = 90,
  code = 3,
  length = 0.05
)
text(STIMUL * 2, 0.67, "ATP", pos = 4)

dev.off()


cat("\n==========\n==========\n\nsessionInfo()\n\n")

sessionInfo()

sink()

