#
# R script for analyzing Ca2+ imaging data from full channel (3 subunits)
#
# please cite:
# Jüschke C et al. CACNB3 defects are associated with infantile idiopathic nystagmus.
# Brain communications. 2025
#
# Figure 4, panels A-D
#
#


library(lmerTest)
library(emmeans)
library(reshape2)
library(stringr)
library(pracma)


#today <- gsub("-", "", Sys.Date())
today <- "fig4"

# folder with the experimental data
folder <- "data/"

# subfolder/file which lists the ROIs to exclude from analysis
# based on visual inspection of 340 nm and 380 nm channels
EXCLUDE <- "filter/fullchannel.csv"


# constant values defining the measurments
MAX_TIME = 550           # duration of measurement 2*300 = 600 sek = 10 min
MAX_ROIs = 225           # max. number of ROIs to consider

UNSTIM  = 1:90           # range for determining resting level
STIMUL  = 115            # time of potassium chloride (KCl, 55 mM) stimulation
STISTOP = 270            # end of KCl stimulation
KCL_END = 400            # time when KCl response is off (empirical)

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
  # files containing the Calcium imaging data, i.e. ration F340/F380
  files <- list.files(paste0(folder, subfolder),
                      pattern = "*.csv",
                      recursive = F)
  
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
  list(dat = data_array[, , c("MOCK", "WT", "MUT"), ], met = meta_data)
}

# read all data
subfolder <- "Data for Analysis_fullchannel_340_380/"
data_array_ratio <- read_data_array(folder, subfolder)


# read list of ROIs that are manually excluded based on visual inspection of 340 nm and 380 nm channels
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

sink(paste0(today, "_output_full.txt"))

cat("Data Summary\n")
meta_data

Ns <- length(unique(meta_data$exp_date))
cat("\nNumber of experimental days: ", Ns, "\n")
cat("\nNumber of ROIs (in total): ", sum(meta_data$roi_num), "\n")
cat("\nNumber of ROIs (split by condition): \n")
tapply(meta_data$roi_num, meta_data$exp_type, "sum")[c("MOCK","WT","MUT")]


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
  (peak_start > STIMUL) &
  (peak_maxt > STIMUL) &
  (peak_end > 420)             # ATP peak (confirms cell responsiveness at the end of the experiment)
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
# analysis should now only include the KCl peak but not the ATP peak:
# i.e. STIMUL:KCL_END
# (the ATP peak could be analysed separately)
#
#

unstimulated_mean <- apply(data_clean[UNSTIM, , , ], c(2, 3, 4), mean, na.rm = T)
unstimulated_std <- apply(data_clean[UNSTIM, , , ], c(2, 3, 4), sd, na.rm = T)

data_clean_substr <- sweep(data_clean, STATS = unstimulated_mean, MARGIN = c(2, 3, 4))


#
# resting Calcium
#
cat("\nResting Calcium\n")

long_data <- melt(unstimulated_mean)
long_data$experiment <- as.factor(long_data$experiment)
model <- lmer(value ~ condition + (1| experiment), data = long_data)
(mmres <- summary(emmeans(model, pairwise ~ condition)))
val <- mmres$emmeans$emmean
names(val) <- mmres$emmeans$condition
val_SE <- mmres$emmeans$SE
names(val_SE) <- mmres$emmeans$condition
signif <- mmres$contrasts$p.value
names(signif) <- mmres$contrasts$contrast
signif

#
# Figure 4B
#
cairo_pdf(
  file = paste0(today, "_Ca_KCl_resting_nMM.pdf"),
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
  ylim = c(0.25, 0.35),
  ylab = "F340/F380",
  col = c(colmock, colwt, colmut),
  main = expression("resting Ca"^{2 * "+"})
)
arrows(
  bar,
  val + val_SE,
  bar,
  val - val_SE,
  lwd = 1.5,
  angle = 90,
  code = 3,
  length = 0.05
)

p <- 0.34
lines(c(bar[1, 1], bar[3, 1]), c(p, p), lwd = 2)
text((bar[1, 1] + bar[3, 1]) / 2, p + 0.005, "n.s.", cex = 0.8)
p <- 0.33
lines(c(bar[1, 1] + 0.1, bar[2, 1] - 0.1), c(p, p), lwd = 2)
text((bar[1, 1] + bar[2, 1]) / 2, p + 0.005, "***", cex = 1.2)
lines(c(bar[2, 1] + 0.1, bar[3, 1] - 0.1), c(p, p), lwd = 2)
text((bar[2, 1] + bar[3, 1]) / 2, p + 0.005, "***", cex = 1.2)

dev.off()



#
# peak height
#
cat("\nCalcium peak amplitude\n")
peak_ampl <- apply(data_clean_substr[STIMUL:KCL_END, , , ], c(2, 3, 4), max, na.rm = T)
peak_ampl[peak_ampl <= 0] <- NA

long_data <- melt(peak_ampl)
long_data$experiment <- as.factor(long_data$experiment)
model <- lmer(value ~ condition + (1| experiment), data = long_data)
(mmres <- summary(emmeans(model, pairwise ~ condition)))
val <- mmres$emmeans$emmean
names(val) <- mmres$emmeans$condition
val_SE <- mmres$emmeans$SE
names(val_SE) <- mmres$emmeans$condition
signif <- mmres$contrasts$p.value
names(signif) <- mmres$contrasts$contrast
signif

#
# Figure 4C
#
cairo_pdf(
  file = paste0(today, "_Ca_KCl_amplitude_nMM.pdf"),
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
  ylim = c(0, 0.24),
  ylab = "ΔF340/F380",
  col = c(colmock, colwt, colmut),
  main = expression("peak amplitude")
)
arrows(
  bar,
  val + val_SE,
  bar,
  val - val_SE,
  lwd = 1.5,
  angle = 90,
  code = 3,
  length = 0.05
)

p <- 0.21
lines(c(bar[1, 1], bar[3, 1]), c(p, p), lwd = 2)
text((bar[1, 1] + bar[3, 1]) / 2, p + 0.01, "***", cex = 1.2)
p <- 0.18
lines(c(bar[1, 1] + 0.1, bar[2, 1] - 0.1), c(p, p), lwd = 2)
text((bar[1, 1] + bar[2, 1]) / 2, p + 0.01, "***", cex = 1.2)
lines(c(bar[2, 1] + 0.1, bar[3, 1] - 0.1), c(p, p), lwd = 2)
text((bar[2, 1] + bar[3, 1]) / 2, p + 0.01, "***", cex = 1.2)

dev.off()


# peak volume
cat("\nCalcium peak volume\n")
peak_vol <- apply(data_clean_substr[STIMUL:KCL_END, , , ], c(2, 3, 4), sum, na.rm = T) * 2
peak_vol[peak_vol <= 0] <- NA

long_data <- melt(peak_vol)
long_data$experiment <- as.factor(long_data$experiment)
model <- lmer(value ~ condition + (1| experiment), data = long_data)
(mmres <- summary(emmeans(model, pairwise ~ condition)))
val <- mmres$emmeans$emmean
names(val) <- mmres$emmeans$condition
val_SE <- mmres$emmeans$SE
names(val_SE) <- mmres$emmeans$condition
signif <- mmres$contrasts$p.value
names(signif) <- mmres$contrasts$contrast
signif

#
# Figure 4D
#
cairo_pdf(
  file = paste0(today, "_Ca_KCl_volume_nMM.pdf"),
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
  ylim = c(0, 40),
  ylab = "ΔF340/F380*s",
  col = c(colmock, colwt, colmut),
  main = expression("area under the curve")
)
arrows(
  bar,
  val + val_SE,
  bar,
  val - val_SE,
  lwd = 1.5,
  angle = 90,
  code = 3,
  length = 0.05
)

p <- 38
lines(c(bar[1, 1], bar[3, 1]), c(p, p), lwd = 2)
text((bar[1, 1] + bar[3, 1]) / 2, p + 1, "*", cex = 1.2)
p <- 35
lines(c(bar[1, 1] + 0.1, bar[2, 1] - 0.1), c(p, p), lwd = 2)
text((bar[1, 1] + bar[2, 1]) / 2, p + 1, "***", cex = 1.2)
lines(c(bar[2, 1] + 0.1, bar[3, 1] - 0.1), c(p, p), lwd = 2)
text((bar[2, 1] + bar[3, 1]) / 2, p + 1, "***", cex = 1.2)

dev.off()


# minimal peak volume/AUC
llimit <- 5
p_na <- which(peak_vol < llimit, arr.ind=T)
for (i in 1:nrow(p_na))
  data_clean[, p_na[i,"roi"], p_na[i,"condition"], p_na[i,"experiment"]] <- NA

#
# Figure 4A
#
# region for plotting
STATIME = 1
ENDTIME = 420
cairo_pdf(
  file = paste0(today, "_Ca_KCl_curves_n.pdf"),
  width = 5,
  height = 4
)
par(
  mar = c(5, 3.5, 2, 1),
  mgp = c(2.4, 0.8, 0),
  las = 1,
  family = "sans"
)
p_wt_n <- apply(data_clean[,,"WT",], c(1,3), mean, na.rm=T)
plot(
  STATIME:ENDTIME * 2,
  rowMeans(p_wt_n, na.rm=T)[STATIME:ENDTIME],
  xlab = "time (s)",
  ylab = "F340/F380",
  main = paste0("Summary (AUC > ", llimit, ")"),
  las = 1,
  ylim = c(0.28, 0.50),
  type = "l",
  col = colwt,
  lwd = 3
)
# sem_wt_n <- apply(p_wt_n, 1, sd)/sqrt(7)
# polygon(c(STATIME:ENDTIME * 2, rev(STATIME:ENDTIME * 2)),
#         c(rowMeans(p_wt_n, na.rm=T)[STATIME:ENDTIME] + sem_wt_n[STATIME:ENDTIME],
#           rev(rowMeans(p_wt_n, na.rm=T)[STATIME:ENDTIME] - sem_wt_n[STATIME:ENDTIME])),
#         col = "#FF00002F", border=NA)
p_mut_n <- apply(data_clean[,,"MUT",], c(1,3), mean, na.rm=T)
points(
  STATIME:ENDTIME * 2,
  rowMeans(p_mut_n, na.rm=T)[STATIME:ENDTIME],
  pch = ".",
  type = "l",
  col = colmut,
  lwd = 3
)
# sem_mut_n <- apply(p_mut_n, 1, sd)/sqrt(7)
# polygon(c(STATIME:ENDTIME * 2, rev(STATIME:ENDTIME * 2)),
#         c(rowMeans(p_mut_n, na.rm=T)[STATIME:ENDTIME] + sem_mut_n[STATIME:ENDTIME],
#           rev(rowMeans(p_mut_n, na.rm=T)[STATIME:ENDTIME] - sem_mut_n[STATIME:ENDTIME])),
#         col = "#0000FF2F", border=NA)
p_mock_n <- apply(data_clean[,,"MOCK",], c(1,3), mean, na.rm=T)
points(
  STATIME:ENDTIME * 2,
  rowMeans(p_mock_n, na.rm=T)[STATIME:ENDTIME],
  pch = ".",
  type = "l",
  col = colmock,
  lwd = 3
)
# sem_mock_n <- apply(p_mock_n, 1, sd)/sqrt(7)
# polygon(c(STATIME:ENDTIME * 2, rev(STATIME:ENDTIME * 2)),
#         c(rowMeans(p_mock_n, na.rm=T)[STATIME:ENDTIME] + sem_mock_n[STATIME:ENDTIME],
#           rev(rowMeans(p_mock_n, na.rm=T)[STATIME:ENDTIME] - sem_mock_n[STATIME:ENDTIME])),
#         col = "#A9A9A92F", border=NA)
legend(
  "topleft",
  paste0(
    c(
      "Cav2.2 + Cavα2δ1 + mock (",
      "Cav2.2 + Cavα2δ1 + Cavβ3 wt (",
      "Cav2.2 + Cavα2δ1 + Cavβ3 mut ("
    ),
    c(length(which(peak_vol[, "MOCK", ] > llimit)), 
      length(which(peak_vol[, "WT", ] > llimit)), 
      length(which(peak_vol[, "MUT", ] > llimit))),
    rep(paste0("/", Ns, ")"), 3)
  ),
  fill = c(colmock, colwt, colmut),
  bty = "n"
)

arrows(
  STIMUL * 2,
  0.42,
  STISTOP * 2,
  0.42,
  lwd = 1.5,
  angle = 90,
  code = 3,
  length = 0.05
)
text(STIMUL * 2, 0.43, "KCl", pos = 4)

dev.off()

cat("\n==========\n==========\n\nsessionInfo()\n\n")

sessionInfo()

sink()

