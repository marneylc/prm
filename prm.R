#!/usr/bin Rscript

# require packages
require(xcms)
require(reshape2)
require(ggplot2)
require(magrittr)
require(data.table)

# define targets in a csv file with columns: targetRT, rtwin, targetMZ, mzwin
# targetRT: retention time of the target compound in seconds
# rtwin: retention time window to search for the target compound in seconds
# targetMZ: m/z of the target compound
# mzwin: m/z window to search for the target compound
# analytes <- read.csv("analytes.csv", header = TRUE, sep = ",", row.names = 1)
main <- function() {
  analytes <- read.csv("analytes.csv", header = TRUE, sep = ",", row.names = 1)
  files <- list.files('.', full.names = TRUE, pattern = '.mzML', recursive = FALSE)
  source("hrms.R")
  peakareas <- PRM(files, analytes)
  write.csv(peakareas, file = "peakareas.csv")
}

# Get the peak areas for all compounds in all files
PRM <- function(files,analytes) {
  for (f in 1:length(files)) {
    print(files[f])
    raw_data <- readMSData(files[f], msLevel = 2, mode = "onDisk") # load data file
    analytes[,files[f]] <- NA
    
    ##############################################################################
    ## running prm without defined mzs
    # in order to determine the prms to use for each compound
    # uses get_spec function and saves the prm_traces table, figures, and spectra
    # usually only needed to run on one file
    ##############################################################################
    # for (i in 1:dim(analytes)[1]) {
    #   sumspec <- get_spec(raw_data,
    #     analyte = rownames(analytes)[i],
    #     target_rt = as.numeric(analytes$targetRT[i]),
    #     rt_win = as.numeric(analytes$rtwin[i]), 
    #     target_mz = as.numeric(analytes$targetMZ[i]), 
    #     mz_win = as.numeric(analytes$mzwin[i]))
    #   peakarea <- sum(sumspec$intensity)
    #   analytes[i,files[f]] <- peakarea
    # }
    # 
    ##############################################################################
    ## running prm with defined mzs
    for (i in 1:dim(analytes)[1]) {
      sumspec <- get_spec_definedmzs(raw_data,
        analyte = rownames(analytes)[i],
        target_rt = as.numeric(analytes$targetRT[i]),
        rt_win = as.numeric(analytes$rtwin[i]), 
        target_mz = as.numeric(analytes$targetMZ[i]), 
        mz_win = as.numeric(analytes$mzwin[i]))
      peakarea <- sum(sumspec$intensity)
      analytes[i,files[f]] <- peakarea
    }
  }
  return(analytes)
}

## For get_spec testing
# a <- 1
# analyte <- rownames(analytes)[a]
# target_rt = as.numeric(analytes$targetRT[a])
# rt_win = as.numeric(analytes$rtwin[a])
# target_mz = as.numeric(analytes$targetMZ[a])
# mz_win = as.numeric(analytes$mzwin[a])

# mzs are defined in msp files and the top 10 are selected for prm
get_spec_definedmzs <- function(raw_data, analyte, target_rt, rt_win, target_mz, mz_win) {
  if (!dir.exists("prm_deviations")) {
    dir.create("prm_deviations")
  }
  print(paste0("Collecting data for: ", analyte, " in ", row.names(raw_data@phenoData)))
  rt_range <- c(target_rt - rt_win, target_rt + rt_win)
  ms2_scans <- rtime(raw_data)[which(rtime(raw_data) > rt_range[1] & rtime(raw_data) < rt_range[2])]
  ms2_scans_index <- names(ms2_scans)
  mz_range <- c(target_mz - mz_win, target_mz + mz_win)
  specs <- list()
  spec_index <- list()
  for (i in 1:length(ms2_scans_index)) {
    sp <- raw_data[[ms2_scans_index[i]]]
    if (sp@precursorMz > mz_range[1] & sp@precursorMz < mz_range[2]) {
      spec_index <- append(spec_index, ms2_scans_index[i])
      specs <- append(specs, sp)
    }
  }
  if (length(specs) == 0) {
    return(data.frame(mz = 0, intensity = 0))
  }
  sumspec <- meanMzInts(specs, mzd = 0.01, intensityFun = max)

  x <- sumspec@mz
  y <- sumspec@intensity
  if (length(x) == length(y)) {
    sumspec_4plot <- data.frame(mz = x, intensity = y)
  } else {
    sumspec_4plot <- data.frame(mz = x, intensity = y[match(x, y)])
  }

  ## plotting tools 
  # g <- ggplot(sumspec_4plot, aes(x = mz, y = intensity)) + geom_line(
  #   color = "black", 
  #   size = 1, 
  #   alpha = 0.8) + 
  #   theme_classic() + 
  #   theme(panel.grid.major = element_blank(), 
  #     panel.grid.minor = element_blank(), 
  #     panel.border = element_blank(), 
  #     panel.background = element_blank()) + 
  #   labs(x = "m/z", y = "Intensity") + 
  #   ggtitle(analyte) + 
  #   theme(plot.title = element_text(hjust = 0.5)) +
  #   scale_x_continuous(limits = c(round(min(x), digits=0), target_mz), breaks = seq(round(min(x), digits=0), target_mz, 100)) 

  # # for testing
  # ggsave("testing.png", plot = g, width = 4, height = 4)

  # define the mz values to use when running prm for all samples
  ref_spec <- read.table(paste0("spectra/", analyte, ".msp"), header=TRUE)
  ref_spec <- data.table(ref_spec)
  setkey(ref_spec, intensity)
  # get the top 10 intensities from reference spectra
  top10_ref <- ref_spec[order(intensity, decreasing = TRUE),][1:10,]
  prmmzs <- round(top10_ref$mz, digits = 4)

  # from sumspec_4plot and prmmzs, get the signal for each mz value that is closes to prmmzs
  # uses hrms.R tool set
  targets <- data.frame(name=prmmzs,mz=prmmzs) # names could be SMILES key in the future
  targets <- data.table(targets)
  spectra <- data.table(sumspec_4plot)
  colnames(spectra) <- c("mz", "V1")
  hwidth <- 0.01
  warnings <- NULL
  target_signals <- peaktable(targets, spectra)

  # save the target_signals table under directory ./prm_deviations
  write.csv(target_signals, file = paste0("prm_deviations/",analyte,"_", row.names(raw_data@phenoData),"_prm_deviations.csv"))

  # return the final peak areas for each prm in a spec variable
  spec <- data.frame(mz = target_signals$targets.name, intensity = target_signals$signal)
  return(spec) # return spectrum with only prm integration values
}

# Define a function to get the peak area for a given compound
get_spec <- function(raw_data, analyte, target_rt, rt_win, target_mz, mz_win) {
  if (!dir.exists("plots")) {
    dir.create("plots")
  }
  if (!dir.exists("spectra")) {
    dir.create("spectra")
  }
  if (!dir.exists("traces")) {
    dir.create("traces")
  }
  print(paste0("Collecting data for: ", analyte))
  rt_range <- c(target_rt - rt_win, target_rt + rt_win)
  ms2_scans <- rtime(raw_data)[which(rtime(raw_data) > rt_range[1] & rtime(raw_data) < rt_range[2])]
  ms2_scans_index <- names(ms2_scans)
  mz_range <- c(target_mz - mz_win, target_mz + mz_win)
  specs <- list()
  spec_index <- list()
  for (i in 1:length(ms2_scans_index)) {
    sp <- raw_data[[ms2_scans_index[i]]]
    if (sp@precursorMz > mz_range[1] & sp@precursorMz < mz_range[2]) {
      spec_index <- append(spec_index, ms2_scans_index[i])
      specs <- append(specs, sp)
    }
  }
  if (length(specs) == 0) {
    return(data.frame(mz = 0, intensity = 0))
  }
  sumspec <- meanMzInts(specs, mzd = 0.01, intensityFun = max)

  x <- sumspec@mz
  y <- sumspec@intensity
  if (length(x) == length(y)) {
    sumspec_4plot <- data.frame(mz = x, intensity = y)
  } else {
    sumspec_4plot <- data.frame(mz = x, intensity = y[match(x, y)])
  }

  ## plotting tools 
  g <- ggplot(sumspec_4plot, aes(x = mz, y = intensity)) + geom_line(
    color = "black", 
    size = 1, 
    alpha = 0.8) + 
    theme_classic() + 
    theme(panel.grid.major = element_blank(), 
      panel.grid.minor = element_blank(), 
      panel.border = element_blank(), 
      panel.background = element_blank()) + 
    labs(x = "m/z", y = "Intensity") + 
    ggtitle(analyte) + 
    theme(plot.title = element_text(hjust = 0.5)) +
    scale_x_continuous(limits = c(round(min(x), digits=0), target_mz), breaks = seq(round(min(x), digits=0), target_mz, 100)) 

  # save the plot to a file with the same name as the csv file
  ggsave(paste0("plots/", analyte,"_MS2.png"), plot = g, width = 4, height = 4)

  # save the sumspec table under directory /spectra with the same name as the plot with the extension .msp
  write.table(sumspec, file = paste0("spectra/", analyte, ".msp"), sep = "\t", row.names = FALSE, col.names = c("mz", "intensity"))

  # pick top ten intensities automatically for mz values from sumspec
  # this is used when running to choose the ms2 values to use for prm
  top10 <- sumspec_4plot[order(sumspec_4plot$intensity, decreasing = TRUE),][1:10,]
  top10$mz <- round(top10$mz, digits = 4)


  # from specs variable, get the rt, mz, and signal for the top ten mz values
  trace_rts <- list()
  for (i in 1:length(specs)) {
    trace_rts[[i]] <- specs[[i]]@rt
  }

  # prepare spectra from each scan
  trace_specs <- list()
  for (i in 1:length(specs)) {
    trace_specs[[i]] <- data.frame(specs[[i]]@mz, specs[[i]]@intensity)
    colnames(trace_specs[[i]]) <- c("mz", "intensity")
    trace_specs[[i]]$mz <- round(trace_specs[[i]]$mz, digits = 4)
    trace_specs[[i]] <- aggregate(trace_specs[[i]]$intensity, by = list(trace_specs[[i]]$mz), FUN = sum)
    colnames(trace_specs[[i]]) <- c("mz", "intensity")
  }

  trace_target_signals <- list()
  for (i in 1:length(trace_specs)) {
    spectra <- data.table(trace_specs[[i]])
    targets <- data.frame(name=top10$mz,top10) # names could be SMILES key in the future
    targets <- data.table(targets)
    colnames(spectra) <- c("mz", "V1")
    hwidth <- 0.01
    warnings <- NULL
    trace_target_signals[[i]] <- peaktable(targets, spectra)
  }

  # concatenate each trace_target_signals into a data from of columns targets.name, rt, signal
  prm_traces <- data.frame(targets$name)
  for (i in 1:length(trace_target_signals)) {
    prm_traces[,i+1] <- data.frame(trace_target_signals[[i]]$signal)
  }

  colnames(prm_traces) <- c("prms",seq(1,length(trace_rts),1))
  prm_traces <- t(prm_traces)
  colnames(prm_traces) <- prm_traces[1,]
  prm_traces <- prm_traces[2:dim(prm_traces)[1],]
  prm_traces <- data.frame(as.numeric(trace_rts),prm_traces)
  # save the prm_traces table under directory /prm_xic
  write.csv(prm_traces, file = paste0("traces/",analyte,"_prm_trace.csv"))

  spec <- data.frame(mz = sumspec@mz, intensity = sumspec@intensity)
  return(spec) # return full spectrum
}

# Run the main function if not interactive
if (!interactive()) {
  main()
}