#! Rscript

##########################################
# Signal processing script for Galaxy
##########################################
#Startup log
sink("startup_log.txt")

pkgs=c("RSpectra", "batch", "MASS", "Matrix", "pracma", "signal")

for(pkg in pkgs) {
  suppressPackageStartupMessages( stopifnot( library(pkg, quietly=TRUE, logical.return=TRUE, character.only=TRUE)))
  cat(pkg,"\t",as.character(packageVersion(pkg)),"\n",sep="")
}

listArguments <- parseCommandArgs(evaluate = FALSE) #interpretation of arguments given in command line as an R list of objects
sink()

#Redirect all stdout to the log file
filename_log <- listArguments[["output_log"]]
sink(filename_log)



# ----- PACKAGE -----
cat("\tPACKAGE INFO\n")
sessionInfo()

source_local <- function(fname) {
  argv <- commandArgs(trailingOnly = FALSE)
  base_dir <- dirname(substring(argv[grep("--file=", argv)], 8))
  source(paste(base_dir, fname, sep="/"))
}

# Load functions
source_local("signal_processing.R")
print("Initial loading successful")



cat('\n\nRunning signal processing routines\n');
options(warn=-1);
#remove rgl warningh5
options(rgl.useNULL = TRUE);
model_type <- "Vibrational Spectral signal processing" ## module name

cat("\nStart of the '", model_type, "' Galaxy module call: ",
    format(Sys.time(), "%a %d %b %Y %X"), "\n", sep="")

data_mat <- read.csv(listArguments[["dataMatrix_in"]], header = TRUE, row.names = NULL)
meta_mat <- read.csv(listArguments[["sampleMetadata_in"]], header = TRUE, row.names = NULL)

col_ids <- colnames(data_mat)


xaxis <- as.numeric(data_mat[, 1])
data_mat <- t(as.matrix(data_mat[, -1]))
col_ids <- col_ids[-1]

data_work <- data_mat
# colnames_work <- colnames(data_work)
# rownames_work <- rownames(data_work)
ns <- dim(data_work)[1]
nv <- dim(data_work)[2]

if (is.null(xaxis)) xaxis <- 1:nv

if (listArguments[["baseline"]] == "yes"){
  smooth_param <- as.numeric(listArguments[["baseline_param"]])
  for (i in 1:ns){
    z = arpls(data_work[i, ], lamda = smooth_param)
    data_work[i, ] = data_work[i, ] - z
  }
#  colnames(data_work) <- as.character(xaxis)
  cat('\nBaseline correction performed with smooth parameter = ', smooth_param, '\n')
}

if (listArguments[["smooth"]] == "yes"){
  sg_win <- as.numeric(listArguments$sm_win)
  sg_order <- as.numeric(listArguments$sm_order)
  sg_m <- as.numeric(listArguments$sm_derv)
  for (i in 1:ns)
    data_work[i, ] <- sgolayfilt(data_work[i, ], p=sg_order, n=sg_win, m = sg_m)
#  colnames(data_work) <- as.character(xaxis)
  cat('\nSignal smoothing/derivative performed using Savizky-Golay filter with window width = ', sg_win, 
      ', polynomial order = ', sg_order, ' and k = ', sg_m, '. .\n')
}

if (listArguments[["norm_name"]] != "none"){
  norm_method <- listArguments[["norm_name"]]
  data_work <- switch (norm_method,
    snv = snv(data_work),
    tot = normalise_tot(data_work),
    emsc = emsc(data_work, p = listArguments[["norm_param"]])
  )
  #data_work <- as.data.frame(data_work, col.names = as.characters(xaxis))
  cat('\nSpectra were normalized using ', norm_method, 'with the parameter of ', listArguments[["norm_param"]], '.\n')
}

if (listArguments[["mos"]] == "yes") {
  mos_scores <- mos(as.matrix(data_work))
  ms <- mos_scores$ms
  meta_names <- colnames(meta_mat)
  meta_mat <- cbind(meta_mat, ms)
  colnames(meta_mat) <- c(meta_names, "morphological_scores")
#  colnames(data_work) <- c(as.character(xaxis), "mos")
  cat('\nMorphological scores are calculated and saved in the meta data. \n')
}

##saving
filename_data <- listArguments[["output_data"]]
filename_meta <- listArguments[["output_meta"]]
filename_figures <- listArguments[["file_figures"]]
if (is.null(xaxis)){
  data_to_save <- cbind(xaxis, t(data_work))
  colnames(data_to_save) <- c('xaxis', col_ids)
} else {
  data_to_save <- t(data_work)
  colnames(data_to_save) <- col_ids
}

write.table(data_to_save, 
            file = filename_data,
            quote = FALSE,
            sep = ",",
            row.names = FALSE
  )
cat("\n Processed data saved. \n")

write.table(meta_mat,
            file = filename_meta,
            quote = FALSE,
            sep = ",",
            row.names = FALSE)
cat("\n Processed meta information saved. \n")

pdf(filename_figures, onefile = TRUE)


data4plot <- data_work

plot(xaxis, data4plot[1, ], type = 'l', xlab = "xaxis", ylab = "processed spectra")
for (i in 2:ns) points(xaxis, data4plot[i, ], type = 'l')

dev.off()


tryCatch({
  save.image(file="processed_spc.RData");
}, warning = function(w) {
  print(paste("Warning: ", w));
}, error = function(err) {
  stop(paste("ERROR saving result RData object:", err));
});

## ending
##-------

cat("\nEnd of the '", model_type, "' Galaxy module call: ",
    format(Sys.time(), "%a %d %b %Y %X"), "\n", sep = "")

sink()

rm(list = ls())
