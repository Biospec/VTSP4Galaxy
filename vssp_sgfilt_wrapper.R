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
model_type <- "Vibrational Spectral signal processing - Savizky-Golay process" ## module name

cat("\nStart of the '", model_type, "' Galaxy module call: ",
    format(Sys.time(), "%a %d %b %Y %X"), "\n", sep="")

data_mat <- read.csv(listArguments[["dataMatrix_in"]], header = TRUE, row.names = NULL)
# meta_mat <- read.csv(listArguments[["sampleMetadata_in"]], header = TRUE, row.names = NULL)

col_ids <- colnames(data_mat)

if (col_ids[1] == 'xaxis') {
  xaxis <- as.numeric(data_mat[, 1])
  data_mat <- t(as.matrix(data_mat[, -1]))
  col_ids <- col_ids[-1]
} else {
  data_mat <- t(as.matrix(data_mat))
  xaxis <- NULL
} 

data_work <- data_mat
# colnames_work <- colnames(data_work)
# rownames_work <- rownames(data_work)
ns <- dim(data_work)[1]
nv <- dim(data_work)[2]

sg_win <- as.numeric(listArguments$sm_win)
sg_order <- as.numeric(listArguments$sm_order)
sg_m <- as.numeric(listArguments$sm_derv)
for (i in 1:ns)
  data_work[i, ] <- sgolayfilt(data_work[i, ], p=sg_order, n=sg_win, m = sg_m)
#  colnames(data_work) <- as.character(xaxis)
cat('\nSignal smoothing/derivative performed using Savizky-Golay filter with window width = ', sg_win, 
    ', polynomial order = ', sg_order, ' and k = ', sg_m, '. .\n')


filename_data <- listArguments[["output_data"]]
filename_figures <- listArguments[["file_figures"]]
if (!is.null(xaxis)){
  data_to_save <- cbind(xaxis, t(data_work))
  colnames(data_to_save) <- c('xaxis', col_ids)
} else {
  data_to_save <- t(data_work)
  cat("\nThe number of variables in final data is", dim(data_to_save)[1], ", the number of samples in final data is ", dim(data_to_save)[2], ".\n")
  colnames(data_to_save) <- col_ids
}

write.table(data_to_save, 
            file = filename_data,
            quote = FALSE,
            sep = ",",
            row.names = FALSE
)
cat("\n Processed data saved. \n")

# write.table(meta_mat,
#             file = filename_meta,
#             quote = FALSE,
#             sep = ",",
#             row.names = FALSE)
# cat("\n Processed meta information saved. \n")

pdf(filename_figures, onefile = TRUE)


data4plot <- data_work
if (is.null(xaxis)){
  xaxis <- 1:nv
  plot(xaxis, data4plot[1, ], type = 'l', xlab = "xaxis or variable id.", ylab = "processed spectra")
  for (i in 2:ns) points(xaxis, data4plot[i, ], type = 'l')
  
} else {
  plot(xaxis, data4plot[1, ], type = 'l', xlab = "xaxis or variable id.", ylab = "processed spectra")
  for (i in 2:ns) points(xaxis, data4plot[i, ], type = 'l')
}

dev.off()

## ending
##-------

cat("\nEnd of the '", model_type, "' Galaxy module call: ",
    format(Sys.time(), "%a %d %b %Y %X"), "\n", sep = "")

sink()

rm(list = ls())
