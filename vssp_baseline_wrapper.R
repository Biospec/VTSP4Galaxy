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
model_type <- "Vibrational Spectral signal processing - Baseline correction" ## module name

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

smooth_param <- as.numeric(listArguments[["baseline_param"]])
# for (i in 1:ns){
#   z = arpls(data_work[i, ], lamda = smooth_param)
#   data_work[i, ] = data_work[i, ] - z
# }

z = apply(data_work, 1, arpls, lamda = smooth_param)
data_work = data_work - t(z)

cat("\n The spectra were baseline corrected using arPLS with smooth parameter of ", listArguments[["baseline_param"]], '.\n')

filename_data <- listArguments[["output_data"]]
# filename_meta <- listArguments[["output_meta"]]
filename_figures <- listArguments[["file_figures"]]
if (!is.null(xaxis)){
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

# write.table(meta_mat,
#             file = filename_meta,
#             quote = FALSE,
#             sep = ",",
#             row.names = FALSE)
# cat("\n Processed meta information saved. \n")

pdf(filename_figures, onefile = TRUE)


data4plot <- data_work

plot(xaxis, data4plot[1, ], type = 'l', xlab = "xaxis or variable id.", ylab = "processed spectra")
for (i in 2:ns) points(xaxis, data4plot[i, ], type = 'l')

dev.off()

## ending
##-------

cat("\nEnd of the '", model_type, "' Galaxy module call: ",
    format(Sys.time(), "%a %d %b %Y %X"), "\n", sep = "")

sink()

rm(list = ls())
