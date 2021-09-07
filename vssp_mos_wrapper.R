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
model_type <- "Vibrational Spectral signal processing - Morphological scores calculation" ## module name

cat("\nStart of the '", model_type, "' Galaxy module call: ",
    format(Sys.time(), "%a %d %b %Y %X"), "\n", sep="")

data_mat <- read.csv(listArguments[["dataMatrix_in"]], header = TRUE, row.names = NULL)
meta_mat <- read.csv(listArguments[["sampleMetadata_in"]], header = TRUE, row.names = NULL)

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

mos_scores <- mos(data_work)
ms <- mos_scores$ms
meta_names <- colnames(meta_mat)
meta_mat <- cbind(meta_mat, ms)
colnames(meta_mat) <- c(meta_names, "morphological_scores")
#  colnames(data_work) <- c(as.character(xaxis), "mos")
cat('\nMorphological scores are calculated and saved in the meta data. \n')

filename_meta <- listArguments[["output_meta"]]
filename_figures <- listArguments[["file_figures"]]

write.table(meta_mat,
            file = filename_meta,
            quote = FALSE,
            sep = ",",
            row.names = FALSE)
cat("\n Processed meta information saved. \n")

pdf(filename_figures, onefile = TRUE)


plot(1:ns, ms, type = 'h', xlab = "Sample id.", ylab = "Morphological scores")

cat("\n Bar chart of mos saved\n")

dev.off()


## ending
##-------

cat("\nEnd of the '", model_type, "' Galaxy module call: ",
    format(Sys.time(), "%a %d %b %Y %X"), "\n", sep = "")

sink()

rm(list = ls())
