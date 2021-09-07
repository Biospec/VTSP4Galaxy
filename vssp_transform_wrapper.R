#! Rscript

##########################################
# PC-DFA script for Galaxy
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


data_mat <- read.csv(listArguments[["dataMatrix_in"]], 
                     header = TRUE,
                     row.names = NULL)
col_ids <- colnames(data_mat)
if ("xaxis" %in% col_ids){
  xaxis_id <- which(col_ids == "xaxis")
  xaxis <- data_mat[, xaxis_id]
  data_mat <- data_mat[, -xaxis_id]
  col_ids <- col_ids[-xaxis_id]
} else{
  xaxis <- NULL
}

data_work <- t(as.matrix(data_mat))

ns <- dim(data_work)[1]
nv <- dim(data_work)[2]
cat("`\nThe imported data matrix has", ns, "samples and", nv, "variables.\n")

if (listArguments[["trans_name"]] != "none"){
  trans_type <- listArguments[['trans_name']]
  min_val = min(data_work)
  data_work <- switch(trans_type,
                      log = {if (min_val < 0){
                        cat("\nNegative number detected, log transformation cannot be applied,skip!")
                        data_work
                        } else if(min_val == 0){
                        cat("\nZero detected in the data matrix, an offset of 1 added\n")
                        log(data_work + 1)
                        } else log(data_work)},
                      asinh = { if (min_val <0) {
                        cat("\nNegative number detected, log transformation cannot be applied,skip!")
                        data_work
                        } else asinh(data_work)},
                      sqrt = { if (min_val <0) {
                        cat("\nNegative number detected, log transformation cannot be applied,skip!")
                        data_work
                        } else sqrt(data_work)}
                      )
  trans_name <- switch(trans_type,
                       log = "Natural Logarithm transformation",
                       asinh = "Inverse hyperbolic sine transformation",
                       sqrt = "Square root transformation"
                       )
  cat('\n', trans_name, ' transformation has been applied to the data.\n')
}

##saving
filename_data <- listArguments[["output_data"]]
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
cat("\nTransformed data saved. \n")
