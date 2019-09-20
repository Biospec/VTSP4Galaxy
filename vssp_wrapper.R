#! Rscript

##########################################
# Signal processing script for Galaxy
##########################################
#Startup log
sink("startup_log.txt")

pkgs=c("RSpectra", "batch", "MASS", "Matrix", "pracma")

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

data_mat <- t(read.csv(listArguments[["dataMatrix_in"]], 
                                 header = TRUE,
                                 row.names = 1))
data_work <- data_mat
ns <- dim(data_work)[1]
nv <- dim(data_work)[2]

xaxis <- as.numeric(colnames(data_mat))
if (any(is.na(xaxis))) xaxis <- 1:nv

if (listArguments[["baseline"]] == "yes"){
  smooth_param <- as.numeric(listArguments[["baseline_param"]])
  for (i in 1:ns){
    z = arpls(data_work[i, ], lamda = smooth_param)
    data_work[i, ] = data_work[i, ] - z
  }
}

if (listArguments[["smooth"]] == "yes"){
  win_width <- as.numeric(listArguments$sm_param)
  data_work <- gaussian_sm(data_work, window_width = win_width)
}

if (listArguments[["norm_name"]] != "none"){
  norm_method <- listArguments[["norm_name"]]
  data_work <- switch (norm_method,
    snv = snv(data_work),
    tot = normalise_tot(data_work),
    emsc = emsc(data_work, p = listArguments[["norm_param"]])
  )
}

if (listArguments[["mos"]] == "yes") {
  mos_scores <- mos(as.matrix(data_work))
  ms <- mos_scores$ms
  data_work <- cbind(data_work, ms)
  colnames(data_work)[nv+1] <- "mos"
  data_work <- as.data.frame(data_work)
}

##saving
filename_data <- listArguments[["output_data"]]
filename_figures <- listArguments[["file_figures"]]

write.table(data_work, 
            file = filename_data,
            quote = FALSE,
            sep = ",")

pdf(filename_figures, onefile = TRUE)
if (length(which(colnames(data_work) == "mos")) == 1)
  data4plot <- data_work[, -which(colnames(data_work) == "mos")] else
  data4plot <- data_work
cat('colnames are :', xaxis, '\n')
cat('nv is ', nv, '\n')
plot(1:nv, data4plot[1, ], type = 'l')
for (i in 2:ns) points(1:nv, data4plot[i, ], type = 'l')

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
