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
source_local("modelling_tools.R")
print("Initial loading successful")



cat('\n\nRunning PC-DFA routines\n');
options(warn=-1);
#remove rgl warningh5
options(rgl.useNULL = TRUE);
model_type <- "PCDFA" ## module name

cat("\nStart of the '", model_type, "' Galaxy module call: ",
    format(Sys.time(), "%a %d %b %Y %X"), "\n", sep="")

data_mat <- t(read.csv(listArguments[["dataMatrix_in"]], 
                                 header = TRUE,
                                 row.names = 1))
meta <- read.csv(listArguments[["sampleMetadata_in"]])
data_work <- data_mat
meta_work <- meta
if(listArguments[["test_flag"]] !="no"){
  data_test_mat <-t(read.csv(listArguments[["test_data"]],
                             header = TRUE,
                             row.names = 1))
  meta_test <- read.csv(listArguments[["test_meta"]])
  cat("\nTest set data has been loaded successfully.")
}

ns <- dim(data_work)[1]
nv <- dim(data_work)[2]

xaxis <- as.numeric(colnames(data_mat))
if (any(is.na(xaxis))) xaxis <- 1:nv

if (listArguments[["mos_flag"]] == "yes"){
  mos_threshold <- listArguments[["mos_thresh"]]
  if ("mos" %in% colnames(data_work)) {
    mos_now <- data_work[,'mos']
    null_idx <- which(mos_now < mos_threshold)
    data_work <- data_work[, -null_idx]
    meta_work <- meta[, -null_idx]
    cat('\n', length(null_idx), 'samples with low MOS have been removed. \n')
    cat('The indices of removed samples are: ', null_idx)} 
  else 
      cat('\n Although MOS threshold has been set there is no MOS in the data, MOS filter step skipped!\n')
}

if (listArguments[["test_flag"]] != "no"){
  if ("mos" %in% colnames(data_test_mat)){
    mos_now <- data_test_mat[,'mos']
    null_idx <- which(mos_now < mos_threshold)
    data_test_mat <- data_test_mat[, -null_idx]
    meta_test <- meta_test[, -null_idx]
    cat('\n', 'Additional ', length(null_idx), " samples with low MOS in the test set has been removed. \n")
    cat('The indices of removed test samples are: ', null_idx)
  }
  else if (listArguments[["mos_flag"]] == "yes")
    cat('\n Although MOS threshold has been set there is no MOS in the test set, MOS filter step skipped!\n')
}

if (listArguments[["trans_type"]] != "none"){
  trans_type <- listArguments[['trans_type']]
  data_work <- switch(trans_type,
                      log = log(data_work),
                      asinh = asinh(data_work),
                      sqrt = sqrt(data_work))
  cat('\n', trans_type, 'transformation has been applied to the data.\n')
  if (listArguments[["test_flag"]] !="none")
    data_test_mat  <- switch(trans_type,
                             log = log(data_test_mat), 
                             asinh = asinh(data_test_mat),
                             sqrt = sqrt(data_test_mat))
  cat('\n', trans_type, 'transformation has been applied to the test set.\n')
}
#  data_work <- as.data.frame(data_work)


if (listArguments[["scaling_type"]] != "none"){
  scaling_method <- listArguments[["scaling_type"]]
  data_scaled <- switch (scaling_method,
    auto = auto_sc(data_work),
    pareto = pareto_sc(data_work)
  )
  data_work <- data_scaled$data
  scale_factor <- data_scaled$scale_factor
  data_centre <- data_scaled$data_centre
  colnames(data_work) <- xaxis
  cat('\n' , scaling_method, 'scaling has been applied to the data .\n')
} else {
  scale_factor <- NULL
  data_centre <- NULL
}

no_pcs <- listArguments[['no_pcs']]
no_dfs <- listArguments[['no_dfs']]
cat('\n', no_pcs, 'principal components have been extracted and used for discriminant function analysis. \n')
if (dim(meta_work)[2] > 1)
  class_label <- meta_work[, 1] else class_label <-meta_work
if (is.data.frame(class_label)){
  class_label <- as.vector(class_label[, 1])
} else if(is.list(class_label)){
  class_label <- as.vector(unlist(class_label))
} else{
  class_label <- as.vector(class_label)
}
unique_class <- unique(class_label)
no_class <- length(unique_class)
class_no <- matrix(0, ns, 1)
for (i in 1:no_class) class_no[which(class_label == unique_class[i] )] <- i
pcdfa_model <- pcdfa(data_work, class_no, no_pcs, no_dfs, xmean = data_centre)
cat('\n', no_dfs, 'discriminant functions have been extracted. \n')

dfa_scores <- pcdfa_model$scores
if(listArguments[["test_flag"]] !="no"){
  test_scores <- pcdfa_pred(data_test_mat, model = pcdfa_model, scaling_factor = scale_factor)
  if (dim(meta_test)[2] > 1)
    label_test <- meta_test[, 1] else
    label_test <- meta_test
} else {
  test_scores <- NULL
  label_test <- NULL
  }

##saving
filename_data <- listArguments[["output_data"]]
filename_figures <- listArguments[["file_figures"]]
# colnames(data_to_save) <- row.names(data_work)

write.table("PC-DFA scores of the training set: \n", 
            file = filename_data,
            sep = ",",
            col.names = FALSE)
write.table(dfa_scores, 
            file = filename_data,
            quote = FALSE,
            sep = ",",
            append = TRUE
  )
write.table("\nPC-DFA loadings: \n",
            file = filename_data, 
            sep = ",",
            append = TRUE)
write.table(pcdfa_model$loadings,
            file = filename_data,
            sep = ",",
            append = TRUE)

if (listArguments[["test_flag"]] != "no"){
  write.table("PC-DFA scores of the test set \n",
              file = filename_data,
              sep = ",",
              col.names=FALSE,
              append = TRUE)
  write.table(test_scores, 
              file = filename_data,
              sep = ",",
              append = TRUE)
}

pdf(filename_figures, onefile = TRUE)
plot.pcdfa(dfa_scores, test_scores, meta_work, label_test)
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
