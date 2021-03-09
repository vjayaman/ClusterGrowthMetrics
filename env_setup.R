#! /usr/bin/env Rscript

# 1) Rscript env_setup.R <output data directory> <time point 1 dataset> <time point 2 dataset> 
#   - to install packages and set up required directory structure
# 2) Rscript nawc.R -i data/tp1-clusters_for_nawc.tsv -o data
#   - can use this script to identify the heights to plug into the following:
# 3) Rscript datacollection.R -a "data/timepoint1_data.csv" -b "data/timepoint2_data.csv" -x "5,10,15,20"

# This should be run first, to make sure the required packages are installed
msg <- file("logfile_env.txt", open="wt")
sink(msg, type="message")

message(paste0("For reporting an issue, see https://github.com/vjayaman/ClusterGrowthMetrics/issues.\n"))

required_packages <- c("tibble", "magrittr", "dplyr", "reshape2", "scales", "progress", 
                       "stringr", "ggplot2", "plotly", "optparse", "methods", "R6", "testit")

not_installed <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]

install.packages(not_installed, quiet = TRUE)

# Testing packages were installed correctly:
x <- lapply(required_packages, require, character.only = TRUE)
names(x) <- required_packages

dir.create("outputs", showWarnings = FALSE)

if (all(unlist(x))) {
  cat("\nEnvironment set up successful.\n")
}else {
  if ("plotly" %in% names(which(x == FALSE))) {
    message("Linux libraries missing. Try: \n")
    message("   $ sudo apt-get update")
    message("   $ sudo apt-get install libcurl4-openssl-dev")
    message("   $ sudo apt-get install libssl-dev\n")
    message("Then run env_setup.R again.")
  }
  
  cat("\nNot all packages were installed successfully. Please see logfile_env.txt for details.")  
}

