#! /usr/bin/env Rscript
require("optparse")
option_list <- list(
  make_option(c("-a", "--tp1"), metavar = "file", default = NULL, help = "Time point 1 file name"),
  make_option(c("-b", "--tp2"), metavar = "file", default = NULL, help = "Time point 2 file name")
)

arg <- parse_args(OptionParser(option_list=option_list))

# 1) Rscript env_setup.R <output data directory> <time point 1 dataset> <time point 2 dataset> 
#   - to install packages and set up required directory structure
# 2) Rscript nawc.R -i data/tp1-clusters_for_nawc.tsv -o data
#   - can use this script to identify the heights to plug into the following:
# 3) Rscript datacollection.R -a "data/timepoint1_data.csv" -b "data/timepoint2_data.csv" -x "5,10,15,20"

# This should be run first, to make sure the required packages are installed
msg <- file("logfile_env.txt", open="wt")
sink(msg, type="message")

pb <- txtProgressBar(min = 0, max = 4, initial = 1, style = 3)

message(paste0("For reporting an issue, see https://github.com/vjayaman/ClusterGrowthMetrics/issues.\n"))

required_packages <- c("tibble", "magrittr", "dplyr", "reshape2", "scales", "progress", "ggplot2", 
                       "plotly", "flexdashboard", "RColorBrewer", "stringr",  "optparse", "shiny", "DT")

not_installed <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]
setTxtProgressBar(pb, 2)

install.packages(not_installed, quiet = TRUE)
setTxtProgressBar(pb, 4)

# Testing packages were installed correctly:
x <- lapply(required_packages, require, character.only = TRUE)
names(x) <- required_packages

dir.create("outputs", showWarnings = FALSE)
dir.create("data", showWarnings = FALSE)

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

# OPTIONAL: DATA PREPARATION
library(tibble); library(magrittr)

convertAndSave <- function(ip, op) {
  df <- read.csv(ip, stringsAsFactors = FALSE, sep = ",", numerals = "no.loss") %>% as_tibble()
  m1 <- ncol(df)-2
  df %>% set_colnames(c("isolate", 0:m1)) %>% 
    write.table(., op, row.names = FALSE, quote = FALSE, sep = "\t")
}

convertAndSave(ip = arg$tp1, op = "data/tp1.csv")
convertAndSave(ip = arg$tp2, op = "data/tp2.csv")

cat("\nFormatted datasets for use in tracking scripts.\n")
close(pb)
