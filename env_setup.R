#! /usr/bin/env Rscript
input_args = commandArgs(trailingOnly = TRUE)

# 1) Rscript env_setup.R - to install packages and set up required directory structure
# 2) Rscript data_prep.R: on the files in data/, including data/european-t1_clusters.csv
# 3) Rscript nawc.R -i data/european-t1-clusters_for_nawc.tsv -o data
# 4) Rscript datacollection.R data/timepoint1_data.csv data/timepoint2_data.csv
# 5) Rscript result_files.R -a data/timepoint1_data.csv -b data/timepoint2_data.csv -t "nawc"

# This should be run first, to make sure the required packages are installed
msg <- file("logfile_env.txt", open="wt")
sink(msg, type="message")

pb <- txtProgressBar(min = 0, max = 4, initial = 1, style = 3)

message(paste0("For reporting an issue, see https://github.com/vjayaman/ClusterGrowthMetrics/issues.\n"))

required_packages <- c("tibble", "magrittr", "dplyr", "reshape2", "scales", "progress", 
                       "ggplot2", "plotly", "flexdashboard", "RColorBrewer")

not_installed <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]
setTxtProgressBar(pb, 2)

install.packages(not_installed, quiet = TRUE)
setTxtProgressBar(pb, 4)

# Testing packages were installed correctly:
x <- lapply(required_packages, require, character.only = TRUE)
dir.create("outputs", showWarnings = FALSE)
dir.create("data")

cat("\nEnvironment set up successful.\n")
# sink()

# OPTIONAL: DATA PREPARATION
library(tibble); library(magrittr)

convertAndSave <- function(datadir, ip, op) {
  file.path(datadir, ip) %>% 
    read.csv(file = ., stringsAsFactors = FALSE, numerals = "no.loss") %>% 
    as_tibble() %>% 
    set_colnames(c("isolate", 1:(ncol(.)-1))) %>% 
    write.table(., file.path(datadir, op), row.names = FALSE, quote = FALSE, sep = "\t")  
}

# convertAndSave(datadir = "data", ip = "european-t1_clusters.csv", op = "european-t1-clusters_for_nawc.tsv")
# convertAndSave(datadir = "data", ip = "european-t1_clusters.csv", op = "timepoint1.csv")
# convertAndSave(datadir = "data", ip = "european-t2_clusters.csv", op = "timepoint2.csv")

convertAndSave(datadir = input_args[1], ip = input_args[2], op = "tp1-clusters_for_nawc.tsv")
convertAndSave(datadir = input_args[1], ip = input_args[2], op = "timepoint1.csv")
convertAndSave(datadir = input_args[1], ip = input_args[3], op = "timepoint2.csv")

cat("\nFormatted datasets for use in tracking scripts.\n")
close(pb)