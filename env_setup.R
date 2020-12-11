#! /usr/bin/env Rscript

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

cat("\nEnvironment set up successful.\n")
# sink()
