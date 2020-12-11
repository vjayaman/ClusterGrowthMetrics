ClusterGrowthMetrics


**WIP:** do not download as is (in testing stage)

To run the cluster growth metrics script from a **terminal**, navigate to the directory, then run *Rscript env_setup.R*. After the relevant packages and data directories are setup, run *Rscript cluster_metrics.R <time point 1 file path> <time point 2 file path>*. Input the full file paths, down to the file name extensions and with no quotations. They should be csvs (for now, allowing other file formats in the future is a WIP). A progress bar will appear, and the metrics outputs will be saved to a newly generated outputs directory. Check the logfiles if anything stops abruptly, or if you do not see a success message (after both environment setup and metrics generation).

To collect the cluster growth metrics from **within R**, first navigate to the cloned directory, then run *source("env_setup.R")*. After the environment is set up, open *datacollection.R* and update the filename input variables if ("time1_raw" and "time2_raw"), to hold the paths of the datasets you want to use as input (they should hold time point 1 and time point 2 data, respectively). Update the same two variables in the *preparingmetrics.R* file, as well as the "running_in_rstudio" parameter if you want progress bars. Then run *> source("datacollection.R")* and when that is done, run *> source("preparingmetrics.R")*. By default all output files will be saved to an *outputs* directory, which will be created if it does not already exist. To change this, see the "saveData()" function in functions/processing_functions.R.
