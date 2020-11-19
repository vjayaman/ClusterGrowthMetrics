ClusterGrowthMetrics


**WIP:** do not download as is (in testing stage)

To run the cluster growth metrics script from a **terminal**, navigate to the directory, then run *Rscript env_setup.R*. After the relevant packages and data directories are setup, run *Rscript cluster_metrics.R <time point 1 file path> <time point 2 file path>*. Input the full file paths, down to the file name extensions and with no quotations. They should be csvs (for now, allowing other file formats is a WIP). A progress bar should appear, and the metrics outputs will be saved to a newly generated outputs directory. Check the logfiles if anything stops abruptly, or if you do not see a success message (after both environment setup and metrics generation).

To collect the cluster growth metrics from **within R**, first navigate to the clone directory, the run *source("env_setup.R")*. After the environment is set up, open *metrics_rmarkdown_version.Rmd* and edit the parameters at the top of the file, "data1" and "data2", to hold the paths of the datasets you want to use as input (they should hold time point 1 and time point 2 data, respectively). Then knit the document, making sure you have write permissions enabled. The document will create an "outputs" directory and save the metrics there. Note that this method is in development.
