ClusterGrowthMetrics

**Method 1:**

To run the cluster growth metrics script from a **terminal**, navigate to the directory, then run *Rscript env_setup.R*. Check "logfile_env.txt" if there are any issues, or if the script stops without a success message. 

After the relevant packages and data directories are setup, run *Rscript datacollection.R <time point 1 file path> <time point 2 file path>*. Input the full file paths, down to the file name extensions (no quotations necessary). They should be csvs (for now, allowing other file formats in the future is a WIP). Progress bars will appear, and tracked cluster information will be saved to a series of height-specific files in a newly generated "outputs" directory. Check "logfile_datacollection.txt" if anything stops abruptly, or if you do not see a success message after the script stops running.

The next (and last) step is to run *Rscript preparingmetrics.R <time point 1 file path> <time point 2 file path>*. Follow the same formatting, and note that this part of the process takes just a few minutes to run (at most). Check "logfile_prepmetrics.txt" if there are any issues or if the script stops without a success message.


**Method 2:**

To collect the cluster growth metrics from **within R**, first navigate to the cloned directory, then run *source("env_setup.R")*. After the environment is set up, open *datacollection.R* and update the filename input variables if ("time1_raw" and "time2_raw"), to hold the paths of the datasets you want to use as input (they should hold time point 1 and time point 2 data, respectively). Update the same two variables in the *preparingmetrics.R* file, as well as the "running_in_rstudio" parameter if you want progress bars. Then run *> source("datacollection.R")* and when that is done, run *> source("preparingmetrics.R")*. By default all output files will be saved to an *outputs* directory, which will be created if it does not already exist. To change this, see the "saveData()" function in functions/processing_functions.R.
