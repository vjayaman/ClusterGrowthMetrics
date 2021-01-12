ClusterGrowthMetrics
=====================================

Generates metrics for analyzing cluster growth across two time points, using isolates and cluster assignments at a series of heights at each time point. We start with a threshold in the time point 1 dataset, then track each cluster to find the "first" cluster at time point 2 that contains all of its isolates.

Prerequisites
-------------

```R
install.packages(c("tibble", "magrittr", "dplyr", "reshape2", "scales", "progress", 
                       "ggplot2", "plotly", "flexdashboard", "RColorBrewer"))
```

Usage
----------

Clusters can be generated from [GrapeTree](https://github.com/achtman-lab/GrapeTree) via [this clustering script](https://github.com/dorbarker/grapetree_cluster)  

Generate clusters for two time point datasets (denoted TP1 and TP2) where TP2 represents the latter dataset, which must contain all the isolates present at TP1 as well as some additional isolates. 

There are two methods for running this script, method 1 is in command line and method 2 is in an R console window (e.g. in RStudio).


Method 1
----------

To run the cluster growth metrics script from a **terminal**, navigate to the directory, then run: 
```sh
Rscript env_setup.R
```
Check `logfile_env.txt` if there are any issues, or if the script stops without a success message. 

After the relevant packages and data directories are set up, run the following. Input the full file paths, down to the file name extensions (no quotations necessary). They should be csvs (for now, allowing other file formats in the future is a WIP). 
```sh
Rscript datacollection.R <time point 1 file path> <time point 2 file path>
```
Progress bars will appear, and tracked cluster information will be saved to a series of height-specific files in a newly generated "outputs" directory. Check `logfile_datacollection.txt` if anything stops abruptly, or if you do not see a success message after the script stops running.

The next (and last) step is to run: 
```sh
Rscript preparingmetrics.R <time point 1 file path> <time point 2 file path>
```
Follow the same formatting, and note that this part of the process takes just a few minutes to run (at most). Check `logfile_prepmetrics.txt` if there are any issues or if the script stops without a success message.


Method 2
----------
To collect the cluster growth metrics from **within R**, first navigate to the cloned directory, then run
```r
source("env_setup.R")
```
After the environment is set up, open `datacollection.R` and update the filename input variables ("time1_raw" and "time2_raw"), to hold the paths of the datasets you want to use as input (they should hold time point 1 and time point 2 data, respectively). Update the same two variables in the `preparingmetrics.R` file. Then run: 
```r
source("datacollection.R")
```
and when that is done, run:
```r
source("preparingmetrics.R")
```
By default all output files will be saved to an *outputs* directory, which will be created if it does not already exist. To change this, see the `saveData()` function in *functions/processing_functions.R*.



