**Production branch (main) under maintenance**

ClusterGrowthMetrics
=====================================

Generates metrics for analyzing cluster growth across two time points, using isolates and cluster assignments at a series of heights at each time point. We start with a threshold in the time point 1 dataset, then track each cluster to find the "first" cluster at time point 2 that contains all of its isolates.

Prerequisites
-------------

```R
install.packages(c("tibble", "magrittr", "dplyr", "reshape2", "scales", "progress", 
                   "stringr", "ggplot2", "plotly", "optparse", "methods"))
```

Usage
----------

Clusters can be generated from [GrapeTree](https://github.com/achtman-lab/GrapeTree) via [this clustering script](https://github.com/dorbarker/grapetree_cluster). (How the clusters are generated is a domain-specific question). 

Generate clusters for two time point datasets (denoted TP1 and TP2) where TP2 represents the latter dataset, which must contain all the isolates present at TP1 as well as some additional isolates. 

There are two methods for running this script, method 1 uses command line and method 2 uses R (e.g. in RStudio).


Method 1
----------

To run the cluster growth metrics script from a **terminal**, navigate to the directory, then run: 
```sh
Rscript env_setup.R *<input data directory>* *<time point 1 dataset>* *<time point 2 dataset>*

e.g. $ Rscript env_setup.R "Desktop/inputs/" "tp1.csv" "tp2.csv"
```
Check `logfile_env.txt` if there are any issues, or if the script stops without a success message. 

After the relevant packages and data directories are set up, run the following. Input the full file paths, down to the file name extensions (no quotations necessary). They should be csvs (for now, allowing other file formats in the future is a WIP). 
```sh
Rscript datacollection.R -a *<time point 1 file path>* -b *<time point 2 file path>* -x *<comma-delimited list of thresholds to run>*

e.g. $ Rscript datacollection.R -a "data/tp1.csv" -b "data/tp2.csv" -x "5,10,25,30"
```
Progress bars will appear, and tracked cluster information will be saved to a file in a newly generated "outputs" directory. Check `logfile_datacollection.txt` if anything stops abruptly, or if you do not see a success message after the script stops running.

The next (optional) step is to run: 
```sh
Rscript tests/sampled_testing.R *<time point 1 data>* *<time point 2 data>* *<decimal (percent of clusters to sample for testing at each threshold)>*

e.g. $ Rscript tests/sampled_testing.R "data/tp1.csv" "data/tp2.csv" 0.25
```

Method 2
----------
To collect the cluster growth metrics from **within R**, first navigate to the cloned directory, then run
```r
source("env_setup.R")
```
After the environment is set up, open `datacollection.R` and update the filename input variables ("time1\_raw" and "time2\_raw"), to hold the paths of the datasets you want to use as input (they should hold time point 1 and time point 2 data, respectively).  Then run: 

**(Working on documenting this)**

By default all output files will be saved to an *outputs* directory, which will be created if it does not already exist. To change this, see the `saveData()` function in *functions/processing_functions.R*.


