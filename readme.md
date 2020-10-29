ClusterGrowthMetrics

To run the cluster growth metrics script *metrics_all_clusters.R* from **within R**, first run *env_setup.R* in a new R project, and then move the data files you would like to use into the newly created "data" directory. Edit the filename variables near the top of the *metrics_all_clusters.R* script, and then run *source("metrics_all_clusters.R")*.  


To run the cluster growth metrics script from a **terminal**, navigate to the directory, then run *Rscript env_setup.R*. Move the data files you would like to use as input into the newly generated "data" directory. Run *Rscript cluster_metrics.R*. You will be prompted for each input file you would like to use (the time point 1 and 2 datasets). Input as the filename with the extension, no quotations. They should be csvs (for now, allowing other file formats is a WIP). For example: *tp1_dataset.csv*, then hit *Enter*. A progress bar should appear, and the metrics outputs will be saved to the data directory. 


