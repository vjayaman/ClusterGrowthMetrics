stopwatch <- rep(0,2)
stopwatch[1] <- Sys.time()

source("global.R")
source("tabulating_functions.R")

## FOR USER: edit these variables as needed (move all data files to the 'data' directory)
tp1_filename <- "data/3692_2020-04-15_thresholds.csv"

# Then, using the same height for comparison, there are 215 multi-strain clusters at TP1. So the 
# number of multi-strain clusters grew by `r ((t2_cl_over_one - t1_cl_over_one)/t2_cl_over_one) %>% 
# scales::percent()`.

## Cluster metric generation

synthetic_filename <- "data/synthetic_tp1.csv"
if (file.exists(synthetic_filename)) {
  time1 <- read.csv(file = synthetic_filename, stringsAsFactors = FALSE, numerals = "no.loss") %>% as_tibble()
}else {
  stop(paste0("File at \'", synthetic_filename, "\' not found."))
}
X <- 1

if (file.exists(tp1_filename)) {
  time2 <- read.csv(file = tp1_filename, stringsAsFactors = FALSE, numerals = "no.loss") %>% as_tibble()
}else {
  stop(paste0("File at \'", tp1_filename, "\' not found."))
}
Y <- 2

# isolates at TP2 and not TP1
novel_isolates <- setdiff(time2$isolate, time1$isolate)
# isolates at TP1
isolates_t1 <- time1$isolate
# isolates at TP2
isolates_t2 <- time2$isolate

# Dataframes of the clusters at each time point and their associated sizes (easier for data retrieval later on).

sizes_for_tpX <- clusterSizes(time1, "1") %>% 
  set_colnames(c("TP1_h", "TP1_cl", "TP1_cl_size", "ID"))
sizes_for_tpY <- clusterSizes(time2, "2") %>% 
  set_colnames(c("TP2_h", "TP2_cl", "TP2_cl_size", "ID"))

# Reorganizing the data so the time points can be merged simply

# Melted time1 data: || isolate | height | cluster | id (height_cluster) ||
df1 <- melt(time1, id = "isolate") %>% as_tibble() %>% set_colnames(c("isolate", "height", "cluster"))
df1$height %<>% as.character()
df1$id <- paste0(df1$height, "-", df1$cluster)
df1_sizes <- df1 %>% set_colnames(c("Isolate", "TP1_h", "TP1_cl", "ID")) %>% 
  left_join(., sizes_for_tpX, by = c("TP1_h", "TP1_cl", "ID"))

# Melted time2 data: || isolate | height | cluster | id (height_cluster) ||
df2 <- melt(time2, id = "isolate") %>% as_tibble() %>% set_colnames(c("isolate", "height", "cluster"))
df2$height %<>% as.character()
df2$id <- paste0(df2$height, "-", df2$cluster)
df2_sizes <- df2 %>% set_colnames(c("Isolate", "TP2_h", "TP2_cl", "ID")) %>% 
  left_join(., sizes_for_tpY, by = c("TP2_h", "TP2_cl", "ID"))

# The foundational frame of cluster data:
#   - sizes and membership at each time point, 
#   - as well as the proportional increase from TP1 to TP2, in decimal form

df_all <- right_join(df1_sizes, df2_sizes, by = c("Isolate", "ID"))

all_clusters <- df_all %>% dplyr::select(Isolate, TP1_h, TP1_cl, TP2_h, TP2_cl, TP1_cl_size, TP2_cl_size)
all_clusters$TP1_cl_size[is.na(all_clusters$TP1_cl_size)] <- 0

all_clusters$prop_inc <- NA
all_clusters$prop_inc[all_clusters$TP1_cl_size == 0] <- 1
inds <- which(all_clusters$TP1_cl_size != 0)
all_clusters$prop_inc[inds] <- (all_clusters$TP2_cl_size[inds] - all_clusters$TP1_cl_size[inds]) / all_clusters$TP1_cl_size[inds]

## Adaptive threshold development

# We can plot what a couple of different threshold models can look like. 
# This part will be heavily altered depending on what the actual data looks like. 
# Currently being used to model synthetic data, and what we might expect significant 
# size change to look like. 

x <- c(1, 25, 50, 100, 150)
y <- c(300, 150, 75, 30, 15)

# http://r-statistics.co/Loess-Regression-With-R.html
loessMod1 <- loess(y ~ x, span = 1)
smoothed1 <- predict(loessMod1)

loessMod5 <- loess(y ~ x, span = 5)
smoothed5 <- predict(loessMod5)

lm_df <- tibble(x, y, smoothed1, smoothed5) %>%
  set_colnames(c("x", "Preset values", "Local regression (span 1)", "Local regression (span 5)")) %>%
  melt(id = "x") %>% as_tibble() %>% set_colnames(c("Cluster size", "Function", "Growth"))

# {ggplot(lm_df, aes(x = `Cluster size`, y = `Growth`, color = `Function`)) + geom_point() + 
#     theme(legend.position = "bottom") + scale_y_continuous(labels = scales::percent) + 
#     ylab("Growth (%)")} %>% ggplotly()

# The model selection step, and using it to determine what the adaptive threshold 
# function values would be for each input of initial cluster size. Anything that 
# needs to be extrapolated is set to a pre-determined plateau value of percent increase. 

model_used <- loessMod1
# predict.lm(model_used, data.frame(x = 30)) # note there are different types of prediction methods

predicted_y <- predict(model_used, newdata = all_clusters$TP1_cl_size)
na_predicted <- all_clusters$TP1_cl_size[which(is.na(predicted_y))]

if (all(na_predicted > max(lm_df$`Cluster size`))) {
  predicted_y[is.na(predicted_y)] <- 15
}

sizes_predicted <- predicted_y %>% round() %>% bind_cols(all_clusters, predicted = .)
sizes_predicted$predicted <- sizes_predicted$predicted*0.01

## Fold change
# Ranking the data by the actual proportional change over the adaptive threshold requirement, 
# to see by how much each cluster exceeds the growth prediction. The data is then saved to 
# the current working directory. 

sizes_predicted$fold_change <- sizes_predicted$prop_inc / sizes_predicted$predicted
sizes_predicted$fold_change %<>% round(., digits = 3)
sizes_predicted %<>% arrange(., -fold_change)

sizes_predicted$predicted %<>% scales::percent()
sizes_predicted$prop_inc %<>% scales::percent()

metrics <- sizes_predicted %>% 
  set_colnames(c("Isolate", "TP1 height", "TP1 cluster", "TP2 height", "TP2 cluster", 
                 "TP1 cluster size", "TP2 cluster size", "Proportional growth", 
                 "Adaptive threshold", "Fold change (Growth/Threshold)"))

# write.csv(metrics, "all_clusters_table.csv", row.names = FALSE)
# write.csv(metrics[1:10,], "cluster_metrics_first_ten.csv", row.names = FALSE)

stopwatch[2] <- Sys.time()

paste0("It takes ", round((stopwatch[2]-stopwatch[1])/60, digits = 2), " minutes to run.") %>% print()

