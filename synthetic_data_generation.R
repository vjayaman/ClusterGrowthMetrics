
# SUMMARY: this is the script for generating a time point 1 dataset (with a pre-determined rate of growth), 
#   given a time point 2 dataset 
#   - the time point 2 dataset should have one column "isolate" with labels for the genomes followed by 
#   a series of columns with cluster assignments at each threshold, and should be able to be read by:
#     read.csv(tp2_filename, stringsAsFactors = FALSE, numerals = "no.loss") %>% as_tibble()

stopwatch <- rep(0,2)
stopwatch[1] <- Sys.time()

source("global.R")
source("tabulating_functions.R")

## FOR USER: edit these variables as needed
tp2_filename <- "data/testing/3692_2020-04-15_thresholds.csv"

## Synthetic dataset development: 

if (file.exists(tp2_filename)) {
  tp2 <- read.csv(tp2_filename, stringsAsFactors = FALSE, numerals = "no.loss") %>% as_tibble()
}else {
  stop(paste0("File at \'", tp2_filename, "\' not found."))
}

melted_tp2 <- melt(tp2, id = "isolate") %>% as_tibble()
melted_tp2$h <- melted_tp2$variable %>% as.character() %>% gsub("h_", "", .) %>% as.integer()
colnames(melted_tp2) <- c("isolate", "height", "cluster", "num_h")

heights <- melted_tp2$num_h %>% unique()

# clusters and sizes for all heights for TP2
melted_tp2$height %<>% as.character()
cs_sizes <- lapply(colnames(tp2)[-1], function(h) {
  tp2[,h] %>% table() %>% as.data.frame() %>% set_colnames(c("cluster", "size")) %>% 
    as_tibble() %>% bind_cols(height = h, .) %>% return()
}) %>% bind_rows() %>% as_tibble()

cs_sizes$cluster %<>% as.integer()
# | isolate | character height | cluster | number height | cluster size |

# Looking at first ~300 heights to identify a useful height for the dropped genome analysis.
# Plotted along these heights to find a noticeable rise in cluster size before a plateau.

hvals <- colnames(tp2)[2:301]
isolate_select <- melted_tp2 %>% filter(height %in% hvals)
cs_select <- cs_sizes %>% filter(height %in% hvals)

### Plotting cluster numbers and sizes
# Zoom in along the x-axis to see in more detail where the plateaus are.

# plot_data <- left_join(isolate_select, cs_select, by = c("height", "cluster"))
# plot_data$height <- factor(plot_data$height, levels = unique(plot_data$height))
# # plot_data %<>% filter(num_h %in% 150:190)
# 
# ggplot(plot_data, aes(x = num_h, y = size, group = height)) + geom_boxplot() +
#   geom_vline(xintercept = 177, linetype = "dashed") + ylab("Cluster sizes") +
#   theme(axis.text.x = element_text(angle = 45)) + xlab("Heights") +
#   ggtitle("Close up of section with height used in synthetic TP1 generation")

# Height h_177 was selected, with 493 non-singleton clusters.
# TP2 clusters and sizes
h <- "h_177"
# clusters at TP2 for height 177, for all isolates
t2_height_select <- tp2 %>% dplyr::select(isolate, all_of(h))
sizes <- t2_height_select %>% pull(h) %>% table() %>% as.data.frame() %>% 
  as_tibble() %>% set_colnames(c("cluster", "size"))
sizes$cluster %<>% as.integer()
multi_strain <- sizes %>% filter(size > 1)
t2_cl_over_one <- multi_strain %>% nrow()

# Within those multistrain clusters, 60% of those isolates were randomly selected (1188 such isolates). 
# We then filter the TP2 dataset to exclude these isolates, this generates the TP1 dataset. Hence, 
# only 40% of the isolates found in multi-strain clusters at TP2 were originally found in the TP1 dataset.
x <- 0.6
tp2_ms <- t2_height_select %>% set_colnames(c("isolate", "cluster")) %>% 
  right_join(., multi_strain, by = "cluster")

isolates_to_rm <- sample(tp2_ms$isolate, round(x*nrow(tp2_ms)))
tp1 <- tp2 %>% dplyr::filter(!(isolate %in% isolates_to_rm))

# TP2 clusters and sizes
t1_height_select <- tp1 %>% dplyr::select(isolate, all_of(h))
t1_sizes <- t1_height_select %>% pull(h) %>% table() %>% 
  as.data.frame() %>% as_tibble() %>% set_colnames(c("cluster", "size"))
t1_cl_over_one <- t1_sizes %>% filter(size > 1) %>% nrow()

num_ms <- list("ms_t1" = t1_cl_over_one, "ms_t2" = t2_cl_over_one)

print("Note, writing over existing \'synthetic_tp1.csv\'")
write.csv(tp1, "synthetic_tp1.csv", row.names = FALSE)
