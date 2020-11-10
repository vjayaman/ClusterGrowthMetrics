# Useful functions for data reading and table processing, used in 
# all of the scripts in this project

# x <- c("tibble", "magrittr", "dplyr", "ggplot2", "reshape2","plotly", "tidyverse", 
#        "DT", "ggpubr", "rcompanion", "fitdistrplus", "logspline", "gamlss",
#        "goftest", "kableExtra", "EnvStats", "pracma", "e1071", "kernlab", "caret", 
#        "LiblineaR", "gridExtra", "transformr", "gifski", "gganimate", "fitdistrplus", 
#        "flexdashboard", "R.utils")
# lapply(x, require, character.only = TRUE)

x <- c("tibble", "magrittr", "dplyr", "reshape2", "scales")
lapply(x, require, character.only = TRUE)

get <- .Primitive("[[")

separateText <- function(txt, sp, i) {
  strsplit(txt, split = sp) %>% unlist() %>% extract2(i) %>% return()
}

# given the defining filename in:
# C:\Users\vasen\Documents\Pre-MSc courses\Honours Research Project\SummerProject-2020\, 
# read in the data
readData <- function(filename) {
  filename %>% 
    read.csv(file = ., numerals = "no.loss", stringsAsFactors = FALSE) %>% 
    as_tibble() %>% return()
}

# merge tables of 
# 1) isolates and cluster assignments at all heights, and
# 2) clusters and cluster sizes at all heights
# output: a melted tibble
isolateClusterSize <- function(clusters, sizes) {
  
  df_clusters <- reshape2::melt(clusters, id.vars = "isolate") %>% 
    as_tibble() %>% set_colnames(c("isolate", "height", "cluster"))
  
  df_sizes <- reshape2::melt(sizes, id.vars = "cluster") %>% 
    as_tibble() %>% set_colnames(c("cluster", "height", "size"))
  
  df_both <- left_join(df_clusters, df_sizes, by = c("height", "cluster"))
  df_both$height <- as.character(df_both$height)
  
  return(df_both)
}

clusterSizes <- function(clusters) { #tp_val
  df1 <- clusters %>% dplyr::select(-isolate)
  cs_sizes <- lapply(1:ncol(df1), function(i) {
    
    # if (i %% 200 == 0) {paste0(i, "/", ncol(df1)) %>% print()}
    
    h_x <- colnames(df1)[i]
    
    sizing <- df1[,h_x] %>% table() %>% as.data.frame() %>% as_tibble() %>% 
      set_colnames(c("cluster", "size")) %>% 
      add_column(height = h_x, .before = 1)
    sizing$cluster %<>% as.character() %>% as.numeric()
    
    return(sizing)
  }) %>% bind_rows()
  
  cs_sizes$id <- paste0(cs_sizes$height, "-", cs_sizes$cluster)
  # paste0("Done collecting cluster sizes for time point ", tp_val) %>% print()
  
  return(cs_sizes)
}

sparseCheck <- function(hist_x) {
  df <- hist_x[['counts']] %>% table() %>% as_tibble() %>% set_colnames(c("value","freq"))
  if (0 %in% df$value) {
    denom <- hist_x[['breaks']] %>% length()
    numerator <- df$freq[df$value == 0]
    sparse_val <- (numerator/denom) %>% scales::percent()
  }else {
    sparse_val <- "none"
  }
  return(sparse_val)
}

boxplotSubsetPoints <- function(df) {
  p <- ggplot(df, aes(x = heights, y = size_difference)) + geom_boxplot()
  p_data <- ggplot_build(p) %>% extract2("data")
  p_points <- tibble(heights = df$heights %>% unique(), 
                     median = p_data[[1]]$middle, q3 = p_data[[1]]$upper)
  return(p_points)
}

boxplotSubset <- function(df, title_sec) {
  p_points <- boxplotSubsetPoints(df)
  
  ggplot(df, aes(x = heights, y = size_difference)) + 
    geom_boxplot() + 
    geom_point(data = p_points, aes(x = heights, y = median, color = "red")) + 
    geom_point(data = p_points, aes(x = heights, y = q3, color = "blue")) +
    scale_x_discrete(guide = guide_axis(check.overlap = TRUE)) +
    scale_color_manual(name = "Point color", values = c("red"="red","blue"="blue"), 
                       labels = c("Quantile 3", "Median")) +
    xlab("Heights") + ylab("Change in cluster size") + 
    ggtitle(paste0("Change in cluster across thresholds", title_sec)) %>% return()
}

getUserInput <- function(msg) {
  cat(msg)
  op = readLines(con = "stdin", 1)
}

# checkDirInput <- function() {
#   cat(paste0("Please move any data files you\'d like ",
#              "to use to the newly created \'data\' directory. ",
#              "Hit Enter to continue."))
#   msg = readLines(con = "stdin", 1)
# }
# 
# checkDirInput()







