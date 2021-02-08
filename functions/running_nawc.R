
shannon <- compiler::cmpfun(function(clusters, base=2) {
  # Calculates the Shannon Entropy for a group of clusters uses base 2, returning bits by default
  p_i <- function(i) {
    sum(clusters == i) / length(clusters)
  }
  N <- unique(clusters)
  
  -sum(sapply(N, function(i) p_i(i) * log(p_i(i), base = base)))
})

neighbour_awc <- compiler::cmpfun(function(clusts, i) {
  # Calculates the Adjusted Wallace Coefficent of the ith column versus the i-1th column
  cur <- clusts[ , i]
  if ((i - 1) == 0 || length(unique(cur)) == 1) {
    result <- NA
  } else {
    suppressWarnings({
      result <- adj_wallace(cur, clusts[, i-1]) %>% use_series("Adjusted_Wallace_A_vs_B")
    })
  }
  result
})

singleton_proportion <- function(x) {
  # Determines what proportion of clusters have only a single member genome
  singletons <- sum(table(x) == 1)
  p_singleton <- singletons / length(unique(x))
  
  p_singleton
}

stats_df <- function(clusters) {
  # Calculates various statistics for clusters, and binds them together in a data.frame
  thresholds <- clusters %>% colnames %>% as.integer
  entropy <- sapply(clusters, shannon)
  
  nawc <- clusters %>% seq_along %>% sapply(neighbour_awc, clusts=clusters)
  n_clusts <- sapply(clusters, function(x) length(unique(x)))
  
  p_singletons <- sapply(clusters, singleton_proportion)
  
  data.frame("Threshold" = thresholds, "Neighbour AWC" = nawc, "Shannon (bits)" = entropy, 
             "Number of Clusters" = n_clusts, "Proportion of Singleton Clusters" = p_singletons, 
             check.names = FALSE) %>% return()
}

process <- function(input_file, outdir, delimiter) {
  if (! dir.exists(outdir)) {
    dir.create(outdir, recursive = TRUE)
  }
  delimiter = "\t"
  clusters <- input_file %>% read.table(sep = delimiter, row.names = 1, header = TRUE, check.names = FALSE)
  
  a2 <- cluster_stats <- stats_df(clusters)
  a2$Neighbour.AWC[is.na(a2$Neighbour.AWC)] <- 0
  
  # differences x[i+1] - x[i]
  # find the points that shrank from the left and grew to the right, note the indexing starts from 0
  y <- which(diff(sign(diff(a2$Neighbour.AWC))) > 0) + 1
  pt <- a2[y,] %>% filter(Neighbour.AWC < 0.95)
  
  a2$Selected.for.Analysis <- 0
  a2$Selected.for.Analysis[a2$Threshold %in% pt$Threshold] <- 1
  
  write.table(a2, file=file.path(outdir, "stats.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
}

