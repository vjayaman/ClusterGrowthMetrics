Timedata <- R6Class(
  "Timedata", lock_objects = FALSE, 
  public = list(
    name = NULL, raw = NULL, flagged = NULL, cnames = NULL, 
    coded = NULL, melted = NULL, comps = NULL, status = NULL, 
    
    initialize = function(name, raw, isos, pad_height, pad_cluster) {
      self$name <- name
      self$start()
      self$raw <- raw
      self$coded <- codeIsolates(raw, name, isos, pad_height, pad_cluster)
      self$melted <- meltedIDs(raw, name, pad_height, pad_cluster)
      invisible(self)
    }, 
    start = function() {
      cat(paste0("  Constructing ", toupper(self$name), " data object:\n"))
    }, 
    coded_status = function(nov_code) {
      self$status <- self$coded %>% mutate(status = ifelse(isolate %in% nov_code, "novs", NA))
    }, 
    set_comps = function() {
      self$comps <- self$coded %>% 
        set_colnames(gsub(self$name, "tp", colnames(.))) %>% 
        compsSet(., toupper(self$name), indicate_progress = TRUE)
      invisible(self)
    }, 
    set_cnames = function() {
      self$cnames <- self$coded %>% pull(grep("h", colnames(self$coded), value = TRUE)) %>% unique() %>% sort()
      invisible(self)
    }, 
    flag_clusters = function() {
      self$flagged <- flaggingClusters(self$comps, self$name)
      invisible(self)
    }
  )
)

Heightdata <- R6Class(
  "Heightdata", lock_objects = FALSE, 
  public = list(
    h_before = NULL, h_after = NULL, comps = NULL, changed = tibble(), same = tibble(), 
    tracked = NULL, bef = NULL, aft = tibble(), results = NULL, 
    
    initialize = function(starter, t1_comps, hvals) {
      self$h_after <- self$h_before <- starter
      self$comps <- t1_comps %>% filter(tp1_h == starter) %>% arrange(tp1_h, tp1_cl)
      self$bef <- t1_comps %>% filter(tp1_h == self$h_before) %>% 
        set_colnames(c("h_bef", "cl_bef", "id_bef", "comp", "size_bef"))
      self$results <- vector(mode = "list", length = length(hvals)) %>% set_names(hvals)
      invisible(self)
    }, 
    clust_tracking = function(t2_comps, t2_cnames, t1_coded, t2_coded, indp) {
      self$changed <- trackClusters(self$comps, t2_comps, t2_cnames, t1_coded, t2_coded, indp)
      invisible(self)
    }, 
    post_data = function(t1_comps) {
      self$aft <- t1_comps %>% filter(tp1_h == self$h_after) %>% 
        set_colnames(c("h_aft", "cl_aft", "id_aft", "comp", "size_aft"))
      invisible(self)
    }, 
    unchanged = function() {
      self$same <- noChange(self$aft, self$bef, self$tracked)
    }, 
    update_iteration = function() {
      self$tracked <- bind_rows(self$changed, self$same) %>% arrange(tp1_h, tp1_cl)
      self$results[[self$h_after]] <- self$tracked
      invisible(self)
    }, 
    reset_values = function() {
      self$h_before <- self$h_after
      self$h_after <- NULL
      
      self$bef <- self$aft %>% set_colnames(c("h_bef", "cl_bef", "id_bef", "comp", "size_bef"))
      self$aft <- tibble()
    }
  )
)

