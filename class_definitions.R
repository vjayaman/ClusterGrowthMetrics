Timedata <- R6Class(
  "Timedata", lock_objects = FALSE, 
  public = list(
    name = NULL, raw = NULL, isos = NULL, 
    initialize = function(name, raw, isos, coded, melted, comps) {
      self$name <- name
      self$raw <- raw
      self$isos <- raw %>% pull(isolate)
      self$coded <- codeIsolates(raw, name, isos)
      self$melted <- meltedIDs(raw, name)
      self$start()
    }, 
    start = function() {cat(paste0("  Initialized object for ", toupper(self$name), "\n"))}, 
    set_comps = function(coded_data) {
      self$comps <- compsSet(coded_data, toupper(self$name), indicate_progress = TRUE)
    }
  )
)

Heightdata <- R6Class(
  "Heightdata", lock_objects = FALSE, 
  public = list(
    h_before = NULL, h_after = NULL, comps = NULL, changed = NULL, same = NULL, 
    tracked = NULL, bef = NULL, aft = NULL, 
    
    initialize = function(h_before, t1_comps) {
      self$h_before <- h_before
      self$comps <- t1_comps %>% filter(tp1_h == h_before) %>% arrange(tp1_h, tp1_cl)
    }, 
    clust_tracking = function(t2_comps, t2_cnames, t1_coded, t2_coded, indp) {
      self$changed <- self$comps %>% 
        trackClusters(., t2_comps, t2_cnames, t1_coded, t2_coded, indp) %>% 
        createID(., "tp1", "tp1_h", "tp1_cl")
      invisible(self)
    }, 
    prior_data = function(t1_comps) {
      self$bef <- t1_comps %>% filter(tp1_h == self$h_before) %>% 
        set_colnames(c("h_bef", "cl_bef", "id_bef", "comp", "size_bef"))
    }, 
    post_data = function(t1_comps) {
      self$aft <- t1_comps %>% filter(tp1_h == self$h_after) %>% 
        set_colnames(c("h_aft", "cl_aft", "id_aft", "comp", "size_aft"))
      invisible(self)
    }, 
    unchanged = function() {self$same <- noChange(self$aft, self$bef, self$tracked)}, 
    update_iteration = function() {
      self$h_before <- self$h_after
      self$bef <- self$aft %>% set_colnames(c("h_bef", "cl_bef", "id_bef", "comp", "size_bef"))
    }, 
    update_tracking = function() {
      self$tracked <- self$changed
      invisible(self)
    }, 
    add_flag = function() {
      self$tracked <- self$changed %>% add_column(flag = self$changed$id)
      invisible(self)
    }, 
    saveTempFile = function(t1_coded, hval, op) {
      m1 <- t1_coded$tp1_h %>% max() %>% nchar()
      m2 <- formatC(as.integer(hval), width = m1, format = "d", flag = "0")
      saveRDS(self$tracked, file.path(op, paste0("h", m2, ".Rds")))
    }
  )
)
