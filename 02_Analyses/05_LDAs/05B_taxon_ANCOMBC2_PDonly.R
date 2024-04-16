# LDA taxon ANCOMBC2 PD only

base.dir <- ""
if (getwd() != base.dir) { setwd(base.dir) }
if (!exists("color.scales")) { source("_setup.R") }
redo.var <- "lda"
if (is.null(get.redo.state(redo.var))) {
  assign.redo(redo.var, state = T)
}

lda.data <- readRDS(files$pa_filt_mDs) # generated in 01_required_first.R

ancom.pdonly.taxon.res <- lapply(sources, function(source) {
  lapply(focal.sets$taxon, function(lvl) {
    vals <- track.values(source, lvl, return.val = T)
    element <- list(source, "taxon", lvl)
    mD <- lda.data %>% pluck(!!!element)
    if (is.null(mD)) {
      paste("Element at", vals, "is null. Skipping.") %>% 
        rlang::inform()
      next
    }
    track.values(source, lvl)
    mD <- lda.data %>% pluck(source, "taxon", lvl)
    assign.tbl <- get.assignments(mD)
    names(assign.tbl) %<>% str_replace("OTU", lvl)
    assign.col <- names(assign.tbl)[
      str_which(names(assign.tbl), lvl) - 1
    ]
    keep.cols <- str_subset(
      names(assign.tbl), 
      paste(lvl, assign.col, sep = "|")
    )
    assign.tbl <- assign.tbl[, ..keep.cols] %>% unique()
    names(assign.tbl) %<>%  str_replace(lvl, "Feature") %>% 
      str_replace(assign.col, "Assignment")
    
    ps <- microbData2phyloseq(mD)
    ancom.res.file <- file.path(
      dirs$saved, 
      paste0(
        "ancombc2_", 
        paste(source, lvl, sep = "_"),
        "_PDonly_results_dt.rds"
      )
    )
    res <- redo.if(redo.var, ancom.res.file, {
      ANCOMBC::ancombc2(
        data = ps, 
        fix_formula = diagnosis.col.name, 
        p_adj_method = "BY", 
        group = diagnosis.col.name, 
        prv_cut = prev.cut,
        struc_zero = TRUE, 
        neg_lb = TRUE, 
        pseudo_sens = FALSE,
        verbose = TRUE
      )
    })
    as.data.table(res$res) %>% 
      set_names(str_replace(names(.), "taxon", "Feature")) %>% 
      merge(assign.tbl, all.x = T, by = "Feature") %>%
      .[
        , `:=`(
          Source = source, 
          Feature.type = "Taxonomic", 
          Feature.lvl = lvl
        )
      ] %>%
      setcolorder(c(15:17, 1, 14)) %>% 
      return()
  }) %>% rbindlist()
}) %>% rbindlist()
warnings()


