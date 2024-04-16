# LDA taxon ALDEx2 PD only

base.dir <- ""
if (getwd() != base.dir) { setwd(base.dir) }
if (!exists("color.scales")) { source("_setup.R") }
redo.var <- "lda"
if (is.null(get.redo.state(redo.var))) {
  assign.redo(redo.var, state = T)
}

lda.data <- readRDS(files$pa_filt_mDs) # generated in 01_required_first.R

aldex.pdonly.taxon.res <- lapply(sources, function(source) {
  lapply(focal.sets$taxon, function(lvl) {
    track.values(source, lvl)
    mD <- lda.data %>% pluck(source, "taxon", lvl)
    assign.tbl <- get.assignments(mD)
    names(assign.tbl) %<>% str_replace("OTU", "Taxon")
    assign.col <- names(assign.tbl)[
      str_which(names(assign.tbl), lvl) - 1
    ]
    keep.cols <- str_subset(
      names(assign.tbl), 
      paste(lvl, assign.col, sep = "|")
    )
    assign.tbl <- assign.tbl[, ..keep.cols] %>% unique()
    names(assign.tbl) <- str_replace(
      names(assign.tbl), 
      lvl, 
      "Feature"
    ) %>% 
      str_replace(assign.col, "Assignment")
    
    aldex.mat <- get.abundances(mD) %>% t()
    
    aldex.res.file <- file.path(
      dirs$saved, 
      paste0(
        "aldex_", 
        paste(source, lvl, sep = "_"),
        "_PDonly_results_dt.rds"
      )
    )
    res <- redo.if(redo.var, aldex.res.file, {
      ALDEx2::aldex(
        reads = aldex.mat,
        conditions = get.metadata(mD)[[diagnosis.col.name]], 
        verbose = TRUE, 
        effect = TRUE
      )
    })
    as.data.table(res, keep.rownames = "Feature") %>% 
      merge(assign.tbl, all.x = T, by = "Feature") %>% 
      .[
        , `:=`(
          Source = source, 
          Feature.type = "Taxonomic", 
          Feature.lvl = lvl
        )
      ] %>%
      setcolorder(c(14:16, 1, 13)) %>% 
      return()
  }) %>% rbindlist()
}) %>% rbindlist()
warnings()