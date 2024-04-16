# LDA function ALDEx2 PD only

base.dir <- ""
if (getwd() != base.dir) { setwd(base.dir) }
if (!exists("color.scales")) { source("_setup.R") }
redo.var <- "lda"
if (is.null(get.redo.state(redo.var))) {
  assign.redo(redo.var, state = T)
}

lda.data.file <- file.path(
  dirs$saved, 
  "func_lda_data_list.rds"
)
lda.data <- readRDS(files$pa_filt_mDs) # generated in 01_required_first.R

aldex.pdonly.func.res <- lapply(sources, function(source) {
  lapply(focal.sets$func, function(lvl) {
    track.values(source, lvl)
    mD <- lda.data %>% pluck(source, "func", lvl)
    if (lvl == "Fam") {
      relabs <- feature.relAbunds(mD)
      mD <- keep.features(mD, features = names(relabs[relabs >= abund.cut]))
    }
    aldex.mat <- get.abundances(mD) %>% t()

    aldex.res.file <- file.path(
      dirs$saved,
      paste0(
        "aldex_",
        paste(source, lvl, sep = "_"),
        ifelse(lvl == "Fam", paste0("_abundCut", abund.cut), ""),
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
      .[
        , `:=`(
          Source = source,
          Feature.type = "Functional",
          Feature.lvl = lvl
        )
      ] %>%
      setcolorder(c(13:15, 1)) %>%
      return()
  }) %>% rbindlist()
}) %>% rbindlist()
warnings()
