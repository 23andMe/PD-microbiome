# LDA function ANCOMBC2 PD only

base.dir <- ""
if (getwd() != base.dir) { setwd(base.dir) }
if (!exists("color.scales")) { source("_setup.R") }
redo.var <- "lda"
if (is.null(get.redo.state(redo.var))) {
  assign.redo(redo.var, state = T)
}

lda.data <- readRDS(files$pa_filt_mDs) # generated in 01_required_first.R

ancom.pdonly.func.res <- lapply(sources, function(source) {
  lapply(focal.sets$func, function(lvl) {
    track.values(source, lvl)
    mD <- lda.data %>% pluck(source, "func", lvl)
    if (lvl == "Fam") {
      relabs <- feature.relAbunds(mD)
      mD <- keep.features(mD, features = names(relabs[relabs >= abund.cut]))
    }
    ps <- microbData2phyloseq(mD)
    ancom.res.file <- file.path(
      dirs$saved,
      paste0(
        "ancombc2_",
        paste(source, lvl, sep = "_"),
        ifelse(lvl == "Fam", paste0("_abundCut", abund.cut), ""),
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
        verbose = TRUE,
        n_cl = 60
      )
    })
    as.data.table(res$res) %>%
      set_names(str_replace(names(.), "func", "Feature")) %>%
      .[
        , `:=`(
          Source = source,
          Feature.type = "Functional",
          Feature.lvl = lvl
        )
      ] %>%
      setcolorder(c(14:16, 1)) %>%
      return()
  }) %>% rbindlist()
}) %>% rbindlist()
warnings()
