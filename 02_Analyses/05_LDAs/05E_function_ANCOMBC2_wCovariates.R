# LDA function ANCOMBC2 with covariates

base.dir <- ""
if (getwd() != base.dir) { setwd(base.dir) }
if (!exists("color.scales")) { source("_setup.R") }
redo.var <- "lda"
if (is.null(get.redo.state(redo.var))) {
  assign.redo(redo.var, state = T)
}

lda.data <- readRDS(files$pa_filt_mDs) # generated in 01_required_first.R
covars <- unlist(ctrl.vars) %>% 
  unique() %>% 
  sort()

ancom.covars.func.res <- lapply(sources, function(source) {
  lapply(focal.sets$func, function(lvl) {
    track.values(source, lvl)
    element <- list(source, "func", lvl)
    mD <- lda.data %>% pluck(!!!element)
    if (lvl == "Fam") {
      relabs <- feature.relAbunds(mD)
      mD <- keep.features(mD, features = names(relabs[relabs >= abund.cut]))
    }
    # covars <- ctrl.vars %>% pluck(!!!element)
    all.terms <- c(mD@Sample.col, diagnosis.col.name, covars)
    mD <- get.metadata(mD)[, ..all.terms] %>%
      build.rec(response = diagnosis.col.name) %>%
      prep() %>%
      bake(new_data = get.metadata(mD)[, ..all.terms]) %>%
      as.data.table() %>%
      setkeyv(mD@Sample.col) %>%
      replace.metadata(mD = mD, new.tbl = .)
    ps <- microbData2phyloseq(mD)
    ancom.res.file <- file.path(
      dirs$saved,
      paste0(
        "ancombc2_",
        paste(source, lvl, sep = "_"),
        ifelse(lvl == "Fam", paste0("_abundCut", abund.cut), ""),
        "_wcovars_results_dt.rds"
      )
    )
    rand.frm <- covars %>%
      paste0("(1|", ., ")") %>%
      paste(collapse = " + ")
    res <- redo.if(redo.var, ancom.res.file, {
      ANCOMBC::ancombc2(
        data = ps,
        fix_formula = diagnosis.col.name,
        rand_formula = rand.frm,
        p_adj_method = "BY",
        group = diagnosis.col.name,
        prv_cut = prev.cut,
        struc_zero = TRUE,
        neg_lb = TRUE,
        pseudo_sens = FALSE,
        verbose = TRUE,
        n_cl = para.cores
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
