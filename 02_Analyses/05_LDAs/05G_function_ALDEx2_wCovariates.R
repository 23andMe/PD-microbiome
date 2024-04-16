# LDA function ALDEx2 with covariates

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

aldex.covars.func.res <- lapply(sources, function(source) {
  lapply(focal.sets$func, function(lvl) {
    track.values(source, lvl)
    element <- list(source, "func", lvl)
    mD <- lda.data %>% pluck(!!!element)
    if (lvl == "Fam") {
      relabs <- feature.relAbunds(mD)
      mD <- keep.features(mD, features = names(relabs[relabs >= abund.cut]))
    }
    aldex.mat <- get.abundances(mD) %>% t()
    # covars <- ctrl.vars %>% pluck(!!!element)
    mm.frm <- paste(
      "~", diagnosis.col.name, "+",
      paste(covars, collapse = " + ")
    ) %>% as.formula()
    all.terms <- c(diagnosis.col.name, covars)
    metadata <- get.metadata(mD)[, ..all.terms] %>%
      build.rec(response = diagnosis.col.name) %>%
      prep() %>%
      bake(new_data = get.metadata(mD)[, ..all.terms]) %>%
      as.data.table()
    aldex.conds <- model.matrix(mm.frm, metadata)

    aldex.res.file <- file.path(
      dirs$saved,
      paste0(
        "aldex_",
        paste(source, lvl, sep = "_"),
        ifelse(lvl == "Fam", paste0("_abundCut", abund.cut), ""),
        "_wcovars_results_dt.rds"
      )
    )
    res <- redo.if(redo.var, aldex.res.file, {
      clr <- ALDEx2::aldex.clr(
        reads = aldex.mat,
        conds = aldex.conds,
        verbose = TRUE
      )
      test <- ALDEx2::aldex.glm(clr, aldex.conds)
      effect <- ALDEx2::aldex.glm.effect(clr)
      merge(
        as.data.table(test, keep.rownames = "Feature"),
        as.data.table(
          effect[[as.name(paste0(diagnosis.col.name, "PD case"))]],
          keep.rownames = "Feature"
        ),
        by = "Feature"
      )
    })
    if (!{"data.table" %in% class(res)}) {
      res <- as.data.table(res, keep.rownames = "Feature")
    }
    keep.cols <- str_subset(
      names(res),
      paste0(diagnosis.col.name, "|Feature|Intercept|rab|diff|effect|overlap")
    )
    res[, ..keep.cols] %>%
      .[
        , `:=`(
          Source = source,
          Feature.type = "Functional",
          Feature.lvl = lvl
        )
      ] %>%
      setcolorder(c(19:21, 1)) %>%
      return()
  }) %>% rbindlist()
}) %>% rbindlist()
warnings()
