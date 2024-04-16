# LDA taxon ANCOMBC2 with covariates

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

ancom.covars.taxon.res <- lapply(sources, function(source) {
  lapply(focal.sets$taxon, function(lvl) {
    vals <- track.values(source, lvl, return.val = T)
    element <- list(source, "taxon", lvl)
    mD <- lda.data %>% pluck(!!!element)
    if (is.null(mD)) {
      paste("Element at", vals, "is null. Skipping.") %>% 
        rlang::inform()
      return(NULL)
    }
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
