# random_forests_taxon_import.R
base.dir <- ""
if (getwd() != base.dir) { setwd(base.dir) }
if (!exists("color.scales")) { source("_setup.R") }
redo.var <- "rand.forest"

covars <- unique(unlist(ctrl.vars))

all.rf.mod.files <- list.files(
  path = dirs$saved,
  pattern = "^randForest",
  full.names = TRUE
)

n.cores <- para.cores

if (get.redo.state(redo.var)) {
  library(doParallel)
  cl <- makeForkCluster(n.cores, outfile = "")
  registerDoParallel(cl)
}

permute.split.aucs <- function(rf.res, iters = 99) {
  if ("Full" %in% names(rf.res$Data)) {
    full.data <- rf.res$Data$Full
  } else {
    full.data <- rbind(rf.res$Data$Train, rf.res$Data$Test)
  }
  wf <- workflow() %>% 
    add_recipe(build.rec.rf(full.data, response = diagnosis.col.name)) %>% 
    add_model(rf.spec)
  # results <- list(ROC.AUC = NULL, ROC.curve = NULL)
  para.res <- foreach(i = 1:iters, .verbose = T) %dopar% {
    track.values(i)
    set.seed(seed = i)
    split.data <- full.data %>% 
      rsample::initial_split(strata = diagnosis.col.name, prop = 0.7)
    train.data <- training(split.data)
    test.data <- testing(split.data)
    res <- finalize_workflow(wf, rf.res$Tune.res) %>% 
      fit(train.data)
    test.probs <- predict(
      res, 
      test.data, 
      type = "prob"
    ) %>% 
      bind_cols(obs = test.data[[diagnosis.col.name]]) %>% 
      bind_cols(predict(res, test.data))
    
    list(
      ROC.AUC = {
        roc_auc(
          test.probs,
          obs,
          `.pred_PD case`,
          event_level = "second"
        ) %>%
          use_series(".estimate")
      },
      ROC.curve = {
        roc_curve(
          test.probs,
          obs,
          `.pred_PD case`,
          event_level = "second"
        ) %>%
          as.data.table() %>%
          .[, Iter := i]
      }
    ) %>% return()
  }
  list(
    ROC.AUC = sapply(para.res, extract2, "ROC.AUC"),
    ROC.curve = rbindlist(lapply(para.res, extract2, "ROC.curve"))
  ) %>% return()
}

# iters <- 4
iters <- 99

covar.rf.mods.files <- list.files(
  path = dirs$saved,
  pattern = "randForest_covar",
  full.names = TRUE
)
covar.rfs.aucs.file <- file.path(
  dirs$saved,
  "covarsOnly_rfs_permutedAUCs_dt.rds"
)
covar.rfs.aucs <- redo.if(redo.var, covar.rfs.aucs.file, {
  raw.aucs.file <- str_replace(covar.rfs.aucs.file, "AUCs", "rawAUCs")
  roc.curve.file <- str_replace(covar.rfs.aucs.file, "AUCs", "ROCcurves")
  if (file.exists(raw.aucs.file)) { file.remove(raw.aucs.file)}
  if (file.exists(roc.curve.file)) { file.remove(roc.curve.file)}
  
  lapply(covar.rf.mods.files, function(file) {
    covar.rf.mods <- readRDS(file)
    
    lapply(names(covar.rf.mods), function(idx) {
      track.values(file, idx)
      source <- ifelse(idx %in% sources, idx, "both")
      cmb.m <- ifelse(idx %in% cmb.methods, idx, "only")
      track.values(source, cmb.m)
      res <- permute.split.aucs(rf.res = covar.rf.mods[[idx]], iters = iters)
      
      auc.dt <- data.table(
        Source = source,
        Combine.method = cmb.m,
        Feature.set = "Base",
        Input.data = "Covariates",
        ROC.AUC = res$ROC.AUC
      )
      if (file.exists(raw.aucs.file)) {
        auc.dt <- rbind(readRDS(raw.aucs.file), auc.dt)
      }
      saveRDS(auc.dt, file = raw.aucs.file)
      roc.dt <- res$ROC.curve[
        , `:=`(
          Source = source,
          Combine.method = cmb.m,
          Feature.set = NA,
          Input.data = "Covariates"
        )
      ]
      if (file.exists(roc.curve.file)) {
        roc.dt <- rbind(readRDS(roc.curve.file), roc.dt)
      }
      saveRDS(roc.dt, file = roc.curve.file)
      data.table(
        Source = source,
        Combine.method = cmb.m,
        Feature.set = NA,
        Input.data = "Covariates",
        Min.ROC.AUC = round(min(res$ROC.AUC), 4),
        Mean.ROC.AUC = round(mean(res$ROC.AUC), 4),
        Max.ROC.AUC = round(max(res$ROC.AUC), 4)
      ) %>% return()
    }) %>% rbindlist()
  }) %>% rbindlist()
})

taxon.rf.mods.files <- list.files(
  path = dirs$saved,
  pattern = "randForest_taxon",
  full.names = TRUE
)
taxon.rfs.aucs.file <- file.path(
  dirs$saved,
  "taxon_rfs_permutedAUCs_dt.rds"
)
taxon.rfs.aucs <- redo.if(redo.var, taxon.rfs.aucs.file, {
  raw.aucs.file <- str_replace(taxon.rfs.aucs.file, "AUCs", "rawAUCs")
  roc.curve.file <- str_replace(taxon.rfs.aucs.file, "AUCs", "ROCcurves")
  if (file.exists(raw.aucs.file)) { file.remove(raw.aucs.file)}
  if (file.exists(roc.curve.file)) { file.remove(roc.curve.file)}
  
  lapply(taxon.rf.mods.files, function(file) {
    taxon.rf.mods <- readRDS(file)
    
    lapply(names(taxon.rf.mods), function(idx1) {
      lapply(names(taxon.rf.mods[[idx1]]), function(idx2) {
        track.values(file, idx1, idx2)
        mod.res <- taxon.rf.mods[[idx1]][[idx2]]
        source <- ifelse(idx1 %in% sources, idx1, "both")
        cmb.m <- ifelse(idx1 %in% cmb.methods, idx1, "only")
        in.data <- ifelse(
          idx2 %in% names(agg.methods),
          paste(str_split(idx2, "\\.")[[1]], collapse = " + "),
          idx2
        ) %>%
          str_replace("Taxon", "Strain") %>%
          paste("abundances")
        if (any(covars %in% names(mod.res$Data$Train))) {
          in.data %<>% paste("Covariates &", .)
        }
        track.values(source, cmb.m, in.data)
        res <- permute.split.aucs(mod.res, iters = iters)
        auc.dt <- data.table(
          Source = source,
          Combine.method = cmb.m,
          Feature.set = "Taxonomic",
          Input.data = in.data,
          ROC.AUC = res$ROC.AUC
        )
        if (file.exists(raw.aucs.file)) {
          auc.dt <- rbind(readRDS(raw.aucs.file), auc.dt)
        }
        saveRDS(auc.dt, file = raw.aucs.file)
        roc.dt <- res$ROC.curve[
          , `:=`(
            Source = source,
            Combine.method = cmb.m,
            Feature.set = "Taxonomic",
            Input.data = in.data
          )
        ]
        if (file.exists(roc.curve.file)) {
          roc.dt <- rbind(readRDS(roc.curve.file), roc.dt)
        }
        saveRDS(roc.dt, file = roc.curve.file)
        data.table(
          Source = source,
          Combine.method = cmb.m,
          Feature.set = "Taxonomic",
          Input.data = in.data,
          Min.ROC.AUC = round(min(res$ROC.AUC), 4),
          Mean.ROC.AUC = round(mean(res$ROC.AUC), 4),
          Max.ROC.AUC = round(max(res$ROC.AUC), 4)
        ) %>% return()
      }) %>%
        rbindlist() %>%
        return()
    }) %>% rbindlist()
  }) %>% rbindlist()
})

taxon.rf.mods.files <- list.files(
  path = dirs$saved,
  pattern = "randForest_taxon",
  full.names = T
)
taxon.sig.impt.dt.file <- file.path(
  dirs$saved,
  "rf_sig_impt_features_taxon.rds"
)
taxon.sig.impt.dt <- redo.if(
  redo.var,
  taxon.sig.impt.dt.file,
  {
    best.data <- taxon.rfs.aucs[Mean.ROC.AUC == max(Mean.ROC.AUC)]
    idx1 <- ifelse(
      best.data$Source == "both",
      best.data$Combine.method,
      best.data$Source
    )
    file.pattern <- ifelse(str_detect(best.data$Input.data, "\\+"), "aggregated", "")
    if (best.data$Source == "both") {
      file.pattern %<>% paste0("_combined")
    }
    if (str_detect(best.data$Input.data, "\\&")) {
      file.pattern %<>% paste0("_covar")
    }
    file.pattern %<>% paste0("_list.rds")
    mod.list <- str_subset(taxon.rf.mods.files, file.pattern) %>% readRDS()
    idx2s <- names(mod.list[[idx1]])
    
    if (any(str_detect(idx2s, "\\."))) {
      idx2 <- str_split(idx2s, "\\.") %>%
        lapply(function(x) {
          sapply(x, function(y) {
            str_detect(
              str_replace_all(best.data$Input.data, "Strain", "Taxon"), 
              y
            )
          }) %>% sum()
        }) %>%
        unlist() %>%
        which.max() %>%
        idx2s[.]
    } else {
      idx2 <- str_which(idx2s, best.data$Input.data)
    }
    best.mod.list <- mod.list[[idx1]][[idx2]]
    best.mod <- extract_fit_engine(best.mod.list$Model)
    impt.data <- build.rec.rf(
      best.mod.list$Data$Train, 
      response = diagnosis.col.name
    ) %>%
      prep() %>%
      bake(best.mod.list$Data$Train)
    impt.dt <- ranger::importance_pvalues(
      best.mod,
      method = "altmann",
      formula = as.formula(paste(diagnosis.col.name, "~ .")),
      data = impt.data
    ) %>% as.data.table(keep.rownames = "Feature")
    sig.impt.dt <- impt.dt[pvalue <= 0.05][
      order(importance, decreasing = T)
    ]
    sig.impt.dt[
      , `:=`(
        Source = best.data$Source,
        Feature.type = "taxon",
        RF.data = best.data$Input.data
      )
    ]
    setcolorder(sig.impt.dt, tail(1:ncol(sig.impt.dt), 3))
    sig.impt.dt
  })

func.rf.mods.files <- list.files(
  path = dirs$saved,
  pattern = "randForest_func",
  full.names = TRUE
)
func.rfs.aucs.file <- file.path(
  dirs$saved,
  "func_rfs_permutedAUCs_dt.rds"
)
func.rfs.aucs <- redo.if(redo.var, func.rfs.aucs.file, {
  raw.aucs.file <- str_replace(func.rfs.aucs.file, "AUCs", "rawAUCs")
  roc.curve.file <- str_replace(func.rfs.aucs.file, "AUCs", "ROCcurves")
  if (file.exists(raw.aucs.file)) { file.remove(raw.aucs.file)}
  if (file.exists(roc.curve.file)) { file.remove(roc.curve.file)}
  
  lapply(func.rf.mods.files, function(file) {
    func.rf.mods <- readRDS(file)
    
    lapply(names(func.rf.mods), function(idx1) {
      lapply(names(func.rf.mods[[idx1]]), function(idx2) {
        track.values(file, idx1, idx2)
        mod.res <- func.rf.mods[[idx1]][[idx2]]
        source <- ifelse(idx1 %in% sources, idx1, "both")
        cmb.m <- ifelse(idx1 %in% cmb.methods, idx1, "only")
        in.data <- ifelse(
          idx2 %in% names(agg.methods),
          paste(str_split(idx2, "\\.")[[1]], collapse = " + "),
          idx2
        ) %>%
          paste("abundances")
        if (any(covars %in% names(mod.res$Data$Train))) {
          in.data %<>% paste("Covariates &", .)
        }
        track.values(source, cmb.m, in.data)
        res <- permute.split.aucs(mod.res, iters = iters)
        auc.dt <- data.table(
          Source = source,
          Combine.method = cmb.m,
          Feature.set = "Functional",
          Input.data = in.data,
          ROC.AUC = res$ROC.AUC
        )
        if (file.exists(raw.aucs.file)) {
          auc.dt <- rbind(readRDS(raw.aucs.file), auc.dt)
        }
        saveRDS(auc.dt, file = raw.aucs.file)
        roc.dt <- res$ROC.curve[
          , `:=`(
            Source = source,
            Combine.method = cmb.m,
            Feature.set = "Functional",
            Input.data = in.data
          )
        ]
        if (file.exists(roc.curve.file)) {
          roc.dt <- rbind(readRDS(roc.curve.file), roc.dt)
        }
        saveRDS(roc.dt, file = roc.curve.file)
        data.table(
          Source = source,
          Combine.method = cmb.m,
          Feature.set = "Functional",
          Input.data = in.data,
          Min.ROC.AUC = round(min(res$ROC.AUC), 4),
          Mean.ROC.AUC = round(mean(res$ROC.AUC), 4),
          Max.ROC.AUC = round(max(res$ROC.AUC), 4)
        ) %>% return()
      }) %>%
        rbindlist() %>%
        return()
    }) %>% rbindlist()
  }) %>% rbindlist()
})

if (get.redo.state(redo.var)) { stopCluster(cl) }

aucs.dt <- list.files(
  path = dirs$saved,
  pattern = "^(covarsOnly_rfs|func|taxon)_.+rawAUCs",
  full.names = T
) %>%
  lapply(readRDS) %>%
  rbindlist(fill = T)