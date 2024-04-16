# random forests taxonomic abundances split types split sources

base.dir <- ""
if (getwd() != base.dir) { setwd(base.dir) }
if (!exists("color.scales")) { source("_setup.R") }
redo.var <- "rand.forest"
if (is.null(get.redo.state(redo.var))) {
  assign.redo(redo.var, state = T)
}

if (get.redo.state(redo.var)) {
  library(doParallel)
  cl <- makeForkCluster(para.cores)
  registerDoParallel(cl)
}

### generated in 01_required_first.R
taxon.cmb.rf.data.file <- file.path(
  dirs$saved,
  "taxon_combined_rf_data_list.rds"
)
taxon.cmb.rf.data <- readRDS(taxon.cmb.rf.data.file)
####################################

taxon.agg.cmb.covar.rf.mods.file <- file.path(
  dirs$saved,
  "randForest_taxon_aggregated_combined_covar_list.rds"
)
taxon.agg.cmb.covar.rf.mods <- redo.if(
  redo.var,
  taxon.agg.cmb.covar.rf.mods.file,
  max.save.mb = max.save.size,
  {
    ctrl.vars <- unique(unlist(ctrl.vars))
    repeat.cols <- c(ctrl.vars, diagnosis.col.name)
    lapply(cmb.methods, function(cmb.m) {
      track.values(cmb.m)
      lapply(agg.methods, function(agg.method) {
        set.seed(seed = 0)
        
        train.data <- NULL
        for (agg.m in agg.method) {
          track.values(agg.m)
          new.dt <- taxon.cmb.rf.data %>% 
            pluck(cmb.m, agg.m, "Train")
          if (is.null(train.data)) {
            train.data <- new.dt
          } else {
            train.data <- merge(
              train.data, 
              new.dt[
                , .SD, .SDcols = names(new.dt)[!{names(new.dt) %in% repeat.cols}]
              ], 
              by = internalID.col.name
            )
          }
        }
        train.data[, (internalID.col.name) := NULL]
        test.data <- NULL
        for (agg.m in agg.method) {
          track.values(agg.m)
          new.dt <- taxon.cmb.rf.data %>% 
            pluck(cmb.m, agg.m, "Test")
          if (is.null(test.data)) {
            test.data <- new.dt
          } else {
            test.data <- merge(
              test.data, 
              new.dt[
                , .SD, .SDcols = names(new.dt)[!{names(new.dt) %in% repeat.cols}]
              ], 
              by = internalID.col.name
            )
          }
        }
        test.data[, (internalID.col.name) := NULL]
        
        wf <- workflow() %>%
          add_recipe(build.rec.rf(train.data, response = diagnosis.col.name)) %>%
          add_model(rf.spec)
        
        max.mtry <- {ncol(train.data) - 1} %>% sqrt() %>% floor()
        t.grid <- tune::tune_grid(
          wf,
          resamples = rsample::vfold_cv(
            train.data,
            v = 5,
            strata = diagnosis.col.name
          ),
          grid = expand.grid(
            mtry = round(
              seq(from = 1, to = max.mtry, length.out = 6),
              0
            ),
            trees = c(10, 50, 100, 500, 1000, 1500),
            splitrule = c("gini", "extratrees")
          ),
          control = control_grid(
            verbose = TRUE,
            allow_par = TRUE,
            pkgs = c("tidyverse", "tidymodels", "ranger", "data.table")
          ),
          metrics = yardstick::metric_set(roc_auc)
        )
        best.tune <- select_best(t.grid, "roc_auc")
        mod <- finalize_workflow(wf, best.tune) %>%
          fit(train.data)
        test.probs <- predict(
          mod,
          test.data,
          type = "prob"
        ) %>%
          bind_cols(obs = test.data[[diagnosis.col.name]]) %>%
          bind_cols(predict(mod, test.data))
        roc.auc.ci <- pROC::roc(
          test.probs,
          obs,
          `.pred_PD case`
        ) %>% 
          pROC::ci.auc(method = "boot")
        list(
          Model = mod,
          Data = list(
            Full = rbind(train.data, test.data),
            Train = train.data,
            Test = test.data
          ),
          Classification = test.probs,
          Tune.res = best.tune,
          ROC.plot = autoplot(
            roc_curve(
              test.probs,
              obs,
              `.pred_PD case`,
              event_level = "second"
            )
          ),
          ROC.AUC = attributes(roc.auc.ci)$auc,
          ROC.AUC.CIs = roc.auc.ci[1:3]
        ) %>% return()
      })
    })
  })

if (get.redo.state(redo.var)) { stopCluster(cl) }