# random forests taxonomic abundances covariates split type split sources

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
taxon.rf.data.file <- file.path(
  dirs$saved,
  "taxon_rf_data_list.rds"
)
taxon.rf.data <- readRDS(taxon.rf.data.file)
####################################

taxon.covar.rf.mods.file <- file.path(
  dirs$saved,
  "randForest_taxon_covar_list.rds"
)
taxon.covar.rf.mods <- redo.if(
  redo.var,
  taxon.covar.rf.mods.file,
  max.save.mb = max.save.size,
  {
    lapply(sources, function(source) {
      lapply(taxon.lvls, function(lvl) {
        track.values(source, lvl)
        set.seed(seed = 0)
        
        train.data <- taxon.rf.data %>% 
          pluck(source, lvl, "Train")
        train.data[, (internalID.col.name) := NULL]
        test.data <- taxon.rf.data %>% 
          pluck(source, lvl, "Test")
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
