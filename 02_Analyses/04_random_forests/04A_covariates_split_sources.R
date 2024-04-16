# random forests covariates split sources

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
pa.filt.mDs <- readRDS(files$pa_filt_mDs) 
covars <- readRDS(file.path(dirs$saved, "rf_covars_vector.rds"))
####################################

covar.rf.mods.file <- file.path(
  dirs$saved,
  "randForest_covar_list.rds"
)
covar.rf.mods <- redo.if(
  redo.var,
  covar.rf.mods.file,
  max.save.mb = max.save.size,
  {
    lapply(sources, function(source) {
      set.seed(seed = 0)
      
      mod.data <- pa.filt.mDs %>% 
        pluck(source, "taxon", "Taxon") %>%
        get.metadata()
      mod.data <- mod.data[
        , .SD, .SDcols = c(diagnosis.col.name, covars)
      ]
      wf <- workflow() %>%
        add_recipe(build.rec.rf(mod.data, response = diagnosis.col.name)) %>%
        add_model(rf.spec)

      split.data <- mod.data %>%
        rsample::initial_split(strata = diagnosis.col.name, prop = 0.7)
      train.data <- training(split.data)
      test.data <- testing(split.data)

      max.mtry <- { dim(mod.data)[2] -1 } %>% sqrt() %>% floor()
      t.grid <- tune::tune_grid(
        wf,
        resamples = rsample::vfold_cv(
          train.data,
          v = 5,
          strata = diagnosis.col.name
        ),
        grid = expand.grid(
          mtry = round(
            seq(
              from = 1,
              to = max.mtry,
              length.out = min(c(max.mtry, 6))
            ),
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
          Full = mod.data,
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

if (get.redo.state(redo.var)) { stopCluster(cl) }

