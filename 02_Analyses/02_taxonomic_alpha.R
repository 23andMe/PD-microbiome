# taxonomic_alpha.pred.PD.R

base.dir <- ""
if (getwd() != base.dir) { setwd(base.dir) }
if (!exists("color.scales")) { source("_setup.R") }
redo.var <- "tax.alpha"
if (is.null(get.redo.state(redo.var))) {
  assign.redo(redo.var, state = T)
}

focal.alphas <- alphas
mDs <- readRDS(files$rarefied_filt_mDs)
all.covars <- unlist(ctrl.vars) %>% 
  unique() %>% 
  sort()
tax.ft <- "Taxon"
stool.taxa <- mDs$stool$taxon$Taxon@Assignments


glm.spec <- logistic_reg() %>% set_engine("glm")

alpha.pred.PD.glms.file <- file.path(
  dirs$saved,
  "taxa_alpha_pred_PD_glms_list.rds"
)
alpha.pred.PD.glms <- redo.if(
  redo.var, 
  alpha.pred.PD.glms.file, 
  max.save.mb = max.save.size, 
  {
    source <- "saliva"
    element <- list(source, "taxon", tax.ft)
    mD <- pluck(mDs, !!!element)
    # all.covars <- pluck(ctrl.vars, !!!element)
    mod.data <- get.metadata(mD)
    mod.data <- mod.data[
      , .SD, .SDcols = c(diagnosis.col.name, focal.alphas, all.covars)
    ]
    frms <- NULL
    mods <- NULL
    i <- 1
    frms[i] <- paste(
      diagnosis.col.name, "~", 
      "Chao1 + InvSimpson +", 
      paste(all.covars, collapse = " + ")
    )
    mods[[i]] <- workflow() %>% 
      add_recipe(build.rec(mod.data, response = diagnosis.col.name)) %>% 
      add_model(glm.spec, formula = as.formula(frms[i])) %>% 
      fit(mod.data) %>% 
      extract_fit_engine()
    summary(mods[[i]]) %>% print()
    Anova(mods[[i]]) %>% print()
    
    i <- 2
    frms[i] <- paste(
      diagnosis.col.name, "~", 
      "InvSimpson +", 
      paste(all.covars, collapse = " + ")
    )
    mods[[i]] <- workflow() %>% 
      add_recipe(build.rec(mod.data, response = diagnosis.col.name)) %>% 
      add_model(glm.spec, formula = as.formula(frms[i])) %>% 
      fit(mod.data) %>% 
      extract_fit_engine()
    summary(mods[[i]]) %>% print()
    Anova(mods[[i]]) %>% print()
    do.call(anova, c(mods, test = "Chisq")) %>% print()
    
    i <- 3
    frms[i] <- paste(
      diagnosis.col.name, "~", 
      "InvSimpson +", 
      paste(all.covars, collapse = " + "),
      "+", paste0("InvSimpson:", wkly.BMs.col.name)
    )
    mods[[i]] <- workflow() %>% 
      add_recipe(build.rec(mod.data, response = diagnosis.col.name)) %>% 
      add_model(glm.spec, formula = as.formula(frms[i])) %>% 
      fit(mod.data) %>% 
      extract_fit_engine()
    summary(mods[[i]]) %>% print()
    Anova(mods[[i]]) %>% print()
    do.call(anova, c(mods[2:3], test = "Chisq")) %>% print()
    
    i <- 4
    frms[i] <- paste(
      diagnosis.col.name, "~", 
      "InvSimpson +", 
      paste(all.covars, collapse = " + "),
      "+", paste0("InvSimpson:", flossing.col.name)
    )
    mods[[i]] <- workflow() %>% 
      add_recipe(build.rec(mod.data, response = diagnosis.col.name)) %>% 
      add_model(glm.spec, formula = as.formula(frms[i])) %>% 
      fit(mod.data) %>% 
      extract_fit_engine()
    summary(mods[[i]]) %>% print()
    Anova(mods[[i]]) %>% print()
    do.call(anova, c(mods[c(2, 4)], test = "Chisq")) %>% print()
    mod.stats <- new.env()
    mod.stats$brs.plots <- NULL
    mod.stats$hoslem.plots <- NULL
    lapply(seq_along(mods), function(i) {
      mod <- mods[[i]]
      frm <- frms[i]
      mod.stats$brs.plots <- binnedplot_tjfs(
        fitted(mod),
        residuals(mod, type = 'response'),
        main = paste0("Frm: ", str_wrap(frm, width = 50))
      ) %>% 
        list() %>% 
        c(mod.stats$brs.plots,  .)
      hoslem.t <- ResourceSelection::hoslem.test(
        as.integer(mod.data[[diagnosis.col.name]]) - 1, 
        fitted(mod)
      )
      mod.stats$hoslem.plots <- {
        data.table(
          obs_p = hoslem.t$observed[, 2], 
          exp_p = hoslem.t$expected[, 2]
        ) %>% 
          ggplot(aes(x = obs_p, y = exp_p)) +
          geom_point() + 
          geom_smooth() +
          geom_abline(intercept = 0, slope = 1, linewidth = 0.5) +
          labs(subtitle = paste("Frm:", str_wrap(frm, width = 50)))
      } %>% 
        list() %>% 
        c(mod.stats$hoslem.plots, .)
      chisq.dev <- pchisq(
        q = mod$null.deviance - mod$deviance, 
        df = length(coef(mod))
      ) %>% 
        subtract(1, .)
      p.rsq <- { mod$deviance / mod$null.deviance } %>% 
        subtract(1, .)
      data.table(
        Frm = frm,
        Deviance.pval = chisq.dev,
        Hoslem.pval = hoslem.t$p.value,
        Pseudo.Rsq = p.rsq,
        AIC = mod$aic 
      )
    }) %>% 
      rbindlist() %>% 
      print()
    fig.file <- file.path(dirs$plots, "alpha_saliva_glms_binnedResids.png")
    cowplot::plot_grid(plotlist = mod.stats$brs.plots) %>% 
      ggsave(filename = fig.file, width = 12, height = 10, dpi = 300)
    fig.file <- file.path(dirs$plots, "alpha_saliva_glms_hoslemPlots.png")
    cowplot::plot_grid(plotlist = mod.stats$hoslem.plots) %>% 
      ggsave(filename = fig.file, width = 12, height = 10, dpi = 300)
    saliva.mod <- mods[[2]]
    
    source <- "stool"
    element <- list(source, "taxon", tax.ft)
    mD <- pluck(mDs, !!!element)
    # all.covars <- pluck(ctrl.vars, !!!element)
    mod.data <- get.metadata(mD)
    mod.data <- mod.data[
      , .SD, .SDcols = c(diagnosis.col.name, focal.alphas, all.covars)
    ]
    
    frms <- NULL
    mods <- NULL
    i <- 1
    frms[i] <- paste(
      diagnosis.col.name, "~", 
      "Chao1 + InvSimpson +", 
      paste(all.covars, collapse = " + ")
    )
    mods[[i]] <- workflow() %>% 
      add_recipe(build.rec(mod.data, response = diagnosis.col.name)) %>% 
      add_model(glm.spec, formula = as.formula(frms[i])) %>% 
      fit(mod.data) %>% 
      extract_fit_engine()
    summary(mods[[i]]) %>% print()
    Anova(mods[[i]]) %>% print()
    
    i <- 2
    frms[i] <- paste(
      diagnosis.col.name, "~", 
      "Chao1 +", 
      paste(all.covars, collapse = " + ")
    )
    mods[[i]] <- workflow() %>% 
      add_recipe(build.rec(mod.data, response = diagnosis.col.name)) %>% 
      add_model(glm.spec, formula = as.formula(frms[i])) %>% 
      fit(mod.data) %>% 
      extract_fit_engine()
    summary(mods[[i]]) %>% print()
    Anova(mods[[i]]) %>% print()
    do.call(anova, c(mods, test = "Chisq")) %>% print()
    
    i <- 3
    frms[i] <- paste(
      diagnosis.col.name, "~", 
      "Chao1 +", 
      paste(all.covars, collapse = " + "),
      "+", paste0("Chao1:", wkly.BMs.col.name)
    )
    mods[[i]] <- workflow() %>% 
      add_recipe(build.rec(mod.data, response = diagnosis.col.name)) %>% 
      add_model(glm.spec, formula = as.formula(frms[i])) %>% 
      fit(mod.data) %>% 
      extract_fit_engine()
    summary(mods[[i]]) %>% print()
    Anova(mods[[i]]) %>% print()
    do.call(anova, c(mods[2:3], test = "Chisq")) %>% print()
    
    i <- 4
    frms[i] <- paste(
      diagnosis.col.name, "~", 
      "Chao1 +", 
      paste(all.covars, collapse = " + "),
      "+", paste0("Chao1:", age.col.name)
    )
    mods[[i]] <- workflow() %>% 
      add_recipe(build.rec(mod.data, response = diagnosis.col.name)) %>% 
      add_model(glm.spec, formula = as.formula(frms[i])) %>% 
      fit(mod.data) %>% 
      extract_fit_engine()
    summary(mods[[i]]) %>% print()
    Anova(mods[[i]]) %>% print()
    do.call(anova, c(mods[c(2, 4)], test = "Chisq")) %>% print()
    mod.stats <- new.env()
    mod.stats$brs.plots <- NULL
    mod.stats$hoslem.plots <- NULL
    lapply(seq_along(mods), function(i) {
      mod <- mods[[i]]
      frm <- frms[i]
      mod.stats$brs.plots <- binnedplot_tjfs(
        fitted(mod),
        residuals(mod, type = 'response'),
        main = paste0("Frm: ", str_wrap(frm, width = 50))
      ) %>% 
        list() %>% 
        c(mod.stats$brs.plots,  .)
      hoslem.t <- ResourceSelection::hoslem.test(
        as.integer(mod.data[[diagnosis.col.name]]) - 1, 
        fitted(mod)
      )
      mod.stats$hoslem.plots <- {
        data.table(
          obs_p = hoslem.t$observed[, 2], 
          exp_p = hoslem.t$expected[, 2]
        ) %>% 
          ggplot(aes(x = obs_p, y = exp_p)) +
          geom_point() + 
          geom_smooth() +
          geom_abline(intercept = 0, slope = 1, linewidth = 0.5) +
          labs(subtitle = paste("Frm:", str_wrap(frm, width = 50)))
      } %>% 
        list() %>% 
        c(mod.stats$hoslem.plots, .)
      chisq.dev <- pchisq(
        q = mod$null.deviance - mod$deviance, 
        df = length(coef(mod))
      ) %>% 
        subtract(1, .)
      p.rsq <- { mod$deviance / mod$null.deviance } %>% 
        subtract(1, .)
      data.table(
        Frm = frm,
        Deviance.pval = chisq.dev,
        Hoslem.pval = hoslem.t$p.value,
        Pseudo.Rsq = p.rsq,
        AIC = mod$aic 
      )
    }) %>% 
      rbindlist() %>%
      print()
    fig.file <- file.path(dirs$plots, "alpha_stool_glms_binnedResids.png")
    cowplot::plot_grid(plotlist = mod.stats$brs.plots) %>% 
      ggsave(filename = fig.file, width = 12, height = 10, dpi = 300)
    fig.file <- file.path(dirs$plots, "alpha_stool_glms_hoslemPlots.png")
    cowplot::plot_grid(plotlist = mod.stats$hoslem.plots) %>% 
      ggsave(filename = fig.file, width = 12, height = 10, dpi = 300)
    stool.mod <- mods[[3]]
    list(saliva = saliva.mod, stool = stool.mod)
  })
###
{
  source <- "saliva"
  element <- list(source, "taxon", tax.ft)
  mD <- pluck(mDs, !!!element)
  # all.covars <- pluck(ctrl.vars, !!!element)
  mod.data <- get.metadata(mD) %>% 
    merge(batches.dt, by = externalID.col.name)
  mod.data <- mod.data[
    , .SD, .SDcols = c(diagnosis.col.name, focal.alphas, all.covars, batchID.col.name)
  ]
  frms <- NULL
  mods <- NULL
  i <- 1
  frms[i] <- paste(
    diagnosis.col.name, "~", 
    "Chao1 + InvSimpson +", batchID.col.name, "+", 
    paste(all.covars, collapse = " + ")
  )
  mods[[i]] <- workflow() %>% 
    add_recipe(build.rec(mod.data, response = diagnosis.col.name)) %>% 
    add_model(glm.spec, formula = as.formula(frms[i])) %>% 
    fit(mod.data) %>% 
    extract_fit_engine()
  summary(mods[[i]]) %>% print()
  Anova(mods[[i]]) %>% print()
  
  i <- 2
  frms[i] <- paste(
    diagnosis.col.name, "~", 
    "InvSimpson +", batchID.col.name, "+", 
    paste(all.covars, collapse = " + ")
  )
  mods[[i]] <- workflow() %>% 
    add_recipe(build.rec(mod.data, response = diagnosis.col.name)) %>% 
    add_model(glm.spec, formula = as.formula(frms[i])) %>% 
    fit(mod.data) %>% 
    extract_fit_engine()
  summary(mods[[i]]) %>% print()
  Anova(mods[[i]]) %>% print()
  do.call(anova, c(mods, test = "Chisq")) %>% print()
  
  i <- 3
  frms[i] <- paste(
    diagnosis.col.name, "~", 
    "InvSimpson +", batchID.col.name, "+", 
    paste(all.covars, collapse = " + "),
    "+", paste0("InvSimpson:", wkly.BMs.col.name)
  )
  mods[[i]] <- workflow() %>% 
    add_recipe(build.rec(mod.data, response = diagnosis.col.name)) %>% 
    add_model(glm.spec, formula = as.formula(frms[i])) %>% 
    fit(mod.data) %>% 
    extract_fit_engine()
  summary(mods[[i]]) %>% print()
  Anova(mods[[i]]) %>% print()
  do.call(anova, c(mods[2:3], test = "Chisq")) %>% print()
  
  i <- 4
  frms[i] <- paste(
    diagnosis.col.name, "~", 
    "InvSimpson +", batchID.col.name, "+", 
    paste(all.covars, collapse = " + "),
    "+", paste0("InvSimpson:", flossing.col.name)
  )
  mods[[i]] <- workflow() %>% 
    add_recipe(build.rec(mod.data, response = diagnosis.col.name)) %>% 
    add_model(glm.spec, formula = as.formula(frms[i])) %>% 
    fit(mod.data) %>% 
    extract_fit_engine()
  summary(mods[[i]]) %>% print()
  Anova(mods[[i]]) %>% print()
  do.call(anova, c(mods[c(2, 4)], test = "Chisq")) %>% print()
  mod.stats <- new.env()
  mod.stats$brs.plots <- NULL
  mod.stats$hoslem.plots <- NULL
  lapply(seq_along(mods), function(i) {
    mod <- mods[[i]]
    frm <- frms[i]
    mod.stats$brs.plots <- binnedplot_tjfs(
      fitted(mod),
      residuals(mod, type = 'response'),
      main = paste0("Frm: ", str_wrap(frm, width = 50))
    ) %>% 
      list() %>% 
      c(mod.stats$brs.plots,  .)
    hoslem.t <- ResourceSelection::hoslem.test(
      as.integer(mod.data[[diagnosis.col.name]]) - 1, 
      fitted(mod)
    )
    mod.stats$hoslem.plots <- {
      data.table(
        obs_p = hoslem.t$observed[, 2], 
        exp_p = hoslem.t$expected[, 2]
      ) %>% 
        ggplot(aes(x = obs_p, y = exp_p)) +
        geom_point() + 
        geom_smooth() +
        geom_abline(intercept = 0, slope = 1, linewidth = 0.5) +
        labs(subtitle = paste("Frm:", str_wrap(frm, width = 50)))
    } %>% 
      list() %>% 
      c(mod.stats$hoslem.plots, .)
    chisq.dev <- pchisq(
      q = mod$null.deviance - mod$deviance, 
      df = length(coef(mod))
    ) %>% 
      subtract(1, .)
    p.rsq <- { mod$deviance / mod$null.deviance } %>% 
      subtract(1, .)
    data.table(
      Frm = frm,
      Deviance.pval = chisq.dev,
      Hoslem.pval = hoslem.t$p.value,
      Pseudo.Rsq = p.rsq,
      AIC = mod$aic 
    )
  }) %>% 
    rbindlist() %>% 
    print()
  fig.file <- file.path(dirs$plots, "alpha_saliva_glms_binnedResids.png")
  cowplot::plot_grid(plotlist = mod.stats$brs.plots) %>% 
    ggsave(filename = fig.file, width = 12, height = 10, dpi = 300)
  fig.file <- file.path(dirs$plots, "alpha_saliva_glms_hoslemPlots.png")
  cowplot::plot_grid(plotlist = mod.stats$hoslem.plots) %>% 
    ggsave(filename = fig.file, width = 12, height = 10, dpi = 300)
  saliva.mod <- mods[[2]]
  
  source <- "stool"
  element <- list(source, "taxon", tax.ft)
  mD <- pluck(mDs, !!!element)
  # all.covars <- pluck(ctrl.vars, !!!element)
  mod.data <- get.metadata(mD) %>% 
    merge(batches.dt, by.x = externalID.col.name, by.y = "Sample.ID")
  mod.data <- mod.data[
    , .SD, .SDcols = c(diagnosis.col.name, focal.alphas, all.covars, batchID.col.name)
  ]
  
  frms <- NULL
  mods <- NULL
  i <- 1
  frms[i] <- paste(
    diagnosis.col.name, "~", 
    "Chao1 + InvSimpson +", batchID.col.name, "+", 
    paste(all.covars, collapse = " + ")
  )
  mods[[i]] <- workflow() %>% 
    add_recipe(build.rec(mod.data, response = diagnosis.col.name)) %>% 
    add_model(glm.spec, formula = as.formula(frms[i])) %>% 
    fit(mod.data) %>% 
    extract_fit_engine()
  summary(mods[[i]]) %>% print()
  Anova(mods[[i]]) %>% print()
  
  i <- 2
  frms[i] <- paste(
    diagnosis.col.name, "~", 
    "Chao1 +", batchID.col.name, "+", 
    paste(all.covars, collapse = " + ")
  )
  mods[[i]] <- workflow() %>% 
    add_recipe(build.rec(mod.data, response = diagnosis.col.name)) %>% 
    add_model(glm.spec, formula = as.formula(frms[i])) %>% 
    fit(mod.data) %>% 
    extract_fit_engine()
  summary(mods[[i]]) %>% print()
  Anova(mods[[i]]) %>% print()
  do.call(anova, c(mods, test = "Chisq")) %>% print()
  
  i <- 3
  frms[i] <- paste(
    diagnosis.col.name, "~", 
    "Chao1 +", batchID.col.name, "+", 
    paste(all.covars, collapse = " + "),
    "+", paste0("Chao1:", wkly.BMs.col.name)
  )
  mods[[i]] <- workflow() %>% 
    add_recipe(build.rec(mod.data, response = diagnosis.col.name)) %>% 
    add_model(glm.spec, formula = as.formula(frms[i])) %>% 
    fit(mod.data) %>% 
    extract_fit_engine()
  summary(mods[[i]]) %>% print()
  Anova(mods[[i]]) %>% print()
  do.call(anova, c(mods[2:3], test = "Chisq")) %>% print()
  
  i <- 4
  frms[i] <- paste(
    diagnosis.col.name, "~", 
    "Chao1 +", batchID.col.name, "+", 
    paste(all.covars, collapse = " + "),
    "+", paste0("Chao1:", age.col.name)
  )
  mods[[i]] <- workflow() %>% 
    add_recipe(build.rec(mod.data, response = diagnosis.col.name)) %>% 
    add_model(glm.spec, formula = as.formula(frms[i])) %>% 
    fit(mod.data) %>% 
    extract_fit_engine()
  summary(mods[[i]]) %>% print()
  Anova(mods[[i]]) %>% print()
  do.call(anova, c(mods[c(2, 4)], test = "Chisq")) %>% print()
  mod.stats <- new.env()
  mod.stats$brs.plots <- NULL
  mod.stats$hoslem.plots <- NULL
  lapply(seq_along(mods), function(i) {
    mod <- mods[[i]]
    frm <- frms[i]
    mod.stats$brs.plots <- binnedplot_tjfs(
      fitted(mod),
      residuals(mod, type = 'response'),
      main = paste0("Frm: ", str_wrap(frm, width = 50))
    ) %>% 
      list() %>% 
      c(mod.stats$brs.plots,  .)
    hoslem.t <- ResourceSelection::hoslem.test(
      as.integer(mod.data[[diagnosis.col.name]]) - 1, 
      fitted(mod)
    )
    mod.stats$hoslem.plots <- {
      data.table(
        obs_p = hoslem.t$observed[, 2], 
        exp_p = hoslem.t$expected[, 2]
      ) %>% 
        ggplot(aes(x = obs_p, y = exp_p)) +
        geom_point() + 
        geom_smooth() +
        geom_abline(intercept = 0, slope = 1, linewidth = 0.5) +
        labs(subtitle = paste("Frm:", str_wrap(frm, width = 50)))
    } %>% 
      list() %>% 
      c(mod.stats$hoslem.plots, .)
    chisq.dev <- pchisq(
      q = mod$null.deviance - mod$deviance, 
      df = length(coef(mod))
    ) %>% 
      subtract(1, .)
    p.rsq <- { mod$deviance / mod$null.deviance } %>% 
      subtract(1, .)
    data.table(
      Frm = frm,
      Deviance.pval = chisq.dev,
      Hoslem.pval = hoslem.t$p.value,
      Pseudo.Rsq = p.rsq,
      AIC = mod$aic 
    )
  }) %>% 
    rbindlist() %>%
    print()
  fig.file <- file.path(dirs$plots, "alpha_stool_glms_binnedResids.png")
  cowplot::plot_grid(plotlist = mod.stats$brs.plots) %>% 
    ggsave(filename = fig.file, width = 12, height = 10, dpi = 300)
  fig.file <- file.path(dirs$plots, "alpha_stool_glms_hoslemPlots.png")
  cowplot::plot_grid(plotlist = mod.stats$hoslem.plots) %>% 
    ggsave(filename = fig.file, width = 12, height = 10, dpi = 300)
  stool.mod <- mods[[3]]
  list(saliva = saliva.mod, stool = stool.mod)
}
###

alpha.cors.file <- file.path(
  dirs$saved,
  "taxa_alpha_cors.rds"
)
alpha.cors <- redo.if(
  redo.var,
  alpha.cors.file,
  max.save.mb = max.save.size,
  {
    lapply(sources, function(source) {
      lapply(focal.sets$taxon, function(ft) {
        track.values(source, ft)
        element <- list(source, "taxon", ft)
        mD <- pluck(mDs, !!!element)
        # ctrls <- pluck(ctrl.vars, !!!element)
        mod.data <- get.metadata(mD)
        mod.data <- mod.data[
          , .SD, .SDcols = c(diagnosis.col.name, focal.alphas, all.covars)
        ]
        lapply(1:{length(alphas) - 1}, function(i) {
          lapply({i + 1}:length(alphas), function(j) {
            mod <- lm(as.formula(paste(alphas[i], "~", alphas[j])), data = mod.data)
            pval <- tidy(mod)[2, 5]
            rsq <- summary(mod)$adj.r.squared
            p <- ggplot(mod.data, aes(x = !!sym(alphas[i]), y = !!sym(alphas[j]))) +
              geom_point() +
              labs(
                title = paste("Source:", source),
                subtitle = paste0("Adj. R-sq: ", round(rsq, 2), "; p-value: ", round(pval, 3))
              )
            if (pval <= 0.05) {
              p <- p + geom_smooth(method = "lm", formula = y ~ x)
            }
            return(p)
          })
        }) %>% unlist(recursive = F)
      }) %>% unlist(recursive = F)
    }) %>%
      unlist(recursive = F) %>%
      cowplot::plot_grid(plotlist = .)
  }
)
# 
nonpara.pd.alpha.res.file <- file.path(dirs$saved, "pd_alpha_nonparametric_dt.rds")
nonpara.pd.alpha.res <- redo.if(redo.var, nonpara.pd.alpha.res.file, {
  lapply(sources, function(source) {
    lapply(focal.sets$taxon, function(ft) {
      element <- list(source, "taxon", ft)
      mD <- pluck(mDs, !!!element)
      alpha <- ifelse(source == "saliva", "InvSimpson", "Chao1")
      track.values(source, ft, alpha)
      mod.data <- get.metadata(mD)
      frm <- paste(alpha, "~", diagnosis.col.name)
      wt <- wilcox.test(as.formula(frm), data = mod.data)
      data.table(
        Source = source,
        Feature = ft,
        Alpha = alpha,
        W = wt$statistic,
        P.value = wt$p.value,
        Sig = ifelse(wt$p.value <= 0.05, "*", "")
      ) %>% return()
      
    }) %>% rbindlist()
  }) %>% rbindlist()
})
nonpara.pd.alpha.res[, Adj.pval := p.adjust(P.value, method = "fdr")]

element <- list("stool", "taxon", "Taxon")
mD <- pluck(mDs, !!!element)
alpha <- "Chao1"
mod.data <- get.metadata(mD)

frm <- paste(alpha, "~", diagnosis.col.name, "*", wkly.BMs.col.name)
mod <- aov(
  as.formula(frm), 
  data = mod.data[get(wkly.BMs.col.name) != "3-5"]
)
summary(mod)
Anova(mod)
TukeyHSD(mod, which = paste0(diagnosis.col.name, ":", wkly.BMs.col.name))

dts <- lapply(sources, function(source) {
  dt <- get.metadata(mDs[[source]]$taxon$Taxon)
  new.alphas <- paste(source, focal.alphas, sep = ".")
  for (alpha in focal.alphas) {
    dt[, (str_subset(new.alphas, alpha)) := get(alpha)]
  }
  dt[, .SD, .SDcols = c(internalID.col.name, diagnosis.col.name, new.alphas)] %>%
    return()
})
cmb.dts <-merge(dts[[1]], dts[[2]], by  = c(internalID.col.name, diagnosis.col.name)) %>%
  return()

ft <- "Taxon"
med.vars <- c(
  # PD meds variable names
)
med.grps <- c(
  # PD meds group names
)
med.dict <- data.table(
  Name = med.vars,
  Group = med.grps
)

pd.meds.res <- lapply(sources, function(source) {
  element <- list(source, "taxon", ft)
  mD <- pluck(mDs, !!!element)
  keep.vars <- c(internalID.col.name, med.vars, alphas)
  meds.dt <- get.metadata(mD)[, ..keep.vars] %>% 
    .[complete.cases(.)]
  for (med.grp in med.grps) {
    grp.vars <- med.dict[Name %in% med.vars & Group == med.grp]$Name
    meds.dt[[make.names(med.grp)]] <- rowSums(meds.dt[, ..grp.vars])
  }
  med.grps %<>% make.names()
  meds.dt$N.types.meds <- rowSums(meds.dt[, .SD > 0, .SDcols = med.grps])
  
  p1 <- melt(
    meds.dt, 
    measure.vars = med.grps, 
    variable.name = "Med.type", 
    value.name = "Taking"
  ) %>% 
    melt(measure.vars = alphas, variable.name = "Alpha", value.name = "Score") %>% 
    ggplot(aes(x = Taking, y = Score)) + 
    geom_quasirandom() +
    stat_summary(fun.data = "mean_se", geom = "errorbar", color = "blue") +
    facet_grid(Alpha ~ Med.type, scales = "free") +
    labs(subtitle = source)
  p2 <- melt(
    meds.dt, 
    measure.vars = alphas, 
    variable.name = "Alpha", 
    value.name = "Score"
  ) %>% 
    ggplot(aes(x = N.types.meds, y = Score)) + 
    geom_quasirandom() +
    stat_summary(fun.data = "mean_se", geom = "errorbar", color = "blue") +
    facet_wrap(~ Alpha, scales = "free") +
    labs(subtitle = source)
  kw.dt <- lapply(alphas, function(alpha) {
    grps.dt <- lapply(med.grps, function(med.grp) {
      frm <- paste(alpha, "~", med.grp)
      kruskal.test(as.formula(frm), data = meds.dt) %>% 
        broom::tidy() %>% 
        as.data.table() %>% 
        .[,`:=`(Alpha = alpha, Med.group = med.grp, Source = source)]
    }) %>% rbindlist()
    frm <- paste(alpha, "~ N.types.meds")
    kruskal.test(as.formula(frm), data = meds.dt) %>% 
      broom::tidy() %>% 
      as.data.table() %>% 
      .[,`:=`(Alpha = alpha, Med.group = "N.types.meds", Source = source)] %>% 
      rbind(grps.dt, .) %>% 
      return()
  }) %>% rbindlist()
  list(Plots = list(p1, p2), Mods = kw.dt) %>% 
    return()
})
pd.meds.kws <- lapply(sources, function(source) {
  pd.meds.res[[source]]$Mods
}) %>% rbindlist()
pd.meds.kws[, `:=`(parameter = NULL, method = NULL)]
setcolorder(pd.meds.kws, c(5, 3, 4))
kable(pd.meds.kws, digits = c(0, 0, 0, 2, 3)) %>% 
  kable_styling(full_width = F) %>% 
  collapse_rows(columns = 1:2)


ft <- "Taxon"
batches.dt <- read.table(
  "", # batches file 
  sep = "\t", 
  header = TRUE
) %>% as.data.table()
batches.res <- lapply(sources, function(source) {
  require(effectsize)
  element <- list(source, "taxon", ft)
  mD <- pluck(mDs, !!!element)
  metadata <- get.metadata(mD) %>% 
    merge(batches.dt, by = externalID.col.name, all.x = TRUE)
  lapply(alphas, function(alpha) {
    frm <- paste(alpha, "~", batchID.col.name)
    mod.data <- metadata[, .SD, .SDcols = c(alpha, diagnosis.col.name, batchID.col.name)] %>% 
      .[complete.cases(.)]
    mod.data[, .N, by = c(diagnosis.col.name, batchID.col.name)] %>% 
      dcast(as.formula(paste(batchID.col.name, "~", diagnosis.col.name)))
    lm0 <- lm(as.formula(frm), data = mod.data)
    lm1 <- update(lm0, as.formula(paste(". ~ . *", diagnosis.col.name)))
    effectsize::eta_squared(lm1, alternative = "two")
    mods <- list(
      mod0.lm = summary(lm0),
      mod0.aov = Anova(lm0),
      mod1.lm = summary(lm1),
      mod1.aov = Anova(lm1),
      comp.aov = anova(lm0, lm1)
    )
    sub.title <- paste(
      "PD x Batch is",
      ifelse(tidy(mods$comp.aov)$p.value[2] < 0.05, "significant", "non-significant")
    )

    alpha.p <- ggplot(mod.data, aes(x = !!sym(batchID.col.name), y = !!sym(alpha))) +
      geom_quasirandom() +
      stat_summary(fun.data = "mean_se", geom = "errorbar") +
      labs(title = source, subtitle = sub.title)
    list(Mod.res = mods, Plot = alpha.p) %>% 
      return()
  })
})

lapply(batches.res, function(x) { lapply(x, function(y) {y$Plot})}) %>% 
  unlist(recursive = F) %>% 
  cowplot::plot_grid(plotlist = .)

library(maps)
pd.zip.file <- ""
pd.zip.dt <- readRDS(pd.zip.file)
county.data.file <- ""
states.file <- ""

county.data <- read.csv(county.data.file) %>% as.data.table()
regions.dt <- data.table(
  STATE = c(
    "AK", "HI", "WA", "OR", "CA", "NV", "ID", "MT", "WY", "CO", "UT", "AZ", "NM",
    "ND", "SD", "NE", "KS", "MO", "IA", "MN", "WI", "IL", "IN", "MI", "OH",
    "TX", "OK", "AR", "LA", "MS", "AL", "TN", "KY", "FL", "GA", "SC", "NC", "VA", "WV", "MD", "DE", "DC",
    "PA", "NJ", "CT", "RI", "NY", "MA", "VT", "NH", "ME"
  ),
  Region = c(
    rep("West", 13),
    rep("Midwest", 12),
    rep("South", 17),
    rep("Northeast", 9)
  )
)


pd.region.dt <- merge(
  pd.zip.dt,
  county.data,
  by = "zip.code",
  all.x = T
) %>% 
  merge(regions.dt, by = "state") %>% 
  unique()
saveRDS(pd.region.dt, file.path(dirs$saved, "US_regions_dt.rds"))
keep.cols <- c(internalID.col.name, diagnosis.col.name, "Region", alphas)
lapply(sources, function(source) {
  mod.data <- merge(
    get.metadata(mDs[[source]]$taxon$Taxon), 
    pd.region.dt,
    by = c(internalID.col.name, diagnosis.col.name),
  ) %>% 
    .[, ..keep.cols] %>% 
    .[complete.cases(.)] %>% 
    unique()
  lapply(alphas, function(alpha) {
    track.values(source, alpha)
    frm0 <- paste(alpha, "~", diagnosis.col.name)
    mod0 <- lm(as.formula(frm0), data = mod.data)
    summary(mod0) %>% print()
    mod1 <- update(mod0, . ~ . + Region)
    summary(mod1) %>% print()
    Anova(mod1) %>% print()
    mod2 <- update(mod0, . ~ . * Region)
    summary(mod2) %>% print()
    Anova(mod2) %>% print()
  })
})
