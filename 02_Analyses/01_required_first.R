# candidate_control_covariates.R
base.dir <- ""
if (getwd() != base.dir) { setwd(base.dir) }
if (!exists("color.scales")) { source("_setup.R") }
redo.var <- "req.first"
if (is.null(get.redo.state(redo.var))) {
  assign.redo(redo.var, state = T)
}

raw.split.mDs <- list(
  saliva = list(
    taxon = list(
      Taxon = "" # saliva OTU raw input mD
    ),
    func = list(
      KO = "" # saliva KO raw input mD
    )
  ),
  stool = list(
    taxon = list(
      Taxon = "" # stool OTU raw input mD
    ),
    func = list(
      KO = "" # stool KO raw input mD
    )
  )
)

pa.filt.mDs.file <- file.path(
  dirs$input, 
  paste0(
    "prevCut", prev.cut,
    "_abundCut", abund.cut, 
    "_mDs_list.rds"
  )
)
update.files("pa_filt_mDs", pa.filt.mDs.file)
pa.filt.mDs <- redo.if(
  redo.var, 
  pa.filt.mDs.file, 
  max.save.mb = max.save.size, 
  {
    lapply(sources, function(source) {
      lapply(focal.sets, function(set) {
        lapply(set, function(ft) {
          set.name <- get.set.name(ft)
          track.values(abund.cut, source, ft)
          
          mD <- purrr::pluck(raw.split.mDs, source, set.name, ft)
          if (set.name == "taxon") {
            mD <- filter.features(mD, Kingdom != "NA_Kingdom")
            if (is.null(attributes(mD@Assignments)$sorted)) {
              setkeyv(mD@Assignments, mD@Feature.col)
            }
          }
          
          feat.relab <- feature.relAbunds(mD)
          mD <- keep.features(
            mD,
            feature = names(feat.relab[feat.relab >= abund.cut])
          )
          feat.prev <- feature.prevalences(mD)
          mD <- keep.features(
            mD,
            features = names(feat.prev[feat.prev >= prev.cut])
          ) %>% key.mD.tbls()
          
          rlang::inform("generating rarefaction curve...")
          fig.file <- file.path(
            dirs$plots,
            paste("rarefaction_curve", source, ft, sep = "_") %>% 
              paste0(".png")
          )
          sp <- get.abundances(mD) %>% 
            vegan::specaccum(method = "random")
          
          png(fig.file)
          plot(
            sp, 
            ylab = plurals[ft]$Plural, 
            main = paste0("Abund. cut = ", abund.cut, "; ", source)
          )
          dev.off()
          
          return(mD)
        })
      })
    })
  }
) %>% filter.list(condition = save.space, keep.names = focal.sets)

rar.pa.mDs.file <- str_replace(
  pa.filt.mDs.file, 
  "prevCut", 
  "rarefied_prevCut"
)
update.files("rarefied_filt_mDs", rar.pa.mDs.file)
rar.mDs <- redo.if(redo.var, rar.pa.mDs.file, max.save.mb = max.save.size, {
  lapply(sources, function(source) {
    lapply(focal.sets, function(set) {
      lapply(set, function(ft) {
        set.name <- get.set.name(ft)
        track.values(source, ft)
        pluck(pa.filt.mDs, source, set.name, ft) %>% 
          rarefy(user.seed = seed) %>% 
          alpha.diversity(metrics = alphas) %>% 
          beta.diversity(metrics = betas) %>% 
          key.mD.tbls() %>% 
          return()
      })
    })
  })
})

raw.combined.mDs <- list(
  taxon = list(
    Taxon = "" # combined sources OTU raw input mD
  ),
  func = list(
    KO = "" # combined sources KO raw input mD
  )
)
cmb.mDs.file <- file.path(dirs$saved, "combined_mDs_glommed_paFilt_rar.rds")
update.files("cmb_mDs", cmb.mDs.file)
cmb.mDs <- redo.if(redo.var, cmb.mDs.file, max.save.mb = max.save.size, {
  lapply(focal.sets, function(set) {
    lapply(set, function(ft) {
      set.name <- get.set.name(ft)
      track.values(ft)
      
      mD <- pluck(raw.combined.mDs, set.name, ft)
      if (set.name == "taxon") {
        mD <- filter.features(mD, Kingdom != "NA_Kingdom")
        if (is.null(attributes(mD@Assignments)$sorted)) {
          setkeyv(mD@Assignments, mD@Feature.col)
        }
      }
      
      feat.relab <- feature.relAbunds(mD)
      mD <- keep.features(
        mD,
        feature = names(feat.relab[feat.relab >= abund.cut])
      )
      feat.prev <- feature.prevalences(mD)
      mD <- keep.features(
        mD,
        features = names(feat.prev[feat.prev >= prev.cut])
      ) %>% key.mD.tbls()
      
      rlang::inform("generating rarefaction curve...")
      fig.file <- file.path(
        dirs$plots,
        paste("rarefaction_curve_combined", ft, sep = "_") %>% 
          paste0(".png")
      )
      sp <- get.abundances(mD) %>% 
        specaccum(method = "random")
      
      png(fig.file)
      plot(
        sp, 
        ylab = plurals[ft]$Plural, 
        main = paste0("Abund. cut = ", abund.cut, "; ", source)
      )
      dev.off()
      rarefy(mD, user.seed = seed) %>% 
        alpha.diversity(metrics = alphas) %>% 
        beta.diversity(metrics = betas) %>% 
        key.mD.tbls() %>% 
        return()
    })
  })
})


# TEST MAIN EFFECTS

if (get(redo.var, envir = redo)) {
  candidates <- list(
    stool = c(
      # selected covariates column names
    ),
    saliva = c(
      # selected covariates column names
    )
  )
  
  cat("Generating models...", sep = "\n")
  sources.res <- lapply(sources, function(source) {
    covars <- candidates[[source]] %>% set_names()
    covars.res <- lapply(covars, function(covar) {
      set.res <- lapply(focal.sets, function(set) {
        ft.res <- lapply(set, function(ft) {
          track.values(source, covar, ft)
          set.name <- get.set.name(ft)
          mD <- pluck(rar.mDs, source, set.name, ft)
          if (set.name == "func") {
            alpha.mods <- NULL
          } else {
            alpha.mods <- lapply(alphas, function(alpha) {
              track.values(alpha, prepend = "\t")
              frm <- paste(alpha, "~", covar)
              lm(
                as.formula(frm), 
                data = mD@Metadata, 
                na.action = na.omit
              ) %>% 
                Anova() %>% 
                tidy() %>% 
                as.data.table() %>% 
                set_names(tools::toTitleCase(names(.))) %>% 
                .[Term == covar] %>% 
                .[
                  , `:=`(
                    Sig = ifelse(P.value <= 0.05, "*", ""),
                    Source = source,
                    Feature = ft,
                    Alpha = alpha
                  )
                ] %>% 
                setcolorder(c(1, tail(1:ncol(.), 3))) %>%
                return()
            }) %>% rbindlist()
          }
          frm <- paste("~", covar)
          ords <- microbData::ordinate(
            mD, 
            formula = as.formula(frm), 
            update.mD = FALSE, 
            na.action = na.omit
          )
          beta.mods <- lapply(betas, function(beta) {
            track.values(beta, prepend = "\t")
            ord <- ords[[str_subset(names(ords), beta)]]
            anova(ord, by = "terms") %>% 
              as.data.table(keep.rownames = T) %>% 
              set_names(
                c("Term", "DF", "Sum.of.sqs", "F", "P.value")
              ) %>% 
              .[Term == covar] %>% 
              .[
                , `:=`(
                  Sig = ifelse(P.value <= 0.05, "*", ""),
                  Source = source,
                  Feature = ft,
                  Beta = beta
                )
              ] %>% 
              setcolorder(c(1, tail(1:ncol(.), 3))) %>% 
              return()
          }) %>% rbindlist()
          
          list(Alphas = alpha.mods, Betas = beta.mods) %>% 
            return()
        })
        list(
          Alphas = lapply(ft.res, function(f) { f$Alphas }) %>% 
            rbindlist(),
          Betas = lapply(ft.res, function(f) { f$Betas }) %>% 
            rbindlist()
        ) %>% return()
      })
      list(
        Alphas = lapply(set.res, function(f) { f$Alphas }) %>% 
          rbindlist(),
        Betas = lapply(set.res, function(f) { f$Betas }) %>% 
          rbindlist()
      ) %>% return()
    })
    list(
      Alphas = lapply(covars.res, function(c) { c$Alphas }) %>%
        rbindlist(),
      Betas = lapply(covars.res, function(c) { c$Betas }) %>% 
        rbindlist()
    ) %>% return()
  })
  models.res <- list(
    Alphas = lapply(sources.res, function(s) { s$Alphas }) %>% 
      rbindlist(),
    Betas = lapply(sources.res, function(s) { s$Betas }) %>% 
      rbindlist()
  ) 
  res.tbl0 <- merge(
    dcast(
      data = models.res$Alphas,
      formula = Source + Term ~ Feature + Alpha,
      value.var = "Sig"
    ),
    dcast(
      data = models.res$Betas,
      formula = Source + Term ~ Feature + Beta,
      value.var = "Sig"
    ),
    by = c("Source", "Term")
  )
  res.tbl <- res.tbl0[, lapply(.SD, function(x) { ifelse(is.na(x), "", x) })]
  
  write.table(
    res.tbl,
    file = file.path(dirs$plots, "candidate_ctrl_covars.tsv"),
    sep = "\t",
    quote = TRUE,
    row.names = FALSE
  )
}

# SELECT COVARIATES
if (get(redo.var, envir = redo)) {
  ctrl.vars.file <- file.path(dirs$saved, "control_covariates.rds")
  update.files("ctrl_vars", ctrl.vars.file)
  ctrl.vars <- lapply(sources, function(source) {
    s.tbl <- res.tbl[Source == source]
    lapply(focal.sets, function(set) {
      lapply(set, function(ft) {
        track.values(source, ft)
        cols <- str_subset(names(s.tbl), paste0("^", ft, "_"))
        if (ft %in% feature.sets$func) {
          s.tbl[, ..cols] %>% 
            apply(MARGIN = 1, function(r) any(str_detect(r, "\\*"))) %>% 
            s.tbl$Term[.] %>% 
            return()
        } else {
          a.cols <- cols[1:3]
          b.cols <- cols[4:6]
          c(
            {
              s.tbl[, ..a.cols] %>% 
                apply(
                  MARGIN = 1, 
                  function(r) any(str_detect(r, "\\*"))
                ) %>% 
                s.tbl$Term[.]
            },
            {
              s.tbl[, ..b.cols] %>% 
                apply(
                  MARGIN = 1, 
                  function(r) any(str_detect(r, "\\*"))
                ) %>% 
                s.tbl$Term[.]
            }
          ) %>% 
            table() %>% 
            .[. > 1] %>% 
            names() %>% 
            return()
        }
      })
    })
  })
  saveRDS(ctrl.vars, file = ctrl.vars.file)
} else {
  ctrl.vars <- readRDS(ctrl.vars.file)
}

# PLOTS
if (get(redo.var, envir = redo)) { 
  wd <- 10
  ht <- 4
  dpi <- 300
  to.save <- lapply(sources, function(source) {
    lapply(feature.sets, function(set) {
      lapply(set, function(ft) {
        set.name <- get.set.name(ft)
        element <-  list(source, set.name, ft)
        lapply(pluck(ctrl.vars, !!!element), function(var) {
          track.values(source, ft, var)
          sig.alphas <- models.res$Alphas[
            Source == source &
              Feature == ft &
              Term == var &
              Sig == "*"
          ]$Alpha
          if (length(sig.alphas) > 0) {
            alpha.dt <- melt(
              get.metadata(pluck(rar.mDs, !!!element)),
              measure.vars = sig.alphas, 
              id.var = c(diagnosis.col.name, var),
              value.name = "Score",
              variable.name = "Alpha"
            )
            alpha.dt <- alpha.dt[complete.cases(alpha.dt)]
            if (class(alpha.dt[[var]]) %in% c("character", "factor")) {
              base.plot <- ggplot(
                alpha.dt, aes(x = !!sym(var), y = Score)
              ) +
                facet_wrap(~ Alpha, scales = "free_y")
              alpha.plots <- list(
                { # without PD
                  base.plot +
                    geom_quasirandom() + 
                    stat_summary(
                      fun.data = "mean_cl_boot",
                      geom = "errorbar",
                      color = "red"
                    )
                },
                { # with PD
                  base.plot + 
                    geom_quasirandom(
                      aes(color = !!sym(diagnosis.col.name)), 
                      dodge.width = 1
                    ) + 
                    stat_summary(
                      aes(group = !!sym(diagnosis.col.name)),
                      fun.data = "mean_cl_boot",
                      position = position_dodge(width = 1),
                      geom = "errorbar",
                      color = "black"
                    ) + 
                    color.scales$color$Diagnosis
                }
              ) %>% 
                cowplot::plot_grid(
                  plotlist = .,
                  nrow = ifelse(length(sig.alphas) > 1, 2, 1),
                  align = "hv",
                  axis = "tbl"
                )
            } else {
              base.plot <- ggplot(
                alpha.dt, aes(x = !!sym(var), y = Score)
              ) +
                facet_wrap(~ Alpha, scales = "free_y")
              alpha.plots <- list(
                { # without PD
                  base.plot +
                    geom_point() +
                    geom_smooth(
                      method = "lm", 
                      formula = y ~ x, 
                      color = "red"
                    )
                },
                { # with PD
                  base.plot +
                    geom_point(aes(color = !!sym(diagnosis.col.name))) +
                    geom_smooth(
                      aes(color = !!sym(diagnosis.col.name)), 
                      method = "lm", 
                      formula = y ~ x
                    ) +
                    color.scales$color$Diagnosis
                }
              ) %>% 
                cowplot::plot_grid(
                  plotlist = .,
                  nrow = ifelse(length(sig.alphas) > 1, 2, 1),
                  align = "hv",
                  axis = "tbl"
                )
            }
          } else {
            alpha.plots <- NULL
          }
          sig.betas <- models.res$Betas[
            Source == source &
              Feature == ft &
              Term == var &
              Sig == "*"
          ]$Beta
          if (length(sig.betas) > 0) {
            mD.noPD <- microbData::ordinate(
              mD = pluck(rar.mDs, !!!element), 
              formula = as.formula(paste("~", var)),
              na.action = na.omit
            )
            coords.noPD <- ordination.coords(
              mD.noPD, 
              axis.digits = 3
            ) %>% 
              lapply(function(x) return(x[Beta.metric %in% sig.betas]))
            mD.yesPD <- microbData::ordinate(
              mD = pluck(rar.mDs, !!!element),
              formula = as.formula(paste("~", diagnosis.col.name, "*", var)),
              na.action = na.omit
            )
            coords.yesPD <- ordination.coords(
              mD.yesPD, 
              axis.digits = 3
            ) %>% 
              lapply(function(x) return(x[Beta.metric %in% sig.betas]))
            beta.plots <- list(
              { # without PD
                bp <- ggplot(
                  coords.noPD$Samples, aes(x = Axis1, y = Axis2)
                ) + 
                  geom_point(aes(color = !!sym(var))) + 
                  geom_text(
                    data = coords.noPD$Axis.labs,
                    aes(label = Label, angle = Angle, vjust = Vjust),
                    hjust = 1
                  ) +
                  facet_wrap(~ Beta.metric, scales = "free") +
                  theme(axis.title = element_blank())
                if (
                  class(coords.yesPD$Samples[[var]]) %in% 
                  c("numeric", "integer")
                ) {
                  bp
                } else {
                  bp + stat_ellipse(aes(color = !!sym(var)))
                }
              },
              { # with PD
                if (
                  class(coords.yesPD$Samples[[var]]) %in% 
                  c("numeric", "integer")
                ) {
                  plot.var <- paste0("Binned.", var)
                  quants <- quantile(
                    coords.yesPD$Samples[[var]], 
                    c(1/3, 2/3, 1),
                    na.rm = T
                  ) %>% 
                    ceiling()
                  coords.yesPD$Samples[
                    , (plot.var) := case_when(
                      .SD < quants[1] ~ paste("<", quants[1]),
                      .SD < quants[2] ~ paste("<", quants[2]),
                      .SD < quants[3] ~ paste("<", quants[3])
                    ),
                    .SDcols = var
                  ]
                  
                } else {
                  plot.var <- var
                }
                ggplot(coords.yesPD$Samples, aes(x = Axis1, y = Axis2)) + 
                  geom_point(
                    aes(color = !!sym(diagnosis.col.name), shape = !!sym(plot.var))
                  ) +
                  stat_ellipse(
                    aes(color = !!sym(diagnosis.col.name), linetype = !!sym(plot.var))
                  ) +
                  color.scales$color$Diagnosis + 
                  geom_text(
                    data = coords.yesPD$Axis.labs,
                    aes(label = Label, angle = Angle, vjust = Vjust),
                    hjust = 1
                  ) +
                  facet_wrap(~ Beta.metric, scales = "free") +
                  theme(axis.title = element_blank())
              }
            ) %>% 
              cowplot::plot_grid(
                plotlist = .,
                nrow = ifelse(length(sig.betas) > 1, 2, 1),
                align = "hv",
                axis = "tbl"
              )
          } else {
            beta.plots <- NULL
          }
          plot.file <- file.path(
            dirs$plots, 
            paste0(
              paste("candidate_ctrl", var, source, ft, sep = "_"), 
              ".pdf"
            )
          )
          if (!is.null(alpha.plots) & is.null(beta.plots)) {
            ggsave(
              alpha.plots,
              filename = plot.file,
              width = wd,
              height = ht * 2,
              dpi = dpi
            )
          } else if (is.null(alpha.plots) & !is.null(beta.plots)) {
            ggsave(
              beta.plots,
              filename = plot.file,
              width = wd,
              height = ht * 2,
              dpi = dpi
            )
          } else if (!is.null(alpha.plots) & !is.null(beta.plots)) {
            scales <- c(
              ifelse(length(sig.alphas) > 1, 2, 1.2), 
              ifelse(length(sig.betas) > 1, 2, 1.2)
            )
            cowplot::plot_grid(
              alpha.plots, beta.plots,
              ncol = 1,
              align = "v",
              axis = "l",
              rel_heights = scales
            ) %>% 
              ggsave(
                filename = plot.file,
                width = wd,
                height = sum(ht * scales),
                dpi = dpi
              )
          }
        })
      })
    })
  })
}

covars <- c(
  # Random forest covariates list
)
saveRDS(covars, file = file.path(dirs$saved, "rf_covars_vector.rds"))

set.redo.true(redo.var)

taxon.rf.data.file <- file.path(
  dirs$saved, 
  "taxon_rf_data_list.rds"
)
taxon.rf.data <- redo.if(
  redo.var, 
  taxon.rf.data.file, 
  max.save.mb = max.save.size, 
  {
    lapply(sources, function(source) {
      lapply(taxon.lvls, function(lvl) {
        track.values(source, lvl)
        set.seed(0)
        # covars <- unlist(ctrl.vars[[source]]) %>% unique()
        element <- list(source, "taxon", lvl)
        raw.mD <- pluck(raw.split.mDs, !!!element)
        metadata <- get.metadata(raw.mD)
        split <- rsample::initial_split(
          metadata, 
          strata = diagnosis.col.name, 
          prop = 0.7
        )
        train.meta <- training(split)
        train.mD <- keep.samples(raw.mD, samples = train.meta[[raw.mD@Sample.col]])
        feat.relab <- feature.relAbunds(train.mD)
        train.mD <- keep.features(
          train.mD,
          features = names(feat.relab[feat.relab >= abund.cut])
        )
        feat.prev <- feature.prevalences(train.mD)
        to.keep <- names(feat.prev[feat.prev >= prev.cut])
        
        clr.mD <- keep.features(raw.mD, features = to.keep) %>% 
          key.mD.tbls() %>% 
          center.log.ratio()
        train.abunds <- keep.samples(
          clr.mD, 
          samples = train.meta[[raw.mD@Sample.col]]
        ) %>% 
          get.abundances(as.DT = TRUE)
        keep.cols <- c(raw.mD@Sample.col, smplID.col.name, diagnosis.col.name, covars)
        train.data <- merge(
          train.meta[, .SD, .SDcols = keep.cols],
          train.abunds
        ) %>% 
          set_names(make.names(names(.))) %>% 
          .[, (raw.mD@Sample.col) := NULL]
        
        
        test.meta <- testing(split)
        test.abunds <- keep.samples(
          clr.mD, 
          samples = test.meta[[raw.mD@Sample.col]]
        ) %>% 
          get.abundances(as.DT = TRUE)
        test.data <- merge(
          test.meta[, .SD, .SDcols = keep.cols],
          test.abunds
        ) %>% 
          set_names(make.names(names(.))) %>% 
          .[, (raw.mD@Sample.col) := NULL]
        
        list(Train = train.data, Test = test.data) %>% 
          return()
      }) %>% return()
    })
  }
)

taxon.cmb.rf.data.file <- file.path(
  dirs$saved, 
  "taxon_combined_rf_data_list.rds"
)
taxon.cmb.rf.data <- redo.if(
  redo.var, 
  taxon.cmb.rf.data.file, 
  max.save.mb = max.save.size, 
  {
    lapply(cmb.methods, function(cmb.m) {
      lapply(taxon.lvls, function(lvl) {
        track.values(cmb.m, lvl)
        set.seed(0)
        # covars <- unlist(ctrl.vars) %>% unique()
        
        saliva.mD <- pluck(raw.split.mDs, "saliva", "taxon", lvl)
        saliva.mD@Sample.col <- smplID.col.name
        saliva.mD@Sample.names <- saliva.mD@Metadata[[saliva.mD@Sample.col]]
        rownames(saliva.mD@Abundances) <- saliva.mD@Sample.names
        stool.mD <- raw.split.mDs$stool$taxon[[lvl]]
        stool.mD@Sample.col <- smplID.col.name
        stool.mD@Sample.names <- stool.mD@Metadata[[stool.mD@Sample.col]]
        rownames(stool.mD@Abundances) <- stool.mD@Sample.names
        saliva.meta <- get.metadata(saliva.mD)
        stool.meta <- get.metadata(stool.mD)
        saliva.assign <- get.assignments(saliva.mD)
        stool.assign <- get.assignments(stool.mD)
        
        cmb.meta <- stool.meta[get(smplID.col.name) %in% saliva.meta[[smplID.col.name]]]
        cmb.meta[, sample_id := NULL]
        setkeyv(cmb.meta, smplID.col.name)
        
        cmb.abunds <- switch(
          cmb.m,
          split = {
            saliva.abunds <- get.abundances(saliva.mD)
            stool.abunds <- get.abundances(stool.mD)
            colnames(saliva.abunds) %<>% paste0("saliva__", .)
            colnames(stool.abunds) %<>% paste0("stool__", .)
            cmb.df <- merge(stool.abunds, saliva.abunds, by = 0, all = T)
            cmb.df <- cmb.df[cmb.df$Row.names %in% cmb.meta[[smplID.col.name]], ]
            cmb.abunds <- as.matrix(cmb.df[, -1])
            rownames(cmb.abunds) <- cmb.df$Row.names
            cmb.abunds[cmb.meta[[smplID.col.name]], ]
          },
          lumped = {
            saliva.abunds <- get.abundances(saliva.mD, as.DT = T) %>% 
              setkeyv(smplID.col.name)
            stool.abunds <- get.abundances(stool.mD, as.DT  = T) %>% 
              setkeyv(smplID.col.name)
            saliva.fts <- names(saliva.abunds)[-1]
            stool.fts <- names(stool.abunds)[-1]
            cmb.fts <- unique(c(stool.fts, saliva.fts))
            cmb.dt <- data.table(1)[
              (smplID.col.name) := cmb.meta[[smplID.col.name]]
            ]
            setkeyv(cmb.dt, smplID.col.name)
            for (ft in cmb.fts) {
              keep.cols <- c(smplID.col.name, ft)
              if (ft %in% saliva.fts & ft %in% stool.fts) {
                ft.dt <- merge(
                  stool.abunds[, ..keep.cols], 
                  saliva.abunds[, ..keep.cols]
                )
                ft.dt[is.na(ft.dt)] <- 0
                ft.dt[, (ft) := get(paste0(ft, ".x")) + get(paste0(ft, ".y"))]
                ft.dt <- ft.dt[, ..keep.cols]
              } else if (ft %in% saliva.fts) {
                ft.dt <- saliva.abunds[cmb.meta[[smplID.col.name]], ..keep.cols]
                ft.dt[is.na(ft.dt)] <- 0
              } else {
                ft.dt <- stool.abunds[cmb.meta[[smplID.col.name]], ..keep.cols]
                ft.dt[is.na(ft.dt)] <- 0
              }
              cmb.dt <- cmb.dt[ft.dt, on = smplID.col.name]
            }
            cmb.abunds <- as.matrix(cmb.dt[, -1])
            rownames(cmb.abunds) <- cmb.dt[[smplID.col.name]]
            cmb.abunds[cmb.meta[[smplID.col.name]], ]
          }
        )
        cmb.assign <- switch(
          cmb.m,
          split = {
            saliva.assign[[stool.mD@Feature.col]] %<>% paste0("saliva__", .)
            stool.assign[[stool.mD@Feature.col]] %<>% paste0("stool__", .)
            assign.dt <- rbind(saliva.assign, stool.assign) %>% unique()
            assign.dt[get(stool.mD@Feature.col) %in% colnames(cmb.abunds)]
          },
          lumped = {
            assign.dt <- rbind(saliva.assign, stool.assign) %>% unique()
            assign.dt[get(stool.mD@Feature.col) %in% colnames(cmb.abunds)]
          }
        )
        setkeyv(cmb.assign, stool.mD@Feature.col)
        cmb.mD <- microbData(
          metadata = cmb.meta,
          abundances = cmb.abunds,
          assignments = cmb.assign
        )
        
        split <- rsample::initial_split(
          cmb.meta, 
          strata = diagnosis.col.name, 
          prop = 0.7
        )
        train.meta <- training(split)
        train.mD <- keep.samples(cmb.mD, samples = train.meta[[cmb.mD@Sample.col]])
        
        feat.relab <- feature.relAbunds(train.mD)
        train.mD <- keep.features(
          train.mD,
          features = names(feat.relab[feat.relab >= abund.cut])
        )
        feat.prev <- feature.prevalences(train.mD)
        to.keep <- names(feat.prev[feat.prev >= prev.cut])
        
        clr.mD <- keep.features(cmb.mD, features = to.keep) %>% 
          key.mD.tbls() %>% 
          center.log.ratio()
        train.abunds <- keep.samples(
          clr.mD, 
          samples = train.meta[[cmb.mD@Sample.col]]
        ) %>% 
          get.abundances(as.DT = TRUE)
        keep.cols <- c(cmb.mD@Sample.col, diagnosis.col.name, covars)
        train.data <- merge(
          train.meta[, .SD, .SDcols = keep.cols],
          train.abunds
        ) %>% 
          set_names(make.names(names(.)))
        
        
        test.meta <- testing(split)
        test.abunds <- keep.samples(
          clr.mD, 
          samples = test.meta[[cmb.mD@Sample.col]]
        ) %>% 
          get.abundances(as.DT = TRUE)
        test.data <- merge(
          test.meta[, .SD, .SDcols = keep.cols],
          test.abunds
        ) %>% 
          set_names(make.names(names(.)))
        
        list(Train = train.data, Test = test.data) %>% 
          return()
      }) %>% return()
    })
  }
)

func.rf.data.file <- file.path(dirs$saved, "func_rf_data_list.rds")
func.rf.data <- redo.if(
  redo.var, 
  func.rf.data.file, 
  max.save.mb = max.save.size, 
  {
    lapply(sources, function(source) {
      lapply(func.lvls, function(lvl) {
        track.values(source, lvl)
        set.seed(0)
        covars <- unlist( <- [[source]]) %>% unique()
        
        raw.mD <- pluck(raw.split.mDs, source, "func", lvl)
        metadata <- get.metadata(raw.mD)
        split <- rsample::initial_split(
          metadata, 
          strata = diagnosis.col.name, 
          prop = 0.7
        )
        train.meta <- training(split)
        train.mD <- keep.samples(raw.mD, samples = train.meta[[raw.mD@Sample.col]])
        feat.relab <- feature.relAbunds(train.mD)
        train.mD <- keep.features(
          train.mD,
          features = names(feat.relab[feat.relab >= abund.cut])
        )
        feat.prev <- feature.prevalences(train.mD)
        to.keep <- names(feat.prev[feat.prev >= prev.cut])
        
        clr.mD <- keep.features(raw.mD, features = to.keep) %>% 
          key.mD.tbls() %>% 
          center.log.ratio()
        train.abunds <- keep.samples(
          clr.mD, 
          samples = train.meta[[raw.mD@Sample.col]]
        ) %>% 
          get.abundances(as.DT = TRUE)
        train.data <- merge(
          train.meta[, .SD, .SDcols = c(raw.mD@Sample.col, diagnosis.col.name, covars)],
          train.abunds
        ) %>% 
          set_names(make.names(names(.))) %>% 
          .[, (raw.mD@Sample.col) := NULL]
        
        
        test.meta <- testing(split)
        test.abunds <- keep.samples(
          clr.mD, 
          samples = test.meta[[raw.mD@Sample.col]]
        ) %>% 
          get.abundances(as.DT = TRUE)
        test.data <- merge(
          test.meta[, .SD, .SDcols = c(raw.mD@Sample.col, diagnosis.col.name, covars)],
          test.abunds
        ) %>% 
          set_names(make.names(names(.))) %>% 
          .[, (raw.mD@Sample.col) := NULL]
        list(Train = train.data, Test = test.data) %>% 
          return()
        
      }) %>% return()
    })
  }
)
