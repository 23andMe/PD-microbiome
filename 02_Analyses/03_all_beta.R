# all_beta.R

base.dir <- ""
if (getwd() != base.dir) { setwd(base.dir) }
if (!exists("color.scales")) { source("_setup.R") }
redo.var <- "beta"
if (is.null(get.redo.state(redo.var))) {
  assign.redo(redo.var, state = T)
}


mDs <- readRDS(files$rarefied_filt_mDs) # generated in 01_required_first.R
all.covars <- unlist(ctrl.vars) %>% 
  unique() %>% 
  sort()

ords.file <- file.path(dirs$saved, "best_ordinations.rds")
ords <- redo.if(redo.var, ords.file, max.save.mb = max.save.size, {
  lapply(sources, function(source) {
    lapply(focal.sets, function(set) {
      lapply(set, function(ft) {
        set.name <- get.set.name(ft)
        track.values(source, ft)
        element <- list(source, set.name, ft)
        mD <- pluck(mDs, !!!element)
        # all.covars <- pluck(ctrl.vars, !!!element)
        all.terms <- c(mD@Sample.col, diagnosis.col.name, all.covars)
        mD <- get.metadata(mD)[, ..all.terms] %>% 
          build.rec(response = diagnosis.col.name) %>% 
          prep() %>% 
          bake(new_data = get.metadata(mD)[, ..all.terms]) %>% 
          as.data.table() %>% 
          setkeyv(mD@Sample.col) %>% 
          replace.metadata(mD, new.tbl = .)
        # mD %<>% beta.diversity(metrics = betas)
        frm <- paste("~", diagnosis.col.name, "+", paste(all.covars, collapse = " + "))
        base.mods <- microbData::ordinate(
          mD,
          method = "dbRDA", 
          formula = as.formula(frm),
          update.mD = FALSE,
          sqrt.dist = FALSE
        )
        cat("Model Adj. R-squares: ", sep = "\n")
        mod.name <- sapply(base.mods[-2], function(x) RsquareAdj(x)$adj.r.squared) %T>%
          print() %>% 
          which.max() %>% 
          names()
        curr.mod <- base.mods[[mod.name]]
        terms <- attributes(curr.mod$terms)$term.labels %>% 
          str_subset(paste0(diagnosis.col.name, "|Condition"), negate = T)
        for (term in terms) {
          new.frm <- paste0(". ~ . + ", diagnosis.col.name, ":", term)
          new.mod <- update(curr.mod, as.formula(new.frm))
          test.tbl <- anova(curr.mod, new.mod)
          print(test.tbl)
          if (tidy(test.tbl)$p.value[2] < 0.05) { curr.mod <- new.mod }
        }
        list(
          curr.mod,
          anova(curr.mod, by = "terms")
        ) %>% 
          set_names(mod.name, "ANOVA") %>% 
          return()
      })
    })
  })
})

updrs.scores <- str_subset(names(mDs$stool$taxon$Taxon@Metadata), "\\.score$") %>% 
  set_names()
updrs.scores <- updrs.scores[1]
updrs.ords.file.base <- "_updrs_ords_list.rds"
updrs.ords <- lapply(sources, function(source) {
  lapply(focal.sets, function(set) {
    lapply(set, function(ft) {
      set.name <- get.set.name(ft)
      track.values(source, ft)
      
      ords.file <- file.path(
        dirs$saved,
        paste0(source, "_", ft, updrs.ords.file.base)
      )
      redo.if(redo.var, ords.file, max.save.mb = max.save.size, {
        element <- list(source, set.name, ft)
        mD <- pluck(mDs, !!!element)
        all.terms <- c(mD@Sample.col, updrs.scores)
        keep.smpls <- get.metadata(mD)[get(diagnosis.col.name) == "PD case", ..all.terms] %>% 
          .[complete.cases(.)] %>% 
          pluck(mD@Sample.col)
        mD <- keep.samples(mD, keep.smpls)
        frm <- paste("~", paste(updrs.scores, collapse = " + "))
        microbData::ordinate(
          mD,
          method = "dbRDA", 
          formula = as.formula(frm),
          update.mD = FALSE,
          sqrt.dist = FALSE
        ) %>% return()
      }) %>% return()
    })
  })
})

updrs.perms.file <- file.path(
  dirs$saved,
  "rarAbund_updrs_ord_permanovas_dt.rds"
)
updrs.perms <- redo.if(redo.var, updrs.perms.file, {
  lapply(sources, function(source) {
    lapply(focal.sets, function(set) {
      lapply(set, function(ft) {
        set.name <- get.set.name(ft)
        track.values(source, ft)
        element <- list(source, set.name, ft)
        ord.names <- pluck(updrs.ords, !!!element) %>% 
          names() %>% 
          set_names()
        lapply(ord.names, function(ord.name) {
          track.values(ord.name, prepend = "  ")
          pluck(updrs.ords, !!!element, ord.name) %>% 
            anova(by = "term") %>% 
            broom::tidy() %>% 
            set_names(tools::toTitleCase(names(.))) %>% 
            as.data.table() %>% 
            .[
              , `:=`(
                Sig = ifelse(P.value <= 0.05, "*", ""), 
                Source = source, 
                Type = set.name,
                Feature = ft,
                Beta = str_split(ord.name, "_")[[1]][1]
              )
            ] %>% 
            setcolorder(tail(1:ncol(.), 4)) %>%
            return()
        }) %>% rbindlist()
      }) %>% rbindlist()
    }) %>% rbindlist()
  }) %>% rbindlist()
})

cmb.mDs <- readRDS(files$cmb_mDs)

cmb2.mDs.file <- file.path(
  dirs$saved, 
  "combined_mDs_for_betaDispersion_and_ordCondSmplType.rds"
)
cmb2.mDs <- redo.if(redo.var, cmb2.mDs.file, max.save.mb = max.save.size, {
  lapply(focal.sets, function(set) {
    lapply(set, function(ft) {
      set.name <- get.set.name(ft)
      track.values(set.name, ft)
      filter.samples(cmb.mDs[[set.name]][[ft]], !is.na(get(diagnosis.col.name))) %>%
        beta.diversity(metrics = betas) %>%
        beta.dispersion(group = diagnosis.col.name) %>% 
        microbData::ordinate(
          formula = as.formula(
            paste("~", diagnosis.col.name, "+ Condition(sample_type)")
          )
        )
    })
  })
})

cmb2.ord.dts <- lapply(focal.sets, function(set) {
  lapply(set, function(ft) {
    set.name <- get.set.name(ft)
    track.values(set.name, ft)
    ordination.coords(
      mD = cmb2.mDs[[set.name]][[ft]], 
      constraint.coords = TRUE, 
      axis.digits = 4
    ) %>%
      return()
  })
})

cmb2.aov.dt.file <- file.path(
  dirs$saved,
  "combined_mD_beta_anovas_condSmplType.rds"
)
cmb2.aov.dt <- redo.if(redo.var, cmb2.aov.dt.file, {
  lapply(focal.sets, function(set) {
    lapply(set, function(ft) {
      set.name <- get.set.name(ft)
      track.values(set.name, ft)
      ords <- cmb2.mDs[[set.name]][[ft]]@Other.data$Ordinations
      lapply(1:length(ords), function(i) {
        cat(paste(" ", names(ords)[i]), sep = "\n")
        dt <- anova(ords[[i]], by = "terms") %>% tidy() %>% as.data.table()
        dt[, `:=`(Feature = ft, Beta = str_remove(names(ords)[i], "_dbRDA$"))]
        setcolorder(dt, {ncol(dt) - 1}:ncol(dt))
      }) %>% rbindlist()
    }) %>% rbindlist()
  }) %>% rbindlist()
})

plot.exts <- c(PNG = ".png", PDF = ".pdf")
dispersion.tests.dt <- lapply(sources, function(source) {
  lapply(focal.sets, function(ft) {
    track.values(source, ft)
    set.name <- get.set.name(ft)
    element <- list(source, set.name, ft)
    mD <- pluck(mDs, !!!element)
    mD %<>% beta.dispersion(group = diagnosis.col.name)
    dt <- get.metadata(mD)
    disp.cols <- str_subset(names(dt), paste(betas, collapse = "|"))
    fig.files <- file.path(
      dirs$figs, 
      paste("suppFig_betadisp", source, ft, sep = "_")
    ) %>% 
      fig.file.names(plot.exts)
    {
      melt(dt, measure.vars = disp.cols) %>% 
        .[, variable := str_remove(variable, "\\.dispersion")] %>% 
        ggplot(aes(x = !!sym(diagnosis.col.name), y = value)) + 
        geom_violin(aes(fill = !!sym(diagnosis.col.name)), show.legend = F) +
        color.scales$fill$Diagnosis + 
        stat_summary(fun.data = "mean_cl_boot", geom = "errorbar", color = "black") +
        facet_wrap(~ variable, scales = "free_y") +
        labs(x = "PD diagnosis", y = "Beta-dispersion")
      } %>% 
      mggsave(files = fig.files)
    lapply(betas, function(beta) {
      frm <- paste0("`", beta, ".dispersion` ~ ", diagnosis.col.name)
      res <- wilcox.test(as.formula(frm), data = dt) %>% tidy() %>% as.data.table()
      res[
        , `:=`(
          method = NULL, 
          alternative = NULL, 
          sample_type = tools::toTitleCase(source), 
          Feature = ft, 
          Beta = beta
        )
      ]
      setcolorder(res, {ncol(res) - 2}:ncol(res))
    }) %>% rbindlist()
  }) %>% rbindlist()
}) %>% rbindlist()
dispersion.tests.dt[, Adj.pval := p.adjust(p.value, method = "fdr")]
dispersion.tests.dt[, Sig := ifelse(Adj.pval <= 0.05, "*", "")]
# View(dispersion.tests.dt)

## PD meds in PD cases ONLY
pd.meds.ords.file.base <- "_pd.meds_ords_list.rds"
pd.meds.ords <- lapply(sources, function(source) {
  lapply(focal.sets, function(set) {
    lapply(set, function(ft) {
      set.name <- get.set.name(ft)
      track.values(source, ft)
      
      ords.file <- file.path(
        dirs$saved,
        paste0(source, "_", ft, pd.meds.ords.file.base)
      )
      redo.if(redo.var, ords.file, max.save.mb = max.save.size, {
        element <- list(source, set.name, ft)
        mD <- pluck(mDs, !!!element)
        keep.vars <- c(externalID.col.name, med.vars)
        keep.smpls <- get.metadata(mD)[, ..keep.vars] %>% 
          .[complete.cases(.)] %>% 
          pluck(mD@Sample.col)
        mD <- keep.samples(mD, keep.smpls)
        
        meds.dt <- get.metadata(mD)[, ..keep.vars]
        for (med.grp in med.grps) {
          grp.vars <- med.dict[Name %in% med.vars & Group == med.grp]$Name
          meds.dt[[make.names(med.grp)]] <- rowSums(meds.dt[, ..grp.vars])
        }
        med.grps %<>% make.names()
        meds.dt$N.types.meds <- rowSums(meds.dt[, .SD > 0, .SDcols = med.grps])
        setkeyv(meds.dt, mD@Sample.col)
        mD %<>% replace.metadata(meds.dt)
        
        frm1 <- paste("~", paste(med.grps, collapse = " + "))
        ords1 <- microbData::ordinate(
          mD,
          method = "dbRDA", 
          formula = as.formula(frm1),
          update.mD = FALSE,
          sqrt.dist = FALSE
        )
        frm2 <- "~ N.types.meds"
        ords2 <- microbData::ordinate(
          mD,
          method = "dbRDA", 
          formula = as.formula(frm2),
          update.mD = FALSE,
          sqrt.dist = FALSE
        )
        list(Med.groups = ords1, Number.meds = ords2) %>% 
          return()
      }) %>% return()
    })
  })
})

pd.meds.perms <- lapply(sources, function(source) {
  lapply(focal.sets, function(set) {
    lapply(set, function(ft) {
      set.name <- get.set.name(ft)
      lapply(c("Med.groups", "Number.meds") %>% set_names(), function(ord.type) {
        lapply(betas, function(beta) {
          track.values(source, ft, ord.type, beta)
          element <- list(source, set.name, ft, ord.type, paste0(beta, "_dbRDA"))
          pluck(pd.meds.ords, !!!element) %>% 
            anova(by = "terms") %>% 
            broom::tidy() %>% 
            as.data.table() %>% 
            .[
              , `:=`(
                Beta = beta, 
                Mod.type = ord.type, 
                Feature = ft, 
                Source = source
              )
            ] %>% 
            return()
        }) %>% rbindlist()
      }) %>% rbindlist()
    }) %>% rbindlist()
  }) %>% rbindlist()
}) %>% rbindlist()
pd.meds.perms[
  , Adj.pval := p.adjust(p.value, method = "fdr"), 
  by = c("Feature", "Source")
]
setcolorder(pd.meds.perms, c(9, 8, 6, 7))
pd.meds.perms[, Sig := ifelse(Adj.pval < 0.05, "*", "")]
pd.meds.perms[order(Source, Feature, Beta, Mod.type)] %>% 
  kable(digits = c(rep(0, 6), 2, 2, 3, 3, 0)) %>% 
  kable_styling(full_width = F) %>% 
  collapse_rows(columns = 1:4)

## Dopamine agonists in all participants
pd.DAmeds.ords.file.base <- "_pd.DAmeds_ords_list.rds"
pd.DAmeds.ords <- lapply(sources, function(source) {
  lapply(focal.sets, function(set) {
    lapply(set, function(ft) {
      set.name <- get.set.name(ft)
      track.values(source, ft)
      
      ords.file <- file.path(
        dirs$saved,
        paste0(source, "_", ft, pd.DAmeds.ords.file.base)
      )
      redo.if(redo.var, ords.file, max.save.mb = max.save.size, {
        element <- list(source, set.name, ft)
        mD <- pluck(mDs, !!!element)
        metadata <- get.metadata(mD)
        for (med.var in med.vars) {
          med.name <- str_split(med.var, "\\.")[[1]] %>% tail(1)
          metadata[[med.var]] <- ifelse(
            is.na(metadata[[med.var]]), 
            paste("Not", med.name), 
            metadata[[med.var]]
          )
        }
        keep.vars <- c(externalID.col.name, diagnosis.col.name, med.vars)
        keep.smpls <- metadata[, ..keep.vars] %>% 
          .[complete.cases(.)] %>% 
          pluck(mD@Sample.col)
        mD <- keep.samples(mD, keep.smpls)
        med.grps <- med.dict[Name %in% med.vars]$Group %>% unique() %>% sort()
        
        meds.dt <- metadata[, ..keep.vars]
        for (med.grp in med.grps) {
          grp.vars <- med.dict[Name %in% med.vars & Group == med.grp]$Name
          meds.dt[[make.names(med.grp)]] <- rowSums(meds.dt[, ..grp.vars])
        }
        med.grps %<>% make.names()
        meds.dt$N.types.meds <- rowSums(meds.dt[, .SD > 0, .SDcols = med.grps])
        setkeyv(meds.dt, mD@Sample.col)
        mD %<>% replace.metadata(meds.dt)
        
        frm <- paste("~", diagnosis.col.name, "* dopamine.agonist")
        ords <- microbData::ordinate(
          mD,
          method = "dbRDA", 
          formula = as.formula(frm),
          update.mD = FALSE,
          sqrt.dist = FALSE
        )
        adonises <- lapply(betas, function(beta) {
          paste("mD@Distance.matrices[[beta]] ~", diagnosis.col.name, "* dopamine.agonist") %>% 
            as.formula() %>% 
            adonis2(data = mD@Metadata)
        })
        list(Ords = ords, Adonises = adonises, Metadatas = meds.dt) %>% 
          return()
      }) %>% return()
    })
  })
})

pd.DAmeds.adns <- lapply(sources, function(source) {
  lapply(focal.sets, function(set) {
    lapply(set, function(ft) {
      set.name <- get.set.name(ft)
      lapply(betas, function(beta) {
        track.values(source, ft, beta)
        element <- list(source, set.name, ft, "Adonises", beta)
        pluck(pd.DAmeds.ords, !!!element) %>% 
          broom::tidy() %>% 
          as.data.table() %>% 
          .[, `:=`(Beta = beta, Feature = ft, Source = source)] %>% 
          return()
      }) %>% rbindlist()
    }) %>% rbindlist()
  }) %>% rbindlist()
}) %>% rbindlist()
setcolorder(pd.DAmeds.adns, 9:7)
pd.DAmeds.adns[, Adj.pval := p.adjust(p.value, method = "fdr"), by = c("Source", "Feature")]
pd.DAmeds.adns[, Sig := ifelse(Adj.pval < 0.05, "*", "")]
kable(pd.DAmeds.adns, digits = c(rep(0, 5), 2, 3, 2, 3, 3, 0)) %>% 
  kable_styling(full_width = F) %>% 
  collapse_rows(columns = 1:3)

pd.DAmeds.aovs <- lapply(sources, function(source) {
  lapply(focal.sets, function(set) {
    lapply(set, function(ft) {
      set.name <- get.set.name(ft)
      lapply(betas, function(beta) {
        track.values(source, ft, beta)
        element <- list(source, set.name, ft, "Ords", paste0(beta, "_dbRDA"))
        ord <- pluck(pd.DAmeds.ords, !!!element)
        attributes(ord)
        pluck(pd.DAmeds.ords, !!!element) %>% 
          anova(by = "terms") %>% 
          broom::tidy() %>% 
          as.data.table() %>% 
          .[, `:=`(Beta = beta, Feature = ft, Source = source)] %>% 
          return()
      }) %>% rbindlist()
    }) %>% rbindlist()
  }) %>% rbindlist()
}) %>% rbindlist()
setcolorder(pd.DAmeds.aovs, 8:6)
pd.DAmeds.aovs[, Adj.pval := p.adjust(p.value, method = "fdr"), by = c("Source", "Feature")]
pd.DAmeds.aovs[, Sig := ifelse(Adj.pval < 0.05, "*", "")]
kable(pd.DAmeds.aovs, digits = c(rep(0, 5), 2, 2, 3, 3, 0)) %>% 
  kable_styling(full_width = F) %>% 
  collapse_rows(columns = 1:3)

pd.DAmeds.coords <- lapply(sources, function(source) {
  lapply(focal.sets, function(set) {
    lapply(set, function(ft) {
      set.name <- get.set.name(ft)
      track.values(source, ft)
      element <- list(source, set.name, ft, "Ords")
      ord.list <- pluck(pd.DAmeds.ords, !!!element)
      element[[4]] <- "Metadatas"
      metadata <- pluck(pd.DAmeds.ords, !!!element)
      ordination.coords(
        ord.list = ord.list,
        metadata = metadata,
        axis.digits = 4,
        constraint.coords = TRUE
      ) %>% return()
    })
  })
})
pd.DAmeds.smpls <- lapply(sources, function(source) {
  lapply(focal.sets, function(set) {
    lapply(set, function(ft) {
      set.name <- get.set.name(ft)
      track.values(source, ft)
      element <- list(source, set.name, ft, "Samples")
      pluck(pd.DAmeds.coords, !!!element) %>% 
        .[, `:=`(Source = source, Feature = ft)]
    }) %>% rbindlist()
  }) %>% rbindlist()
}) %>% rbindlist()

list(
  {
    ggplot(pd.DAmeds.smpls, aes(x = Axis1, y = Axis2, color = !!sym(diagnosis.col.name))) + 
      geom_point() +
      stat_ellipse() +
      color.scales$color$Diagnosis +
      ggh4x::facet_nested(Beta.metric ~ Source + Feature, scales = "free")
  },
  {
    ggplot(pd.DAmeds.smpls, aes(x = Axis1, y = Axis2, color = factor(dopamine.agonist))) + 
      geom_point() +
      stat_ellipse() +
      scale_color_manual(
        name = "Num. dopamine agonists", 
        values = c("purple", "green", "orange")
      ) +
      ggh4x::facet_nested(Beta.metric ~ Source + Feature, scales = "free")
  }
) %>% 
  cowplot::plot_grid(plotlist = ., nrow = 1)

### Batch effects
batches.dt <- read.table(
  "", # batches file path
  sep = "\t", 
  header = TRUE
) %>% as.data.table()

batch.ords.file.base <- "_batch_ords_list.rds"
batch.ords <- lapply(sources, function(source) {
  lapply(focal.sets, function(set) {
    lapply(set, function(ft) {
      set.name <- get.set.name(ft)
      track.values(source, ft)
      
      ords.file <- file.path(
        dirs$saved,
        paste0(source, "_", ft, batch.ords.file.base)
      )
      redo.if(redo.var, ords.file, max.save.mb = max.save.size, {
        element <- list(source, set.name, ft)
        mD <- pluck(mDs, !!!element)
        metadata <- get.metadata(mD) %>% 
          merge(batches.dt, by = externalID.col.name, all.x = TRUE)
        
        keep.vars <- c(externalID.col.name, diagnosis.col.name, batchID.col.name)
        keep.smpls <- metadata[, ..keep.vars] %>% 
          .[complete.cases(.)] %>% 
          pluck(mD@Sample.col)
        mD <- keep.samples(mD, keep.smpls)
        metadata <- metadata[keep.smpls]
        mD %<>% replace.metadata(metadata)
        
        frm <- paste("~", diagnosis.col.name,  "*", batchID.col.name)
        ords <- microbData::ordinate(
          mD,
          method = "dbRDA", 
          formula = as.formula(frm),
          update.mD = FALSE,
          sqrt.dist = FALSE
        )
        adonises <- lapply(betas, function(beta) {
          paste(
            "mD@Distance.matrices[[beta]] ~", 
            diagnosis.col.name, "*", batchID.col.name
          ) %>% 
            as.formula() %>% 
            adonis2(data = mD@Metadata)
        })
        list(Ords = ords, Adonises = adonises, Metadatas = metadata) %>% 
          return()
      }) %>% return()
    })
  })
})

batch.adns <- lapply(sources, function(source) {
  lapply(focal.sets, function(set) {
    lapply(set, function(ft) {
      set.name <- get.set.name(ft)
      lapply(betas, function(beta) {
        track.values(source, ft, beta)
        element <- list(source, set.name, ft, "Adonises", beta)
        pluck(batch.ords, !!!element) %>% 
          broom::tidy() %>% 
          as.data.table() %>% 
          .[, `:=`(Beta = beta, Feature = ft, Source = source)] %>% 
          return()
      }) %>% rbindlist()
    }) %>% rbindlist()
  }) %>% rbindlist()
}) %>% rbindlist()
setcolorder(batch.adns, 9:7)
batch.adns[, Adj.pval := p.adjust(p.value, method = "fdr"), by = c("Source", "Feature")]
batch.adns[, Sig := ifelse(Adj.pval < 0.05, "*", "")]
kable(batch.adns, digits = c(rep(0, 5), 2, 3, 2, 3, 3, 0)) %>% 
  kable_styling(full_width = F) %>% 
  collapse_rows(columns = 1:3)

batch.coords <- lapply(sources, function(source) {
  lapply(focal.sets, function(set) {
    lapply(set, function(ft) {
      set.name <- get.set.name(ft)
      track.values(source, ft)
      element <- list(source, set.name, ft, "Ords")
      ord.list <- pluck(batch.ords, !!!element)
      element[[4]] <- "Metadatas"
      metadata <- pluck(batch.ords, !!!element)
      ordination.coords(
        ord.list = ord.list,
        metadata = metadata,
        axis.digits = 4,
        constraint.coords = TRUE
      ) %>% return()
    })
  })
})

batch.smpls <- lapply(sources, function(source) {
  lapply(focal.sets, function(set) {
    lapply(set, function(ft) {
      set.name <- get.set.name(ft)
      track.values(source, ft)
      element <- list(source, set.name, ft, "Samples")
      pluck(batch.coords, !!!element) %>% 
        .[, `:=`(Source = source, Feature = ft)]
    }) %>% rbindlist()
  }) %>% rbindlist()
}) %>% rbindlist()

list(
  {
    ggplot(batch.smpls, aes(x = Axis1, y = Axis2, color = !!sym(diagnosis.col.name))) + 
      geom_point() +
      stat_ellipse() +
      color.scales$color$Diagnosis +
      ggh4x::facet_nested(Beta.metric ~ Source + Feature, scales = "free")
  },
  {
    ggplot(batch.smpls, aes(x = Axis1, y = Axis2, color = !!sym(batchID.col.name))) + 
      geom_point() +
      stat_ellipse() +
      ggh4x::facet_nested(Beta.metric ~ Source + Feature, scales = "free")
  }
) %>% 
  cowplot::plot_grid(plotlist = ., nrow = 1)

###

regions.dt <- readRDS("") # regions dt file path
region.ords.file.base <- "_region_ords_list.rds"
region.ords <- lapply(sources, function(source) {
  lapply(focal.sets, function(set) {
    lapply(set, function(ft) {
      set.name <- get.set.name(ft)
      track.values(source, ft)
      
      ords.file <- file.path(
        dirs$saved,
        paste0(source, "_", ft, region.ords.file.base)
      )
      redo.if(redo.var, ords.file, max.save.mb = max.save.size, {
        element <- list(source, set.name, ft)
        mD <- pluck(mDs, !!!element)
        metadata <- get.metadata(mD) %>% 
          merge(regions.dt, by = c(smplID.col.name, diagnosis.col.name), all.x = TRUE)
        setkeyv(metadata, mD@Sample.col)
        keep.vars <- c(externalID.col.name, diagnosis.col.name, "Region")
        metadata <- metadata[, ..keep.vars] %>% 
          .[complete.cases(.)] %>% 
          unique()
        setkeyv(metadata, mD@Sample.col)
        keep.smpls <- metadata[[mD@Sample.col]]
        mD <- keep.samples(mD, keep.smpls)
        
        mD %<>% replace.metadata(metadata)
        
        frm <- paste("~", diagnosis.col.name, "* Region")
        ords <- microbData::ordinate(
          mD,
          method = "dbRDA", 
          formula = as.formula(frm),
          update.mD = FALSE,
          sqrt.dist = FALSE
        )
        adonises <- lapply(betas, function(beta) {
          paste("mD@Distance.matrices[[beta]] ~", diagnosis.col.name, "* Region") %>%
            as.formula() %>% 
            adonis2(data = mD@Metadata)
        })
        list(Ords = ords, Adonises = adonises, Metadatas = metadata) %>% 
          return()
      }) %>% return()
    })
  })
})

region.adns <- lapply(sources, function(source) {
  lapply(focal.sets, function(set) {
    lapply(set, function(ft) {
      set.name <- get.set.name(ft)
      lapply(betas, function(beta) {
        track.values(source, ft, beta)
        element <- list(source, set.name, ft, "Adonises", beta)
        pluck(region.ords, !!!element) %>% 
          broom::tidy() %>% 
          as.data.table() %>% 
          .[, `:=`(Beta = beta, Feature = ft, Source = source)] %>% 
          return()
      }) %>% rbindlist()
    }) %>% rbindlist()
  }) %>% rbindlist()
}) %>% rbindlist()
setcolorder(region.adns, 9:7)
region.adns[, Adj.pval := p.adjust(p.value, method = "fdr"), by = c("Source", "Feature")]
region.adns[, Sig := ifelse(Adj.pval < 0.05, "*", "")]
kable(region.adns, digits = c(rep(0, 5), 2, 3, 2, 3, 3, 0)) %>% 
  kable_styling(full_width = F) %>% 
  collapse_rows(columns = 1:3)
