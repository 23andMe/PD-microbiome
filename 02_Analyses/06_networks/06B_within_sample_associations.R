# SpiecEasi within sample associations

base.dir <- ""
if (getwd() != base.dir) { setwd(base.dir) }
if (!exists("color.scales")) { source("_setup.R") }
redo.var <- "network"
if (is.null(get.redo.state(redo.var))) {
  assign.redo(redo.var, state = T)
}

mDs <- readRDS(files$pa_filt_mDs)

se.method <- "mb"
sel.crit <- "stars"
reps <- 50
n.lambda <- 20

options(mc.cores = reps)

se.win.src <- lapply(sources, function(source) {
  lapply(focal.sets, function(set) {
    lapply(set, function(ft) {
      set.name <- get.set.name(ft)
      lapply(pd.statuses, function(status) {
        track.values(source, ft, status)
        mD <- mDs %>% pluck(source, set.name, ft)
        meta.dt <- get.metadata(mD)
        abund.mat <- keep.samples(
          mD, 
          samples = mD@Metadata[get(diagnosis.col.name) == status][[externalID.col.name]]
        ) %>% 
          get.abundances()
        if (ft == "Taxon") {
          assign.dt <- get.assignments(mD)
          bad.otus <- assign.dt[Phylum == "Bacteria_Phylum"]$OTU
          abund.mat <- abund.mat[, !{colnames(abund.mat) %in% bad.otus}]
        }
        rownames(abund.mat) <- meta.dt[rownames(abund.mat)][[internalID.col.name]]
        save.file <- file.path(
          dirs$saved, 
          paste(
            "spieceasi_withinSource", 
            source,
            ft,
            str_remove(status, " "),
            paste0("method", toupper(se.method)), 
            paste0("selCrit", toupper(sel.crit)),
            paste0("reps", reps),
            paste0("nlambda", n.lambda),
            "results.rds",
            sep = "_"
          )
        )
        redo.if(redo.var, save.file, max.save.mb = max.save.size, {
          spiec.easi(
            abund.mat, 
            method = se.method,
            sel.criterion = sel.crit,
            pulsar.select = TRUE,
            nlambda = n.lambda,
            pulsar.params = list(rep.num = reps, ncores = reps)
          )
        }) %>%
          return()
      })
    })
  })
})
