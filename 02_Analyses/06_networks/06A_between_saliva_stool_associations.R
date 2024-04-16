# SpiecEasi between saliva stool sample associations

base.dir <- ""
if (getwd() != base.dir) { setwd(base.dir) }
if (!exists("color.scales")) { source("_setup.R") }
redo.var <- "network"
if (is.null(get.redo.state(redo.var))) {
  assign.redo(redo.var, state = T)
}

cmb.mDs.file <- file.path(dirs$saved, "combined_mDs_glommed.rds")
cmb.mDs <- redo.if(redo.var, cmb.mDs.file, max.save.mb = max.save.size, {
  lapply(focal.sets, function(set) {
    lapply(set, function(ft) {
      set.name <- get.set.name(ft)
      track.values(set.name, ft)
      mD <- list.files(
        path = dirs$mds, 
        pattern = paste0(set.name, ".+combined"), 
        full.names = TRUE
      ) %>% readRDS()
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
      )
      mD <- filter.samples(mD, sample.sums(mD) > 0)
      if (is.null(attributes(mD@Metadata)$sorted)) {
        setkeyv(mD@Metadata, mD@Sample.col)
      }
      if (is.null(attributes(mD@Assignments)$sorted)) {
        setkeyv(mD@Assignments, mD@Feature.col)
      }
      return(mD)
    })
  })
})

se.method <- "mb"
sel.crit <- "stars"
reps <- 20
n.lambda <- 20

options(mc.cores = reps)

se.btw.src <- lapply(focal.sets, function(set) {
  lapply(set, function(ft) {
    set.name <- get.set.name(ft)
    mD <- cmb.mDs %>% pluck(set.name, ft)
    lapply(pd.statuses, function(status) {
      track.values(ft, status)
      meta.dt <- get.metadata(mD)
      saliva.dt <- keep.samples(
        mD,
        samples = mD@Metadata[
          get(diagnosis.col.name) == status & sample_type == "Saliva"
        ][[externalID.col.name]]
      ) %>%
        get.abundances(as.DT = T)
      saliva.dt[
        , (internalID.col.name) := meta.dt[get(externalID.col.name)][[internalID.col.name]]
      ]
      saliva.dt[, (externalID.col.name) := NULL]
      setcolorder(saliva.dt, ncol(saliva.dt))
      names(saliva.dt)[-1] <- paste0("saliva__", names(saliva.dt)[-1])
      
      stool.dt <- keep.samples(
        mD,
        samples = mD@Metadata[
          get(diagnosis.col.name) == status & sample_type == "Stool"
        ][[externalID.col.name]]
      ) %>%
        get.abundances(as.DT = T)
      stool.dt[
        , (internalID.col.name) := meta.dt[get(externalID.col.name)][[internalID.col.name]]
      ]
      stool.dt[, (externalID.col.name) := NULL]
      setcolorder(stool.dt, ncol(stool.dt))
      names(stool.dt)[-1] <- paste0("stool__", names(stool.dt)[-1])
      abund.dt <- merge(saliva.dt, stool.dt)
      abund.mat <- as.matrix(abund.dt[, -1])
      if (ft == "Taxon") {
        assign.dt <- get.assignments(mD)
        bad.otus <- assign.dt[Phylum == "Bacteria_Phylum"]$OTU %>% 
          lapply(c("saliva__", "stool__"), paste0, .) %>% 
          unlist()
        abund.mat <- abund.mat[, !{colnames(abund.mat) %in% bad.otus}]
      }
      rownames(abund.mat) <- abund.dt[[internalID.col.name]]
      
      save.file <- file.path(
        dirs$saved, 
        paste(
          "spieceasi_btwnSources", 
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

