# setup.R

base.dir <- "" # SET WORKING DIRECTORY
if (getwd() != base.dir) { setwd(base.dir) }
source(".Rprofile")

knitr::opts_chunk$set(
  echo = FALSE, 
  warning = FALSE, 
  message = FALSE, 
  dpi = 300
)
source("_packages_and_sources.R")
theme_set(theme_cowplot())
my.theme <- theme_update(
  legend.position = "top",
  legend.box.just = "left",
  legend.text = element_text(size = 8),
  legend.title = element_text(size = 10),
  legend.justification = "left",
  legend.key.size = unit(0.7, "line"),
  plot.caption = element_text(hjust = 0, size = 8)
)
options(knitr.kable.NA = "")
options(tidymodels.dark = TRUE)
max.save.size <- 10000
to.run <- c("req.first", "tax.alpha", "beta", "lda", "network", "rand.forest")
to.skip <- c() 
assign.redo(to.run, state = T)
if (length(to.skip) > 0) { set.redo.false(to.skip) }

taxa.otu.combined.mD.file <- file.path(dirs$mds, "") ### SET APPROPRIATE FILE NAME
func.otu.combined.mD.file <- file.path(dirs$mds, "") ### SET APPROPRIATE FILE NAME

internalID.col.name <- "" # name of internal sample IDs column
externalID.col.name <- "" # name of internal sample IDs column
diagnosis.col.name <- "" # name of PD diagnosis column
wkly.BMs.col.name <- "" # name of number weekly BMs column
flossing.col.name <- "" # name of flossing frequency column
age.col.name <- "" # name of current age column
batchID.col.name <- "" # name of batch ID column

prev.cut <- 0.1
abund.cut <- 0.01 / 100

seed <- 42
sources <- c("saliva", "stool") %>% set_names()
taxon.lvls <- c("Taxon", "Species", "Genus") %>% set_names()
func.lvls <- c("KO", "Module", "Pathway") %>% set_names()
feature.sets <- list(taxon = taxon.lvls, func = func.lvls)
focal.sets <- lapply(feature.sets, `[`, 1)

alphas <- c("Chao1", "InvSimpson") %>% set_names()
betas <- c("Sorensen", "Bray-Curtis") %>% set_names()
pd.statuses <- c("Control", "PD case") %>% set_names()
cmb.methods <- c("split", "lumped") %>% set_names()
agg.methods <- list(Taxon.Species = taxon.lvls[1:2], Taxon.Species.Genus = taxon.lvls)
if (file.exists(files$ctrl_vars)) {
  ctrl.vars <- readRDS(files$ctrl_vars)
} else {
  rlang::inform(
    paste(
      "File", files$ctrl_vars, "doesn't exist; make sure to run analysis file 01_required_first.R before proceeding with other analyses"
    )
  )
}

para.cores <- 50

na.color <- "grey50"
age.binwd <- 5
create.scale <- function(
  type = c("color", "fill"), 
  palette.source = c("brewer", "sjplot", "manual"), 
  ...
) {
  switch(
    type,
    color = switch(
      palette.source,
      brewer = scale_color_brewer(...),
      sjplot = scale_color_sjplot(...),
      manual = scale_color_manual(...)
    ),
    fill = switch(
      palette.source,
      brewer = scale_fill_brewer(...),
      sjplot = scale_fill_sjplot(...),
      manual = scale_fill_manual(...)
    )
  )
}

color.scales <- lapply(
  setNames(c("color", "fill"), c("color", "fill")), 
  function(scale.type) {
    list(
      Diagnosis = create.scale(
        type = scale.type, 
        palette.source = "manual", 
        name = "PD diagnosis", 
        values = c("#0072B2", "#CC79A7"), 
        na.value = na.color
      ),
      Floss_freq = create.scale(
        type = scale.type,
        palette.source = "sjplot",
        name = "Floss teeth freq.",
        palette = "ipsum",
        na.value = na.color
      ),
      Bristol_stool_type_4 = create.scale(
        type = scale.type,
        palette.source = "brewer",
        name = "Bristol score 4 freq.",
        palette = "Set2",
        na.value = na.color
      ),
      Bristol_stool_type_6 = create.scale(
        type = scale.type,
        palette.source = "brewer",
        name = "Bristol score 6 freq.",
        palette = "Dark2",
        na.value = na.color
      ),
      Number_weekly_BMs = create.scale(
        type = scale.type,
        palette.source = "sjplot",
        name = "Num. BMs weekly",
        palette = "system",
        na.value = na.color
      ),
      Laxative_past_month = create.scale(
        type = scale.type,
        palette.source = "manual",
        name = "Used laxative past month?",
        values = brewer.pal(4, name = "Set2")[3:4],
        na.value = na.color
      ),
      Age = create.scale(
        type = scale.type,
        palette.source = "brewer",
        name = "Age",
        palette = "Set1",
        na.value = na.color
      )
    ) %>% return()
  })

centr.shapes <- c(3, 4, 8, 2, 6, 0, 5, 1)

plurals <- data.table(
  Singular = c(
    "Taxon",
    "Species",
    "Genus",
    "KO",
    "Module",
    "Pathway"
  ),
  Plural = c(
    "OTUs", 
    "Species", 
    "Genera", 
    "KOs", 
    "Modules", 
    "Pathways"
    )
)
setkey(plurals, Singular)

assignments <- readRDS(taxa.otu.combined.mD.file) %>% 
  get.assignments()
names(assignments) %<>% str_replace("OTU", "Taxon")
families <- sort(unique(assignments$Family))
exclude <- c(
  "gray",
  "grey",
  "white",
  "light",
  "pale",
  "ivory",
  "snow",
  "azure",
  "honeydew",
  "cornsilk",
  "beige",
  "seagreen",
  "lemon",
  "linen",
  "iceblue",
  "seashell",
  "oldlace"
) %>% paste(collapse = "|")
rand.colors <- grDevices::colors()[
  grep(exclude, grDevices::colors(), invert = T)
]
set.seed(seed * 4)
fam.colors <- data.table(
  Family = families,
  Color = sample(rand.colors, size = length(families), replace = F)
)
setkey(fam.colors, Family)

assignments <- readRDS(func.otu.combined.mD.file) %>% 
  get.assignments()
pathways <- kegg.dicts$Pathway[sort(unique(assignments$Pathway))]$Path.name
path.colors <- data.table(
  Pathway = pathways,
  Color = sample(rand.colors, size = length(pathways), replace = F)
)
setkey(path.colors, Pathway)

brite.classes <- sort(unique(kegg.links$ko.class$Brite.class))
set.seed(seed * 6)
brite.class.colors <- data.table(
  Class = brite.classes,
  Color = sample(rand.colors, size = length(brite.classes), replace = F)
)
setkey(brite.class.colors, Class)

mod.path.classes <- c(kegg.links$mod.class$Class, kegg.links$path.class$Class) %>% 
  unique() %>% 
  sort() 
set.seed(seed * 7)
class.colors <- data.table(
  Class = mod.path.classes,
  Color = sample(rand.colors, size = length(mod.path.classes), replace = F)
)
setkey(class.colors, Class)
