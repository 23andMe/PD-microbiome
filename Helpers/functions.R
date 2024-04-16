# functions.R

require(data.table)
require(rlang)
require(magrittr)
require(purrr)
require(stringr)

table <- function(..., useNA = "always") { base::table(..., useNA = useNA) }

fig.file.names <- function(base.name, plot.exts) {
  lapply(plot.exts, function(x) paste0(base.name, x)) %>% 
    return()
}
mggsave <- function(plot, files, ...) {
  for (file in files) {
    ggsave(plot = plot, filename = file, ...)
  }
}

update.files <- function(var, name) {
  save.path <- "_files_env.rds"
  assign(var, name, envir = files)
  saveRDS(object = files, file = save.path)
}

fix.column.class <- function(col, coerce.to = "character") {
  coerce.to <- rlang::arg_match(
    coerce.to, 
    values = c("character", "numeric", "integer", "factor")
  )
  new.col <- sapply(col, function(x) {
    if (is.null(x)) {
      return(NA)
    } else if (is.na(x) | x == "NA") {
      return(NA)
    } else {
      new.x <- as.character(x)
      if (coerce.to != "character") {
        new.x <- as(new.x, coerce.to)
      }
      return(new.x)
    }
  }) %>% unlist()
  return(new.col)
}

build.rec <- function(x, response) {
  paste(response, "~ .") %>% 
    as.formula() %>% 
    recipe(x, formula = .) %>%
    step_string2factor(all_string_predictors()) %>%
    step_impute_knn(all_predictors(), neighbors = 10) %>%
    step_nzv(all_predictors()) %>% 
    step_scale(all_numeric_predictors()) %>% 
    return()
}

build.rec.rf <- function(x, response) {
  paste(response, "~ .") %>% 
    as.formula() %>% 
    recipe(x, formula = .) %>%
    step_string2factor(all_string_predictors()) %>%
    step_impute_knn(all_predictors(), neighbors = 10) %>%
    step_nzv(all_predictors()) %>% 
    return()
}

track.values <- function(..., prepend = NULL, return.val = FALSE) {
  args <- rlang::dots_list(..., .named = TRUE)
  text <- sapply(names(args), function(name) {
    ifelse(
      !is.integer(args[[name]]) & !is.numeric(args[[name]]),
      paste0(name, " <- '", args[[name]], "'"),
      paste(name, "<-", args[[name]])
    ) %>% return()
  }) %>% 
    paste(collapse = "; ") %>% 
    ifelse(is.null(prepend), ., paste(prepend, .)) 
  cat(text, sep = "\n")
  if (return.val) {
    return(text)
  }
} 

feature.relAbunds <- function(mD) {
  feature.sums <- get.abundances(mD) %>% 
    colSums()
  total.sum <- get.abundances(mD) %>% 
    colSums() %>% 
    sum()
  {feature.sums / total.sum} %>% 
    return()
}

feature.prevalences <- function(mD) {
  get.abundances(mD) %>% 
    apply(MARGIN = 1, FUN = function(x) ifelse(x > 0, 1, 0)) %>% 
    rowSums() %>% 
    divide_by(microbData::nsamples(mD)) %>% 
    return()
}

get.set.name <- function(ft) {
  set.name <- sapply(feature.sets, function(x) str_which(x, ft)) %>%
    unlist() %>% 
    names()
}

key.mD.tbls <- function(mD) {
  if (is.null(attributes(mD@Metadata)$sorted)) {
    setkeyv(mD@Metadata, mD@Sample.col)
  }
  if (is.null(attributes(mD@Assignments)$sorted)) {
    setkeyv(mD@Assignments, mD@Feature.col)
  }
  return(mD)
}


