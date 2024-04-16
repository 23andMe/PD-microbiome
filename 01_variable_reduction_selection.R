# variable_reduction_selection.R

source("_setup.R")

files$base_metadata <- "" ### SET NAME OF INPUT METADATA FILE
update.files(
  "working_metadata",
  file.path(dirs$input, "working_metadata.rds")
)
metadata <- readRDS(files$base_metadata)

metadata[, BMI := {weight / 1000} / {{height / 1000}^2}]
metadata[BMI > 100]$BMI <- median(metadata[BMI < 100]$BMI)
metadata <- metadata[
  { between(current_age, 50, 69) | is.na(current_age) } &
    !is.na(get(diagnosis.col.name)) &
    gastroenteritis != "yes"
]

# all of the patterns below match variable names in the MJFF Fox DEN database,
# not 23andMe's
for (col in str_subset(names(metadata), "PDMBStl")) { 
  new.col <- paste0(col, "_reduced")
  lvls <- c("never", "sometimes", "frequently")
  if (str_detect(col, "Snake|Blob")) { lvls %<>% rev() }
  metadata[
    , (new.col) := ifelse(
      get(col) %in% c(
        "3-6 times per week", 
        "Daily", 
        "2-3 times per day", 
        "More than 3 times per day"
      ),
      "frequently",
      ifelse(get(col) == "Less than 3 times per week", "sometimes", "never")
    ) %>% 
      factor(levels = lvls)
  ]
  cat(new.col, sep = "\n")
  print(table(metadata[[new.col]], useNA = "always"))
}
metadata[
  ,`:=` (
    PDMBFloFreq = ifelse(
      as.character(PDMBFloFreq) %in% c("Once a day", "Twice a day or more"), 
      "Regularly", 
      as.character(PDMBFloFreq)
    ) %>% 
      factor(levels = c("Regularly", "Sometimes", "Rarely", "Never")),
    PDMBBowlMovWk = factor(
      PDMBBowlMovWk, 
      levels = c("Six or more", "Between 3 and 5", "Less than 3")
    ),
    PDMBLaxMonth = factor(PDMBLaxMonth, levels = c("No", "Yes"))
  )
]

saveRDS(metadata, files$working_metadata)

covars <- names(metadata)[
  -str_which(
    names(metadata), 
    "" # remove factors, like PD diagnosis, you don't want to consider
  )
]
covar.test.dt <- lapply(covars, function(covar) {
  res.dt <- data.table(
    Covar = covar, 
    Responses = capture.output(str(metadata[[covar]]))
  ) %>% unique()
  cols <- c(diagnosis.col.name, covar)
  
  if (class(metadata[[covar]]) %in% c("numeric", "integer")) {
    case.types <- na.omit(metadata[, ..cols])[[diagnosis.col.name]] %>% 
      unique()
    
    if (length(case.types) < 2) {
      res.dt[
        , `:=`(
          Test.type = NA, 
          P.val = NA, PD.vs.Ctrl = paste(case.types, "only")
        )
      ]
    } else {
      aov.dt <- wilcox.test(
        as.formula(paste(covar, "~", diagnosis.col.name)), 
        data = metadata, 
        na.action = na.omit
      ) %>% 
        tidy() %>% 
        as.data.table()
      res.dt[
        , `:=`(
          Test.type = "Wilcoxon", 
          P.val = aov.dt$p.value[1],
          PD.vs.Ctrl = ifelse(
            aov.dt$p.value[1] <= 0.05, 
            "significant", 
            "non-significant"
          )
        )
      ]
    }
  } else {
    # chi-sq
    covar.mat <- na.omit(metadata[, ..cols])[, .(N = .N), by = cols] %>% 
      dcast(as.formula(paste(covar, "~", diagnosis.col.name)), value.var = "N") %>% 
      as.matrix(rownames = covar)
    covar.mat[is.na(covar.mat)] <- 0
    
    if (ncol(covar.mat) < 2) {
      res.dt[
        , `:=`(
          Test.type = NA, 
          P.val = NA, 
          PD.vs.Ctrl = paste(colnames(covar.mat), "only")
        )
      ]
    } else {
      chisq.dt <- chisq.test(covar.mat) %>% 
        tidy() %>% 
        as.data.table()
      res.dt[
        , `:=`(
          Test.type = "Chi-squared", 
          P.val = chisq.dt$p.value,
          PD.vs.Ctrl = ifelse(
            chisq.dt$p.value <= 0.05, 
            "significant", 
            "non-significant"
          )
        )
      ]
    }
  }
  return(res.dt)
}) %>% rbindlist()
write.table(
  covar.test.dt[, -3],
  file = file.path(dirs$plots, "covariate_test_results.tsv"), 
  row.names = FALSE,
  sep = "\t",
  quote = TRUE
)

require(ggplot2)
require(ggbeeswarm)
require(cowplot)
theme_set(theme_cowplot())

covar.test.sig <- covar.test.dt[PD.vs.Ctrl == "significant"]
covar.test.sig.file <- file.path(
  dirs$saved, 
  "PDstatus_significant_variables.rds"
)
update.files("status_sig_vars", covar.test.sig.file)
saveRDS(covar.test.sig, covar.test.sig.file)

fig.file <- file.path(dirs$plots, "sig_diff_covars_by_PDstatus.pdf")
if (file.exists(fig.file)) { file.remove(fig.file) }
pdf(
  file = fig.file,
  width = 10,
  height = 10
)
lapply(1:nrow(covar.test.sig), function(r) {
  row <- covar.test.sig[order(P.val)][r]
  cat(r, sep = "\n")
  if (row$Test.type == "Chi-squared") {
    p <- lapply(levels(metadata[[diagnosis.col.name]]), function(status) {
      metadata[
        get(diagnosis.col.name) == status, 
        .(
          Prop = .N / sum(metadata[[diagnosis.col.name]] == status, na.rm = T),
          .N
        ), 
        by = c(diagnosis.col.name, row$Covar)
      ]
    }) %>% 
      rbindlist() %>% 
      ggplot(
        aes(
          x = reorder(!!sym(row$Covar), -N),
          y = Prop
        )
      ) +
      geom_text(
        aes(label = N, group = !!sym(diagnosis.col.name)),
        position = position_dodge(width = 1),
        vjust = -1
      ) +
      geom_col(aes(fill = !!sym(diagnosis.col.name)), position = "dodge") +
      scale_y_continuous(expand = expansion(mult = c(0, 0.1))) + 
      color.scales$fill$Diagnosis +
      labs(x = row$Covar, y = "Proportion within PD status\n(Control/PD case)")
    if (row$Covar != covar.test.sig[order(P.val)][1]$Covar) {
      p <- p + theme(legend.position = "none")
    }
  } else {
    p <- ggplot(metadata, aes(x = !!sym(diagnosis.col.name), y = !!sym(row$Covar))) +
      geom_quasirandom(aes(color = !!sym(diagnosis.col.name))) +
      stat_summary(
        fun.data = "mean_cl_boot",
        color = "black",
        geom = "errorbar"
      ) + 
      color.scales$color$Diagnosis
  }
  print(p + theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1)))
})
dev.off()
