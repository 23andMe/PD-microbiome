# model_specs.R

require(tidymodels)

rf.spec <- rand_forest(
  mode = "classification",
  mtry = tune(), 
  trees = tune(),
  min_n = 1
) %>% 
  set_engine(
    "ranger", 
    importance = "permutation", 
    splitrule = tune()
  ) %>% 
  translate()
