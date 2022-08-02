## orderly::orderly_develop_start(use_draft = "newer")
source("utils.R")
fit <- readRDS("no_isol_gamma_alietal.rds")
tab1 <- fitted_params(fit)
samples_tost <- estimated_TOST_gamma(
  tab1, n = 1e4, extract(fit)
)
obs_data <- readRDS("alietal_data_clean.rds")
common_params <- readRDS("common_params.rds")

best_si <- estimated_SI(
  obs_data, samples_tost$TOST_bestpars,
  mixture = TRUE, isol = FALSE,
  tab1 = tab1, shape_inc = common_params$param_inc$shape,
  scale_inc = common_params$param_inc$scale
)

