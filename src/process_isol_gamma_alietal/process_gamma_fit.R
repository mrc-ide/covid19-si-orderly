## orderly::orderly_develop_start(use_draft = "newer")
source("utils.R")
fit <- readRDS("isol_gamma_alietal.rds")
## tab1 <- fitted_params(fit)
## samples_tost <- estimated_TOST_gamma(
##   tab1, n = 1e4, extract(fit)
## )

tost_pp <- postpred_TOST_gamma(fit, 1e4)
summary_tost <- tidy(summary(tost_pp[[1]]))

obs_data <- readRDS("alietal_data_clean.rds")
common_params <- readRDS("common_params.rds")

best_si <- estimated_SI(
  obs_data, tost_pp[[1]],
  mixture = TRUE, isol = TRUE,
  tab1 = tab1, shape_inc = common_params$param_inc$shape,
  scale_inc = common_params$param_inc$scale
)

saveRDS(best_si, "best_si.rds")
saveRDS(summary_tost, "summary_tost.rds")
