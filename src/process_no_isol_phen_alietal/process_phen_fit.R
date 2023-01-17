## orderly::orderly_develop_start(use_draft = "newer")
source("utils.R")
fit <- readRDS("no_isol_phen_alietal.rds")

tost_pp <- postpred_TOST_phen(fit, taus = seq(-20, 40, 0.1), n = 1e4)

summary_tost <- tidy(summary(tost_pp[[1]]))

obs_data <- readRDS("alietal_data_clean.rds")
common_params <- readRDS("common_params.rds")

best_si <- estimated_SI(
  obs_data, tost_pp[[1]],
  mixture = TRUE, isol = FALSE,
  tab1 = tab1, shape_inc = common_params$param_inc$shape,
  scale_inc = common_params$param_inc$scale
)

saveRDS(best_si, "best_si.rds")
saveRDS(summary_tost, "summary_tost.rds")
