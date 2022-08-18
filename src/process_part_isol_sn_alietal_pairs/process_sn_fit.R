## orderly::orderly_develop_start(use_draft = "newer")
source("utils.R")
fit <- readRDS("part_isol_sn_alietal_pairs.rds")
ggmcmc(ggs(fit))
tab1 <- fitted_params(fit)
samples_tost <- estimated_TOST_sn(
  tab1, n = 1e4
)
summary_tost <- tidy(summary(samples_tost[[1]]))

obs_data <- readRDS("alietal_data_clean.rds")
common_params <- readRDS("common_params.rds")

best_si <- estimated_SI(
  obs_data, samples_tost$TOST_bestpars,
  mixture = TRUE, isol = TRUE,
  tab1 = tab1, shape_inc = common_params$param_inc$shape,
  scale_inc = common_params$param_inc$scale
)

saveRDS(best_si, "best_si.rds")
saveRDS(summary_tost, "summary_tost.rds")
