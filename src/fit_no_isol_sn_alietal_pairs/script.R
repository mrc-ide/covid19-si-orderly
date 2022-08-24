##orderly::orderly_develop_start()
alietal <- readRDS("alietal_data_clean.rds")
alietal$si <- as.numeric(alietal$si, units="days")
common_params <- readRDS("common_params.rds")
fit <- stan(
  file = "scenario3a_mixture_skew_normal.stan",
  data = list(
    N = length(alietal$si),
    si = alietal$si,
    max_shed = common_params$max_shed,
    alpha2 = common_params$param_inc$shape,
    beta2 = common_params$param_inc$scale,
    max_invalid_si = common_params$max_invalid_si,
    min_invalid_si = common_params$min_invalid_si,
    width = common_params$width,
    M = length(common_params$si_vec),
    si_vec = common_params$si_vec
  ),
  seed = common_params$seed,
  chains = common_params$chains, iter = common_params$iter,
  verbose = TRUE
)

saveRDS(fit, "no_isol_sn_alietal.rds")
