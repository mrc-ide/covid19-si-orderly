##orderly::orderly_develop_start()
alietal <- readRDS("alietal_data_clean.rds")
alietal$si <- as.numeric(alietal$si, units="days")
## Incubation period distribution
## https://pubmed.ncbi.nlm.nih.gov/32150748/
param_inc <- list(
  shape = 5.807, scale = 0.948
)

si_vec <- seq(-20, 40, 1)

fit <- stan(
  file = "scenario3a_mixture_gamma.stan",
  data = list(
    N = length(alietal$si),
    si = alietal$si,
    max_shed = 21,
    offset = 3,
    alpha2 = param_inc$shape,
    beta2 = param_inc$scale,
    max_invalid_si = 40,
    min_invalid_si = -20,
    width = 0.5,
    M = length(si_vec),
    si_vec = si_vec
  ),
  seed = 123,
  chains = 2,
  verbose = TRUE
)

saveRDS(fit, "no_isol_gamma_alietal.rds")
