## Incubation period distribution
## https://pubmed.ncbi.nlm.nih.gov/32150748/

common_params <- list(
  param_inc = list(
    shape = 5.807, scale = 0.948
  ),
  offset = 20,
  max_shed = 21,
  max_invalid_si = 40,
  min_invalid_si = -20,
  si_vec = seq(-20, 40, 1),
  ## Step size of integration for Stan.
  width = 0.5,
  chains = 2,
  iter = 10000,
  seed = 123
)

saveRDS(common_params, "common_params.rds")
