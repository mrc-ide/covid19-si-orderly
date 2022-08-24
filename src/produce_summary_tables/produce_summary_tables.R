## orderly::orderly_develop_start(use_draft = 'newer')
infiles <- c(
  `Without Isolation` = "no_isol_gamma_alietal.rds",
  `Without Isolation (Pairs)` = "no_isol_gamma_alietal_pairs.rds",
  `With Isolation` = "isol_gamma_alietal.rds",
  `With Isolation (Pairs)` = "isol_gamma_alietal_pairs.rds",
  `Part Isolation` = "part_isol_gamma_alietal.rds",
  `Part Isolation (Pairs)` = "part_isol_gamma_alietal_pairs.rds"
)

produce_params_table <- function(x) {
  x <- data.frame(x)
  x <- rownames_to_column(x, "parameter")
  x
}

fits <- map(infiles, readRDS)
params <- map(fits, ~ summary(.)$summary)
params <- map_dfr(params, produce_params_table, .id = "Model")
saveRDS(params, "gamma_fits_params.rds")


fits <- gsub("gamma", "sn", infiles) %>%
  map(readRDS)
params <- map(fits, ~ summary(.)$summary)
params <- map_dfr(params, produce_params_table, .id = "Model")
saveRDS(params, "sn_fits_params.rds")

fits <- gsub("gamma", "phen", infiles) %>%
  map(readRDS)
params <- map(fits, ~ summary(.)$summary)
params <- map_dfr(params, produce_params_table, .id = "Model")
saveRDS(params, "phen_fits_params.rds")
