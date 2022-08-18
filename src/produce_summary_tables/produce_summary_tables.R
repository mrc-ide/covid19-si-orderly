## orderly::orderly_develop_start(use_draft = 'newer')
infiles <- c(
  `Without Isolation` = "no_isol_gamma_alietal.rds",
  `With Isolation` = "no_isol_gamma_alietal.rds"
)

fits <- map(infiles, readRDS)
params <- map(fits, ~ summary(.)$summary)
params <- map(
  params, function(x) {
    x <- data.frame(x)
    x <- rownames_to_column(x, "parameter")
    x
  }
)
