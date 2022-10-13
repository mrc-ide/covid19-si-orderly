## orderly::orderly_develop_start(use_draft = "newer")
source("utils.R")
fit <- readRDS("isol_gamma_alietal.rds")
tab1 <- fitted_params(fit)
samples_tost <- estimated_TOST_gamma(
  tab1, n = 1e4, extract(fit)
)

infiles <- list(
  full_dataset = c(`Without Isolation` = "no_isol_gamma_alietal.rds",
                   `With Isolation` = "isol_gamma_alietal.rds"),

  pairs= c(`Without Isolation (Pairs)` = "no_isol_gamma_alietal_pairs.rds",
           `With Isolation (Pairs)` = "isol_gamma_alietal_pairs.rds"),

  part_isol = c(`Part Isolation` = "part_isol_gamma_alietal.rds",
                `Part Isolation (Pairs)` = "part_isol_gamma_alietal_pairs.rds")
)


palette <- list(
  full_dataset = c(`Without Isolation` = "#E69F00",
                   `With Isolation` = "#56B4E9"),

  pairs= c(`Without Isolation (Pairs)` = "#009E73",
           `With Isolation (Pairs)` = "#CC79A7"),

  part_isol = c(`Part Isolation` = "#0072B2",
                `Part Isolation (Pairs)` = "#D55E00")
)
alpha <- 0.3
## Re-generate the TOST distribution here as otherwise it takes up
## too much space filling up my hard disk.
plots <- imap(
  infiles, function(files, model) {
    browser()
    fits <- map(files, readRDS)
    tab1 <- map(fits, fitted_params)
    samples_tost <- map(
      tab1, function(x) {
        estimated_TOST_gamma(
          x, n = 1e4, extract(x)
        )
      }
    )
    sub_models <- names(samples_tost)
    fill1 <- palette[[model]][[sub_models[1]]]
    fill2 <- palette[[model]][[sub_models[2]]]

    p <- ggplot() +
      geom_density(
        aes(samples_tost[[1]][[1]]), fill = fill1, alpha = alpha, col = NA
      ) +
      geom_density(
        aes(samples_tost[[2]][[1]]), fill = fill2, alpha = alpha, col = NA
      ) +
      theme_bw()

  })


