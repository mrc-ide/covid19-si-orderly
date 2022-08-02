##orderly::orderly_develop_start()

fits <- unzip("s3_stanfits.zip")
param_grid <- readRDS("param_grid.rds")
params <- readRDS("params.rds")
simulated_data <- readRDS("simulated_data.rds")
fits <- map(
  seq_along(param_grid$params_inc),
  function(index1) {
    map(
      seq_along(simulated_data[[1]]),
      function(index2) {
        infile <- glue(
          "stanfits/param_{index1}_sim_{index2}_fit.rds"
        )
        fit <- readRDS(infile)
      }
    )
  }
)

table1 <- map2(
  simulated_data,
  fits,
  function(sim_data, fits) {
    map2(
      sim_data, fits, function(x, fit) {
        list(
          true_si_mean = mean(x$si),
          realised_si_mean = mean(x$si[x$type == "valid filtered"]),
          tab1 = fitted_params(fit)
        )
      }
    )
  }
)

samples_tost <- map2(
  table1, fits,
  function(tab1, fit) {
    map2(
      tab1, fit, function(x, y) {
        estimated_TOST_gamma(
          x$tab1, n = 1e3, rstan::extract(y), offset = 0
        )
      }
    )
  }
)


best_si <- map2(
  samples_tost,
  param_grid$params_inc,
  function(tost, param_inc) {
    param_inc <- params[[param_inc]]
    param_inc <- gamma_mucv2shapescale(
      mu = param_inc$mean_inc,
      cv = param_inc$sd_inc / param_inc$mean_inc
    )
    inc <- rgamma(1e3, shape = param_inc$shape, scale = param_inc$scale)
    map(
      tost, function(x) {
        out <- x$TOST_bestpars + inc
        data.frame(
          mean_si = mean(out), sd_si = sd(out),
          median_si = quantile(out, 0.5),
          low_si = quantile(out, 0.025),
          high_si = quantile(out, 0.975)
        )
      }
    )
  }
)

obs <- imap_dfr(
  table1,
  function(tab1, index) {
    map_dfr(tab1, function(x) {
      data.frame(
        realised_si_mean = x$realised_si_mean,
        true_si_mean = x$true_si_mean
      )
    }, .id = "sim")
  }, .id = "param_set"
)

fitted <- imap_dfr(
  best_si, function(x, index) bind_rows(x, .id = "sim"),
  .id = "param_set"
)

saveRDS(obs, "simulated_si.rds")
saveRDS(fitted, "fitted_si_s3.rds")

## out <- cbind(obs, fitted[, -c(1, 2)])
## out$param_set <- glue("Parameter set {out$param_set}")
## out <- pivot_longer(
##   out,
##   cols = c("realised_si_mean", "mean_si")
## )

## x <- group_by(out, param_set) %>%
##   summarise(true_si_mean = mean(true_si_mean)) %>%
##   ungroup()

## p <- ggplot() +
##   geom_boxplot(
##     data = out,
##     aes(param_set, value, fill = name), position = position_dodge(width = 1),
##     alpha = 0.5
##   ) +
##   geom_point(
##     data = x, aes(param_set, true_si_mean), shape = 4, size = 3,
##     inherit.aes = FALSE
##   ) +
##   ylab("Mean serial interval") +
##   scale_fill_manual(
##     values = c(realised_si_mean = "#999999", mean_si = "#56B4E9"),
##     labels = c("Realised", "Fitted"),
##     name = "Mean serial interval"
##   ) +
##   theme_classic() +
##     theme(
##       legend.position = "top",
##       legend.title = element_blank(),
##       axis.title.x = element_blank()
##     )

## ggsave("simulated_s3.png", p)
