## orderly::orderly_develop_start()
map_into_interval <- function(x_1, y_1, x_2, y_2) {
  slope <- (x_2 - y_2) / (x_1 - y_1)
  intercept <- x_2 - x_1 * slope
  f <- function(x) slope * x + intercept
  f
}

simulate_si <- function(tost, params_inc, params_nu, nsim) {
  t_2 <- stats::rgamma(
    n = nsim, shape = params_inc$shape, rate = 1 / params_inc$scale
  )

  out <- data.frame(t_1 = tost, t_2 = t_2, si = tost + t_2)

  if (! is.null(params_nu)) {
    out$nu <- stats::rgamma(
      n = nsim, shape = params_nu$shape,
      scale = params_nu$scale
    )
  } else out$nu <- NULL
  out
}

## offset is positive so to go left you subtract.
simulate_tost_gamma <- function(params_inf, min_inf, max_inf, nsim = 50) {
  stats::rgamma(
    n = nsim, scale = params_inf$scale, shape = params_inf$shape
  ) - min_inf

}

## Simulate lots more than you need because lots
## will be filtered out
dir.create("figures")
nsim <- 50000
nsim_post_filter <- 100
## For each parameter combination simulate
## multiple data sets
ndatasets <- 50
max_shed <- 21
params <- list(
  inf_par1 = list(mean_inf = 4, sd_inf = 2), ## short
  inf_par2 = list(mean_inf = 8, sd_inf = 2), ## long
  inc_par1 = list(mean_inc = 2, sd_inc = 1), ## short
  inc_par2 = list(mean_inc = 6, sd_inc = 2), ## long
  iso_par1 = list(mean_iso = 2, sd_iso = 1), ## short
  iso_par2 = list(mean_iso = 6, sd_iso = 3), ## long
  offset1 = 0,
  offset2 = 1,
  offset3 = 2,
  pinvalid1 = 0,
  pinvalid2 = 0.05,
  pinvalid3 = 0.2
)

param_grid <- expand.grid(
  params_inf = c("inf_par1", "inf_par2"),
  params_inc = c("inc_par1", "inc_par2"),
  params_iso = c("iso_par1"), offset = params$offset1,
  stringsAsFactors = FALSE
)

saveRDS(params, "params.rds")
saveRDS(param_grid, "param_grid.rds")
saveRDS(max_shed, "max_shed.rds")

valid_unfiltered <- pmap(
  param_grid,
  function(params_inf, params_inc, params_iso, offset) {
    param_inf <- params[[params_inf]]
    param_inc <- params[[params_inc]]
    param_iso <- params[[params_iso]]
    param_inf <- gamma_mucv2shapescale(
      param_inf$mean_inf, param_inf$sd_inf / param_inf$mean_inf
    )
    param_inc <- gamma_mucv2shapescale(
      mu = param_inc$mean_inc,
      cv = param_inc$sd_inc / param_inc$mean_inc
    )
    param_iso <- gamma_mucv2shapescale(
      mu = param_iso$mean_iso,
      cv = param_iso$sd_iso/ param_iso$mean_iso
    )
    map(
      seq_len(ndatasets), function(i) {
        valid_t1 <- simulate_tost_gamma(
          param_inf, offset, max_shed, nsim
        )
        simulate_si(
          valid_t1, param_inc, param_iso, nsim
        )
      }
    )
  }
)

valid_filtered <- map_depth(
  valid_unfiltered, 2, function(x) {
    out <- x[x$t_1 <= x$nu, ]
    if (nrow(out) < nsim_post_filter) {
      warning(
        "Number of filtered rows is ",
        nrow(out)
      )
      out <- out
    } else {
      idx <- sample(
        1:nrow(out), size = nsim_post_filter
      )
      out <- out[idx, ]
    }
    out
 }
)


all_data <- pmap(
  list(
    x = valid_unfiltered,
    y = valid_filtered
  ),
  function(x, y) {
    map2(x, y, function(unf, fil) {
      rbind(
        data.frame(
          type = "all valid", si = unf$si, nu = unf$nu
        ),
        data.frame(
          type = "valid filtered", si = fil$si, nu = fil$nu
        )
      )
    })
  }
)
saveRDS(all_data, "simulated_data.rds")

palette <- c(
 "all valid" = "#999999",
 "valid filtered" = "#E69F00"
)

iwalk(
  all_data,
  function(sims, index1) {
    iwalk(
      sims, function(df, index2) {
        p <- ggplot(df, aes(si, fill = type)) +
          geom_density(alpha = 0.4, col = NA, adjust = 2) +
          scale_fill_manual(
            values = palette,
            breaks = c("all valid", "valid filtered"),
            labels = c("All possible SIs",
                       "Infection of infectee before isolation of primary"),
            guide = guide_legend(nrow = 2)
          ) +
          theme_minimal() +
          theme(legend.position = "top", legend.title = element_blank()) +
          xlab("Serial Interval") +
          ylab("Density")
        ggsave(
          glue("figures/{index1}_{index2}.pdf"), p
        )
      }
    )
  }
)
