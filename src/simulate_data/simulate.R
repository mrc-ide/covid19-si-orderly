map_into_interval <- function(x_1, y_1, x_2, y_2) {
  slope <- (x_2 - y_2) / (x_1 - y_1)
  intercept <- x_2 - x_1 * slope
  f <- function(x) slope * x + intercept
  f
}

simulate_si <- function(tost, params_inc, params_nu ) {
  t_2 <- stats::rgamma(
    n = nsim, shape = params_inc$shape, rate = 1 / params_inc$scale
  )

  out <- data.frame(t_1 = tost, t_2 = t_2, si = t_1 + t_2)

  if (! is.null(params_nu)) {
    out$nu <- stats::rgamma(
      n = nsim, shape = params_iso$shape, rate = 1 / params_iso$scale
      )
  } else out$nu <- NULL
  out
}

simulate_tost_gamma <- function(params_inf, min_inf, max_inf, nsim = 50) {

  t_1 <- stats::rgamma(
    n = nsim, shape1 = params_inf$shape1, shape2 = params_inf$shape2
  )
  ## Map the possible infection times into interval
  ## (min_inf, max_inf)
  f <- map_into_interval(0, 1, min_inf, max_inf)
  f(t_1)
}


nsim <- 5000
params <- list(
  inf_par1 = list(mean_inf = 2, sd_inf = 1), ## short
  inf_par2 = list(mean_inf = 8, sd_inf = 2), ## long
  inc_par1 = list(mean_inc = 2, sd_inc = 1), ## short
  inc_par2 = list(mean_inc = 6, sd_inc = 2), ## long
  iso_par1 = list(mean_iso = 2, sd_iso = 2), ## short
  iso_par2 = list(mean_iso = 10, sd_iso = 3), ## long
  offset1 = 0,
  offset2 = -1,
  offset3 = -2,
  pinvalid1 = 0,
  pinvalid2 = 0.05,
  pinvalid3 = 0.2,
)

param_grid <- expand.grid(
  params_inf = c("inf_par1", "inf_par2"),
  params_inc = c("inc_par1", "inc_par2"),
  params_iso = c("iso_par1"), offset = params$offset2,
  stringsAsFactors = FALSE
)

valid_unfiltered <- pmap(
  param_grid,
  function(params_inf, params_inc, params_iso, offset) {
    params_inf <- params[[params_inf]]
    params_inc <- params[[params_inc]]

    params_inf <- gamma_mucv2shapescale(
      params_inf$mean_inf, params_inf$sd_inf / params_inf$mean_inf
    )
    params_inc <- gamma_mucv2shapescale(
      mu = params_inc[[1]], cv = params_inc[[2]] / params_inc[[1]]
    )
    valid_t1 <- simulate_tost_gamma(params_inf, offset, max_shed, nsim)
    valid_si <- simulate_si(
      valid_t1, params_inc, params_iso, nsim
    )
    valid_si
  }
)

valid_filtered <- map(
  valid_unfiltered, function(simulated_si) {
    simulated_si[simulated_si$t_1 <= simulated_si$nu, ]
 }
)

invalid_remapped <- map(
  valid_unfiltered,
  function(valid) {
    df <- invalid_si
    max_si <- max(valid$si)
    df$si <- (max_si - min_si) * df$si + min_si
    df
  }
)

mixed <- pmap(
  list(
    valid = valid_filtered,
    invalid = invalid_remapped,
    params_pinv = param_grid$params_pinv
  ),
  function(valid, invalid, params_pinv) {
    pinvalid <- params[[params_pinv]]
    toss <- runif(nrow(valid))
    rbind(
      valid[toss > pinvalid , c("si", "nu")],
      invalid[toss <= pinvalid ,c("si", "nu")]
    )
  }
)


with_recall_bias <- map2(
  mixed,
  param_grid$params_recall,
  function(df, recall_true) {
    recall_true <- params[[recall_true]]
    df$p_si <- exp(abs(df$si - df$nu) * -recall_true)
    idx <- sample(
      nrow(df), nrow(df), replace = TRUE, prob = df$p_si
    )
    df[idx, ]
  }
)

all_data <- pmap(
  list(
    x = valid_unfiltered,
    y = valid_filtered,
    a = mixed,
    b = with_recall_bias
  ),
  function(x, y, a, b) {
    rbind(
      data.frame(type = "all valid", si = x$si, nu = x$nu),
      data.frame(type = "valid filtered", si = y$si, nu = y$nu),
      data.frame(type = "mixed", si = a$si, nu = a$nu),
      data.frame(type = "mixed and recall", si = b$si, nu = b$nu)
    )
  }
)

iwalk(
  all_data,
  function(df, i) saveRDS(df, glue::glue("data/simulated_{i}.rds"))
)

palette <- c(
 "all valid" = "#999999",
 "valid filtered" = "#E69F00",
  "mixed" = "#56B4E9",
  "mixed and recall" = "#CC79A7"
)

iwalk(
  all_data,
  function(df, index) {
    p <- ggplot(df, aes(si, fill = type)) +
      geom_density(alpha = 0.4, col = NA) +
      scale_fill_manual(
        values = palette,
        breaks = c("all valid", "valid filtered",
                   "mixed", "mixed and recall"),
        labels = c("All possible SIs",
                   "Infection of infectee before isolation of primary",
                   "Mixture distribution",
                   "Mixture distribution sampled with recall bias"),
        guide = guide_legend(nrow = 2)
      ) +
      theme_minimal() +
      theme(legend.position = "top", legend.title = element_blank()) +
      xlab("Serial Interval") +
      ylab("Density")

    ggsave(glue::glue("figures/{index}.pdf"), p)
  }
)

iwalk(
  all_data,
  function(df, index) {
    p <- ggplot(df, aes(nu, si, col = type)) +
      geom_point() +
      scale_color_manual(
        values = palette,
        breaks = c("all valid", "valid filtered",
                   "mixed", "mixed and recall"),
        labels = c("All possible SIs",
                   "Infection of infectee before isolation of primary",
                   "Mixture distribution",
                   "Mixture distribution sampled with recall bias"),
        guide = guide_legend(nrow = 2)
      ) +
      theme_minimal() +
      theme(legend.position = "top", legend.title = element_blank()) +
      xlab("Delay from symptom onset to isolation") +
      ylab("Serial Interval")

    ggsave(glue("figures/nu_vs_si_{index}.pdf"), p)
  }
)
