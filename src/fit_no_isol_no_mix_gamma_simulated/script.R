## orderly::orderly_develop_start(use_draft = "newer")
params <- readRDS("params.rds")
param_grid <- readRDS("param_grid.rds")
simulated_data <- readRDS("simulated_data.rds")
max_shed <- readRDS("max_shed.rds")
max_valid_si <- 40

width <- 0.1

pwalk(
  list(
    params_inc = param_grid$params_inc,
    offset = param_grid$offset,
    sim_data = simulated_data,
    index1 = seq_along(param_grid$params_inc)
  ),
  function(params_inc, offset, sim_data, index1) {
    param_inc <- params[[params_inc]]
    param_inc <- gamma_mucv2shapescale(
      mu = param_inc$mean_inc,
      cv = param_inc$sd_inc / param_inc$mean_inc
    )
    si_vec <- seq(offset + 0.5, max_valid_si, 1)
    out <- imap(sim_data, function(x, index2) {
      x <- x[x$type == "valid filtered", ]
      x <- arrange(x, nu)
      fit <- stan(
        file = "scenario3a_gamma.stan",
        data = list(
          N = length(sim_data$si),
          si = sim_data$si,
          nu = sim_data$nu,
          max_shed = max_shed,
          offset1 = offset,
          alpha2 = param_inc[["shape"]],
          beta2 = 1 / param_inc[["scale"]],
          max_valid_si = max_valid_si,
          min_valid_si = offset,
          width = width,
          M = length(si_vec),
          si_vec = si_vec,
          first_valid_nu = 1
      ),
      chains = 2, iter = 2000,
      seed = 42,
      verbose = FALSE
      )
      ## control = list(adapt_delta = 0.99)
    }
    )
    outfile <- glue(
      "param_{index1}_sim_{index2}.rds"
    )
    saveRDS(out, outfile)
  }
)



