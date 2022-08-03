fitted_params <- function(fit, digits = 2) {
  tab1 <- as.data.frame(rstan::summary(fit)$summary)
  fit <- rstan::extract(fit)

  best_idx <- which.max(fit[["lp__"]])
  best <- unlist(map(fit, function(x) x[[best_idx]]))
  tab1 <- add_column(tab1, best = best, .after = 1)

  round(tab1, digits)
}

# functions
# density function for NF distribution
dnf <- function(t, a = 0.5, b = 0.5, c = 0.1, tmax = 0) {
  numerator <- a + b
  growing <- a * exp(b * (t - tmax))
  failing <- b * exp(-a * (t - tmax))
  (numerator / (growing + failing))^c
}

# random sample from NF distribution
rnf <- function(n, taus = seq(-20, 40, 0.1), a = 0.5, b = 0.5, c = 0.1, tmax = 0) {
  p <- map_dbl(taus, function(t) dnf(t, a, b, c, tmax))
  sample(taus, n, replace = TRUE, prob = p)
}

estimated_TOST_phen <- function(tab1, taus, n) {
  # TOST
  best <- tab1$best
  names(best) <- rownames(tab1)
  TOST_bestpars <- rnf(
    n = n, taus = taus, a = best["a"], b = best["b"],
    c = best["c"], tmax = best["tmax"]
  )

  list(
    TOST_bestpars = TOST_bestpars
  )
}

estimated_TOST_sn <- function(tab1, n) {
  # TOST
  best <- tab1$best
  names(best) <- rownames(tab1)
  TOST_bestpars <- sn::rsn(
    n = n, xi = best["a"], omega = best["b"], alpha = best["c"]
  )

  list(
    TOST_bestpars = TOST_bestpars
  )
}

## Returns TOST for
## (a) parameters with maximum posterior likelihood
## (b) at the mean of the posterior distribution
## (c) a distribution of distributions of TOST ; 1
## distribution for each parameter in the posterior
## distributions of xi, omega, and alpha  sampled jointly.
## tab1 is the output from fitted_params function
estimated_TOST_gamma <- function(tab1, n, fit, offset = 20) {
  # TOST
  best <- tab1$best
  names(best) <- rownames(tab1)
  TOST_bestpars <- rgamma(
    n = n, shape = best[["alpha1"]], rate = best[["beta1"]]
  ) - offset

  list(
    TOST_bestpars = TOST_bestpars
  )
}


nice_model_name <- function(model) {
  suffix <- case_when(
    grepl("scenario3", model) ~ "BASELINE",
    grepl("scenario4", model) ~ "BASELINE + ISOL"
  )
  suffix1 <- case_when(
    grepl("mixture", model) ~ " + MIX",
    TRUE ~ ""
  )
  suffix2 <- case_when(
    grepl("recall", model) ~ " + RECALL",
    TRUE ~ ""
  )
  paste0(suffix, suffix1, suffix2)
}

### extract fitted parameters
fitted_params <- function(fit, digits = 2) {
  tab1 <- as.data.frame(rstan::summary(fit)$summary)
  fit <- rstan::extract(fit)

  best_idx <- which.max(fit[["lp__"]])
  best <- unlist(map(fit, function(x) x[[best_idx]]))
  tab1 <- add_column(tab1, best = best, .after = 1)

  round(tab1, digits)
}


unconditional_si <- function(inf_samples, shape_inc, scale_inc) {
  inc_samples <- rgamma(
    shape = shape_inc,
    scale = scale_inc,
    n = length(inf_samples)
  )
  list(
    inf = inf_samples,
    si = inc_samples + inf_samples
  )
}

## Return a named list where the names are the nus
## in the observe data set. This is so that
## (a) we can use s3 + recall etc i.e. make use of
## nu even when we don't condition on isolation,
## and (b) we have a consistent output format
## across all functions
unconditional_si_all <- function(obs, inf_samples, shape_inc, scale_inc) {

  x <- unconditional_si(inf_samples, shape_inc, scale_inc)
  nu_freq <- janitor::tabyl(obs$nu)
  colnames(nu_freq)[1] <- "nu"
  si_given_nu <- map(nu_freq$nu, function(nu) x)
  names(si_given_nu) <- nu_freq$nu
  si_given_nu
}

## un_si is the output of unconditional_si
##
conditional_si <- function(un_si, nu, nsim) {
  inf_samples <- un_si[[1]]
  idx <- which(inf_samples <= nu)
  filtered <- inf_samples[idx]
  filtered_si <- un_si[[2]][idx]

  if (length(filtered) == 0) return(NULL)
  inf_samples <- sample(filtered, nsim, replace = TRUE)
  si <- sample(filtered_si, nsim, replace = TRUE)
  list(inf = inf_samples, si = si)
}

## obs is observed data with column nu
## processed_fit is the output of process_beta_fit
conditional_si_all <- function(un_si, nsim) {
  si_given_nu <- imap(
    un_si, function(si, nu) {
      conditional_si(si, as.numeric(nu), nsim = nsim)
    }
  )
  keep(si_given_nu, function(x) !is.null(x[["si"]]))
}

## si_given_nu is the output of beta_fit_si_all, a list of conditional
## SIs for each nu in the data
## Returns a named list
conditional_si_pooled <- function(obs, si_given_nu, nsim) {
  ## Under some -ve nus, we have 0 possible SIs, we wil filter them out
  ## here.
  si_given_nu <- map(si_given_nu, function(x) x[["si"]])

  ## Renormalise weights over this smaller set of nus
  nu_freq <- janitor::tabyl(obs$nu[obs$nu %in% as.numeric(names(si_given_nu))])
  ## Number of samples to be drawn is nsim * percent
  nu_freq$nu_weights <- ceiling(
    nu_freq$percent * nsim
  )
  ## This can give you more or less than nsim
  ## SIs in all. You can resample when pooling.
  map2(
    si_given_nu, nu_freq$nu_weights, function(x, size) {
      sample(x, size, replace = TRUE)
    }
  )
}

leaky_si <- function(un_si, nu, pleak, nsim) {
  inf_samples <- un_si[[1]]

  ## Split inf_times into before and
  ## after nu. And then from the before
  ## group, select 1 - pleak samples
  ## and from the after group, select
  ## pleak samples
  before_nu <- which(inf_samples <= nu)
  after_nu <- which(inf_samples > nu)
  ## Do it this way rather than sampling the vector
  ## directly so that we can use the index to select
  ## SIs
  leaky <- sample(
    length(after_nu), size = floor(nsim * pleak),
    replace = TRUE
  )
  not_leaky <- sample(
    length(before_nu), size = floor(nsim * (1 - pleak)),
    replace = TRUE
  )
  filtered <- c(
    inf_samples[before_nu][not_leaky],
    inf_samples[after_nu][leaky]
  )
  filtered_si <- c(
    un_si[[2]][before_nu][not_leaky],
    un_si[[2]][after_nu][leaky]
  )
  if (length(filtered) == 0) return(NULL)
  ## TODO Sample these again to ensure length is
  ## nsim
  list(inf = inf_samples, si = filtered_si)

}

leaky_si_all <- function(un_si, pleak, nsim) {
  si_given_nu <- imap(
    un_si, function(si, nu) {
      leaky_si(si, as.numeric(nu), pleak, nsim = nsim)
    }
  )
  keep(si_given_nu, function(x) !is.null(x[["si"]]))
}


## predicted observed SIs (under assumed biases)
## currently assumes recall and isolation biases only affect valid SIs
##
## This function will then apply relevant biases to
estimated_SI <- function(obs, inf_times, mixture,
                         isol, tab1, nsim = 1e4, tmin = -20,
                         shape_inc = 5.807, scale_inc = 0.948) {


  ## First simulate unconditional SI and then apply biases
  ## This returns a named list where the names are the distinct
  ## nus in obs, EVEN THOUGH no conditioning is done on nu at this
  ## point. All elements of this list are the same.
  un_si <- unconditional_si_all(obs, inf_times, shape_inc, scale_inc)
  ## Now pool. Pooling involves sampling from
  ## nu-conditional SI. However, in case of unconditional
  ## SIs, it makes no difference because all elements of the
  ## list are identical.
  pooled_unc <- conditional_si_pooled(obs, un_si, nsim)
  pooled_unc <- unname(unlist(pooled_unc))

  # with isolation bias
  if (isol) {
    with_iso <- conditional_si_all(un_si, nsim)
    pooled_c <- conditional_si_pooled(obs, with_iso, nsim)
    pooled_c <- unname(unlist(pooled_c))
  }

  ## Mixture
  invalid_si <- runif(nsim, tmin, 40)
  pinvalid <- ifelse(mixture, tab1["pinvalid", "best"], 0)
  toss <- runif(nsim)
  valid <- which(toss >= pinvalid)

  mixed_unc <- c(pooled_unc[valid], invalid_si[-valid])
  if (isol) mixed_c <- c(pooled_c[valid], invalid_si[-valid])

  out <- list(
    unconditional = NULL,
    mixed_unc = NULL,
    mixed_c = NULL
  )
  ## Across all scenarions,
  ## we want the unconditional SI
  ## distribution as the true underlying
  ## distribution.
  out[["unconditional"]] <- pooled_unc

  if (isol) {
    ## Isolation is only going to fitted
    ## with mixture, so that is the only
    ## thing we care about.
    out[["mixed_c"]] <- mixed_c
  } else {
    ## When not accounting for isolation
    ## we only consider mixture model
    out[["mixed_unc"]] <- mixed_unc
  }
  out
}

