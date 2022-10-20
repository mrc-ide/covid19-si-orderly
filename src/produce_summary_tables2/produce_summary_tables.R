## orderly::orderly_develop_start(use_draft = 'newer')
infiles <- c(
  `Without Isolation` = "no_isol_gamma_alietal_si.rds",
  `Without Isolation (Pairs)` = "no_isol_gamma_alietal_si_pairs.rds",
  `With Isolation` = "isol_gamma_alietal_si.rds",
  `With Isolation (Pairs)` = "isol_gamma_alietal_si_pairs.rds",
  `Part Isolation` = "part_isol_gamma_alietal_si.rds",
  `Part Isolation (Pairs)` = "part_isol_gamma_alietal_si_pairs.rds"
)


## Given a vector, returns prettily formatted string
pretty_summary <- function(vec, probs = c(0.025, 0.5, 0.975), digits = 1) {
  out <- quantile(vec, probs = probs)
  out <- lapply(
    out, function(x) format(round(x, digits), nsmall = digits)
  )
  glue(
    "{out[[2]]} ({out[[1]]},{out[[3]]})"
  )
}

produce_si_table <- function(si) {
  out <- map(si, function(x) {
    if (length(x) == 0) NULL
    else data.frame(median_95_cri = pretty_summary(x))
  })

  out <- keep(out, ~ ! is.null(.))

  bind_rows(out, .id = "SI Type")
}

si_summary <- map(infiles, readRDS) %>%
  map_dfr(produce_si_table, .id = "Model") %>%
  spread(`SI Type`, median_95_cri)


saveRDS(si_summary, "gamma_fits_si.rds")


si_summary <- gsub("gamma", "sn", infiles) %>%
  map(readRDS) %>%
  map_dfr(produce_si_table, .id = "Model") %>%
  spread(`SI Type`, median_95_cri)


saveRDS(si_summary, "sn_fits_si.rds")


si_summary <- gsub("gamma", "phen", infiles) %>%
  map(readRDS) %>%
  map_dfr(produce_si_table, .id = "Model") %>%
  spread(`SI Type`, median_95_cri)


saveRDS(si_summary, "phen_fits_si.rds")
