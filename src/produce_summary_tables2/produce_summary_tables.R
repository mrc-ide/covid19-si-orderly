## orderly::orderly_develop_start(use_draft = 'newer')
infiles <- c(
  `Without Isolation` = "no_isol_gamma_alietal_si.rds",
  `Without Isolation (Pairs)` = "no_isol_gamma_alietal_si_pairs.rds",
  `With Isolation` = "isol_gamma_alietal_si.rds",
  `With Isolation (Pairs)` = "isol_gamma_alietal_si_pairs.rds",
  `Part Isolation` = "part_isol_gamma_alietal_si.rds",
  `Part Isolation (Pairs)` = "part_isol_gamma_alietal_si_pairs.rds"
)

produce_si_table <- function(si) {
  out <- map(si, function(x) {
    if (length(x) == 0) NULL
    else tidy(summary(x))
  })

  out <- keep(out, ~ ! is.null(.))

  bind_rows(out, .id = "SI Type")
}

si_summary <- map(infiles, readRDS) %>%
  map_dfr(produce_si_table, .id = "Model")


saveRDS(si_summary, "gamma_fits_si.rds")


si_summary <- gsub("gamma", "sn", infiles) %>%
  map(readRDS) %>%
  map_dfr(produce_si_table, .id = "Model")


saveRDS(si_summary, "sn_fits_si.rds")


si_summary <- gsub("gamma", "phen", infiles) %>%
  map(readRDS) %>%
  map_dfr(produce_si_table, .id = "Model")


saveRDS(si_summary, "phen_fits_si.rds")

tost_summary <- gsub("gamma", "sn", infiles) %>%
  gsub("si", "tost", .) %>%
  map_dfr(readRDS, .id = "Model")

saveRDS(tost_summary, "sn_fits_tost.rds")


tost_summary <- gsub("gamma", "phen", infiles) %>%
  gsub("si", "tost", .) %>%
  map_dfr(readRDS, .id = "Model")

saveRDS(tost_summary, "phen_fits_tost.rds")
