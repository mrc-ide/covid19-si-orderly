## orderly::orderly_develop_start(use_draft = 'newer')
pretty_table <- function(x) {
  x$`Realised SI` <- case_when(
    is.na(x$mixed_c) ~ x$mixed_unc,
    is.na(x$mixed_unc) ~ x$mixed_c
  )
  x <- rename(x, `True SI` = "unconditional", TOST = "Median (95% CrI)")
  select(
    x, Model, `Realised SI`, `True SI`, TOST
  )
}

roworder <- c("With Isolation", "Without Isolation",
              "With Isolation (Pairs)",  "Without Isolation (Pairs)",
              "Part Isolation", "Part Isolation (Pairs)")

infiles <- c("gamma_si.rds", "gamma_tost.rds")
tabs <- map(infiles, readRDS)
gamma_tab <- left_join(tabs[[1]], tabs[[2]]) |> pretty_table()
gamma_tab <- gamma_tab[match(roworder, gamma_tab$Model), ]

gtsave(gt(gamma_tab), "gamma_summary.docx")


infiles <- c("sn_si.rds", "sn_tost.rds")
tabs <- map(infiles, readRDS)
sn_tab <- left_join(tabs[[1]], tabs[[2]]) |> pretty_table()
sn_tab <- sn_tab[match(roworder, sn_tab$Model), ]

gtsave(gt(sn_tab), "sn_summary.docx")


infiles <- c("phen_si.rds", "phen_tost.rds")
tabs <- map(infiles, readRDS)
phen_tab <- left_join(tabs[[1]], tabs[[2]]) |> pretty_table()
phen_tab <- phen_tab[match(roworder, phen_tab$Model), ]

gtsave(gt(phen_tab), "phen_summary.docx")
