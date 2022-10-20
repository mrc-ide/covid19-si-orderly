## orderly::orderly_develop_start(use_draft = "newer")
source("utils.R")
tost_summary <- function(samples_tost) {
    ## Prepare table for annotation
  df <- data.frame(
    Model = names(samples_tost),
    median = c(quantile(samples_tost[[1]][[1]], 0.5),
               quantile(samples_tost[[2]][[1]], 0.5)),
    low = c(quantile(samples_tost[[1]][[1]], 0.025),
            quantile(samples_tost[[2]][[1]], 0.025)),
    high = c(quantile(samples_tost[[1]][[1]], 0.975),
             quantile(samples_tost[[2]][[1]], 0.975)))

  f <- function(x) format(round(x, 1), nsmall = 1)

  df$`Median (95% CrI)` <- glue(
    "{f(df$median)} ({f(df$low)} - {f(df$high)})"
  )

  df
}

tost_summary_figure <- function(tost, distr) {
  iwalk(tost, function(samples_tost, model) {
    sub_models <- names(samples_tost)
    fill1 <- palette[[model]][[sub_models[1]]]
    fill2 <- palette[[model]][[sub_models[2]]]

    df <- tost_summary(samples_tost)
    medians <- df$median

    p <- ggplot() +
      geom_density(
        aes(samples_tost[[1]][[1]], fill = fill1),
        alpha = alpha, col = NA
      ) +
      geom_density(
        aes(samples_tost[[2]][[1]], fill = fill2),
        alpha = alpha, col = NA
      ) +
      geom_vline(
        aes(xintercept = medians[[1]]), col = fill1, linetype = "dashed",
        size = linesize
      ) +
      geom_vline(
        aes(xintercept = medians[[2]]), col = fill2, linetype = "dashed",
        size = linesize
      ) +
      geom_table_npc(
        data = df, label = list(df), npcx = 0.95, npcy = 0.95,
        table.theme = ttheme_default(base_size = 8)
      ) +
      xlab("TOST (days)") +
      scale_fill_identity(
        labels = sub_models, breaks = c(fill1, fill2), guide = "legend"
      ) +
      theme_bw() +
      theme(
        text = element_text(size = 20 / .pt),
        legend.title = element_blank(),
        legend.position = "top"
      )

    ggsave(glue("{distr}_{model}_tost.png"), p)

  })

}

palette <- list(
  full_dataset = c(`Without Isolation` = "#E69F00",
                   `With Isolation` = "#56B4E9"),

  pairs= c(`Without Isolation (Pairs)` = "#009E73",
           `With Isolation (Pairs)` = "#CC79A7"),

  part_isol = c(`Part Isolation` = "#0072B2",
                `Part Isolation (Pairs)` = "#D55E00")
)


alpha <- 0.3
linesize <- 1.1
infiles <- list(
  full_dataset = c(`Without Isolation` = "no_isol_gamma_alietal.rds",
                   `With Isolation` = "isol_gamma_alietal.rds"),

  pairs= c(`Without Isolation (Pairs)` = "no_isol_gamma_alietal_pairs.rds",
           `With Isolation (Pairs)` = "isol_gamma_alietal_pairs.rds"),

  part_isol = c(`Part Isolation` = "part_isol_gamma_alietal.rds",
                `Part Isolation (Pairs)` = "part_isol_gamma_alietal_pairs.rds")
)

## Re-generate the TOST distribution here as otherwise it takes up
## too much space filling up my hard disk.
tost <- imap(
  infiles, function(files, model) {
    fits <- map(files, readRDS)
    tab1 <- map(fits, fitted_params)
    map(
      tab1, function(x) {
        estimated_TOST_gamma(
          x, n = 1e4, extract(x)
        )
      }
    )
  }
)

out <- tost_summary_figure(tost, "gamma")


tost <- map(infiles, ~ gsub("gamma", "sn", .)) %>%
  imap(function(files, model) {
    fits <- map(files, readRDS)
    tab1 <- map(fits, fitted_params)
    map(
      tab1, function(x) {
        estimated_TOST_sn(x, n = 1e4)
      }
    )
  }
)

tost_summary_figure(tost, "sn")


tost <- map(infiles, ~ gsub("gamma", "phen", .)) %>%
  imap(function(files, model) {
    fits <- map(files, readRDS)
    tab1 <- map(fits, fitted_params)
    map(
      tab1, function(x) {
        estimated_TOST_phen(x, taus = seq(-20, 40, 0.1), n = 1e4)
      }
    )
  }
)

tost_summary_figure(tost, "phen")
