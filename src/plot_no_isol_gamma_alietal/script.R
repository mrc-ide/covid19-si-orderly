## orderly::orderly_develop_start(use_draft = "newer")
obs_data <- readRDS("alietal_data_clean.rds")
best_si <- readRDS("best_si.rds")

plot_both_SIs <- function(SI1, SI2, data) {

  p <- ggplot() +
    geom_histogram(
      data = data, aes(si, y = ..density.., fill = "gray77"),
      alpha = 0.8,
      binwidth = 1
    ) +
    geom_density(aes(SI1, fill = "red"),
      alpha = 0.3, colour = NA
    ) +
    geom_density(aes(SI2, fill = "blue"),
      alpha = 0.3, colour = NA
    ) +
    geom_vline(
      xintercept = mean(data$si), col = "gray77", linetype = "dashed"
    ) +
    geom_vline(
      xintercept = mean(SI1), col = "red", linetype = "dashed"
    ) +
    geom_vline(
      xintercept = mean(SI2), col = "blue", linetype = "dashed"
    ) +
    scale_fill_identity(
      guide = "legend",
      labels = c("Data", "Posterior - True", "Posterior - Expected"),
      breaks = c("gray77", "red", "blue")
    ) +
    theme_minimal() +
    xlab("Serial Interval (days)") +
    ## annotate(
    ##   geom = "text", x = mean(SI1), y = 0.055, color = "red",
    ##   label = paste(" Mean:", format(round(mean(SI1), 1), nsmall = 1), " days", sep = ""),
    ##   hjust = -0.1, size = 5
    ## ) +
    ## annotate(
    ##   geom = "text", x = mean(SI2), y = 0.095, color = "blue",
    ##   label = paste(" Mean:", format(round(mean(SI2), 1), nsmall = 1), " days", sep = ""),
    ##   hjust = -0.1, size = 5
    ## ) +
    theme(legend.title = element_blank()) +
    guides(
      fill = guide_legend(override.aes = list(alpha = c(0.8, 0.3, 0.3)))
    ) +
    theme(legend.position = "top")
  p
}

summary_tab <- data.frame(
  SI = c("True", "Expected"),
  Mean = c(mean(best_si$unconditional),
           mean(best_si$mixed_unc)),
  SD = c(sd(best_si$unconditional),
         sd(best_si$mixed_unc))
)

make_pretty <- function(x) {
  format(round(x, 1), nsmall = 1)
}

summary_tab <- mutate_if(
  summary_tab, is.numeric,
  make_pretty
)

dfnpc <- tibble(x = 0.95, y = 0.95, tb = list(summary_tab))


p <- plot_both_SIs(
  best_si$unconditional,
  best_si$mixed_unc,
  obs_data
)


p <- p +
  theme_minimal() +
  geom_table_npc(
    data = dfnpc, aes(npcx = x, npcy = y, label = tb),
    table.theme = ttheme_gtminimal()
  )

ggsave("no_isol_gamma_alietal.png", p)
