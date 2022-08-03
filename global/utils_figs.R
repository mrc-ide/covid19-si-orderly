plot_both_SIs <- function(unconditional, conditional, data) {

  p <- ggplot() +
    geom_histogram(
      data = data, aes(si, y = ..density.., fill = "gray77"),
      alpha = 0.8,
      binwidth = 1
    ) +
    geom_density(aes(unconditional, fill = "red"),
      alpha = 0.3, colour = NA
    ) +
    geom_density(aes(conditional, fill = "blue"),
      alpha = 0.3, colour = NA
    ) +
    geom_vline(
      xintercept = mean(data$si), col = "gray77", linetype = "dashed"
    ) +
    geom_vline(
      xintercept = mean(unconditional), col = "red", linetype = "dashed"
    ) +
    geom_vline(
      xintercept = mean(conditional), col = "blue", linetype = "dashed"
    ) +
    scale_fill_identity(
      guide = "legend",
      labels = c("Data", "Posterior - True", "Posterior - Expected"),
      breaks = c("gray77", "red", "blue")
    ) +
    theme_minimal() +
    xlab("Serial Interval (days)") +
    ## annotate(
    ##   geom = "text", x = mean(unconditional), y = 0.055, color = "red",
    ##   label = paste(" Mean:", format(round(mean(unconditional), 1), nsmall = 1), " days", sep = ""),
    ##   hjust = -0.1, size = 5
    ## ) +
    ## annotate(
    ##   geom = "text", x = mean(conditional), y = 0.095, color = "blue",
    ##   label = paste(" Mean:", format(round(mean(conditional), 1), nsmall = 1), " days", sep = ""),
    ##   hjust = -0.1, size = 5
    ## ) +
    theme(legend.title = element_blank()) +
    guides(
      fill = guide_legend(override.aes = list(alpha = c(0.8, 0.3, 0.3)))
    ) +
    theme(legend.position = "top", legend.title = element_blank())
  p
}

make_pretty <- function(x) {
  format(round(x, 1), nsmall = 1)
}
