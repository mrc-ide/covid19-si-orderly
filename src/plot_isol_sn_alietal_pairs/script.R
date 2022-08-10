## orderly::orderly_develop_start(use_draft = "newer")
source("utils_figs.R")
obs_data <- readRDS("alietal_data_clean.rds")
best_si <- readRDS("best_si.rds")


summary_tab <- data.frame(
  SI = c("True       ", "Expected"),
  Mean = c(mean(best_si$unconditional),
           mean(best_si$mixed_c)),
  SD = c(sd(best_si$unconditional),
         sd(best_si$mixed_c))
)


summary_tab <- mutate_if(
  summary_tab, is.numeric,
  make_pretty
)

dfnpc <- tibble(x = 0.95, y = 0.95, tb = list(summary_tab))


p <- plot_both_SIs(
  best_si$unconditional,
  best_si$mixed_c,
  obs_data
)


p <- p +
  geom_table_npc(
    data = dfnpc, aes(npcx = x, npcy = y, label = tb),
    table.theme = ttheme_gtminimal()
  ) + ggtitle("Serial interval estimated using discrete transmission pairs",
              subtitle = "Infectious period modelled as a Skew normal rv")

ggsave("isol_sn_alietal_pairs.png", p)
