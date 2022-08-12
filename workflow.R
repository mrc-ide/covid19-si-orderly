library(orderly)
use_draft <- 'newer'

orderly_run('fit_no_isol_gamma_alietal_pairs', use_draft = use_draft)

orderly_run('process_no_isol_gamma_alietal_pairs', use_draft = use_draft)
orderly_run('plot_no_isol_gamma_alietal_pairs', use_draft = use_draft)

orderly_run('process_isol_sn_alietal_pairs', use_draft = use_draft)
orderly_run('plot_isol_sn_alietal_pairs', use_draft = use_draft)

orderly_run('process_isol_phen_alietal_pairs', use_draft = use_draft)
orderly_run('plot_isol_phen_alietal_pairs', use_draft = use_draft)


orderly_run('process_isol_gamma_alietal_pairs', use_draft = use_draft)
orderly_run('plot_isol_gamma_alietal_pairs', use_draft = use_draft)

orderly_run('process_part_isol_gamma_alietal', use_draft = use_draft)
orderly_run('plot_part_isol_gamma_alietal', use_draft = use_draft)
