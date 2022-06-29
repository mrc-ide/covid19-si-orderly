##orderly::orderly_develop_start()
alietal <- readRDS("alietal_data_clean.rds")
alietal$si <- as.numeric(alietal$si, units="days")
alietal$nu <- alietal$onset_first_iso
alietal <- alietal[ , c("si", "nu")]
alietal <- alietal[complete.cases(alietal), ]
