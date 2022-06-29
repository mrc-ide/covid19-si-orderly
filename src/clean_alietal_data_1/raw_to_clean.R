# orderly::orderly_develop_start()
# Here we take the data made available by Ali et al.
# https://github.com/PDGLin/COVID19_EffSerialInterval_NPI
# and explore it and subset it for our purposes, producing to two outputs
# 1. A data frame with one row per transmission pair that is used in fitting and data exploration
# 2. A linelist with one row per case used in data exploration

###################################
# TRANSMISSION PAIR DATA CLEANING #
###################################

data <- read.csv("alietal_data.csv", stringsAsFactors = FALSE) # csv downloaded from github

# change "unknown" --> NA
data[data == "unknown"] <- NA

# create new columns to indicate local/imported/NA
data <- data %>%
  mutate(infector_returned_fromOtherCity = infector_returnDate_fromOtherCity)%>%
  mutate(infectee_returned_fromOtherCity = infectee_returnDate_fromOtherCity)

# replace local/? with NA in date of return from other city columns
data$infectee_returnDate_fromOtherCity[data$infectee_returnDate_fromOtherCity == "local"] <- NA
data$infectee_returnDate_fromOtherCity[data$infectee_returnDate_fromOtherCity == "?"] <- NA
data$infector_returnDate_fromOtherCity[data$infector_returnDate_fromOtherCity == "local"] <- NA

# change dates to date format
date_index <- grep("Date|HospitalVisit", x = names(data))
data[ , date_index] <- lapply(data[ , date_index], dmy)

# subset to exclude pairs with unknown onset date
si_data <- data%>%filter(!is.na(infector_onsetDate))%>%
  filter(!is.na(infectee_onsetDate))%>%
  mutate(si = infectee_onsetDate - infector_onsetDate) # create si column

# how many index cases have 2 dates of isolation?
si_data <- si_data%>%
  mutate(before_symp = !is.na(infector_isolateDate_beforeSymptom))%>%
  mutate(after_symp = !is.na(infector_isolateDate_afterSymptom))%>%
  mutate(both_dates = before_symp + after_symp)

# what is the difference between the two dates?
si_data <- si_data %>%
  mutate(onset_iso_b4 = infector_isolateDate_beforeSymptom - infector_onsetDate)%>%
  mutate(onset_iso_aft = infector_isolateDate_afterSymptom - infector_onsetDate)%>%
  mutate(iso_diff = infector_isolateDate_afterSymptom - infector_isolateDate_beforeSymptom)

# create a column using the earliest date of isolation available
si_data <- si_data %>%
  mutate(onset_first_iso = pmap_dbl(list(onset_iso_b4, onset_iso_aft), min, na.rm = TRUE))
#when both values were NA it resulted in Inf
si_data$onset_first_iso[!is.finite(si_data$onset_first_iso)] <- NA

# what is the difference between date of first hospital visit and isolation after symptoms
si_data <- si_data%>%
  mutate(isohosp_diff = infector_isolateDate_afterSymptom - infector_firstHospitalVisit)


###################
# MAKE A LINELIST #
###################

# for epicontacts we need a linelist as well as a list of contacts

# separate the columns pertaining to infector and infectee so we can make long
infectee_ind <- grep("infectee",
                     x = names(si_data))

infector_ind <- grep("infector",
                     x = names(si_data))

infectee_ind <- infectee_ind[c(-1,-15)] # we don't want the relationship bc not symmetrical
infector_ind <- infector_ind[c(-14,-15)]

infectee_data <- si_data[,infectee_ind]
infector_data <- si_data[,infector_ind]

#get rid of the infectee or infector in column names so we can make long
ll_names <- c("id", "age", "sex", "reportCity",
              "infectCity", "onsetDate", "returnDate_fromOtherCity",
              "isolateDate_beforeSymptom", "firstHospitalVisit",
              "isolateDate_afterSymptom", "labConfirmDate", "disclosureDate",
              "severity", "returned_fromOtherCity")

names(infectee_data) <- ll_names
names(infector_data) <- ll_names

infectors <- infector_data$id # get list of ids of infectors
infectees <- infectee_data$id # get list of ids of infectees

# some people may be both infectees and infectors so through removal of duplicates
# we may lose info on who is an infector so retain here
linelist <- rbind(infectee_data, infector_data)
linelist_unique <- linelist[!duplicated(linelist$id), ]
linelist_unique$infector <- linelist_unique$id %in% infectors
linelist_unique$infectee <- linelist_unique$id %in% infectees

# extract contacts
contacts <- si_data%>%
  rename(from = infector_id, to = infectee_id)%>%
  mutate(si = infectee_onsetDate - infector_onsetDate)%>%
  select(from,to, si, isHousehold)

alietal_epicontacts <- make_epicontacts(linelist = linelist_unique,
                                        contacts = contacts,
                                        directed = TRUE)
# extract cluster information for each individual case
linelist_clusters <- get_clusters(alietal_epicontacts,
                                  output = "data.frame",
                                  member_col = "cluster_member",
                                  size_col = "cluster_size",
                                  override = FALSE)

# add cluster information to the transmission pair data frame
si_data <- si_data%>%
  left_join(x = si_data, y = linelist_clusters, by = c("infector_id" = "id"))

############################
# SAVE CLEAN/WRANGLED DATA #
############################
si_data$si <- as.numeric(si_data$si, units="days")
si_data$nu <- si_data$onset_first_iso
si_data <- si_data[! is.na(si_data$nu), ]
si_data <- arrange(si_data, nu)
saveRDS(file = "alietal_data_clean.rds", si_data) # save transmission pair data for fitting
saveRDS(file = "alietal_linelist_clean.rds", linelist_unique) # save linelist for data exploration
saveRDS(file = "alietal_contacts_clean.rds", contacts) # save contacts for data exploration

## S3S4 and Discrete pairs
data_discrete_pairs <- filter(si_data, cluster_size == 2)
saveRDS(data_discrete_pairs, "alietal_pairs.rds")

# whole data set with nu replacement for s3 type pairs
data_s3_s4mix <- si_data %>%
  mutate(
    infector_first_isolateDate = pmap_dbl(
      list(infector_isolateDate_beforeSymptom,
           infector_isolateDate_afterSymptom), min, na.rm = TRUE
    ),
    infectee_first_isolateDate = pmap_dbl(
      list(infectee_isolateDate_beforeSymptom,
           infectee_isolateDate_afterSymptom), min, na.rm = TRUE
    )
  )

## when a value is NA it results in Inf
data_s3_s4mix$infector_first_isolateDate[! is.finite(data_s3_s4mix$infector_first_isolateDate)] <- NA
data_s3_s4mix$infector_first_isolateDate <- as.Date(data_s3_s4mix$infector_first_isolateDate, origin = "1970-01-01")
## when both values are NA it results in Inf
data_s3_s4mix$infectee_first_isolateDate[! is.finite(data_s3_s4mix$infectee_first_isolateDate)] <- NA
data_s3_s4mix$infectee_first_isolateDate <- as.Date(data_s3_s4mix$infectee_first_isolateDate, origin = "1970-01-01")

for(i in 1:(length(data_s3_s4mix$infector_first_isolateDate))){
  if(!is.na(data_s3_s4mix$infectee_first_isolateDate[i])){
    if(data_s3_s4mix$infector_first_isolateDate[i] >= data_s3_s4mix$infectee_onsetDate[i] && data_s3_s4mix$infector_first_isolateDate[i] >= data_s3_s4mix$infectee_first_isolateDate[i]){
      if(data_s3_s4mix$infectee_first_isolateDate[i] >= data_s3_s4mix$infectee_onsetDate[i]){
        data_s3_s4mix$nu[i] <- 41
      }
    }
  }
}
data_s3_s4mix <- arrange(data_s3_s4mix, nu)
saveRDS(data_s3_s4mix, "alietal_part_isol.rds")
## discrete pairs with nu replacement for s3 type pairs
data_discrete_pairs_s3_s4mix <- filter(data_s3_s4mix, cluster_size == 2)
data_discrete_pairs_s3_s4mix <- arrange(data_discrete_pairs_s3_s4mix, nu)
saveRDS(data_s3_s4mix, "alietal_part_isol_pairs.rds")
