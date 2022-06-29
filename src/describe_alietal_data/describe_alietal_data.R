# This task completes basic descriptive analysis of the Ali et al. data
# to produces figures 1-6 for the supplementary material

# read in the data
tp_data <- readRDS("alietal_data_clean.rds")
linelist <- readRDS("alietal_linelist_clean.rds")
contacts <- readRDS("alietal_contacts_clean.rds")

tp_data <- tp_data %>% filter(onset_first_iso != "NA")

# Figure 1A - the serial interval distribution for all pairs

si_median <- median(tp_data$si)
si_lowquant <- quantile(tp_data$si, probs = 0.25)
si_hiquant <- quantile(tp_data$si, probs = 0.75)
si_mean <- mean(tp_data$si)

pS1A <- ggplot(tp_data)+
  geom_histogram(aes(x = as.numeric(si)), binwidth = 1,
                 col = "white", fill = "grey")+
  xlab(" ")+
  theme_minimal()+
  theme(text = element_text(size = 18),
        plot.title = element_text(hjust = 0.5, size = 15, face = "italic"),
        axis.text.x = element_blank())+
  geom_vline(xintercept = si_median,
             lwd = 1)+
  geom_vline(xintercept = si_lowquant,
             lty = 2,
             lwd = 1)+
  geom_vline(xintercept = si_hiquant,
             lty = 2,
             lwd = 1)+
  geom_vline(xintercept = si_mean,
             lwd = 1, col = "firebrick")+
  ggtitle("all pairs (including those in clusters)")

# Figure 1B - the serial interval distribution for discrete pairs

tp_data_d <- tp_data %>% filter(cluster_size == 2)
si_median <- median(tp_data_d$si)
si_lowquant <- quantile(tp_data_d$si, probs = 0.25)
si_hiquant <- quantile(tp_data_d$si, probs = 0.75)
si_mean <- mean(tp_data_d$si)

pS1B <- ggplot(tp_data_d)+
  geom_histogram(aes(x = as.numeric(si)), binwidth = 1,
                 col = "white", fill = "grey")+
  xlab("Serial Interval (days)")+
  theme_minimal()+
  theme(text = element_text(size = 18),
        plot.title = element_text(hjust = 0.5, size = 15, face = "italic"))+
  geom_vline(xintercept = si_median,
             lwd = 1)+
  geom_vline(xintercept = si_lowquant,
             lty = 2,
             lwd = 1)+
  geom_vline(xintercept = si_hiquant,
             lty = 2,
             lwd = 1)+
  geom_vline(xintercept = si_mean,
             lwd = 1, col = "firebrick")+
  ggtitle("discrete pairs only")


# Figure S1C - distribution of cluster size

cluster_df <- as.data.frame(table(tp_data$cluster_size), stringsAsFactors = FALSE)
cluster_df <- cluster_df %>% 
  mutate(frequency = as.numeric(Freq)/(as.numeric(Var1) - 1))
cluster_dummy <- data.frame(Var1 = 11:21, Freq = 0, frequency = 0)
cluster_df <- rbind(cluster_df, cluster_dummy)
cluster_df$Var1 <- factor(cluster_df$Var1, levels = 1:22)

pS1C <- ggplot(cluster_df)+
  geom_bar(aes(x = Var1, y = frequency), stat="identity")+
  theme_minimal()+
  xlab("cluster size")+
  ylab("frequency")+
  theme(text = element_text(size = 18),
        axis.text = element_text(size = 12))

# Figure S1D - cluster size and SI variability

# find group mean SI for each cluster size
gm <- tp_data %>% 
  group_by(cluster_size) %>% 
  summarise(si = mean(si))

pS1D <- ggplot(data = tp_data, aes(x = cluster_size, y = si))+
  geom_point(alpha = 0.3, size = 3, shape = 21, fill = "black", colour = "transparent")+
  geom_point(data = gm, shape = 17, size = 3, colour = "firebrick")+
  theme_minimal()+
  xlab("cluster size")+
  ylab("Serial Interval (days)")+
  theme(text = element_text(size = 18))

pS1 <- plot_grid(pS1A, pS1B, pS1C, pS1D, ncol = 2, nrow = 2, labels = c("A", "C", "B", "D"), 
          align = "h", byrow = F)

ggsave(filename = "pS1.png", pS1)


# Figure s2 - SI for household versus community transmission pairs

household.labs <- c("community", "household")
names(household.labs) <- c("no", "yes")

tp_data_median <- tp_data %>%
  group_by(isHousehold) %>%
  summarise(median.si = median(si))

tp_data_lowquant <- tp_data %>%
  group_by(isHousehold) %>%
  summarise(lowquant.si = quantile(si, probs = 0.25))

tp_data_hiquant <- tp_data %>%
  group_by(isHousehold) %>%
  summarise(hiquant.si = quantile(si, probs = 0.75))

tp_data_mean <- tp_data %>%
  group_by(isHousehold) %>%
  summarise(mean.si = mean(si))

pS2 <- ggplot(data = tp_data)+
  geom_histogram(aes(x = as.numeric(si), 
                     group = isHousehold), 
                 binwidth = 1,
                 col = "white", fill = "grey")+
  facet_wrap(~isHousehold, 
             labeller = labeller(isHousehold = household.labs),
             nrow = 2)+
  theme_minimal()+
  xlab("Serial Interval (days)")+
  geom_vline(data = tp_data_median,
             mapping = aes(xintercept = median.si),
             lwd = 1)+
  geom_vline(data = tp_data_lowquant,
             mapping = aes(xintercept = lowquant.si),
             lty = 2,
             lwd = 1)+
  geom_vline(data = tp_data_hiquant,
             mapping = aes(xintercept = hiquant.si), 
             lty = 2,
             lwd = 1)+
  geom_vline(data = tp_data_mean,
             mapping = aes(xintercept = mean.si),
             lwd = 1, col = "firebrick")+
  theme(text = element_text(size = 20),
        strip.text = element_text(face = "italic",
                                  size = 18))

ggsave(filename = "pS2.png", pS2)


# Figure S3 - delay from symptom onset to isolation (primary case)

pS3A <- ggplot(data = tp_data,
               aes(y = as.numeric(onset_first_iso), x = as.numeric(si)))+
  geom_point(size = 3, alpha = 0.5, shape = 21, fill = "black", colour = "transparent")+
  geom_smooth(method = "lm", colour = "firebrick")+
  theme_minimal()+
  xlab("Serial Interval (days)")+
  ylab("Delay from symptom onset to \nisolation of the primary case (days)")+
  theme(text = element_text(size = 20))

cor.test(method = "pearson", 
         x = as.numeric(tp_data$onset_first_iso), 
         y = as.numeric(tp_data$si), 
         use = "pairwise.complete.obs")


pS3B <- ggplot(tp_data)+
  geom_histogram(aes(x = as.numeric(onset_first_iso)), binwidth = 1,
                 col = "white", fill = "grey")+
  xlab("Delay from symptom onset to\nisolation of the primary case (days)")+
  theme_minimal()+
  theme(text = element_text(size = 18))

pS3 <- plot_grid(rel_widths = c(1, 1), 
             pS3A, pS3B, 
             ncol = 2, 
             labels = c("A", "B"))

ggsave(filename = "pS3.png", pS3)

# Figure S4 correlation between SI and delay to importation

tp_data <- tp_data%>%
  mutate(onset_imp = infector_returnDate_fromOtherCity - infector_onsetDate)

pS4A <- ggplot(data = tp_data,
               aes(y = as.numeric(onset_imp), x = as.numeric(si)))+
  geom_point(size = 3, alpha = 0.5, shape = 21, fill = "black", colour = "transparent")+
  geom_smooth(method = "lm", colour = "firebrick")+
  theme_minimal()+
  xlab("Serial Interval (days)")+
  ylab("Delay from symptom onset to\narrival in the area (days)")+
  theme(text = element_text(size = 20))


cor.test(method = "pearson", 
         x = as.numeric(tp_data$onset_imp), 
         y = as.numeric(tp_data$si), 
         use = "pairwise.complete.obs")
# comparing observed SI distributions

for(i in 1:length(tp_data$infector_returned_fromOtherCity)){
  
  if(tp_data$infector_returned_fromOtherCity[i] != "local" & 
     !is.na(tp_data$infector_returned_fromOtherCity[i])){
    tp_data$infector_returned_fromOtherCity[i] <- "imported"
  } }
tp_data$infector_returned_fromOtherCity[is.na(
  tp_data$infector_returned_fromOtherCity)] <- "unknown"
table_imp <- as.data.frame(table(tp_data$infector_returned_fromOtherCity))

sum(tp_data$onset_imp>0, na.rm=TRUE)
pS4B <- ggplot(table_imp)+
  geom_bar(aes(x = Var1, y = Freq), stat = "identity")+
  theme_minimal()+
  xlab("primary case")+
  ylab("frequency")+
  theme(text = element_text(size = 18))
  
pS4 <- plot_grid(rel_widths = c(2, 1), 
                 pS4A, pS4B, 
                 ncol = 2, 
                 labels = c("A", "B"))

ggsave(filename = "pS4.png", pS4)

# Figure 3 - potential extra investigations into importation time

# SI v nu but colour points by importation status of primary case
ggplot(data = tp_data,
       aes(y = as.numeric(onset_first_iso), x = as.numeric(si)))+
  geom_point(aes(fill = infector_returned_fromOtherCity), size = 3, alpha = 0.5, shape = 21, colour = "transparent")+
  geom_smooth(method = "lm", colour = "firebrick")+
  theme_minimal()+
  xlab("Serial Interval (days)")+
  ylab("Delay from symptom onset to \nisolation of the primary case (days)")+
  theme(text = element_text(size = 20))

# difference in the difference between isolation delay and the si
# by importation status of the primary case
ggplot(data = tp_data)+
  geom_boxplot(aes(x = infector_returned_fromOtherCity,
                   y = onset_first_iso - si,
                   fill = infector_returned_fromOtherCity))+
  theme_minimal()

# difference in isolation delay by importation status of primary case
ggplot(data = tp_data)+
  geom_boxplot(aes(x = infector_returned_fromOtherCity,
                   y = onset_first_iso,
                   fill = infector_returned_fromOtherCity))+
  theme_minimal()

# correlation between delay to isolation and delay to importation
ggplot(data = tp_data,
       aes(y = as.numeric(onset_first_iso), x = as.numeric(onset_imp)))+
  geom_point(aes(fill = infector_returned_fromOtherCity), size = 3, alpha = 0.5, shape = 21, colour = "transparent")+
  geom_smooth(method = "lm", colour = "firebrick")+
  theme_minimal()+
  xlab("Delay from symptom onset to return (days)")+
  ylab("Delay from symptom onset to \nisolation of the primary case (days)")+
  theme(text = element_text(size = 20))+
  geom_abline(slope = 1, intercept = 0, col = "firebrick", lty = 2)

cor.test(method = "pearson", 
         x = as.numeric(tp_data$onset_imp), 
         y = as.numeric(tp_data$onset_first_iso), 
         use = "pairwise.complete.obs")

# correlation between diff between si and delay to isolation and delay to importation
ggplot(data = tp_data,
       aes(y = as.numeric(onset_first_iso - si), x = as.numeric(onset_imp)))+
  geom_point(aes(fill = infector_returned_fromOtherCity), size = 3, alpha = 0.5, shape = 21, colour = "transparent")+
  geom_smooth(method = "lm", colour = "firebrick")+
  theme_minimal()+
  xlab("Delay from symptom onset to return (days)")+
  ylab("Difference between delay to isolation\nand symptom onset of secondary")+
  theme(text = element_text(size = 20))
cor.test(method = "pearson", 
         x = as.numeric(tp_data$onset_imp), 
         y = as.numeric(tp_data$onset_first_iso - tp_data$si), 
         use = "pairwise.complete.obs")

# correlation between date of symptom onset and delay to isolation
ggplot(data = tp_data,
       aes(y = as.numeric(onset_first_iso), x = infector_onsetDate))+
  geom_point(aes(fill = infector_returned_fromOtherCity), size = 3, alpha = 0.5, shape = 21, colour = "transparent")+
  geom_smooth(method = "lm", colour = "firebrick")+
  theme_minimal()+
  xlab("Date of symptom onset of primary case")+
  ylab("Delay to isolation of primary")+
  theme(text = element_text(size = 20))
cor.test(method = "pearson", 
         x = as.numeric(tp_data$infector_onsetDate), 
         y = as.numeric(tp_data$onset_first_iso), 
         use = "pairwise.complete.obs")

# correlation between date of symptom onset and delay to importation
ggplot(data = tp_data,
       aes(y = as.numeric(onset_imp), x = infector_onsetDate))+
  geom_point(aes(fill = infector_returned_fromOtherCity), size = 3, alpha = 0.5, shape = 21, colour = "transparent")+
  geom_smooth(method = "lm", colour = "firebrick")+
  theme_minimal()+
  xlab("Date of symptom onset of primary case")+
  ylab("Delay from symptom onset to importation")+
  theme(text = element_text(size = 20))
cor.test(method = "pearson", 
         x = as.numeric(tp_data$infector_onsetDate), 
         y = as.numeric(tp_data$onset_imp), 
         use = "pairwise.complete.obs")

# tp_data_median_imp <- tp_data %>%
#   group_by(infector_returned_fromOtherCity) %>%
#   summarise(median.si = median(si))
# 
# tp_data_lowquant_imp <- tp_data %>%
#   group_by(infector_returned_fromOtherCity) %>%
#   summarise(lowquant.si = quantile(si, probs = 0.25))
# 
# tp_data_hiquant_imp <- tp_data %>%
#   group_by(infector_returned_fromOtherCity) %>%
#   summarise(hiquant.si = quantile(si, probs = 0.75))
# 
# tp_data_mean_imp <- tp_data %>%
#   group_by(infector_returned_fromOtherCity) %>%
#   summarise(mean.si = mean(si))
# 
# ggplot(data = tp_data)+
#   geom_histogram(aes(x = as.numeric(si), 
#                      group = infector_returned_fromOtherCity), 
#                  binwidth = 1,
#                  col = "white", fill = "grey")+
#   facet_wrap(~infector_returned_fromOtherCity,
#              nrow = 3)+
#   theme_minimal()+
#   xlab("Serial Interval (days)")+
#   geom_vline(data = tp_data_median_imp,
#              mapping = aes(xintercept = median.si), lwd = 1)+
#   geom_vline(data = tp_data_lowquant_imp,
#              mapping = aes(xintercept = lowquant.si), lty = 2, lwd = 1)+
#   geom_vline(data = tp_data_hiquant_imp,
#              mapping = aes(xintercept = hiquant.si), lty = 2, lwd = 1)+
#   geom_vline(data = tp_data_mean_imp,
#              mapping = aes(xintercept = mean.si), col = "firebrick", lwd = 1)+
#   theme(text = element_text(size = 20))
# 
# 
# # Figure S4 - SI for discrete versus clustered pairs
# tp_data <- tp_data %>% mutate(discrete = cluster_size==2)
# 
# discrete.labs <- c("discrete pair", "part of cluster")
# names(discrete.labs) <- c("TRUE", "FALSE")
# 
# tp_data_median_disc <- tp_data %>%
#   filter(discrete == TRUE) %>%
#   group_by(discrete) %>%
#   summarise(median.si = median(si))
# 
# tp_data_lowquant_disc <- tp_data %>%
#   group_by(discrete) %>%
#   summarise(lowquant.si = quantile(si, probs = 0.25))
# 
# tp_data_hiquant_disc <- tp_data %>%
#   group_by(discrete) %>%
#   summarise(hiquant.si = quantile(si, probs = 0.75))
# 
# tp_data_mean_disc <- tp_data %>%
#   group_by(discrete) %>%
#   summarise(mean.si = mean(si))
# 
# ggplot(data = tp_data)+
#   geom_histogram(aes(x = as.numeric(si), 
#                      group = discrete), 
#                  binwidth = 1,
#                  col = "white", fill = "grey")+
#   facet_wrap(~discrete, 
#              labeller = labeller(discrete = discrete.labs),
#              nrow = 2)+
#   theme_minimal()+
#   xlab("Serial Interval (days)")+
#   geom_vline(data = tp_data_median_disc,
#              mapping = aes(xintercept = median.si),
#              lwd = 1)+
#   geom_vline(data = tp_data_lowquant_disc,
#              mapping = aes(xintercept = lowquant.si),
#              lty = 2,
#              lwd = 1)+
#   geom_vline(data = tp_data_hiquant_disc,
#              mapping = aes(xintercept = hiquant.si), 
#              lty = 2,
#              lwd = 1)+
#   geom_vline(data = tp_data_mean_disc,
#              mapping = aes(xintercept = mean.si),
#              lwd = 1, col = "firebrick")+
#   theme(text = element_text(size = 20))


