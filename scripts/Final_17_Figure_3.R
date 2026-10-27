#### 0. prep ####
rm(list = ls())
gc()

par(mfrow = c(1,1))

requiredPackages = c("dplyr", "tidyr", "data.table", "stringr")
lib.path <- .libPaths()
for(p in requiredPackages){
  if(!require(p,character.only = TRUE)) install.packages(p, lib = lib.path)
  library(p,character.only = TRUE)
}
n_sim <- 10000; set.seed(123)

#---------------------------prepare data for figure------------------------
UNI_indexed     <- readRDS("outputs/UNI_indexed.rds")
coef_df_LRI     <- readRDS("outputs/coef_df_LRI.rds")
coef_df_Fever   <- readRDS("outputs/coef_df_Fever.rds")
coef_df_LRI_DDD <- readRDS("outputs/coef_df_LRI_DDD.rds")
coef_df_Fever_DDD <- readRDS("outputs/coef_df_Fever_DDD.rds")

UNI_indexed_Abaccess <- UNI_indexed %>% filter(Abaccess == 1)
UNI_indexed_Abdata <- UNI_indexed %>% filter(ab.data == 1)
UNI_indexed_both <- UNI_indexed %>% filter(Abaccess == 1, ab.data == 1)

UNI_indexed_0 <- UNI_indexed[UNI_indexed$Endemic == 0,]
UNI_indexed_1 <- UNI_indexed[UNI_indexed$Endemic == 1,]

X_pred <- model.matrix(~ PC1 + PC2, data = UNI_indexed)
X_pred_0 <- model.matrix(~ PC1 + PC2, data = UNI_indexed_0)
X_pred_1 <- model.matrix(~ PC1 + PC2, data = UNI_indexed_1)

X_pred_LRI <- X_pred[, names(coef_df_LRI), drop = FALSE]
X_pred_Fever_0 <- X_pred_0[, names(coef_df_Fever), drop = FALSE]
X_pred_Fever_1 <- X_pred_1[, names(coef_df_Fever), drop = FALSE]

X_pred_LRI_DDD <- X_pred[, names(coef_df_LRI_DDD), drop = FALSE]
X_pred_Fever_DDD_0 <- X_pred_0[, names(coef_df_Fever_DDD), drop = FALSE]
X_pred_Fever_DDD_1 <- X_pred_1[, names(coef_df_Fever_DDD), drop = FALSE]

LRI_matrix <- as.matrix(X_pred_LRI) %*% t(as.matrix(coef_df_LRI))
Fever_matrix_0 <- as.matrix(X_pred_Fever_0) %*% t(as.matrix(coef_df_Fever))
Fever_matrix_1 <- as.matrix(X_pred_Fever_1) %*% t(as.matrix(coef_df_Fever))

LRI_DDD_matrix <- as.matrix(X_pred_LRI_DDD) %*% t(as.matrix(coef_df_LRI_DDD))
Fever_DDD_matrix_0 <- as.matrix(X_pred_Fever_DDD_0) %*% t(as.matrix(coef_df_Fever_DDD))
Fever_DDD_matrix_1 <- as.matrix(X_pred_Fever_DDD_1) %*% t(as.matrix(coef_df_Fever_DDD))

# convert to antibiotic treatment probability using logit link
AbLRI_sim <- plogis(LRI_matrix)
AbFever_sim_0 <- plogis(Fever_matrix_0)
AbFever_sim_1 <- plogis(Fever_matrix_1)

# convert to expected DDD using exp link function
LRI_DDD_sim <- exp(LRI_DDD_matrix)
Fever_DDD_sim_0 <- exp(Fever_DDD_matrix_0)
Fever_DDD_sim_1 <- exp(Fever_DDD_matrix_1)

# Flatten the matrices into vectors (column-wise)
AbLRI_long <- as.vector(AbLRI_sim)
LRI_DDD_long <- as.vector(LRI_DDD_sim)
AbFever_long_0 <- as.vector(AbFever_sim_0)
AbFever_long_1 <- as.vector(AbFever_sim_1)
Fever_DDD_long_0 <- as.vector(Fever_DDD_sim_0)
Fever_DDD_long_1 <- as.vector(Fever_DDD_sim_1)

# Create matching indices
Year_index <- rep(UNI_indexed$Year, times = ncol(AbLRI_sim))
Country_index <- rep(UNI_indexed$Country, times = ncol(AbLRI_sim))
GDP_index <- rep(UNI_indexed$GDP, times = ncol(AbLRI_sim))
Malaria_index <- rep(UNI_indexed$Endemic, times = ncol(AbLRI_sim))

Year_index_0 <- rep(UNI_indexed_0$Year, times = ncol(AbFever_sim_0))
Country_index_0 <- rep(UNI_indexed_0$Country, times = ncol(AbFever_sim_0))
GDP_index_0 <- rep(UNI_indexed_0$GDP, times = ncol(AbFever_sim_0))
Malaria_index_0 <- rep(UNI_indexed_0$Endemic, times = ncol(AbFever_sim_0))

Year_index_1 <- rep(UNI_indexed_1$Year, times = ncol(AbFever_sim_1))
Country_index_1 <- rep(UNI_indexed_1$Country, times = ncol(AbFever_sim_1))
GDP_index_1 <- rep(UNI_indexed_1$GDP, times = ncol(AbFever_sim_1))
Malaria_index_1 <- rep(UNI_indexed_1$Endemic, times = ncol(AbFever_sim_1))

# Combine into dataframe
sim_df_LRI <- data.frame(
  Country = Country_index,
  Year = Year_index,
  GDP = GDP_index,
  Malaria = Malaria_index,
  AbLRI = AbLRI_long,
  lack_AbLRI = 1-AbLRI_long,
  LRI_DDD = LRI_DDD_long
)

sim_df_Fever_0 <- data.frame(
  Country = Country_index_0,
  Year = Year_index_0,
  GDP = GDP_index_0,
  Malaria = Malaria_index_0,
  AbFever = AbFever_long_0,
  lack_AbFever = 1-AbFever_long_0,
  Fever_DDD = Fever_DDD_long_0
)

sim_df_Fever_1 <- data.frame(
  Country = Country_index_1,
  Year = Year_index_1,
  GDP = GDP_index_1,
  Malaria = Malaria_index_1,
  AbFever = AbFever_long_1,
  lack_AbFever = 1-AbFever_long_1,
  Fever_DDD = Fever_DDD_long_1
)

sim_df_LRI_summary <- sim_df_LRI %>%
  group_by(Country, Year) %>%
  summarise(
    GDP = median(GDP, na.rm = TRUE),
    AbLRI_median = median(AbLRI, na.rm = TRUE),
    AbLRI_median_lack = median(lack_AbLRI, na.rm = TRUE),
    LRI_DDD_median = median(LRI_DDD, na.rm = TRUE),
    .groups = "drop")

sim_df_Fever_summary_0 <- sim_df_Fever_0 %>%
  group_by(Country, Year) %>%
  summarise(
    GDP = median(GDP, na.rm = TRUE),
    AbFever_median = median(AbFever, na.rm = TRUE),
    AbFever_median_lack = median(lack_AbFever, na.rm = TRUE),
    Fever_DDD_median = median(Fever_DDD, na.rm = TRUE),
    .groups = "drop")

sim_df_Fever_summary_1 <- sim_df_Fever_1 %>%
  group_by(Country, Year) %>%
  summarise(
    GDP = median(GDP, na.rm = TRUE),
    AbFever_median = median(AbFever, na.rm = TRUE),
    AbFever_median_lack = median(lack_AbFever, na.rm = TRUE),
    Fever_DDD_median = median(Fever_DDD, na.rm = TRUE),
    .groups = "drop")

sampled_df_ordered_LRI <- sim_df_LRI %>% arrange(LRI_DDD) %>%
  mutate(LRI_DDD_rank = cut(LRI_DDD, breaks = quantile(LRI_DDD, probs = seq(0, 1, 0.1), na.rm = TRUE), include.lowest = TRUE))

sampled_df_ordered_Fever_0 <- sim_df_Fever_0 %>% arrange(Fever_DDD) %>%
  mutate(Fever_DDD_rank = cut(Fever_DDD, breaks = quantile(Fever_DDD, probs = seq(0, 1, 0.1), na.rm = TRUE), include.lowest = TRUE))

sampled_df_ordered_Fever_1 <- sim_df_Fever_1 %>% arrange(Fever_DDD) %>%
  mutate(Fever_DDD_rank = cut(Fever_DDD, breaks = quantile(Fever_DDD, probs = seq(0, 1, 0.1), na.rm = TRUE), include.lowest = TRUE))

group_df_LRI <- sampled_df_ordered_LRI %>% group_by(LRI_DDD_rank) %>% 
  summarize(
    med_AbLRI = median(AbLRI),
    lower_AbLRI = quantile(AbLRI, 0.025),
    upper_AbLRI = quantile(AbLRI, 0.975),
    lack_med_AbLRI = median(lack_AbLRI),
    lack_lower_AbLRI = quantile(lack_AbLRI, 0.025),
    lack_upper_AbLRI = quantile(lack_AbLRI, 0.975),
    seq_DDD = median(LRI_DDD),
    .groups = "drop") %>%  
  distinct() %>% na.omit()

group_df_Fever_0 <- sampled_df_ordered_Fever_0 %>% group_by(Fever_DDD_rank) %>% 
  summarize(
    med_AbFever = median(AbFever),
    lower_AbFever = quantile(AbFever, 0.025),
    upper_AbFever = quantile(AbFever, 0.975),
    lack_med_AbFever = median(lack_AbFever),
    lack_lower_AbFever = quantile(lack_AbFever, 0.025),
    lack_upper_AbFever = quantile(lack_AbFever, 0.975),
    seq_DDD = median(Fever_DDD),
    .groups = "drop") %>%  
  distinct() %>% na.omit()

group_df_Fever_1 <- sampled_df_ordered_Fever_1 %>% group_by(Fever_DDD_rank) %>% 
  summarize(
    med_AbFever = median(AbFever),
    lower_AbFever = quantile(AbFever, 0.025),
    upper_AbFever = quantile(AbFever, 0.975),
    lack_med_AbFever = median(lack_AbFever),
    lack_lower_AbFever = quantile(lack_AbFever, 0.025),
    lack_upper_AbFever = quantile(lack_AbFever, 0.975),
    seq_DDD = median(Fever_DDD),
    .groups = "drop") %>%  
  distinct() %>% na.omit()

group_df_LRI$ideal <- NA
for (i in 1:10){
  intercept <- group_df_LRI$med_AbLRI[1]
  slope <- (group_df_LRI$med_AbLRI[2] - group_df_LRI$med_AbLRI[1]) / (group_df_LRI$seq_DDD[2] - group_df_LRI$seq_DDD[1])
  group_df_LRI$ideal[i] <- intercept + slope * (group_df_LRI$seq_DDD[i] - group_df_LRI$seq_DDD[1])
}
group_df_LRI$ideal <- 1-group_df_LRI$ideal

group_df_Fever_0$ideal <- NA
for (i in 1:10){
  intercept <- group_df_Fever_0$med_AbFever[1]
  slope <- (group_df_Fever_0$med_AbFever[2] - group_df_Fever_0$med_AbFever[1]) / (group_df_Fever_0$seq_DDD[2] - group_df_Fever_0$seq_DDD[1])
  group_df_Fever_0$ideal[i] <- intercept + slope * (group_df_Fever_0$seq_DDD[i] - group_df_Fever_0$seq_DDD[1])
}
group_df_Fever_0$ideal <- 1-group_df_Fever_0$ideal

group_df_Fever_1$ideal <- NA
for (i in 1:10){
  intercept <- group_df_Fever_1$med_AbFever[1]
  slope <- (group_df_Fever_1$med_AbFever[2] - group_df_Fever_1$med_AbFever[1]) / (group_df_Fever_1$seq_DDD[2] - group_df_Fever_1$seq_DDD[1])
  group_df_Fever_1$ideal[i] <- intercept + slope * (group_df_Fever_1$seq_DDD[i] - group_df_Fever_1$seq_DDD[1])
}
group_df_Fever_1$ideal <- 1-group_df_Fever_1$ideal

LRI_gdp_plot   <- readRDS("outputs/global_map_df/economic/LRI_GDP.rds")
Fever_gdp_plot <- readRDS("outputs/global_map_df/economic/Fever_GDP.rds")

UNI_indexed <- readRDS("outputs/UNI_indexed.rds")
health_exp <- readRDS("outputs/health_expenditure.rds")
UNI_indexed <- UNI_indexed %>% inner_join(health_exp, by = c("Country", "Year"))

UNI_indexed_Abaccess <- UNI_indexed %>% filter(Abaccess == 1)
UNI_indexed_Abdata <- UNI_indexed %>% filter(ab.data == 1)
UNI_indexed_both <- UNI_indexed %>% filter(Abaccess == 1, ab.data == 1)

LRI_health_exp_plot   <- readRDS("outputs/global_map_df/economic/LRI_health_exp.rds")
Fever_health_exp_plot <- readRDS("outputs/global_map_df/economic/Fever_health_exp.rds")

intercept_LRI_ideal <- group_df_LRI$med_AbLRI[1]
slope_LRI_ideal <- (group_df_LRI$med_AbLRI[2] - group_df_LRI$med_AbLRI[1]) / (group_df_LRI$seq_DDD[2] - group_df_LRI$seq_DDD[1])
xmax_LRI <- (1 - intercept_LRI_ideal)/slope_LRI_ideal + group_df_LRI$seq_DDD[1]

intercept_LRI_obs <- group_df_LRI$med_AbLRI[9]
slope_LRI_obs <- (group_df_LRI$med_AbLRI[10] - group_df_LRI$med_AbLRI[9]) / (group_df_LRI$seq_DDD[10] - group_df_LRI$seq_DDD[9])
med_AbLRI_at_xmax <- intercept_LRI_obs + slope_LRI_obs * (xmax_LRI - group_df_LRI$seq_DDD[8])

group_df_LRI_2 <- group_df_LRI[1:8,]
group_df_LRI_2[9,] <- NA
group_df_LRI_2$seq_DDD[9] <- xmax_LRI
group_df_LRI_2$ideal[9] <- 0
group_df_LRI_2$lack_med_AbLRI[9] <- 1-med_AbLRI_at_xmax

intercept_Fever_0_ideal <- group_df_Fever_0$med_AbFever[1]
slope_Fever_0_ideal <- (group_df_Fever_0$med_AbFever[2] - group_df_Fever_0$med_AbFever[1]) / (group_df_Fever_0$seq_DDD[2] - group_df_Fever_0$seq_DDD[1])
xmax_Fever_0 <- (1 - intercept_Fever_0_ideal)/slope_Fever_0_ideal + group_df_Fever_0$seq_DDD[1]

intercept_Fever_0_obs <- group_df_Fever_0$med_AbFever[9]
slope_Fever_0_obs <- (group_df_Fever_0$med_AbFever[10] - group_df_Fever_0$med_AbFever[9]) / (group_df_Fever_0$seq_DDD[10] - group_df_Fever_0$seq_DDD[9])
med_AbFever_at_xmax_0 <- intercept_Fever_0_obs + slope_Fever_0_obs * (xmax_Fever_0 - group_df_Fever_0$seq_DDD[9])

group_df_Fever_0_2 <- group_df_Fever_0
group_df_Fever_1_2 <- group_df_Fever_1

write.csv(group_df_LRI, "outputs/plot_df/Figure3_LRI.csv")
write.csv(group_df_Fever_0, "outputs/plot_df/Figure3_Fever0.csv")
write.csv(group_df_Fever_1, "outputs/plot_df/Figure3_Fever1.csv")
write.csv(group_df_LRI_2, "outputs/plot_df/Figure3_LRI_2.csv")
write.csv(group_df_Fever_0_2, "outputs/plot_df/Figure3_Fever0_2.csv")
write.csv(group_df_Fever_1_2, "outputs/plot_df/Figure3_Fever1_2.csv")



# ============================================================
# Export source data for Figure 3
# ============================================================
dir.create("outputs/plot_objects", recursive = TRUE, showWarnings = FALSE)


# ------------------------------------------------------------
# A. LRTI: antibiotic access by GDP
# ------------------------------------------------------------

fig3_A <- LRI_gdp_plot %>%
  transmute(
    Panel = "A",
    Syndrome = "LRTI",
    Malaria_status = "Not applicable",
    X_variable = "GDP per capita (PPP)",
    X_value = GDP,
    Median = med_AbLRI,
    Lower = lower_AbLRI,
    Upper = upper_AbLRI,
    Projected = NA_real_
  )


# ------------------------------------------------------------
# B. AFI: antibiotic access by GDP
# ------------------------------------------------------------

fig3_B <- Fever_gdp_plot %>%
  transmute(
    Panel = "B",
    Syndrome = "AFI",
    Malaria_status = "All",
    X_variable = "GDP per capita (PPP)",
    X_value = GDP,
    Median = med_AbFever,
    Lower = lower_AbFever,
    Upper = upper_AbFever,
    Projected = NA_real_
  )


# ------------------------------------------------------------
# C. LRTI: antibiotic access by health expenditure
# ------------------------------------------------------------

fig3_C <- LRI_health_exp_plot %>%
  transmute(
    Panel = "C",
    Syndrome = "LRTI",
    Malaria_status = "Not applicable",
    X_variable = "Health expenditure per capita (PPP)",
    X_value = health_exp,
    Median = med_AbLRI,
    Lower = lower_AbLRI,
    Upper = upper_AbLRI,
    Projected = NA_real_
  )


# ------------------------------------------------------------
# D. AFI: antibiotic access by health expenditure
# ------------------------------------------------------------

fig3_D <- Fever_health_exp_plot %>%
  transmute(
    Panel = "D",
    Syndrome = "AFI",
    Malaria_status = "All",
    X_variable = "Health expenditure per capita (PPP)",
    X_value = health_exp,
    Median = med_AbFever,
    Lower = lower_AbFever,
    Upper = upper_AbFever,
    Projected = NA_real_
  )


# ------------------------------------------------------------
# E. LRTI: treatment probability vs DDD per episode
# ------------------------------------------------------------

fig3_E <- group_df_LRI %>%
  transmute(
    Panel = "E",
    Syndrome = "LRTI",
    Malaria_status = "Not applicable",
    X_variable = "DDD per episode",
    X_value = seq_DDD,
    Median = med_AbLRI,
    Lower = lower_AbLRI,
    Upper = upper_AbLRI,
    Projected = NA_real_
  )


# ------------------------------------------------------------
# F. AFI: malaria non-endemic
# ------------------------------------------------------------

fig3_F <- group_df_Fever_0 %>%
  transmute(
    Panel = "F",
    Syndrome = "AFI",
    Malaria_status = "Malaria non-endemic",
    X_variable = "DDD per episode",
    X_value = seq_DDD,
    Median = med_AbFever,
    Lower = lower_AbFever,
    Upper = upper_AbFever,
    Projected = NA_real_
  )


# ------------------------------------------------------------
# G. AFI: malaria endemic
# ------------------------------------------------------------

fig3_G <- group_df_Fever_1 %>%
  transmute(
    Panel = "G",
    Syndrome = "AFI",
    Malaria_status = "Malaria endemic",
    X_variable = "DDD per episode",
    X_value = seq_DDD,
    Median = med_AbFever,
    Lower = lower_AbFever,
    Upper = upper_AbFever,
    Projected = NA_real_
  )


# ------------------------------------------------------------
# H. LRTI: observed vs projected access
# ------------------------------------------------------------

fig3_H <- group_df_LRI_2 %>%
  transmute(
    Panel = "H",
    Syndrome = "LRTI",
    Malaria_status = "Not applicable",
    X_variable = "DDD per episode",
    X_value = seq_DDD,
    Median = 1 - lack_med_AbLRI,
    Lower = NA_real_,
    Upper = NA_real_,
    Projected = 1 - ideal
  )


# ------------------------------------------------------------
# I. AFI: malaria non-endemic, observed vs projected
# ------------------------------------------------------------

fig3_I <- group_df_Fever_0_2 %>%
  transmute(
    Panel = "I",
    Syndrome = "AFI",
    Malaria_status = "Malaria non-endemic",
    X_variable = "DDD per episode",
    X_value = seq_DDD,
    Median = 1 - lack_med_AbFever,
    Lower = NA_real_,
    Upper = NA_real_,
    Projected = 1 - ideal
  )


# ------------------------------------------------------------
# J. AFI: malaria endemic, observed vs projected
# ------------------------------------------------------------

fig3_J <- group_df_Fever_1_2 %>%
  transmute(
    Panel = "J",
    Syndrome = "AFI",
    Malaria_status = "Malaria endemic",
    X_variable = "DDD per episode",
    X_value = seq_DDD,
    Median = 1 - lack_med_AbFever,
    Lower = NA_real_,
    Upper = NA_real_,
    Projected = 1 - ideal
  )


# ------------------------------------------------------------
# Combine all panels and save
# ------------------------------------------------------------

Figure_3_source_data <- bind_rows(
  fig3_A,
  fig3_B,
  fig3_C,
  fig3_D,
  fig3_E,
  fig3_F,
  fig3_G,
  fig3_H,
  fig3_I,
  fig3_J
) %>%
  arrange(Panel, X_value)


write.csv(
  Figure_3_source_data,
  "outputs/plot_objects/Figure_3_source_data.csv",
  row.names = FALSE
)





# ------------------------1. P(Antibiotic treatment|infection) ~ Antibiotic consumption/incidence-----------------

pdf("figures/Figure_3.pdf", width = 7.5, height = 7)


layout(matrix(c(
                1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                0, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
                2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
                0, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9,10,10,10,10,10,10,10,10,10,10,10,11,11,11,11,11,11,11,11,11,11,11,
                3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3,
                0,12,12,12,12,12,12,12,12,12,12,12,13,13,13,13,13,13,13,13,13,13,13,14,14,14,14,14,14,14,14,14,14,14,15,15,15,15,15,15,15,15,15,15,15), 
              nrow = 6, byrow = TRUE), 
   widths = c(1.2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1), 
   heights = c(0.1,1,0.1,1,0.1,1))

#### 1. (title 1) ####
par(mar = c(0, 0, 0, 0))
plot.new()
text(x = -0.035, y = 0.5, labels = "Antibiotic access by GDP and health expenditure", 
     adj = 0, cex = 0.9, font = 2, xpd = TRUE)

#### 2. (title 2)  ####
par(mar = c(0, 0, 0, 0))
plot.new()
text(x = -0.035, y = 0.5, labels = "Antibiotic access and antibiotic consumption relative to incidence", 
     adj = 0, cex = 0.9, font = 2, xpd = TRUE)

#### 3. (title 3)  ####
par(mar = c(0, 0, 0, 0))
plot.new()
text(x = -0.035, y = 0.5, labels = "Observed and expected antibiotic access, and antibiotic consumption relative to incidence", 
     adj = 0, cex = 0.9, font = 2, xpd = TRUE)


#### 4 (figure A) ####

par(mar = c(3, 2.5, 1.5, 1))
par(mgp = c(3, 0.3, 0))     
par(tck = -0.025)
par(lwd = 1)

xmin <- log10(1000)         ;     xmax <- log10(30000)       
ymin <- 0 ; ymax <- 1

plot(med_AbLRI ~ GDP, 
     data = LRI_gdp_plot,
     xlab = "", 
     ylab = "",
     pch = 16, col = rgb(0.75, 0.75, 0.75, 0),
     xlim = c(xmin,xmax),
     ylim = c(ymin,ymax),
     lwd = "1",
     axes=F,ann=F,
     main = "")


xaxis <- log10(c(1000,10000,30000))
xticks <- log10(c(1000, 10000))
xlabels <- round(c(1000, 2000, 10000, 20000, 30000), digits = 0)
minor_xlabels <- as.vector(c(sapply(c(1000), function(v) v * 2:9),sapply(c(10000), function(v) v * 2:3)))

axis(1, at = xaxis, labels = FALSE, las = 2, pos = ymin, tck = 0)
axis(1, at = xticks, labels = FALSE, las = 2, pos = ymin)
axis(1, at = log10(minor_xlabels), labels = FALSE, las = 2, pos = ymin, tck = -0.015, lwd = 1)
text(x = log10(xlabels), y = ymin-(ymax-ymin)*0.065,  
     labels = xlabels/1000, srt = 0, adj = 0.5, xpd = TRUE, cex = 0.65)

axis(2,  at = seq(ymin, ymax, by = (ymax - ymin)/4), labels = FALSE, las = 2, pos = xmin)
text(y = seq(ymin, ymax, by = (ymax-ymin)/4), x = xmin-0.04*(xmax-xmin),  
     labels = format(round(seq(ymin, ymax, by = (ymax-ymin)/4), digits = 2), nsmall = 1), srt = 0, adj = 1, xpd = TRUE, cex = 0.65)

lines(med_AbLRI ~ log10(GDP), data = LRI_gdp_plot, col = "#0072B2", lwd = 1.5)

# Add shaded area for confidence interval
polygon(c(log10(LRI_gdp_plot$GDP), rev(log10(LRI_gdp_plot$GDP))),
        c(LRI_gdp_plot$lower_AbLRI, rev(LRI_gdp_plot$upper_AbLRI)),
        col = rgb(0/255, 114/255, 178/255, 0.25), border = NA)

mtext(expression("Antibiotic treatment probability"), side = 2, line = 2, cex = 0.50)  
mtext("GDP per capita (thousands, PPP)", side = 1, line = 1, cex = 0.50)
mtext("A. LRTI", at = xmin - 0.33*(xmax-xmin), side = 3, line = 0.5, adj = 0, cex = 0.50, font = 1)




#### 5 (figure B) ####

xmin <- log10(1000)         ;     xmax <- log10(30000)       
ymin <- 0 ; ymax <- 1

plot(med_AbFever ~ GDP, 
     data = Fever_gdp_plot,
     xlab = "", 
     ylab = "",
     pch = 16, col = rgb(0.75, 0.75, 0.75, 0),
     xlim = c(xmin,xmax),
     ylim = c(ymin,ymax),
     lwd = "1",
     axes=F,ann=F,
     main = "")

xaxis <- log10(c(1000,10000,30000))
xticks <- log10(c(1000, 10000))
xlabels <- round(c(1000, 2000, 10000, 20000, 30000), digits = 0)
minor_xlabels <- as.vector(c(sapply(c(1000), function(v) v * 2:9),sapply(c(10000), function(v) v * 2:3)))

axis(1, at = xaxis, labels = FALSE, las = 2, pos = ymin, tck = 0)
axis(1, at = xticks, labels = FALSE, las = 2, pos = ymin)
axis(1, at = log10(minor_xlabels), labels = FALSE, las = 2, pos = ymin, tck = -0.015, lwd = 1)
text(x = log10(xlabels), y = ymin-(ymax-ymin)*0.065,  
     labels = xlabels/1000, srt = 0, adj = 0.5, xpd = TRUE, cex = 0.65)


axis(2,  at = seq(ymin, ymax, by = (ymax - ymin)/4), labels = FALSE, las = 2, pos = xmin)
text(y = seq(ymin, ymax, by = (ymax-ymin)/4), x = xmin-0.04*(xmax-xmin),  
     labels = format(round(seq(ymin, ymax, by = (ymax-ymin)/4), digits = 2), nsmall = 1), srt = 0, adj = 1, xpd = TRUE, cex = 0.65)

lines(med_AbFever ~ log10(GDP), data = Fever_gdp_plot, col = "#0072B2", lwd = 1.5)

# Add shaded area for confidence interval
polygon(c(log10(Fever_gdp_plot$GDP), rev(log10(Fever_gdp_plot$GDP))),
        c(Fever_gdp_plot$lower_AbFever, rev(Fever_gdp_plot$upper_AbFever)),
        col = rgb(0/255, 114/255, 178/255, 0.25), border = NA)

mtext("GDP per capita (thousands, PPP)", side = 1, line = 1, cex = 0.50)
mtext("B. AFI", at = xmin - 0.33*(xmax-xmin), side = 3, line = 0.5, adj = 0, cex = 0.50, font = 1)


#### 6 (figure C) ####

par(mar = c(3, 2.5, 1.5, 1))  # Bottom, Left, Top, Right
par(mgp = c(3, 0.3, 0))     # Axis label and tick distance
par(tck = -0.025)
par(lwd = 1)

xmin <- log10(50)         ;     xmax <- log10(2000)       
ymin <- 0 ; ymax <- 1

plot(med_AbLRI ~ health_exp, 
     data = LRI_health_exp_plot,
     xlab = "", 
     ylab = "",
     pch = 16, col = rgb(0.75, 0.75, 0.75, 0),
     xlim = c(xmin,xmax),
     ylim = c(ymin,ymax),
     lwd = "1",
     axes=F,ann=F,
     main = "")

xaxis <- log10(c(50,100,1000, 2000))
xticks <- log10(c(100, 1000))
xlabels <- round(c(50, 100, 200, 1000, 2000), digits = 0)
minor_xlabels <- as.vector(c(sapply(c(10), function(v) v * 5:9),sapply(c(100), function(v) v * 2:9),sapply(c(1000), function(v) v * 2)))

axis(1, at = xaxis, labels = FALSE, las = 2, pos = ymin, tck = 0)
axis(1, at = xticks, labels = FALSE, las = 2, pos = ymin)
axis(1, at = log10(minor_xlabels), labels = FALSE, las = 2, pos = ymin, tck = -0.015, lwd = 1)
text(x = log10(xlabels), y = ymin-(ymax-ymin)*0.065,  
     labels = xlabels, srt = 0, adj = 0.5, xpd = TRUE, cex = 0.65)

axis(2,  at = seq(ymin, ymax, by = (ymax - ymin)/4), labels = FALSE, las = 2, pos = xmin)
text(y = seq(ymin, ymax, by = (ymax-ymin)/4), x = xmin-0.04*(xmax-xmin),  
     labels = format(round(seq(ymin, ymax, by = (ymax-ymin)/4), digits = 2), nsmall = 1), srt = 0, adj = 1, xpd = TRUE, cex = 0.65)

lines(med_AbLRI ~ log10(health_exp), data = LRI_health_exp_plot, col = "#0072B2", lwd = 1.5)

# Add shaded area for confidence interval
polygon(c(log10(LRI_health_exp_plot$health_exp), rev(log10(LRI_health_exp_plot$health_exp))),
        c(LRI_health_exp_plot$lower_AbLRI, rev(LRI_health_exp_plot$upper_AbLRI)),
        col = rgb(0/255, 114/255, 178/255, 0.25), border = NA)

mtext("Health expenditure per capita (PPP)", side = 1, line = 1, cex = 0.50)
mtext("C. LRTI", at = xmin - 0.33*(xmax-xmin), side = 3, line = 0.5, adj = 0, cex = 0.50, font = 1)

#### 7 (figure D) ####
xmin <- log10(50)         ;     xmax <- log10(2000)       
ymin <- 0 ; ymax <- 1

plot(med_AbFever ~ health_exp, 
     data = Fever_health_exp_plot,
     xlab = "", 
     ylab = "",
     pch = 16, col = rgb(0.75, 0.75, 0.75, 0),
     xlim = c(xmin,xmax),
     ylim = c(ymin,ymax),
     lwd = "1",
     axes=F,ann=F,
     main = "")

xaxis <- log10(c(50,100,1000, 2000))
xticks <- log10(c(100, 1000))
xlabels <- round(c(50, 100, 200, 1000, 2000), digits = 0)
minor_xlabels <- as.vector(c(sapply(c(10), function(v) v * 5:9),sapply(c(100), function(v) v * 2:9),sapply(c(1000), function(v) v * 2)))

axis(1, at = xaxis, labels = FALSE, las = 2, pos = ymin, tck = 0)
axis(1, at = xticks, labels = FALSE, las = 2, pos = ymin)
axis(1, at = log10(minor_xlabels), labels = FALSE, las = 2, pos = ymin, tck = -0.015, lwd = 1)
text(x = log10(xlabels), y = ymin-(ymax-ymin)*0.065,  
     labels = xlabels, srt = 0, adj = 0.5, xpd = TRUE, cex = 0.65)


axis(2,  at = seq(ymin, ymax, by = (ymax - ymin)/4), labels = FALSE, las = 2, pos = xmin)
text(y = seq(ymin, ymax, by = (ymax-ymin)/4), x = xmin-0.04*(xmax-xmin),  
     labels = format(round(seq(ymin, ymax, by = (ymax-ymin)/4), digits = 2), nsmall = 1), srt = 0, adj = 1, xpd = TRUE, cex = 0.65)

lines(med_AbFever ~ log10(health_exp), data = Fever_health_exp_plot, col = "#0072B2", lwd = 1.5)

# Add shaded area for confidence interval
polygon(c(log10(Fever_health_exp_plot$health_exp), rev(log10(Fever_health_exp_plot$health_exp))),
        c(Fever_health_exp_plot$lower_AbFever, rev(Fever_health_exp_plot$upper_AbFever)),
        col = rgb(0/255, 114/255, 178/255, 0.25), border = NA)

mtext("Health expenditure per capita (PPP)", side = 1, line = 1, cex = 0.50)
mtext("D. AFI", at = xmin - 0.33*(xmax-xmin), side = 3, line = 0.5, adj = 0, cex = 0.50, font = 1)


#### 8. (Fig E)  ####
xmin <- log10(10)         ;     xmax <- log10(300)       
ymin <- 0 ; ymax <- 1

par(mar = c(3, 2.5, 1.5, 1))
par(mgp = c(3, 0.3, 0))     
par(tck = -0.025)
par(lwd = 1)

# plot(type = "n", plot_data$GDP_Capita_PPP, plot_data$AbARI,
plot(med_AbLRI ~ seq_DDD, 
     data = group_df_LRI,
     xlab = "", 
     ylab = "",
     pch = 16, col = rgb(0.75, 0.75, 0.75, 0),
     xlim = c(xmin,xmax),
     ylim = c(ymin,ymax),
     lwd = "1",
     axes=F,ann=F,
     main = "")

# title(main = "Probability of treatment for LRI ~ Antibiotics consumption per incidence", cex.main = 1.5)

xaxis <- log10(c(10,300))
xticks <- log10(c(10, 100))
xlabels <- round(c(10, 20, 100, 200, 300), digits = 0)
minor_xlabels <- as.vector(c(sapply(c(10), function(v) v * 2:9), sapply(c(100), function(v) v * 2:3)))

axis(1, at = xaxis, labels = FALSE, las = 2, pos = ymin, tck = 0)
axis(1, at = xticks, labels = FALSE, las = 2, pos = ymin)
axis(1, at = log10(minor_xlabels), labels = FALSE, las = 2, pos = ymin, tck = -0.015, lwd = 1)
text(x = log10(xlabels), y = ymin-(ymax-ymin)*0.065,  
     labels = xlabels, srt = 0, adj = 0.5, xpd = TRUE, cex = 0.65)

axis(2,  at = seq(ymin, ymax, by = (ymax - ymin)/4), labels = FALSE, las = 2, pos = xmin)
text(y = seq(ymin, ymax, by = (ymax-ymin)/4), x = xmin-0.04*(xmax-xmin),  
     labels = format(round(seq(ymin, ymax, by = (ymax-ymin)/4), digits = 2), nsmall = 1), srt = 0, adj = 1, xpd = TRUE, cex = 0.65)



lines(med_AbLRI ~ log10(seq_DDD), data = group_df_LRI, col = "#0072B2", lwd = 1.5)

# Add shaded area for confidence interval
polygon(c(log10(group_df_LRI$seq_DDD), rev(log10(group_df_LRI$seq_DDD))),
        c(group_df_LRI$lower_AbLRI, rev(group_df_LRI$upper_AbLRI)),
        col = rgb(0/255, 114/255, 178/255, 0.25), border = NA)

mtext("E. LRTI", at = xmin - 0.3*(xmax-xmin), side = 3, line = 0.5, adj = 0, cex = 0.5, font = 1)
mtext("DDD per episode", side = 1, line = 1, cex = 0.5)
mtext(expression("Antibiotic treatment probability"), side = 2, line = 2, cex = 0.5)  



#### 9. (Fig F) ####
xmin <- log10(1)         ;     xmax <- log10(7)       
ymin <- 0 ; ymax <- 1

plot(med_AbFever ~ seq_DDD, 
     data = group_df_Fever_0,
     xlab = "", 
     ylab = "",
     pch = 16, col = rgb(0.75, 0.75, 0.75, 0),
     xlim = c(xmin,xmax),
     ylim = c(ymin,ymax),
     lwd = "1",
     axes=F,ann=F,
     main = "")

xaxis <- log10(c(1,7))
xticks <- log10(c(1))
xlabels <- round(c(1, 2, 7), digits = 0)
minor_xlabels <- as.vector(sapply(c(1), function(v) v * 2:7))

axis(1, at = xaxis, labels = FALSE, las = 2, pos = ymin, tck = 0)
axis(1, at = xticks, labels = FALSE, las = 2, pos = ymin)
axis(1, at = log10(minor_xlabels), labels = FALSE, las = 2, pos = ymin, tck = -0.015, lwd = 1)
text(x = log10(xlabels), y = ymin-(ymax-ymin)*0.065,  
     labels = xlabels, srt = 0, adj = 0.5, xpd = TRUE, cex = 0.65)

axis(2,  at = seq(ymin, ymax, by = (ymax - ymin)/4), labels = FALSE, las = 2, pos = xmin)
text(y = seq(ymin, ymax, by = (ymax-ymin)/4), x = xmin-0.04*(xmax-xmin),  
     labels = format(round(seq(ymin, ymax, by = (ymax-ymin)/4), digits = 2), nsmall = 1), srt = 0, adj = 1, xpd = TRUE, cex = 0.65)


lines(med_AbFever ~ log10(seq_DDD), data = group_df_Fever_0, col = "#0072B2", lwd = 1.5)

# Add shaded area for confidence interval
polygon(c(log10(group_df_Fever_0$seq_DDD), rev(log10(group_df_Fever_0$seq_DDD))),
        c(group_df_Fever_0$lower_AbFever, rev(group_df_Fever_0$upper_AbFever)),
        col = rgb(0/255, 114/255, 178/255, 0.25), border = NA)

mtext("F. AFI (malaria non-endemic countries)", at = xmin - 0.3*(xmax-xmin), side = 3, line = 0.5, adj = 0, cex = 0.5, font = 1)
mtext("DDD per episode", side = 1, line = 1, cex = 0.5)# mtext(expression("Probability of antibotics treatment"), side = 2, line = 2, cex = 0.75)  





#### 10. (Fig G) ####
xmin <- log10(0.9)         ;     xmax <- log10(4)       
ymin <- 0 ; ymax <- 1

plot(med_AbFever ~ seq_DDD, 
     data = group_df_Fever_1,
     xlab = "", 
     ylab = "",
     pch = 16, col = rgb(0.75, 0.75, 0.75, 0),
     xlim = c(xmin,xmax),
     ylim = c(ymin,ymax),
     lwd = "1",
     axes=F,ann=F,
     main = "")

xticks <- log10(c(1))
xlabels <- c(0.9, 1, 2, 4)
minor_xlabels <- as.vector(c(0.9,1,2,3,4))

axis(1, at = xticks, labels = FALSE, las = 2, pos = ymin)
axis(1, at = log10(minor_xlabels), labels = FALSE, las = 2, pos = ymin, tck = -0.015, lwd = 1)
text(x = log10(xlabels), y = ymin-(ymax-ymin)*0.065,  
     labels = xlabels, srt = 0, adj = 0.5, xpd = TRUE, cex = 0.65)

axis(2,  at = seq(ymin, ymax, by = (ymax - ymin)/4), labels = FALSE, las = 2, pos = xmin)
text(y = seq(ymin, ymax, by = (ymax-ymin)/4), x = xmin-0.04*(xmax-xmin),  
     labels = format(round(seq(ymin, ymax, by = (ymax-ymin)/4), digits = 2), nsmall = 1), srt = 0, adj = 1, xpd = TRUE, cex = 0.65)


lines(med_AbFever ~ log10(seq_DDD), data = group_df_Fever_1, col = "#0072B2", lwd = 1.5)

# Add shaded area for confidence interval
polygon(c(log10(group_df_Fever_1$seq_DDD), rev(log10(group_df_Fever_1$seq_DDD))),
        c(group_df_Fever_1$lower_AbFever, rev(group_df_Fever_1$upper_AbFever)),
        col = rgb(0/255, 114/255, 178/255, 0.25), border = NA)

mtext("G. AFI (malaria-endemic countries)", at = xmin - 0.3*(xmax-xmin), side = 3, line = 0.5, adj = 0, cex = 0.5, font = 1)
mtext("DDD per episode", side = 1, line = 1, cex = 0.5)

#### 11. (legend) ####
par(mar = c(0, 0, 0, 0))
plot.new()

#### 12. (Fig H)  ####
xmin <- log10(10)         ;     xmax <- log10(300)       
ymin <- 0 ; ymax <- 1

par(mar = c(3, 2.5, 1.5, 1)) 
par(mgp = c(3, 0.3, 0))     
par(tck = -0.025)
par(lwd = 1)

plot(lack_med_AbLRI ~ seq_DDD, 
     data = group_df_LRI_2,
     xlab = "", 
     ylab = "",
     pch = 16, col = rgb(0.75, 0.75, 0.75, 0),
     xlim = c(xmin,xmax),
     ylim = c(ymin,ymax),
     lwd = "1",
     axes=F,ann=F,
     main = "")

xaxis <- log10(c(10,200))
xticks <- log10(c(10, 100))
xlabels <- round(c(10, 20, 100, 200), digits = 0)
minor_xlabels <- as.vector(c(sapply(c(10), function(v) v * 2:9), sapply(c(100), function(v) v * 2)))

axis(1, at = xaxis, labels = FALSE, las = 2, pos = ymin, tck = 0)
axis(1, at = xticks, labels = FALSE, las = 2, pos = ymin)
axis(1, at = log10(minor_xlabels), labels = FALSE, las = 2, pos = ymin, tck = -0.015, lwd = 1)
text(x = log10(xlabels), y = ymin-(ymax-ymin)*0.065,  
     labels = xlabels, srt = 0, adj = 0.5, xpd = TRUE, cex = 0.65)

axis(2,  at = seq(ymin, ymax, by = (ymax - ymin)/4), labels = FALSE, las = 2, pos = xmin)
text(y = seq(ymin, ymax, by = (ymax-ymin)/4), x = xmin-0.04*(xmax-xmin),  
     labels = format(round(seq(ymin, ymax, by = (ymax-ymin)/4), digits = 2), nsmall = 1), srt = 0, adj = 1, xpd = TRUE, cex = 0.65)

lines((1-lack_med_AbLRI) ~ log10(seq_DDD), data = group_df_LRI_2, col = "#0072B2", lwd = 1.5)
lines((1-ideal) ~ log10(seq_DDD), data = group_df_LRI_2, col = "darkblue", lwd = 1.5)

polygon(c(log10(group_df_LRI_2$seq_DDD), rev(log10(group_df_LRI_2$seq_DDD))),
        c((1-group_df_LRI_2$ideal), rev((1-group_df_LRI_2$lack_med_AbLRI))),
        col = rgb(0/255, 0/255, 139/255, 0.25), border = NA)

# Probability of antibotics treatment by antibiotic consumption (DDD) per incidence
mtext("H. LRTI", at = xmin - 0.3*(xmax-xmin), side = 3, line = 0.5, adj = 0, cex = 0.5, font = 1)
mtext("DDD per episode", side = 1, line = 1, cex = 0.5)
mtext(expression("Antibiotic treatment probability"), side = 2, line = 2, cex = 0.5)  






#### 13. (Fig I) ####
xmin <- log10(1)         ;     xmax <- log10(7)       
ymin <- 0 ; ymax <- 1

# plot(type = "n", plot_data$GDP_Capita_PPP, plot_data$AbARI,
plot(med_AbFever ~ seq_DDD, 
     data = group_df_Fever_0_2,
     xlab = "", 
     ylab = "",
     pch = 16, col = rgb(0.75, 0.75, 0.75, 0),
     xlim = c(xmin,xmax),
     ylim = c(ymin,ymax),
     lwd = "1", 
     axes=F,ann=F,
     main = "")

# title(main = "Probability of treatment for Fever ~ Antibiotics consumption per incidence", cex.main = 1.5)

xaxis <- log10(c(1,7))
xticks <- log10(c(1))
xlabels <- round(c(1, 2, 7), digits = 0)
minor_xlabels <- as.vector(sapply(c(1), function(v) v * 2:7))

axis(1, at = xaxis, labels = FALSE, las = 2, pos = ymin, tck = 0)
axis(1, at = xticks, labels = FALSE, las = 2, pos = ymin)
axis(1, at = log10(minor_xlabels), labels = FALSE, las = 2, pos = ymin, tck = -0.015, lwd = 1)
text(x = log10(xlabels), y = ymin-(ymax-ymin)*0.065,  
     labels = xlabels, srt = 0, adj = 0.5, xpd = TRUE, cex = 0.65)

axis(2,  at = seq(ymin, ymax, by = (ymax - ymin)/4), labels = FALSE, las = 2, pos = xmin)
text(y = seq(ymin, ymax, by = (ymax-ymin)/4), x = xmin-0.04*(xmax-xmin),  
     labels = format(round(seq(ymin, ymax, by = (ymax-ymin)/4), digits = 2), nsmall = 1), srt = 0, adj = 1, xpd = TRUE, cex = 0.65)

lines((1-lack_med_AbFever) ~ log10(seq_DDD), data = group_df_Fever_0_2, col = "#0072B2", lwd = 1.5)
lines((1-ideal) ~ log10(seq_DDD), data = group_df_Fever_0_2, col = "darkblue", lwd = 1.5)

# Add shaded area for confidence interval
polygon(c(log10(group_df_Fever_0_2$seq_DDD), rev(log10(group_df_Fever_0_2$seq_DDD))),
        c((1-group_df_Fever_0_2$ideal), rev((1-group_df_Fever_0_2$lack_med_AbFever))),
        col = rgb(0/255, 0/255, 139/255, 0.25), border = NA)
 
mtext("I. AFI (malaria non-endemic countries)", at = xmin - 0.3*(xmax-xmin), side = 3, line = 0.5, adj = 0, cex = 0.5, font = 1)
mtext("DDD per episode", side = 1, line = 1, cex = 0.5)# mtext(expression("Probability of antibotics treatment"), side = 2, line = 2, cex = 0.75)  





#### 104 (Fig J) ####
xmin <- log10(0.9)         ;     xmax <- log10(4)       
ymin <- 0 ; ymax <- 1

# plot(type = "n", plot_data$GDP_Capita_PPP, plot_data$AbARI,
plot(med_AbFever ~ seq_DDD, 
     data = group_df_Fever_1,
     xlab = "", 
     ylab = "",
     pch = 16, col = rgb(0.75, 0.75, 0.75, 0),
     xlim = c(xmin,xmax),
     ylim = c(ymin,ymax),
     lwd = "1",
     axes=F,ann=F,
     main = "")

# title(main = "Probability of treatment for Fever ~ Antibiotics consumption per incidence", cex.main = 1.5)

xticks <- log10(c(1))
xlabels <- c(0.9, 1, 2, 4)
minor_xlabels <- as.vector(c(0.9,1,2,3,4))

axis(1, at = xticks, labels = FALSE, las = 2, pos = ymin)
axis(1, at = log10(minor_xlabels), labels = FALSE, las = 2, pos = ymin, tck = -0.015, lwd = 1)
text(x = log10(xlabels), y = ymin-(ymax-ymin)*0.065,  
     labels = xlabels, srt = 0, adj = 0.5, xpd = TRUE, cex = 0.65)

axis(2,  at = seq(ymin, ymax, by = (ymax - ymin)/4), labels = FALSE, las = 2, pos = xmin)
text(y = seq(ymin, ymax, by = (ymax-ymin)/4), x = xmin-0.04*(xmax-xmin),  
     labels = format(round(seq(ymin, ymax, by = (ymax-ymin)/4), digits = 2), nsmall = 1), srt = 0, adj = 1, xpd = TRUE, cex = 0.65)

lines((1-lack_med_AbFever) ~ log10(seq_DDD), data = group_df_Fever_1_2, col = "#0072B2", lwd = 1.5)
lines((1-ideal) ~ log10(seq_DDD), data = group_df_Fever_1_2, col = "darkblue", lwd = 1.5)

# Add shaded area for confidence interval
polygon(c(log10(group_df_Fever_1_2$seq_DDD), rev(log10(group_df_Fever_1_2$seq_DDD))),
        c((1-group_df_Fever_1_2$ideal), rev((1-group_df_Fever_1_2$lack_med_AbFever))),
        col = rgb(0/255, 0/255, 139/255, 0.25), border = NA)

mtext("J. AFI (malaria-endemic countries)", at = xmin - 0.3*(xmax-xmin), side = 3, line = 0.5, adj = 0, cex = 0.5, font = 1)
mtext("DDD per episode", side = 1, line = 1, cex = 0.5)
# mtext(expression("Probability of antibotics treatment"), side = 2, line = 2, cex = 0.75)  


#### 15. (legend) ####
par(mar = c(0, 0, 0, 0))
plot.new()

legend_labels <- c(
  "Projected value from lowest consumption deciles", 
  "Observed value", 
  "Gap between observed and expected values"
)
legend_colors <- c( "darkblue","#0072B2", "darkblue")
n <- length(legend_labels)


# Coordinates for the custom legend
x0 <- 0.0
x1 <- 0.07
y_start <- 0.61
line_spacing <- 0.07


for (i in seq_len(n)) {
  y <- y_start - (i - 1) * line_spacing
  
  if (i == 3) {
    # Polygon only (no line)
    polygon(
      x = c(x0-0.01, x1+0.01, x1+0.01, x0-0.01),
      y = c(y - 0.02, y - 0.02, y + 0.02, y + 0.02),
      col = rgb(0/255, 0/255, 139/255, 0.25),
      border = NA
    )
  } else {
    # Line for the first two
    segments(x0, y, x1, y, col = legend_colors[i], lwd = 1.5)
  }
  
  # Label
  text(x = x1 + 0.04, y = y, labels = legend_labels[i], adj = 0, cex = 0.56)
}

dev.off()



