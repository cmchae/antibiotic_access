#### 0. prep ####
rm(list = ls())
gc()
par(mfrow = c(1,1))

requiredPackages = c("dplyr", "tidyr", "data.table", "stringr", "MASS")
lib.path <- .libPaths()
for(p in requiredPackages){
  if(!require(p,character.only = TRUE)) install.packages(p, lib = lib.path)
  library(p,character.only = TRUE)
}

n_sim <- 10000
set.seed(123)

#---------------------------prepare data for figure------------------------
UNI_indexed <- readRDS("outputs/UNI_indexed.rds")
coef_df_LRI <- readRDS("outputs/coef_df_LRI.rds")
coef_df_Fever <- readRDS("outputs/coef_df_Fever.rds")
coef_df_LRI_DDD <- readRDS("outputs/coef_df_LRI_DDD.rds")
coef_df_Fever_DDD <- readRDS("outputs/coef_df_Fever_DDD.rds")

health_exp <- readRDS("outputs/health_expenditure.rds")
UNI_indexed <- UNI_indexed %>% inner_join(health_exp, by = c("Country", "Year"))

UNI_indexed_Abaccess <- UNI_indexed %>% filter(Abaccess == 1)
UNI_indexed_Abdata <- UNI_indexed %>% filter(ab.data == 1)
UNI_indexed_both <- UNI_indexed %>% filter(Abaccess == 1, ab.data == 1)

expected_coef_names <- c("(Intercept)", "PC1", "PC2")

stopifnot(
  identical(names(coef_df_LRI), expected_coef_names),
  identical(names(coef_df_Fever), expected_coef_names),
  identical(names(coef_df_LRI_DDD), expected_coef_names),
  identical(names(coef_df_Fever_DDD), expected_coef_names)
)

X_pred <- model.matrix(
  ~ PC1 + PC2,
  data = UNI_indexed
)

# Match design matrix columns to each coefficient matrix
X_pred_LRI <- X_pred[, names(coef_df_LRI), drop = FALSE]
X_pred_Fever <- X_pred[, names(coef_df_Fever), drop = FALSE]

X_pred_LRI_DDD <- X_pred[, names(coef_df_LRI_DDD), drop = FALSE]
X_pred_Fever_DDD <- X_pred[, names(coef_df_Fever_DDD), drop = FALSE]


# Linear predictors
LRI_matrix <- as.matrix(X_pred_LRI) %*%
  t(as.matrix(coef_df_LRI))

Fever_matrix <- as.matrix(X_pred_Fever) %*%
  t(as.matrix(coef_df_Fever))

LRI_DDD_matrix <- as.matrix(X_pred_LRI_DDD) %*%
  t(as.matrix(coef_df_LRI_DDD))

Fever_DDD_matrix <- as.matrix(X_pred_Fever_DDD) %*%
  t(as.matrix(coef_df_Fever_DDD))


# Antibiotic-treatment models: logistic link
AbLRI_sim <- plogis(LRI_matrix)
AbFever_sim <- plogis(Fever_matrix)


# DDD models: log-linear models
LRI_DDD_sim <- exp(LRI_DDD_matrix)
Fever_DDD_sim <- exp(Fever_DDD_matrix)

# Flatten the matrices into vectors (column-wise)
AbLRI_long <- as.vector(AbLRI_sim)
LRI_DDD_long <- as.vector(LRI_DDD_sim)
AbFever_long <- as.vector(AbFever_sim)
Fever_DDD_long <- as.vector(Fever_DDD_sim)

# Create matching indices
Year_index <- rep(UNI_indexed$Year, times = ncol(AbLRI_sim))
Country_index <- rep(UNI_indexed$Country, times = ncol(AbLRI_sim))
GDP_index <- rep(UNI_indexed$GDP, times = ncol(AbLRI_sim))
Malaria_index <- rep(UNI_indexed$Endemic, times = ncol(AbLRI_sim))
Health_exp_index <- rep(UNI_indexed$Current_health_expenditure_per_capita_PPP, times = ncol(AbLRI_sim))

# Combine into dataframe
sim_df_LRI <- data.frame(
  Country = Country_index,
  Year = Year_index,
  GDP = GDP_index,
  Malaria = Malaria_index,
  AbLRI = AbLRI_long,
  LRI_DDD = LRI_DDD_long,
  Health_exp = Health_exp_index
)

sim_df_Fever <- data.frame(
  Country = Country_index,
  Year = Year_index,
  GDP = GDP_index,
  Malaria = Malaria_index,
  AbFever = AbFever_long,
  Fever_DDD = Fever_DDD_long,
  Health_exp = Health_exp_index
)

sim_df_LRI_summary <- sim_df_LRI %>%
  group_by(Country, Year) %>%
  summarise(
    GDP = median(GDP, na.rm = TRUE),
    AbLRI_median = median(AbLRI, na.rm = TRUE),
    LRI_DDD_median = median(LRI_DDD, na.rm = TRUE),
    .groups = "drop")

sim_df_Fever_summary <- sim_df_Fever %>%
  group_by(Country, Year) %>%
  summarise(
    GDP = median(GDP, na.rm = TRUE),
    AbFever_median = median(AbFever, na.rm = TRUE),
    Fever_DDD_median = median(Fever_DDD, na.rm = TRUE),
    .groups = "drop")

sampled_df_ordered_LRI <- sim_df_LRI %>% arrange(LRI_DDD) %>%
  mutate(LRI_DDD_rank = cut(LRI_DDD,breaks = quantile(LRI_DDD, probs = seq(0, 1, 0.1), na.rm = TRUE),include.lowest = TRUE))

sampled_df_ordered_Fever <- sim_df_Fever %>% arrange(Fever_DDD) %>%
  mutate(Fever_DDD_rank = cut(Fever_DDD,breaks = quantile(Fever_DDD, probs = seq(0, 1, 0.1), na.rm = TRUE),include.lowest = TRUE))

group_df_LRI <- sampled_df_ordered_LRI %>%
  group_by(LRI_DDD_rank) %>%
  summarise(
    med_AbLRI   = median(AbLRI, na.rm = TRUE),
    lower_AbLRI = quantile(AbLRI, 0.025, na.rm = TRUE),
    upper_AbLRI = quantile(AbLRI, 0.975, na.rm = TRUE),
    seq_DDD     = median(LRI_DDD, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  na.omit()

group_df_Fever <- sampled_df_ordered_Fever %>%
  group_by(Fever_DDD_rank) %>%
  summarise(
    med_AbFever   = median(AbFever, na.rm = TRUE),
    lower_AbFever = quantile(AbFever, 0.025, na.rm = TRUE),
    upper_AbFever = quantile(AbFever, 0.975, na.rm = TRUE),
    seq_DDD       = median(Fever_DDD, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  na.omit()



  # ------------------------1. P(Antibiotic treatment|infection) ~ Antibiotic consumption/incidence-----------------
  sim_df_LRI_every_1000 <- sim_df_LRI[seq(1, nrow(sim_df_LRI), by = 10000), ]
  year_range <- range(sim_df_LRI_every_1000$Year)
  sim_df_Fever_every_1000 <- sim_df_Fever[seq(1, nrow(sim_df_Fever), by = 10000), ]

  layout(matrix(c(1,2), nrow = 1, byrow = TRUE))
  xmin <- min(group_df_LRI$seq_DDD); xmax <- max(group_df_LRI$seq_DDD)

  ymin <- 0 ; ymax <- 1
  
  par(oma = c(0, 0, 0, 0))
  par(mar = c(2.5, 2.5, 2, 1))  # Bottom, Left, Top, Right
  par(mgp = c(3, 0.3, 0))     # Axis label and tick distance
  par(tck = -0.025)
  par(lwd = 1)
  
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
  
  axis(1, at = seq(xmin, xmax, by = (xmax - xmin)/4), labels = FALSE, las = 2, pos = ymin)
  text(x = seq(xmin, xmax, by = (xmax - xmin)/4), y = ymin-0.08*(ymax-ymin),  
       labels = format(round(seq(xmin, xmax, by = ((xmax - xmin)/4)), digits = 0), nsmall = 0), xpd = TRUE, cex = 0.65)
  
  axis(2,  at = seq(ymin, ymax, by = (ymax - ymin)/4), labels = FALSE, las = 2, pos = xmin)
  text(y = seq(ymin, ymax, by = (ymax-ymin)/4), x = xmin-0.04*(xmax-xmin),  
       labels = format(round(seq(ymin, ymax, by = (ymax-ymin)/4), digits = 2), nsmall = 1), srt = 0, adj = 1, xpd = TRUE, cex = 0.65)

  lines(med_AbLRI ~ seq_DDD, data = group_df_LRI, col = "#0072B2", lwd = 1.5)
  
  # Add shaded area for confidence interval
  polygon(c(group_df_LRI$seq_DDD, rev(group_df_LRI$seq_DDD)),
          c(group_df_LRI$lower_AbLRI, rev(group_df_LRI$upper_AbLRI)),
          col = rgb(0/255, 114/255, 178/255, 0.25), border = NA)
  
  # Probability of antibotics treatment by antibiotic consumption (DDD) per incidence
  mtext("A. Lower Respiratory Infection", at = xmin - 0.2*(xmax-xmin), side = 3, line = 1, adj = 0, cex = 0.65, font = 2)
  mtext("Antibiotic consumption (DDD) per incidence", side = 1, line = 1.5, cex = 0.65)
  mtext(expression("Probability of antibiotic treatment"), side = 2, line = 1.5, cex = 0.65)  
  
  xmin <-  min(group_df_Fever$seq_DDD); xmax <- max(group_df_Fever$seq_DDD)
  
  plot(med_AbFever ~ seq_DDD, 
       data = group_df_Fever,
       xlab = "", 
       ylab = "",
       pch = 16, col = rgb(0.75, 0.75, 0.75, 0),
       xlim = c(xmin,xmax),
       ylim = c(ymin,ymax),
       lwd = "1",
       axes=F,ann=F,
       main = "")
  
  axis(1, at = seq(xmin, xmax, by = (xmax - xmin)/4), labels = FALSE, las = 2, pos = ymin)
  text(x = seq(xmin, xmax, by = (xmax - xmin)/4), y = ymin-0.08*(ymax-ymin),  
       labels = format(round(seq(xmin, xmax, by = ((xmax - xmin)/4)), digits = 1), nsmall = 0), xpd = TRUE, cex = 0.65)
  
  axis(2, at = seq(ymin, ymax, by = (ymax - ymin)/4), labels = FALSE, las = 2, pos = xmin)
  text(y = seq(ymin, ymax, by = (ymax-ymin)/4), x = xmin-0.04*(xmax-xmin),  
       labels = format(round(seq(ymin, ymax, by = (ymax-ymin)/4), digits = 2), nsmall = 1), srt = 0, adj = 1, xpd = TRUE, cex = 0.65)
  
  lines(med_AbFever ~ seq_DDD, data = group_df_Fever, col = "#0072B2", lwd = 1.5)
  
  # Add shaded area for confidence interval
  polygon(c(group_df_Fever$seq_DDD, rev(group_df_Fever$seq_DDD)),
          c(group_df_Fever$lower_AbFever, rev(group_df_Fever$upper_AbFever)),
          col = rgb(0/255, 114/255, 178/255, 0.25), border = NA)

  mtext("B. Febrile Illness without Lower Respiratory Infection", at = xmin - 0.2*(xmax-xmin), side = 3, line = 1, adj = 0, cex = 0.65, font = 2)
  mtext("Antibiotic consumption (DDD) per incidence", side = 1, line = 1.5, cex = 0.65)

  # ------------------------2.Antibiotic consumption/incidence ~ GDP --------------------
  LRI_DDD_sim_long <- as.vector(LRI_DDD_sim)
  LRI_DDD_gdp <- as.data.frame(UNI_indexed$GDP)
  
  for (i in 1:nrow(LRI_DDD_gdp)){
    LRI_DDD_gdp[i,2] <- quantile(LRI_DDD_sim[i,1:n_sim], probs = 0.5)
    LRI_DDD_gdp[i,3] <- quantile(LRI_DDD_sim[i,1:n_sim], probs = 0.025)
    LRI_DDD_gdp[i,4] <- quantile(LRI_DDD_sim[i,1:n_sim], probs = 0.975)
  }
  
  colnames(LRI_DDD_gdp) = c("GDP","med","lower","upper")

  # Create and process the merged dataset
  LRI_DDD_gdp_long <- rep(LRI_DDD_gdp$GDP, n_sim)
  LRI_DDD_gdp_merged <- data.frame(LRI_DDD = LRI_DDD_sim_long, GDP = LRI_DDD_gdp_long) %>%
    arrange(GDP) %>% mutate(GDP_rank = cut(GDP,breaks = quantile(GDP, probs = seq(0, 1, 0.1), na.rm = TRUE),include.lowest = TRUE))
  
  LRI_DDD_gdp_plot <- LRI_DDD_gdp_merged %>% group_by(GDP_rank) %>%
    summarize(
      med_DDD = median(LRI_DDD),
      lower_DDD = quantile(LRI_DDD, 0.025),
      upper_DDD = quantile(LRI_DDD, 0.975),
      GDP = mean(GDP),
      .groups = "drop")
  
  
  Fever_DDD_gdp <- as.data.frame(UNI_indexed$GDP)
  for (i in 1:nrow(Fever_DDD_gdp)){
    Fever_DDD_gdp[i,2] <- quantile(Fever_DDD_sim[i,1:n_sim], probs = 0.5)
    Fever_DDD_gdp[i,3] <- quantile(Fever_DDD_sim[i,1:n_sim], probs = 0.025)
    Fever_DDD_gdp[i,4] <- quantile(Fever_DDD_sim[i,1:n_sim], probs = 0.975)
  }
  
  colnames(Fever_DDD_gdp) = c("GDP","med","lower","upper")

  Fever_DDD_gdp <- Fever_DDD_gdp %>% arrange(Fever_DDD_gdp$GDP)
  
  Fever_DDD_sim_long <- as.vector(Fever_DDD_sim)
  Fever_DDD_gdp_long <- rep(Fever_DDD_gdp$GDP, n_sim)
  
  # Create and process the merged dataset
  Fever_DDD_gdp_merged <- data.frame(Fever_DDD = Fever_DDD_sim_long, GDP = Fever_DDD_gdp_long) %>%
    arrange(GDP) %>% mutate(GDP_rank = cut(GDP,breaks = quantile(GDP, probs = seq(0, 1, 0.1), na.rm = TRUE),include.lowest = TRUE))
  
  Fever_DDD_gdp_plot <- Fever_DDD_gdp_merged %>% group_by(GDP_rank) %>%
    summarize(
      med_DDD = median(Fever_DDD),
      lower_DDD = quantile(Fever_DDD, 0.025),
      upper_DDD = quantile(Fever_DDD, 0.975),
      GDP = mean(GDP),
      .groups = "drop")
  
  
  

  #-----------------------P(Antibiotic treatment|infection) ~ GDP ----------------------
  # 3.1 AbLRI ~ GDP -------------------------
  AbLRI_gdp <- as.data.frame(UNI_indexed$GDP)
  AbLRI_health_exp <- as.data.frame(UNI_indexed$Current_health_expenditure_per_capita_PPP)
  AbLRI_sim_long <- as.vector(AbLRI_sim)
  AbLRI_gdp_long <- as.numeric(rep(AbLRI_gdp$'UNI_indexed$GDP', n_sim))
  AbLRI_health_exp_long <- as.numeric(rep(AbLRI_health_exp$'UNI_indexed$Current_health_expenditure_per_capita_PPP', n_sim))
    
  # Create and process the merged dataset
  AbLRI_gdp_merged <- data.frame(AbLRI = AbLRI_sim_long, GDP = AbLRI_gdp_long) %>% 
    arrange(GDP) %>%
    mutate(GDP_rank = cut(GDP,breaks = quantile(GDP, probs = seq(0, 1, 0.1), na.rm = TRUE),include.lowest = TRUE))
  
  AbLRI_health_exp_merged <- data.frame(AbLRI = AbLRI_sim_long, health_exp = AbLRI_health_exp_long) %>% 
    arrange(health_exp) %>%
    mutate(health_exp_rank = cut(health_exp, breaks = quantile(health_exp, probs = seq(0, 1, 0.1), na.rm = TRUE),include.lowest = TRUE)) 
  
  AbLRI_gdp_plot <- AbLRI_gdp_merged %>% group_by(GDP_rank) %>% 
    summarize(
      med_AbLRI = median(AbLRI),
      lower_AbLRI = quantile(AbLRI, 0.025),
      upper_AbLRI = quantile(AbLRI, 0.975),
      GDP = mean(GDP),
      .groups = "drop")
  
  AbLRI_health_exp_plot <- AbLRI_health_exp_merged %>% group_by(health_exp_rank) %>% 
    summarize(
      med_AbLRI = median(AbLRI),
      lower_AbLRI = quantile(AbLRI, 0.025),
      upper_AbLRI = quantile(AbLRI, 0.975),
      health_exp = mean(health_exp),
      .groups = "drop")
  
  AbFever_gdp <- as.data.frame(UNI_indexed$GDP)
  AbFever_health_exp <- as.data.frame(UNI_indexed$Current_health_expenditure_per_capita_PPP)

  AbFever_sim_long <- as.vector(AbFever_sim)
  AbFever_gdp_long <- as.numeric(rep(AbFever_gdp$'UNI_indexed$GDP', n_sim))
  AbFever_health_exp_long <- as.numeric(rep(AbFever_health_exp$'UNI_indexed$Current_health_expenditure_per_capita_PPP', n_sim))
  
  # Create and process the merged dataset
  AbFever_gdp_merged <- data.frame(AbFever = AbFever_sim_long, GDP = AbFever_gdp_long) %>% 
    arrange(GDP) %>%
    mutate(GDP_rank = cut(GDP,breaks = quantile(GDP, probs = seq(0, 1, 0.1), na.rm = TRUE),include.lowest = TRUE))

  AbFever_health_exp_merged <- data.frame(AbFever = AbFever_sim_long, health_exp = AbFever_health_exp_long) %>% 
    arrange(health_exp) %>%
    mutate(health_exp_rank = cut(health_exp, breaks = quantile(health_exp, probs = seq(0, 1, 0.1), na.rm = TRUE),include.lowest = TRUE)) 
  
  AbFever_gdp_plot <- AbFever_gdp_merged %>% group_by(GDP_rank) %>% 
    summarize(
      med_AbFever = median(AbFever),
      lower_AbFever = quantile(AbFever, 0.025),
      upper_AbFever = quantile(AbFever, 0.975),
      GDP = mean(GDP),
      .groups = "drop")
  
  AbFever_health_exp_plot <- AbFever_health_exp_merged %>% group_by(health_exp_rank) %>% 
    summarize(
      med_AbFever = median(AbFever),
      lower_AbFever = quantile(AbFever, 0.025),
      upper_AbFever = quantile(AbFever, 0.975),
      health_exp = mean(health_exp),
      .groups = "drop")
  
  dir.create(
    "outputs/global_map_df/economic",
    recursive = TRUE,
    showWarnings = FALSE
  )
  
  saveRDS(AbLRI_gdp_plot, "outputs/global_map_df/economic/LRI_GDP.rds")
  saveRDS(AbFever_gdp_plot, "outputs/global_map_df/economic/Fever_GDP.rds")
  saveRDS(AbLRI_health_exp_plot, "outputs/global_map_df/economic/LRI_health_exp.rds")
  saveRDS(AbFever_health_exp_plot, "outputs/global_map_df/economic/Fever_health_exp.rds")
  
  