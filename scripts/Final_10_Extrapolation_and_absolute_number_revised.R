#### 0. prep

rm(list = ls())
gc()

requiredPackages = c("dplyr", "tidyr", "data.table", "stringr", "MASS", "naniar")
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
coef_df_all_DDD <- readRDS("outputs/coef_df_all_DDD.rds")

# confirm all models use intercept + PC1 + PC2
stopifnot(
  identical(names(coef_df_LRI), c("(Intercept)", "PC1", "PC2")),
  identical(names(coef_df_Fever), c("(Intercept)", "PC1", "PC2")),
  identical(names(coef_df_LRI_DDD), c("(Intercept)", "PC1", "PC2")),
  identical(names(coef_df_Fever_DDD), c("(Intercept)", "PC1", "PC2")),
  identical(names(coef_df_all_DDD), c("(Intercept)", "PC1", "PC2"))
)

UNI_indexed_Abaccess <- UNI_indexed %>% filter(Abaccess == 1)
UNI_indexed_Abdata <- UNI_indexed %>% filter(ab.data == 1)
UNI_indexed_both <- UNI_indexed %>% filter(Abaccess == 1, ab.data == 1)

with(UNI_indexed, table(Country, Year))

# Antibiotic treatment: PC1-PC2
X_pred <- model.matrix(~ PC1 + PC2, data = UNI_indexed)
X_pred_LRI <- X_pred[, names(coef_df_LRI), drop = FALSE]
X_pred_Fever <- X_pred[, names(coef_df_Fever), drop = FALSE]

# DDD: PC1-PC2
X_pred_DDD <- model.matrix(~ PC1 + PC2, data = UNI_indexed)
X_pred_LRI_DDD <- X_pred_DDD[, names(coef_df_LRI_DDD), drop = FALSE]
X_pred_Fever_DDD <- X_pred_DDD[, names(coef_df_Fever_DDD), drop = FALSE]
X_pred_all_DDD <- X_pred_DDD[, names(coef_df_all_DDD), drop = FALSE]

LRI_matrix <- as.matrix(X_pred_LRI) %*% t(as.matrix(coef_df_LRI))
Fever_matrix <- as.matrix(X_pred_Fever) %*% t(as.matrix(coef_df_Fever))
LRI_DDD_matrix <- as.matrix(X_pred_LRI_DDD) %*% t(as.matrix(coef_df_LRI_DDD))
Fever_DDD_matrix <- as.matrix(X_pred_Fever_DDD) %*% t(as.matrix(coef_df_Fever_DDD))
all_DDD_matrix <- as.matrix(X_pred_all_DDD) %*% t(as.matrix(coef_df_all_DDD))

# Convert logit-scale predictions to antibiotic treatment probabilities
AbLRI_sim <- plogis(LRI_matrix)  # Dimensions: [country-years x 10000 simulations]
AbFever_sim <- plogis(Fever_matrix)  # Dimensions: [country-years x 10000 simulations]

# Convert log-scale predictions to expected DDD
LRI_DDD_sim <- exp(LRI_DDD_matrix)  # Dimensions: [country-years x 10000 simulations]
Fever_DDD_sim <- exp(Fever_DDD_matrix)  # Dimensions: [country-years x 10000 simulations]
all_DDD_sim <- exp(all_DDD_matrix)  # Dimensions: [country-years x 10000 simulations]

epsilon_min <- 1e-6
epsilon_max <- 1e-6 + 1e-9

LRI_DDD_sim[LRI_DDD_sim <= 0] <- runif(sum(LRI_DDD_sim <= 0), min = epsilon_min, max = epsilon_max)
Fever_DDD_sim[Fever_DDD_sim <= 0] <- runif(sum(Fever_DDD_sim <= 0), min = epsilon_min, max = epsilon_max)
all_DDD_sim[all_DDD_sim <= 0] <- runif(sum(all_DDD_sim <= 0), min = epsilon_min, max = epsilon_max)

# Flatten the matrices into vectors (column-wise)

AbLRI_long <- as.vector(AbLRI_sim)
LRI_DDD_long <- as.vector(LRI_DDD_sim)
AbFever_long <- as.vector(AbFever_sim)
Fever_DDD_long <- as.vector(Fever_DDD_sim)
all_DDD_long <- as.vector(all_DDD_sim)

# Create matching indices

Year_index <- rep(UNI_indexed$Year, times = ncol(AbLRI_sim))
Country_index <- rep(UNI_indexed$Country, times = ncol(AbLRI_sim))
GDP_index <- rep(UNI_indexed$GDP, times = ncol(AbLRI_sim))
Income_index <- rep(UNI_indexed$Income, times = ncol(AbLRI_sim))
Malaria_index <- rep(UNI_indexed$Endemic, times = ncol(AbLRI_sim))
GBD_index <- rep(UNI_indexed$gbd_region, times = ncol(AbLRI_sim))
Population <- rep(UNI_indexed$Population, times = ncol(AbLRI_sim))

# Combine into dataframe

sim_df_LRI <- data.frame(
  Country = Country_index,
  Year = Year_index,
  GDP = GDP_index,
  Income = Income_index,
  Malaria = Malaria_index,
  GBD_region = GBD_index,
  AbLRI = AbLRI_long,
  Population = Population,
  sim = rep(1:n_sim, each = nrow(AbLRI_sim))
)

# Combine into dataframe

sim_df_Fever <- data.frame(
  Country = Country_index,
  Year = Year_index,
  GDP = GDP_index,
  Income = Income_index,
  Malaria = Malaria_index,
  GBD_region = GBD_index,
  AbFever = AbFever_long,
  Population = Population,
  sim = rep(1:n_sim, each = nrow(AbFever_sim))
)

AbLRI_matrix <- sim_df_LRI %>% pivot_wider(names_from = sim, values_from = AbLRI)
AbFever_matrix <- sim_df_Fever %>% pivot_wider(names_from = sim, values_from = AbFever)

Country_Year_Antibiotics <- data.frame(
  Country = AbLRI_matrix$Country,
  Year = AbLRI_matrix$Year)

U5_matrix <- readRDS("outputs/matrix.result.ver.1.rds")

Country_Year_Incidence <- data.frame(
  Country = U5_matrix$Country,
  Year = U5_matrix$Year)

a <- unique(Country_Year_Antibiotics$Country)
b <- unique(Country_Year_Incidence$Country)
setdiff(a,b)
setdiff(b,a)

country_mapping <- c(
  "Bolivia" = "Bolivia (Plurinational State of)",
  "Congo, Dem. Rep." = "Democratic Republic of the Congo",
  "Congo, Rep." = "Congo",
  "Cote d'Ivoire" = "Côte d'Ivoire",
  "Egypt, Arab Rep." = "Egypt",
  "Gambia, The" = "Gambia",
  "Iran, Islamic Rep." = "Iran (Islamic Republic of)",
  "Kyrgyz Republic" = "Kyrgyzstan",
  "Lao PDR" = "Lao People's Democratic Republic",
  "Moldova" = "Republic of Moldova",
  "St. Vincent and the Grenadines" = "Saint Vincent and the Grenadines",
  "Tanzania" = "United Republic of Tanzania",
  "Turkiye" = "Türkiye",
  "Venezuela, RB" = "Venezuela (Bolivarian Republic of)",
  "West Bank and Gaza" = "Palestine",
  "Yemen, Rep." = "Yemen"  # No match in your `b`, might need manual check
)

Country_Year_Antibiotics <- Country_Year_Antibiotics %>% mutate(Country = ifelse(Country %in% names(country_mapping), country_mapping[Country], Country))
AbLRI_matrix <- AbLRI_matrix %>% mutate(Country = ifelse(Country %in% names(country_mapping), country_mapping[Country], Country))
AbFever_matrix <- AbFever_matrix %>% mutate(Country = ifelse(Country %in% names(country_mapping), country_mapping[Country], Country))

Country_Year_meta <- inner_join(Country_Year_Antibiotics, Country_Year_Incidence, by = c("Country", "Year"))

Incidence <- Country_Year_meta %>% inner_join(U5_matrix, by = c("Country","Year"))
Incidence_LRI_matrix <- Incidence[,3:(n_sim+2)] + Incidence[,(2*n_sim+3):(3*n_sim+2)]
Incidence_Fever_matrix <- Incidence[,(n_sim+3):(2*n_sim+2)]
Incidence_all_matrix <- Incidence[,3:(n_sim+2)] + Incidence[,(n_sim+3):(2*n_sim+2)] + Incidence[,(2*n_sim+3):(3*n_sim+2)]

Antibiotics_LRI <- Country_Year_meta %>% inner_join(AbLRI_matrix, by = c("Country","Year"))
Antibiotics_LRI_matrix <- Antibiotics_LRI[,8:(n_sim+7)]
Antibiotics_Fever <- Country_Year_meta %>% inner_join(AbFever_matrix, by = c("Country","Year"))
Antibiotics_Fever_matrix <- Antibiotics_Fever[,8:(n_sim+7)]

dim(Incidence_LRI_matrix)
dim(Antibiotics_LRI_matrix)
treated_LRI_matrix <- Incidence_LRI_matrix * Antibiotics_LRI_matrix
untreated_LRI_matrix <- Incidence_LRI_matrix - treated_LRI_matrix

dim(Incidence_Fever_matrix)
dim(Antibiotics_Fever_matrix)
treated_Fever_matrix <- Incidence_Fever_matrix * Antibiotics_Fever_matrix
untreated_Fever_matrix <- Incidence_Fever_matrix - treated_Fever_matrix


### for gini index
dir.create(
  "outputs/global_map_df/gini",
  recursive = TRUE,
  showWarnings = FALSE
)

saveRDS(Country_Year_meta,             "outputs/global_map_df/gini/Country_Year_meta.rds")
saveRDS(Antibiotics_LRI[,c(1,2,3,5)], "outputs/global_map_df/gini/Antibiotics_LRI.rds")
saveRDS(treated_LRI_matrix,            "outputs/global_map_df/gini/treated_LRI_matrix.rds")
saveRDS(Incidence_LRI_matrix,          "outputs/global_map_df/gini/Incidence_LRI_matrix.rds")
saveRDS(Antibiotics_Fever[,c(1,2,3,5)], "outputs/global_map_df/gini/Antibiotics_Fever.rds")
saveRDS(treated_Fever_matrix,          "outputs/global_map_df/gini/treated_Fever_matrix.rds")
saveRDS(Incidence_Fever_matrix,        "outputs/global_map_df/gini/Incidence_Fever_matrix.rds")

#-----------------------1.1 by country and year (treatment)------------------------
make_summary <- function(meta, income, treated, untreated, total) {
  
  # Compute proportion treated
  
  prop_treated <- treated / (treated + untreated)
  inc_rate_total <- (treated + untreated) / income$Population
  
  data.frame(
    Country = meta$Country,
    Year = meta$Year,
    Income = income$Income,
    Malaria = income$Malaria,
    GBD_region = income$GBD_region,
    Population = income$Population,
    
    Treated = sprintf("%s (%s  %s)",
                      formatC(apply(treated, 1, median, na.rm = TRUE), format = "f", big.mark = ",", digits = 0),
                      formatC(apply(treated, 1, quantile, probs = 0.025, na.rm = TRUE), format = "f", big.mark = ",", digits = 0),
                      formatC(apply(treated, 1, quantile, probs = 0.975, na.rm = TRUE), format = "f", big.mark = ",", digits = 0)),
    
    Untreated = sprintf("%s (%s  %s)",
                        formatC(apply(untreated, 1, median, na.rm = TRUE), format = "f", big.mark = ",", digits = 0),
                        formatC(apply(untreated, 1, quantile, probs = 0.025, na.rm = TRUE), format = "f", big.mark = ",", digits = 0),
                        formatC(apply(untreated, 1, quantile, probs = 0.975, na.rm = TRUE), format = "f", big.mark = ",", digits = 0)),
    
    Total = sprintf("%s (%s  %s)",
                    formatC(apply(total, 1, median, na.rm = TRUE), format = "f", big.mark = ",", digits = 0),
                    formatC(apply(total, 1, quantile, probs = 0.025, na.rm = TRUE), format = "f", big.mark = ",", digits = 0),
                    formatC(apply(total, 1, quantile, probs = 0.975, na.rm = TRUE), format = "f", big.mark = ",", digits = 0)),
    
    Incidence_rate = sprintf("%s (%s  %s)",
                             formatC(apply(inc_rate_total, 1, median, na.rm = TRUE), format = "f", big.mark = ",", digits = 3),
                             formatC(apply(inc_rate_total, 1, quantile, probs = 0.025, na.rm = TRUE), format = "f", big.mark = ",", digits = 3),
                             formatC(apply(inc_rate_total, 1, quantile, probs = 0.975, na.rm = TRUE), format = "f", big.mark = ",", digits = 3)),
    
    Treated_Proportion = sprintf("%.3f (%.3f, %.3f)",
                                 apply(prop_treated, 1, median, na.rm = TRUE),
                                 apply(prop_treated, 1, quantile, probs = 0.025, na.rm = TRUE),
                                 apply(prop_treated, 1, quantile, probs = 0.975, na.rm = TRUE))
  )
}

make_df <- function(meta, income, treated, untreated, total) {
  
  # Compute proportion treated
  
  prop_treated <- treated / (treated + untreated)
  inc_rate_total <- (treated + untreated) / income$Population
  
  data.frame(
    Country = meta$Country,
    Year = meta$Year,
    GDP = income$GDP,
    Malaria = income$Malaria,
    Income = income$Income,
    
    Treated_median   = apply(treated, 1, median, na.rm = TRUE),
    Treated_lower    = apply(treated, 1, quantile, probs = 0.025, na.rm = TRUE),
    Treated_upper    = apply(treated, 1, quantile, probs = 0.975, na.rm = TRUE),
    
    Untreated_median = apply(untreated, 1, median, na.rm = TRUE),
    Untreated_lower  = apply(untreated, 1, quantile, probs = 0.025, na.rm = TRUE),
    Untreated_upper  = apply(untreated, 1, quantile, probs = 0.975, na.rm = TRUE),
    
    Total_median     = apply(total, 1, median, na.rm = TRUE),
    Total_lower      = apply(total, 1, quantile, probs = 0.025, na.rm = TRUE),
    Total_upper      = apply(total, 1, quantile, probs = 0.975, na.rm = TRUE),
    
    Incidence_rate_median = apply(inc_rate_total, 1, median, na.rm = TRUE),
    Incidence_rate_lower  = apply(inc_rate_total, 1, quantile, probs = 0.025, na.rm = TRUE),
    Incidence_rate_upper  = apply(inc_rate_total, 1, quantile, probs = 0.975, na.rm = TRUE),
    
    Proportion_median = apply(prop_treated, 1, median, na.rm = TRUE),
    Proportion_lower  = apply(prop_treated, 1, quantile, probs = 0.025, na.rm = TRUE),
    Proportion_upper  = apply(prop_treated, 1, quantile, probs = 0.975, na.rm = TRUE)
  )
}

LRI_summary <- make_summary(Country_Year_meta, Antibiotics_LRI, treated_LRI_matrix, untreated_LRI_matrix, Incidence_LRI_matrix)
Fever_summary <- make_summary(Country_Year_meta, Antibiotics_Fever, treated_Fever_matrix, untreated_Fever_matrix, Incidence_Fever_matrix)

LRI_df <- make_df(Country_Year_meta, Antibiotics_LRI, treated_LRI_matrix, untreated_LRI_matrix, Incidence_LRI_matrix)
Fever_df <- make_df(Country_Year_meta, Antibiotics_Fever, treated_Fever_matrix, untreated_Fever_matrix, Incidence_Fever_matrix)

Fever_df_0 <- Fever_df[Fever_df$Malaria == 0,]
Fever_df_1 <- Fever_df[Fever_df$Malaria == 1,]

years_to_save <- 2003:2019

#### 1.1 CSV file

dir.create(
  "outputs/table/by_country_year/probability",
  recursive = TRUE,
  showWarnings = FALSE
)

for (yr in years_to_save) {
  write.csv(LRI_summary[LRI_summary$Year == yr, ], paste0("outputs/table/by_country_year/probability/LRI_", yr, ".csv"), row.names = FALSE)
  write.csv(Fever_summary[Fever_summary$Year == yr, ], paste0("outputs/table/by_country_year/probability/Fever_", yr, ".csv"), row.names = FALSE)
}

for (yr in years_to_save) {
  write.csv(LRI_df[LRI_df$Year == yr, ], paste0("outputs/global_map_df/1.1_Untreated_number_", yr, "_LRI.csv"), row.names = FALSE)
  write.csv(Fever_df_0[Fever_df_0$Year == yr, ], paste0("outputs/global_map_df/1.1_Untreated_number_0_", yr, "_Fever.csv"), row.names = FALSE)
  write.csv(Fever_df_1[Fever_df_1$Year == yr, ], paste0("outputs/global_map_df/1.1_Untreated_number_1_", yr, "_Fever.csv"), row.names = FALSE)
}








#----------------------------1.2 by year (treatment) ------------------------------

years <- 2003:2019
summary_list_LRI <- summary_list_LRI <- list()


#---LRI---
for (yr in years) {
  # Filter matrices by year
  treated_LRI <- cbind(Country_Year_meta, treated_LRI_matrix) %>%
    filter(Year == yr) %>%
    dplyr::select(-Country, -Year)
  
  untreated_LRI <- cbind(Country_Year_meta, untreated_LRI_matrix) %>%
    filter(Year == yr) %>%
    dplyr::select(-Country, -Year)
  
  total_LRI <- cbind(Country_Year_meta, Incidence_LRI_matrix) %>%
    filter(Year == yr) %>%
    dplyr::select(-Country, -Year)
  
  # Sum over countries for each simulation
  sum_treated_LRI   <- colSums(treated_LRI)
  sum_untreated_LRI <- colSums(untreated_LRI)
  sum_total_LRI     <- colSums(total_LRI)
  # --- population summed over countries (fixed, no sim dimension) ---
  sum_population_LRI <- sum(
    Antibiotics_LRI$Population[Country_Year_meta$Year == yr],
    na.rm = TRUE)
  incidence_LRI <- sum_total_LRI / sum_population_LRI
  # --- proportion treated ---
  proportion_treated_LRI <- sum_treated_LRI / 
    (sum_treated_LRI + sum_untreated_LRI)
  
  # Store results as long-format table
  summary_list_LRI[[as.character(yr)]] <- data.frame(
    Year = rep(yr, 15),
    Condition = rep(c("Treated LRI", "Untreated LRI", "Total LRI", "Incidence LRI", "Proportion LRI"), each = 3),
    Quantile = rep(c("Median", "2.5%", "97.5%"), times = 5),
    Value = c(
      quantile(sum_treated_LRI, probs = 0.50, na.rm = TRUE),
      quantile(sum_treated_LRI, probs = 0.025, na.rm = TRUE),
      quantile(sum_treated_LRI, probs = 0.975, na.rm = TRUE),
      
      quantile(sum_untreated_LRI, probs = 0.50, na.rm = TRUE),
      quantile(sum_untreated_LRI, probs = 0.025, na.rm = TRUE),
      quantile(sum_untreated_LRI, probs = 0.975, na.rm = TRUE),
      
      quantile(sum_total_LRI, probs = 0.50, na.rm = TRUE),
      quantile(sum_total_LRI, probs = 0.025, na.rm = TRUE),
      quantile(sum_total_LRI, probs = 0.975, na.rm = TRUE),
      
      quantile(incidence_LRI, probs = 0.50, na.rm = TRUE),
      quantile(incidence_LRI, probs = 0.025, na.rm = TRUE),
      quantile(incidence_LRI, probs = 0.975, na.rm = TRUE),
      
      quantile(proportion_treated_LRI, probs = 0.50, na.rm = TRUE),
      quantile(proportion_treated_LRI, probs = 0.025, na.rm = TRUE),
      quantile(proportion_treated_LRI, probs = 0.975, na.rm = TRUE)
    )
  )
}

# Combine all sublists
summary_all_years_LRI <- do.call(rbind, summary_list_LRI)

# Combine all years into one data frame
summary_wide_LRI <- summary_all_years_LRI %>%
  mutate(Column = paste0(gsub(" ", "_", Condition), "_", Quantile)) %>%
  dplyr::select(Year, Column, Value) %>%
  pivot_wider(names_from = Column, values_from = Value)

summary_wide_LRI_combined <- summary_wide_LRI %>%
  mutate(
    Treated_LRI = sprintf("%s (%s  %s)", 
                          formatC(Treated_LRI_Median, format = "f", digits = 0, big.mark = ","),
                          formatC(`Treated_LRI_2.5%`, format = "f", digits = 0, big.mark = ","),
                          formatC(`Treated_LRI_97.5%`, format = "f", digits = 0, big.mark = ",")),
    
    Untreated_LRI = sprintf("%s (%s  %s)", 
                            formatC(Untreated_LRI_Median, format = "f", digits = 0, big.mark = ","),
                            formatC(`Untreated_LRI_2.5%`, format = "f", digits = 0, big.mark = ","),
                            formatC(`Untreated_LRI_97.5%`, format = "f", digits = 0, big.mark = ",")),
    
    Total_LRI = sprintf("%s (%s  %s)", 
                        formatC(Total_LRI_Median, format = "f", digits = 0, big.mark = ","),
                        formatC(`Total_LRI_2.5%`, format = "f", digits = 0, big.mark = ","),
                        formatC(`Total_LRI_97.5%`, format = "f", digits = 0, big.mark = ",")),
    
    Incidence_LRI = sprintf("%.2f (%.2f  %.2f)",
      Incidence_LRI_Median,
      `Incidence_LRI_2.5%`,
      `Incidence_LRI_97.5%`
    ),
    
    Proportion_LRI = sprintf("%.3f (%.3f, %.3f)", 
                             Proportion_LRI_Median, 
                             `Proportion_LRI_2.5%`, 
                             `Proportion_LRI_97.5%`)
  ) %>%
  dplyr::select(Year, Treated_LRI, Untreated_LRI, Total_LRI, Incidence_LRI, Proportion_LRI)

LRI_year_proportion <- summary_wide_LRI %>%
  dplyr::select(Year, Proportion_LRI_Median, `Proportion_LRI_2.5%`, `Proportion_LRI_97.5%`) %>%
  dplyr::rename(Proportion_Median = Proportion_LRI_Median, Proportion_Lower = `Proportion_LRI_2.5%`, Proportion_Upper = `Proportion_LRI_97.5%`)




#---Fever---
years <- 2003:2019
summary_list_Fever <- list()

for (yr in years) {
  # Filter matrices by year
  treated_Fever <- cbind(Country_Year_meta, treated_Fever_matrix) %>%
    filter(Year == yr) %>%
    dplyr::select(-Country, -Year)
  
  untreated_Fever <- cbind(Country_Year_meta, untreated_Fever_matrix) %>%
    filter(Year == yr) %>%
    dplyr::select(-Country, -Year)
  
  total_Fever <- cbind(Country_Year_meta, Incidence_Fever_matrix) %>%
    filter(Year == yr) %>%
    dplyr::select(-Country, -Year)
  
  # Sum over countries for each simulation
  sum_treated_Fever <- colSums(treated_Fever)
  sum_untreated_Fever <- colSums(untreated_Fever)
  sum_total_Fever <- colSums(total_Fever)
  
  sum_population_Fever <- sum(
    Antibiotics_Fever$Population[Country_Year_meta$Year == yr],
    na.rm = TRUE)
  incidence_Fever <- sum_total_Fever / sum_population_Fever
  
  proportion_treated_Fever <- sum_treated_Fever / (sum_treated_Fever + sum_untreated_Fever)
  
  # Store results as long-format table
  summary_list_Fever[[as.character(yr)]] <- data.frame(
    Year = rep(yr, 15),
    Condition = rep(c("Treated Fever", "Untreated Fever", "Total Fever", "Incidence Fever","Proportion Fever"), each = 3),
    Quantile = rep(c("Median", "2.5%", "97.5%"), times = 5),
    Value = c(
      quantile(sum_treated_Fever, probs = 0.50, na.rm = TRUE),
      quantile(sum_treated_Fever, probs = 0.025, na.rm = TRUE),
      quantile(sum_treated_Fever, probs = 0.975, na.rm = TRUE),
      
      quantile(sum_untreated_Fever, probs = 0.50, na.rm = TRUE),
      quantile(sum_untreated_Fever, probs = 0.025, na.rm = TRUE),
      quantile(sum_untreated_Fever, probs = 0.975, na.rm = TRUE),
      
      quantile(sum_total_Fever, probs = 0.50, na.rm = TRUE),
      quantile(sum_total_Fever, probs = 0.025, na.rm = TRUE),
      quantile(sum_total_Fever, probs = 0.975, na.rm = TRUE),
      
      quantile(incidence_Fever, probs = 0.50, na.rm = TRUE),
      quantile(incidence_Fever, probs = 0.025, na.rm = TRUE),
      quantile(incidence_Fever, probs = 0.975, na.rm = TRUE),
      
      quantile(proportion_treated_Fever, probs = 0.50, na.rm = TRUE),
      quantile(proportion_treated_Fever, probs = 0.025, na.rm = TRUE),
      quantile(proportion_treated_Fever, probs = 0.975, na.rm = TRUE)
    )
  )
}

# Combine all sublists
summary_all_years_Fever <- do.call(rbind, summary_list_Fever)

# Combine all years into one data frame
summary_wide_Fever <- summary_all_years_Fever %>%
  mutate(Column = paste0(gsub(" ", "_", Condition), "_", Quantile)) %>%
  dplyr::select(Year, Column, Value) %>%
  pivot_wider(names_from = Column, values_from = Value)

summary_wide_Fever_combined <- summary_wide_Fever %>%
  mutate(
    Treated_Fever = sprintf("%s (%s  %s)", 
                          formatC(Treated_Fever_Median, format = "f", digits = 0, big.mark = ","),
                          formatC(`Treated_Fever_2.5%`, format = "f", digits = 0, big.mark = ","),
                          formatC(`Treated_Fever_97.5%`, format = "f", digits = 0, big.mark = ",")),
    
    Untreated_Fever = sprintf("%s (%s  %s)", 
                            formatC(Untreated_Fever_Median, format = "f", digits = 0, big.mark = ","),
                            formatC(`Untreated_Fever_2.5%`, format = "f", digits = 0, big.mark = ","),
                            formatC(`Untreated_Fever_97.5%`, format = "f", digits = 0, big.mark = ",")),
    
    Total_Fever = sprintf("%s (%s  %s)", 
                        formatC(Total_Fever_Median, format = "f", digits = 0, big.mark = ","),
                        formatC(`Total_Fever_2.5%`, format = "f", digits = 0, big.mark = ","),
                        formatC(`Total_Fever_97.5%`, format = "f", digits = 0, big.mark = ",")),
    
    Incidence_Fever = sprintf("%.2f (%.2f  %.2f)",
                            Incidence_Fever_Median,
                            `Incidence_Fever_2.5%`,
                            `Incidence_Fever_97.5%`
    ),
    Proportion_Fever = sprintf("%.3f (%.3f, %.3f)", 
                             Proportion_Fever_Median, 
                             `Proportion_Fever_2.5%`, 
                             `Proportion_Fever_97.5%`)
  ) %>%
  dplyr::select(Year, Treated_Fever, Untreated_Fever, Total_Fever, Incidence_Fever, Proportion_Fever)


Fever_year_proportion <- summary_wide_Fever %>%
  dplyr::select(Year, Proportion_Fever_Median, `Proportion_Fever_2.5%`, `Proportion_Fever_97.5%`) %>%
  dplyr::rename(Proportion_Median = Proportion_Fever_Median, Proportion_Lower = `Proportion_Fever_2.5%`, Proportion_Upper = `Proportion_Fever_97.5%`)


colnames(summary_wide_LRI_combined) <- c("Year", "Treated_LRI", "Untreated_LRI", "Total_LRI", "Incidence_LRI", "Proportion_treated_LRI")
colnames(summary_wide_Fever_combined) <- c("Year", "Treated_Fever_without_LRI", "Untreated_Fever_without_LRI", "Total_Fever_without_LRI", "Incidence_Fever_without_LRI", "Proportion_treated_Fever_without_LRI")

dir.create(
  "outputs/plot_df",
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  "outputs/table/by_year/probability",
  recursive = TRUE,
  showWarnings = FALSE
)

write.csv(LRI_year_proportion   , "outputs/plot_df/1.3_plot_Ab_LRI_year_proportion.csv")
write.csv(Fever_year_proportion   , "outputs/plot_df/1.3_plot_Ab_Fever_year_proportion.csv")

#### 1.2 CSV file ####
write.csv(summary_wide_LRI_combined, "outputs/table/by_year/probability/LRI.csv")
write.csv(summary_wide_Fever_combined, "outputs/table/by_year/probability/Fever_without_LRI.csv")







#----------------------------1.3 by Malaria and year (treatment) ----------------------------

#------LRI-------
Country_Year_Malaria_meta <- Country_Year_meta %>% inner_join(AbLRI_matrix[,c(1,2,5)], by = c("Country", "Year"))
Malaria_level <- c(0,1)
years <- 2003:2019
summary_list_malaria_LRI <- list()
for (mal in Malaria_level) {
  for (yr in years) {
    
    treated_LRI <- cbind(Country_Year_Malaria_meta, treated_LRI_matrix) %>%
      filter(Malaria == mal, Year == yr) %>%
      dplyr::select(-Country, -Year, -Malaria)
    
    untreated_LRI <- cbind(Country_Year_Malaria_meta, untreated_LRI_matrix) %>%
      filter(Malaria == mal, Year == yr) %>%
      dplyr::select(-Country, -Year, -Malaria)
    
    total_LRI <- cbind(Country_Year_Malaria_meta, Incidence_LRI_matrix) %>%
      filter(Malaria == mal, Year == yr) %>%
      dplyr::select(-Country, -Year, -Malaria)
    
    # sums per simulation
    sum_treated   <- colSums(treated_LRI, na.rm = TRUE)
    sum_untreated <- colSums(untreated_LRI, na.rm = TRUE)
    sum_total     <- colSums(total_LRI, na.rm = TRUE)
    
    # population (fixed, stratified)
    sum_population_LRI <- sum(
      Antibiotics_LRI$Population[
        Antibiotics_LRI$Year == yr &
          Antibiotics_LRI$Malaria == mal
      ],
      na.rm = TRUE
    )
    
    # incidence per simulation
    incidence_LRI <- sum_total / sum_population_LRI
    
    # proportion treated
    proportion_treated <- sum_treated / (sum_treated + sum_untreated)
    
    summary_list_malaria_LRI[[paste(mal, yr, sep = "_")]] <- data.frame(
      Malaria = rep(mal, 15),
      Year = rep(yr, 15),
      Condition = rep(
        c("Treated LRI",
          "Untreated LRI",
          "Total LRI",
          "Incidence LRI",
          "Proportion LRI"),
        each = 3
      ),
      Quantile = rep(c("Median", "2.5%", "97.5%"), times = 5),
      Value = c(
        quantile(sum_treated, 0.50),
        quantile(sum_treated, 0.025),
        quantile(sum_treated, 0.975),
        
        quantile(sum_untreated, 0.50),
        quantile(sum_untreated, 0.025),
        quantile(sum_untreated, 0.975),
        
        quantile(sum_total, 0.50),
        quantile(sum_total, 0.025),
        quantile(sum_total, 0.975),
        
        quantile(incidence_LRI, 0.50),
        quantile(incidence_LRI, 0.025),
        quantile(incidence_LRI, 0.975),
        
        quantile(proportion_treated, 0.50),
        quantile(proportion_treated, 0.025),
        quantile(proportion_treated, 0.975)
      )
    )
  }
}

# Combine all sublists
summary_all_years_LRI <- do.call(rbind, summary_list_malaria_LRI)

# Convert to wide format
summary_wide_LRI <- summary_all_years_LRI %>%
  mutate(Column = paste0(gsub(" ", "_", Condition), "_", Quantile)) %>%
  dplyr::select(Year, Malaria, Column, Value) %>%
  pivot_wider(names_from = Column, values_from = Value)


summary_wide_LRI_by_malaria <- summary_wide_LRI %>%
  mutate(
    Treated_LRI = sprintf("%s (%s  %s)", 
                          formatC(Treated_LRI_Median, format = "f", digits = 0, big.mark = ","),
                          formatC(`Treated_LRI_2.5%`, format = "f", digits = 0, big.mark = ","),
                          formatC(`Treated_LRI_97.5%`, format = "f", digits = 0, big.mark = ",")),
    
    Untreated_LRI = sprintf("%s (%s  %s)", 
                            formatC(Untreated_LRI_Median, format = "f", digits = 0, big.mark = ","),
                            formatC(`Untreated_LRI_2.5%`, format = "f", digits = 0, big.mark = ","),
                            formatC(`Untreated_LRI_97.5%`, format = "f", digits = 0, big.mark = ",")),
    
    Total_LRI = sprintf("%s (%s  %s)", 
                        formatC(Total_LRI_Median, format = "f", digits = 0, big.mark = ","),
                        formatC(`Total_LRI_2.5%`, format = "f", digits = 0, big.mark = ","),
                        formatC(`Total_LRI_97.5%`, format = "f", digits = 0, big.mark = ",")),

    Incidence_LRI = sprintf("%.2f (%.2f  %.2f)",
                            Incidence_LRI_Median,
                            `Incidence_LRI_2.5%`,
                            `Incidence_LRI_97.5%`
    ),
    
    Proportion_LRI = sprintf("%.3f (%.3f, %.3f)", 
                             Proportion_LRI_Median, 
                             `Proportion_LRI_2.5%`, 
                             `Proportion_LRI_97.5%`)
  ) %>%
  dplyr::select(Year, Malaria, Treated_LRI, Untreated_LRI, Total_LRI, Incidence_LRI, Proportion_LRI)


LRI_year_proportion_malaria_0 <- summary_wide_LRI %>%
  filter(Malaria == 0) %>%
  dplyr::select(Year, Malaria, Proportion_LRI_Median, `Proportion_LRI_2.5%`, `Proportion_LRI_97.5%`) %>%
  dplyr::rename(Proportion_Median = Proportion_LRI_Median, Proportion_Lower = `Proportion_LRI_2.5%`, Proportion_Upper = `Proportion_LRI_97.5%`)

LRI_year_proportion_malaria_1 <- summary_wide_LRI %>%
  filter(Malaria == 1) %>%
  dplyr::select(Year, Malaria, Proportion_LRI_Median, `Proportion_LRI_2.5%`, `Proportion_LRI_97.5%`) %>%
  dplyr::rename(Proportion_Median = Proportion_LRI_Median, Proportion_Lower = `Proportion_LRI_2.5%`, Proportion_Upper = `Proportion_LRI_97.5%`)








#------Fever-------
Country_Year_Malaria_meta <- Country_Year_meta %>% inner_join(AbFever_matrix[,c(1,2,5)], by = c("Country", "Year"))
Malaria_level <- c(0,1)
years <- 2003:2019
summary_list_malaria_Fever <- list()
for (mal in Malaria_level) {
  for (yr in years) {
    
    treated_Fever <- cbind(Country_Year_Malaria_meta, treated_Fever_matrix) %>%
      filter(Malaria == mal, Year == yr) %>%
      dplyr::select(-Country, -Year, -Malaria)
    
    untreated_Fever <- cbind(Country_Year_Malaria_meta, untreated_Fever_matrix) %>%
      filter(Malaria == mal, Year == yr) %>%
      dplyr::select(-Country, -Year, -Malaria)
    
    total_Fever <- cbind(Country_Year_Malaria_meta, Incidence_Fever_matrix) %>%
      filter(Malaria == mal, Year == yr) %>%
      dplyr::select(-Country, -Year, -Malaria)
    
    # sums per simulation
    sum_treated   <- colSums(treated_Fever, na.rm = TRUE)
    sum_untreated <- colSums(untreated_Fever, na.rm = TRUE)
    sum_total     <- colSums(total_Fever, na.rm = TRUE)
    
    # population (fixed, stratified)
    sum_population_Fever <- sum(
      Antibiotics_Fever$Population[
        Antibiotics_Fever$Year == yr &
          Antibiotics_Fever$Malaria == mal
      ],
      na.rm = TRUE
    )
    
    # incidence per simulation
    incidence_Fever <- sum_total / sum_population_Fever
    
    # proportion treated
    proportion_treated <- sum_treated / (sum_treated + sum_untreated)
    
    summary_list_malaria_Fever[[paste(mal, yr, sep = "_")]] <- data.frame(
      Malaria = rep(mal, 15),
      Year = rep(yr, 15),
      Condition = rep(
        c("Treated Fever",
          "Untreated Fever",
          "Total Fever",
          "Incidence Fever",
          "Proportion Fever"),
        each = 3
      ),
      Quantile = rep(c("Median", "2.5%", "97.5%"), times = 5),
      Value = c(
        quantile(sum_treated, 0.50),
        quantile(sum_treated, 0.025),
        quantile(sum_treated, 0.975),
        
        quantile(sum_untreated, 0.50),
        quantile(sum_untreated, 0.025),
        quantile(sum_untreated, 0.975),
        
        quantile(sum_total, 0.50),
        quantile(sum_total, 0.025),
        quantile(sum_total, 0.975),
        
        quantile(incidence_Fever, 0.50),
        quantile(incidence_Fever, 0.025),
        quantile(incidence_Fever, 0.975),
        
        quantile(proportion_treated, 0.50),
        quantile(proportion_treated, 0.025),
        quantile(proportion_treated, 0.975)
      )
    )
  }
}

# Combine all sublists
summary_all_years_Fever <- do.call(rbind, summary_list_malaria_Fever)

# Convert to wide format
summary_wide_Fever <- summary_all_years_Fever %>%
  mutate(Column = paste0(gsub(" ", "_", Condition), "_", Quantile)) %>%
  dplyr::select(Year, Malaria, Column, Value) %>%
  pivot_wider(names_from = Column, values_from = Value)


summary_wide_Fever_by_malaria <- summary_wide_Fever %>%
  mutate(
    Treated_Fever = sprintf("%s (%s  %s)", 
                          formatC(Treated_Fever_Median, format = "f", digits = 0, big.mark = ","),
                          formatC(`Treated_Fever_2.5%`, format = "f", digits = 0, big.mark = ","),
                          formatC(`Treated_Fever_97.5%`, format = "f", digits = 0, big.mark = ",")),
    
    Untreated_Fever = sprintf("%s (%s  %s)", 
                            formatC(Untreated_Fever_Median, format = "f", digits = 0, big.mark = ","),
                            formatC(`Untreated_Fever_2.5%`, format = "f", digits = 0, big.mark = ","),
                            formatC(`Untreated_Fever_97.5%`, format = "f", digits = 0, big.mark = ",")),
    
    Total_Fever = sprintf("%s (%s  %s)", 
                        formatC(Total_Fever_Median, format = "f", digits = 0, big.mark = ","),
                        formatC(`Total_Fever_2.5%`, format = "f", digits = 0, big.mark = ","),
                        formatC(`Total_Fever_97.5%`, format = "f", digits = 0, big.mark = ",")),
    
    Incidence_Fever = sprintf("%.2f (%.2f  %.2f)",
                            Incidence_Fever_Median,
                            `Incidence_Fever_2.5%`,
                            `Incidence_Fever_97.5%`
    ),
    
    Proportion_Fever = sprintf("%.3f (%.3f, %.3f)", 
                             Proportion_Fever_Median, 
                             `Proportion_Fever_2.5%`, 
                             `Proportion_Fever_97.5%`)
  ) %>%
  dplyr::select(Year, Malaria, Treated_Fever, Untreated_Fever, Total_Fever, Incidence_Fever, Proportion_Fever)


Fever_year_proportion_malaria_0 <- summary_wide_Fever %>%
  filter(Malaria == 0) %>%
  dplyr::select(Year, Malaria, Proportion_Fever_Median, `Proportion_Fever_2.5%`, `Proportion_Fever_97.5%`) %>%
  dplyr::rename(Proportion_Median = Proportion_Fever_Median, Proportion_Lower = `Proportion_Fever_2.5%`, Proportion_Upper = `Proportion_Fever_97.5%`)

Fever_year_proportion_malaria_1 <- summary_wide_Fever %>%
  filter(Malaria == 1) %>%
  dplyr::select(Year, Malaria, Proportion_Fever_Median, `Proportion_Fever_2.5%`, `Proportion_Fever_97.5%`) %>%
  dplyr::rename(Proportion_Median = Proportion_Fever_Median, Proportion_Lower = `Proportion_Fever_2.5%`, Proportion_Upper = `Proportion_Fever_97.5%`)


colnames(summary_wide_LRI_by_malaria) <- c("Year", "Malaria", "Treated_LRI", "Untreated_LRI", "Total_LRI", "Incidence_LRI", "Proportion_treated_LRI")
colnames(summary_wide_Fever_by_malaria) <- c("Year", "Malaria", "Treated_Fever_without_LRI", "Untreated_Fever_without_LRI", "Total_Fever_without_LRI", "Incidence_Fever", "Proportion_treated_Fever_without_LRI")


dir.create(
  "outputs/table/by_year_malaria/probability",
  recursive = TRUE,
  showWarnings = FALSE
)

# #### 1.4 CSV file ####
write.csv(summary_wide_LRI_by_malaria, "outputs/table/by_year_malaria/probability/LRI.csv")
write.csv(summary_wide_Fever_by_malaria, "outputs/table/by_year_malaria/probability/Fever_wihtout_LRI.csv")
Malaria_distribution <- as.matrix(with(Country_Year_Malaria_meta, addmargins(table(Year, Malaria), margin = c(1,2))))
write.csv(Malaria_distribution, "outputs/table/by_year_malaria/Country_number_by_Malaria_Year.csv")

# save the plot data
write.csv(LRI_year_proportion_malaria_0  , "outputs/plot_df/1.4_plot_Ab_LRI_year_proportion_malaria_0.csv")
write.csv(LRI_year_proportion_malaria_1  , "outputs/plot_df/1.4_plot_Ab_LRI_year_proportion_malaria_1.csv")
write.csv(Fever_year_proportion_malaria_0, "outputs/plot_df/1.4_plot_Ab_Fever_year_proportion_malaria_0.csv")
write.csv(Fever_year_proportion_malaria_1, "outputs/plot_df/1.4_plot_Ab_Fever_year_proportion_malaria_1.csv")










#----------------------------1.4 by Income and year (treatment) ----------------------------
first_sim <- AbLRI_matrix[,c(1,2,4)]
income_dist <- pivot_wider(first_sim, names_from = Year, values_from = Income)


first5_years <- as.character(2003:2007)
last5_years  <- as.character(2015:2019)

income_dist_summary <- income_dist %>%
  rowwise() %>%
  mutate(
    First5 = {
      first_vals <- c_across(all_of(first5_years))
      tab <- table(first_vals)
      if(length(tab) == 0) NA_character_ else names(sort(tab, decreasing = TRUE))[1]
    },
    Last5 = {
      last_vals <- c_across(all_of(last5_years))
      tab <- table(last_vals)
      if(length(tab) == 0) NA_character_ else names(sort(tab, decreasing = TRUE))[1]
    }
  ) %>%
  ungroup()

income_dist_summary[which_na(income_dist_summary$First5),]
income_dist_summary[income_dist_summary$Country == "Djibouti",]$First5 <- "LM" 
AbLRI_matrix_2 <- AbLRI_matrix %>% inner_join(income_dist_summary[,c("Country","First5","Last5")], by = "Country") 

AbFever_matrix_2 <- AbFever_matrix %>% inner_join(income_dist_summary[,c("Country","First5","Last5")], by = "Country") 
rm(income_dist, income_dist_summary, first_sim, first5_years, last5_years)

Country_Year_info <- Country_Year_meta %>% inner_join(AbLRI_matrix_2[,c(1:7,10008:10009)], by = c("Country", "Year"))
Country_Year_info <- Country_Year_info[,c("Country","Year","GDP","Last5","Malaria","GBD_region","Population")]
Country_Year_info_2 <- Country_Year_info %>% dplyr::select(Country, Last5, Malaria, GBD_region, Population) %>% filter(!duplicated(Country))
write.csv(Country_Year_info, "outputs/draft_df/Suppl_3_108LMICs.csv")
write.csv(Country_Year_info_2, "outputs/draft_df/Suppl_3_2_108LMICs.csv")


#------LRI-------
Country_Year_Income_meta <- Country_Year_meta %>% inner_join(AbLRI_matrix_2[,c("Country","Year","Last5")], by = c("Country", "Year"))
Antibiotics_LRI_2 <- Antibiotics_LRI %>% inner_join(AbLRI_matrix_2[,c("Country","Year","Last5")], by = c("Country", "Year"))
Income_level <- c('L','LM','UM')
years <- 2003:2019
summary_list_LRI <- list()

for (inc in Income_level) {
  for (yr in years) {
    # Filter and select simulation matrices
    treated_LRI <- cbind(Country_Year_Income_meta, treated_LRI_matrix) %>%
      filter(Last5 == inc, Year == yr) %>%
      dplyr::select(-Country, -Year, -Last5)
    
    untreated_LRI <- cbind(Country_Year_Income_meta, untreated_LRI_matrix) %>%
      filter(Last5 == inc, Year == yr) %>%
      dplyr::select(-Country, -Year, -Last5)
    
    total_LRI <- cbind(Country_Year_Income_meta, Incidence_LRI_matrix) %>%
      filter(Last5 == inc, Year == yr) %>%
      dplyr::select(-Country, -Year, -Last5)
    
    # Sum over countries for each simulation
    sum_treated <- colSums(treated_LRI, na.rm = TRUE)
    sum_untreated <- colSums(untreated_LRI, na.rm = TRUE)
    sum_total <- colSums(total_LRI, na.rm = TRUE)
    # population (fixed, stratified)
    sum_population_LRI <- sum(
      Antibiotics_LRI_2$Population[
        Antibiotics_LRI_2$Year == yr &
          Antibiotics_LRI_2$Last5 == inc
      ],
      na.rm = TRUE
    )
    
    # incidence per simulation
    incidence_LRI <- sum_total / sum_population_LRI
    
    # proportion treated
    proportion_treated <- sum_treated / (sum_treated + sum_untreated)
    
    # Store results as long-format table
    summary_list_LRI[[paste(inc, yr, sep = "_")]] <- data.frame(
      Last5 = inc,
      Year = rep(yr, 15),
      Condition = rep(c("Treated LRI", "Untreated LRI", "Total LRI", "Incidence LRI", "Proportion LRI"), each = 3),
      Quantile = rep(c("Median", "2.5%", "97.5%"), times = 5),
      Value = c(
        quantile(sum_treated, probs = 0.50, na.rm = TRUE),
        quantile(sum_treated, probs = 0.025, na.rm = TRUE),
        quantile(sum_treated, probs = 0.975, na.rm = TRUE),
        
        quantile(sum_untreated, probs = 0.50, na.rm = TRUE),
        quantile(sum_untreated, probs = 0.025, na.rm = TRUE),
        quantile(sum_untreated, probs = 0.975, na.rm = TRUE),
        
        quantile(sum_total, probs = 0.50, na.rm = TRUE),
        quantile(sum_total, probs = 0.025, na.rm = TRUE),
        quantile(sum_total, probs = 0.975, na.rm = TRUE),
        
        quantile(incidence_LRI, 0.50),
        quantile(incidence_LRI, 0.025),
        quantile(incidence_LRI, 0.975),
        
        quantile(proportion_treated, probs = 0.50, na.rm = TRUE),
        quantile(proportion_treated, probs = 0.025, na.rm = TRUE),
        quantile(proportion_treated, probs = 0.975, na.rm = TRUE)
      )
    )
  }
}

# Combine all sublists
summary_all_years_LRI <- do.call(rbind, summary_list_LRI)

# Convert to wide format
summary_wide_LRI <- summary_all_years_LRI %>%
  mutate(Column = paste0(gsub(" ", "_", Condition), "_", Quantile)) %>%
  dplyr::select(Year, Last5, Column, Value) %>%
  pivot_wider(names_from = Column, values_from = Value)


summary_wide_LRI_by_income <- summary_wide_LRI %>%
  mutate(
    Treated_LRI = sprintf("%s (%s  %s)", 
                          formatC(Treated_LRI_Median, format = "f", digits = 0, big.mark = ","),
                          formatC(`Treated_LRI_2.5%`, format = "f", digits = 0, big.mark = ","),
                          formatC(`Treated_LRI_97.5%`, format = "f", digits = 0, big.mark = ",")),
    
    Untreated_LRI = sprintf("%s (%s  %s)", 
                            formatC(Untreated_LRI_Median, format = "f", digits = 0, big.mark = ","),
                            formatC(`Untreated_LRI_2.5%`, format = "f", digits = 0, big.mark = ","),
                            formatC(`Untreated_LRI_97.5%`, format = "f", digits = 0, big.mark = ",")),
    
    Total_LRI = sprintf("%s (%s  %s)", 
                        formatC(Total_LRI_Median, format = "f", digits = 0, big.mark = ","),
                        formatC(`Total_LRI_2.5%`, format = "f", digits = 0, big.mark = ","),
                        formatC(`Total_LRI_97.5%`, format = "f", digits = 0, big.mark = ",")),
    
    Incidence_LRI = sprintf("%.2f (%.2f  %.2f)",
                            Incidence_LRI_Median,
                            `Incidence_LRI_2.5%`,
                            `Incidence_LRI_97.5%`),
    
    Proportion_LRI = sprintf("%.3f (%.3f, %.3f)", 
                             Proportion_LRI_Median, 
                             `Proportion_LRI_2.5%`, 
                             `Proportion_LRI_97.5%`)
  ) %>%
  dplyr::select(Year, Last5, Treated_LRI, Untreated_LRI, Total_LRI, Incidence_LRI, Proportion_LRI)


LRI_year_proportion_L <- summary_wide_LRI %>%
  filter(Last5 == 'L') %>%
  dplyr::select(Year, Last5, Proportion_LRI_Median, `Proportion_LRI_2.5%`, `Proportion_LRI_97.5%`) %>%
  dplyr::rename(Proportion_Median = Proportion_LRI_Median, Proportion_Lower = `Proportion_LRI_2.5%`, Proportion_Upper = `Proportion_LRI_97.5%`)

LRI_year_proportion_LM <- summary_wide_LRI %>%
  filter(Last5 == 'LM') %>%
  dplyr::select(Year, Last5, Proportion_LRI_Median, `Proportion_LRI_2.5%`, `Proportion_LRI_97.5%`) %>%
  dplyr::rename(Proportion_Median = Proportion_LRI_Median, Proportion_Lower = `Proportion_LRI_2.5%`, Proportion_Upper = `Proportion_LRI_97.5%`)

LRI_year_proportion_UM <- summary_wide_LRI %>%
  filter(Last5 == 'UM') %>%
  dplyr::select(Year, Last5, Proportion_LRI_Median, `Proportion_LRI_2.5%`, `Proportion_LRI_97.5%`) %>%
  dplyr::rename(Proportion_Median = Proportion_LRI_Median, Proportion_Lower = `Proportion_LRI_2.5%`, Proportion_Upper = `Proportion_LRI_97.5%`)




#------Fever-------
Country_Year_Income_meta <- Country_Year_meta %>% inner_join(AbFever_matrix_2[,c("Country","Year","Malaria","Last5")], by = c("Country", "Year"))
Antibiotics_Fever_2 <- Antibiotics_Fever %>% inner_join(AbFever_matrix_2[,c("Country","Year","Last5")], by = c("Country", "Year"))
Income_level <- c('L','LM','UM')
Malaria_level <- c(0,1)
years <- 2003:2019
summary_list_Fever <- list()

for (mal in Malaria_level) {
  for (inc in Income_level) {
    for (yr in years) {
      # Filter and select simulation matrices
      treated_Fever <- cbind(Country_Year_Income_meta, treated_Fever_matrix) %>%
        filter(Malaria == mal, Last5 == inc, Year == yr) %>%
        dplyr::select(-Country, -Year, -Malaria, -Last5)
      
      untreated_Fever <- cbind(Country_Year_Income_meta, untreated_Fever_matrix) %>%
        filter(Malaria == mal, Last5 == inc, Year == yr) %>%
        dplyr::select(-Country, -Year, -Malaria, -Last5)
      
      total_Fever <- cbind(Country_Year_Income_meta, Incidence_Fever_matrix) %>%
        filter(Malaria == mal, Last5 == inc, Year == yr) %>%
        dplyr::select(-Country, -Year, -Malaria, -Last5)
      
      # Sum over countries for each simulation
      sum_treated <- colSums(treated_Fever, na.rm = TRUE)
      sum_untreated <- colSums(untreated_Fever, na.rm = TRUE)
      sum_total <- colSums(total_Fever, na.rm = TRUE)
      # population (fixed, stratified)
      sum_population_Fever <- sum(
        Antibiotics_Fever_2$Population[
          Antibiotics_Fever_2$Year == yr &
            Antibiotics_Fever_2$Malaria == mal &
              Antibiotics_Fever_2$Last5 == inc
        ],
        na.rm = TRUE
      )
      
      # incidence per simulation
      incidence_Fever <- sum_total / sum_population_Fever
      
      # proportion treated
      proportion_treated <- sum_treated / (sum_treated + sum_untreated)
      
      # Store results as long-format table
      summary_list_Fever[[paste(mal, inc, yr, sep = "_")]] <- data.frame(
        Malaria = mal,
        Last5 = inc,
        Year = rep(yr, 15),
        Condition = rep(c("Treated Fever", "Untreated Fever", "Total Fever", "Incidence Fever", "Proportion Fever"), each = 3),
        Quantile = rep(c("Median", "2.5%", "97.5%"), times = 5),
        Value = c(
          quantile(sum_treated, probs = 0.50, na.rm = TRUE),
          quantile(sum_treated, probs = 0.025, na.rm = TRUE),
          quantile(sum_treated, probs = 0.975, na.rm = TRUE),
          
          quantile(sum_untreated, probs = 0.50, na.rm = TRUE),
          quantile(sum_untreated, probs = 0.025, na.rm = TRUE),
          quantile(sum_untreated, probs = 0.975, na.rm = TRUE),
          
          quantile(sum_total, probs = 0.50, na.rm = TRUE),
          quantile(sum_total, probs = 0.025, na.rm = TRUE),
          quantile(sum_total, probs = 0.975, na.rm = TRUE),
          
          quantile(incidence_Fever, 0.50),
          quantile(incidence_Fever, 0.025),
          quantile(incidence_Fever, 0.975),
          
          quantile(proportion_treated, probs = 0.50, na.rm = TRUE),
          quantile(proportion_treated, probs = 0.025, na.rm = TRUE),
          quantile(proportion_treated, probs = 0.975, na.rm = TRUE)
        )
      )
    }
  }
}  

# Combine all sublists
summary_all_years_Fever <- do.call(rbind, summary_list_Fever)

# Convert to wide format
summary_wide_Fever <- summary_all_years_Fever %>%
  mutate(Column = paste0(gsub(" ", "_", Condition), "_", Quantile)) %>%
  dplyr::select(Year, Malaria, Last5, Column, Value) %>%
  pivot_wider(names_from = Column, values_from = Value)


summary_wide_Fever_by_income <- summary_wide_Fever %>%
  mutate(
    Treated_Fever = sprintf("%s (%s  %s)", 
                            formatC(Treated_Fever_Median, format = "f", digits = 0, big.mark = ","),
                            formatC(`Treated_Fever_2.5%`, format = "f", digits = 0, big.mark = ","),
                            formatC(`Treated_Fever_97.5%`, format = "f", digits = 0, big.mark = ",")),
    
    Untreated_Fever = sprintf("%s (%s  %s)", 
                              formatC(Untreated_Fever_Median, format = "f", digits = 0, big.mark = ","),
                              formatC(`Untreated_Fever_2.5%`, format = "f", digits = 0, big.mark = ","),
                              formatC(`Untreated_Fever_97.5%`, format = "f", digits = 0, big.mark = ",")),
    
    Total_Fever = sprintf("%s (%s  %s)", 
                          formatC(Total_Fever_Median, format = "f", digits = 0, big.mark = ","),
                          formatC(`Total_Fever_2.5%`, format = "f", digits = 0, big.mark = ","),
                          formatC(`Total_Fever_97.5%`, format = "f", digits = 0, big.mark = ",")),
    
    Incidence_Fever = sprintf("%.2f (%.2f  %.2f)",
                            Incidence_Fever_Median,
                            `Incidence_Fever_2.5%`,
                            `Incidence_Fever_97.5%`),
    
    Proportion_Fever = sprintf("%.3f (%.3f, %.3f)", 
                               Proportion_Fever_Median, 
                               `Proportion_Fever_2.5%`, 
                               `Proportion_Fever_97.5%`)
  ) %>%
  dplyr::select(Year, Malaria, Last5, Treated_Fever, Untreated_Fever, Total_Fever, Incidence_Fever, Proportion_Fever)


Fever_year_proportion_mal0_L <- summary_wide_Fever %>%
  filter(Malaria == 0, Last5 == 'L') %>%  
  
  dplyr::select(Year, Last5, Proportion_Fever_Median, `Proportion_Fever_2.5%`, `Proportion_Fever_97.5%`) %>%
  dplyr::rename(Proportion_Median = Proportion_Fever_Median, Proportion_Lower = `Proportion_Fever_2.5%`, Proportion_Upper = `Proportion_Fever_97.5%`)

Fever_year_proportion_mal0_LM <- summary_wide_Fever %>%
  filter(Malaria == 0, Last5 == 'LM') %>%
  dplyr::select(Year, Last5, Proportion_Fever_Median, `Proportion_Fever_2.5%`, `Proportion_Fever_97.5%`) %>%
  dplyr::rename(Proportion_Median = Proportion_Fever_Median, Proportion_Lower = `Proportion_Fever_2.5%`, Proportion_Upper = `Proportion_Fever_97.5%`)

Fever_year_proportion_mal0_UM <- summary_wide_Fever %>%
  filter(Malaria == 0, Last5 == 'UM') %>%
  dplyr::select(Year, Last5, Proportion_Fever_Median, `Proportion_Fever_2.5%`, `Proportion_Fever_97.5%`) %>%
  dplyr::rename(Proportion_Median = Proportion_Fever_Median, Proportion_Lower = `Proportion_Fever_2.5%`, Proportion_Upper = `Proportion_Fever_97.5%`)


Fever_year_proportion_mal1_L <- summary_wide_Fever %>%
  filter(Malaria == 1, Last5 == 'L') %>%
  dplyr::select(Year, Last5, Proportion_Fever_Median, `Proportion_Fever_2.5%`, `Proportion_Fever_97.5%`) %>%
  dplyr::rename(Proportion_Median = Proportion_Fever_Median, Proportion_Lower = `Proportion_Fever_2.5%`, Proportion_Upper = `Proportion_Fever_97.5%`)

Fever_year_proportion_mal1_LM <- summary_wide_Fever %>%
  filter(Malaria == 1, Last5 == 'LM') %>%
  dplyr::select(Year, Last5, Proportion_Fever_Median, `Proportion_Fever_2.5%`, `Proportion_Fever_97.5%`) %>%
  dplyr::rename(Proportion_Median = Proportion_Fever_Median, Proportion_Lower = `Proportion_Fever_2.5%`, Proportion_Upper = `Proportion_Fever_97.5%`)

Fever_year_proportion_mal1_UM <- summary_wide_Fever %>%
  filter(Malaria == 1, Last5 == 'UM') %>%
  dplyr::select(Year, Last5, Proportion_Fever_Median, `Proportion_Fever_2.5%`, `Proportion_Fever_97.5%`) %>%
  dplyr::rename(Proportion_Median = Proportion_Fever_Median, Proportion_Lower = `Proportion_Fever_2.5%`, Proportion_Upper = `Proportion_Fever_97.5%`)

dir.create(
  "outputs/table/by_year_income/probability",
  recursive = TRUE,
  showWarnings = FALSE
)

# #### 1.3 CSV file ####
write.csv(summary_wide_LRI_by_income, "outputs/table/by_year_income/probability/LRI.csv")
write.csv(summary_wide_Fever_by_income, "outputs/table/by_year_income/probability/Fever_wihtout_LRI.csv")
Income_distribution <- as.matrix(with(Country_Year_Income_meta, addmargins(table(Year, Last5), margin = c(1,2))))
write.csv(Income_distribution, "outputs/table/by_year_income/Country_number_by_Income_Year.csv")

# save the plot data
write.csv(LRI_year_proportion_L   , "outputs/plot_df/1.3_plot_Ab_LRI_year_proportion_L.csv")
write.csv(LRI_year_proportion_LM  , "outputs/plot_df/1.3_plot_Ab_LRI_year_proportion_LM.csv")
write.csv(LRI_year_proportion_UM  , "outputs/plot_df/1.3_plot_Ab_LRI_year_proportion_UM.csv")
write.csv(Fever_year_proportion_mal0_L , "outputs/plot_df/1.3_plot_Ab_Fever_year_proportion_mal0_L.csv")
write.csv(Fever_year_proportion_mal0_LM, "outputs/plot_df/1.3_plot_Ab_Fever_year_proportion_mal0_LM.csv")
write.csv(Fever_year_proportion_mal0_UM, "outputs/plot_df/1.3_plot_Ab_Fever_year_proportion_mal0_UM.csv")
write.csv(Fever_year_proportion_mal1_L , "outputs/plot_df/1.3_plot_Ab_Fever_year_proportion_mal1_L.csv")
write.csv(Fever_year_proportion_mal1_LM, "outputs/plot_df/1.3_plot_Ab_Fever_year_proportion_mal1_LM.csv")
write.csv(Fever_year_proportion_mal1_UM, "outputs/plot_df/1.3_plot_Ab_Fever_year_proportion_mal1_UM.csv")











#----------------------------1.5 GBD_region (treatment) ----------------------------
#------LRI-------
Country_Year_GBD_region_meta <- Country_Year_meta %>% inner_join(AbLRI_matrix[,c(1,2,6)], by = c("Country", "Year"))
GBD_region_level <- c("North Africa and Middle East", "Central Europe, Eastern Europe, and Central Asia", "Sub-Saharan Africa", "High-income", "South Asia", "Latin America and Caribbean", "Southeast Asia, East Asia, and Oceania")

years <- 2003:2019
summary_list_GBD_region_LRI <- list()

for (gbd in GBD_region_level) {
  for (yr in years) {
    # Filter and select simulation matrices
    treated_LRI <- cbind(Country_Year_GBD_region_meta, treated_LRI_matrix) %>%
      filter(GBD_region == gbd, Year == yr) %>%
      dplyr::select(-Country, -Year, -GBD_region)
    
    untreated_LRI <- cbind(Country_Year_GBD_region_meta, untreated_LRI_matrix) %>%
      filter(GBD_region == gbd, Year == yr) %>%
      dplyr::select(-Country, -Year, -GBD_region)
    
    total_LRI <- cbind(Country_Year_GBD_region_meta, Incidence_LRI_matrix) %>%
      filter(GBD_region == gbd, Year == yr) %>%
      dplyr::select(-Country, -Year, -GBD_region)
    
    # Sum over countries for each simulation
    sum_treated <- colSums(treated_LRI, na.rm = TRUE)
    sum_untreated <- colSums(untreated_LRI, na.rm = TRUE)
    sum_total <- colSums(total_LRI, na.rm = TRUE)
    # population (fixed, stratified)
    sum_population_LRI <- sum(
      Antibiotics_LRI$Population[
        Antibiotics_LRI$Year == yr &
          Antibiotics_LRI$GBD_region == gbd
      ],
      na.rm = TRUE
    )
    
    # incidence per simulation
    incidence_LRI <- sum_total / sum_population_LRI

    proportion_treated <- sum_treated / (sum_treated + sum_untreated)
    
    # Store results as long-format table
    summary_list_GBD_region_LRI[[paste(gbd, yr, sep = "_")]] <- data.frame(
      GBD_region = gbd,
      Year = rep(yr, 15),
      Condition = rep(c("Treated LRI", "Untreated LRI", "Total LRI", "Incidence LRI", "Proportion LRI"), each = 3),
      Quantile = rep(c("Median", "2.5%", "97.5%"), times = 5),
      Value = c(
        quantile(sum_treated, probs = 0.50, na.rm = TRUE),
        quantile(sum_treated, probs = 0.025, na.rm = TRUE),
        quantile(sum_treated, probs = 0.975, na.rm = TRUE),
        
        quantile(sum_untreated, probs = 0.50, na.rm = TRUE),
        quantile(sum_untreated, probs = 0.025, na.rm = TRUE),
        quantile(sum_untreated, probs = 0.975, na.rm = TRUE),
        
        quantile(sum_total, probs = 0.50, na.rm = TRUE),
        quantile(sum_total, probs = 0.025, na.rm = TRUE),
        quantile(sum_total, probs = 0.975, na.rm = TRUE),
        
        quantile(incidence_LRI, 0.50),
        quantile(incidence_LRI, 0.025),
        quantile(incidence_LRI, 0.975),
        
        quantile(proportion_treated, probs = 0.50, na.rm = TRUE),
        quantile(proportion_treated, probs = 0.025, na.rm = TRUE),
        quantile(proportion_treated, probs = 0.975, na.rm = TRUE)
      )
    )
  }
}

# Combine all sublists
summary_all_years_LRI <- do.call(rbind, summary_list_GBD_region_LRI)

# Convert to wide format
summary_wide_LRI <- summary_all_years_LRI %>%
  mutate(Column = paste0(gsub(" ", "_", Condition), "_", Quantile)) %>%
  dplyr::select(Year, GBD_region, Column, Value) %>%
  pivot_wider(names_from = Column, values_from = Value)


summary_wide_LRI_by_GBD_region <- summary_wide_LRI %>%
  mutate(
    Treated_LRI = sprintf("%s (%s  %s)", 
                          formatC(Treated_LRI_Median, format = "f", digits = 0, big.mark = ","),
                          formatC(`Treated_LRI_2.5%`, format = "f", digits = 0, big.mark = ","),
                          formatC(`Treated_LRI_97.5%`, format = "f", digits = 0, big.mark = ",")),
    
    Untreated_LRI = sprintf("%s (%s  %s)", 
                            formatC(Untreated_LRI_Median, format = "f", digits = 0, big.mark = ","),
                            formatC(`Untreated_LRI_2.5%`, format = "f", digits = 0, big.mark = ","),
                            formatC(`Untreated_LRI_97.5%`, format = "f", digits = 0, big.mark = ",")),
    
    Total_LRI = sprintf("%s (%s  %s)", 
                        formatC(Total_LRI_Median, format = "f", digits = 0, big.mark = ","),
                        formatC(`Total_LRI_2.5%`, format = "f", digits = 0, big.mark = ","),
                        formatC(`Total_LRI_97.5%`, format = "f", digits = 0, big.mark = ",")),
    
    Incidence_LRI = sprintf("%.2f (%.2f  %.2f)",
                            Incidence_LRI_Median,
                            `Incidence_LRI_2.5%`,
                            `Incidence_LRI_97.5%`),
    
    Proportion_LRI = sprintf("%.3f (%.3f, %.3f)", 
                             Proportion_LRI_Median, 
                             `Proportion_LRI_2.5%`, 
                             `Proportion_LRI_97.5%`)
  ) %>%
  dplyr::select(Year, GBD_region, Treated_LRI, Untreated_LRI, Total_LRI, Incidence_LRI, Proportion_LRI)


for (i in 1:length(GBD_region_level)) {
  temp <- summary_wide_LRI %>%
    filter(GBD_region == GBD_region_level[i]) %>%
    dplyr::select(Year, GBD_region, Proportion_LRI_Median, `Proportion_LRI_2.5%`, `Proportion_LRI_97.5%`) %>%
    rename(
      Proportion_Median = Proportion_LRI_Median,
      Proportion_Lower = `Proportion_LRI_2.5%`,
      Proportion_Upper = `Proportion_LRI_97.5%`
    )
  obj_name <- paste0("LRI_year_proportion_GBD_region_", i)
  assign(obj_name, temp)
}

for (i in 1:7) {
  obj_name <- paste0("LRI_year_proportion_GBD_region_", i)
  path_LRI <- paste0("outputs/plot_df/1.5_plot_Ab_LRI_year_proportion_GBD_region_", i,".csv")
    
  write.csv(get(obj_name), path_LRI, row.names = FALSE)
}






#------Fever-------
Country_Year_GBD_region_meta <- Country_Year_meta %>% inner_join(AbFever_matrix[,c("Country","Year","Malaria","GBD_region")], by = c("Country", "Year"))
GBD_region_level <- c("North Africa and Middle East", "Central Europe, Eastern Europe, and Central Asia", "Sub-Saharan Africa", "High-income", "South Asia", "Latin America and Caribbean", "Southeast Asia, East Asia, and Oceania")
Malaria_level <- c(0,1)
years <- 2003:2019
summary_list_GBD_region_Fever <- list()

for (mal in Malaria_level) {
  for (gbd in GBD_region_level) {
    for (yr in years) {
      # Filter and select simulation matrices
      treated_Fever <- cbind(Country_Year_GBD_region_meta, treated_Fever_matrix) %>%
        filter(Malaria == mal, GBD_region == gbd, Year == yr) %>%
        dplyr::select(-Country, -Year, -Malaria, -GBD_region)
      
      untreated_Fever <- cbind(Country_Year_GBD_region_meta, untreated_Fever_matrix) %>%
        filter(Malaria == mal, GBD_region == gbd, Year == yr) %>%
        dplyr::select(-Country, -Year, -Malaria, -GBD_region)
      
      total_Fever <- cbind(Country_Year_GBD_region_meta, Incidence_Fever_matrix) %>%
        filter(Malaria == mal, GBD_region == gbd, Year == yr) %>%
        dplyr::select(-Country, -Year, -Malaria, -GBD_region)
      
      # Sum over countries for each simulation
      sum_treated <- colSums(treated_Fever, na.rm = TRUE)
      sum_untreated <- colSums(untreated_Fever, na.rm = TRUE)
      sum_total <- colSums(total_Fever, na.rm = TRUE)
      # population (fixed, stratified)
      sum_population_Fever <- sum(
        Antibiotics_Fever$Population[
          Antibiotics_Fever$Malaria == mal &
            Antibiotics_Fever$Year == yr &
              Antibiotics_Fever$GBD_region == gbd
        ],
        na.rm = TRUE
      )
      
      # incidence per simulation
      incidence_Fever <- sum_total / sum_population_Fever
      
      proportion_treated <- sum_treated / (sum_treated + sum_untreated)
      
      # Store results as long-format table
      summary_list_GBD_region_Fever[[paste(mal, gbd, yr, sep = "_")]] <- data.frame(
        Malaria = mal,
        GBD_region = gbd,
        Year = rep(yr, 15),
        Condition = rep(c("Treated Fever", "Untreated Fever", "Total Fever", "Incidence Fever", "Proportion Fever"), each = 3),
        Quantile = rep(c("Median", "2.5%", "97.5%"), times = 5),
        Value = c(
          quantile(sum_treated, probs = 0.50, na.rm = TRUE),
          quantile(sum_treated, probs = 0.025, na.rm = TRUE),
          quantile(sum_treated, probs = 0.975, na.rm = TRUE),
          
          quantile(sum_untreated, probs = 0.50, na.rm = TRUE),
          quantile(sum_untreated, probs = 0.025, na.rm = TRUE),
          quantile(sum_untreated, probs = 0.975, na.rm = TRUE),
          
          quantile(sum_total, probs = 0.50, na.rm = TRUE),
          quantile(sum_total, probs = 0.025, na.rm = TRUE),
          quantile(sum_total, probs = 0.975, na.rm = TRUE),
          
          quantile(incidence_Fever, 0.50, na.rm = TRUE),
          quantile(incidence_Fever, 0.025, na.rm = TRUE),
          quantile(incidence_Fever, 0.975, na.rm = TRUE),
          
          quantile(proportion_treated, probs = 0.50, na.rm = TRUE),
          quantile(proportion_treated, probs = 0.025, na.rm = TRUE),
          quantile(proportion_treated, probs = 0.975, na.rm = TRUE)
        )
      )
    }
  }
}  

# Combine all sublists
summary_all_years_Fever <- do.call(rbind, summary_list_GBD_region_Fever)

# Convert to wide format
summary_wide_Fever <- summary_all_years_Fever %>%
  mutate(Column = paste0(gsub(" ", "_", Condition), "_", Quantile)) %>%
  dplyr::select(Year, Malaria, GBD_region, Column, Value) %>%
  pivot_wider(names_from = Column, values_from = Value)


summary_wide_Fever_by_GBD_region <- summary_wide_Fever %>%
  mutate(
    Treated_Fever = sprintf("%s (%s  %s)", 
                            formatC(Treated_Fever_Median, format = "f", digits = 0, big.mark = ","),
                            formatC(`Treated_Fever_2.5%`, format = "f", digits = 0, big.mark = ","),
                            formatC(`Treated_Fever_97.5%`, format = "f", digits = 0, big.mark = ",")),
    
    Untreated_Fever = sprintf("%s (%s  %s)", 
                              formatC(Untreated_Fever_Median, format = "f", digits = 0, big.mark = ","),
                              formatC(`Untreated_Fever_2.5%`, format = "f", digits = 0, big.mark = ","),
                              formatC(`Untreated_Fever_97.5%`, format = "f", digits = 0, big.mark = ",")),
    
    Total_Fever = sprintf("%s (%s  %s)", 
                          formatC(Total_Fever_Median, format = "f", digits = 0, big.mark = ","),
                          formatC(`Total_Fever_2.5%`, format = "f", digits = 0, big.mark = ","),
                          formatC(`Total_Fever_97.5%`, format = "f", digits = 0, big.mark = ",")),
 
    Incidence_Fever = sprintf("%.2f (%.2f  %.2f)",
                            Incidence_Fever_Median,
                            `Incidence_Fever_2.5%`,
                            `Incidence_Fever_97.5%`),
    
    Proportion_Fever = sprintf("%.3f (%.3f, %.3f)", 
                               Proportion_Fever_Median, 
                               `Proportion_Fever_2.5%`, 
                               `Proportion_Fever_97.5%`)
  ) %>%
  dplyr::select(Year, Malaria, GBD_region, Treated_Fever, Untreated_Fever, Total_Fever, Incidence_Fever, Proportion_Fever)

dir.create(
  "outputs/table/by_year_GBD_region/probability",
  recursive = TRUE,
  showWarnings = FALSE
)

write.csv(summary_wide_LRI_by_GBD_region, "outputs/table/by_year_GBD_region/probability/LRI.csv")
write.csv(summary_wide_Fever_by_GBD_region, "outputs/table/by_year_GBD_region/probability/Fever_wihtout_LRI.csv")


for (mal in Malaria_level){
for (i in 1:length(GBD_region_level)) {
  temp <- summary_wide_Fever %>%
    filter(Malaria == mal, GBD_region == GBD_region_level[i]) %>%
    dplyr::select(Year, GBD_region, Proportion_Fever_Median, `Proportion_Fever_2.5%`, `Proportion_Fever_97.5%`) %>%
    rename(
      Proportion_Median = Proportion_Fever_Median,
      Proportion_Lower = `Proportion_Fever_2.5%`,
      Proportion_Upper = `Proportion_Fever_97.5%`
    )
  obj_name <- paste0("Fever_year_proportion_GBD_region_", mal, "_", i)
  assign(obj_name, temp)
  }
}



for (mal in Malaria_level) {
  for (i in 1:7) {
  obj_Fever <- paste0("Fever_year_proportion_GBD_region_", mal,"_", i)
  path_Fever <- paste0("outputs/plot_df/1.5_plot_Ab_Fever_year_proportion_GBD_region_", mal,"_", i,".csv")

  write.csv(get(obj_Fever), path_Fever, row.names = FALSE)
  }
}  



# #### 1.5 CSV file ####
GBD_region_distribution <- as.matrix(with(Country_Year_GBD_region_meta, addmargins(table(Year, GBD_region), margin = c(1,2))))
write.csv(GBD_region_distribution, "outputs/table/by_year_GBD_region/Country_number_by_GBD_region_Year.csv")











#___________________________________________________-----------------------

#-----------------------2.1 by country and year (DDD)------------------------

# Create matching indices
Year_index <- rep(UNI_indexed$Year, times = ncol(AbLRI_sim))
Country_index <- rep(UNI_indexed$Country, times = ncol(AbLRI_sim))
GDP_index <- rep(UNI_indexed$GDP, times = ncol(AbLRI_sim))
Income_index <- rep(UNI_indexed$Income, times = ncol(AbLRI_sim))


# Combine into dataframe
sim_df_LRI_DDD <- data.frame(
  Country = UNI_indexed$Country,
  Year = UNI_indexed$Year,
  GDP = GDP_index,
  Income = Income_index,
  LRI_DDD = LRI_DDD_long,
  sim = rep(1:10000, each = nrow(LRI_DDD_sim))
)

# Combine into dataframe
sim_df_Fever_DDD <- data.frame(
  Country = UNI_indexed$Country,
  Year = UNI_indexed$Year,
  GDP = GDP_index,
  Income = Income_index,
  Fever_DDD = Fever_DDD_long,
  sim = rep(1:10000, each = nrow(Fever_DDD_sim))
)

# Combine into dataframe
sim_df_all_DDD <- data.frame(
  Country = UNI_indexed$Country,
  Year = UNI_indexed$Year,
  GDP = GDP_index,
  Income = Income_index,
  all_DDD = all_DDD_long,
  sim = rep(1:10000, each = nrow(all_DDD_sim))
)

LRI_DDD_matrix <- sim_df_LRI_DDD %>% pivot_wider(names_from = sim, values_from = LRI_DDD)
Fever_DDD_matrix <- sim_df_Fever_DDD %>% pivot_wider(names_from = sim, values_from = Fever_DDD)
all_DDD_matrix <- sim_df_all_DDD %>% pivot_wider(names_from = sim, values_from = all_DDD)

LRI_DDD_matrix <- LRI_DDD_matrix %>% mutate(Country = ifelse(Country %in% names(country_mapping), country_mapping[Country], Country))
Fever_DDD_matrix <- Fever_DDD_matrix %>% mutate(Country = ifelse(Country %in% names(country_mapping), country_mapping[Country], Country))
all_DDD_matrix <- all_DDD_matrix %>% mutate(Country = ifelse(Country %in% names(country_mapping), country_mapping[Country], Country))

DDD_LRI <- Country_Year_meta %>% inner_join(LRI_DDD_matrix, by = c("Country","Year"))
DDD_LRI_matrix <- DDD_LRI[,5:(n_sim+4)]
DDD_Fever <- Country_Year_meta %>% inner_join(Fever_DDD_matrix, by = c("Country","Year"))
DDD_Fever_matrix <- DDD_Fever[,5:(n_sim+4)]
DDD_all <- Country_Year_meta %>% inner_join(all_DDD_matrix, by = c("Country","Year"))
DDD_all_matrix <- DDD_all[,5:(n_sim+4)]

total_DDD_LRI_matrix              <- Incidence_LRI_matrix * DDD_LRI_matrix
total_DDD_Fever_matrix            <- Incidence_Fever_matrix * DDD_Fever_matrix
total_DDD_all_matrix              <- Incidence_all_matrix * DDD_all_matrix
treated_LRI_matrix                <- Incidence_LRI_matrix*Antibiotics_LRI_matrix
treated_Fever_matrix              <- Incidence_Fever_matrix*Antibiotics_Fever_matrix
treated_all_matrix                <- Incidence_LRI_matrix*Antibiotics_LRI_matrix + Incidence_Fever_matrix*Antibiotics_Fever_matrix

DDD_per_treated_LRI               <- DDD_LRI_matrix/Antibiotics_LRI_matrix
DDD_per_treated_Fever             <- DDD_Fever_matrix/Antibiotics_Fever_matrix
DDD_per_treated_all               <- total_DDD_all_matrix /treated_all_matrix 

make_summary_2 <- function(meta, income, total, per_treated) {

  data.frame(
    Year = meta$Year,
    Country = meta$Country,
    Income = income$Income,
    
    Total = sprintf("%s (%s  %s)",
                      formatC(apply(total, 1, median, na.rm = TRUE), format = "f", big.mark = ",", digits = 0),
                      formatC(apply(total, 1, quantile, probs = 0.025, na.rm = TRUE), format = "f", big.mark = ",", digits = 0),
                      formatC(apply(total, 1, quantile, probs = 0.975, na.rm = TRUE), format = "f", big.mark = ",", digits = 0)),
    
    Per_treated = sprintf("%s (%s  %s)",
                        formatC(apply(per_treated, 1, median, na.rm = TRUE), format = "f", big.mark = ",", digits = 2),
                        formatC(apply(per_treated, 1, quantile, probs = 0.025, na.rm = TRUE), format = "f", big.mark = ",", digits = 2),
                        formatC(apply(per_treated, 1, quantile, probs = 0.975, na.rm = TRUE), format = "f", big.mark = ",", digits = 2))
  )
}

make_df_2 <- function(meta, income, total, per_treated) {
  data.frame(
    Year = meta$Year,
    Country = meta$Country,
    Income = income$Income,
    
    Total_median = apply(total, 1, median, na.rm = TRUE),
    Total_lower  = apply(total, 1, quantile, probs = 0.025, na.rm = TRUE),
    Total_upper  = apply(total, 1, quantile, probs = 0.975, na.rm = TRUE),
    
    Per_treated_median = apply(per_treated, 1, median, na.rm = TRUE),
    Per_treated_lower  = apply(per_treated, 1, quantile, probs = 0.025, na.rm = TRUE),
    Per_treated_upper  = apply(per_treated, 1, quantile, probs = 0.975, na.rm = TRUE)
  )
}


LRI_DDD_summary          <- make_summary_2(Country_Year_meta, DDD_LRI, total_DDD_LRI_matrix              , DDD_per_treated_LRI         )
Fever_DDD_summary        <- make_summary_2(Country_Year_meta, DDD_LRI, total_DDD_Fever_matrix            , DDD_per_treated_Fever       )
all_DDD_summary          <- make_summary_2(Country_Year_meta, DDD_LRI, total_DDD_all_matrix              , DDD_per_treated_all)

LRI_DDD_df          <- make_df_2(Country_Year_meta, DDD_LRI, total_DDD_LRI_matrix              , DDD_per_treated_LRI         )
Fever_DDD_df        <- make_df_2(Country_Year_meta, DDD_LRI, total_DDD_Fever_matrix            , DDD_per_treated_Fever       )
all_DDD_df          <- make_df_2(Country_Year_meta, DDD_LRI, total_DDD_all_matrix              , DDD_per_treated_all)

years_to_save <- 2003:2019

DDD_summary <- data.frame(
  'Year' = LRI_DDD_summary$Year,
  'Country' = LRI_DDD_summary$Country,
  'Income' = LRI_DDD_summary$Income,
  'DDD_per_treated_LRI' = LRI_DDD_summary$Per_treated,
  'DDD_per_treated_Fever_without_LRI' = Fever_DDD_summary$Per_treated,
  'DDD_per_treated_LRI_or_Fever' = all_DDD_summary$Per_treated 
)

dir.create(
  "outputs/table/by_country_year",
  recursive = TRUE,
  showWarnings = FALSE
)

#### 2.1 CSV file ####
for (yr in years_to_save) {
  write.csv(DDD_summary[DDD_summary$Year == yr, ], paste0("outputs/table/by_country_year/DDD_per_treated_", yr, ".csv"), row.names = FALSE)
}

dir.create(
  "outputs/global_map_df",
  recursive = TRUE,
  showWarnings = FALSE
)

#### 2.1 CSV file ####
for (yr in years_to_save) {
  write.csv(LRI_DDD_df[LRI_DDD_df$Year == yr, ], paste0("outputs/global_map_df/2.1_DDD_LRI_", yr, ".csv"), row.names = FALSE)
  write.csv(Fever_DDD_df[Fever_DDD_df$Year == yr, ], paste0("outputs/global_map_df/2.1_DDD_Fever_", yr, ".csv"), row.names = FALSE)
  write.csv(all_DDD_df[all_DDD_df$Year == yr, ], paste0("outputs/global_map_df/2.1_DDD_all_", yr, ".csv"), row.names = FALSE)
}





#----------------------------2.2 by year (DDD) ------------------------------

years <- 2003:2019
summary_list_DDD_all <- summary_list_DDD_Fever <- summary_list_DDD_LRI <- list()

# ----- LRI------
for (yr in years) {
  # Filter matrices by year
  total_DDD_LRI_bind <- cbind(Country_Year_meta, total_DDD_LRI_matrix) %>%
    filter(Year == yr) %>%
    dplyr::select(-Country, -Year)
  
  treated_LRI_bind <- cbind(Country_Year_meta, treated_LRI_matrix) %>%
    filter(Year == yr) %>%
    dplyr::select(-Country, -Year)

  # Sum over countries for each simulation
  total_DDD_LRI <- colSums(total_DDD_LRI_bind)
  DDD_per_treated_LRI <- colSums(total_DDD_LRI_bind)/colSums(treated_LRI_bind)

  # Store results as long-format table
  summary_list_DDD_LRI[[as.character(yr)]] <- data.frame(
    Year = rep(yr, 6),
    Condition = rep(c("Total_DDD_LRI", "DDD_per_treated_LRI"), each = 3),
    Quantile = rep(c("Median", "2.5%", "97.5%"), times = 2),
    Value = c(
      quantile(total_DDD_LRI, probs = 0.50, na.rm = TRUE),
      quantile(total_DDD_LRI, probs = 0.025, na.rm = TRUE),
      quantile(total_DDD_LRI, probs = 0.975, na.rm = TRUE),
      
      quantile(DDD_per_treated_LRI, probs = 0.50, na.rm = TRUE),
      quantile(DDD_per_treated_LRI, probs = 0.025, na.rm = TRUE),
      quantile(DDD_per_treated_LRI, probs = 0.975, na.rm = TRUE)
    )
  )
}

# Combine all sublists
summary_all_years_DDD_LRI <- do.call(rbind, summary_list_DDD_LRI)

# Combine all years into one data frame
summary_wide_DDD_LRI <- summary_all_years_DDD_LRI %>%
  mutate(Column = paste0(gsub(" ", "_", Condition), "_", Quantile)) %>%
  dplyr::select(Year, Column, Value) %>%
  pivot_wider(names_from = Column, values_from = Value)

summary_wide_DDD_LRI_combined <- summary_wide_DDD_LRI %>%
  mutate(
    Total_DDD_LRI = sprintf("%s (%s  %s)", 
                          formatC(Total_DDD_LRI_Median, format = "f", digits = 0, big.mark = ","),
                          formatC(`Total_DDD_LRI_2.5%`, format = "f", digits = 0, big.mark = ","),
                          formatC(`Total_DDD_LRI_97.5%`, format = "f", digits = 0, big.mark = ",")),
    
    DDD_per_treated_LRI = sprintf("%s (%s  %s)", 
                            formatC(DDD_per_treated_LRI_Median, format = "f", digits = 2, big.mark = ","),
                            formatC(`DDD_per_treated_LRI_2.5%`, format = "f", digits = 2, big.mark = ","),
                            formatC(`DDD_per_treated_LRI_97.5%`, format = "f", digits = 2, big.mark = ","))
  ) %>%
  dplyr::select(Year, Total_DDD_LRI, DDD_per_treated_LRI)

LRI_year_DDD_per_treated <- summary_wide_DDD_LRI %>%
  dplyr::select(Year, DDD_per_treated_LRI_Median, `DDD_per_treated_LRI_2.5%`, `DDD_per_treated_LRI_97.5%`) %>%
  dplyr::rename(DDD_per_treated_LRI = DDD_per_treated_LRI_Median, DDD_per_treated_LRI_lower = `DDD_per_treated_LRI_2.5%`, DDD_per_treated_LRI_upper = `DDD_per_treated_LRI_97.5%`)



# ----- Fever------
for (yr in years) {
  # Filter matrices by year
  total_DDD_Fever_bind <- cbind(Country_Year_meta, total_DDD_Fever_matrix) %>%
    filter(Year == yr) %>%
    dplyr::select(-Country, -Year)
  
  treated_Fever_bind <- cbind(Country_Year_meta, treated_Fever_matrix) %>%
    filter(Year == yr) %>%
    dplyr::select(-Country, -Year)
  
  # Sum over countries for each simulation
  total_DDD_Fever <- colSums(total_DDD_Fever_bind)
  DDD_per_treated_Fever <- colSums(total_DDD_Fever_bind)/colSums(treated_Fever_bind)
  
  # Store results as long-format table
  summary_list_DDD_Fever[[as.character(yr)]] <- data.frame(
    Year = rep(yr, 6),
    Condition = rep(c("Total_DDD_Fever", "DDD_per_treated_Fever"), each = 3),
    Quantile = rep(c("Median", "2.5%", "97.5%"), times = 2),
    Value = c(
      quantile(total_DDD_Fever, probs = 0.50, na.rm = TRUE),
      quantile(total_DDD_Fever, probs = 0.025, na.rm = TRUE),
      quantile(total_DDD_Fever, probs = 0.975, na.rm = TRUE),
      
      quantile(DDD_per_treated_Fever, probs = 0.50, na.rm = TRUE),
      quantile(DDD_per_treated_Fever, probs = 0.025, na.rm = TRUE),
      quantile(DDD_per_treated_Fever, probs = 0.975, na.rm = TRUE)
    )
  )
}

# Combine all sublists
summary_all_years_DDD_Fever <- do.call(rbind, summary_list_DDD_Fever)

# Combine all years into one data frame
summary_wide_DDD_Fever <- summary_all_years_DDD_Fever %>%
  mutate(Column = paste0(gsub(" ", "_", Condition), "_", Quantile)) %>%
  dplyr::select(Year, Column, Value) %>%
  pivot_wider(names_from = Column, values_from = Value)

summary_wide_DDD_Fever_combined <- summary_wide_DDD_Fever %>%
  mutate(
    Total_DDD_Fever = sprintf("%s (%s  %s)", 
                            formatC(Total_DDD_Fever_Median, format = "f", digits = 0, big.mark = ","),
                            formatC(`Total_DDD_Fever_2.5%`, format = "f", digits = 0, big.mark = ","),
                            formatC(`Total_DDD_Fever_97.5%`, format = "f", digits = 0, big.mark = ",")),
    
    DDD_per_treated_Fever = sprintf("%s (%s  %s)", 
                                  formatC(DDD_per_treated_Fever_Median, format = "f", digits = 2, big.mark = ","),
                                  formatC(`DDD_per_treated_Fever_2.5%`, format = "f", digits = 2, big.mark = ","),
                                  formatC(`DDD_per_treated_Fever_97.5%`, format = "f", digits = 2, big.mark = ","))
  ) %>%
  dplyr::select(Year, Total_DDD_Fever, DDD_per_treated_Fever)

Fever_year_DDD_per_treated <- summary_wide_DDD_Fever %>%
  dplyr::select(Year, DDD_per_treated_Fever_Median, `DDD_per_treated_Fever_2.5%`, `DDD_per_treated_Fever_97.5%`) %>%
  dplyr::rename(DDD_per_treated_Fever = DDD_per_treated_Fever_Median, DDD_per_treated_Fever_lower = `DDD_per_treated_Fever_2.5%`, DDD_per_treated_Fever_upper = `DDD_per_treated_Fever_97.5%`)



# ----- all------
for (yr in years) {
  # Filter matrices by year
  total_DDD_all_bind <- cbind(Country_Year_meta, total_DDD_all_matrix) %>%
    filter(Year == yr) %>%
    dplyr::select(-Country, -Year)
  
  treated_all_bind <- cbind(Country_Year_meta, treated_all_matrix) %>%
    filter(Year == yr) %>%
    dplyr::select(-Country, -Year)
  
  # Sum over countries for each simulation
  total_DDD_all <- colSums(total_DDD_all_bind)
  DDD_per_treated_all <- colSums(total_DDD_all_bind) / colSums(treated_all_bind)
  
  # Store results as long-format table
  summary_list_DDD_all[[as.character(yr)]] <- data.frame(
    Year = rep(yr, 6),
    Condition = rep(c("Total_DDD_all", "DDD_per_treated_all"), each = 3),
    Quantile = rep(c("Median", "2.5%", "97.5%"), times = 2),
    Value = c(
      quantile(total_DDD_all, probs = 0.50, na.rm = TRUE),
      quantile(total_DDD_all, probs = 0.025, na.rm = TRUE),
      quantile(total_DDD_all, probs = 0.975, na.rm = TRUE),
      
      quantile(DDD_per_treated_all, probs = 0.50, na.rm = TRUE),
      quantile(DDD_per_treated_all, probs = 0.025, na.rm = TRUE),
      quantile(DDD_per_treated_all, probs = 0.975, na.rm = TRUE)
    )
  )
}

# Combine all sublists
summary_all_years_DDD_all <- do.call(rbind, summary_list_DDD_all)

# Combine all years into one data frame
summary_wide_DDD_all <- summary_all_years_DDD_all %>%
  mutate(Column = paste0(gsub(" ", "_", Condition), "_", Quantile)) %>%
  dplyr::select(Year, Column, Value) %>%
  pivot_wider(names_from = Column, values_from = Value)

summary_wide_DDD_all_combined <- summary_wide_DDD_all %>%
  mutate(
    Total_DDD_all = sprintf("%s (%s  %s)", 
                            formatC(Total_DDD_all_Median, format = "f", digits = 0, big.mark = ","),
                            formatC(`Total_DDD_all_2.5%`, format = "f", digits = 0, big.mark = ","),
                            formatC(`Total_DDD_all_97.5%`, format = "f", digits = 0, big.mark = ",")),
    
    DDD_per_treated_all = sprintf("%s (%s  %s)", 
                                  formatC(DDD_per_treated_all_Median, format = "f", digits = 2, big.mark = ","),
                                  formatC(`DDD_per_treated_all_2.5%`, format = "f", digits = 2, big.mark = ","),
                                  formatC(`DDD_per_treated_all_97.5%`, format = "f", digits = 2, big.mark = ","))
  ) %>%
  dplyr::select(Year, Total_DDD_all, DDD_per_treated_all)

all_year_DDD_per_treated <- summary_wide_DDD_all %>%
  dplyr::select(Year, DDD_per_treated_all_Median, `DDD_per_treated_all_2.5%`, `DDD_per_treated_all_97.5%`) %>%
  dplyr::rename(DDD_per_treated_all = DDD_per_treated_all_Median, DDD_per_treated_all_lower = `DDD_per_treated_all_2.5%`, DDD_per_treated_all_upper = `DDD_per_treated_all_97.5%`)

DDD_by_year <- data.frame(
  'Year' = summary_wide_DDD_LRI_combined$Year,
  'DDD_per_treated_LRI' = summary_wide_DDD_LRI_combined$DDD_per_treated_LRI, 
  'DDD_per_treated_Fever_without_LRI' = summary_wide_DDD_Fever_combined$DDD_per_treated_Fever,
  'DDD_per_treated_LRI_or_Fever' = summary_wide_DDD_all_combined$DDD_per_treated_all 
)

dir.create(
  "outputs/table/by_year/",
  recursive = TRUE,
  showWarnings = FALSE
)

#### 2.2 CSV file ####
write.csv(DDD_by_year, "outputs/table/by_year/DDD_per_treated_by_year.csv", row.names = FALSE)





#----------------------------2.2 by year (DDD increase) ------------------------------
total_DDD_LRI_2003 <- cbind(Country_Year_meta, total_DDD_LRI_matrix) %>%
  filter(Year == 2003) %>%
  dplyr::select(-Country, -Year)

incidence_LRI_2003 <- cbind(Country_Year_meta, Incidence_LRI_matrix) %>%
  filter(Year == 2003) %>%
  dplyr::select(-Country, -Year)

DDD_per_LRI_2003 <- colSums(total_DDD_LRI_2003)/colSums(incidence_LRI_2003)

total_DDD_LRI_2019 <- cbind(Country_Year_meta, total_DDD_LRI_matrix) %>%
  filter(Year == 2019) %>%
  dplyr::select(-Country, -Year)

incidence_LRI_2019 <- cbind(Country_Year_meta, Incidence_LRI_matrix) %>%
  filter(Year == 2019) %>%
  dplyr::select(-Country, -Year)

DDD_per_LRI_2019 <- colSums(total_DDD_LRI_2019)/colSums(incidence_LRI_2019)

ratio_LRI <- DDD_per_LRI_2019/DDD_per_LRI_2003

round(quantile(ratio_LRI, probs = 0.50, na.rm = TRUE), digits = 3)
round(quantile(ratio_LRI, probs = 0.025, na.rm = TRUE), digits = 3)
round(quantile(ratio_LRI, probs = 0.975, na.rm = TRUE), digits = 3)

total_DDD_Fever_2003 <- cbind(Country_Year_meta, total_DDD_Fever_matrix) %>%
  filter(Year == 2003) %>%
  dplyr::select(-Country, -Year)

incidence_Fever_2003 <- cbind(Country_Year_meta, Incidence_Fever_matrix) %>%
  filter(Year == 2003) %>%
  dplyr::select(-Country, -Year)

DDD_per_Fever_2003 <- colSums(total_DDD_Fever_2003)/colSums(incidence_Fever_2003)

total_DDD_Fever_2019 <- cbind(Country_Year_meta, total_DDD_Fever_matrix) %>%
  filter(Year == 2019) %>%
  dplyr::select(-Country, -Year)

incidence_Fever_2019 <- cbind(Country_Year_meta, Incidence_Fever_matrix) %>%
  filter(Year == 2019) %>%
  dplyr::select(-Country, -Year)

DDD_per_Fever_2019 <- colSums(total_DDD_Fever_2019)/colSums(incidence_Fever_2019)

ratio_Fever <- DDD_per_Fever_2019/DDD_per_Fever_2003

round(quantile(ratio_Fever, probs = 0.50, na.rm = TRUE), digits = 3)
round(quantile(ratio_Fever, probs = 0.025, na.rm = TRUE), digits = 3)
round(quantile(ratio_Fever, probs = 0.975, na.rm = TRUE), digits = 3)

   



#----------------------------2.3 by Income and year (treatment) ----------------------------

# -------LRI---------
Country_Year_Income_meta <- Country_Year_meta %>% inner_join(LRI_DDD_matrix[,c(1,2,4)], by = c("Country", "Year"))
Income_level <- c('L','LM','UM')
years <- 2003:2019
summary_list_DDD_LRI <- list()



for (inc in Income_level) {
  for (yr in years) {
    # Filter and select simulation matrices
    total_DDD_LRI_bind <- cbind(Country_Year_Income_meta, total_DDD_LRI_matrix) %>%
      filter(Income == inc, Year == yr) %>%
      dplyr::select(-Country, -Year, -Income)
    
    treated_LRI_bind <- cbind(Country_Year_Income_meta, treated_LRI_matrix) %>%
      filter(Income == inc, Year == yr) %>%
      dplyr::select(-Country, -Year, -Income)
    
    # Sum over countries for each simulation
    total_DDD_LRI <- colSums(total_DDD_LRI_bind)
    DDD_per_treated_LRI <- colSums(total_DDD_LRI_bind)/colSums(treated_LRI_bind)
    
    # Store results as long-format table
    summary_list_DDD_LRI[[paste(inc, yr, sep = "_")]] <- data.frame(
      Income = inc,
      Year = rep(yr, 6),
      Condition = rep(c("Total_DDD_LRI", "DDD_per_treated_LRI"), each = 3),
      Quantile = rep(c("Median", "2.5%", "97.5%"), times = 2),
      Value = c(
        quantile(total_DDD_LRI, probs = 0.50, na.rm = TRUE),
        quantile(total_DDD_LRI, probs = 0.025, na.rm = TRUE),
        quantile(total_DDD_LRI, probs = 0.975, na.rm = TRUE),
        
        quantile(DDD_per_treated_LRI, probs = 0.50, na.rm = TRUE),
        quantile(DDD_per_treated_LRI, probs = 0.025, na.rm = TRUE),
        quantile(DDD_per_treated_LRI, probs = 0.975, na.rm = TRUE)
      )
    )
  }
}

# Combine all sublists
summary_all_years_DDD_LRI <- do.call(rbind, summary_list_DDD_LRI)

# Combine all years into one data frame
summary_wide_DDD_LRI <- summary_all_years_DDD_LRI %>%
  mutate(Column = paste0(gsub(" ", "_", Condition), "_", Quantile)) %>%
  dplyr::select(Year, Income, Column, Value) %>%
  pivot_wider(names_from = Column, values_from = Value)


summary_wide_DDD_LRI_combined <- summary_wide_DDD_LRI %>%
  mutate(
    Total_DDD_LRI = sprintf("%s (%s  %s)", 
                            formatC(Total_DDD_LRI_Median, format = "f", digits = 0, big.mark = ","),
                            formatC(`Total_DDD_LRI_2.5%`, format = "f", digits = 0, big.mark = ","),
                            formatC(`Total_DDD_LRI_97.5%`, format = "f", digits = 0, big.mark = ",")),
    
    DDD_per_treated_LRI = sprintf("%s (%s  %s)", 
                                  formatC(DDD_per_treated_LRI_Median, format = "f", digits = 2, big.mark = ","),
                                  formatC(`DDD_per_treated_LRI_2.5%`, format = "f", digits = 2, big.mark = ","),
                                  formatC(`DDD_per_treated_LRI_97.5%`, format = "f", digits = 2, big.mark = ","))
  ) %>%
  dplyr::select(Year, Income, Total_DDD_LRI, DDD_per_treated_LRI)

LRI_year_DDD_per_treated_L <- summary_wide_DDD_LRI %>%
  filter(Income == 'L') %>%
  dplyr::select(Year, Income, DDD_per_treated_LRI_Median, `DDD_per_treated_LRI_2.5%`, `DDD_per_treated_LRI_97.5%`) %>%
  dplyr::rename(DDD_per_treated_Median = DDD_per_treated_LRI_Median, DDD_per_treated_LRI_Lower = `DDD_per_treated_LRI_2.5%`, DDD_per_treated_LRI_Upper = `DDD_per_treated_LRI_97.5%`)

LRI_year_DDD_per_treated_LM <- summary_wide_DDD_LRI %>%
  filter(Income == 'LM') %>%
  dplyr::select(Year, Income, DDD_per_treated_LRI_Median, `DDD_per_treated_LRI_2.5%`, `DDD_per_treated_LRI_97.5%`) %>%
  dplyr::rename(DDD_per_treated_Median = DDD_per_treated_LRI_Median, DDD_per_treated_LRI_Lower = `DDD_per_treated_LRI_2.5%`, DDD_per_treated_LRI_Upper = `DDD_per_treated_LRI_97.5%`)

LRI_year_DDD_per_treated_UM <- summary_wide_DDD_LRI %>%
  filter(Income == 'UM') %>%
  dplyr::select(Year, Income, DDD_per_treated_LRI_Median, `DDD_per_treated_LRI_2.5%`, `DDD_per_treated_LRI_97.5%`) %>%
  dplyr::rename(DDD_per_treated_Median = DDD_per_treated_LRI_Median, DDD_per_treated_LRI_Lower = `DDD_per_treated_LRI_2.5%`, DDD_per_treated_LRI_Upper = `DDD_per_treated_LRI_97.5%`)



# -------Fever---------
Country_Year_Income_meta <- Country_Year_meta %>% inner_join(Fever_DDD_matrix[,c(1,2,4)], by = c("Country", "Year"))
Income_level <- c('L','LM','UM')
years <- 2003:2019
summary_list_DDD_Fever <- list()

for (inc in Income_level) {
  for (yr in years) {
    # Filter and select simulation matrices
    total_DDD_Fever_bind <- cbind(Country_Year_Income_meta, total_DDD_Fever_matrix) %>%
      filter(Income == inc, Year == yr) %>%
      dplyr::select(-Country, -Year, -Income)
    
    treated_Fever_bind <- cbind(Country_Year_Income_meta, treated_Fever_matrix) %>%
      filter(Income == inc, Year == yr) %>%
      dplyr::select(-Country, -Year, -Income)
    
    # Sum over countries for each simulation
    total_DDD_Fever <- colSums(total_DDD_Fever_bind)
    DDD_per_treated_Fever <- colSums(total_DDD_Fever_bind)/colSums(treated_Fever_bind)
    
    # Store results as long-format table
    summary_list_DDD_Fever[[paste(inc, yr, sep = "_")]] <- data.frame(
      Income = inc,
      Year = rep(yr, 6),
      Condition = rep(c("Total_DDD_Fever", "DDD_per_treated_Fever"), each = 3),
      Quantile = rep(c("Median", "2.5%", "97.5%"), times = 2),
      Value = c(
        quantile(total_DDD_Fever, probs = 0.50, na.rm = TRUE),
        quantile(total_DDD_Fever, probs = 0.025, na.rm = TRUE),
        quantile(total_DDD_Fever, probs = 0.975, na.rm = TRUE),
        
        quantile(DDD_per_treated_Fever, probs = 0.50, na.rm = TRUE),
        quantile(DDD_per_treated_Fever, probs = 0.025, na.rm = TRUE),
        quantile(DDD_per_treated_Fever, probs = 0.975, na.rm = TRUE)
      )
    )
  }
}

# Combine all sublists
summary_all_years_DDD_Fever <- do.call(rbind, summary_list_DDD_Fever)

# Combine all years into one data frame
summary_wide_DDD_Fever <- summary_all_years_DDD_Fever %>%
  mutate(Column = paste0(gsub(" ", "_", Condition), "_", Quantile)) %>%
  dplyr::select(Year, Income, Column, Value) %>%
  pivot_wider(names_from = Column, values_from = Value)


summary_wide_DDD_Fever_combined <- summary_wide_DDD_Fever %>%
  mutate(
    Total_DDD_Fever = sprintf("%s (%s  %s)", 
                            formatC(Total_DDD_Fever_Median, format = "f", digits = 0, big.mark = ","),
                            formatC(`Total_DDD_Fever_2.5%`, format = "f", digits = 0, big.mark = ","),
                            formatC(`Total_DDD_Fever_97.5%`, format = "f", digits = 0, big.mark = ",")),
    
    DDD_per_treated_Fever = sprintf("%s (%s  %s)", 
                                  formatC(DDD_per_treated_Fever_Median, format = "f", digits = 2, big.mark = ","),
                                  formatC(`DDD_per_treated_Fever_2.5%`, format = "f", digits = 2, big.mark = ","),
                                  formatC(`DDD_per_treated_Fever_97.5%`, format = "f", digits = 2, big.mark = ","))
  ) %>%
  dplyr::select(Year, Income, Total_DDD_Fever, DDD_per_treated_Fever)



Fever_year_DDD_per_treated_L <- summary_wide_DDD_Fever %>%
  filter(Income == 'L') %>%
  dplyr::select(Year, Income, DDD_per_treated_Fever_Median, `DDD_per_treated_Fever_2.5%`, `DDD_per_treated_Fever_97.5%`) %>%
  dplyr::rename(DDD_per_treated_Median = DDD_per_treated_Fever_Median, DDD_per_treated_Fever_Lower = `DDD_per_treated_Fever_2.5%`, DDD_per_treated_Fever_Upper = `DDD_per_treated_Fever_97.5%`)

Fever_year_DDD_per_treated_LM <- summary_wide_DDD_Fever %>%
  filter(Income == 'LM') %>%
  dplyr::select(Year, Income, DDD_per_treated_Fever_Median, `DDD_per_treated_Fever_2.5%`, `DDD_per_treated_Fever_97.5%`) %>%
  dplyr::rename(DDD_per_treated_Median = DDD_per_treated_Fever_Median, DDD_per_treated_Fever_Lower = `DDD_per_treated_Fever_2.5%`, DDD_per_treated_Fever_Upper = `DDD_per_treated_Fever_97.5%`)

Fever_year_DDD_per_treated_UM <- summary_wide_DDD_Fever %>%
  filter(Income == 'UM') %>%
  dplyr::select(Year, Income, DDD_per_treated_Fever_Median, `DDD_per_treated_Fever_2.5%`, `DDD_per_treated_Fever_97.5%`) %>%
  dplyr::rename(DDD_per_treated_Median = DDD_per_treated_Fever_Median, DDD_per_treated_Fever_Lower = `DDD_per_treated_Fever_2.5%`, DDD_per_treated_Fever_Upper = `DDD_per_treated_Fever_97.5%`)



# -------all---------
Country_Year_Income_meta <- Country_Year_meta %>% inner_join(all_DDD_matrix[,c(1,2,4)], by = c("Country", "Year"))
Income_level <- c('L','LM','UM')
years <- 2003:2019
summary_list_DDD_all <- list()

for (inc in Income_level) {
  for (yr in years) {
    # Filter and select simulation matrices
    total_DDD_all_bind <- cbind(Country_Year_Income_meta, total_DDD_all_matrix) %>%
      filter(Income == inc, Year == yr) %>%
      dplyr::select(-Country, -Year, -Income)
    
    treated_all_bind <- cbind(Country_Year_Income_meta, treated_all_matrix) %>%
      filter(Income == inc, Year == yr) %>%
      dplyr::select(-Country, -Year, -Income)
    
    # Sum over countries for each simulation
    total_DDD_all <- colSums(total_DDD_all_bind)
    DDD_per_treated_all <- colSums(total_DDD_all_bind)/colSums(treated_all_bind)
    
    # Store results as long-format table
    summary_list_DDD_all[[paste(inc, yr, sep = "_")]] <- data.frame(
      Income = inc,
      Year = rep(yr, 6),
      Condition = rep(c("Total_DDD_all", "DDD_per_treated_all"), each = 3),
      Quantile = rep(c("Median", "2.5%", "97.5%"), times = 2),
      Value = c(
        quantile(total_DDD_all, probs = 0.50, na.rm = TRUE),
        quantile(total_DDD_all, probs = 0.025, na.rm = TRUE),
        quantile(total_DDD_all, probs = 0.975, na.rm = TRUE),
        
        quantile(DDD_per_treated_all, probs = 0.50, na.rm = TRUE),
        quantile(DDD_per_treated_all, probs = 0.025, na.rm = TRUE),
        quantile(DDD_per_treated_all, probs = 0.975, na.rm = TRUE)
      )
    )
  }
}

# Combine all sublists
summary_all_years_DDD_all <- do.call(rbind, summary_list_DDD_all)

# Combine all years into one data frame
summary_wide_DDD_all <- summary_all_years_DDD_all %>%
  mutate(Column = paste0(gsub(" ", "_", Condition), "_", Quantile)) %>%
  dplyr::select(Year, Income, Column, Value) %>%
  pivot_wider(names_from = Column, values_from = Value)


summary_wide_DDD_all_combined <- summary_wide_DDD_all %>%
  mutate(
    Total_DDD_all = sprintf("%s (%s  %s)", 
                            formatC(Total_DDD_all_Median, format = "f", digits = 0, big.mark = ","),
                            formatC(`Total_DDD_all_2.5%`, format = "f", digits = 0, big.mark = ","),
                            formatC(`Total_DDD_all_97.5%`, format = "f", digits = 0, big.mark = ",")),
    
    DDD_per_treated_all = sprintf("%s (%s  %s)", 
                                  formatC(DDD_per_treated_all_Median, format = "f", digits = 2, big.mark = ","),
                                  formatC(`DDD_per_treated_all_2.5%`, format = "f", digits = 2, big.mark = ","),
                                  formatC(`DDD_per_treated_all_97.5%`, format = "f", digits = 2, big.mark = ","))
  ) %>%
  dplyr::select(Year, Income, Total_DDD_all, DDD_per_treated_all)

all_year_DDD_per_treated_L <- summary_wide_DDD_all %>%
  filter(Income == 'L') %>%
  dplyr::select(Year, Income, DDD_per_treated_all_Median, `DDD_per_treated_all_2.5%`, `DDD_per_treated_all_97.5%`) %>%
  dplyr::rename(DDD_per_treated_Median = DDD_per_treated_all_Median, DDD_per_treated_all_Lower = `DDD_per_treated_all_2.5%`, DDD_per_treated_all_Upper = `DDD_per_treated_all_97.5%`)

all_year_DDD_per_treated_LM <- summary_wide_DDD_all %>%
  filter(Income == 'LM') %>%
  dplyr::select(Year, Income, DDD_per_treated_all_Median, `DDD_per_treated_all_2.5%`, `DDD_per_treated_all_97.5%`) %>%
  dplyr::rename(DDD_per_treated_Median = DDD_per_treated_all_Median, DDD_per_treated_all_Lower = `DDD_per_treated_all_2.5%`, DDD_per_treated_all_Upper = `DDD_per_treated_all_97.5%`)

all_year_DDD_per_treated_UM <- summary_wide_DDD_all %>%
  filter(Income == 'UM') %>%
  dplyr::select(Year, Income, DDD_per_treated_all_Median, `DDD_per_treated_all_2.5%`, `DDD_per_treated_all_97.5%`) %>%
  dplyr::rename(DDD_per_treated_Median = DDD_per_treated_all_Median, DDD_per_treated_all_Lower = `DDD_per_treated_all_2.5%`, DDD_per_treated_all_Upper = `DDD_per_treated_all_97.5%`)

DDD_by_year_income <- data.frame(
  'Year' = summary_wide_DDD_LRI_combined$Year,
  'Income' = summary_wide_DDD_LRI_combined$Income,
  'DDD_per_treated_LRI' = summary_wide_DDD_LRI_combined$DDD_per_treated_LRI, 
  'DDD_per_treated_Fever_without_LRI' = summary_wide_DDD_Fever_combined$DDD_per_treated_Fever,
  'DDD_per_treated_LRI_or_Fever' = summary_wide_DDD_all_combined$DDD_per_treated_all 
)

dir.create(
  "outputs/table/by_year_income/",
  recursive = TRUE,
  showWarnings = FALSE
)

#### 2.2 CSV file ####
write.csv(DDD_by_year_income, "outputs/table/by_year_income/DDD_per_treated_by_year_income.csv", row.names = FALSE)


# save the plot data
write.csv(LRI_year_DDD_per_treated     , "outputs/plot_df/2.3_plot_Ab_LRI_year_DDD_per_treated.csv")
write.csv(LRI_year_DDD_per_treated_L   , "outputs/plot_df/2.3_plot_Ab_LRI_year_DDD_per_treated_L.csv")
write.csv(LRI_year_DDD_per_treated_LM  , "outputs/plot_df/2.3_plot_Ab_LRI_year_DDD_per_treated_LM.csv")
write.csv(LRI_year_DDD_per_treated_UM  , "outputs/plot_df/2.3_plot_Ab_LRI_year_DDD_per_treated_UM.csv")
write.csv(Fever_year_DDD_per_treated   , "outputs/plot_df/2.3_plot_Ab_Fever_year_DDD_per_treated.csv")
write.csv(Fever_year_DDD_per_treated_L , "outputs/plot_df/2.3_plot_Ab_Fever_year_DDD_per_treated_L.csv")
write.csv(Fever_year_DDD_per_treated_LM, "outputs/plot_df/2.3_plot_Ab_Fever_year_DDD_per_treated_LM.csv")
write.csv(Fever_year_DDD_per_treated_UM, "outputs/plot_df/2.3_plot_Ab_Fever_year_DDD_per_treated_UM.csv")
write.csv(all_year_DDD_per_treated     , "outputs/plot_df/2.3_plot_Ab_all_year_DDD_per_treated.csv")
write.csv(all_year_DDD_per_treated_L   , "outputs/plot_df/2.3_plot_Ab_all_year_DDD_per_treated_L.csv")
write.csv(all_year_DDD_per_treated_LM  , "outputs/plot_df/2.3_plot_Ab_all_year_DDD_per_treated_LM.csv")
write.csv(all_year_DDD_per_treated_UM  , "outputs/plot_df/2.3_plot_Ab_all_year_DDD_per_treated_UM.csv")

#