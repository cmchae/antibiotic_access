##### 0 #####
rm(list = ls())
gc()

requiredPackages = c("dplyr", "tidyr", "tibble", "data.table", "stringr", "MASS")
lib.path <- .libPaths()
for(p in requiredPackages){
  if(!require(p,character.only = TRUE)) install.packages(p, lib = lib.path)
  library(p,character.only = TRUE)
}



##### 1 #####
rm(list = ls())
gc()
n_sim <- 10000; set.seed(123)

#---------------------------prepare data for figure------------------------
Ab_num <- 4
UNI_indexed <- readRDS("outputs/UNI_indexed.rds")
UNI_indexed_Abaccess <- UNI_indexed %>% filter(Abaccess == 1)
UNI_indexed_Abdata <- UNI_indexed %>% filter(ab.data == 1)
UNI_indexed_both <- UNI_indexed %>% filter(Abaccess == 1, ab.data == 1)

# outputs generated from the current PC1-PC2 models
p_best_LRI <- readRDS("outputs/plot_df/p_best/p_LRI4.rds")
LRI_DDD_ratio_sim <- readRDS("outputs/plot_df/p_best/ratio_LRI4.rds")
LRI_Ab_probability.matrix <- list()

for (i in 1:Ab_num) {
  LRI_ratio_single_Ab <- LRI_DDD_ratio_sim[,(n_sim*(i-1)+1):(n_sim*i)]
  LRI_Ab_probability.matrix[[i]] <- p_best_LRI*LRI_ratio_single_Ab
}

Ab.names <- c("BSP","NSP","Cephalosporins","Aminoglycosides")
names(LRI_Ab_probability.matrix) <- Ab.names

for (i in 1:Ab_num) {
  ab_name <- Ab.names[i]
  AbLRI_long <- as.vector(LRI_Ab_probability.matrix[[i]])
  assign(paste0("AbLRI_long_", ab_name), AbLRI_long)
}

# Create matching indices
Year_index <- rep(UNI_indexed$Year, times = n_sim)
Country_index <- rep(UNI_indexed$Country, times = n_sim)
GDP_index <- rep(UNI_indexed$GDP, times = n_sim)
Income_index <- rep(UNI_indexed$Income, times = n_sim)
Malaria_index <- rep(UNI_indexed$Endemic, times = n_sim)
GBD_index <- rep(UNI_indexed$gbd_region, times = n_sim)

# Combine into dataframe
sim_df_LRI <- tibble(
  Country = Country_index,
  Year = Year_index,
  GDP = GDP_index,
  Income = Income_index,
  Malaria = Malaria_index,
  GBD_region = GBD_index,
  AbLRI_BSP = AbLRI_long_BSP,
  AbLRI_NSP = AbLRI_long_NSP,
  AbLRI_Cephalosporins = AbLRI_long_Cephalosporins,
  AbLRI_Aminoglycosides = AbLRI_long_Aminoglycosides,
  sim = rep(1:n_sim, each = nrow(UNI_indexed))
)

rm(p_best_LRI, LRI_DDD_ratio_sim)
rm_list <- c(paste0("AbLRI_long_", Ab.names))
rm(list = rm_list, envir = .GlobalEnv)

Country_Year_Antibiotics <- tibble(
  Country = sim_df_LRI$Country,
  Year = sim_df_LRI$Year,
  sim = sim_df_LRI$sim,
  GDP = sim_df_LRI$GDP,
  Income = sim_df_LRI$Income,
  Malaria = sim_df_LRI$Malaria,
  GBD_region = sim_df_LRI$GBD_region
)

U5_matrix <- readRDS("outputs/matrix.result.ver.1.rds")

U5_matrix_meta <- U5_matrix[,1:2]
U5_matrix_LRI <- U5_matrix[,3:(n_sim+2)] + U5_matrix[,(2*n_sim+3):(3*n_sim+2)]
U5_matrix_LRI_index <- as.vector(as.matrix(U5_matrix_LRI))
U5_Year_index <- rep(U5_matrix_meta$Year, times = n_sim)
U5_Country_index <- rep(U5_matrix_meta$Country, times = n_sim)

sim_incidence <- tibble(
  Country = U5_Country_index,
  Year = U5_Year_index,
  Incidence_LRI = U5_matrix_LRI_index,
  sim = rep(1:n_sim, each = nrow(U5_matrix_meta))
)

Country_Year_Incidence <- tibble(
  Country = sim_incidence$Country,
  Year = sim_incidence$Year,
  sim = sim_incidence$sim
)

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
sim_df_LRI <- sim_df_LRI %>% mutate(Country = ifelse(Country %in% names(country_mapping), country_mapping[Country], Country))

Country_Year_meta <- inner_join(Country_Year_Incidence, Country_Year_Antibiotics, by = c("Country", "Year", "sim"))

rm(U5_matrix, U5_matrix_LRI, U5_matrix_meta, 
   UNI_indexed, UNI_indexed_Abaccess, UNI_indexed_Abdata, UNI_indexed_both,
   LRI_ratio_single_Ab, LRI_Ab_probability.matrix,
   U5_Country_index, U5_matrix_LRI_index, U5_Year_index, 
   Year_index, Malaria_index, Income_index, GBD_index, GDP_index, Country_index, AbLRI_long)
gc()

###### 2 ######
Antibiotics_LRI <- Country_Year_meta %>% inner_join(sim_df_LRI, by = c("Country","Year","sim","GDP","Income","Malaria","GBD_region"))
rm(sim_df_LRI)

colnames(sim_incidence)
colnames(Antibiotics_LRI)

rm(Country_Year_Incidence, Country_Year_Antibiotics)
gc()

###### 3 ########
treated_LRI_matrix <- tibble(matrix(NA_real_, nrow = nrow(sim_incidence), ncol = Ab_num))
for (i in 1:Ab_num) {
  treated_LRI_matrix[,i] <- Antibiotics_LRI[, 7 + i] * sim_incidence$Incidence_LRI
}

proportion_LRI_matrix <- tibble(matrix(NA_real_, nrow = nrow(sim_incidence), ncol = Ab_num))
for (i in 1:Ab_num) {
  proportion_LRI_matrix[,i] <- treated_LRI_matrix[,i] / sim_incidence$Incidence_LRI
}

colnames(treated_LRI_matrix) <- gsub("AbLRI_", "Treated_LRI_", colnames(Antibiotics_LRI[, 8:11]))
colnames(proportion_LRI_matrix) <- gsub("AbLRI_", "proportion_LRI_", colnames(Antibiotics_LRI[, 8:11]))

treated_LRI_matrix <- cbind(Country_Year_meta, treated_LRI_matrix)
proportion_LRI_matrix <- cbind(Country_Year_meta, proportion_LRI_matrix)

rm(Antibiotics_LRI)
gc()

saveRDS(Country_Year_meta,"outputs/Country_Year_meta4.rds")
saveRDS(sim_incidence,"outputs/sim_incidence4.rds")
saveRDS(treated_LRI_matrix,"outputs/treated_LRI_matrix4.rds")
saveRDS(proportion_LRI_matrix, "outputs/proportion_LRI_matrix4.rds")

rm(list = ls())
gc()



##### 0 #####
rm(list = ls())
gc()

requiredPackages = c("dplyr", "tidyr", "tibble", "data.table", "stringr", "MASS")
lib.path <- .libPaths()
for(p in requiredPackages){
  if(!require(p,character.only = TRUE)) install.packages(p, lib = lib.path)
  library(p,character.only = TRUE)
}

##### 1 #####
n_sim <- 10000; set.seed(123)
Country_Year_meta <- readRDS("outputs/Country_Year_meta4.rds")
sim_incidence <- readRDS("outputs/sim_incidence4.rds")
treated_LRI_matrix <- readRDS("outputs/treated_LRI_matrix4.rds")
proportion_LRI_matrix <- readRDS("outputs/proportion_LRI_matrix4.rds")
digit_number <- c(3,4,3,3)

#----------------------------1.1 by country and year (treatment) ------------------------------
### LRI ###
treated_LRI_result <- treated_LRI_matrix %>%
  group_by(Country, Year, GDP, Income, Malaria, GBD_region) %>%
  summarise(
    Treated_LRI_BSP_med = median(Treated_LRI_BSP, na.rm = TRUE),
    Treated_LRI_BSP_lower = quantile(Treated_LRI_BSP, 0.025, na.rm = TRUE),
    Treated_LRI_BSP_upper = quantile(Treated_LRI_BSP, 0.975, na.rm = TRUE),
    
    Treated_LRI_NSP_med = median(Treated_LRI_NSP, na.rm = TRUE),
    Treated_LRI_NSP_lower = quantile(Treated_LRI_NSP, 0.025, na.rm = TRUE),
    Treated_LRI_NSP_upper = quantile(Treated_LRI_NSP, 0.975, na.rm = TRUE),
    
    Treated_LRI_Cephalosporins_med = median(Treated_LRI_Cephalosporins, na.rm = TRUE),
    Treated_LRI_Cephalosporins_lower = quantile(Treated_LRI_Cephalosporins, 0.025, na.rm = TRUE),
    Treated_LRI_Cephalosporins_upper = quantile(Treated_LRI_Cephalosporins, 0.975, na.rm = TRUE),
   
    Treated_LRI_Aminoglycosides_med = median(Treated_LRI_Aminoglycosides, na.rm = TRUE),
    Treated_LRI_Aminoglycosides_lower = quantile(Treated_LRI_Aminoglycosides, 0.025, na.rm = TRUE),
    Treated_LRI_Aminoglycosides_upper = quantile(Treated_LRI_Aminoglycosides, 0.975, na.rm = TRUE),
    
    .groups = "drop"
  )


abx_names <- c(
  "Treated_LRI_BSP", "Treated_LRI_NSP", "Treated_LRI_Cephalosporins", "Treated_LRI_Aminoglycosides")

for (i in 1:length(abx_names)) {
  med_col <- paste0(abx_names[i], "_med")
  lwr_col <- paste0(abx_names[i], "_lower")
  upr_col <- paste0(abx_names[i], "_upper")
  fmt_col <- paste0(abx_names[i])
  
  treated_LRI_result[[fmt_col]] <- sprintf(
    "%s (%s, %s)",
    formatC(treated_LRI_result[[med_col]], format = "f", big.mark = " ", digits = 0),
    formatC(treated_LRI_result[[lwr_col]], format = "f", big.mark = " ", digits = 0),
    formatC(treated_LRI_result[[upr_col]], format = "f", big.mark = " ", digits = 0)
  )
}

treated_LRI_result_summary <- treated_LRI_result[, c(1:6, 19:22)]
treated_LRI_result_df <- treated_LRI_result[, c(1:18)]

rm(treated_LRI_result)
gc()


proportion_LRI_result <- proportion_LRI_matrix %>%
  group_by(Country, Year, GDP, Income, Malaria, GBD_region) %>%
  summarise(
    proportion_LRI_BSP_med = median(proportion_LRI_BSP, na.rm = TRUE),
    proportion_LRI_BSP_lower = quantile(proportion_LRI_BSP, 0.025, na.rm = TRUE),
    proportion_LRI_BSP_upper = quantile(proportion_LRI_BSP, 0.975, na.rm = TRUE),
    
    proportion_LRI_NSP_med = median(proportion_LRI_NSP, na.rm = TRUE),
    proportion_LRI_NSP_lower = quantile(proportion_LRI_NSP, 0.025, na.rm = TRUE),
    proportion_LRI_NSP_upper = quantile(proportion_LRI_NSP, 0.975, na.rm = TRUE),
    
    proportion_LRI_Cephalosporins_med = median(proportion_LRI_Cephalosporins, na.rm = TRUE),
    proportion_LRI_Cephalosporins_lower = quantile(proportion_LRI_Cephalosporins, 0.025, na.rm = TRUE),
    proportion_LRI_Cephalosporins_upper = quantile(proportion_LRI_Cephalosporins, 0.975, na.rm = TRUE),
    
    proportion_LRI_Aminoglycosides_med = median(proportion_LRI_Aminoglycosides, na.rm = TRUE),
    proportion_LRI_Aminoglycosides_lower = quantile(proportion_LRI_Aminoglycosides, 0.025, na.rm = TRUE),
    proportion_LRI_Aminoglycosides_upper = quantile(proportion_LRI_Aminoglycosides, 0.975, na.rm = TRUE),
    
    .groups = "drop"
  )


abx_names <- c(
  "proportion_LRI_BSP", "proportion_LRI_NSP", "proportion_LRI_Cephalosporins", "proportion_LRI_Aminoglycosides"
)

for (i in 1:length(abx_names)) {
  med_col <- paste0(abx_names[i], "_med")
  lwr_col <- paste0(abx_names[i], "_lower")
  upr_col <- paste0(abx_names[i], "_upper")
  fmt_col <- paste0(abx_names[i])
  
  proportion_LRI_result[[fmt_col]] <- sprintf(
    "%s (%s, %s)",
    formatC(proportion_LRI_result[[med_col]], format = "f", big.mark = " ", digit_number[i]),
    formatC(proportion_LRI_result[[lwr_col]], format = "f", big.mark = " ", digit_number[i]),
    formatC(proportion_LRI_result[[upr_col]], format = "f", big.mark = " ", digit_number[i])
  )
}

proportion_LRI_result_summary <- proportion_LRI_result[, c(1:6, 19:22)]
proportion_LRI_result_df <- proportion_LRI_result[, c(1:18)]
rm(proportion_LRI_result)
gc()

LRI_result_summary <- cbind(treated_LRI_result_summary, proportion_LRI_result_summary[7:10])
years_to_save <- 2003:2019


dir.create(
  "outputs/table/by_class_country_year/probability",
  recursive = TRUE,
  showWarnings = FALSE
)

#### 1.1 CSV file ####
for (yr in years_to_save) {
  write.csv(LRI_result_summary[LRI_result_summary$Year == yr, ], paste0("outputs/table/by_class_country_year/probability/LRI_", yr, "4.csv"), row.names = FALSE)
}

rm(LRI_result_summary, treated_LRI_result, treated_LRI_result_df, treated_LRI_result_summary, 
   proportion_LRI_matrix, proportion_LRI_result, proportion_LRI_result_df, proportion_LRI_result_summary)
gc()

treated_LRI_matrix <- cbind(treated_LRI_matrix, sim_incidence[,3])



#----------------------------1.2 by year (treatment) ------------------------------
#---LRI---
treated_LRI_sim <- treated_LRI_matrix %>% 
  group_by(Year, sim) %>%
  summarise(
    Treated_LRI_BSP = sum(Treated_LRI_BSP, na.rm = TRUE),
    Treated_LRI_NSP = sum(Treated_LRI_NSP, na.rm = TRUE),
    Treated_LRI_Cephalosporins = sum(Treated_LRI_Cephalosporins, na.rm = TRUE),
    Treated_LRI_Aminoglycosides = sum(Treated_LRI_Aminoglycosides, na.rm = TRUE),
    Incidence = sum(Incidence_LRI, na.rm = TRUE),
    .groups = "drop"
  )

treated_LRI_result <- treated_LRI_sim %>%
  group_by(Year) %>%
  summarise(
    
    proportion_LRI_BSP_med = median(Treated_LRI_BSP/Incidence, na.rm = TRUE),
    proportion_LRI_BSP_lower = quantile(Treated_LRI_BSP/Incidence, 0.025, na.rm = TRUE),
    proportion_LRI_BSP_upper = quantile(Treated_LRI_BSP/Incidence, 0.975, na.rm = TRUE),
    
    proportion_LRI_NSP_med = median(Treated_LRI_NSP/Incidence, na.rm = TRUE),
    proportion_LRI_NSP_lower = quantile(Treated_LRI_NSP/Incidence, 0.025, na.rm = TRUE),
    proportion_LRI_NSP_upper = quantile(Treated_LRI_NSP/Incidence, 0.975, na.rm = TRUE),
    
    proportion_LRI_Cephalosporins_med = median(Treated_LRI_Cephalosporins/Incidence, na.rm = TRUE),
    proportion_LRI_Cephalosporins_lower = quantile(Treated_LRI_Cephalosporins/Incidence, 0.025, na.rm = TRUE),
    proportion_LRI_Cephalosporins_upper = quantile(Treated_LRI_Cephalosporins/Incidence, 0.975, na.rm = TRUE),
    
    proportion_LRI_Aminoglycosides_med = median(Treated_LRI_Aminoglycosides/Incidence, na.rm = TRUE),
    proportion_LRI_Aminoglycosides_lower = quantile(Treated_LRI_Aminoglycosides/Incidence, 0.025, na.rm = TRUE),
    proportion_LRI_Aminoglycosides_upper = quantile(Treated_LRI_Aminoglycosides/Incidence, 0.975, na.rm = TRUE),

    .groups = "drop"
  )


abx_names <- c(
  "proportion_LRI_BSP", "proportion_LRI_NSP", "proportion_LRI_Cephalosporins","proportion_LRI_Aminoglycosides"
)

for (i in 1:length(abx_names)) {
  med_col <- paste0(abx_names[i], "_med")
  lwr_col <- paste0(abx_names[i], "_lower")
  upr_col <- paste0(abx_names[i], "_upper")
  fmt_col <- paste0(abx_names[i])
  
  treated_LRI_result[[fmt_col]] <- sprintf(
    "%s (%s, %s)",
    formatC(treated_LRI_result[[med_col]], format = "f", big.mark = " ", digits = digit_number[i]),
    formatC(treated_LRI_result[[lwr_col]], format = "f", big.mark = " ", digits = digit_number[i]),
    formatC(treated_LRI_result[[upr_col]], format = "f", big.mark = " ", digits = digit_number[i])
  )
}

treated_LRI_summary <- treated_LRI_result[,c(1,14:17)]
treated_LRI_df <- treated_LRI_result[,c(1:13)]

dir.create(
  "outputs/table/by_class_year/probability",
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  "outputs/by_class_df/by_year/probability",
  recursive = TRUE,
  showWarnings = FALSE
)

#### 1.2 CSV file ####
write.csv(treated_LRI_summary, "outputs/table/by_class_year/probability/LRI4.csv")
write.csv(treated_LRI_df, "outputs/by_class_df/by_year/probability/LRI4.csv")
rm(treated_LRI_sim, treated_LRI_result, treated_LRI_summary, treated_LRI_df)



#----------------------------1.3 by year and malaria (treatment) ------------------------------
#---LRI---
treated_LRI_sim <- treated_LRI_matrix %>% 
  group_by(Year, Malaria, sim) %>%
  summarise(
    Treated_LRI_BSP = sum(Treated_LRI_BSP, na.rm = TRUE),
    Treated_LRI_NSP = sum(Treated_LRI_NSP, na.rm = TRUE),
    Treated_LRI_Cephalosporins = sum(Treated_LRI_Cephalosporins, na.rm = TRUE),
    Treated_LRI_Aminoglycosides = sum(Treated_LRI_Aminoglycosides, na.rm = TRUE),
    Incidence = sum(Incidence_LRI, na.rm = TRUE),
    .groups = "drop"
  )

treated_LRI_result <- treated_LRI_sim %>%
  group_by(Year, Malaria) %>%
  summarise(
    
    proportion_LRI_BSP_med = median(Treated_LRI_BSP/Incidence, na.rm = TRUE),
    proportion_LRI_BSP_lower = quantile(Treated_LRI_BSP/Incidence, 0.025, na.rm = TRUE),
    proportion_LRI_BSP_upper = quantile(Treated_LRI_BSP/Incidence, 0.975, na.rm = TRUE),
    
    proportion_LRI_NSP_med = median(Treated_LRI_NSP/Incidence, na.rm = TRUE),
    proportion_LRI_NSP_lower = quantile(Treated_LRI_NSP/Incidence, 0.025, na.rm = TRUE),
    proportion_LRI_NSP_upper = quantile(Treated_LRI_NSP/Incidence, 0.975, na.rm = TRUE),
    
    proportion_LRI_Cephalosporins_med = median(Treated_LRI_Cephalosporins/Incidence, na.rm = TRUE),
    proportion_LRI_Cephalosporins_lower = quantile(Treated_LRI_Cephalosporins/Incidence, 0.025, na.rm = TRUE),
    proportion_LRI_Cephalosporins_upper = quantile(Treated_LRI_Cephalosporins/Incidence, 0.975, na.rm = TRUE),
    
    proportion_LRI_Aminoglycosides_med = median(Treated_LRI_Aminoglycosides/Incidence, na.rm = TRUE),
    proportion_LRI_Aminoglycosides_lower = quantile(Treated_LRI_Aminoglycosides/Incidence, 0.025, na.rm = TRUE),
    proportion_LRI_Aminoglycosides_upper = quantile(Treated_LRI_Aminoglycosides/Incidence, 0.975, na.rm = TRUE),
    
    .groups = "drop"
  )

abx_names <- c(
  "proportion_LRI_BSP", "proportion_LRI_NSP", "proportion_LRI_Cephalosporins", "proportion_LRI_Aminoglycosides")

for (i in 1:length(abx_names)) {
  med_col <- paste0(abx_names[i], "_med")
  lwr_col <- paste0(abx_names[i], "_lower")
  upr_col <- paste0(abx_names[i], "_upper")
  fmt_col <- paste0(abx_names[i])
  
  treated_LRI_result[[fmt_col]] <- sprintf(
    "%s (%s, %s)",
    formatC(treated_LRI_result[[med_col]], format = "f", big.mark = " ", digits = digit_number[i]),
    formatC(treated_LRI_result[[lwr_col]], format = "f", big.mark = " ", digits = digit_number[i]),
    formatC(treated_LRI_result[[upr_col]], format = "f", big.mark = " ", digits = digit_number[i])
  )
}

treated_LRI_summary <- treated_LRI_result[,c(1:2,15:18)]
treated_LRI_df <- treated_LRI_result[,c(1:14)]

dir.create(
  "outputs/table/by_class_malaria_year/probability",
  recursive = TRUE,
  showWarnings = FALSE
)


#### 1.3 CSV file ####
write.csv(treated_LRI_summary, "outputs/table/by_class_malaria_year/probability/LRI4.csv")

Malaria_name <- unique(treated_LRI_df$Malaria)

for (mal in Malaria_name) {
  write.csv(treated_LRI_summary[treated_LRI_summary$Malaria == mal, ], paste0("outputs/table/by_class_malaria_year/probability/LRI_", mal, "4.csv"), row.names = FALSE)
}

dir.create(
  "outputs/by_class_df/by_malaria_year/probability",
  recursive = TRUE,
  showWarnings = FALSE
)

write.csv(treated_LRI_df, "outputs/by_class_df/by_malaria_year/probability/LRI_all4.csv")

for (mal in Malaria_name) {
  write.csv(treated_LRI_df[treated_LRI_df$Malaria == mal, ], paste0("outputs/by_class_df/by_malaria_year/probability/LRI_", mal, "4.csv"), row.names = FALSE)
}

rm(treated_LRI_sim, treated_LRI_result, treated_LRI_summary, treated_LRI_df)





#----------------------------1.4 by year and GBD region (treatment) ------------------------------

#---LRI---
treated_LRI_sim <- treated_LRI_matrix %>% 
  group_by(Year, GBD_region, sim) %>%
  summarise(
    Treated_LRI_BSP = sum(Treated_LRI_BSP, na.rm = TRUE),
    Treated_LRI_NSP = sum(Treated_LRI_NSP, na.rm = TRUE),
    Treated_LRI_Cephalosporins = sum(Treated_LRI_Cephalosporins, na.rm = TRUE),
    Treated_LRI_Aminoglycosides = sum(Treated_LRI_Aminoglycosides, na.rm = TRUE),
    Incidence = sum(Incidence_LRI, na.rm = TRUE),
    .groups = "drop"
  )

treated_LRI_result <- treated_LRI_sim %>%
  group_by(Year, GBD_region) %>%
  summarise(
    
    proportion_LRI_BSP_med = median(Treated_LRI_BSP/Incidence, na.rm = TRUE),
    proportion_LRI_BSP_lower = quantile(Treated_LRI_BSP/Incidence, 0.025, na.rm = TRUE),
    proportion_LRI_BSP_upper = quantile(Treated_LRI_BSP/Incidence, 0.975, na.rm = TRUE),
    
    proportion_LRI_NSP_med = median(Treated_LRI_NSP/Incidence, na.rm = TRUE),
    proportion_LRI_NSP_lower = quantile(Treated_LRI_NSP/Incidence, 0.025, na.rm = TRUE),
    proportion_LRI_NSP_upper = quantile(Treated_LRI_NSP/Incidence, 0.975, na.rm = TRUE),
    
    proportion_LRI_Cephalosporins_med = median(Treated_LRI_Cephalosporins/Incidence, na.rm = TRUE),
    proportion_LRI_Cephalosporins_lower = quantile(Treated_LRI_Cephalosporins/Incidence, 0.025, na.rm = TRUE),
    proportion_LRI_Cephalosporins_upper = quantile(Treated_LRI_Cephalosporins/Incidence, 0.975, na.rm = TRUE),
    
    proportion_LRI_Aminoglycosides_med = median(Treated_LRI_Aminoglycosides/Incidence, na.rm = TRUE),
    proportion_LRI_Aminoglycosides_lower = quantile(Treated_LRI_Aminoglycosides/Incidence, 0.025, na.rm = TRUE),
    proportion_LRI_Aminoglycosides_upper = quantile(Treated_LRI_Aminoglycosides/Incidence, 0.975, na.rm = TRUE),

    .groups = "drop"
  )


abx_names <- c(
  "proportion_LRI_BSP", "proportion_LRI_NSP", "proportion_LRI_Cephalosporins", "proportion_LRI_Aminoglycosides")

for (i in 1:length(abx_names)) {
  med_col <- paste0(abx_names[i], "_med")
  lwr_col <- paste0(abx_names[i], "_lower")
  upr_col <- paste0(abx_names[i], "_upper")
  fmt_col <- paste0(abx_names[i])
  
  treated_LRI_result[[fmt_col]] <- sprintf(
    "%s (%s, %s)",
    formatC(treated_LRI_result[[med_col]], format = "f", big.mark = " ", digits = digit_number[i]),
    formatC(treated_LRI_result[[lwr_col]], format = "f", big.mark = " ", digits = digit_number[i]),
    formatC(treated_LRI_result[[upr_col]], format = "f", big.mark = " ", digits = digit_number[i])
  )
}

treated_LRI_summary <- treated_LRI_result[,c(1:2,15:18)]
treated_LRI_df <- treated_LRI_result[,c(1:14)]



#### 1.4 CSV file ####
dir.create(
  "outputs/table/by_class_gbd_region_year/probability",
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  "outputs/by_class_df/by_gbd_region_year/probability",
  recursive = TRUE,
  showWarnings = FALSE
)

write.csv(treated_LRI_summary, "outputs/table/by_class_gbd_region_year/probability/LRI4.csv")
write.csv(treated_LRI_df, "outputs/by_class_df/by_gbd_region_year/probability/LRI_all_Region4.csv")

GBD_region_name <- unique(treated_LRI_df$GBD_region)
for (gbd in GBD_region_name) {
  write.csv(treated_LRI_summary[treated_LRI_summary$GBD_region == gbd, ], paste0("outputs/table/by_class_gbd_region_year/probability/LRI_", gbd, "4.csv"), row.names = FALSE)
}

GBD_region_name <- unique(treated_LRI_df$GBD_region)
for (gbd in GBD_region_name) {
  write.csv(treated_LRI_df[treated_LRI_df$GBD_region == gbd, ], paste0("outputs/by_class_df/by_gbd_region_year/probability/LRI_", gbd, "4.csv"), row.names = FALSE)
}

rm(treated_LRI_sim, treated_LRI_result, treated_LRI_summary, treated_LRI_df)






#----------------------------1.5 by year and Income (treatment-LRI) ------------------------------


first_sim <- treated_LRI_matrix[1:(nrow(treated_LRI_matrix)/10000),]
first_sim <- first_sim[,c(1,2,5)]
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
treated_LRI_matrix <- treated_LRI_matrix %>% inner_join(income_dist_summary[,c("Country","First5","Last5")], by = "Country") 
rm(income_dist, income_dist_summary, first_sim, first5_years, last5_years)




#---LRI---
treated_LRI_sim <- treated_LRI_matrix %>% 
  group_by(Year, First5, sim) %>%
  summarise(
    Treated_LRI_BSP = sum(Treated_LRI_BSP, na.rm = TRUE),
    Treated_LRI_NSP = sum(Treated_LRI_NSP, na.rm = TRUE),
    Treated_LRI_Cephalosporins = sum(Treated_LRI_Cephalosporins, na.rm = TRUE),
    Treated_LRI_Aminoglycosides = sum(Treated_LRI_Aminoglycosides, na.rm = TRUE),
    Incidence = sum(Incidence_LRI, na.rm = TRUE),
    .groups = "drop"
  )

treated_LRI_result <- treated_LRI_sim %>%
  group_by(Year, First5) %>%
  summarise(
    
    proportion_LRI_BSP_med = median(Treated_LRI_BSP/Incidence, na.rm = TRUE),
    proportion_LRI_BSP_lower = quantile(Treated_LRI_BSP/Incidence, 0.025, na.rm = TRUE),
    proportion_LRI_BSP_upper = quantile(Treated_LRI_BSP/Incidence, 0.975, na.rm = TRUE),
    
    proportion_LRI_NSP_med = median(Treated_LRI_NSP/Incidence, na.rm = TRUE),
    proportion_LRI_NSP_lower = quantile(Treated_LRI_NSP/Incidence, 0.025, na.rm = TRUE),
    proportion_LRI_NSP_upper = quantile(Treated_LRI_NSP/Incidence, 0.975, na.rm = TRUE),
    
    proportion_LRI_Cephalosporins_med = median(Treated_LRI_Cephalosporins/Incidence, na.rm = TRUE),
    proportion_LRI_Cephalosporins_lower = quantile(Treated_LRI_Cephalosporins/Incidence, 0.025, na.rm = TRUE),
    proportion_LRI_Cephalosporins_upper = quantile(Treated_LRI_Cephalosporins/Incidence, 0.975, na.rm = TRUE),
   
    proportion_LRI_Aminoglycosides_med = median(Treated_LRI_Aminoglycosides/Incidence, na.rm = TRUE),
    proportion_LRI_Aminoglycosides_lower = quantile(Treated_LRI_Aminoglycosides/Incidence, 0.025, na.rm = TRUE),
    proportion_LRI_Aminoglycosides_upper = quantile(Treated_LRI_Aminoglycosides/Incidence, 0.975, na.rm = TRUE),

    .groups = "drop"
  )



abx_names <- c(
  "proportion_LRI_BSP", "proportion_LRI_NSP", "proportion_LRI_Cephalosporins", "proportion_LRI_Aminoglycosides")

for (i in 1:length(abx_names)) {
  med_col <- paste0(abx_names[i], "_med")
  lwr_col <- paste0(abx_names[i], "_lower")
  upr_col <- paste0(abx_names[i], "_upper")
  fmt_col <- paste0(abx_names[i])
  
  treated_LRI_result[[fmt_col]] <- sprintf(
    "%s (%s, %s)",
    formatC(treated_LRI_result[[med_col]], format = "f", big.mark = " ", digits = digit_number[i]),
    formatC(treated_LRI_result[[lwr_col]], format = "f", big.mark = " ", digits = digit_number[i]),
    formatC(treated_LRI_result[[upr_col]], format = "f", big.mark = " ", digits = digit_number[i])
  )
}

treated_LRI_summary_first5 <- treated_LRI_result[,c(1:2,15:18)]
treated_LRI_df_first5 <- treated_LRI_result[,c(1:14)]



#### 1.5.1 CSV file ####
dir.create(
  "outputs/table/by_class_income_year/probability",
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  "outputs/by_class_df/by_income_year/probability",
  recursive = TRUE,
  showWarnings = FALSE
)

write.csv(treated_LRI_summary_first5, "outputs/table/by_class_income_year/probability/LRI4_first5.csv")
write.csv(treated_LRI_df_first5, "outputs/by_class_df/by_income_year/probability/LRI4_first5.csv")

Income <- unique(treated_LRI_df_first5$First5)
for (inc in Income) {
  write.csv(treated_LRI_summary_first5[treated_LRI_summary_first5$First5 == inc, ], paste0("outputs/table/by_class_income_year/probability/LRI4_first5_", inc, ".csv"), row.names = FALSE)
}

dir.create(
  "outputs/by_class_df/by_income_year/probability",
  recursive = TRUE,
  showWarnings = FALSE
)

Income <- unique(treated_LRI_df_first5$First5)
for (inc in Income) {
  write.csv(treated_LRI_df_first5[treated_LRI_df_first5$First5 == inc, ], paste0("outputs/by_class_df/by_income_year/probability/LRI4_first5_", inc, ".csv"), row.names = FALSE)
}

rm(treated_LRI_result, treated_LRI_df_first5, treated_LRI_summary_first5, treated_LRI_sim)



#---LRI---
treated_LRI_sim <- treated_LRI_matrix %>% 
  group_by(Year, Last5, sim) %>%
  summarise(
    Treated_LRI_BSP = sum(Treated_LRI_BSP, na.rm = TRUE),
    Treated_LRI_NSP = sum(Treated_LRI_NSP, na.rm = TRUE),
    Treated_LRI_Cephalosporins = sum(Treated_LRI_Cephalosporins, na.rm = TRUE),
    Treated_LRI_Aminoglycosides = sum(Treated_LRI_Aminoglycosides, na.rm = TRUE),
    Incidence = sum(Incidence_LRI, na.rm = TRUE),
    .groups = "drop"
  )

treated_LRI_result <- treated_LRI_sim %>%
  group_by(Year, Last5) %>%
  summarise(
    
    proportion_LRI_BSP_med = median(Treated_LRI_BSP/Incidence, na.rm = TRUE),
    proportion_LRI_BSP_lower = quantile(Treated_LRI_BSP/Incidence, 0.025, na.rm = TRUE),
    proportion_LRI_BSP_upper = quantile(Treated_LRI_BSP/Incidence, 0.975, na.rm = TRUE),
    
    proportion_LRI_NSP_med = median(Treated_LRI_NSP/Incidence, na.rm = TRUE),
    proportion_LRI_NSP_lower = quantile(Treated_LRI_NSP/Incidence, 0.025, na.rm = TRUE),
    proportion_LRI_NSP_upper = quantile(Treated_LRI_NSP/Incidence, 0.975, na.rm = TRUE),
    
    proportion_LRI_Cephalosporins_med = median(Treated_LRI_Cephalosporins/Incidence, na.rm = TRUE),
    proportion_LRI_Cephalosporins_lower = quantile(Treated_LRI_Cephalosporins/Incidence, 0.025, na.rm = TRUE),
    proportion_LRI_Cephalosporins_upper = quantile(Treated_LRI_Cephalosporins/Incidence, 0.975, na.rm = TRUE),

    proportion_LRI_Aminoglycosides_med = median(Treated_LRI_Aminoglycosides/Incidence, na.rm = TRUE),
    proportion_LRI_Aminoglycosides_lower = quantile(Treated_LRI_Aminoglycosides/Incidence, 0.025, na.rm = TRUE),
    proportion_LRI_Aminoglycosides_upper = quantile(Treated_LRI_Aminoglycosides/Incidence, 0.975, na.rm = TRUE),
   
    .groups = "drop"
  )


abx_names <- c(
  "proportion_LRI_BSP", "proportion_LRI_NSP", "proportion_LRI_Cephalosporins", "proportion_LRI_Aminoglycosides")

for (i in 1:length(abx_names)) {
  med_col <- paste0(abx_names[i], "_med")
  lwr_col <- paste0(abx_names[i], "_lower")
  upr_col <- paste0(abx_names[i], "_upper")
  fmt_col <- paste0(abx_names[i])
  
  treated_LRI_result[[fmt_col]] <- sprintf(
    "%s (%s, %s)",
    formatC(treated_LRI_result[[med_col]], format = "f", big.mark = " ", digits = digit_number[i]),
    formatC(treated_LRI_result[[lwr_col]], format = "f", big.mark = " ", digits = digit_number[i]),
    formatC(treated_LRI_result[[upr_col]], format = "f", big.mark = " ", digits = digit_number[i])
  )
}

treated_LRI_summary_last5 <- treated_LRI_result[,c(1:2,15:18)]
treated_LRI_df_last5 <- treated_LRI_result[,c(1:14)]




#### 1.5.2 CSV file ####
write.csv(treated_LRI_summary_last5, "outputs/table/by_class_income_year/probability/LRI4_last5.csv")
write.csv(treated_LRI_df_last5, "outputs/by_class_df/by_income_year/probability/LRI4_last5.csv")

Income <- unique(treated_LRI_df_last5$Last5)
for (inc in Income) {
  write.csv(treated_LRI_summary_last5[treated_LRI_summary_last5$Last5 == inc, ], paste0("outputs/table/by_class_income_year/probability/LRI4_last5_", inc, ".csv"), row.names = FALSE)
}

Income <- unique(treated_LRI_df_last5$Last5)
for (inc in Income) {
  write.csv(treated_LRI_df_last5[treated_LRI_df_last5$Last5 == inc, ], paste0("outputs/by_class_df/by_income_year/probability/LRI4_last5_", inc, ".csv"), row.names = FALSE)
}

rm(treated_LRI_result, treated_LRI_df_last5, treated_LRI_summary_last5, treated_LRI_sim)
