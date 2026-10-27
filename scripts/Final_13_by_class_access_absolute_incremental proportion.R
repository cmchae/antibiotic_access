##### 0 #####
rm(list = ls())
gc()

# Optimized strictly for dataframe manipulation
requiredPackages = c("dplyr", "tidyr", "tibble")
lib.path <- .libPaths()
for(p in requiredPackages){
  if(!require(p,character.only = TRUE)) install.packages(p, lib = lib.path)
  library(p,character.only = TRUE)
}


##### 1 #####
n_sim <- 10000
set.seed(123)

# confirm upstream antibiotic-treatment model uses PC1-PC2
coef_df_LRI <- readRDS("outputs/coef_df_LRI.rds")

stopifnot(
  identical(
    names(coef_df_LRI),
    c("(Intercept)", "PC1", "PC2")
  )
)

rm(coef_df_LRI)

# The following objects were generated using the current PC1-PC2 models
Country_Year_meta <- readRDS("outputs/Country_Year_meta4.rds")
sim_incidence <- readRDS("outputs/sim_incidence4.rds")
treated_LRI_matrix <- readRDS("outputs/treated_LRI_matrix4.rds")
proportion_LRI_matrix <- readRDS("outputs/proportion_LRI_matrix4.rds")

treated_LRI_matrix$Treated_all <- rowSums(treated_LRI_matrix[,8:11])
proportion_LRI_matrix$Proportion_LRI_all <- rowSums(proportion_LRI_matrix[,8:11])

head(treated_LRI_matrix)
unique(treated_LRI_matrix$Year)
unique(treated_LRI_matrix$sim)
digit_number <- c(3,4,3,3)


#----------------------------1.2 by year (treatment) ------------------------------
# Incremental antibiotic share (LRI)

baseline_year <- 2003

#--- 1. Aggregate by Year × sim ---
prop_LRI_sim <- proportion_LRI_matrix %>%
  group_by(Year, sim) %>%
  summarise(
    BSP   = mean(proportion_LRI_BSP, na.rm = TRUE),
    NSP   = mean(proportion_LRI_NSP, na.rm = TRUE),
    CEPH  = mean(proportion_LRI_Cephalosporins, na.rm = TRUE),
    AMINO = mean(proportion_LRI_Aminoglycosides, na.rm = TRUE),
    ALL   = mean(Proportion_LRI_all, na.rm = TRUE),
    .groups = "drop"
  )

#--- 2. Extract baseline per simulation ---
baseline_2003 <- prop_LRI_sim %>%
  filter(Year == baseline_year) %>%
  dplyr::select(
    sim,
    BSP_2003   = BSP,
    NSP_2003   = NSP,
    CEPH_2003  = CEPH,
    AMINO_2003 = AMINO,
    ALL_2003   = ALL
  )

#--- 3. Compute incremental shares per simulation ---
prop_LRI_incremental <- prop_LRI_sim %>%
  left_join(baseline_2003, by = "sim") %>%
  mutate(
    d_ALL = ALL - ALL_2003
  ) %>%
  mutate(
    inc_BSP   = (BSP   - BSP_2003)   / d_ALL,
    inc_NSP   = (NSP   - NSP_2003)   / d_ALL,
    inc_CEPH  = (CEPH  - CEPH_2003)  / d_ALL,
    inc_AMINO = (AMINO - AMINO_2003) / d_ALL
  ) %>%
  # protect against division by ~0
  mutate(
    across(starts_with("inc_"),
           ~ ifelse(abs(d_ALL) < 1e-6, NA, .))
  )

#--- 4. Summarise across simulations by year ---
treated_LRI_result <- prop_LRI_incremental %>%
  group_by(Year) %>%
  summarise(
    inc_LRI_BSP_med   = median(inc_BSP, na.rm = TRUE),
    inc_LRI_BSP_lower = quantile(inc_BSP, 0.025, na.rm = TRUE),
    inc_LRI_BSP_upper = quantile(inc_BSP, 0.975, na.rm = TRUE),
    
    inc_LRI_NSP_med   = median(inc_NSP, na.rm = TRUE),
    inc_LRI_NSP_lower = quantile(inc_NSP, 0.025, na.rm = TRUE),
    inc_LRI_NSP_upper = quantile(inc_NSP, 0.975, na.rm = TRUE),
    
    inc_LRI_Cephalosporins_med   = median(inc_CEPH, na.rm = TRUE),
    inc_LRI_Cephalosporins_lower = quantile(inc_CEPH, 0.025, na.rm = TRUE),
    inc_LRI_Cephalosporins_upper = quantile(inc_CEPH, 0.975, na.rm = TRUE),
    
    inc_LRI_Aminoglycosides_med   = median(inc_AMINO, na.rm = TRUE),
    inc_LRI_Aminoglycosides_lower = quantile(inc_AMINO, 0.025, na.rm = TRUE),
    inc_LRI_Aminoglycosides_upper = quantile(inc_AMINO, 0.975, na.rm = TRUE),
    
    .groups = "drop"
  )

#--- 5. Format for table output ---
abx_names <- c(
  "inc_LRI_BSP",
  "inc_LRI_NSP",
  "inc_LRI_Cephalosporins",
  "inc_LRI_Aminoglycosides"
)

for (i in seq_along(abx_names)) {
  med_col <- paste0(abx_names[i], "_med")
  lwr_col <- paste0(abx_names[i], "_lower")
  upr_col <- paste0(abx_names[i], "_upper")
  fmt_col <- abx_names[i]
  
  treated_LRI_result[[fmt_col]] <- sprintf(
    "%s (%s, %s)",
    formatC(treated_LRI_result[[med_col]], format = "f",
            digits = digit_number[i]),
    formatC(treated_LRI_result[[lwr_col]], format = "f",
            digits = digit_number[i]),
    formatC(treated_LRI_result[[upr_col]], format = "f",
            digits = digit_number[i])
  )
}

#--- 6. Final outputs ---
treated_LRI_summary <- treated_LRI_result %>%
  filter(Year != baseline_year) %>%
  dplyr::select(Year, all_of(abx_names))

treated_LRI_df <- treated_LRI_result %>%
  filter(Year != baseline_year)

#--- 7. Write CSVs ---
write.csv(
  treated_LRI_summary,
  "outputs/table/by_class_year/probability/LRI_incremental_share_prop.csv",
  row.names = FALSE
)

write.csv(
  treated_LRI_df,
  "outputs/by_class_df/by_year/probability/LRI_incremental_share_df_prop.csv",
  row.names = FALSE
)

#--- 8. Clean up ---
rm(
  prop_LRI_sim,
  prop_LRI_incremental,
  baseline_2003,
  treated_LRI_result,
  treated_LRI_summary,
  treated_LRI_df
)



#----------------------------1.4 by year and GBD region (treatment) ------------------------------

baseline_year <- 2003

prop_LRI_sim <- proportion_LRI_matrix %>%
  group_by(Year, GBD_region, sim) %>%
  summarise(
    BSP   = mean(proportion_LRI_BSP, na.rm = TRUE),
    NSP   = mean(proportion_LRI_NSP, na.rm = TRUE),
    CEPH  = mean(proportion_LRI_Cephalosporins, na.rm = TRUE),
    AMINO = mean(proportion_LRI_Aminoglycosides, na.rm = TRUE),
    ALL   = mean(Proportion_LRI_all, na.rm = TRUE),
    .groups = "drop"
  )

baseline_2003 <- prop_LRI_sim %>%
  filter(Year == baseline_year) %>%
  dplyr::select(
    GBD_region, sim,
    BSP_2003   = BSP,
    NSP_2003   = NSP,
    CEPH_2003  = CEPH,
    AMINO_2003 = AMINO,
    ALL_2003   = ALL
  )

prop_LRI_incremental <- prop_LRI_sim %>%
  left_join(baseline_2003, by = c("GBD_region", "sim")) %>%
  mutate(
    d_ALL = ALL - ALL_2003
  ) %>%
  mutate(
    inc_BSP   = (BSP   - BSP_2003)   / d_ALL,
    inc_NSP   = (NSP   - NSP_2003)   / d_ALL,
    inc_CEPH  = (CEPH  - CEPH_2003)  / d_ALL,
    inc_AMINO = (AMINO - AMINO_2003) / d_ALL
  ) %>%
  mutate(
    across(starts_with("inc_"),
           ~ ifelse(abs(d_ALL) < 1e-6, NA, .))
  )

treated_LRI_result <- prop_LRI_incremental %>%
  group_by(Year, GBD_region) %>%
  summarise(
    inc_LRI_BSP_med   = median(inc_BSP, na.rm = TRUE),
    inc_LRI_BSP_lower = quantile(inc_BSP, 0.025, na.rm = TRUE),
    inc_LRI_BSP_upper = quantile(inc_BSP, 0.975, na.rm = TRUE),
    
    inc_LRI_NSP_med   = median(inc_NSP, na.rm = TRUE),
    inc_LRI_NSP_lower = quantile(inc_NSP, 0.025, na.rm = TRUE),
    inc_LRI_NSP_upper = quantile(inc_NSP, 0.975, na.rm = TRUE),
    
    inc_LRI_Cephalosporins_med   = median(inc_CEPH, na.rm = TRUE),
    inc_LRI_Cephalosporins_lower = quantile(inc_CEPH, 0.025, na.rm = TRUE),
    inc_LRI_Cephalosporins_upper = quantile(inc_CEPH, 0.975, na.rm = TRUE),
    
    inc_LRI_Aminoglycosides_med   = median(inc_AMINO, na.rm = TRUE),
    inc_LRI_Aminoglycosides_lower = quantile(inc_AMINO, 0.025, na.rm = TRUE),
    inc_LRI_Aminoglycosides_upper = quantile(inc_AMINO, 0.975, na.rm = TRUE),
    
    .groups = "drop"
  )

abx_names <- c(
  "inc_LRI_BSP",
  "inc_LRI_NSP",
  "inc_LRI_Cephalosporins",
  "inc_LRI_Aminoglycosides"
)

for (i in seq_along(abx_names)) {
  med_col <- paste0(abx_names[i], "_med")
  lwr_col <- paste0(abx_names[i], "_lower")
  upr_col <- paste0(abx_names[i], "_upper")
  fmt_col <- abx_names[i]
  
  treated_LRI_result[[fmt_col]] <- sprintf(
    "%s (%s, %s)",
    formatC(treated_LRI_result[[med_col]], format = "f",
            digits = digit_number[i]),
    formatC(treated_LRI_result[[lwr_col]], format = "f",
            digits = digit_number[i]),
    formatC(treated_LRI_result[[upr_col]], format = "f",
            digits = digit_number[i])
  )
}

treated_LRI_summary <- treated_LRI_result %>%
  filter(Year != baseline_year) %>%
  dplyr::select(Year, GBD_region, all_of(abx_names))

treated_LRI_df <- treated_LRI_result %>%
  filter(Year != baseline_year)


#### 1.4 CSV file ####
write.csv(
  treated_LRI_summary,
  "outputs/table/by_class_gbd_region_year/probability/LRI_incremental_prop.csv",
  row.names = FALSE
)

write.csv(
  treated_LRI_df,
  "outputs/by_class_df/by_gbd_region_year/probability/LRI_incremental_df_prop.csv",
  row.names = FALSE
)

GBD_region_name <- unique(treated_LRI_df$GBD_region)
for (gbd in GBD_region_name) {
  write.csv(treated_LRI_summary[treated_LRI_summary$GBD_region == gbd, ], paste0("outputs/table/by_class_gbd_region_year/probability/LRI_prop_", gbd, "4.csv"), row.names = FALSE)
}

GBD_region_name <- unique(treated_LRI_df$GBD_region)
for (gbd in GBD_region_name) {
  write.csv(treated_LRI_df[treated_LRI_df$GBD_region == gbd, ], paste0("outputs/by_class_df/by_gbd_region_year/probability/LRI_prop_", gbd, "4.csv"), row.names = FALSE)
}

rm(
  prop_LRI_sim,
  prop_LRI_incremental,
  baseline_2003,
  treated_LRI_result,
  treated_LRI_summary,
  treated_LRI_df
)



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
baseline_year <- 2003

prop_LRI_sim <- proportion_LRI_matrix %>%
  inner_join(
    treated_LRI_matrix %>% distinct(Country, Last5),
    by = "Country"
  ) %>%
  group_by(Year, Last5, sim) %>%
  summarise(
    BSP   = mean(proportion_LRI_BSP, na.rm = TRUE),
    NSP   = mean(proportion_LRI_NSP, na.rm = TRUE),
    CEPH  = mean(proportion_LRI_Cephalosporins, na.rm = TRUE),
    AMINO = mean(proportion_LRI_Aminoglycosides, na.rm = TRUE),
    ALL   = mean(Proportion_LRI_all, na.rm = TRUE),
    .groups = "drop"
  )

baseline_2003 <- prop_LRI_sim %>%
  filter(Year == baseline_year) %>%
  dplyr::select(
    Last5, sim,
    BSP_2003   = BSP,
    NSP_2003   = NSP,
    CEPH_2003  = CEPH,
    AMINO_2003 = AMINO,
    ALL_2003   = ALL
  )

prop_LRI_incremental <- prop_LRI_sim %>%
  left_join(baseline_2003, by = c("Last5", "sim")) %>%
  mutate(
    d_ALL = ALL - ALL_2003,
    inc_BSP   = (BSP   - BSP_2003)   / d_ALL,
    inc_NSP   = (NSP   - NSP_2003)   / d_ALL,
    inc_CEPH  = (CEPH  - CEPH_2003)  / d_ALL,
    inc_AMINO = (AMINO - AMINO_2003) / d_ALL
  ) %>%
  mutate(
    across(starts_with("inc_"),
           ~ ifelse(abs(d_ALL) < 1e-6, NA, .))
  )

treated_LRI_result <- prop_LRI_incremental %>%
  group_by(Year, Last5) %>%
  summarise(
    inc_LRI_BSP_med   = median(inc_BSP, na.rm = TRUE),
    inc_LRI_BSP_lower = quantile(inc_BSP, 0.025, na.rm = TRUE),
    inc_LRI_BSP_upper = quantile(inc_BSP, 0.975, na.rm = TRUE),
    
    inc_LRI_NSP_med   = median(inc_NSP, na.rm = TRUE),
    inc_LRI_NSP_lower = quantile(inc_NSP, 0.025, na.rm = TRUE),
    inc_LRI_NSP_upper = quantile(inc_NSP, 0.975, na.rm = TRUE),
    
    inc_LRI_Cephalosporins_med   = median(inc_CEPH, na.rm = TRUE),
    inc_LRI_Cephalosporins_lower = quantile(inc_CEPH, 0.025, na.rm = TRUE),
    inc_LRI_Cephalosporins_upper = quantile(inc_CEPH, 0.975, na.rm = TRUE),
    
    inc_LRI_Aminoglycosides_med   = median(inc_AMINO, na.rm = TRUE),
    inc_LRI_Aminoglycosides_lower = quantile(inc_AMINO, 0.025, na.rm = TRUE),
    inc_LRI_Aminoglycosides_upper = quantile(inc_AMINO, 0.975, na.rm = TRUE),
    
    .groups = "drop"
  )

abx_names <- c(
  "inc_LRI_BSP",
  "inc_LRI_NSP",
  "inc_LRI_Cephalosporins",
  "inc_LRI_Aminoglycosides"
)

for (i in seq_along(abx_names)) {
  med_col <- paste0(abx_names[i], "_med")
  lwr_col <- paste0(abx_names[i], "_lower")
  upr_col <- paste0(abx_names[i], "_upper")
  fmt_col <- abx_names[i]
  
  treated_LRI_result[[fmt_col]] <- sprintf(
    "%s (%s, %s)",
    formatC(treated_LRI_result[[med_col]], format = "f",
            digits = digit_number[i]),
    formatC(treated_LRI_result[[lwr_col]], format = "f",
            digits = digit_number[i]),
    formatC(treated_LRI_result[[upr_col]], format = "f",
            digits = digit_number[i])
  )
}

treated_LRI_summary_last5 <- treated_LRI_result %>%
  filter(Year != baseline_year) %>%
  dplyr::select(Year, Last5, all_of(abx_names))

treated_LRI_df_last5 <- treated_LRI_result %>%
  filter(Year != baseline_year)


#### 1.5.2 CSV file ####
write.csv(
  treated_LRI_summary_last5,
  "outputs/table/by_class_income_year/probability/LRI_incremental_last5_prop.csv",
  row.names = FALSE
)

write.csv(
  treated_LRI_df_last5,
  "outputs/by_class_df/by_income_year/probability/LRI_incremental_last5_df_prop.csv",
  row.names = FALSE
)

Income <- unique(treated_LRI_df_last5$Last5)
for (inc in Income) {
  write.csv(treated_LRI_summary_last5[treated_LRI_summary_last5$Last5 == inc, ], paste0("outputs/table/by_class_income_year/probability/LRI4_last5_prop_", inc, ".csv"), row.names = FALSE)
}

Income <- unique(treated_LRI_df_last5$Last5)
for (inc in Income) {
  write.csv(treated_LRI_df_last5[treated_LRI_df_last5$Last5 == inc, ], paste0("outputs/by_class_df/by_income_year/probability/LRI4_last5_prop_", inc, ".csv"), row.names = FALSE)
}

rm(
  prop_LRI_sim,
  prop_LRI_incremental,
  baseline_2003,
  treated_LRI_result,
  treated_LRI_df_last5,
  treated_LRI_summary_last5
)

