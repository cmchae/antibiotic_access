rm(list = ls())
gc()

options(scipen = 999, digits = 4)

requiredPackages <- c("dplyr", "tidyr", "glmmTMB")
lib.path <- .libPaths()
for(p in requiredPackages){
  if(!require(p,character.only = TRUE)) install.packages(p, lib = lib.path)
  library(p,character.only = TRUE)
}

n_sim <- 10000

# ----------- Data from literautre review dataset ------------
prop.dat <- read.csv("data/ver.2_res.fev.prop.csv", header = TRUE)
UNI <- readRDS("outputs/pc_gdp_ten.rds")
country_mapping <- c("Yemen" = "Yemen, Rep.", "Iran" = "Iran, Islamic Rep.", "Turkey" = "Turkiye")  
prop.dat <- prop.dat %>% mutate(Country = ifelse(Country %in% names(country_mapping), country_mapping[Country], Country))
prop.data <- prop.dat %>% inner_join(UNI, by = c("Country","Year"))
prop.data <- prop.data %>% distinct() %>% na.omit()
UNI.data <- UNI[,1:3]

# Malaria proportion :IHME data (GBD2021, adjust the proportion of malaria in 2021 to 2022-2024) - https://vizhub.healthdata.org/gbd-results/
mal.prop <- read.csv("data/malaria_proportion.csv")
mal.prop <- mal.prop[,-1]
colnames(mal.prop) <- c("Country", "Age", "Year", "mal.prop", "RI.prop")

country_mapping_mal <- c("Taiwan (Province of China)" = "Taiwan",   "Lao People's Democratic Republic" = "Lao PDR", "Kyrgyzstan" = "Kyrgyz Republic", "Democratic People's Republic of Korea" = "Korea (Democratic People's Rep. of)", "Micronesia (Federated States of)" = "Micronesia, Fed. Sts.", "Republic of Moldova" = "Moldova", "Republic of Korea" = "Korea, Rep.", "Slovakia" = "Slovak Republic", "Bahamas" = "Bahamas, The", "Venezuela (Bolivarian Republic of)" = "Venezuela, RB", "United States of America" = "United States", "Saint Vincent and the Grenadines" = "St. Vincent and the Grenadines", "Bolivia (Plurinational State of)" = "Bolivia", "Saint Lucia" = "St. Lucia", "Iran (Islamic Republic of)" = "Iran, Islamic Rep.", "Palestine" = "West Bank and Gaza", "Egypt" = "Egypt, Arab Rep.", "Türkiye" = "Turkiye", "Yemen" = "Yemen, Rep.", "Congo" = "Congo, Rep.", "Democratic Republic of the Congo" = "Congo, Dem. Rep.", "United Republic of Tanzania" = "Tanzania", "Côte d'Ivoire" = "Cote d'Ivoire", "Gambia" = "Gambia, The", "Saint Kitts and Nevis" = "St. Kitts and Nevis", "United States Virgin Islands" = "Virgin Islands (U.S.)"
)

mal.prop <- mal.prop %>% mutate(Country = ifelse(Country %in% names(country_mapping_mal), country_mapping_mal[Country], Country))

mal.prop_all.ages <- mal.prop[mal.prop$Age == "All ages", ]

prop.data <- prop.data %>% inner_join(mal.prop_all.ages, by = c("Country" = "Country", "Year" = "Year"))
UNI.data <- UNI.data %>% inner_join(mal.prop_all.ages, by = c("Country" = "Country", "Year" = "Year"))



#Temperature: World bank - Climate Change Knowledge Portal, https://climateknowledgeportal.worldbank.org/
Temperature <- read.csv("data/Temperature_2024.csv")
Temperature <- Temperature %>% pivot_longer(cols = starts_with("X"), names_to = "year_month", values_to = "temperature")  # Convert wide month columns into long format
Temperature <- Temperature %>%mutate(year_month = sub("^X", "", year_month)) %>%  separate(year_month, into = c("Year", "Month"), sep = "\\.") # Remove leading "X" and split into Year and Month
Temperature <- Temperature %>% mutate(Year = as.integer(Year), Month = as.integer(Month)) # Convert types
Temperature <- Temperature %>% group_by(name, Year, Month) %>% summarise(temperature = mean(temperature), .groups = "drop") # merging if there are multiple station in one country
Temperature <- Temperature %>% pivot_wider(names_from = Month, names_prefix = "Month", values_from = temperature) # Re-spread long data to have Month columns
Temperature <- Temperature %>%  arrange(name, Year) %>% rename(Country = name) # Rearrange columns: Country, Year, Month1–Month12
Temperature$Temperature <- apply(Temperature[, 3:14], 1, mean) # Mean of monthly temperature
Temperature$Temp_maximum <- apply(Temperature[, 3:14], 1, max) # Maximum of monthly temperature
plot(Temperature$Temperature, Temperature$Temp_maximum)
Temperature <- Temperature[, c("Country", "Year", "Temperature")]

country_mapping_Temp <- c("Aruba (Neth.)" = "Aruba", "The Bahamas" = "Bahamas, The", "Bermuda (U.K.)" = "Bermuda", "Cayman Islands (U.K.)" = "Cayman Islands", "Democratic Republic of Congo" = "Congo, Dem. Rep.", "Congo" = "Congo, Rep.", "C\x99te d'Ivoire" = "Cote d'Ivoire", "Cura\x8dao (Neth.)" = "Curacao", "Czech Republic" = "Czechia", "Arab Republic of Egypt" = "Egypt, Arab Rep.", "Faroe Islands (Den.)" = "Faroe Islands", "The Gambia" = "Gambia, The", "Greenland (Den.)" = "Greenland", "Hong Kong, SAR" = "Hong Kong SAR, China", "Islamic Republic of Iran" = "Iran, Islamic Rep.", "Republic of Korea" = "Korea, Rep.", "Lao People's Democratic Republic" = "Lao PDR", "Macau, SAR" = "Macao SAR, China", "Federated States of Micronesia" = "Micronesia, Fed. Sts.", "Puerto Rico (U.S.)" = "Puerto Rico", "S\x8bo Tom\x8e and Pr\x92ncipe" = "Sao Tome and Principe", "Sint Maarten (Neth.)" = "Sint Maarten (Dutch part)", "Saint Kitts and Nevis" = "St. Kitts and Nevis", "Saint Lucia" = "St. Lucia", "Saint Vincent and the Grenadines" = "St. Vincent and the Grenadines", "T\x9frkiye" = "Turkiye", "Turks and Caicos Islands (U.K.)" = "Turks and Caicos Islands", "United States of America" = "United States", "R. B. de Venezuela" = "Venezuela, RB", "Vietnam" = "Viet Nam", "United States Virgin Islands (U.S.)" = "Virgin Islands (U.S.)", "Republic of Yemen" = "Yemen, Rep."
)

Temperature <- Temperature %>% mutate(Country = ifelse(Country %in% names(country_mapping_Temp), country_mapping_Temp[Country], Country))
prop.data <- prop.data %>% inner_join(Temperature, by = c("Country" = "Country", "Year" = "Year"))
UNI.data <- UNI.data %>% inner_join(Temperature, by = c("Country" = "Country", "Year" = "Year"))



#Precipitation: Our world in data - Annual precipitation, 1940 to 2024, https://ourworldindata.org/grapher/average-precipitation-per-year?tab=table
Precipitation <- read.csv("data/Precipitation.csv")
Precipitation <- as.data.frame(Precipitation)
Precipitation$Year <- as.numeric(Precipitation$Year) 

country_mapping_Prec <- c("Bahamas" = "Bahamas, The", "Brunei" = "Brunei Darussalam", "Cape Verde" = "Cabo Verde", "Democratic Republic of Congo" = "Congo, Dem. Rep.", "Congo" = "Congo, Rep.", "Curacao" = "Curacao", "Dominica" = "Dominica", "Egypt" = "Egypt, Arab Rep.", "Gambia" = "Gambia, The", "Grenada" = "Grenada", "Hong Kong" = "Hong Kong SAR, China", "Iran" = "Iran, Islamic Rep.", "South Korea" = "Korea, Rep.", "Kyrgyzstan" = "Kyrgyz Republic", "Laos" = "Lao PDR", "Macau" = "Macao SAR, China", "Micronesia" = "Micronesia, Fed. Sts.", "Russia" = "Russian Federation", "Slovakia" = "Slovak Republic", "Saint Vincent and the Grenadines" = "St. Vincent and the Grenadines", "Saint Kitts and Nevis" = "St. Kitts and Nevis", "Saint Lucia" = "St. Lucia", "Syria" = "Syrian Arab Republic", "East Timor" = "Timor-Leste", "Turkey" = "Turkiye", "Turks and Caicos Islands" = "Turks and Caicos Islands", "United States Virgin Islands" = "Virgin Islands (U.S.)", "Venezuela" = "Venezuela, RB", "Vietnam" = "Viet Nam", "Yemen" = "Yemen, Rep.", "Palestine" = "West Bank and Gaza"
)

Precipitation <- Precipitation %>% mutate(Country = ifelse(Country %in% names(country_mapping_Prec), country_mapping_Prec[Country], Country))
prop.data <- prop.data %>% inner_join(Precipitation, by = c("Country" = "Country", "Year" = "Year"))
UNI.data <- UNI.data %>% inner_join(Precipitation, by = c("Country" = "Country", "Year" = "Year"))



### RI/Fever:
set.seed(123); coef_prop <- list(); vcov_prop <- list()

prop.data$Study_ID <- as.factor(prop.data$Study_ID)

for (i in 1:n_sim) {
  sim_data <- prop.data %>%
    mutate(Proportion_sim = pmin(pmax(rnorm(n(), Proportion, sqrt(Proportion * (1 - Proportion) / Sample.size)), 1e-5), 1 - 1e-5))
  mod_sim <- glmmTMB(Proportion_sim ~ Age1 + Age2 + Age3 + w.o.malaria + Year + GDP + Temperature + Precipitation + mal.prop + (1|Study_ID), data = sim_data, family = beta_family(link = "logit"))
  
  if (!is.null(mod_sim)) {
    coef_prop[[i]] <- fixef(mod_sim)$cond
    vcov_prop[[i]] <- vcov(mod_sim)$cond
    
    if (i %% (n_sim / 100) == 0) {
      progress <- round(i / n_sim * 100)
      current_time <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
      print(paste0(current_time, " - ", progress, "% completed"))
    }
  }
}


new.data.UNI <- UNI.data %>% dplyr::select(c("Year", "GDP","Temperature","Precipitation","mal.prop"))
new.data.UNI = cbind(1, 0, 0, 0, new.data.UNI) #, UNI.data[,3:30] (UNI index)
colnames(new.data.UNI)[1:4] <- c("Age1", "Age2", "Age3", "w.o.malaria")
new.data.UNI = cbind(`(Intercept)` = 1, new.data.UNI)

X <- as.matrix(new.data.UNI)
n_sim <- length(coef_prop); n_obs <- nrow(X)
pred.prop <- matrix(NA, nrow = n_obs, ncol = n_sim)

for (i in 1:n_sim) {
  beta_i <- unlist(coef_prop[[i]])
  beta_i <- beta_i[match(colnames(X), names(beta_i))]
  pred.prop[, i] <- X %*% beta_i
  if (i %% (n_sim / 100) == 0) {
    progress <- round(i / n_sim * 100)
    current_time <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
    print(paste0(current_time, " - ", progress, "% completed"))
  }
}

pred.prop <- plogis(pred.prop)

pred.prop.med <- apply(pred.prop, 1, median)
pred.prop.lower <- apply(pred.prop, 1, quantile, probs = 0.025)
pred.prop.upper <- apply(pred.prop, 1, quantile, probs = 0.975)

# Step 6: Combine into summary table
pred.summary <- data.frame(
  Country = UNI.data$Country,
  Year = UNI.data$Year,
  Median = round(pred.prop.med, 3),
  Lower = round(pred.prop.lower, 3),
  Upper = round(pred.prop.upper, 3)
)

prop.RI.Fever <- data.frame(
  Country = UNI.data$Country,
  Year = UNI.data$Year,
  RI_Fever = pred.prop)

saveRDS(prop.RI.Fever, "outputs/prop.RI.Fever.RE.rds")
