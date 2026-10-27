rm(list = ls())
gc()

requiredPackages = c("dplyr","tidyr", "MASS")
lib.path <- .libPaths()
for(p in requiredPackages){
  if(!require(p,character.only = TRUE)) install.packages(p, lib = lib.path)
  library(p,character.only = TRUE)
}
n_sim <- 10000

analysis_data <- readRDS("outputs/analysis_data_prop.rds")
UNI <- readRDS("outputs/pc_gdp_ten.rds")

prop.dat <- analysis_data %>% dplyr::select(c(Country, Year, Fever_LRI, Fever_URI, Fever_LRI_SE, Fever_URI_SE, n_Fever_LRI, n_Fever_URI, n_LRI, n_URI))
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

mal.prop_under5 <- mal.prop[mal.prop$Age == "<5 years", ]
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
Temperature_filtered <- Temperature %>% filter(Country %in% b)
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



### Fever/URI: Logistic binomial --------------
set.seed(123)
mod_sim_URI_logit <- glm(
  cbind(n_Fever_URI,n_URI-n_Fever_URI) ~ Year + GDP + Temperature + Precipitation + mal.prop,
               family = binomial(link="logit"), data = prop.data)

coef_URI <- coef(mod_sim_URI_logit); vcov_URI <- vcov(mod_sim_URI_logit)

new.data.UNI <- UNI.data %>% dplyr::select(c("Year", "GDP", "Temperature", "Precipitation", "mal.prop"))
X <- cbind(1, as.matrix(new.data.UNI))  # Add intercept manually
n_obs <- nrow(X)
beta_sim <- MASS::mvrnorm(n = n_sim, mu = coef_URI, Sigma = vcov_URI)  
logit_p_sim <- X %*% t(beta_sim)
pred.prop.URI <- plogis(logit_p_sim)

pred.prop.med_URI <- apply(pred.prop.URI, 1, median)
pred.prop.lower_URI <- apply(pred.prop.URI, 1, quantile, probs = 0.025)
pred.prop.upper_URI <- apply(pred.prop.URI, 1, quantile, probs = 0.975)

prop.med_URI <- round(pred.prop.med_URI, 3)
prop.lower_URI <- round(pred.prop.lower_URI, 3)
prop.upper_URI <- round(pred.prop.upper_URI, 3)



### Fever/LRI
mod_sim_LRI_logit <- glm(
  cbind(n_Fever_LRI,n_LRI-n_Fever_LRI) ~ Year + GDP + Temperature + Precipitation + mal.prop,
  family = binomial(link="logit"), data = prop.data)

coef_LRI <- coef(mod_sim_LRI_logit); vcov_LRI <- vcov(mod_sim_LRI_logit)

new.data.UNI <- UNI.data %>% dplyr::select(c("Year", "GDP", "Temperature", "Precipitation", "mal.prop"))
X <- cbind(1, as.matrix(new.data.UNI))  # Add intercept manually
n_obs <- nrow(X)
beta_sim <- MASS::mvrnorm(n = n_sim, mu = coef_LRI, Sigma = vcov_LRI)  
logit_p_sim <- X %*% t(beta_sim)
pred.prop.LRI <- plogis(logit_p_sim)

pred.prop.med_LRI <- apply(pred.prop.LRI, 1, median)
pred.prop.lower_LRI <- apply(pred.prop.LRI, 1, quantile, probs = 0.025)
pred.prop.upper_LRI <- apply(pred.prop.LRI, 1, quantile, probs = 0.975)

prop.med_LRI <- round(pred.prop.med_LRI, 3)
prop.lower_LRI <- round(pred.prop.lower_LRI, 3)
prop.upper_LRI <- round(pred.prop.upper_LRI, 3)

sum(residuals(mod_sim_LRI_logit, type = "pearson")^2) /
  df.residual(mod_sim_LRI_logit)
sum(residuals(mod_sim_URI_logit, type = "pearson")^2) /
  df.residual(mod_sim_URI_logit)



### Result--------------
par(mfrow = c(2, 1))

hist(unique(analysis_data$Fever_LRI), xlim = c(0, 1), ylim = c(0, 15), breaks = seq(0,1,by=0.05), main = "P(Fever|LRI)", xlab = "")
hist(unique(analysis_data$Fever_URI), xlim = c(0, 1), ylim = c(0, 15), breaks = seq(0,1,by=0.05), main = "P(Fever|URI)", xlab = "")

hist(prop.med_URI, xlim = c(0, 1), breaks = seq(0,1,by=0.05), ylim = c(0, 1000), main = "P(Fever|URI)", xlab = "")
hist(prop.med_LRI, xlim = c(0, 1), breaks = seq(0,1,by=0.05), ylim = c(0, 1000), main = "P(Fever|LRI)", xlab = "")

par(mfrow = c(1, 1))

prop.Fever.RI <- data.frame(
  Country = UNI.data$Country,
  Year = UNI.data$Year,
  Fever_LRI = pred.prop.LRI,
  Fever_URI = pred.prop.URI)

saveRDS(prop.Fever.RI, "outputs/prop.Fever.RI_logit.rds")
