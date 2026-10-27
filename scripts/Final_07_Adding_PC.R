rm(list = ls())
gc()

requiredPackages = c("dplyr", "readr", "data.table")
lib.path <- .libPaths()
for(p in requiredPackages){
  if(!require(p,character.only = TRUE)) install.packages(p, lib = lib.path)
  library(p,character.only = TRUE)
}

options(scipen = 999, digits = 4)
set.seed(123)
UNI <- readRDS("outputs/pc_gdp_ten.rds")

ab.data_LRI      <- readRDS("outputs/ab.data_LRI.rds")
ab.data_Fever    <- readRDS("outputs/ab.data_Fever_only.rds")
ab.data_all    <- readRDS("outputs/ab.data_all.rds")
analysis_data <- readRDS("outputs/analysis_data_prop.rds")

ab.data_LRI <- ab.data_LRI %>% dplyr::select(-GDP)
ab.data_Fever <- ab.data_Fever %>% dplyr::select(-GDP)
ab.data_all <- ab.data_all %>% dplyr::select(-GDP)

analysis_data <- analysis_data %>%
  dplyr::select(c("Antibiotics", "ARI", 
                  "ChildSex","ChildAgeMonth", 
                  "CountryYear", "Country", "Year", "AbLRI", "AbLRI_SE", "AbFever.nonARI", "AbFever.nonARI_SE", 
                  "n_AbLRI", "n_LRI", "n_Ab_Fever.nonARI", "n_Fever.nonARI"))

analysis_data <- na.omit(analysis_data)
analysis_data <- analysis_data %>% mutate(PersonTime = 1/26) 
analysis_data$Year <- as.integer(analysis_data$Year)

UNI_subset <- UNI
analysis_data$Country[analysis_data$Country == "C\x99te d'Ivoire"] <- "Cote d'Ivoire"
country_mapping <- c("Korea, Rep."="Republic of Korea", "United States"="United States of America", "Egypt, Arab Rep."="Egypt", "Slovak Republic"="Slovakia", "Turkiye"="Türkiye")
UNI_subset <- UNI_subset %>% mutate(Country = ifelse(Country %in% names(country_mapping), country_mapping[Country], Country))
analysis_data$Country[analysis_data$Country == "Egypt, Arab Rep."] <- "Egypt"
analysis_data$Year <- as.integer(analysis_data$Year)

# ab.data mainly does not have data in african countries
ab.data_LRI          <- ab.data_LRI           %>% inner_join(UNI_subset, by = c("Country","Year"))
ab.data_Fever        <- ab.data_Fever         %>% inner_join(UNI_subset, by = c("Country","Year"))
ab.data_all          <- ab.data_all           %>% inner_join(UNI_subset, by = c("Country","Year"))
AbLRI.data           <- analysis_data         %>% inner_join(UNI_subset, by = c("Country","Year"))

ab.data_LRI <- na.omit(ab.data_LRI)
ab.data_Fever <- na.omit(ab.data_Fever)
ab.data_all <- na.omit(ab.data_all)

AbLRI.data <- na.omit(AbLRI.data)

years <- sort(unique(UNI$Year))
colors <- rainbow(length(years))  
year_colors <- colors[as.numeric(factor(UNI$Year, levels = years))]
plot(UNI$PC1,UNI$PC2, col = UNI$Year)

AbLRI.data_PC <- AbLRI.data %>% dplyr::select(PC1, PC2) %>% dplyr::filter(!duplicated(.))
points(AbLRI.data_PC$PC1, AbLRI.data_PC$PC2, col = 'red', pch = 19)
points(ab.data_LRI$PC1, ab.data_LRI$PC2, col = 'blue', pch = 19)

legend("topright", legend = years, col = colors, pch = 19, cex = 0.6, title = "Year")

saveRDS(AbLRI.data, "outputs/AbLRI.data.UNI.rds")
saveRDS(ab.data_LRI, "outputs/ab.data_LRI.UNI.rds")
saveRDS(ab.data_Fever, "outputs/ab.data_Fever.UNI.rds")
saveRDS(ab.data_all, "outputs/ab.data_all.UNI.rds")
write.csv(UNI, "outputs/UNI.data.csv")