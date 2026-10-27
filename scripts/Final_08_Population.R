rm(list = ls())
gc()

requiredPackages = c("dplyr", "readr", "data.table")
lib.path <- .libPaths()
for(p in requiredPackages){
  if(!require(p,character.only = TRUE)) install.packages(p, lib = lib.path)
  library(p,character.only = TRUE)
}


Population <- read.csv("data/population_under_5.csv")
UNI <- read.csv("outputs/UNI.data.csv")
Population <- Population %>% dplyr::select(c(Entity,Year,Population))
colnames(Population) <- c("Country","Year","Population")

head(Population)
head(UNI)

a <- unique(UNI$Country)
b <- unique(Population$Country)

setdiff(a,b)
setdiff(b,a)

country_mapping <- c("Democratic Republic of Congo" ="Congo, Dem. Rep.", "Congo"="Congo, Rep.", "Egypt"="Egypt, Arab Rep.", "Gambia"="Gambia, The", "Iran"="Iran, Islamic Rep.", "Kyrgyzstan"="Kyrgyz Republic", "Laos"="Lao PDR", "Russia"="Russian Federation", "East Timor"="Timor-Leste", "Turkey"="Turkiye", "Vietnam"="Viet Nam")
Population <- Population %>% mutate(Country = ifelse(Country %in% names(country_mapping), country_mapping[Country], Country))

UNI <- UNI %>% inner_join(Population, by = c("Country","Year"))
head(UNI)
sum(is.na(UNI))

write.csv(UNI, "outputs/UNI.data.csv")
saveRDS(UNI, "outputs/pc_gdp.rds")
