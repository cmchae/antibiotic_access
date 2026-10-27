rm(list = ls())
gc()

requiredPackages = c("rdhs", "dplyr", "tidyr", "gt", "zoo", "haven", "data.table", "labelled")

lib.path <- .libPaths()
for(p in requiredPackages){
  if(!require(p,character.only = TRUE)) install.packages(p, lib = lib.path)
  library(p,character.only = TRUE)
}

options(scipen = 999, digits = 2)

datasets_dir <- "data/datasets/"
local_files <- list.files(datasets_dir, pattern = "\\.rds$", full.names = TRUE)

# Variables to keep
keep_variables <- c(
  "caseid", "v000", "v001", "v007", "v008", "bord", "b3", "b4", "b5", "b7", "b8", 
  "h22", "h31", "h31b", "h31c", "h37i", "h37j", "h37n", "h37o", "h37p",
  "v025", "v106", "v190", "v481",
  "h37a", "h37b", "h37c", "h37d", "h37da", "h37e", "h37aa", "h37ab", "h37f", "h37g", "h37h"
)

# Function to load and filter a dataset
load_and_filter <- function(file) {
  dat <- as.data.table(readRDS(file))
  dat <- dat[, intersect(keep_variables, names(dat)), with = FALSE]
  dat$filename <- tools::file_path_sans_ext(basename(file))
  dat
}

# Load all datasets into a named list
extract <- lapply(local_files, load_and_filter)
names(extract) <- sapply(local_files, function(x) tools::file_path_sans_ext(basename(x)))

for (i in 1:length(extract)) {
  year_val <- min(unique(extract[[i]]$v007), na.rm = TRUE)
  file_label <- names(extract)[i]
  if (year_val > 2019) {
    cat("i =", i, "filename =", file_label, "Year =", year_val, "\n")
  }
}

datasets_after_2019 <- c("BDKR81FL", "BFKR81FL", "CIKR81FL", "GHKR8AFL", "GHKR8BFL", "JOKR81FL", "KEKR8CFL", "KHKR82FL", "LSKR81FL", "MDKR81FL", "MZKR81FL", "PHKR82FL","SNKR8RFL","TZKR82FL")
extract <- extract[!names(extract) %in% datasets_after_2019]

# assigning correct survey year (Afghanistan and Nepal: assigning survey starting year)
extract$AFKR71FL$v007 <- 2015
extract$NPKR61FL$v007 <- 2011
extract$NPKR7HFL$v007 <- 2016
extract$NPKR82FL$v007 <- 2022

# excluding DHS ver.4
version.4 <- names(extract)[substr(names(extract), 5, 5) == "4"]
length(version.4)
extract <- extract[!names(extract) %in% version.4]

## datasets without variables
column_names <- lapply(extract, colnames); length(keep_variables); sum(lapply(column_names, length)<24) 
datasets_to_remove <- c("AMKR54FL", "COKR53FL","EGKR51FL","ETKR51FL", "GNKR52FL", "RWKR53FL") 
extract <- extract[!names(extract) %in% datasets_to_remove]

## datasets without respiratory symptoms information (h31~h31c)
columns_to_check <- c("h22", "h31", "h31b", "h31c")
column_names <- lapply(extract, colnames)
datasets_with_all_columns <- sapply(column_names, function(cols) {all(columns_to_check %in% cols)})
which(!datasets_with_all_columns)

# Exclude the datasets with only NA in antibiotic variables
datasets_all_na_h37i_h37j <- c()
for (dataset_name in names(extract)) {
  df <- extract[[dataset_name]]
  if (all(c("h37i", "h37j") %in% names(df))) {
    if (all(is.na(df$h37i)) && all(is.na(df$h37j))) {
      datasets_all_na_h37i_h37j <- c(datasets_all_na_h37i_h37j, dataset_name)
}}}
length(datasets_all_na_h37i_h37j)
datasets_with_all_NA <- c("AOKR71FL", "BDKR51FL", "BDKR72FL", "BDKR7RFL", "BDKR81FL", 
                        "BFKR62FL", "BJKR51FL", "CDKR51FL", "CDKR61FL", "CGKR51FL", 
                        "CGKR61FL", "CIKR62FL", "CMKR61FL", "COKR61FL", "COKR72FL", 
                        "ETKR81FL", "GAKR61FL", "GHKR5AFL", "GHKR72FL", "GMKR61FL", 
                        "GMKR81FL", "HTKR71FL", "IDKR51FL", "JOKR61FL", "KHKR51FL", 
                        "KMKR61FL", "LSKR61FL", "MDKR51FL", "MDKR81FL", "MLKR53FL", 
                        "MLKR7AFL", "MMKR71FL", "MRKR71FL", "MWKR61FL", "NGKR53FL", 
                        "NGKR6AFL", "NGKR7BFL", "NIKR51FL", "NIKR61FL", "PGKR71FL", 
                        "RWKR5AFL", "RWKR61FL", "RWKR70FL", "SLKR51FL", "SLKR61FL", 
                        "SLKR7AFL", "SNKR61FL", "SNKR6DFL", "SNKR71FL", "SNKR7AFL", 
                        "SNKR7IFL", "SNKR81FL", "SNKR8BFL", "STKR51FL", "TDKR71FL", 
                        "TGKR61FL", "TRKR51FL", "TRKR62FL", "TZKR63FL", "TZKR7BFL", 
                        "UAKR51FL", "UGKR7BFL", "ZMKR61FL", "ZWKR52FL", "ZWKR62FL", 
                        "ZWKR72FL")
extract <- extract[!names(extract) %in% datasets_with_all_NA]
length(extract)

# # Antibiotic demand (h37)
# Treated with antibiotics: h37i(pills), h37j(inection), h37(f,g,n,o,p: country specific antibiotic)
h37_cs <- c( #Country specific antibiotics: label 
  "Amoxicillin taken for fever", "Beta lactum taken for fever",
  "Amoxicillin taken for fever/cough","Cotramoxazole (CS) taken for fever/cough", 
  "CS other taken for fever/cough: Pen VK", "Macrolides for fever/cough", 
  "Cotrimoxazole taken for fever/cough", "CS other for fever/cough: Erithromycin", 
  "Bactrim taken for fever/cough", "Aminophyllin for fever/cough", 
  "Septrin taken for fever/cough", "Erythromycine(CS) taken for fever/cough",
  "Beta lactum Tab/Syrup taken for fever/cough", "Ampicilin taken for fever/cough", 
  "Amoxicillin(CS) taken for fever/cough", "Quinolone for fever/cough", 
  "Macrolides Tablets taken for fever/cough", "Quinolone Syrup taken for fever/cough", 
  "Other antibiotic taken for fever/cough", "Antibiotic", "Amoxycillin taken for fever/cough")

# finding antimalaria agent from country specific malarial agent
h37_mal_lab <- c("NA - CS antimalarial taken for fever", "Other antimalarial taken for fever", "CS antimalarial taken for fever/cough", "Other antimalarial taken for fever/cough", "NA - Other antimalarial taken for fever", "NA - CS antimalarial taken for fever/cough", "NA - Other antimalarial taken for fever/cough", "CS antimalarial taken for fever", "CS antimalarial taken for fever/cough (Delagil)", "CS antimalarial taken for fever/cough (Premakhin)", "Artemether-Lumefantrine taken for fever/cough", "Mefloquine (CS) taken for fever/cough", "Artesunate/artinate (CS) taken for fever/cough", "Primaquine antimalarial taken for fever/cough", "CS antimalarial taken for f-NA", "NA - NA-CS antimalarial taken for fever/cough", "malarine: antimalarial taken for fever/cough", "A+M 2,3,4: antimalarial taken for fever/cough", "Other antimalarial taken fo-NA", "Unknown antimalarial taken for fever/cough", "NA-CS antimalarial taken for fever", "NA-CS antimalarial taken for fever/cough")

h37_vars <- c("h37f", "h37g", "h37h", "h37n", "h37o", "h37p")
h37_mal <- c("h37f", "h37g", "h37h")

for (dataset_name in names(extract)) {
  df <- extract[[dataset_name]]
  for (var in h37_mal) {
    if (var %in% names(df)) {
      new_var <- paste0(var, "_mal")   
      df[[new_var]] <- df[[var]]       
      var_label(df[[new_var]]) <- paste0(var_label(df[[var]]), " (mal)") 
    }
  }
  extract[[dataset_name]] <- df
}

for (dataset_name in names(extract)) {
  df <- extract[[dataset_name]]
  for (var in h37_vars) {
    if (var %in% names(df)) {
      var_lbl <- var_label(df[[var]])
      if (is.null(var_lbl) || (!(var_lbl %in% h37_cs))) {
        df[[var]] <- NA }}}
  extract[[dataset_name]] <- df
}

h37_new_mal <- c("h37f_mal", "h37g_mal", "h37h_mal")

for (dataset_name in names(extract)) {
  df <- extract[[dataset_name]]
  for (var in h37_new_mal) {
    if (var %in% names(df)) {
      var_lbl <- var_label(df[[var]])
      if (is.null(var_lbl) || (!(var_lbl %in% h37_mal_lab))) {
        df[[var]] <- NA }}}
  extract[[dataset_name]] <- df
}

presence_matrix <-      
  matrix(0, nrow = length(extract), ncol = length(h37_vars), dimnames = list(names(extract), h37_vars))
for (dataset_name in names(extract)) {
  df <- extract[[dataset_name]]
  for (var in h37_vars) {
    if (var %in% names(df) && any(!is.na(df[[var]]))) {
      presence_matrix[dataset_name, var] <- 1
}}}
as.data.frame(presence_matrix)


# Combining the datasets into one dataframe and align the variables with order by A-Z
combined_data <- rbindlist(extract, fill = TRUE, ignore.attr = TRUE)
combined_data <- combined_data %>% dplyr::select(order(names(.))) 


# making overall antibiotic treatment variable
# have at least 1 -> 1, only have 0 and NA -> 0, only have NA -> NA  
determine_value <- function(x) {                     
  if (any(x == 1, na.rm = TRUE)) {return(1)}         
  else if (any(x == 0, na.rm = TRUE)) {return(0)}    
  else {return(NA)}}                                 

h37g_labels <- purrr::imap_dfr(
  extract,
  ~ tibble(
    filename = .y,
    h37g_present = "h37g" %in% names(.x),
    h37g_label = if ("h37g" %in% names(.x)) {
      as.character(labelled::var_label(.x$h37g))
    } else {
      NA_character_}))
print(h37g_labels, n = Inf)

combined_data$Ab_treatment <- 
  apply(combined_data[, c("h37i", "h37j", "h37g", "h37n", "h37o", "h37p")], 1, determine_value)
combined_data$mal_treatment <- 
  apply(combined_data[, c("h37a", "h37b", "h37c", "h37d", "h37e", "h37f", "h37da", "h37aa", "h37ab", "h37f_mal", "h37g_mal", "h37h_mal")], 1, determine_value) 
combined_data <- combined_data %>% dplyr::select(-starts_with("h37"))

# DHS survey version from filename
combined_data <- combined_data %>% mutate(survey_v = substr(filename, 5,5))

# Country name
country_name <- c(AL5 = "Albania", AF7 = "Afghanistan", BO5 = "Bolivia", NP6 = "Nepal", NP7 = "Nepal", NP8 = "Nepal", 
                     AL7 = "Albania", AM6 = "Armenia", AM7 = "Armenia", AZ5 = "Azerbaijan", BD5 = "Bangladesh", 
                     BD6 = "Bangladesh", BD7 = "Bangladesh", BD8 = "Bangladesh", BF8 = "Burkina Faso", BJ6 = "Benin", 
                     BJ7 = "Benin", CI8 = "Cote d'Ivoire", CM7 = "Cameroon", EG5 = "Egypt", EG6 = "Egypt", ET6 = "Ethiopia", 
                     ET7 = "Ethiopia", GA7 = "Gabon", GH8 = "Ghana", GN6 = "Guinea", GN7 = "Guinea", GY5 = "Guyana", 
                     HN5 = "Honduras", HN6 = "Honduras", HT5 = "Haiti", HT6 = "Haiti",  ID6 = "Indonesia", ID7 = "Indonesia",  
                     JO5 = "Jordan", JO6 = "Jordan", JO7 = "Jordan", JO8 = "Jordan",  KE5 = "Kenya", KE6 = "Kenya", KE8 = "Kenya",  
                     KH5 = "Cambodia", KH6 = "Cambodia", KH8 = "Cambodia",  LS6 = "Lesotho", LS8 = "Lesotho",  ML6 = "Mali",
                     MV5 = "Maldives", MV7 = "Maldives",  MW7 = "Malawi",  MZ6 = "Mozambique", MZ8 = "Mozambique",
                     NG5 = "Nigeria",  PE5 = "Peru", PE6 = "Peru",  PH5 = "Philippines", PH6 = "Philippines", 
                     PH7 = "Philippines", PH8 = "Philippines",  PK5 = "Pakistan", PK6 = "Pakistan", PK7 = "Pakistan",
                     RW7 = "Rwanda",  SN7 = "Senegal", SN8 = "Senegal",  SZ5 = "Eswatini",  TJ6 = "Tajikistan", 
                     TJ7 = "Tajikistan",  TZ8 = "Tanzania",  UG5 = "Uganda", UG6 = "Uganda",  YE6 = "Yemen",
                     ZA7 = "South Africa",  ZM5 = "Zambia", ZM7 = "Zambia",  ZW5 = "Zimbabwe")


combined_data$Country <- country_name[combined_data$v000]
combined_data_unique <- combined_data[!duplicated(combined_data[, c("filename", "Country")]), c("filename", "Country")]
length(unique(combined_data$Country)); length(unique(combined_data$v000))
sum(with(combined_data, table(Country, v007))>0) 


# Mutate "I don't know (6, 8, 9)" value in h22 (fever) and h31 (respiratory symptoms) to NA
combined_data[, c("h22","h31","h31b","h31c")] <- lapply(combined_data[, c("h22","h31","h31b","h31c")], function(x) ifelse(x  %in% c(6, 8, 9), NA, x))
combined_data$h31 <- ifelse(combined_data$h31 == 2, 1, combined_data$h31)
combined_data <- combined_data %>%
  rename(sy_fever = h22, sy_cough = h31, sy_rapid_breath = h31b, sy_cpain_rnose = h31c)


## b section
combined_data %>% dplyr::select(starts_with("b")) %>% colnames()
b_columns <- combined_data %>% dplyr::select(starts_with("b"))
b_column_tables <-
  lapply(b_columns, function(column) {table(column, useNA = "always")})
b_column_tables
combined_data <- combined_data %>% 
  rename(ch_bth_order = bord, ch_bth_month = b3, ch_sex = b4, ch_alive = b5, ch_dth_month = b7, ch_age_yr = b8)
colnames(combined_data)


## v section
combined_data %>% dplyr::select(starts_with("v")) %>% colnames()
v_columns <- combined_data %>% dplyr::select(starts_with("v"))
v_column_tables <- 
  lapply(v_columns, function(column) {table(column, useNA = "always")})
v_column_tables
combined_data <- 
  combined_data %>% rename(
  survey_code = v000,
  cluster_code = v001,
  Year = v007,
  int_month = v008)


# household section
combined_data <- #mutate household status names
  combined_data %>%  rename(
  hh_location = v025,   
  hh_educ = v106, 
  hh_wealth = v190,     
  hh_insurance = v481)  
hh_columns <- combined_data %>% dplyr::select(starts_with("hh"))
hh_column_tables <- lapply(hh_columns, function(column) {table(column, useNA = "always")})
hh_column_tables
combined_data[, c("hh_educ","hh_insurance")] <-
  lapply(combined_data[, c("hh_educ","hh_insurance")], function(x) ifelse(x  %in% c(8, 9), NA, x))


# check the duplication
combined_data <- 
  combined_data %>% mutate(ch_id = paste(survey_code, caseid, ch_bth_order, sep = "_"))
any(duplicated(combined_data$ch_id)) 
combined_data <- combined_data %>% dplyr::select(-c("caseid", "ch_bth_order"))

#age in month
combined_data$ch_age_month <- combined_data$int_month - combined_data$ch_bth_month

## adding GBD regions
gbd_region <- read.csv("data/gbd_region.csv")
gbd_region[gbd_region$Country == "Côte d'Ivoire", ] <- "Cote d'Ivoire"
gbd_region[gbd_region$Country == "Arab Republic of Egypt", ] <- "Egypt"
gbd_region[gbd_region$Country == "Republic of Yemen", ] <- "Yemen"
combined_data <- combined_data %>% left_join(gbd_region, by = "Country")

##selecting dataset
analysis_data <- combined_data %>%
  dplyr::select(c("Country", "Year", "ch_id", "ch_sex", "ch_age_month","ch_alive","ch_dth_month", "survey_code", "survey_v", "gbd_region", "hh_location",	"cluster_code", "int_month", "hh_educ",	"hh_wealth",	"hh_insurance", "sy_fever",	"sy_cough",	"sy_rapid_breath",	"sy_cpain_rnose", "Ab_treatment", "mal_treatment"))

# Remove rows where sy_cough or sy_fever is NA
with(analysis_data, table(sy_cough,sy_fever, useNA = 'always'))
analysis_data <- analysis_data[!is.na(analysis_data$sy_cough) & !is.na(analysis_data$sy_fever), ]

# Check the distribution of symptoms
with(analysis_data, table(sy_cough, useNA='always')); with(analysis_data, table(sy_fever, useNA='always'))
with(analysis_data, table(sy_rapid_breath, useNA='always')); with(analysis_data, table(sy_cpain_rnose, useNA='always'))
with(analysis_data, table(sy_cough, sy_rapid_breath, useNA='always'))
with(analysis_data, table(sy_cough, sy_cpain_rnose, useNA='always'))
with(analysis_data, table(sy_rapid_breath, sy_cpain_rnose, useNA='always'))
with(analysis_data, table(Ab_treatment, useNA = 'always'))
with(analysis_data, table(sy_fever, Ab_treatment, mal_treatment), margin = 1)

# ARI variable
analysis_data$ARI <- NA
# Step 1: Assign 3 if rapid breath and chest pain
analysis_data$ARI[(analysis_data$sy_cpain_rnose %in% c(1,3) & analysis_data$sy_rapid_breath == 1)] <- 3
# Step 2: Assign 2 if rapid breath without chest pain
analysis_data$ARI[is.na(analysis_data$ARI) & (analysis_data$sy_rapid_breath == 1)] <- 2
# Step 3: Assign 1 if only with fever or cough
analysis_data$ARI[is.na(analysis_data$ARI) & (analysis_data$sy_cou == 1 | analysis_data$sy_fev == 1)] <- 1
# Step 4: Assign 0 only if without fever and cough
analysis_data$ARI[is.na(analysis_data$ARI) & (analysis_data$sy_cou == 0 | analysis_data$sy_fev == 0)] <- 0
with(analysis_data, table(ARI, useNA = 'always'))
with(analysis_data, table(ARI, sy_cough, useNA = 'always'))
with(analysis_data, table(ARI, sy_fever, useNA = 'always'))
with(analysis_data, table(sy_rapid_breath, sy_fever, useNA = 'always'))
with(analysis_data, table(ARI, Ab_treatment, useNA = 'always'))
with(analysis_data, table(sy_rapid_breath, sy_cough, survey_v, useNA = 'always'))

# ARI ~ Fever table by survey version
df <- as.data.frame(with(analysis_data, table(sy_fever, ARI, survey_v)))
clean_df <- df %>%  pivot_wider(names_from = ARI, values_from = Freq, values_fill = 0) %>% arrange(survey_v, sy_fever) %>%
  rename(`Fever` = sy_fever, `Survey Version` = survey_v,`ARI 0` = `0`,`ARI 1` = `1`,`ARI 2` = `2`,`ARI 3` = `3`)
table_fever <- clean_df %>%  gt(groupname_col = "Survey Version") %>% tab_header(title = "ARI by Fever Status and Survey Round") %>%
  cols_label( `Fever` = "Fever Status", `ARI 0` = "ARI = 0", `ARI 1` = "ARI = 1", `ARI 2` = "ARI = 2", `ARI 3` = "ARI = 3")
table_fever

# ARI ~ Cough table by survey version
df <- as.data.frame(with(analysis_data, table(sy_cough, ARI, survey_v)))
clean_df <- df %>%  pivot_wider(names_from = ARI, values_from = Freq, values_fill = 0) %>% arrange(survey_v, sy_cough) %>%
  rename(`Cough` = sy_cough, `Survey Version` = survey_v,`ARI 0` = `0`,`ARI 1` = `1`,`ARI 2` = `2`,`ARI 3` = `3`)
table_cough <- clean_df %>%  gt(groupname_col = "Survey Version") %>% tab_header(title = "ARI by Cough Status and Survey Round") %>%
  cols_label( `Cough` = "Cough Status", `ARI 0` = "ARI = 0", `ARI 1` = "ARI = 1", `ARI 2` = "ARI = 2", `ARI 3` = "ARI = 3")
table_cough

# Cough ~ Fever table by survey version
df <- as.data.frame(with(analysis_data, table(sy_fever, sy_cough, survey_v)))
clean_df <- df %>%  pivot_wider(names_from = sy_cough, values_from = Freq, values_fill = 0) %>% arrange(survey_v, sy_fever) %>%
  rename(`Fever` = sy_fever, `Survey Round` = survey_v,`sy_cough 0` = `0`,`sy_cough 1` = `1`)
table_fever_cough <- clean_df %>%  gt(groupname_col = "Survey Round") %>% tab_header(title = "Cough by Fever Status and Survey Round") %>%
  cols_label( `Fever` = "Fever Status", `sy_cough 0` = "Cough = 0", `sy_cough 1` = "Cough = 1")
table_fever_cough

# make Country_Year variable
with(analysis_data, table(Year, Country, useNA ='always'))
analysis_data$CountryYear <- paste(analysis_data$Country, analysis_data$Year, sep = "_")
unique(analysis_data$CountryYear)

# merging and excluding the CountryYear with very small number of Antibiotics treatment and ARI data.
with(analysis_data,length(unique(Country))); with(analysis_data,length(unique(CountryYear))) 
with(analysis_data, table(Ab_treatment, ARI, CountryYear))

with(analysis_data[analysis_data$CountryYear == "Albania_2008"], table(Ab_treatment, ARI, CountryYear))
with(analysis_data[analysis_data$CountryYear == "Albania_2018"], table(Ab_treatment, ARI, CountryYear))
with(analysis_data[analysis_data$CountryYear == "Cameroon_2019"], table(Ab_treatment, ARI, CountryYear))
with(analysis_data[analysis_data$CountryYear == "Eswatini_2007"], table(Ab_treatment, ARI, CountryYear))
with(analysis_data[analysis_data$CountryYear == "Jordan_2018"], table(Ab_treatment, ARI, CountryYear))
with(analysis_data[analysis_data$CountryYear == "Mali_2013"], table(Ab_treatment, ARI, CountryYear))


analysis_data$CountryYear[analysis_data$CountryYear == "Albania_2008"] <- "Albania_2009" 
analysis_data$CountryYear[analysis_data$CountryYear == "Albania_2018"] <- "Albania_2017" 
analysis_data$CountryYear[analysis_data$CountryYear == "Cameroon_2019"] <- "Cameroon_2018" 
analysis_data$CountryYear[analysis_data$CountryYear == "Eswatini_2007"] <- "Eswatini_2006" 
analysis_data$CountryYear[analysis_data$CountryYear == "Gabon_2020"] <- "Gabon_2019" 
analysis_data$CountryYear[analysis_data$CountryYear == "Gabon_2021"] <- "Gabon_2019" 
analysis_data$CountryYear[analysis_data$CountryYear == "Jordan_2018"] <- "Jordan_2017" 
analysis_data$CountryYear[analysis_data$CountryYear == "Mali_2013"] <- "Mali_2012"
analysis_data$CountryYear[analysis_data$CountryYear == "Rwanda_2020"] <- "Rwanda_2019" 

analysis_data$Year <- as.numeric(sub(".*_(\\d{4})$", "\\1", analysis_data$CountryYear)) 

analysis_data <- analysis_data[!analysis_data$CountryYear %in% 
                                 c("Armenia_2015", "Armenia_2016", "Cambodia_2011", 
                                   "Maldives_2009", "Maldives_2016", "Maldives_2017", "Peru_2003", 
                                   "Peru_2004", "Peru_2005", "Tajikistan_2017", "Zambia_2019"), ]

with(analysis_data,length(unique(Country))); with(analysis_data,length(unique(CountryYear))) 
nrow(distinct(analysis_data[,c("Country","Year")]))


# add GDP variable
analysis_data[analysis_data$Country == "Egypt",]$Country <- "Egypt, Arab Rep."
analysis_data[analysis_data$Country == "Yemen",]$Country <- "Yemen, Rep."

GDP.PPP <- read.csv("data/GDP_PPP.csv")
colnames(GDP.PPP)[-1] <- gsub("^X", "", colnames(GDP.PPP)[-1])
GDP.PPP <- GDP.PPP %>% pivot_longer(cols = -Country, names_to = "Year", values_to = "GDP")
GDP.PPP$Year <- as.numeric(GDP.PPP$Year)
a <- GDP.PPP$Country
b <- analysis_data$Country
setdiff(b,a)
analysis_data <- analysis_data %>% inner_join(GDP.PPP, by = c("Country" = "Country", "Year" = "Year"))
unique(analysis_data$GDP)

colnames(analysis_data) <-c("Country", "Year", "ChildId", "ChildSex", "ChildAgeMonth", "ChildAlive", "ChildDeathMonth", "SurveyId", "SurveyVersion", "GBDRegion", "UrbanRural", "Cluster", "InterviewMonth", "Education", "Wealth", "Insurance", "Fever", "Cough", "RapidBreath", "ChestpainRunnynose", "Antibiotics", "Antimalaria", "ARI", "CountryYear", "GDP")

analysis_data
dir.create("outputs", recursive = TRUE, showWarnings = FALSE)
saveRDS(analysis_data, file ="outputs/analysis_data.rds")

rm(list = ls())
gc()
