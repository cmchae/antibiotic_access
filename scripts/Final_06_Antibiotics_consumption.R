rm(list = ls())
gc()

requiredPackages = c("dplyr", "tidyr")
lib.path <- .libPaths()
for(p in requiredPackages){
  if(!require(p,character.only = TRUE)) install.packages(p, lib = lib.path)
  library(p,character.only = TRUE)
}

options(scipen = 100, digits = 4)

n_sim <- 10000
ab.data <- read.csv("data/DDD_country_year_class.csv")
ab.data <- ab.data[,-3]

# matching country name
ab.data$Country <- trimws(ab.data$Country, which = "right")
country_mapping = c("Korea, Rep." = "Republic of Korea", "Czech Republic" = "Czechia", "Russia" = "Russian Federation", "United States" = "United States of America", "Venezuela" = "Venezuela (Bolivarian Republic of)", "Vietnam" = "Viet Nam")
ab.data <- ab.data %>% mutate(Country = ifelse(Country %in% names(country_mapping), country_mapping[Country], Country))

# Select included antibiotics
ab.data <- ab.data %>% group_by(Country, Year, Antimicrobial) %>% mutate(DDD_sum = sum(DDD, na.rm = TRUE))
ab.data <- ab.data[, c("Country", "Year", "Antimicrobial", "DDD_sum")]
ab.data <- ab.data %>% pivot_wider(names_from = Antimicrobial, values_from = DDD_sum)

colnames(ab.data) <- c("Country", "Year", "Aminoglycosides", "BSP", "Carbapenems",     
                       "Cephalosporins" , "Chloramphenicols" , "Fluoroquinolones" ,  "Glycopeptides"    , "Glycylcyclines"   ,          
                       "Lipopeptides"   , "Macrolides"       , "Monobactams"      ,  "NSP"              ,            "Others",          
                       "Oxazolidinones" , "Phosphonics"      , "Polymyxins"       ,  "Tetracyclines"    , "Trimethoprim"     ,          
                       "Rifampicin")

ab.data <- ab.data %>% dplyr::select(c(Country, Year, BSP, NSP, Cephalosporins, Macrolides, Aminoglycosides, Carbapenems, Glycopeptides, Trimethoprim))

U5_prop <- readRDS("outputs/U5_prop.rds")
ab.data <- ab.data %>% inner_join(U5_prop, by = c("Country", "Year"))
ab.data[,3:10] <- ab.data[,3:10]*ab.data$U5_prop
ab.data

nrow(ab.data)
ab.data <- na.omit(ab.data)
nrow(ab.data)

for(i in 3:10){
  print(min(ab.data[,i]))
}

ab.data <- ab.data %>% filter(Carbapenems != 0) %>% filter(Glycopeptides != 0) %>% filter(NSP != 0)
nrow(ab.data)

U5_matrix <- readRDS("outputs/matrix.result.ver.1.rds")
U5_matrix <- ab.data[,1:2] %>% inner_join(U5_matrix, by = c("Country", "Year"))
U5_LRI_matrix <- U5_matrix[,3:(n_sim+2)] + U5_matrix[,(2*n_sim+3):(3*n_sim++2)]
U5_fever_only_matrix <- U5_matrix[,(n_sim+3):(2*n_sim+2)]
U5_all_matrix <- U5_matrix[,3:(n_sim+2)] + U5_matrix[,(n_sim+3):(2*n_sim+2)] + U5_matrix[,(2*n_sim+3):(3*n_sim++2)]

n_abx <- 8
U5_DDD_LRI <- U5_DDD_fever_only <- U5_DDD_all <- matrix(NA, nrow = nrow(ab.data), ncol = n_sim * n_abx)
U5_LRI_matrix <- apply(U5_LRI_matrix, 2, as.numeric)
U5_fever_only_matrix <- apply(U5_fever_only_matrix, 2, as.numeric)
U5_all_matrix <- apply(U5_all_matrix, 2, as.numeric)

# Loop to compute
for (i in 1:n_abx) {
  abx_col <- as.numeric(ab.data[[i + 2]]) 
  for (j in 1:n_sim) {
    col_index <- (i - 1) * n_sim + j
    U5_DDD_LRI[, col_index] <- abx_col / U5_LRI_matrix[, j]
    U5_DDD_fever_only[, col_index] <- abx_col / U5_fever_only_matrix[, j]
    U5_DDD_all[, col_index] <- abx_col / U5_all_matrix[, j]
    }
}

abx_names <- colnames(ab.data)[3:10]
col_names <- unlist(lapply(abx_names, function(name) paste0(name, "_sim", 1:n_sim)))
colnames(U5_DDD_LRI) <- col_names
U5_DDD_LRI_df <- cbind(ab.data[, 1:2], as.data.frame(U5_DDD_LRI))
U5_DDD_fever_only_df <- cbind(ab.data[, 1:2], as.data.frame(U5_DDD_fever_only))
U5_DDD_all_df <- cbind(ab.data[, 1:2], as.data.frame(U5_DDD_all))

chunk_size <- 10000 
total_cols <- ncol(U5_DDD_LRI) 
n_chunks <- total_cols / chunk_size 

LRI_median_list <- fever_only_median_list <- all_median_list <- 
LRI_lower_list  <- fever_only_lower_list  <- all_lower_list  <- 
LRI_upper_list  <- fever_only_upper_list  <- all_upper_list  <- list()

chunk_size <- n_sim; total_cols <- ncol(U5_DDD_LRI); n_chunks <- total_cols / chunk_size

# Loop through each chunk
for (i in 1:n_chunks) {
  cols <- ((i - 1) * chunk_size + 1):(i * chunk_size)
  
  LRI_median_list[[i]] <- apply(U5_DDD_LRI[, cols], 1, median)
  fever_only_median_list[[i]] <- apply(U5_DDD_fever_only[, cols], 1, median)
  all_median_list[[i]] <- apply(U5_DDD_all[, cols], 1, median)
  
  LRI_lower_list[[i]] <- apply(U5_DDD_LRI[, cols], 1, quantile, probs = 0.025)
  fever_only_lower_list[[i]] <- apply(U5_DDD_fever_only[, cols], 1, quantile, probs = 0.025)
  all_lower_list[[i]] <- apply(U5_DDD_all[, cols], 1, quantile, probs = 0.025)
  
  LRI_upper_list[[i]] <- apply(U5_DDD_LRI[, cols], 1, quantile, probs = 0.975)
  fever_only_upper_list[[i]] <- apply(U5_DDD_fever_only[, cols], 1, quantile, probs = 0.975)
  all_upper_list[[i]] <- apply(U5_DDD_all[, cols], 1, quantile, probs = 0.975)
}

#assigning antibiotics name
ab_names <- c("BSP", "NSP", "Cephalosporins", "Macrolides", 
              "Aminoglycosides", "Carbapenems", "Glycopeptides", "Trimethoprim")
names(LRI_median_list) <- names(LRI_lower_list) <- names(LRI_upper_list) <- 
  names(fever_only_median_list) <- names(fever_only_lower_list) <- names(fever_only_upper_list) <- 
  names(all_median_list) <- names(all_lower_list) <- names(all_upper_list) <- ab_names

ab.data_all <- ab.data_Fever_only <- ab.data_LRI <- ab.data_CountryYear <- data.frame(
  Country = ab.data$Country,
  Year = ab.data$Year)

for (ab in ab_names) {
  ab.data_LRI[[paste0("U5_DDD_LRI_", ab)]]         <- LRI_median_list[[ab]]
  ab.data_LRI[[paste0("U5_DDD_LRI_", ab, "_lower")]] <- LRI_lower_list[[ab]]
  ab.data_LRI[[paste0("U5_DDD_LRI_", ab, "_upper")]] <- LRI_upper_list[[ab]]
  
  ab.data_Fever_only[[paste0("U5_DDD_fever_", ab)]]         <- fever_only_median_list[[ab]]
  ab.data_Fever_only[[paste0("U5_DDD_fever_", ab, "_lower")]] <- fever_only_lower_list[[ab]]
  ab.data_Fever_only[[paste0("U5_DDD_fever_", ab, "_upper")]] <- fever_only_upper_list[[ab]]

  ab.data_all[[paste0("U5_DDD_fever_", ab)]]         <- all_median_list[[ab]]
  ab.data_all[[paste0("U5_DDD_fever_", ab, "_lower")]] <- all_lower_list[[ab]]
  ab.data_all[[paste0("U5_DDD_fever_", ab, "_upper")]] <- all_upper_list[[ab]]
  
  print(i)
}

GDP_PPP <- read.csv("data/GDP_PPP.csv")
GDP_PPP <- GDP_PPP[,1:22]
colnames(GDP_PPP) <- c("Country","2001",  "2002",  "2003",  "2004",  "2005",  "2006",  "2007",  "2008",  "2009",  "2010",  "2011",  "2012","2013", "2014",  "2015",  "2016",  "2017",  "2018",  "2019",  "2020",  "2021")
GDP_PPP <- GDP_PPP %>% pivot_longer(cols = -Country, names_to = "Year", values_to = "GDP")
GDP_PPP <- na.omit(GDP_PPP)

country_mapping2 = c("Egypt" = "Egypt, Arab Rep.", "Republic of Korea" = "Korea, Rep.", "Slovakia" = "Slovak Republic", "Türkiye" ="Turkiye", "United States of America" = "United States")
GDP_PPP <- GDP_PPP %>% mutate(Country = ifelse(Country %in% names(country_mapping), country_mapping[Country], Country))
GDP_PPP$Year <- as.numeric(GDP_PPP$Year)
ab.data_LRI <- ab.data_LRI %>%left_join(GDP_PPP, by = c("Country","Year"))
ab.data_Fever_only <- ab.data_Fever_only %>%left_join(GDP_PPP, by = c("Country","Year"))
ab.data_all <- ab.data_all %>%left_join(GDP_PPP, by = c("Country","Year"))


U5_DDD_LRI       <- cbind(ab.data_CountryYear,U5_DDD_LRI)      
U5_DDD_Fever_only     <- cbind(ab.data_CountryYear,U5_DDD_fever_only) 
U5_DDD_all     <- cbind(ab.data_CountryYear,U5_DDD_all) 


saveRDS(U5_DDD_LRI, "outputs/U5_DDD_LRI_matrix.rds")
saveRDS(U5_DDD_Fever_only, "outputs/U5_DDD_Fever_only_matrix.rds")
saveRDS(U5_DDD_all, "outputs/U5_DDD_all_matrix.rds")

saveRDS(ab.data_LRI, "outputs/ab.data_LRI.rds")
saveRDS(ab.data_Fever_only, "outputs/ab.data_Fever_only.rds")
saveRDS(ab.data_all, "outputs/ab.data_all.rds")