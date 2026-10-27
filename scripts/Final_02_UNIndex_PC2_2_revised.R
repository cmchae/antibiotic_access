rm(list = ls())
gc()

requiredPackages = c("data.table", "dplyr", "tidyr", "tibble", "missMDA", "loadings")
lib.path <- .libPaths()
for(p in requiredPackages){
  if(!require(p,character.only = TRUE)) install.packages(p, lib = lib.path)
  library(p,character.only = TRUE)
}

options(scipen = 999)

# -----------------prepare UNIndex dataset -------------------
# merging datsets
df.1 <- fread("data/UNIndex/1to2.csv"  , select = c("Series Name","Country Name", "2001 [YR2001]", "2002 [YR2002]"), col.names = c("Series", "Country", "2001", "2002"))
df.2 <- fread("data/UNIndex/3to5.csv"  , select = c("Series Name","Country Name", "2003 [YR2003]", "2004 [YR2004]", "2005 [YR2005]"), col.names = c("Series", "Country", "2003", "2004", "2005"))
df.3 <- fread("data/UNIndex/6to8.csv"  , select = c("Series Name","Country Name", "2006 [YR2006]", "2007 [YR2007]", "2008 [YR2008]"), col.names = c("Series", "Country", "2006", "2007", "2008"))
df.4 <- fread("data/UNIndex/9to11.csv" , select = c("Series Name","Country Name", "2009 [YR2009]", "2010 [YR2010]", "2011 [YR2011]"), col.names = c("Series", "Country", "2009", "2010", "2011"))
df.5 <- fread("data/UNIndex/12to13.csv", select = c("Series Name","Country Name", "2012 [YR2012]", "2013 [YR2013]"), col.names = c("Series", "Country", "2012", "2013"))
df.6 <- fread("data/UNIndex/14to16.csv", select = c("Series Name","Country Name", "2014 [YR2014]", "2015 [YR2015]", "2016 [YR2016]"), col.names = c("Series", "Country", "2014", "2015", "2016"))
df.7 <- fread("data/UNIndex/17to19.csv", select = c("Series Name","Country Name", "2017 [YR2017]", "2018 [YR2018]", "2019 [YR2019]"), col.names = c("Series", "Country", "2017", "2018", "2019"))
df.8 <- fread("data/UNIndex/20to22.csv", select = c("Series Name","Country Name", "2020 [YR2020]", "2021 [YR2021]", "2022 [YR2022]"), col.names = c("Series", "Country", "2020", "2021", "2022"))

df_list <- list(df.1, df.2, df.3, df.4, df.5, df.6, df.7, df.8)
df_list <- lapply(df_list, function(df) df[1:(.N - 5)])

first_two_cols <- lapply(df_list, function(df) df[, .(Series, Country)])
sapply(2:8, function(i) all.equal(first_two_cols[[1]], first_two_cols[[i]])) # check whether the Country series are same.
dfs.years <- lapply(df_list, function(df) df[, -c("Series", "Country"), with = FALSE])
df.combined <- cbind(df_list[[1]][,1:2], do.call(cbind, dfs.years))
dim(df.combined); colnames(df.combined)
df.long <- df.combined %>% pivot_longer(cols = -c(Series, Country), names_to = "Year", values_to = "Value")
df.long <- df.long %>% filter(Series != "")
df.ready <- df.long %>% pivot_wider(names_from = Series, values_from = Value)
df.ready <- df.ready[, -((ncol(df.ready)-1):ncol(df.ready))] #(5874,330)
dim(df.ready) # (5874,328)
unique(df.ready$Country)

non_countries <- c("Africa Eastern and Southern", "Africa Western and Central", "Arab World",
                   "Caribbean small states", "Central Europe and the Baltics", "East Asia & Pacific",
                   "East Asia & Pacific (excluding high income)", "East Asia & Pacific (IDA & IBRD countries)",
                   "Early-demographic dividend", "Euro area", "Europe & Central Asia",
                   "Europe & Central Asia (excluding high income)", "Europe & Central Asia (IDA & IBRD countries)",
                   "European Union", "Fragile and conflict affected situations",
                   "Heavily indebted poor countries (HIPC)", "High income", "IBRD only",
                   "IDA & IBRD total", "IDA blend", "IDA only", "IDA total",
                   "Late-demographic dividend", "Latin America & Caribbean",
                   "Latin America & Caribbean (excluding high income)",
                   "Latin America & the Caribbean (IDA & IBRD countries)",
                   "Least developed countries: UN classification", "Low & middle income", "Low income",
                   "Lower middle income", "Middle East & North Africa",
                   "Middle East & North Africa (excluding high income)",
                   "Middle East & North Africa (IDA & IBRD countries)", "Middle income",
                   "North America", "Not classified", "OECD members", "Other small states",
                   "Pacific island small states", "Post-demographic dividend", "Pre-demographic dividend",
                   "Small states", "South Asia", "South Asia (IDA & IBRD)", "Sub-Saharan Africa",
                   "Sub-Saharan Africa (excluding high income)", "Sub-Saharan Africa (IDA & IBRD countries)",
                   "Upper middle income", "World")
df.ready <- df.ready %>% filter(!(Country %in% non_countries))  # exclude non countries
df.ready <- df.ready %>% filter(df.ready$Year < 2020)
df.ready <- df.ready %>% filter(df.ready$Year > 2002)

Income <- read.csv("data/Income_level.csv")
colnames(Income)[-1] <- gsub("^X", "", colnames(Income)[-1])
Income <- Income %>% pivot_longer(cols = 2:24, names_to = "Year", values_to = "level")  # Convert wide month 
country_mapping <- c("C\x99te d'Ivoire" = "Cote d'Ivoire", "T\x9frkiye" = "Turkiye", "Korea, Dem. Rep." = "Korea, Dem. People's Rep.")
Income <- Income %>% mutate(Country = ifelse(Country %in% names(country_mapping), country_mapping[Country], Country))
colnames(Income)[3] <- "Income" 
Income$Year <- as.numeric(Income$Year)
df.ready$Year <- as.numeric(df.ready$Year)

df.ready <- df.ready %>% inner_join(Income, by = c("Country" = "Country", "Year" = "Year"))
df.ready <- df.ready %>% filter(Income %in% c("L", "LM", "UM", "R"))
df.ready[df.ready$Income == "R",]$Income <- "H"

a1 <- as.data.frame(with(df.ready, table(Country, Year)))
b1 <- a1 %>% group_by(Country) %>% summarise(Country_sum = mean(Freq), .groups = "drop")
c1 <- b1[b1$Country_sum > 0.3,]$Country
df.ready <- df.ready %>% filter(Country %in% c1)
with(df.ready, table(Country, Year))
df.ready

a <- df.ready[,c(1,2,38)]
colnames(a) <- c("Country", "Year", "Current_health_expenditure_per_capita_PPP")
saveRDS(a, "outputs/health_expenditure.rds")


# adding GDP variable 
GDP.PPP <- read.csv("data/GDP_PPP.csv") # Yemen and Venezuela from IMF data (GDP per capita international dollar)
colnames(GDP.PPP)[-1] <- gsub("^X", "", colnames(GDP.PPP)[-1])
GDP.PPP <- GDP.PPP %>% pivot_longer(cols = -Country, names_to = "Year", values_to = "GDP")
GDP.PPP$Year <- as.numeric(GDP.PPP$Year)
df.ready.gdp <- df.ready %>% inner_join(GDP.PPP, by = c("Country" = "Country", "Year" = "Year"))
unique(df.ready.gdp[!complete.cases(df.ready.gdp), ]$Country)
df.ready.gdp <- df.ready.gdp[!is.na(df.ready.gdp$GDP),]

meta_cols <- c("Country", "Year", "GDP", "Income")
df.pca <- df.ready.gdp %>% mutate(across(-all_of(meta_cols), ~ na_if(., "..")   )) %>%
  mutate(across(-all_of(meta_cols), ~ as.numeric(.), .names = "{.col}")) 

par(mfrow = c(1,1))
hist(round(colSums(is.na(df.pca))/nrow(df.pca), digits = 3),)
# Overall missingness
na_prop <- colSums(is.na(df.pca)) / nrow(df.pca)
df.pca.ready <- df.pca[, na_prop <= 0.3]
dim(df.pca.ready) #   <0.3:(2344,187)
df.numeric <- df.pca.ready[, !(names(df.pca.ready) %in% c("Country", "Year", "GDP", "Income"))]

# Multiple imputation using PCA
set.seed(23)
mi.pca <- MIPCA(df.numeric, ncp = 5, nboot = 5)
imputed.list <- mi.pca$res.MI

length(imputed.list)
sapply(imputed.list, nrow)
sapply(imputed.list, ncol)

# Convert completed datasets to numeric matrices
imputed.list <- lapply(imputed.list, function(x) {
  x <- as.matrix(x)
  storage.mode(x) <- "double"
  x
})

stopifnot(all(sapply(imputed.list, nrow) == nrow(imputed.list[[1]])))
stopifnot(all(sapply(imputed.list, ncol) == ncol(imputed.list[[1]])))
stopifnot(all(sapply(imputed.list, function(x) identical(colnames(x), colnames(imputed.list[[1]])))))


# Mean of the five completed datasets
imputed.array <- array(
  NA_real_,
  dim = c(
    nrow(imputed.list[[1]]),
    ncol(imputed.list[[1]]),
    length(imputed.list)
  ),
  dimnames = list(
    NULL,
    colnames(imputed.list[[1]]),
    paste0("Imp", seq_along(imputed.list))
  )
)

for (m in seq_along(imputed.list)) {
  imputed.array[, , m] <- imputed.list[[m]]
}

imputed.mean <- apply(imputed.array, c(1, 2), mean)
storage.mode(imputed.mean) <- "double"



# Define a common PCA coordinate system using the mean completed dataset


pca.common <- prcomp(
  imputed.mean,
  center = TRUE,
  scale. = TRUE
)

summary(pca.common)
screeplot(
  pca.common,
  main = "Screeplot from mean-imputation reference PCA",
  type = "lines"
)

# Retain PC1-PC5
common_rotation <- pca.common$rotation[, 1:5, drop = FALSE]
colnames(common_rotation) <- paste0("PC", 1:5)

common_center <- pca.common$center
common_scale <- pca.common$scale



# Match PC signs to PCA from the first completed dataset
# Sign has no statistical consequence but maintains comparability


pca.first <- prcomp(
  imputed.list[[1]],
  center = TRUE,
  scale. = TRUE
)

first_rotation <- pca.first$rotation[, 1:5, drop = FALSE]

for (j in 1:5) {
  sgn_j <- sign(
    sum(common_rotation[, j] * first_rotation[, j], na.rm = TRUE)
  )
  
  if (is.na(sgn_j) || sgn_j == 0) {
    sgn_j <- 1
  }
  
  common_rotation[, j] <- common_rotation[, j] * sgn_j
}

round(crossprod(first_rotation, common_rotation), 3)



# Project each completed dataset onto the SAME PCA axes


meta_pca <- df.pca.ready[, c("Country", "Year", "GDP", "Income")] %>%
  mutate(.row = row_number())

pc_scores.list <- vector("list", length(imputed.list))
df.pc.mi.list <- vector("list", length(imputed.list))

for (m in seq_along(imputed.list)) {
  
  dat_m <- imputed.list[[m]]
  dat_m_standardized <- sweep(dat_m, 2, common_center, FUN = "-")
  dat_m_standardized <- sweep(dat_m_standardized, 2, common_scale, FUN = "/")
  
  pc_scores_m <- dat_m_standardized %*% common_rotation
  colnames(pc_scores_m) <- paste0("PC", 1:5)
  
  pc_scores.list[[m]] <- pc_scores_m
  
  df.pc.mi.list[[m]] <- bind_cols(
    tibble(.imp = rep(m, nrow(meta_pca))),
    meta_pca,
    as_tibble(pc_scores_m)
  )
}

df.pc.mi <- bind_rows(df.pc.mi.list)

# Average PC1-PC5 across the five completed datasets
df.pc.mean <- df.pc.mi %>%
  group_by(.row, Country, Year, GDP, Income) %>%
  summarise(
    PC1 = mean(PC1),
    PC2 = mean(PC2),
    PC3 = mean(PC3),
    PC4 = mean(PC4),
    PC5 = mean(PC5),
    .groups = "drop"
  ) %>%
  arrange(.row)

df.pc <- df.pc.mean %>%
  dplyr::select(-.row)

# Validation:
# Averaged projected scores should equal scores obtained directly from
# the mean completed dataset
imputed.mean.standardized <- sweep(imputed.mean, 2, common_center, FUN = "-")
imputed.mean.standardized <- sweep(imputed.mean.standardized, 2, common_scale, FUN = "/")

pc_scores.mean.direct <- imputed.mean.standardized %*% common_rotation
colnames(pc_scores.mean.direct) <- paste0("PC", 1:5)

max(abs(as.matrix(df.pc[, paste0("PC", 1:5)]) - pc_scores.mean.direct))

# Basic checks
unique(df.pc$Country)
head(df.pc)
str(df.pc)

par(mfrow = c(1, 2))
hist(df.pc$PC1, main = "PC1", xlab = "PC1")
hist(df.pc$PC2, main = "PC2", xlab = "PC2")
par(mfrow = c(1, 1))

# PC scores over calendar year
years <- sort(unique(df.pc$Year))
colors <- rainbow(length(years))
year_colors <- colors[as.numeric(factor(df.pc$Year, levels = years))]

par(mfrow = c(1, 2))
plot(
  df.pc$PC1,
  df.pc$Year,
  col = year_colors,
  xlab = "PC1",
  ylab = "Year"
)
plot(
  df.pc$PC2,
  df.pc$Year,
  col = year_colors,
  xlab = "PC2",
  ylab = "Year"
)
par(mfrow = c(1, 1))



# Common-axis loadings, correlations, and nominal p-values


stopifnot(nrow(imputed.mean) == nrow(df.pc.mean))
stopifnot(identical(colnames(imputed.mean), rownames(common_rotation)))
stopifnot(identical(df.pc.mean$.row, seq_len(nrow(imputed.mean))))

rvals <- matrix(
  NA_real_,
  nrow = ncol(imputed.mean),
  ncol = 2
)

pvals <- matrix(
  NA_real_,
  nrow = ncol(imputed.mean),
  ncol = 2
)

rownames(rvals) <- rownames(pvals) <- colnames(imputed.mean)
colnames(rvals) <- colnames(pvals) <- c("PC1", "PC2")

for (i in 1:ncol(imputed.mean)) {
  for (j in 1:2) {
    
    test <- cor.test(
      imputed.mean[, i],
      df.pc.mean[[paste0("PC", j)]],
      method = "pearson"
    )
    
    rvals[i, j] <- unname(test$estimate)
    pvals[i, j] <- test$p.value
  }
}

loadings <- common_rotation

result_all <- data.frame(
  Variable = rownames(loadings),
  Loading_PC1 = loadings[, 1],
  Loading_PC2 = loadings[, 2],
  Correlation_PC1 = rvals[, 1],
  Correlation_PC2 = rvals[, 2],
  Pvalue_PC1 = pvals[, 1],
  Pvalue_PC2 = pvals[, 2]
)

result_PC1 <- result_all %>%
  transmute(
    Variable,
    Loading_PC1 = round(Loading_PC1, 3),
    Correlation_PC1 = round(Correlation_PC1, 3),
    Pvalue_PC1 = format.pval(
      Pvalue_PC1,
      digits = 4,
      eps = 0.0001
    )
  ) %>%
  arrange(desc(abs(Loading_PC1)))

result_PC2 <- result_all %>%
  transmute(
    Variable,
    Loading_PC2 = round(Loading_PC2, 3),
    Correlation_PC2 = round(Correlation_PC2, 3),
    Pvalue_PC2 = format.pval(
      Pvalue_PC2,
      digits = 4,
      eps = 0.0001
    )
  ) %>%
  arrange(desc(abs(Loading_PC2)))

write.csv(result_PC1, "outputs/result_PC1_mean_reference.csv", row.names = FALSE)
write.csv(result_PC2, "outputs/result_PC2_mean_reference.csv", row.names = FALSE)


# small island countries ~500,000 population, civil war

exclude_countries_small_island <- c("Cabo Verde", "Dominica", "Grenada", "Kiribati", "Maldives", "Marshall Islands", "Micronesia, Fed. Sts.", "Palau", "Samoa", "Seychelles", "St. Kitts and Nevis", "St. Lucia", "St. Vincent and the Grenadines", "Tonga", "Tuvalu", "Vanuatu")
df.pc <- df.pc %>% filter(!Country %in% exclude_countries_small_island)

exclude_countries_war <- c("Kosovo", "Libya", "Somalia", "Syrian Arab Republic", "West Bank and Gaza")
df.pc <- df.pc %>% filter(!Country %in% exclude_countries_war)

exclude_countries_last_3yr_high_income <- c("Chile", "Croatia", "Hungary", "Latvia", "Lithuania", "Panama", "Poland", "Uruguay")
df.pc <- df.pc %>% filter(!Country %in% exclude_countries_last_3yr_high_income)

exclude_countries_outboundary <- c("Montenegro","Zimbabwe","Equatorial Guinea")
df.pc <- df.pc %>% filter(!Country %in% exclude_countries_outboundary)

plot(df.pc$Year, df.pc$PC1, col = df.pc$Year)
df.pc[which.min(df.pc$PC1),]


#--------- adding malarial endemic countries -----------

mal.endemic.df <- read.csv("data/WHO2024_malarial_endemic.csv")


country_mapping <- c(
  "Comoros " = "Comoros",
  "Congo" = "Congo, Rep.",
  "C\x99te d\xd5Ivoire" = "Cote d'Ivoire",
  "Democratic Republic of the Congo" = "Congo, Dem. Rep.",
  "Gambia" = "Gambia, The",
  "Niger " = "Niger",
  "United Republic of Tanzania" = "Tanzania",
  "Iran (Islamic Republic of)" = "Iran, Islamic Rep.",
  "Bolivia (Plurinational State of)" = "Bolivia",
  "Venezuela (Bolivarian Republic of)" = "Venezuela",
  "Democratic People\xd5s Republic of Korea" = "Korea, Dem. Rep.",
  "Lao People\xd5s Democratic Republic " = "Lao PDR",
  "Republic of Korea" = "Korea, Rep."
)

mal.endemic.df <- mal.endemic.df %>% mutate(Country = ifelse(Country %in% names(country_mapping), country_mapping[Country], Country))
df.pc <- df.pc %>% left_join(mal.endemic.df, by = c("Country" = "Country"))
df.pc <- df.pc %>% mutate(Endemic = ifelse(is.na(Endemic), 0, Endemic))


# ------- adding GBD region -------
gbd_region <- read.csv("data/gbd_region.csv")
country_mapping <- c("Arab Republic of Egypt" = "Egypt, Arab Rep.", "Islamic Republic of Iran" = "Iran, Islamic Rep.", "Turkey" = "Turkiye", "Lao People's Democratic Republic" = "Lao PDR", "Vietnam" = "Viet Nam", "Cape Verde" = "Cabo Verde", "Côte d'Ivoire" = "Cote d'Ivoire", "The Gambia" = "Gambia, The", "Congo" = "Congo, Rep.", "Democratic Republic of Congo" = "Congo, Dem. Rep.", "Federated States of Micronesia" = "Micronesia, Fed. Sts.", "Saint Kitts and Nevis" = "St. Kitts and Nevis", "Saint Lucia" = "St. Lucia", "Saint Vincent and the Grenadines" = "St. Vincent and the Grenadines")
gbd_region <- gbd_region %>% mutate(Country = ifelse(Country %in% names(country_mapping), country_mapping[Country], Country))
df.pc <- df.pc %>% inner_join(gbd_region, by = c("Country"))


# adding population --------------
Population <- read.csv("data/population_under_5.csv") %>%
  dplyr::select(Entity, Year, Population)
colnames(Population) <- c("Country", "Year", "Population")
head(Population)
head(df.pc)

country_mapping <- c(
  "Democratic Republic of Congo" = "Congo, Dem. Rep.",
  "Congo" = "Congo, Rep.",
  "Egypt" = "Egypt, Arab Rep.",
  "Gambia" = "Gambia, The",
  "Iran" = "Iran, Islamic Rep.",
  "Kyrgyzstan" = "Kyrgyz Republic",
  "Laos" = "Lao PDR",
  "Russia" = "Russian Federation",
  "East Timor" = "Timor-Leste",
  "Turkey" = "Turkiye",
  "Vietnam" = "Viet Nam",
  "Yemen" = "Yemen, Rep."
)

Population <- Population %>%
  mutate(Country = ifelse(Country %in% names(country_mapping), country_mapping[Country], Country))

df.pc <- df.pc %>%
  inner_join(Population, by = c("Country", "Year"))

head(df.pc)
sum(is.na(df.pc))

# Save one common PC dataset for all downstream models
# Final extrapolation models use PC1-PC2 only.
saveRDS(df.pc, "outputs/pc_gdp_ten.rds") 
saveRDS(df.pc, "outputs/pc_gdp_mean_reference.rds")

eigvals <- pca.common$sdev^2
propvar <- eigvals / sum(eigvals)
pcs <- 1:10
propvar <- propvar[1:10]

# ============================================================
# Export source data for Figure S1
# ============================================================

dir.create( "outputs/plot_objects", recursive = TRUE, showWarnings = FALSE)

Figure_S1_source_data <- data.frame(
  Principal_component = pcs,
  Eigenvalue = eigvals[1:10],
  Proportion_variance_explained = propvar
)

write.csv(
  Figure_S1_source_data,
  "outputs/plot_objects/Figure_S1_source_data.csv",
  row.names = FALSE
)

### supplemental plot

dir.create(
  "figures",
  recursive = TRUE,
  showWarnings = FALSE
)

pdf("figures/Figure_S1.pdf", width = 7.5, height = 4)

par(mar = c(2.5, 2.5, 1, 1))
par(mgp = c(3, 0.3, 0))
par(tck = -0.025)
par(lwd = 1)

# Eigenvalues and proportion of variance
xmin <- 1
xmax <- 10
ymin <- 0
ymax <- 0.3

plot(pcs, propvar, type = "n",
     xlab = "", ylab = "",
     xlim = c(xmin, xmax), ylim = c(ymin, ymax),
     axes = FALSE, ann = FALSE, main = "")

axis(1, at = pcs, labels = FALSE, pos = ymin)
text(x = pcs,
     y = ymin - 0.05*(ymax - ymin),
     labels = pcs,
     xpd = TRUE, cex = 0.65)

axis(2, at = seq(0, ymax, length.out = 4), labels = FALSE, pos = xmin)
text(y = seq(0, ymax, length.out = 4),
     x = xmin - 0.02*(xmax - xmin),
     labels = round(seq(0, ymax, length.out = 4), 2),
     srt = 0, adj = 1, xpd = TRUE, cex = 0.65)

lines(pcs, propvar, col = "#0072B2", lwd = 1.5)
points(pcs, propvar, pch = 16, cex = 0.7, col = "#0072B2")

mtext("Principal components", side = 1, line = 1.5, cex = 0.8)
mtext("Proportion of variance explained", side = 2, line = 1.5, cex = 0.8)

dev.off()