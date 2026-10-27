# ------------------8 antibiotic classes: proportion of LRTI antibiotic consumption ------------------
#### 0. prep ####
rm(list = ls())
gc()

requiredPackages = c("dplyr", "tidyr", "MASS")
lib.path <- .libPaths()
for(p in requiredPackages){
  if(!require(p,character.only = TRUE)) install.packages(p, lib = lib.path)
  library(p,character.only = TRUE)
}

options(scipen = 100, digits = 6)

n_sim <- 10000
set.seed(123)

years <- 2003:2019

Ab.names <- c(
  "BSP", "NSP", "Cephalosporins", "Macrolides",
  "Aminoglycosides", "Carbapenems", "Glycopeptides", "Trimethoprim"
)

Ab.colnames_LRI <- paste0("U5_DDD_LRI_", Ab.names)

income_levels <- c("L", "LM", "UM")

gbd_region_levels <- c(
  "Southeast Asia, East Asia, and Oceania",
  "Central Europe, Eastern Europe, and Central Asia",
  "Latin America and Caribbean",
  "South Asia",
  "North Africa and Middle East",
  "Sub-Saharan Africa"
)

dir.create("outputs/antibiotic_class", recursive = TRUE, showWarnings = FALSE)
dir.create("outputs/plot_df", recursive = TRUE, showWarnings = FALSE)

# ==============================================================================
# PART 1. Refit current PC1 + PC2 DDD models for all 8 classes
# ==============================================================================

# TRUE first time; FALSE after current 8-class coefficient RDS is generated
rerun_8class_DDD_models <- TRUE

if(rerun_8class_DDD_models){
  
  ab.data_LRI <- readRDS("outputs/ab.data_LRI.UNI.rds")
  matrix.LRI <- readRDS("outputs/U5_DDD_LRI_matrix.rds")
  
  stopifnot(
    all(c("Country", "Year", "PC1", "PC2", Ab.colnames_LRI) %in% names(ab.data_LRI))
  )
  
  matrix.meta <- ab.data_LRI %>% dplyr::select(Country, Year)
  matrix.LRI <- matrix.meta %>% inner_join(matrix.LRI, by = c("Country", "Year"))
  
  stopifnot(
    nrow(matrix.LRI) == nrow(ab.data_LRI),
    ncol(matrix.LRI) == 2 + 8 * n_sim
  )
  
  LRI.DDD_sim <- matrix.LRI %>%
    dplyr::select(-Country, -Year) %>%
    as.data.frame()
  
  coef_matrix_LRI_DDD <- matrix(
    NA_real_,
    nrow = n_sim,
    ncol = 8 * 3
  )
  
  colnames(coef_matrix_LRI_DDD) <- unlist(
    lapply(
      Ab.colnames_LRI,
      function(x) paste0(x, "_", c("(Intercept)", "PC1", "PC2"))
    )
  )
  
  for (j in 1:8) {
    
    for (i in 1:n_sim) {
      
      ab.data_LRI$LRI.DDD_sim <- LRI.DDD_sim[, n_sim * (j - 1) + i]
      
      mod <- glm(
        log(LRI.DDD_sim) ~ PC1 + PC2,
        data = ab.data_LRI
      )
      
      coef_mean <- coef(mod)
      coef_vcov <- vcov(mod)
      
      coef_matrix_LRI_DDD[
        i,
        3 * (j - 1) + (1:3)
      ] <- MASS::mvrnorm(
        1,
        mu = coef_mean,
        Sigma = coef_vcov
      )
      
      if (i %% 100 == 0) cat("Completed simulation:", Ab.names[j], i, "\n")
    }
  }
  
  coef_df_LRI_DDD <- as.data.frame(coef_matrix_LRI_DDD)
  
  stopifnot(
    nrow(coef_df_LRI_DDD) == n_sim,
    ncol(coef_df_LRI_DDD) == 24,
    !any(grepl("PC3", names(coef_df_LRI_DDD)))
  )
  
  saveRDS(
    coef_df_LRI_DDD,
    "outputs/antibiotic_class/coef_df_LRI_DDD_by_8class.rds"
  )
  
  rm(
    ab.data_LRI,
    matrix.LRI,
    matrix.meta,
    LRI.DDD_sim,
    coef_matrix_LRI_DDD,
    coef_df_LRI_DDD,
    mod
  )
  
  gc()
}


# ==============================================================================
# PART 2. Generate proportion of antibiotic consumption by 8 classes
# ==============================================================================

UNI_indexed <- readRDS("outputs/UNI_indexed.rds")
coef_df_LRI_DDD <- readRDS("outputs/antibiotic_class/coef_df_LRI_DDD_by_8class.rds")
U5_matrix <- readRDS("outputs/matrix.result.ver.2.rds")

UNI_indexed <- UNI_indexed %>% filter(Year %in% years)

stopifnot(
  all(c("Country", "Year", "Income", "gbd_region", "PC1", "PC2") %in% names(UNI_indexed)),
  nrow(coef_df_LRI_DDD) == n_sim,
  ncol(coef_df_LRI_DDD) == 24,
  ncol(U5_matrix) >= 2 + n_sim,
  !any(duplicated(UNI_indexed[,c("Country","Year")]))
)


# ---------------- define income group from modal income level in 2015-2019 ----------------
income_dist <- UNI_indexed %>%
  dplyr::select(Country, Year, Income) %>%
  pivot_wider(names_from = Year, values_from = Income)

last5_years <- as.character(2015:2019)

income_dist_summary <- income_dist %>%
  rowwise() %>%
  mutate(
    Last5 = {
      last_vals <- c_across(all_of(last5_years))
      tab <- table(last_vals)
      if(length(tab) == 0) NA_character_ else names(sort(tab, decreasing = TRUE))[1]
    }
  ) %>%
  ungroup()

UNI_indexed <- UNI_indexed %>%
  left_join(
    income_dist_summary %>% dplyr::select(Country, Last5),
    by = "Country"
  )

stopifnot(!any(is.na(UNI_indexed$Last5)))


# ---------------- align current LRTI incidence matrix ----------------

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
  "Yemen, Rep." = "Yemen"
)

Country_match <- as.character(UNI_indexed$Country)
map_index <- Country_match %in% names(country_mapping)
Country_match[map_index] <- unname(country_mapping[Country_match[map_index]])

key_UNI <- paste(Country_match, UNI_indexed$Year, sep = "___")
key_U5 <- paste(as.character(U5_matrix$Country), U5_matrix$Year, sep = "___")

U5_index <- match(key_UNI, key_U5)

if(any(is.na(U5_index))){
  print(
    data.frame(
      Country_UNI = UNI_indexed$Country[is.na(U5_index)],
      Country_match = Country_match[is.na(U5_index)],
      Year = UNI_indexed$Year[is.na(U5_index)]
    )
  )
  stop("Some UNI country-years could not be matched to outputs/matrix.result.ver.2.rds.")
}

U5_matrix <- U5_matrix[U5_index,,drop = FALSE]

Incidence_LRI_matrix <- as.matrix(
  U5_matrix[,3:(n_sim+2),drop = FALSE]
)

storage.mode(Incidence_LRI_matrix) <- "double"

stopifnot(
  nrow(Incidence_LRI_matrix) == nrow(UNI_indexed),
  ncol(Incidence_LRI_matrix) == n_sim
)


# ---------------- predict DDD per LRTI and calculate total consumption ----------------

X_pred <- model.matrix(~ PC1 + PC2, data = UNI_indexed)

stopifnot(
  identical(colnames(X_pred), c("(Intercept)", "PC1", "PC2"))
)

Total_DDD_LRI_class <- vector("list", length(Ab.names))
names(Total_DDD_LRI_class) <- Ab.names

for(i in 1:8){
  
  coef_cols <- 3 * (i - 1) + (1:3)
  
  LRI_DDD_sim <- exp(
    as.matrix(X_pred) %*%
      t(as.matrix(coef_df_LRI_DDD[,coef_cols,drop = FALSE]))
  )
  
  Total_DDD_LRI_class[[i]] <- Incidence_LRI_matrix * LRI_DDD_sim
  
  rm(LRI_DDD_sim)
  gc()
}

rm(U5_matrix, Incidence_LRI_matrix)


# ---------------- aggregate by income/GBD region and calculate class shares ----------------

make_share_summary <- function(group_var, group_levels){
  
  summary_list <- list()
  draw_list <- list()
  row_i <- 1
  
  for(group_i in group_levels){
    
    for(year_i in years){
      
      idx <- which(
        as.character(UNI_indexed[[group_var]]) == group_i &
          UNI_indexed$Year == year_i
      )
      
      if(length(idx) == 0) next
      
      class_total <- matrix(
        NA_real_,
        nrow = n_sim,
        ncol = 8,
        dimnames = list(NULL, Ab.names)
      )
      
      for(j in 1:8){
        class_total[,j] <- colSums(
          Total_DDD_LRI_class[[j]][idx,,drop = FALSE],
          na.rm = TRUE
        )
      }
      
      total_all_classes <- rowSums(class_total, na.rm = TRUE)
      
      class_share <- sweep(
        class_total,
        1,
        total_all_classes,
        FUN = "/"
      )
      
      class_share[!is.finite(class_share)] <- NA_real_
      
      result_row <- data.frame(
        Group = group_i,
        Year = year_i,
        stringsAsFactors = FALSE
      )
      
      for(j in 1:8){
        
        ab_i <- Ab.names[j]
        
        result_row[[paste0(ab_i, "_Median")]] <- median(
          class_share[,j],
          na.rm = TRUE
        )
        
        result_row[[paste0(ab_i, "_Q025")]] <- quantile(
          class_share[,j],
          0.025,
          na.rm = TRUE,
          names = FALSE
        )
        
        result_row[[paste0(ab_i, "_Q975")]] <- quantile(
          class_share[,j],
          0.975,
          na.rm = TRUE,
          names = FALSE
        )
      }
      
      summary_list[[row_i]] <- result_row
      draw_list[[paste(group_i, year_i, sep = "__")]] <- class_share
      
      row_i <- row_i + 1
    }
  }
  
  summary_result <- bind_rows(summary_list)
  names(summary_result)[1] <- group_var
  
  list(
    summary = summary_result,
    draws = draw_list
  )
}


income_result <- make_share_summary("Last5", income_levels)
gbd_result <- make_share_summary("gbd_region", gbd_region_levels)

LRI_income_original <- income_result$summary
LRI_gbd_region_original <- gbd_result$summary

LRI_income_draws <- income_result$draws
LRI_gbd_region_draws <- gbd_result$draws


# ---------------- order rows ----------------

LRI_income_original <- LRI_income_original %>%
  mutate(Last5 = factor(Last5, levels = income_levels)) %>%
  arrange(Last5, Year) %>%
  mutate(Last5 = as.character(Last5))

LRI_gbd_region_original <- LRI_gbd_region_original %>%
  mutate(gbd_region = factor(gbd_region, levels = gbd_region_levels)) %>%
  arrange(gbd_region, Year) %>%
  mutate(gbd_region = as.character(gbd_region))


# ---------------- save ORIGINAL proportions ----------------
saveRDS(
  LRI_income_original,
  "outputs/plot_df/Figure5_8class_income.rds"
)

saveRDS(
  LRI_gbd_region_original,
  "outputs/plot_df/Figure5_8class_gbd_region.rds"
)

write.csv(
  LRI_income_original,
  "outputs/plot_df/Figure5_8class_income.csv",
  row.names = FALSE
)

write.csv(
  LRI_gbd_region_original,
  "outputs/plot_df/Figure5_8class_gbd_region.csv",
  row.names = FALSE
)


# ---------------- log10-transformed proportions ----------------
LRI_income <- LRI_income_original %>%
  mutate(across(matches("_(Median|Q025|Q975)$"), log10))

LRI_gbd_region <- LRI_gbd_region_original %>%
  mutate(across(matches("_(Median|Q025|Q975)$"), log10))

saveRDS(
  LRI_income,
  "outputs/plot_df/Figure5_8class_income_log10.rds"
)

saveRDS(
  LRI_gbd_region,
  "outputs/plot_df/Figure5_8class_gbd_region_log10.rds"
)

write.csv(
  LRI_income,
  "outputs/plot_df/Figure5_8class_income_log10.csv",
  row.names = FALSE
)

write.csv(
  LRI_gbd_region,
  "outputs/plot_df/Figure5_8class_gbd_region_log10.csv",
  row.names = FALSE
)


# ---------------- objects compatible with your previous plot code ----------------

LRI_income_L  <- LRI_income %>% filter(Last5 == "L")
LRI_income_LM <- LRI_income %>% filter(Last5 == "LM")
LRI_income_UM <- LRI_income %>% filter(Last5 == "UM")

LRI_gbd_region_1 <- LRI_gbd_region %>%
  filter(gbd_region == "North Africa and Middle East")

LRI_gbd_region_2 <- LRI_gbd_region %>%
  filter(gbd_region == "Central Europe, Eastern Europe, and Central Asia")

LRI_gbd_region_3 <- LRI_gbd_region %>%
  filter(gbd_region == "Sub-Saharan Africa")

LRI_gbd_region_5 <- LRI_gbd_region %>%
  filter(gbd_region == "South Asia")

LRI_gbd_region_6 <- LRI_gbd_region %>%
  filter(gbd_region == "Latin America and Caribbean")

LRI_gbd_region_7 <- LRI_gbd_region %>%
  filter(gbd_region == "Southeast Asia, East Asia, and Oceania")


# ---------------- final checks ----------------

head(LRI_income_original)
head(LRI_income)

head(LRI_gbd_region_original)
head(LRI_gbd_region)

dim(LRI_income_original)
dim(LRI_gbd_region_original)



#### Figure S3: 8 antibiotic classes using current main-analysis outputs ####

rm(list = ls())
gc()

requiredPackages = c("dplyr")
lib.path <- .libPaths()
for(p in requiredPackages){
  if(!require(p,character.only = TRUE)) install.packages(p, lib = lib.path)
  library(p,character.only = TRUE)
}

options(scipen = 100, digits = 6)

dir.create("figures", recursive = TRUE, showWarnings = FALSE)

# Read the log10-transformed proportion summaries generated by the upstream 8-class code.
LRI_income <- readRDS("outputs/plot_df/Figure5_8class_income_log10.rds")
LRI_gbd_region <- readRDS("outputs/plot_df/Figure5_8class_gbd_region_log10.rds")

stopifnot(
  all(c("Last5", "Year") %in% names(LRI_income)),
  all(c("gbd_region", "Year") %in% names(LRI_gbd_region))
)

# Income groups
LRI_income_L  <- LRI_income %>% filter(Last5 == "L")
LRI_income_LM <- LRI_income %>% filter(Last5 == "LM")
LRI_income_UM <- LRI_income %>% filter(Last5 == "UM")

# GBD regions: retain numbering used in the previous Figure S3 code
LRI_gbd_region_1 <- LRI_gbd_region %>%
  filter(gbd_region == "North Africa and Middle East")

LRI_gbd_region_2 <- LRI_gbd_region %>%
  filter(gbd_region == "Central Europe, Eastern Europe, and Central Asia")

LRI_gbd_region_3 <- LRI_gbd_region %>%
  filter(gbd_region == "Sub-Saharan Africa")

LRI_gbd_region_5 <- LRI_gbd_region %>%
  filter(gbd_region == "South Asia")

LRI_gbd_region_6 <- LRI_gbd_region %>%
  filter(gbd_region == "Latin America and Caribbean")

LRI_gbd_region_7 <- LRI_gbd_region %>%
  filter(gbd_region == "Southeast Asia, East Asia, and Oceania")



# ============================================================
# Export source data for Figure S3
# ============================================================

dir.create(
  "outputs/plot_objects",
  recursive = TRUE,
  showWarnings = FALSE
)

# Antibiotic classes and corresponding panels
abx_classes <- c(
  "BSP",
  "Aminoglycosides",
  "Cephalosporins",
  "Macrolides",
  "Carbapenems",
  "Glycopeptides",
  "NSP",
  "Trimethoprim"
)

income_panels <- setNames(
  LETTERS[1:8],
  abx_classes
)

gbd_panels <- setNames(
  LETTERS[9:16],
  abx_classes
)


# ------------------------------------------------------------
# Panels A-H: proportion of antibiotic consumption by income
# ------------------------------------------------------------

income_labels <- c(
  "L"  = "Low income",
  "LM" = "Lower middle income",
  "UM" = "Upper middle income"
)

Figure_S3_income <- bind_rows(
  lapply(abx_classes, function(abx) {
    
    LRI_income %>%
      transmute(
        Panel = income_panels[[abx]],
        Antibiotic_class = abx,
        Stratification = "Income",
        Group = income_labels[Last5],
        Year = Year,
        
        Log10_Median = .data[[paste0(abx, "_Median")]],
        Log10_Lower = .data[[paste0(abx, "_Q025")]],
        Log10_Upper = .data[[paste0(abx, "_Q975")]],
        
        Median = 10^.data[[paste0(abx, "_Median")]],
        Lower = 10^.data[[paste0(abx, "_Q025")]],
        Upper = 10^.data[[paste0(abx, "_Q975")]]
      )
  })
)


# ------------------------------------------------------------
# Panels I-P: proportion of antibiotic consumption by GBD region
# ------------------------------------------------------------

Figure_S3_GBD <- bind_rows(
  lapply(abx_classes, function(abx) {
    
    LRI_gbd_region %>%
      transmute(
        Panel = gbd_panels[[abx]],
        Antibiotic_class = abx,
        Stratification = "GBD region",
        Group = gbd_region,
        Year = Year,
        
        Log10_Median = .data[[paste0(abx, "_Median")]],
        Log10_Lower = .data[[paste0(abx, "_Q025")]],
        Log10_Upper = .data[[paste0(abx, "_Q975")]],
        
        Median = 10^.data[[paste0(abx, "_Median")]],
        Lower = 10^.data[[paste0(abx, "_Q025")]],
        Upper = 10^.data[[paste0(abx, "_Q975")]]
      )
  })
)


# ------------------------------------------------------------
# Combine all panels
# ------------------------------------------------------------

Figure_S3_source_data <- bind_rows(
  Figure_S3_income,
  Figure_S3_GBD
) %>%
  arrange(Panel, Group, Year)


# ------------------------------------------------------------
# Save
# ------------------------------------------------------------

write.csv(
  Figure_S3_source_data,
  "outputs/plot_objects/Figure_S3_source_data.csv",
  row.names = FALSE
)

pdf("figures/Figure_S3.pdf", width = 7.5, height = 7)

layout(matrix(c(1,  1,  1,  1,  1,  1,
                0,  3,  4,  5,  6,  11,
                0,  7,  8,  9,  10, 11,
                0,  0,  0,  0,  0,  0,
                2,  2,  2,  2,  2,  2,
                0,  12, 13, 14, 15, 20, 
                0,  16, 17, 18, 19, 20), nrow = 7, byrow = TRUE), widths = c(0.15 ,1, 1, 1, 1, 1.05), heights = c(0.1,1,1,0.1,0.1,1,1))


#### 1. (title 1) ####

par(mar = c(0, 0, 0, 0))
plot.new()
text(x = -0.035, y = 0.5, labels = "Proportion of antibiotic consumption by class and income", 
     adj = 0, cex = 0.7, font = 2, xpd = TRUE)


#### 2. (title 2) ####

par(mar = c(0, 0, 0, 0))
plot.new()
text(x = -0.035, y = 0.5, labels = "Proportion of antibiotic consumption by class and GBD region", 
     adj = 0, cex = 0.7, font = 2, xpd = TRUE)


#### 3. (Fig.A) #####
par(mar = c(2.5, 2.5, 2, 1))  
par(mgp = c(3, 0.3, 0))    
par(tck = -0.025)
par(lwd = 1)


xmin <- 2003; xmax <- 2019
ymin <- log10(0.1) ; ymax <- log10(1)


plot(LRI_income$BSP_Median ~ LRI_income$Year, 
     data =  LRI_income,
     xlab = "", 
     ylab = "",
     pch = 16, col = rgb(0.8, 0.8, 0.8, 0),
     xlim = c(xmin,xmax),
     ylim = c(ymin,ymax),
     lwd = "1",
     axes=F,ann=F,
     main = "")



yticks <- log10(c(0.1,1))
ylabels <- c(0.1,1)    
minor_labels <- as.vector(sapply(0.1, function(v) v * 1:9))

axis(1, at = seq(xmin, xmax, by = 4), labels = FALSE, las = 2, pos = ymin)
text(x = seq(xmin, xmax, by = 4), y = ymin-(ymax-ymin)*0.08,  
     labels = format(round(seq(xmin, xmax, by = (xmax-xmin)/6*0+4), digits = 0), nsmall = 0), xpd = TRUE, cex = 0.65)

axis(2,  at = yticks, labels = FALSE, las = 2, pos = xmin)
axis(2, at = log10(minor_labels), labels = FALSE, las = 2, pos = xmin, tck = -0.015, lwd = 0.7)
text(y = yticks, x = xmin-0.04*(xmax-xmin),  
     labels = parse(text = paste0("10^", log10(ylabels))), srt = 0, adj = 1, xpd = TRUE, cex = 0.65)
text(y = log10(c(0.2, 0.3)),
     x = xmin - 0.04 * (xmax - xmin),
     labels = c(expression(2 %*% 10^-1), 
                expression(3 %*% 10^-1)),
     srt = 0, adj = 1, xpd = TRUE, cex = 0.5)  

lines(BSP_Median ~ Year, data = LRI_income_L,  col = "#CC79A7", lwd = 1.5)
lines(BSP_Median ~ Year, data = LRI_income_LM, col = "#009E73", lwd = 1.5)
lines(BSP_Median ~ Year, data = LRI_income_UM, col = "#0072B2", lwd = 1.5)

# Add shaded area for confidence interval
polygon(c(LRI_income_L$Year, rev(LRI_income_L$Year)),
        c(LRI_income_L$BSP_Q025, rev(LRI_income_L$BSP_Q975)),
        col = rgb(204/255, 121/255, 167/255, 0.25), border = NA)

polygon(c(LRI_income_LM$Year, rev(LRI_income_LM$Year)),
        c(LRI_income_LM$BSP_Q025, rev(LRI_income_LM$BSP_Q975)),
        col = rgb(0/255, 158/255, 115/255, 0.25), border = NA)

polygon(c(LRI_income_UM$Year, rev(LRI_income_UM$Year)),
        c(LRI_income_UM$BSP_Q025, rev(LRI_income_UM$BSP_Q975)),
        col = rgb(0/255, 114/255, 178/255, 0.25), border = NA)
mtext(expression("Proportion of antibiotic \n consumption by class*"), side = 2, line = 2, cex = 0.5)  
mtext("A. BSP", at = xmin - 0.54*(xmax-xmin), side = 3, line = 0.5, adj = 0, cex = 0.5, font = 1)


#### 4. (Fig.B) #####
xmin <- 2003; xmax <- 2019
ymin <- -3 ; ymax <- -1


plot(LRI_income$Aminoglycosides_Median ~ LRI_income$Year, 
     data =  LRI_income,
     xlab = "", 
     ylab = "",
     pch = 16, col = rgb(0.8, 0.8, 0.8, 0),
     xlim = c(xmin,xmax),
     ylim = c(ymin,ymax),
     lwd = "1",
     axes=F,ann=F,
     main = "")

# Example: log10 scale
yticks <- log10(c(1e-3, 1e-2, 1e-1))   
ylabels <- c(1e-3, 1e-2, 1e-1)         
minor_labels <- as.vector(sapply(c(1e-3, 1e-2), function(v) v * 2:9))


axis(1, at = seq(xmin, xmax, by = 4), labels = FALSE, las = 2, pos = ymin)
text(x = seq(xmin, xmax, by = 4), y = ymin-(ymax-ymin)*0.08,  
     labels = format(round(seq(xmin, xmax, by = (xmax-xmin)/6*0+4), digits = 0), nsmall = 0), xpd = TRUE, cex = 0.65)

axis(2,  at = yticks, labels = FALSE, las = 2, pos = xmin)
axis(2, at = log10(minor_labels), labels = FALSE, las = 2, pos = xmin, tck = -0.015, lwd = 0.7)
text(y = yticks, x = xmin-0.04*(xmax-xmin),  
     labels = parse(text = paste0("10^", log10(ylabels))), srt = 0, adj = 1, xpd = TRUE, cex = 0.65)


lines(Aminoglycosides_Median ~ Year, data = LRI_income_L,  col = "#CC79A7", lwd = 1.5)
lines(Aminoglycosides_Median ~ Year, data = LRI_income_LM, col = "#009E73", lwd = 1.5)
lines(Aminoglycosides_Median ~ Year, data = LRI_income_UM, col = "#0072B2", lwd = 1.5)

# Add shaded area for confidence interval
polygon(c(LRI_income_L$Year, rev(LRI_income_L$Year)),
        c(LRI_income_L$Aminoglycosides_Q025, rev(LRI_income_L$Aminoglycosides_Q975)),
        col = rgb(204/255, 121/255, 167/255, 0.25), border = NA)

polygon(c(LRI_income_LM$Year, rev(LRI_income_LM$Year)),
        c(LRI_income_LM$Aminoglycosides_Q025, rev(LRI_income_LM$Aminoglycosides_Q975)),
        col = rgb(0/255, 158/255, 115/255, 0.25), border = NA)

polygon(c(LRI_income_UM$Year, rev(LRI_income_UM$Year)),
        c(LRI_income_UM$Aminoglycosides_Q025, rev(LRI_income_UM$Aminoglycosides_Q975)),
        col = rgb(0/255, 114/255, 178/255, 0.25), border = NA)
mtext("B. Aminoglycosides", at = xmin - 0.54*(xmax-xmin), side = 3, line = 0.5, adj = 0, cex = 0.5, font = 1)



#### 5. (Fig.C) #####
xmin <- 2003; xmax <- 2019
ymin <- -2 ; ymax <- -0


plot(LRI_income$Cephalosporins_Median ~ LRI_income$Year, 
     data =  LRI_income,
     xlab = "", 
     ylab = "",
     pch = 16, col = rgb(0.8, 0.8, 0.8, 0),
     xlim = c(xmin,xmax),
     ylim = c(ymin,ymax),
     lwd = "1",
     axes=F,ann=F,
     main = "")

yticks <- log10(c(1e-2, 1e-1, 1e-0))   
ylabels <- c(1e-2, 1e-1, 1e-0)         
minor_labels <- as.vector(sapply(c(1e-2, 1e-1), function(v) v * 2:9))


axis(1, at = seq(xmin, xmax, by = 4), labels = FALSE, las = 2, pos = ymin)
text(x = seq(xmin, xmax, by = 4), y = ymin-(ymax-ymin)*0.08,  
     labels = format(round(seq(xmin, xmax, by = (xmax-xmin)/6*0+4), digits = 0), nsmall = 0), xpd = TRUE, cex = 0.65)

axis(2,  at = yticks, labels = FALSE, las = 2, pos = xmin)
axis(2, at = log10(minor_labels), labels = FALSE, las = 2, pos = xmin, tck = -0.015, lwd = 0.7)
text(y = yticks, x = xmin-0.04*(xmax-xmin),  
     labels = parse(text = paste0("10^", log10(ylabels))), srt = 0, adj = 1, xpd = TRUE, cex = 0.65)

lines(Cephalosporins_Median ~ Year, data = LRI_income_L,  col = "#CC79A7", lwd = 1.5)
lines(Cephalosporins_Median ~ Year, data = LRI_income_LM, col = "#009E73", lwd = 1.5)
lines(Cephalosporins_Median ~ Year, data = LRI_income_UM, col = "#0072B2", lwd = 1.5)

# Add shaded area for confidence interval
polygon(c(LRI_income_L$Year, rev(LRI_income_L$Year)),
        c(LRI_income_L$Cephalosporins_Q025, rev(LRI_income_L$Cephalosporins_Q975)),
        col = rgb(204/255, 121/255, 167/255, 0.25), border = NA)

polygon(c(LRI_income_LM$Year, rev(LRI_income_LM$Year)),
        c(LRI_income_LM$Cephalosporins_Q025, rev(LRI_income_LM$Cephalosporins_Q975)),
        col = rgb(0/255, 158/255, 115/255, 0.25), border = NA)

polygon(c(LRI_income_UM$Year, rev(LRI_income_UM$Year)),
        c(LRI_income_UM$Cephalosporins_Q025, rev(LRI_income_UM$Cephalosporins_Q975)),
        col = rgb(0/255, 114/255, 178/255, 0.25), border = NA)
mtext("C. Cephalosporins", at = xmin - 0.54*(xmax-xmin), side = 3, line = 0.5, adj = 0, cex = 0.5, font = 1)



#### 6. (Fig.D) #####
xmin <- 2003; xmax <- 2019
ymin <- -2 ; ymax <- -0


plot(LRI_income$Macrolides_Median ~ LRI_income$Year, 
     data =  LRI_income,
     xlab = "", 
     ylab = "",
     pch = 16, col = rgb(0.8, 0.8, 0.8, 0),
     xlim = c(xmin,xmax),
     ylim = c(ymin,ymax),
     lwd = "1",
     axes=F,ann=F,
     main = "")

yticks <- log10(c(1e-2, 1e-1, 1e-0))   
ylabels <- c(1e-2, 1e-1, 1e-0)         
minor_labels <- as.vector(sapply(c(1e-2, 1e-1), function(v) v * 2:9))


axis(1, at = seq(xmin, xmax, by = 4), labels = FALSE, las = 2, pos = ymin)
text(x = seq(xmin, xmax, by = 4), y = ymin-(ymax-ymin)*0.08,  
     labels = format(round(seq(xmin, xmax, by = (xmax-xmin)/6*0+4), digits = 0), nsmall = 0), xpd = TRUE, cex = 0.65)

axis(2,  at = yticks, labels = FALSE, las = 2, pos = xmin)
axis(2, at = log10(minor_labels), labels = FALSE, las = 2, pos = xmin, tck = -0.015, lwd = 0.7)
text(y = yticks, x = xmin-0.04*(xmax-xmin),  
     labels = parse(text = paste0("10^", log10(ylabels))), srt = 0, adj = 1, xpd = TRUE, cex = 0.65)

lines(Macrolides_Median ~ Year, data = LRI_income_L,  col = "#CC79A7", lwd = 1.5)
lines(Macrolides_Median ~ Year, data = LRI_income_LM, col = "#009E73", lwd = 1.5)
lines(Macrolides_Median ~ Year, data = LRI_income_UM, col = "#0072B2", lwd = 1.5)

# Add shaded area for confidence interval
polygon(c(LRI_income_L$Year, rev(LRI_income_L$Year)),
        c(LRI_income_L$Macrolides_Q025, rev(LRI_income_L$Macrolides_Q975)),
        col = rgb(204/255, 121/255, 167/255, 0.25), border = NA)

polygon(c(LRI_income_LM$Year, rev(LRI_income_LM$Year)),
        c(LRI_income_LM$Macrolides_Q025, rev(LRI_income_LM$Macrolides_Q975)),
        col = rgb(0/255, 158/255, 115/255, 0.25), border = NA)

polygon(c(LRI_income_UM$Year, rev(LRI_income_UM$Year)),
        c(LRI_income_UM$Macrolides_Q025, rev(LRI_income_UM$Macrolides_Q975)),
        col = rgb(0/255, 114/255, 178/255, 0.25), border = NA)
mtext("D. Macrolides", at = xmin - 0.54*(xmax-xmin), side = 3, line = 0.5, adj = 0, cex = 0.5, font = 1)



#### 7. (Fig.E) #####
xmin <- 2003; xmax <- 2019
ymin <- -5 ; ymax <- -2

plot(LRI_income$Carbapenems_Median ~ LRI_income$Year, 
     data =  LRI_income,
     xlab = "", 
     ylab = "",
     pch = 16, col = rgb(0.8, 0.8, 0.8, 0),
     xlim = c(xmin,xmax),
     ylim = c(ymin,ymax),
     lwd = "1",
     axes=F,ann=F,
     main = "")

yticks <- log10(c(1e-5, 1e-4, 1e-3, 1e-2))   
yticks_label <- log10(c(1e-5, 1e-4, 1e-3, 1e-2))   
ylabels <- c(1e-5, 1e-4, 1e-3, 1e-2)         
minor_labels <- as.vector(sapply(c(1e-5, 1e-4, 1e-3), function(v) v * 2:9))


axis(1, at = seq(xmin, xmax, by = 4), labels = FALSE, las = 2, pos = ymin)
text(x = seq(xmin, xmax, by = 4), y = ymin-(ymax-ymin)*0.08,  
     labels = format(round(seq(xmin, xmax, by = (xmax-xmin)/6*0+4), digits = 0), nsmall = 0), xpd = TRUE, cex = 0.65)

axis(2,  at = yticks, labels = FALSE, las = 2, pos = xmin)
axis(2, at = log10(minor_labels), labels = FALSE, las = 2, pos = xmin, tck = -0.015, lwd = 0.7)
text(y = yticks_label, x = xmin-0.04*(xmax-xmin),  
     labels = parse(text = paste0("10^", log10(ylabels))), srt = 0, adj = 1, xpd = TRUE, cex = 0.65)

lines(Carbapenems_Median ~ Year, data = LRI_income_L,  col = "#CC79A7", lwd = 1.5)
lines(Carbapenems_Median ~ Year, data = LRI_income_LM, col = "#009E73", lwd = 1.5)
lines(Carbapenems_Median ~ Year, data = LRI_income_UM, col = "#0072B2", lwd = 1.5)

# Add shaded area for confidence interval
polygon(c(LRI_income_L$Year, rev(LRI_income_L$Year)),
        c(LRI_income_L$Carbapenems_Q025, rev(LRI_income_L$Carbapenems_Q975)),
        col = rgb(204/255, 121/255, 167/255, 0.25), border = NA)

polygon(c(LRI_income_LM$Year, rev(LRI_income_LM$Year)),
        c(LRI_income_LM$Carbapenems_Q025, rev(LRI_income_LM$Carbapenems_Q975)),
        col = rgb(0/255, 158/255, 115/255, 0.25), border = NA)

polygon(c(LRI_income_UM$Year, rev(LRI_income_UM$Year)),
        c(LRI_income_UM$Carbapenems_Q025, rev(LRI_income_UM$Carbapenems_Q975)),
        col = rgb(0/255, 114/255, 178/255, 0.25), border = NA)
mtext("Year", side = 1, line = 1, cex = 0.5)
mtext("E. Carbapenems", at = xmin - 0.54*(xmax-xmin), side = 3, line = 0.5, adj = 0, cex = 0.5, font = 1)
mtext(expression("Proportion of antibiotic \n consumption by class*"), side = 2, line = 2, cex = 0.5)  



#### 8. (Fig.F) #####
xmin <- 2003; xmax <- 2019
ymin <- -6 ; ymax <- -2


plot(LRI_income$Glycopeptides_Median ~ LRI_income$Year, 
     data =  LRI_income,
     xlab = "", 
     ylab = "",
     pch = 16, col = rgb(0.8, 0.8, 0.8, 0),
     xlim = c(xmin,xmax),
     ylim = c(ymin,ymax),
     lwd = "1",
     axes=F,ann=F,
     main = "")



yticks <- log10(c(1e-6, 1e-5, 1e-4, 1e-3, 1e-2))   
yticks_label <- log10(c(1e-6, 1e-5, 1e-4, 1e-3, 1e-2))         
ylabels <- c(1e-6, 1e-5, 1e-4, 1e-3, 1e-2)         
minor_labels <- as.vector(sapply(c(1e-6, 1e-5, 1e-4, 1e-3), function(v) v * 2:9))


axis(1, at = seq(xmin, xmax, by = 4), labels = FALSE, las = 2, pos = ymin)
text(x = seq(xmin, xmax, by = 4), y = ymin-(ymax-ymin)*0.08,  
     labels = format(round(seq(xmin, xmax, by = (xmax-xmin)/6*0+4), digits = 0), nsmall = 0), xpd = TRUE, cex = 0.65)

axis(2,  at = yticks, labels = FALSE, las = 2, pos = xmin)
axis(2, at = log10(minor_labels), labels = FALSE, las = 2, pos = xmin, tck = -0.015, lwd = 0.7)
text(y = yticks_label, x = xmin-0.04*(xmax-xmin),  
     labels = parse(text = paste0("10^", log10(ylabels))), srt = 0, adj = 1, xpd = TRUE, cex = 0.65)

lines(Glycopeptides_Median ~ Year, data = LRI_income_L,  col = "#CC79A7", lwd = 1.5)
lines(Glycopeptides_Median ~ Year, data = LRI_income_LM, col = "#009E73", lwd = 1.5)
lines(Glycopeptides_Median ~ Year, data = LRI_income_UM, col = "#0072B2", lwd = 1.5)

# Add shaded area for confidence interval
polygon(c(LRI_income_L$Year, rev(LRI_income_L$Year)),
        c(LRI_income_L$Glycopeptides_Q025, rev(LRI_income_L$Glycopeptides_Q975)),
        col = rgb(204/255, 121/255, 167/255, 0.25), border = NA)

polygon(c(LRI_income_LM$Year, rev(LRI_income_LM$Year)),
        c(LRI_income_LM$Glycopeptides_Q025, rev(LRI_income_LM$Glycopeptides_Q975)),
        col = rgb(0/255, 158/255, 115/255, 0.25), border = NA)

polygon(c(LRI_income_UM$Year, rev(LRI_income_UM$Year)),
        c(LRI_income_UM$Glycopeptides_Q025, rev(LRI_income_UM$Glycopeptides_Q975)),
        col = rgb(0/255, 114/255, 178/255, 0.25), border = NA)
mtext("Year", side = 1, line = 1, cex = 0.5)
mtext("F. Glycopeptides", at = xmin - 0.54*(xmax-xmin), side = 3, line = 0.5, adj = 0, cex = 0.5, font = 1)



#### 9. (Fig.G) #####
xmin <- 2003; xmax <- 2019
ymin <- -2 ; ymax <- -1


plot(LRI_income$NSP_Median ~ LRI_income$Year, 
     data =  LRI_income,
     xlab = "", 
     ylab = "",
     pch = 16, col = rgb(0.8, 0.8, 0.8, 0),
     xlim = c(xmin,xmax),
     ylim = c(ymin,ymax),
     lwd = "1",
     axes=F,ann=F,
     main = "")



yticks <- log10(c(1e-2, 1e-1))   
ylabels <- c(1e-2, 1e-1)         
minor_labels <- as.vector(sapply(c(1e-2), function(v) v * 2:9))


axis(1, at = seq(xmin, xmax, by = 4), labels = FALSE, las = 2, pos = ymin)
text(x = seq(xmin, xmax, by = 4), y = ymin-(ymax-ymin)*0.08,  
     labels = format(round(seq(xmin, xmax, by = (xmax-xmin)/6*0+4), digits = 0), nsmall = 0), xpd = TRUE, cex = 0.65)

axis(2,  at = yticks, labels = FALSE, las = 2, pos = xmin)
axis(2, at = log10(minor_labels), labels = FALSE, las = 2, pos = xmin, tck = -0.015, lwd = 0.7)
text(y = yticks, x = xmin-0.04*(xmax-xmin),  
     labels = parse(text = paste0("10^", log10(ylabels))), srt = 0, adj = 1, xpd = TRUE, cex = 0.65)

lines(NSP_Median ~ Year, data = LRI_income_L,  col = "#CC79A7", lwd = 1.5)
lines(NSP_Median ~ Year, data = LRI_income_LM, col = "#009E73", lwd = 1.5)
lines(NSP_Median ~ Year, data = LRI_income_UM, col = "#0072B2", lwd = 1.5)

# Add shaded area for confidence interval
polygon(c(LRI_income_L$Year, rev(LRI_income_L$Year)),
        c(LRI_income_L$NSP_Q025, rev(LRI_income_L$NSP_Q975)),
        col = rgb(204/255, 121/255, 167/255, 0.25), border = NA)

polygon(c(LRI_income_LM$Year, rev(LRI_income_LM$Year)),
        c(LRI_income_LM$NSP_Q025, rev(LRI_income_LM$NSP_Q975)),
        col = rgb(0/255, 158/255, 115/255, 0.25), border = NA)

polygon(c(LRI_income_UM$Year, rev(LRI_income_UM$Year)),
        c(LRI_income_UM$NSP_Q025, rev(LRI_income_UM$NSP_Q975)),
        col = rgb(0/255, 114/255, 178/255, 0.25), border = NA)
mtext("G. NSP", at = xmin - 0.54*(xmax-xmin), side = 3, line = 0.5, adj = 0, cex = 0.5, font = 1)
mtext("Year", side = 1, line = 1, cex = 0.5)



#### 10. (Fig.H) #####
xmin <- 2003; xmax <- 2019
ymin <- -2 ; ymax <- -0


plot(LRI_income$Trimethoprim_Median ~ LRI_income$Year, 
     data =  LRI_income,
     xlab = "", 
     ylab = "",
     pch = 16, col = rgb(0.8, 0.8, 0.8, 0),
     xlim = c(xmin,xmax),
     ylim = c(ymin,ymax),
     lwd = "1",
     axes=F,ann=F,
     main = "")



yticks <- log10(c(1e-2, 1e-1, 1e-0))   
ylabels <- c(1e-2, 1e-1, 1e-0)         
minor_labels <- as.vector(sapply(c(1e-2, 1e-1), function(v) v * 2:9))


axis(1, at = seq(xmin, xmax, by = 4), labels = FALSE, las = 2, pos = ymin)
text(x = seq(xmin, xmax, by = 4), y = ymin-(ymax-ymin)*0.08,  
     labels = format(round(seq(xmin, xmax, by = (xmax-xmin)/6*0+4), digits = 0), nsmall = 0), xpd = TRUE, cex = 0.65)

axis(2,  at = yticks, labels = FALSE, las = 2, pos = xmin)
axis(2, at = log10(minor_labels), labels = FALSE, las = 2, pos = xmin, tck = -0.015, lwd = 0.7)
text(y = yticks, x = xmin-0.04*(xmax-xmin),  
     labels = parse(text = paste0("10^", log10(ylabels))), srt = 0, adj = 1, xpd = TRUE, cex = 0.65)

lines(Trimethoprim_Median ~ Year, data = LRI_income_L,  col = "#CC79A7", lwd = 1.5)
lines(Trimethoprim_Median ~ Year, data = LRI_income_LM, col = "#009E73", lwd = 1.5)
lines(Trimethoprim_Median ~ Year, data = LRI_income_UM, col = "#0072B2", lwd = 1.5)

# Add shaded area for confidence interval
polygon(c(LRI_income_L$Year, rev(LRI_income_L$Year)),
        c(LRI_income_L$Trimethoprim_Q025, rev(LRI_income_L$Trimethoprim_Q975)),
        col = rgb(204/255, 121/255, 167/255, 0.25), border = NA)

polygon(c(LRI_income_LM$Year, rev(LRI_income_LM$Year)),
        c(LRI_income_LM$Trimethoprim_Q025, rev(LRI_income_LM$Trimethoprim_Q975)),
        col = rgb(0/255, 158/255, 115/255, 0.25), border = NA)

polygon(c(LRI_income_UM$Year, rev(LRI_income_UM$Year)),
        c(LRI_income_UM$Trimethoprim_Q025, rev(LRI_income_UM$Trimethoprim_Q975)),
        col = rgb(0/255, 114/255, 178/255, 0.25), border = NA)
mtext("Year", side = 1, line = 1, cex = 0.5)
mtext("H. Trimethoprim", at = xmin - 0.54*(xmax-xmin), side = 3, line = 0.5, adj = 0, cex = 0.5, font = 1)



#### 10. (legend) ####
par(mar = c(0, 0, 0, 0))
plot.new()

legend_labels <- c("Upper Middle Income Countries", "Lower Middle Income Countries", "Low Income Countries")
legend_colors <- c("#0072B2", "#009E73", "#CC79A7")
n <- length(legend_labels)

# Coordinates for the custom legend
x0 <- 0.0  
x1 <- 0.07 
y_start <- 0.6
line_spacing <- 0.07

# Draw each line, label, and shaded polygon
for (i in seq_len(n)) {
  y <- y_start - (i - 1) * line_spacing
  

  polygon(
    x = c(x0-0.01, x1+0.01, x1+0.01, x0-0.01),
    y = c(y - 0.02, y - 0.02, y + 0.02, y + 0.02),
    col = adjustcolor(legend_colors[i], alpha.f = 0.25),
    border = NA
  )
  
  
  segments(x0, y, x1, y, col = legend_colors[i], lwd = 1.5)
  

  text(x = x1 + 0.04, y = y, labels = legend_labels[i], adj = 0, cex = 0.55)
}


#### 11. (Fig.I) #####
par(mar = c(2.5, 2.5, 2, 1))  
par(mgp = c(3, 0.3, 0))    
par(tck = -0.025)
par(lwd = 1)


xmin <- 2003; xmax <- 2019
ymin <- -1 ; ymax <- 0


plot(LRI_gbd_region$BSP_Median ~ LRI_gbd_region$Year, 
     data =  LRI_gbd_region,
     xlab = "", 
     ylab = "",
     pch = 16, col = rgb(0.8, 0.8, 0.8, 0),
     xlim = c(xmin,xmax),
     ylim = c(ymin,ymax),
     lwd = "1",
     axes=F,ann=F,
     main = "")



yticks <- log10(c(0.1,1))
ylabels <- c(0.1,1)    
minor_labels <- as.vector(sapply(0.1, function(v) v * 1:9))

axis(1, at = seq(xmin, xmax, by = 4), labels = FALSE, las = 2, pos = ymin)
text(x = seq(xmin, xmax, by = 4), y = ymin-(ymax-ymin)*0.08,  
     labels = format(round(seq(xmin, xmax, by = (xmax-xmin)/6*0+4), digits = 0), nsmall = 0), xpd = TRUE, cex = 0.65)

axis(2,  at = yticks, labels = FALSE, las = 2, pos = xmin)
axis(2, at = log10(minor_labels), labels = FALSE, las = 2, pos = xmin, tck = -0.015, lwd = 0.7)
text(y = yticks, x = xmin-0.04*(xmax-xmin),  
     labels = parse(text = paste0("10^", log10(ylabels))), srt = 0, adj = 1, xpd = TRUE, cex = 0.65)

text(y = log10(c(0.2, 0.3)),
     x = xmin - 0.04 * (xmax - xmin),
     labels = c(expression(2 %*% 10^-1), 
                expression(3 %*% 10^-1)),
     srt = 0, adj = 1, xpd = TRUE, cex = 0.5)  

lines(BSP_Median ~ Year, data = LRI_gbd_region_2, col = "#CC79A7", lwd = 1.5)
polygon(c(LRI_gbd_region_2$Year, rev(LRI_gbd_region_2$Year)),
        c(LRI_gbd_region_2$BSP_Q025, rev(LRI_gbd_region_2$BSP_Q975)),
        col = rgb(204/255, 121/255, 167/255, 0.25), border = NA)

lines(BSP_Median ~ Year, data = LRI_gbd_region_7, col = "#D55E00", lwd = 1.5)
polygon(c(LRI_gbd_region_7$Year, rev(LRI_gbd_region_7$Year)),
        c(LRI_gbd_region_7$BSP_Q025, rev(LRI_gbd_region_7$BSP_Q975)),
        col = rgb(213/255, 94/255, 0/255, 0.25), border = NA)

lines(BSP_Median ~ Year, data = LRI_gbd_region_6, col = "#E69F00", lwd = 1.5)
polygon(c(LRI_gbd_region_6$Year, rev(LRI_gbd_region_6$Year)),
        c(LRI_gbd_region_6$BSP_Q025, rev(LRI_gbd_region_6$BSP_Q975)),
        col = rgb(230/255, 159/255, 0/255, 0.25), border = NA)

lines(BSP_Median ~ Year, data = LRI_gbd_region_5, col = "#A7B97F", lwd = 1.5)
polygon(c(LRI_gbd_region_5$Year, rev(LRI_gbd_region_5$Year)),
        c(LRI_gbd_region_5$BSP_Q025, rev(LRI_gbd_region_5$BSP_Q975)),
        col = rgb(167/255, 185/255, 127/255, 0.25), border = NA)

lines(BSP_Median ~ Year, data = LRI_gbd_region_1, col = "#009E73", lwd = 1.5)
polygon(c(LRI_gbd_region_1$Year, rev(LRI_gbd_region_1$Year)),
        c(LRI_gbd_region_1$BSP_Q025, rev(LRI_gbd_region_1$BSP_Q975)),
        col = rgb(0/255, 158/255, 115/255, 0.25), border = NA)

lines(BSP_Median ~ Year, data = LRI_gbd_region_3, col = "#0072B2", lwd = 1.5)
polygon(c(LRI_gbd_region_3$Year, rev(LRI_gbd_region_3$Year)),
        c(LRI_gbd_region_3$BSP_Q025, rev(LRI_gbd_region_3$BSP_Q975)),
        col = rgb(0/255, 114/255, 178/255, 0.25), border = NA)

mtext(expression("Proportion of antibiotic \n consumption by class*"), side = 2, line = 2, cex = 0.5)  
mtext("I. BSP", at = xmin - 0.54*(xmax-xmin), side = 3, line = 0.5, adj = 0, cex = 0.5, font = 1)




#### 12. (Fig.J) #####


xmin <- 2003; xmax <- 2019
ymin <- -3 ; ymax <- -1


plot(LRI_gbd_region$Aminoglycosides_Median ~ LRI_gbd_region$Year, 
     data =  LRI_gbd_region,
     xlab = "", 
     ylab = "",
     pch = 16, col = rgb(0.8, 0.8, 0.8, 0),
     xlim = c(xmin,xmax),
     ylim = c(ymin,ymax),
     lwd = "1",
     axes=F,ann=F,
     main = "")



# Example: log10 scale
yticks <- log10(c(1e-3, 1e-2, 1e-1))   
ylabels <- c(1e-3, 1e-2, 1e-1)         
minor_labels <- as.vector(sapply(c(1e-3, 1e-2), function(v) v * 2:9))

axis(1, at = seq(xmin, xmax, by = 4), labels = FALSE, las = 2, pos = ymin)
text(x = seq(xmin, xmax, by = 4), y = ymin-(ymax-ymin)*0.08,  
     labels = format(round(seq(xmin, xmax, by = (xmax-xmin)/6*0+4), digits = 0), nsmall = 0), xpd = TRUE, cex = 0.65)

axis(2,  at = yticks, labels = FALSE, las = 2, pos = xmin)
axis(2, at = log10(minor_labels), labels = FALSE, las = 2, pos = xmin, tck = -0.015, lwd = 0.7)
text(y = yticks, x = xmin-0.04*(xmax-xmin),  
     labels = parse(text = paste0("10^", log10(ylabels))), srt = 0, adj = 1, xpd = TRUE, cex = 0.65)

lines(Aminoglycosides_Median ~ Year, data = LRI_gbd_region_2, col = "#CC79A7", lwd = 1.5)
polygon(c(LRI_gbd_region_2$Year, rev(LRI_gbd_region_2$Year)),
        c(LRI_gbd_region_2$Aminoglycosides_Q025, rev(LRI_gbd_region_2$Aminoglycosides_Q975)),
        col = rgb(204/255, 121/255, 167/255, 0.25), border = NA)

lines(Aminoglycosides_Median ~ Year, data = LRI_gbd_region_7, col = "#D55E00", lwd = 1.5)
polygon(c(LRI_gbd_region_7$Year, rev(LRI_gbd_region_7$Year)),
        c(LRI_gbd_region_7$Aminoglycosides_Q025, rev(LRI_gbd_region_7$Aminoglycosides_Q975)),
        col = rgb(213/255, 94/255, 0/255, 0.25), border = NA)

lines(Aminoglycosides_Median ~ Year, data = LRI_gbd_region_6, col = "#E69F00", lwd = 1.5)
polygon(c(LRI_gbd_region_6$Year, rev(LRI_gbd_region_6$Year)),
        c(LRI_gbd_region_6$Aminoglycosides_Q025, rev(LRI_gbd_region_6$Aminoglycosides_Q975)),
        col = rgb(230/255, 159/255, 0/255, 0.25), border = NA)

lines(Aminoglycosides_Median ~ Year, data = LRI_gbd_region_5, col = "#A7B97F", lwd = 1.5)
polygon(c(LRI_gbd_region_5$Year, rev(LRI_gbd_region_5$Year)),
        c(LRI_gbd_region_5$Aminoglycosides_Q025, rev(LRI_gbd_region_5$Aminoglycosides_Q975)),
        col = rgb(167/255, 185/255, 127/255, 0.25), border = NA)

lines(Aminoglycosides_Median ~ Year, data = LRI_gbd_region_1, col = "#009E73", lwd = 1.5)
polygon(c(LRI_gbd_region_1$Year, rev(LRI_gbd_region_1$Year)),
        c(LRI_gbd_region_1$Aminoglycosides_Q025, rev(LRI_gbd_region_1$Aminoglycosides_Q975)),
        col = rgb(0/255, 158/255, 115/255, 0.25), border = NA)

lines(Aminoglycosides_Median ~ Year, data = LRI_gbd_region_3, col = "#0072B2", lwd = 1.5)
polygon(c(LRI_gbd_region_3$Year, rev(LRI_gbd_region_3$Year)),
        c(LRI_gbd_region_3$Aminoglycosides_Q025, rev(LRI_gbd_region_3$Aminoglycosides_Q975)),
        col = rgb(0/255, 114/255, 178/255, 0.25), border = NA)
mtext("J. Aminoglycosides", at = xmin - 0.54*(xmax-xmin), side = 3, line = 0.5, adj = 0, cex = 0.5, font = 1)


#### 13. (Fig.K) #####

xmin <- 2003; xmax <- 2019
ymin <- -2 ; ymax <- -0

plot(LRI_gbd_region$Cephalosporins_Median ~ LRI_gbd_region$Year, 
     data =  LRI_gbd_region,
     xlab = "", 
     ylab = "",
     pch = 16, col = rgb(0.8, 0.8, 0.8, 0),
     xlim = c(xmin,xmax),
     ylim = c(ymin,ymax),
     lwd = "1",
     axes=F,ann=F,
     main = "")

yticks <- log10(c(1e-2, 1e-1, 1e-0))   
ylabels <- c(1e-2, 1e-1, 1e-0)         
minor_labels <- as.vector(sapply(c(1e-2, 1e-1), function(v) v * 2:9))


axis(1, at = seq(xmin, xmax, by = 4), labels = FALSE, las = 2, pos = ymin)
text(x = seq(xmin, xmax, by = 4), y = ymin-(ymax-ymin)*0.08,  
     labels = format(round(seq(xmin, xmax, by = (xmax-xmin)/6*0+4), digits = 0), nsmall = 0), xpd = TRUE, cex = 0.65)

axis(2,  at = yticks, labels = FALSE, las = 2, pos = xmin)
axis(2, at = log10(minor_labels), labels = FALSE, las = 2, pos = xmin, tck = -0.015, lwd = 0.7)
text(y = yticks, x = xmin-0.04*(xmax-xmin),  
     labels = parse(text = paste0("10^", log10(ylabels))), srt = 0, adj = 1, xpd = TRUE, cex = 0.65)

lines(Cephalosporins_Median ~ Year, data = LRI_gbd_region_2, col = "#CC79A7", lwd = 1.5)
polygon(c(LRI_gbd_region_2$Year, rev(LRI_gbd_region_2$Year)),
        c(LRI_gbd_region_2$Cephalosporins_Q025, rev(LRI_gbd_region_2$Cephalosporins_Q975)),
        col = rgb(204/255, 121/255, 167/255, 0.25), border = NA)

lines(Cephalosporins_Median ~ Year, data = LRI_gbd_region_7, col = "#D55E00", lwd = 1.5)
polygon(c(LRI_gbd_region_7$Year, rev(LRI_gbd_region_7$Year)),
        c(LRI_gbd_region_7$Cephalosporins_Q025, rev(LRI_gbd_region_7$Cephalosporins_Q975)),
        col = rgb(213/255, 94/255, 0/255, 0.25), border = NA)

lines(Cephalosporins_Median ~ Year, data = LRI_gbd_region_6, col = "#E69F00", lwd = 1.5)
polygon(c(LRI_gbd_region_6$Year, rev(LRI_gbd_region_6$Year)),
        c(LRI_gbd_region_6$Cephalosporins_Q025, rev(LRI_gbd_region_6$Cephalosporins_Q975)),
        col = rgb(230/255, 159/255, 0/255, 0.25), border = NA)

lines(Cephalosporins_Median ~ Year, data = LRI_gbd_region_5, col = "#A7B97F", lwd = 1.5)
polygon(c(LRI_gbd_region_5$Year, rev(LRI_gbd_region_5$Year)),
        c(LRI_gbd_region_5$Cephalosporins_Q025, rev(LRI_gbd_region_5$Cephalosporins_Q975)),
        col = rgb(167/255, 185/255, 127/255, 0.25), border = NA)

lines(Cephalosporins_Median ~ Year, data = LRI_gbd_region_1, col = "#009E73", lwd = 1.5)
polygon(c(LRI_gbd_region_1$Year, rev(LRI_gbd_region_1$Year)),
        c(LRI_gbd_region_1$Cephalosporins_Q025, rev(LRI_gbd_region_1$Cephalosporins_Q975)),
        col = rgb(0/255, 158/255, 115/255, 0.25), border = NA)

lines(Cephalosporins_Median ~ Year, data = LRI_gbd_region_3, col = "#0072B2", lwd = 1.5)
polygon(c(LRI_gbd_region_3$Year, rev(LRI_gbd_region_3$Year)),
        c(LRI_gbd_region_3$Cephalosporins_Q025, rev(LRI_gbd_region_3$Cephalosporins_Q975)),
        col = rgb(0/255, 114/255, 178/255, 0.25), border = NA)
mtext("K. Cephalosporins", at = xmin - 0.54*(xmax-xmin), side = 3, line = 0.5, adj = 0, cex = 0.5, font = 1)



#### 14. (Fig.L) #####

xmin <- 2003; xmax <- 2019
ymin <- -2 ; ymax <- -0

plot(LRI_gbd_region$Macrolides_Median ~ LRI_gbd_region$Year, 
     data =  LRI_gbd_region,
     xlab = "", 
     ylab = "",
     pch = 16, col = rgb(0.8, 0.8, 0.8, 0),
     xlim = c(xmin,xmax),
     ylim = c(ymin,ymax),
     lwd = "1",
     axes=F,ann=F,
     main = "")

yticks <- log10(c(1e-2, 1e-1, 1e-0))   
ylabels <- c(1e-2, 1e-1, 1e-0)         
minor_labels <- as.vector(sapply(c(1e-2, 1e-1), function(v) v * 2:9))


axis(1, at = seq(xmin, xmax, by = 4), labels = FALSE, las = 2, pos = ymin)
text(x = seq(xmin, xmax, by = 4), y = ymin-(ymax-ymin)*0.08,  
     labels = format(round(seq(xmin, xmax, by = (xmax-xmin)/6*0+4), digits = 0), nsmall = 0), xpd = TRUE, cex = 0.65)

axis(2,  at = yticks, labels = FALSE, las = 2, pos = xmin)
axis(2, at = log10(minor_labels), labels = FALSE, las = 2, pos = xmin, tck = -0.015, lwd = 0.7)
text(y = yticks, x = xmin-0.04*(xmax-xmin),  
     labels = parse(text = paste0("10^", log10(ylabels))), srt = 0, adj = 1, xpd = TRUE, cex = 0.65)

lines(Macrolides_Median ~ Year, data = LRI_gbd_region_2, col = "#CC79A7", lwd = 1.5)
polygon(c(LRI_gbd_region_2$Year, rev(LRI_gbd_region_2$Year)),
        c(LRI_gbd_region_2$Macrolides_Q025, rev(LRI_gbd_region_2$Macrolides_Q975)),
        col = rgb(204/255, 121/255, 167/255, 0.25), border = NA)

lines(Macrolides_Median ~ Year, data = LRI_gbd_region_7, col = "#D55E00", lwd = 1.5)
polygon(c(LRI_gbd_region_7$Year, rev(LRI_gbd_region_7$Year)),
        c(LRI_gbd_region_7$Macrolides_Q025, rev(LRI_gbd_region_7$Macrolides_Q975)),
        col = rgb(213/255, 94/255, 0/255, 0.25), border = NA)

lines(Macrolides_Median ~ Year, data = LRI_gbd_region_6, col = "#E69F00", lwd = 1.5)
polygon(c(LRI_gbd_region_6$Year, rev(LRI_gbd_region_6$Year)),
        c(LRI_gbd_region_6$Macrolides_Q025, rev(LRI_gbd_region_6$Macrolides_Q975)),
        col = rgb(230/255, 159/255, 0/255, 0.25), border = NA)

lines(Macrolides_Median ~ Year, data = LRI_gbd_region_5, col = "#A7B97F", lwd = 1.5)
polygon(c(LRI_gbd_region_5$Year, rev(LRI_gbd_region_5$Year)),
        c(LRI_gbd_region_5$Macrolides_Q025, rev(LRI_gbd_region_5$Macrolides_Q975)),
        col = rgb(167/255, 185/255, 127/255, 0.25), border = NA)

lines(Macrolides_Median ~ Year, data = LRI_gbd_region_1, col = "#009E73", lwd = 1.5)
polygon(c(LRI_gbd_region_1$Year, rev(LRI_gbd_region_1$Year)),
        c(LRI_gbd_region_1$Macrolides_Q025, rev(LRI_gbd_region_1$Macrolides_Q975)),
        col = rgb(0/255, 158/255, 115/255, 0.25), border = NA)

lines(Macrolides_Median ~ Year, data = LRI_gbd_region_3, col = "#0072B2", lwd = 1.5)
polygon(c(LRI_gbd_region_3$Year, rev(LRI_gbd_region_3$Year)),
        c(LRI_gbd_region_3$Macrolides_Q025, rev(LRI_gbd_region_3$Macrolides_Q975)),
        col = rgb(0/255, 114/255, 178/255, 0.25), border = NA)
mtext("L. Macrolides", at = xmin - 0.54*(xmax-xmin), side = 3, line = 0.5, adj = 0, cex = 0.5, font = 1)






#### 15. (Fig.M) #####
xmin <- 2003; xmax <- 2019
ymin <- -5 ; ymax <- -2

plot(LRI_gbd_region$Carbapenems_Median ~ LRI_gbd_region$Year, 
     data =  LRI_gbd_region,
     xlab = "", 
     ylab = "",
     pch = 16, col = rgb(0.8, 0.8, 0.8, 0),
     xlim = c(xmin,xmax),
     ylim = c(ymin,ymax),
     lwd = "1",
     axes=F,ann=F,
     main = "")

yticks <- log10(c(1e-5, 1e-4, 1e-3, 1e-2))   
yticks_label <- log10(c(1e-5, 1e-4, 1e-3, 1e-2))   
ylabels <- c(1e-5, 1e-4, 1e-3, 1e-2)         
minor_labels <- as.vector(sapply(c(1e-5, 1e-4, 1e-3), function(v) v * 2:9))


axis(1, at = seq(xmin, xmax, by = 4), labels = FALSE, las = 2, pos = ymin)
text(x = seq(xmin, xmax, by = 4), y = ymin-(ymax-ymin)*0.08,  
     labels = format(round(seq(xmin, xmax, by = (xmax-xmin)/6*0+4), digits = 0), nsmall = 0), xpd = TRUE, cex = 0.65)

axis(2,  at = yticks, labels = FALSE, las = 2, pos = xmin)
axis(2, at = log10(minor_labels), labels = FALSE, las = 2, pos = xmin, tck = -0.015, lwd = 0.7)
text(y = yticks_label, x = xmin-0.04*(xmax-xmin),  
     labels = parse(text = paste0("10^", log10(ylabels))), srt = 0, adj = 1, xpd = TRUE, cex = 0.65)


lines(Carbapenems_Median ~ Year, data = LRI_gbd_region_2, col = "#CC79A7", lwd = 1.5)
polygon(c(LRI_gbd_region_2$Year, rev(LRI_gbd_region_2$Year)),
        c(LRI_gbd_region_2$Carbapenems_Q025, rev(LRI_gbd_region_2$Carbapenems_Q975)),
        col = rgb(204/255, 121/255, 167/255, 0.25), border = NA)

lines(Carbapenems_Median ~ Year, data = LRI_gbd_region_7, col = "#D55E00", lwd = 1.5)
polygon(c(LRI_gbd_region_7$Year, rev(LRI_gbd_region_7$Year)),
        c(LRI_gbd_region_7$Carbapenems_Q025, rev(LRI_gbd_region_7$Carbapenems_Q975)),
        col = rgb(213/255, 94/255, 0/255, 0.25), border = NA)

lines(Carbapenems_Median ~ Year, data = LRI_gbd_region_6, col = "#E69F00", lwd = 1.5)
polygon(c(LRI_gbd_region_6$Year, rev(LRI_gbd_region_6$Year)),
        c(LRI_gbd_region_6$Carbapenems_Q025, rev(LRI_gbd_region_6$Carbapenems_Q975)),
        col = rgb(230/255, 159/255, 0/255, 0.25), border = NA)

lines(Carbapenems_Median ~ Year, data = LRI_gbd_region_5, col = "#A7B97F", lwd = 1.5)
polygon(c(LRI_gbd_region_5$Year, rev(LRI_gbd_region_5$Year)),
        c(LRI_gbd_region_5$Carbapenems_Q025, rev(LRI_gbd_region_5$Carbapenems_Q975)),
        col = rgb(167/255, 185/255, 127/255, 0.25), border = NA)

lines(Carbapenems_Median ~ Year, data = LRI_gbd_region_1, col = "#009E73", lwd = 1.5)
polygon(c(LRI_gbd_region_1$Year, rev(LRI_gbd_region_1$Year)),
        c(LRI_gbd_region_1$Carbapenems_Q025, rev(LRI_gbd_region_1$Carbapenems_Q975)),
        col = rgb(0/255, 158/255, 115/255, 0.25), border = NA)

lines(Carbapenems_Median ~ Year, data = LRI_gbd_region_3, col = "#0072B2", lwd = 1.5)
polygon(c(LRI_gbd_region_3$Year, rev(LRI_gbd_region_3$Year)),
        c(LRI_gbd_region_3$Carbapenems_Q025, rev(LRI_gbd_region_3$Carbapenems_Q975)),
        col = rgb(0/255, 114/255, 178/255, 0.25), border = NA)
mtext("Year", side = 1, line = 1, cex = 0.5)
mtext("M. Carbapenems", at = xmin - 0.54*(xmax-xmin), side = 3, line = 0.5, adj = 0, cex = 0.5, font = 1)
mtext(expression("Proportion of antibiotic \n consumption by class*"), side = 2, line = 2, cex = 0.5)  




#### 16. (Fig.N) #####
xmin <- 2003; xmax <- 2019
ymin <- -6 ; ymax <- -2

plot(LRI_gbd_region$Glycopeptides_Median ~ LRI_gbd_region$Year, 
     data =  LRI_gbd_region,
     xlab = "", 
     ylab = "",
     pch = 16, col = rgb(0.8, 0.8, 0.8, 0),
     xlim = c(xmin,xmax),
     ylim = c(ymin,ymax),
     lwd = "1",
     axes=F,ann=F,
     main = "")

yticks <- log10(c(1e-6, 1e-5, 1e-4, 1e-3, 1e-2))   
yticks_label <- log10(c(1e-6, 1e-5, 1e-4, 1e-3, 1e-2))   
ylabels <- c(1e-6, 1e-5, 1e-4, 1e-3, 1e-2)         
minor_labels <- as.vector(sapply(c(1e-6, 1e-5, 1e-4, 1e-3), function(v) v * 2:9))


axis(1, at = seq(xmin, xmax, by = 4), labels = FALSE, las = 2, pos = ymin)
text(x = seq(xmin, xmax, by = 4), y = ymin-(ymax-ymin)*0.08,  
     labels = format(round(seq(xmin, xmax, by = (xmax-xmin)/6*0+4), digits = 0), nsmall = 0), xpd = TRUE, cex = 0.65)

axis(2,  at = yticks, labels = FALSE, las = 2, pos = xmin)
axis(2, at = log10(minor_labels), labels = FALSE, las = 2, pos = xmin, tck = -0.015, lwd = 0.7)
text(y = yticks_label, x = xmin-0.04*(xmax-xmin),  
     labels = parse(text = paste0("10^", log10(ylabels))), srt = 0, adj = 1, xpd = TRUE, cex = 0.65)

# Add shaded area for confidence interval
lines(Glycopeptides_Median ~ Year, data = LRI_gbd_region_2, col = "#CC79A7", lwd = 1.5)
polygon(c(LRI_gbd_region_2$Year, rev(LRI_gbd_region_2$Year)),
        c(LRI_gbd_region_2$Glycopeptides_Q025, rev(LRI_gbd_region_2$Glycopeptides_Q975)),
        col = rgb(204/255, 121/255, 167/255, 0.25), border = NA)

lines(Glycopeptides_Median ~ Year, data = LRI_gbd_region_7, col = "#D55E00", lwd = 1.5)
polygon(c(LRI_gbd_region_7$Year, rev(LRI_gbd_region_7$Year)),
        c(LRI_gbd_region_7$Glycopeptides_Q025, rev(LRI_gbd_region_7$Glycopeptides_Q975)),
        col = rgb(213/255, 94/255, 0/255, 0.25), border = NA)

lines(Glycopeptides_Median ~ Year, data = LRI_gbd_region_6, col = "#E69F00", lwd = 1.5)
polygon(c(LRI_gbd_region_6$Year, rev(LRI_gbd_region_6$Year)),
        c(LRI_gbd_region_6$Glycopeptides_Q025, rev(LRI_gbd_region_6$Glycopeptides_Q975)),
        col = rgb(230/255, 159/255, 0/255, 0.25), border = NA)

lines(Glycopeptides_Median ~ Year, data = LRI_gbd_region_5, col = "#A7B97F", lwd = 1.5)
polygon(c(LRI_gbd_region_5$Year, rev(LRI_gbd_region_5$Year)),
        c(LRI_gbd_region_5$Glycopeptides_Q025, rev(LRI_gbd_region_5$Glycopeptides_Q975)),
        col = rgb(167/255, 185/255, 127/255, 0.25), border = NA)

lines(Glycopeptides_Median ~ Year, data = LRI_gbd_region_1, col = "#009E73", lwd = 1.5)
polygon(c(LRI_gbd_region_1$Year, rev(LRI_gbd_region_1$Year)),
        c(LRI_gbd_region_1$Glycopeptides_Q025, rev(LRI_gbd_region_1$Glycopeptides_Q975)),
        col = rgb(0/255, 158/255, 115/255, 0.25), border = NA)

lines(Glycopeptides_Median ~ Year, data = LRI_gbd_region_3, col = "#0072B2", lwd = 1.5)
polygon(c(LRI_gbd_region_3$Year, rev(LRI_gbd_region_3$Year)),
        c(LRI_gbd_region_3$Glycopeptides_Q025, rev(LRI_gbd_region_3$Glycopeptides_Q975)),
        col = rgb(0/255, 114/255, 178/255, 0.25), border = NA)
mtext("Year", side = 1, line = 1, cex = 0.5)
mtext("N. Glycopeptides", at = xmin - 0.54*(xmax-xmin), side = 3, line = 0.5, adj = 0, cex = 0.5, font = 1)


#### 11. (Fig.O) #####

xmin <- 2003; xmax <- 2019
ymin <- -2 ; ymax <- -1

plot(LRI_gbd_region$NSP_Median ~ LRI_gbd_region$Year, 
     data =  LRI_gbd_region,
     xlab = "", 
     ylab = "",
     pch = 16, col = rgb(0.8, 0.8, 0.8, 0),
     xlim = c(xmin,xmax),
     ylim = c(ymin,ymax),
     lwd = "1",
     axes=F,ann=F,
     main = "")

yticks <- log10(c(1e-2, 1e-1))   
ylabels <- c(1e-2, 1e-1)         
minor_labels <- as.vector(sapply(c(1e-2, 1e-1), function(v) v * 2:9))


axis(1, at = seq(xmin, xmax, by = 4), labels = FALSE, las = 2, pos = ymin)
text(x = seq(xmin, xmax, by = 4), y = ymin-(ymax-ymin)*0.08,  
     labels = format(round(seq(xmin, xmax, by = (xmax-xmin)/6*0+4), digits = 0), nsmall = 0), xpd = TRUE, cex = 0.65)

axis(2,  at = yticks, labels = FALSE, las = 2, pos = xmin)
axis(2, at = log10(minor_labels), labels = FALSE, las = 2, pos = xmin, tck = -0.015, lwd = 0.7)
text(y = yticks, x = xmin-0.04*(xmax-xmin),  
     labels = parse(text = paste0("10^", log10(ylabels))), srt = 0, adj = 1, xpd = TRUE, cex = 0.65)

lines(NSP_Median ~ Year, data = LRI_gbd_region_2, col = "#CC79A7", lwd = 1.5)
polygon(c(LRI_gbd_region_2$Year, rev(LRI_gbd_region_2$Year)),
        c(LRI_gbd_region_2$NSP_Q025, rev(LRI_gbd_region_2$NSP_Q975)),
        col = rgb(204/255, 121/255, 167/255, 0.25), border = NA)

lines(NSP_Median ~ Year, data = LRI_gbd_region_7, col = "#D55E00", lwd = 1.5)
polygon(c(LRI_gbd_region_7$Year, rev(LRI_gbd_region_7$Year)),
        c(LRI_gbd_region_7$NSP_Q025, rev(LRI_gbd_region_7$NSP_Q975)),
        col = rgb(213/255, 94/255, 0/255, 0.25), border = NA)

lines(NSP_Median ~ Year, data = LRI_gbd_region_6, col = "#E69F00", lwd = 1.5)
polygon(c(LRI_gbd_region_6$Year, rev(LRI_gbd_region_6$Year)),
        c(LRI_gbd_region_6$NSP_Q025, rev(LRI_gbd_region_6$NSP_Q975)),
        col = rgb(230/255, 159/255, 0/255, 0.25), border = NA)

lines(NSP_Median ~ Year, data = LRI_gbd_region_5, col = "#A7B97F", lwd = 1.5)
polygon(c(LRI_gbd_region_5$Year, rev(LRI_gbd_region_5$Year)),
        c(LRI_gbd_region_5$NSP_Q025, rev(LRI_gbd_region_5$NSP_Q975)),
        col = rgb(167/255, 185/255, 127/255, 0.25), border = NA)

lines(NSP_Median ~ Year, data = LRI_gbd_region_1, col = "#009E73", lwd = 1.5)
polygon(c(LRI_gbd_region_1$Year, rev(LRI_gbd_region_1$Year)),
        c(LRI_gbd_region_1$NSP_Q025, rev(LRI_gbd_region_1$NSP_Q975)),
        col = rgb(0/255, 158/255, 115/255, 0.25), border = NA)

lines(NSP_Median ~ Year, data = LRI_gbd_region_3, col = "#0072B2", lwd = 1.5)
polygon(c(LRI_gbd_region_3$Year, rev(LRI_gbd_region_3$Year)),
        c(LRI_gbd_region_3$NSP_Q025, rev(LRI_gbd_region_3$NSP_Q975)),
        col = rgb(0/255, 114/255, 178/255, 0.25), border = NA)
mtext("O. NSP", at = xmin - 0.54*(xmax-xmin), side = 3, line = 0.5, adj = 0, cex = 0.5, font = 1)
mtext("Year", side = 1, line = 1, cex = 0.5)




#### 18. (Fig.P) #####

xmin <- 2003; xmax <- 2019
ymin <- -2 ; ymax <- -0

plot(LRI_gbd_region$Trimethoprim_Median ~ LRI_gbd_region$Year, 
     data =  LRI_gbd_region,
     xlab = "", 
     ylab = "",
     pch = 16, col = rgb(0.8, 0.8, 0.8, 0),
     xlim = c(xmin,xmax),
     ylim = c(ymin,ymax),
     lwd = "1",
     axes=F,ann=F,
     main = "")

yticks <- log10(c(1e-2, 1e-1, 1e-0))   
ylabels <- c(1e-2, 1e-1, 1e-0)         
minor_labels <- as.vector(sapply(c(1e-2, 1e-1), function(v) v * 2:9))


axis(1, at = seq(xmin, xmax, by = 4), labels = FALSE, las = 2, pos = ymin)
text(x = seq(xmin, xmax, by = 4), y = ymin-(ymax-ymin)*0.08,  
     labels = format(round(seq(xmin, xmax, by = (xmax-xmin)/6*0+4), digits = 0), nsmall = 0), xpd = TRUE, cex = 0.65)

axis(2,  at = yticks, labels = FALSE, las = 2, pos = xmin)
axis(2, at = log10(minor_labels), labels = FALSE, las = 2, pos = xmin, tck = -0.015, lwd = 0.7)
text(y = yticks, x = xmin-0.04*(xmax-xmin),  
     labels = parse(text = paste0("10^", log10(ylabels))), srt = 0, adj = 1, xpd = TRUE, cex = 0.65)

# Add shaded area for confidence interval
lines(Trimethoprim_Median ~ Year, data = LRI_gbd_region_2, col = "#CC79A7", lwd = 1.5)
polygon(c(LRI_gbd_region_2$Year, rev(LRI_gbd_region_2$Year)),
        c(LRI_gbd_region_2$Trimethoprim_Q025, rev(LRI_gbd_region_2$Trimethoprim_Q975)),
        col = rgb(204/255, 121/255, 167/255, 0.25), border = NA)

lines(Trimethoprim_Median ~ Year, data = LRI_gbd_region_7, col = "#D55E00", lwd = 1.5)
polygon(c(LRI_gbd_region_7$Year, rev(LRI_gbd_region_7$Year)),
        c(LRI_gbd_region_7$Trimethoprim_Q025, rev(LRI_gbd_region_7$Trimethoprim_Q975)),
        col = rgb(213/255, 94/255, 0/255, 0.25), border = NA)

lines(Trimethoprim_Median ~ Year, data = LRI_gbd_region_6, col = "#E69F00", lwd = 1.5)
polygon(c(LRI_gbd_region_6$Year, rev(LRI_gbd_region_6$Year)),
        c(LRI_gbd_region_6$Trimethoprim_Q025, rev(LRI_gbd_region_6$Trimethoprim_Q975)),
        col = rgb(230/255, 159/255, 0/255, 0.25), border = NA)

lines(Trimethoprim_Median ~ Year, data = LRI_gbd_region_5, col = "#A7B97F", lwd = 1.5)
polygon(c(LRI_gbd_region_5$Year, rev(LRI_gbd_region_5$Year)),
        c(LRI_gbd_region_5$Trimethoprim_Q025, rev(LRI_gbd_region_5$Trimethoprim_Q975)),
        col = rgb(167/255, 185/255, 127/255, 0.25), border = NA)

lines(Trimethoprim_Median ~ Year, data = LRI_gbd_region_1, col = "#009E73", lwd = 1.5)
polygon(c(LRI_gbd_region_1$Year, rev(LRI_gbd_region_1$Year)),
        c(LRI_gbd_region_1$Trimethoprim_Q025, rev(LRI_gbd_region_1$Trimethoprim_Q975)),
        col = rgb(0/255, 158/255, 115/255, 0.25), border = NA)

lines(Trimethoprim_Median ~ Year, data = LRI_gbd_region_3, col = "#0072B2", lwd = 1.5)
polygon(c(LRI_gbd_region_3$Year, rev(LRI_gbd_region_3$Year)),
        c(LRI_gbd_region_3$Trimethoprim_Q025, rev(LRI_gbd_region_3$Trimethoprim_Q975)),
        col = rgb(0/255, 114/255, 178/255, 0.25), border = NA)
mtext("Year", side = 1, line = 1, cex = 0.5)
mtext("P. Trimethoprim", at = xmin - 0.54*(xmax-xmin), side = 3, line = 0.5, adj = 0, cex = 0.5, font = 1)

#### 10. (legend) ####
par(mar = c(0, 0, 0, 0))
plot.new()

legend_labels <- c("Southeast Asia, East Asia, and Oceania",
                   "Central Europe, Eastern Europe, and Central Asia",
                   "Latin America and Caribbean",
                   "South Asia",
                   "North Africa and Middle East",
                   "Sub-Saharan Africa")
legend_colors <- c("#D55E00",
                   "#CC79A7",
                   "#E69F00",
                   "#A7B97F",
                   "#009E73",
                   "#0072B2")
n <- length(legend_labels)

# Coordinates for the custom legend
x0 <- 0.0  
x1 <- 0.07 
y_start <- 0.7
line_spacing <- 0.07

# Draw each line, label, and shaded polygon
for (i in seq_len(n)) {
  y <- y_start - (i - 1) * line_spacing
  
  polygon(
    x = c(x0-0.01, x1+0.01, x1+0.01, x0-0.01),
    y = c(y - 0.02, y - 0.02, y + 0.02, y + 0.02),
    col = adjustcolor(legend_colors[i], alpha.f = 0.25),
    border = NA
  )
  
  segments(x0, y, x1, y, col = legend_colors[i], lwd = 1.5)
  
  text(x = x1 + 0.04, y = y, labels = legend_labels[i], adj = 0, cex = 0.55)
}

dev.off()

#