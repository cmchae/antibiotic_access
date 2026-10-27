#### 0. prep ####
rm(list = ls())
gc()

requiredPackages = c("dplyr")
lib.path <- .libPaths()
for(p in requiredPackages){
  if(!require(p,character.only = TRUE)) install.packages(p, lib = lib.path)
  library(p,character.only = TRUE)
}

options(scipen = 100, digits = 6) # 지수표시 안되게
n_sim <- 10000
set.seed(123)

coef_df_LRI <- readRDS("outputs/coef_df_LRI.rds")
coef_df_LRI_DDD_by_4class <- readRDS("outputs/antibiotic_class/coef_df_LRI_DDD_by_4class.rds")

stopifnot(
  identical(
    names(coef_df_LRI),
    c("(Intercept)", "PC1", "PC2")
  ),
  ncol(coef_df_LRI_DDD_by_4class) == 12,
  !any(grepl("PC3", names(coef_df_LRI_DDD_by_4class)))
)

rm(
  coef_df_LRI,
  coef_df_LRI_DDD_by_4class
)


# plot data below were generated from the current PC1-PC2 models

gbd_region <- c("North Africa and Middle East","Central Europe, Eastern Europe, and Central Asia","Sub-Saharan Africa","South Asia","Latin America and Caribbean","Southeast Asia, East Asia, and Oceania")

for (i in 1:6) { # GBD_region
  path_LRI <- paste0("outputs/by_class_df/by_gbd_region_year/probability/LRI_", gbd_region[i], "4.csv")
  obj_LRI   <- paste0("LRI_year_proportion_GBD_region_", i)
  assign(obj_LRI,   read.csv(path_LRI,   stringsAsFactors = FALSE))
}

for (i in 1:6) { # GBD_region
  path_LRI_prop <- paste0("outputs/by_class_df/by_gbd_region_year/probability/LRI_prop_", gbd_region[i], "4.csv")
  obj_LRI_prop   <- paste0("LRI_year_proportion_GBD_region_prop_", i)
  assign(obj_LRI_prop,   read.csv(path_LRI_prop,   stringsAsFactors = FALSE))
}

Income <- c("L", "LM", "UM")

for (inc in Income) { # Income
  path_LRI_last   <- paste0("outputs/by_class_df/by_income_year/probability/LRI4_last5_", inc, ".csv")
  obj_LRI_last    <- paste0("LRI_year_proportion_income_last_", inc)
  assign(obj_LRI_last,   read.csv(path_LRI_last,   stringsAsFactors = FALSE))
}

for (inc in Income) { # Income
  path_LRI_last_prop   <- paste0("outputs/by_class_df/by_income_year/probability/LRI4_last5_prop_", inc, ".csv")
  obj_LRI_last_prop    <- paste0("LRI_year_proportion_income_last_prop_", inc)
  assign(obj_LRI_last_prop,   read.csv(path_LRI_last_prop,   stringsAsFactors = FALSE))
}

abx_names <- c(
  "BSP", "NSP", "Cephalosporins", "Aminoglycosides")

LRI_year_proportion_prop <- read.csv("outputs/by_class_df/by_year/probability/LRI_incremental_share_df_prop.csv")

LRI_incr_2019       <- LRI_year_proportion_prop               [LRI_year_proportion_prop               $Year == 2019,]
LRI_incr_2019_L     <- LRI_year_proportion_income_last_prop_L [LRI_year_proportion_income_last_prop_L $Year == 2019,]
LRI_incr_2019_LM    <- LRI_year_proportion_income_last_prop_LM[LRI_year_proportion_income_last_prop_LM$Year == 2019,]
LRI_incr_2019_UM    <- LRI_year_proportion_income_last_prop_UM[LRI_year_proportion_income_last_prop_UM$Year == 2019,]
LRI_incr_2019_GBD_1 <- LRI_year_proportion_GBD_region_prop_1  [LRI_year_proportion_GBD_region_prop_1  $Year == 2019,]
LRI_incr_2019_GBD_2 <- LRI_year_proportion_GBD_region_prop_2  [LRI_year_proportion_GBD_region_prop_2  $Year == 2019,]
LRI_incr_2019_GBD_3 <- LRI_year_proportion_GBD_region_prop_3  [LRI_year_proportion_GBD_region_prop_3  $Year == 2019,]
LRI_incr_2019_GBD_4 <- LRI_year_proportion_GBD_region_prop_4  [LRI_year_proportion_GBD_region_prop_4  $Year == 2019,]
LRI_incr_2019_GBD_5 <- LRI_year_proportion_GBD_region_prop_5  [LRI_year_proportion_GBD_region_prop_5  $Year == 2019,]
LRI_incr_2019_GBD_6 <- LRI_year_proportion_GBD_region_prop_6  [LRI_year_proportion_GBD_region_prop_6  $Year == 2019,]

colnames(LRI_incr_2019_L)    [2] <- "Category"
colnames(LRI_incr_2019_LM   )[2] <- "Category"
colnames(LRI_incr_2019_UM   )[2] <- "Category"
colnames(LRI_incr_2019_GBD_1)[2] <- "Category"
colnames(LRI_incr_2019_GBD_2)[2] <- "Category"
colnames(LRI_incr_2019_GBD_3)[2] <- "Category"
colnames(LRI_incr_2019_GBD_4)[2] <- "Category"
colnames(LRI_incr_2019_GBD_5)[2] <- "Category"
colnames(LRI_incr_2019_GBD_6)[2] <- "Category"

LRI_incr_2019$Category <- "Overall"  
LRI_incr_2019 <- LRI_incr_2019 %>% relocate(Category, .after = Year)

incr_df <- rbind(LRI_incr_2019, LRI_incr_2019_L, LRI_incr_2019_LM, LRI_incr_2019_UM, LRI_incr_2019_GBD_6, LRI_incr_2019_GBD_2, LRI_incr_2019_GBD_5, LRI_incr_2019_GBD_4, LRI_incr_2019_GBD_1, LRI_incr_2019_GBD_3)

incr_df$Category <- factor(
  incr_df$Category,
  levels = c( "Overall", "L", "LM", "UM", "Sub-Saharan Africa", "South Asia", "Southeast Asia, East Asia, and Oceania", "North Africa and Middle East", "Latin America and Caribbean", "Central Europe, Eastern Europe, and Central Asia"))
xpos <- seq_along(incr_df$Category)

incr_df$col <- NA
incr_df$col <- c("#0072B2", "#CC79A7", "#009E73", "#0072B2", "#D55E00", "#CC79A7", "#E69F00", "#A7B97F", "#009E73", "#0072B2")



write.csv(incr_df, "outputs/plot_df/increment_proportion.csv")





# ============================================================
# Export source data for Figure 4
# ============================================================

dir.create(
  "outputs/plot_objects",
  recursive = TRUE,
  showWarnings = FALSE
)


# ------------------------------------------------------------
# Helper: panels A-D, antibiotic access by income
# ------------------------------------------------------------

make_income_panel <- function(abx, panel) {
  
  income_groups <- c(
    "L"  = "Low income",
    "LM" = "Lower middle income",
    "UM" = "Upper middle income"
  )
  
  bind_rows(
    lapply(names(income_groups), function(inc) {
      
      df <- get(
        paste0("LRI_year_proportion_income_last_", inc)
      )
      
      df %>%
        transmute(
          Panel = panel,
          Measure = "Antibiotic treatment probability",
          Antibiotic_class = abx,
          Stratification = "Income",
          Group = income_groups[[inc]],
          Year = Year,
          Median = .data[[paste0(
            "proportion_LRI_", abx, "_med"
          )]],
          Lower = .data[[paste0(
            "proportion_LRI_", abx, "_lower"
          )]],
          Upper = .data[[paste0(
            "proportion_LRI_", abx, "_upper"
          )]]
        )
    })
  )
}


# Panels A-D
fig4_A <- make_income_panel("BSP", "A")
fig4_B <- make_income_panel("NSP", "B")
fig4_C <- make_income_panel("Cephalosporins", "C")
fig4_D <- make_income_panel("Aminoglycosides", "D")



# ------------------------------------------------------------
# Helper: panels E-H, antibiotic access by GBD region
# ------------------------------------------------------------

GBD_region_names <- c(
  "1" = "North Africa and Middle East",
  "2" = "Central Europe, Eastern Europe, and Central Asia",
  "3" = "Sub-Saharan Africa",
  "4" = "South Asia",
  "5" = "Latin America and Caribbean",
  "6" = "Southeast Asia, East Asia, and Oceania"
)

make_GBD_panel <- function(abx, panel) {
  
  bind_rows(
    lapply(1:6, function(i) {
      
      df <- get(
        paste0("LRI_year_proportion_GBD_region_", i)
      )
      
      df %>%
        transmute(
          Panel = panel,
          Measure = "Antibiotic treatment probability",
          Antibiotic_class = abx,
          Stratification = "GBD region",
          Group = GBD_region_names[as.character(i)],
          Year = Year,
          Median = .data[[paste0(
            "proportion_LRI_", abx, "_med"
          )]],
          Lower = .data[[paste0(
            "proportion_LRI_", abx, "_lower"
          )]],
          Upper = .data[[paste0(
            "proportion_LRI_", abx, "_upper"
          )]]
        )
    })
  )
}


# Panels E-H
fig4_E <- make_GBD_panel("BSP", "E")
fig4_F <- make_GBD_panel("NSP", "F")
fig4_G <- make_GBD_panel("Cephalosporins", "G")
fig4_H <- make_GBD_panel("Aminoglycosides", "H")



# ------------------------------------------------------------
# Helper: panels I-L,
# contribution to increase in antibiotic access in 2019
# ------------------------------------------------------------

make_increment_panel <- function(abx, panel) {
  
  incr_df %>%
    transmute(
      Panel = panel,
      Measure = paste0(
        "Proportion of increase in antibiotic access ",
        "attributable to class"
      ),
      Antibiotic_class = abx,
      
      Stratification = case_when(
        as.character(Category) == "Overall" ~ "Overall",
        as.character(Category) %in% c(
          "L", "LM", "UM"
        ) ~ "Income",
        TRUE ~ "GBD region"
      ),
      
      Group = case_when(
        as.character(Category) == "L" ~ "Low income",
        as.character(Category) == "LM" ~ "Lower middle income",
        as.character(Category) == "UM" ~ "Upper middle income",
        TRUE ~ as.character(Category)
      ),
      
      Year = Year,
      
      Median = .data[[paste0(
        "inc_LRI_", abx, "_med"
      )]],
      Lower = .data[[paste0(
        "inc_LRI_", abx, "_lower"
      )]],
      Upper = .data[[paste0(
        "inc_LRI_", abx, "_upper"
      )]]
    )
}


# Panels I-L
fig4_I <- make_increment_panel("BSP", "I")
fig4_J <- make_increment_panel("NSP", "J")
fig4_K <- make_increment_panel("Cephalosporins", "K")
fig4_L <- make_increment_panel("Aminoglycosides", "L")



# ------------------------------------------------------------
# Combine all Figure 4 panels
# ------------------------------------------------------------

Figure_4_source_data <- bind_rows(
  fig4_A,
  fig4_B,
  fig4_C,
  fig4_D,
  fig4_E,
  fig4_F,
  fig4_G,
  fig4_H,
  fig4_I,
  fig4_J,
  fig4_K,
  fig4_L
) %>%
  arrange(Panel, Group, Year)


# ------------------------------------------------------------
# Save
# ------------------------------------------------------------

write.csv(
  Figure_4_source_data,
  "outputs/plot_objects/Figure_4_source_data.csv",
  row.names = FALSE
)





  
pdf("figures/Figure_4.pdf", width = 7.5, height = 7)

layout(matrix(c(1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,
                0,  4,  4,  4,  4,  5,  5,  5,  5,  6,  6,  6,  6,  7,  7,  7,  7,  8,  8,  8,  8,
                0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
                2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,  2,
                0,  9,  9,  9,  9, 10, 10, 10, 10, 11, 11, 11, 11, 12, 12, 12, 12, 13, 13, 13, 13,
                0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
                3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,
                0, 14, 14, 14, 14, 14, 15, 15, 15, 15, 15, 16, 16, 16, 16, 16, 17, 17, 17, 17, 17 
                ), nrow = 8, byrow = TRUE), 
       widths = c(0.75 ,1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1), 
       heights = c(0.1,1,0.03,0.1,1,0.03,0.1,1.6))



#### 1. (title 1) ####

par(mar = c(0, 0, 0, 0))
plot.new()
text(x = -0.035, y = 0.5, labels = "Antibiotic access by class and income", 
     adj = 0, cex = 0.7, font = 2, xpd = TRUE)


#### 2. (title 2) ####

par(mar = c(0, 0, 0, 0))
plot.new()
text(x = -0.035, y = 0.5, labels = "Antibiotic access by class and GBD region", 
     adj = 0, cex = 0.7, font = 2, xpd = TRUE)

#### 3. (title 3) ####

par(mar = c(0, 0, 0, 0))
plot.new()
text(x = -0.035, y = 0.5, labels = "Class-specific contribution to the increase in antibiotic access, 2003\u20132019", 
     adj = 0, cex = 0.7, font = 2, xpd = TRUE)

#### 5. (Fig.A) #####
par(mar = c(2.5, 2.5, 2, 1))  
par(mgp = c(3, 0.3, 0))     
par(tck = -0.025)
par(lwd = 1)


xmin <- 2003; xmax <- 2019
ymin <- log10(0.001) ; ymax <- log10(1)


plot(proportion_LRI_BSP_med ~ Year, 
     data =  LRI_year_proportion_income_last_L,
     xlab = "", 
     ylab = "",
     pch = 16, col = rgb(0.8, 0.8, 0.8, 0),
     xlim = c(xmin,xmax),
     ylim = c(ymin,ymax),
     lwd = "1",
     axes=F,ann=F,
     main = "")



yticks <- log10(c(0.001, 0.01, 0.1, 1))   
ylabels <- c(0.001, 0.01, 0.1, 1)         
minor_labels <- as.vector(sapply(c(0.001, 0.01, 0.1), function(v) v * 2:9))

axis(1, at = seq(xmin, xmax, by = 4), labels = FALSE, las = 2, pos = ymin)
text(x = seq(xmin, xmax, by = 4), y = ymin-(ymax-ymin)*0.08,  
     labels = format(round(seq(xmin, xmax, by = (xmax-xmin)/6*0+4), digits = 0), nsmall = 0), xpd = TRUE, cex = 0.65)

axis(2,  at = yticks, labels = FALSE, las = 2, pos = xmin)
axis(2, at = log10(minor_labels), labels = FALSE, las = 2, pos = xmin, tck = -0.015, lwd = 0.7)
text(y = yticks, x = xmin-0.04*(xmax-xmin),  
     labels = ylabels, srt = 0, adj = 1, xpd = TRUE, cex = 0.65)

lines(log10(proportion_LRI_BSP_med) ~ Year, data = LRI_year_proportion_income_last_L,  col = "#CC79A7", lwd = 1.5)
lines(log10(proportion_LRI_BSP_med) ~ Year, data = LRI_year_proportion_income_last_LM, col = "#009E73", lwd = 1.5)
lines(log10(proportion_LRI_BSP_med) ~ Year, data = LRI_year_proportion_income_last_UM, col = "#0072B2", lwd = 1.5)



# Add shaded area for confidence interval
polygon(c(LRI_year_proportion_income_last_L$Year, rev(LRI_year_proportion_income_last_L$Year)),
        c(log10(LRI_year_proportion_income_last_L$proportion_LRI_BSP_lower), rev(log10(LRI_year_proportion_income_last_L$proportion_LRI_BSP_upper))),
        col = rgb(204/255, 121/255, 167/255, 0.25), border = NA)

polygon(c(LRI_year_proportion_income_last_LM$Year, rev(LRI_year_proportion_income_last_LM$Year)),
        c(log10(LRI_year_proportion_income_last_LM$proportion_LRI_BSP_lower), rev(log10(LRI_year_proportion_income_last_LM$proportion_LRI_BSP_upper))),
        col = rgb(0/255, 158/255, 115/255, 0.25), border = NA)

polygon(c(LRI_year_proportion_income_last_UM$Year, rev(LRI_year_proportion_income_last_UM$Year)),
        c(log10(LRI_year_proportion_income_last_UM$proportion_LRI_BSP_lower), rev(log10(LRI_year_proportion_income_last_UM$proportion_LRI_BSP_upper))),
        col = rgb(0/255, 114/255, 178/255, 0.25), border = NA)
mtext(expression("Antibiotic treatment probability"), side = 2, line = 2.5, cex = 0.5)  
mtext("Year", side = 1, line = 1, cex = 0.5)
mtext("A. BSP", at = xmin - 0.54*(xmax-xmin), side = 3, line = 0.5, adj = -0.1, cex = 0.5, font = 1)



#### 6. (Fig.B) #####
par(mar = c(2.5, 2.5, 2, 1))
par(mgp = c(3, 0.3, 0))     
par(tck = -0.025)
par(lwd = 1)

xmin <- 2003; xmax <- 2019
ymin <- log10(0.001) ; ymax <- log10(1)

plot(proportion_LRI_NSP_med ~ Year, 
     data =  LRI_year_proportion_income_last_L,
     xlab = "", 
     ylab = "",
     pch = 16, col = rgb(0.8, 0.8, 0.8, 0),
     xlim = c(xmin,xmax),
     ylim = c(ymin,ymax),
     lwd = "1",
     axes=F,ann=F,
     main = "")

yticks <- log10(c(0.001, 0.01, 0.1, 1))   
ylabels <- c(0.001, 0.01, 0.1, 1)         
minor_labels <- as.vector(sapply(c(0.001, 0.01, 0.1), function(v) v * 2:9))

axis(1, at = seq(xmin, xmax, by = 4), labels = FALSE, las = 2, pos = ymin)
text(x = seq(xmin, xmax, by = 4), y = ymin-(ymax-ymin)*0.08,  
     labels = format(round(seq(xmin, xmax, by = (xmax-xmin)/6*0+4), digits = 0), nsmall = 0), xpd = TRUE, cex = 0.65)

axis(2,  at = yticks, labels = FALSE, las = 2, pos = xmin)
axis(2, at = log10(minor_labels), labels = FALSE, las = 2, pos = xmin, tck = -0.015, lwd = 0.7)
text(y = yticks, x = xmin-0.04*(xmax-xmin),  
     labels = ylabels, srt = 0, adj = 1, xpd = TRUE, cex = 0.65)  

lines(log10(proportion_LRI_NSP_med) ~ Year, data = LRI_year_proportion_income_last_L,  col = "#CC79A7", lwd = 1.5)
lines(log10(proportion_LRI_NSP_med) ~ Year, data = LRI_year_proportion_income_last_LM, col = "#009E73", lwd = 1.5)
lines(log10(proportion_LRI_NSP_med) ~ Year, data = LRI_year_proportion_income_last_UM, col = "#0072B2", lwd = 1.5)

# Add shaded area for confidence interval
polygon(c(LRI_year_proportion_income_last_L$Year, rev(LRI_year_proportion_income_last_L$Year)),
        c(log10(LRI_year_proportion_income_last_L$proportion_LRI_NSP_lower), rev(log10(LRI_year_proportion_income_last_L$proportion_LRI_NSP_upper))),
        col = rgb(204/255, 121/255, 167/255, 0.25), border = NA)

polygon(c(LRI_year_proportion_income_last_LM$Year, rev(LRI_year_proportion_income_last_LM$Year)),
        c(log10(LRI_year_proportion_income_last_LM$proportion_LRI_NSP_lower), rev(log10(LRI_year_proportion_income_last_LM$proportion_LRI_NSP_upper))),
        col = rgb(0/255, 158/255, 115/255, 0.25), border = NA)

polygon(c(LRI_year_proportion_income_last_UM$Year, rev(LRI_year_proportion_income_last_UM$Year)),
        c(log10(LRI_year_proportion_income_last_UM$proportion_LRI_NSP_lower), rev(log10(LRI_year_proportion_income_last_UM$proportion_LRI_NSP_upper))),
        col = rgb(0/255, 114/255, 178/255, 0.25), border = NA)
mtext("Year", side = 1, line = 1, cex = 0.5)
mtext("B. NSP", at = xmin - 0.54*(xmax-xmin), side = 3, line = 0.5, adj = -0.1, cex = 0.5, font = 1)



#### 7. (Fig.C) #####
par(mar = c(2.5, 2.5, 2, 1))
par(mgp = c(3, 0.3, 0))     
par(tck = -0.025)
par(lwd = 1)

xmin <- 2003; xmax <- 2019
ymin <- log10(0.001) ; ymax <- log10(1)

plot(proportion_LRI_Cephalosporins_med ~ Year, 
     data =  LRI_year_proportion_income_last_L,
     xlab = "", 
     ylab = "",
     pch = 16, col = rgb(0.8, 0.8, 0.8, 0),
     xlim = c(xmin,xmax),
     ylim = c(ymin,ymax),
     lwd = "1",
     axes=F,ann=F,
     main = "")

yticks <- log10(c(0.001, 0.01, 0.1, 1))   
ylabels <- c(0.001, 0.01, 0.1, 1)         
minor_labels <- as.vector(sapply(c(0.001, 0.01, 0.1), function(v) v * 2:9))

axis(1, at = seq(xmin, xmax, by = 4), labels = FALSE, las = 2, pos = ymin)
text(x = seq(xmin, xmax, by = 4), y = ymin-(ymax-ymin)*0.08,  
     labels = format(round(seq(xmin, xmax, by = (xmax-xmin)/6*0+4), digits = 0), nsmall = 0), xpd = TRUE, cex = 0.65)

axis(2,  at = yticks, labels = FALSE, las = 2, pos = xmin)
axis(2, at = log10(minor_labels), labels = FALSE, las = 2, pos = xmin, tck = -0.015, lwd = 0.7)
text(y = yticks, x = xmin-0.04*(xmax-xmin),  
     labels = ylabels, srt = 0, adj = 1, xpd = TRUE, cex = 0.65)

lines(log10(proportion_LRI_Cephalosporins_med) ~ Year, data = LRI_year_proportion_income_last_L,  col = "#CC79A7", lwd = 1.5)
lines(log10(proportion_LRI_Cephalosporins_med) ~ Year, data = LRI_year_proportion_income_last_LM, col = "#009E73", lwd = 1.5)
lines(log10(proportion_LRI_Cephalosporins_med) ~ Year, data = LRI_year_proportion_income_last_UM, col = "#0072B2", lwd = 1.5)

# Add shaded area for confidence interval
polygon(c(LRI_year_proportion_income_last_L$Year, rev(LRI_year_proportion_income_last_L$Year)),
        c(log10(LRI_year_proportion_income_last_L$proportion_LRI_Cephalosporins_lower), rev(log10(LRI_year_proportion_income_last_L$proportion_LRI_Cephalosporins_upper))),
        col = rgb(204/255, 121/255, 167/255, 0.25), border = NA)

polygon(c(LRI_year_proportion_income_last_LM$Year, rev(LRI_year_proportion_income_last_LM$Year)),
        c(log10(LRI_year_proportion_income_last_LM$proportion_LRI_Cephalosporins_lower), rev(log10(LRI_year_proportion_income_last_LM$proportion_LRI_Cephalosporins_upper))),
        col = rgb(0/255, 158/255, 115/255, 0.25), border = NA)

polygon(c(LRI_year_proportion_income_last_UM$Year, rev(LRI_year_proportion_income_last_UM$Year)),
        c(log10(LRI_year_proportion_income_last_UM$proportion_LRI_Cephalosporins_lower), rev(log10(LRI_year_proportion_income_last_UM$proportion_LRI_Cephalosporins_upper))),
        col = rgb(0/255, 114/255, 178/255, 0.25), border = NA)
mtext("Year", side = 1, line = 1, cex = 0.5)
mtext("C. Cephalosporins", at = xmin - 0.54*(xmax-xmin), side = 3, line = 0.5, adj = -0.1, cex = 0.5, font = 1)


#### 8. (Fig.D) #####
par(mar = c(2.5, 2.5, 2, 1))
par(mgp = c(3, 0.3, 0))     
par(tck = -0.025)
par(lwd = 1)

xmin <- 2003; xmax <- 2019
ymin <- log10(0.001) ; ymax <- log10(1)

plot(proportion_LRI_Aminoglycosides_med ~ Year, 
     data =  LRI_year_proportion_income_last_L,
     xlab = "", 
     ylab = "",
     pch = 16, col = rgb(0.8, 0.8, 0.8, 0),
     xlim = c(xmin,xmax),
     ylim = c(ymin,ymax),
     lwd = "1",
     axes=F,ann=F,
     main = "")

yticks <- log10(c(0.001, 0.01, 0.1, 1))   
ylabels <- c(0.001, 0.01, 0.1, 1)         
minor_labels <- as.vector(sapply(c(0.001, 0.01, 0.1), function(v) v * 2:9))

axis(1, at = seq(xmin, xmax, by = 4), labels = FALSE, las = 2, pos = ymin)
text(x = seq(xmin, xmax, by = 4), y = ymin-(ymax-ymin)*0.08,  
     labels = format(round(seq(xmin, xmax, by = (xmax-xmin)/6*0+4), digits = 0), nsmall = 0), xpd = TRUE, cex = 0.65)

axis(2,  at = yticks, labels = FALSE, las = 2, pos = xmin)
axis(2, at = log10(minor_labels), labels = FALSE, las = 2, pos = xmin, tck = -0.015, lwd = 0.7)
text(y = yticks, x = xmin-0.04*(xmax-xmin),  
     labels = ylabels, srt = 0, adj = 1, xpd = TRUE, cex = 0.65)

lines(log10(proportion_LRI_Aminoglycosides_med) ~ Year, data = LRI_year_proportion_income_last_L,  col = "#CC79A7", lwd = 1.5)
lines(log10(proportion_LRI_Aminoglycosides_med) ~ Year, data = LRI_year_proportion_income_last_LM, col = "#009E73", lwd = 1.5)
lines(log10(proportion_LRI_Aminoglycosides_med) ~ Year, data = LRI_year_proportion_income_last_UM, col = "#0072B2", lwd = 1.5)

# Add shaded area for confidence interval
polygon(c(LRI_year_proportion_income_last_L$Year, rev(LRI_year_proportion_income_last_L$Year)),
        c(log10(LRI_year_proportion_income_last_L$proportion_LRI_Aminoglycosides_lower), rev(log10(LRI_year_proportion_income_last_L$proportion_LRI_Aminoglycosides_upper))),
        col = rgb(204/255, 121/255, 167/255, 0.25), border = NA)

polygon(c(LRI_year_proportion_income_last_LM$Year, rev(LRI_year_proportion_income_last_LM$Year)),
        c(log10(LRI_year_proportion_income_last_LM$proportion_LRI_Aminoglycosides_lower), rev(log10(LRI_year_proportion_income_last_LM$proportion_LRI_Aminoglycosides_upper))),
        col = rgb(0/255, 158/255, 115/255, 0.25), border = NA)

polygon(c(LRI_year_proportion_income_last_UM$Year, rev(LRI_year_proportion_income_last_UM$Year)),
        c(log10(LRI_year_proportion_income_last_UM$proportion_LRI_Aminoglycosides_lower), rev(log10(LRI_year_proportion_income_last_UM$proportion_LRI_Aminoglycosides_upper))),
        col = rgb(0/255, 114/255, 178/255, 0.25), border = NA)

mtext("D. Aminoglycosides", at = xmin - 0.54*(xmax-xmin), side = 3, line = 0.5, adj = -0.1, cex = 0.5, font = 1)


#### 9. (legend) ####
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


#### 10. (Fig.E) #####
par(mar = c(2.5, 2.5, 2, 1))
par(mgp = c(3, 0.3, 0))     
par(tck = -0.025)
par(lwd = 1)

xmin <- 2003; xmax <- 2019
ymin <- log10(0.001) ; ymax <- log10(1)

plot(LRI_year_proportion_GBD_region_1$proportion_LRI_BSP_med ~ LRI_year_proportion_GBD_region_1$Year, 
     data =  LRI_year_proportion_GBD_region_1,
     xlab = "", 
     ylab = "",
     pch = 16, col = rgb(0.8, 0.8, 0.8, 0),
     xlim = c(xmin,xmax),
     ylim = c(ymin,ymax),
     lwd = "1",
     axes=F,ann=F,
     main = "")

yticks <- log10(c(0.001, 0.01, 0.1, 1))   
ylabels <- c(0.001, 0.01, 0.1, 1)         
minor_labels <- as.vector(sapply(c(0.001, 0.01, 0.1), function(v) v * 2:9))

axis(1, at = seq(xmin, xmax, by = 4), labels = FALSE, las = 2, pos = ymin)
text(x = seq(xmin, xmax, by = 4), y = ymin-(ymax-ymin)*0.08,  
     labels = format(round(seq(xmin, xmax, by = (xmax-xmin)/6*0+4), digits = 0), nsmall = 0), xpd = TRUE, cex = 0.65)

axis(2,  at = yticks, labels = FALSE, las = 2, pos = xmin)
axis(2, at = log10(minor_labels), labels = FALSE, las = 2, pos = xmin, tck = -0.015, lwd = 0.7)
text(y = yticks, x = xmin-0.04*(xmax-xmin),  
     labels = ylabels, srt = 0, adj = 1, xpd = TRUE, cex = 0.65)

lines(log10(proportion_LRI_BSP_med) ~ Year, data = LRI_year_proportion_GBD_region_6, col = "#D55E00", lwd = 1.5)
polygon(c(LRI_year_proportion_GBD_region_6$Year, rev(LRI_year_proportion_GBD_region_6$Year)),
        c(log10(LRI_year_proportion_GBD_region_6$proportion_LRI_BSP_lower), rev(log10(LRI_year_proportion_GBD_region_6$proportion_LRI_BSP_upper))),
        col = rgb(204/255, 121/255, 167/255, 0.25), border = NA)

lines(log10(proportion_LRI_BSP_med) ~ Year, data = LRI_year_proportion_GBD_region_2, col = "#CC79A7", lwd = 1.5)
polygon(c(LRI_year_proportion_GBD_region_2$Year, rev(LRI_year_proportion_GBD_region_2$Year)),
        c(log10(LRI_year_proportion_GBD_region_2$proportion_LRI_BSP_lower), rev(log10(LRI_year_proportion_GBD_region_2$proportion_LRI_BSP_upper))),
        col = rgb(213/255, 94/255, 0/255, 0.25), border = NA)

lines(log10(proportion_LRI_BSP_med) ~ Year, data = LRI_year_proportion_GBD_region_5, col = "#E69F00", lwd = 1.5)
polygon(c(LRI_year_proportion_GBD_region_5$Year, rev(LRI_year_proportion_GBD_region_5$Year)),
        c(log10(LRI_year_proportion_GBD_region_5$proportion_LRI_BSP_lower), rev(log10(LRI_year_proportion_GBD_region_5$proportion_LRI_BSP_upper))),
        col = rgb(230/255, 159/255, 0/255, 0.25), border = NA)

lines(log10(proportion_LRI_BSP_med) ~ Year, data = LRI_year_proportion_GBD_region_4, col = "#A7B97F", lwd = 1.5)
polygon(c(LRI_year_proportion_GBD_region_4$Year, rev(LRI_year_proportion_GBD_region_4$Year)),
        c(log10(LRI_year_proportion_GBD_region_4$proportion_LRI_BSP_lower), rev(log10(LRI_year_proportion_GBD_region_4$proportion_LRI_BSP_upper))),
        col = rgb(0/255, 158/255, 115/255, 0.25), border = NA)

lines(log10(proportion_LRI_BSP_med) ~ Year, data = LRI_year_proportion_GBD_region_1, col = "#009E73", lwd = 1.5)
polygon(c(LRI_year_proportion_GBD_region_1$Year, rev(LRI_year_proportion_GBD_region_1$Year)),
        c(log10(LRI_year_proportion_GBD_region_1$proportion_LRI_BSP_lower), rev(log10(LRI_year_proportion_GBD_region_1$proportion_LRI_BSP_upper))),
        col = rgb(0/255, 158/255, 115/255, 0.25), border = NA)

lines(log10(proportion_LRI_BSP_med) ~ Year, data = LRI_year_proportion_GBD_region_3, col = "#0072B2", lwd = 1.5)
polygon(c(LRI_year_proportion_GBD_region_3$Year, rev(LRI_year_proportion_GBD_region_3$Year)),
        c(log10(LRI_year_proportion_GBD_region_3$proportion_LRI_BSP_lower), rev(log10(LRI_year_proportion_GBD_region_3$proportion_LRI_BSP_upper))),
        col = rgb(0/255, 114/255, 178/255, 0.25), border = NA)

mtext(expression("Antibiotic treatment probability"), side = 2, line = 2.5, cex = 0.5)  
mtext("Year", side = 1, line = 1, cex = 0.5)
mtext("E. BSP", at = xmin - 0.54*(xmax-xmin), side = 3, line = 0.5, adj = -0.1, cex = 0.5, font = 1)




















#### 11. (Fig.F) #####
par(mar = c(2.5, 2.5, 2, 1))  
par(mgp = c(3, 0.3, 0))     
par(tck = -0.025)
par(lwd = 1)


xmin <- 2003; xmax <- 2019
ymin <- log10(0.001) ; ymax <- log10(1)


plot(LRI_year_proportion_GBD_region_1$proportion_LRI_NSP_med ~ LRI_year_proportion_GBD_region_1$Year, 
     data =  LRI_year_proportion_GBD_region_1,
     xlab = "", 
     ylab = "",
     pch = 16, col = rgb(0.8, 0.8, 0.8, 0),
     xlim = c(xmin,xmax),
     ylim = c(ymin,ymax),
     lwd = "1",
     axes=F,ann=F,
     main = "")

yticks <- log10(c(0.001, 0.01, 0.1, 1))   
ylabels <- c(0.001, 0.01, 0.1, 1)         
minor_labels <- as.vector(sapply(c(0.001, 0.01, 0.1), function(v) v * 2:9))

axis(1, at = seq(xmin, xmax, by = 4), labels = FALSE, las = 2, pos = ymin)
text(x = seq(xmin, xmax, by = 4), y = ymin-(ymax-ymin)*0.08,  
     labels = format(round(seq(xmin, xmax, by = (xmax-xmin)/6*0+4), digits = 0), nsmall = 0), xpd = TRUE, cex = 0.65)

axis(2,  at = yticks, labels = FALSE, las = 2, pos = xmin)
axis(2, at = log10(minor_labels), labels = FALSE, las = 2, pos = xmin, tck = -0.015, lwd = 0.7)
text(y = yticks, x = xmin-0.04*(xmax-xmin),  
     labels = ylabels, srt = 0, adj = 1, xpd = TRUE, cex = 0.65)

lines(log10(proportion_LRI_NSP_med) ~ Year, data = LRI_year_proportion_GBD_region_6, col = "#D55E00", lwd = 1.5)
polygon(c(LRI_year_proportion_GBD_region_6$Year, rev(LRI_year_proportion_GBD_region_6$Year)),
        c(log10(LRI_year_proportion_GBD_region_6$proportion_LRI_NSP_lower), rev(log10(LRI_year_proportion_GBD_region_6$proportion_LRI_NSP_upper))),
        col = rgb(204/255, 121/255, 167/255, 0.25), border = NA)

lines(log10(proportion_LRI_NSP_med) ~ Year, data = LRI_year_proportion_GBD_region_2, col = "#CC79A7", lwd = 1.5)
polygon(c(LRI_year_proportion_GBD_region_2$Year, rev(LRI_year_proportion_GBD_region_2$Year)),
        c(log10(LRI_year_proportion_GBD_region_2$proportion_LRI_NSP_lower), rev(log10(LRI_year_proportion_GBD_region_2$proportion_LRI_NSP_upper))),
        col = rgb(213/255, 94/255, 0/255, 0.25), border = NA)

lines(log10(proportion_LRI_NSP_med) ~ Year, data = LRI_year_proportion_GBD_region_5, col = "#E69F00", lwd = 1.5)
polygon(c(LRI_year_proportion_GBD_region_5$Year, rev(LRI_year_proportion_GBD_region_5$Year)),
        c(log10(LRI_year_proportion_GBD_region_5$proportion_LRI_NSP_lower), rev(log10(LRI_year_proportion_GBD_region_5$proportion_LRI_NSP_upper))),
        col = rgb(230/255, 159/255, 0/255, 0.25), border = NA)

lines(log10(proportion_LRI_NSP_med) ~ Year, data = LRI_year_proportion_GBD_region_4, col = "#A7B97F", lwd = 1.5)
polygon(c(LRI_year_proportion_GBD_region_4$Year, rev(LRI_year_proportion_GBD_region_4$Year)),
        c(log10(LRI_year_proportion_GBD_region_4$proportion_LRI_NSP_lower), rev(log10(LRI_year_proportion_GBD_region_4$proportion_LRI_NSP_upper))),
        col = rgb(0/255, 158/255, 115/255, 0.25), border = NA)

lines(log10(proportion_LRI_NSP_med) ~ Year, data = LRI_year_proportion_GBD_region_1, col = "#009E73", lwd = 1.5)
polygon(c(LRI_year_proportion_GBD_region_1$Year, rev(LRI_year_proportion_GBD_region_1$Year)),
        c(log10(LRI_year_proportion_GBD_region_1$proportion_LRI_NSP_lower), rev(log10(LRI_year_proportion_GBD_region_1$proportion_LRI_NSP_upper))),
        col = rgb(0/255, 158/255, 115/255, 0.25), border = NA)

lines(log10(proportion_LRI_NSP_med) ~ Year, data = LRI_year_proportion_GBD_region_3, col = "#0072B2", lwd = 1.5)
polygon(c(LRI_year_proportion_GBD_region_3$Year, rev(LRI_year_proportion_GBD_region_3$Year)),
        c(log10(LRI_year_proportion_GBD_region_3$proportion_LRI_NSP_lower), rev(log10(LRI_year_proportion_GBD_region_3$proportion_LRI_NSP_upper))),
        col = rgb(0/255, 114/255, 178/255, 0.25), border = NA)

mtext("Year", side = 1, line = 1, cex = 0.5)
mtext("F. NSP", at = xmin - 0.54*(xmax-xmin), side = 3, line = 0.5, adj = -0.1, cex = 0.5, font = 1)



#### 12. (Fig.G) #####
par(mar = c(2.5, 2.5, 2, 1))  
par(mgp = c(3, 0.3, 0))     
par(tck = -0.025)
par(lwd = 1)


xmin <- 2003; xmax <- 2019
ymin <- log10(0.001) ; ymax <- log10(1)


plot(LRI_year_proportion_GBD_region_1$proportion_LRI_Cephalosporins_med ~ LRI_year_proportion_GBD_region_1$Year, 
     data =  LRI_year_proportion_GBD_region_1,
     xlab = "", 
     ylab = "",
     pch = 16, col = rgb(0.8, 0.8, 0.8, 0),
     xlim = c(xmin,xmax),
     ylim = c(ymin,ymax),
     lwd = "1",
     axes=F,ann=F,
     main = "")

yticks <- log10(c(0.001, 0.01, 0.1, 1))   
ylabels <- c(0.001, 0.01, 0.1, 1)         
minor_labels <- as.vector(sapply(c(0.001, 0.01, 0.1), function(v) v * 2:9))

axis(1, at = seq(xmin, xmax, by = 4), labels = FALSE, las = 2, pos = ymin)
text(x = seq(xmin, xmax, by = 4), y = ymin-(ymax-ymin)*0.08,  
     labels = format(round(seq(xmin, xmax, by = (xmax-xmin)/6*0+4), digits = 0), nsmall = 0), xpd = TRUE, cex = 0.65)

axis(2,  at = yticks, labels = FALSE, las = 2, pos = xmin)
axis(2, at = log10(minor_labels), labels = FALSE, las = 2, pos = xmin, tck = -0.015, lwd = 0.7)
text(y = yticks, x = xmin-0.04*(xmax-xmin),  
     labels = ylabels, srt = 0, adj = 1, xpd = TRUE, cex = 0.65)

lines(log10(proportion_LRI_Cephalosporins_med) ~ Year, data = LRI_year_proportion_GBD_region_6, col = "#D55E00", lwd = 1.5)
polygon(c(LRI_year_proportion_GBD_region_6$Year, rev(LRI_year_proportion_GBD_region_6$Year)),
        c(log10(LRI_year_proportion_GBD_region_6$proportion_LRI_Cephalosporins_lower), rev(log10(LRI_year_proportion_GBD_region_6$proportion_LRI_Cephalosporins_upper))),
        col = rgb(204/255, 121/255, 167/255, 0.25), border = NA)

lines(log10(proportion_LRI_Cephalosporins_med) ~ Year, data = LRI_year_proportion_GBD_region_2, col = "#CC79A7", lwd = 1.5)
polygon(c(LRI_year_proportion_GBD_region_2$Year, rev(LRI_year_proportion_GBD_region_2$Year)),
        c(log10(LRI_year_proportion_GBD_region_2$proportion_LRI_Cephalosporins_lower), rev(log10(LRI_year_proportion_GBD_region_2$proportion_LRI_Cephalosporins_upper))),
        col = rgb(213/255, 94/255, 0/255, 0.25), border = NA)

lines(log10(proportion_LRI_Cephalosporins_med) ~ Year, data = LRI_year_proportion_GBD_region_5, col = "#E69F00", lwd = 1.5)
polygon(c(LRI_year_proportion_GBD_region_5$Year, rev(LRI_year_proportion_GBD_region_5$Year)),
        c(log10(LRI_year_proportion_GBD_region_5$proportion_LRI_Cephalosporins_lower), rev(log10(LRI_year_proportion_GBD_region_5$proportion_LRI_Cephalosporins_upper))),
        col = rgb(230/255, 159/255, 0/255, 0.25), border = NA)

lines(log10(proportion_LRI_Cephalosporins_med) ~ Year, data = LRI_year_proportion_GBD_region_4, col = "#A7B97F", lwd = 1.5)
polygon(c(LRI_year_proportion_GBD_region_4$Year, rev(LRI_year_proportion_GBD_region_4$Year)),
        c(log10(LRI_year_proportion_GBD_region_4$proportion_LRI_Cephalosporins_lower), rev(log10(LRI_year_proportion_GBD_region_4$proportion_LRI_Cephalosporins_upper))),
        col = rgb(0/255, 158/255, 115/255, 0.25), border = NA)

lines(log10(proportion_LRI_Cephalosporins_med) ~ Year, data = LRI_year_proportion_GBD_region_1, col = "#009E73", lwd = 1.5)
polygon(c(LRI_year_proportion_GBD_region_1$Year, rev(LRI_year_proportion_GBD_region_1$Year)),
        c(log10(LRI_year_proportion_GBD_region_1$proportion_LRI_Cephalosporins_lower), rev(log10(LRI_year_proportion_GBD_region_1$proportion_LRI_Cephalosporins_upper))),
        col = rgb(0/255, 158/255, 115/255, 0.25), border = NA)

lines(log10(proportion_LRI_Cephalosporins_med) ~ Year, data = LRI_year_proportion_GBD_region_3, col = "#0072B2", lwd = 1.5)
polygon(c(LRI_year_proportion_GBD_region_3$Year, rev(LRI_year_proportion_GBD_region_3$Year)),
        c(log10(LRI_year_proportion_GBD_region_3$proportion_LRI_Cephalosporins_lower), rev(log10(LRI_year_proportion_GBD_region_3$proportion_LRI_Cephalosporins_upper))),
        col = rgb(0/255, 114/255, 178/255, 0.25), border = NA)

mtext("Year", side = 1, line = 1, cex = 0.5)
mtext("G. Cephalosporins", at = xmin - 0.54*(xmax-xmin), side = 3, line = 0.5, adj = -0.1, cex = 0.5, font = 1)



#### 13. (Fig.H) #####
par(mar = c(2.5, 2.5, 2, 1))  
par(mgp = c(3, 0.3, 0))     
par(tck = -0.025)
par(lwd = 1)


xmin <- 2003; xmax <- 2019
ymin <- log10(0.001) ; ymax <- log10(1)


plot(LRI_year_proportion_GBD_region_1$proportion_LRI_Aminoglycosides_med ~ LRI_year_proportion_GBD_region_1$Year, 
     data =  LRI_year_proportion_GBD_region_1,
     xlab = "", 
     ylab = "",
     pch = 16, col = rgb(0.8, 0.8, 0.8, 0),
     xlim = c(xmin,xmax),
     ylim = c(ymin,ymax),
     lwd = "1",
     axes=F,ann=F,
     main = "")



yticks <- log10(c(0.001, 0.01, 0.1, 1))   
ylabels <- c(0.001, 0.01, 0.1, 1)         
minor_labels <- as.vector(sapply(c(0.001, 0.01, 0.1), function(v) v * 2:9))

axis(1, at = seq(xmin, xmax, by = 4), labels = FALSE, las = 2, pos = ymin)
text(x = seq(xmin, xmax, by = 4), y = ymin-(ymax-ymin)*0.08,  
     labels = format(round(seq(xmin, xmax, by = (xmax-xmin)/6*0+4), digits = 0), nsmall = 0), xpd = TRUE, cex = 0.65)

axis(2,  at = yticks, labels = FALSE, las = 2, pos = xmin)
axis(2, at = log10(minor_labels), labels = FALSE, las = 2, pos = xmin, tck = -0.015, lwd = 0.7)
text(y = yticks, x = xmin-0.04*(xmax-xmin),  
     labels = ylabels, srt = 0, adj = 1, xpd = TRUE, cex = 0.65)

lines(log10(proportion_LRI_Aminoglycosides_med) ~ Year, data = LRI_year_proportion_GBD_region_6, col = "#D55E00", lwd = 1.5)
polygon(c(LRI_year_proportion_GBD_region_6$Year, rev(LRI_year_proportion_GBD_region_6$Year)),
        c(log10(LRI_year_proportion_GBD_region_6$proportion_LRI_Aminoglycosides_lower), rev(log10(LRI_year_proportion_GBD_region_6$proportion_LRI_Aminoglycosides_upper))),
        col = rgb(204/255, 121/255, 167/255, 0.25), border = NA)

lines(log10(proportion_LRI_Aminoglycosides_med) ~ Year, data = LRI_year_proportion_GBD_region_2, col = "#CC79A7", lwd = 1.5)
polygon(c(LRI_year_proportion_GBD_region_2$Year, rev(LRI_year_proportion_GBD_region_2$Year)),
        c(log10(LRI_year_proportion_GBD_region_2$proportion_LRI_Aminoglycosides_lower), rev(log10(LRI_year_proportion_GBD_region_2$proportion_LRI_Aminoglycosides_upper))),
        col = rgb(213/255, 94/255, 0/255, 0.25), border = NA)

lines(log10(proportion_LRI_Aminoglycosides_med) ~ Year, data = LRI_year_proportion_GBD_region_5, col = "#E69F00", lwd = 1.5)
polygon(c(LRI_year_proportion_GBD_region_5$Year, rev(LRI_year_proportion_GBD_region_5$Year)),
        c(log10(LRI_year_proportion_GBD_region_5$proportion_LRI_Aminoglycosides_lower), rev(log10(LRI_year_proportion_GBD_region_5$proportion_LRI_Aminoglycosides_upper))),
        col = rgb(230/255, 159/255, 0/255, 0.25), border = NA)

lines(log10(proportion_LRI_Aminoglycosides_med) ~ Year, data = LRI_year_proportion_GBD_region_4, col = "#A7B97F", lwd = 1.5)
polygon(c(LRI_year_proportion_GBD_region_4$Year, rev(LRI_year_proportion_GBD_region_4$Year)),
        c(log10(LRI_year_proportion_GBD_region_4$proportion_LRI_Aminoglycosides_lower), rev(log10(LRI_year_proportion_GBD_region_4$proportion_LRI_Aminoglycosides_upper))),
        col = rgb(0/255, 158/255, 115/255, 0.25), border = NA)

lines(log10(proportion_LRI_Aminoglycosides_med) ~ Year, data = LRI_year_proportion_GBD_region_1, col = "#009E73", lwd = 1.5)
polygon(c(LRI_year_proportion_GBD_region_1$Year, rev(LRI_year_proportion_GBD_region_1$Year)),
        c(log10(LRI_year_proportion_GBD_region_1$proportion_LRI_Aminoglycosides_lower), rev(log10(LRI_year_proportion_GBD_region_1$proportion_LRI_Aminoglycosides_upper))),
        col = rgb(0/255, 158/255, 115/255, 0.25), border = NA)

lines(log10(proportion_LRI_Aminoglycosides_med) ~ Year, data = LRI_year_proportion_GBD_region_3, col = "#0072B2", lwd = 1.5)
polygon(c(LRI_year_proportion_GBD_region_3$Year, rev(LRI_year_proportion_GBD_region_3$Year)),
        c(log10(LRI_year_proportion_GBD_region_3$proportion_LRI_Aminoglycosides_lower), rev(log10(LRI_year_proportion_GBD_region_3$proportion_LRI_Aminoglycosides_upper))),
        col = rgb(0/255, 114/255, 178/255, 0.25), border = NA)

mtext("Year", side = 1, line = 1, cex = 0.5)
mtext("H. Aminoglycosides", at = xmin - 0.54*(xmax-xmin), side = 3, line = 0.5, adj = -0.1, cex = 0.5, font = 1)



#### 14. (legend) ####
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



#### 15. (Fig.I) #####
par(mar = c(9, 4, 2, 1))  
par(mgp = c(3, 0.3, 0))     
par(tck = -0.025)
par(lwd = 1)


xmin <- 0.5; xmax <- length(xpos) + 0.5
ymin <- log10(0.001) ; ymax <- log10(1)


plot(log10(inc_LRI_BSP_med) ~ xpos, 
     data =  incr_df,
     xlab = "", 
     ylab = "",
     pch = 16, col = rgb(0.8, 0.8, 0.8, 0),
     xlim = c(xmin,xmax),
     ylim = c(ymin,ymax),
     lwd = "1",
     axes=F,ann=F,
     main = "")


axis(1, at = seq(xmin-0.1, xmax, by = 3), labels = FALSE, las = 2, pos = ymin, tck = 0)
axis(1, at = xpos, labels = FALSE, pos = ymin)
text(x = xpos, y = ymin - 0.06 * (ymax - ymin), labels = incr_df$Category, srt = 45, adj = 1, xpd = TRUE, cex = 0.6)

yticks <- log10(c(0.001, 0.01, 0.1, 1))   
ylabels <- c(0.001, 0.01, 0.1, 1)         
minor_labels <- as.vector(sapply(c(0.001, 0.01, 0.1), function(v) v * 2:9))

axis(2,  at = yticks, labels = FALSE, las = 2, pos = xmin)
axis(2, at = log10(minor_labels), labels = FALSE, las = 2, pos = xmin, tck = -0.015, lwd = 0.7)
text(y = yticks, x = xmin-0.04*(xmax-xmin),  
     labels = ylabels, srt = 0, adj = 1, xpd = TRUE, cex = 0.65)

segments(x0 = xpos, x1 = xpos, y0 = log10(incr_df$inc_LRI_BSP_lower), y1 = log10(incr_df$inc_LRI_BSP_upper), col = incr_df$col, lwd = 1.2)
points(xpos, log10(incr_df$inc_LRI_BSP_med), pch = 16, col = incr_df$col, cex = 0.8)

mtext(expression("Proportion"), side = 2, line = 2.5, cex = 0.5)  
mtext("I. BSP", at = xmin - 0.45*(xmax-xmin), side = 3, line = 0.5, adj = 0, cex = 0.5, font = 1)

div_x <- c(1.5, 4.5)
segments(x0 = div_x, y0 = ymin+0.01, 
         x1 = div_x, y1 = ymax, 
         col = "grey70", lwd = 1)



#### 16. (Fig.J) #####
par(mar = c(9, 4, 2, 1))  
par(mgp = c(3, 0.3, 0))     
par(tck = -0.025)
par(lwd = 1)

xmin <- 0.5; xmax <- length(xpos) + 0.5
ymin <- log10(0.001) ; ymax <- log10(1)


plot(log10(inc_LRI_NSP_med) ~ xpos, 
     data =  incr_df,
     xlab = "", 
     ylab = "",
     pch = 16, col = rgb(0.8, 0.8, 0.8, 0),
     xlim = c(xmin,xmax),
     ylim = c(ymin,ymax),
     lwd = "1",
     axes=F,ann=F,
     main = "")

axis(1, at = seq(xmin-0.1, xmax, by = 3), labels = FALSE, las = 2, pos = ymin, tck = 0)
axis(1, at = xpos, labels = FALSE, pos = ymin)
text(x = xpos, y = ymin - 0.06 * (ymax - ymin), labels = incr_df$Category, srt = 45, adj = 1, xpd = TRUE, cex = 0.6)

yticks <- log10(c(0.001, 0.01, 0.1, 1))   
ylabels <- c(0.001, 0.01, 0.1, 1)         
minor_labels <- as.vector(sapply(c(0.001, 0.01, 0.1), function(v) v * 2:9))

axis(2,  at = yticks, labels = FALSE, las = 2, pos = xmin)
axis(2, at = log10(minor_labels), labels = FALSE, las = 2, pos = xmin, tck = -0.015, lwd = 0.7)
text(y = yticks, x = xmin-0.04*(xmax-xmin),  
     labels = ylabels, srt = 0, adj = 1, xpd = TRUE, cex = 0.65)

incr_df$inc_LRI_NSP_lower[c(2,10)] <- 0.001

points(xpos, log10(incr_df$inc_LRI_NSP_med), pch = 16, col = incr_df$col, cex = 0.8)

arrow_indices <- c(2, 10)
segment_indices <- c(1, 3, 4, 5, 6, 7, 8, 9)

segments(x0 = xpos[segment_indices], x1 = xpos[segment_indices], 
         y0 = log10(incr_df$inc_LRI_NSP_lower[segment_indices]), 
         y1 = log10(incr_df$inc_LRI_NSP_upper[segment_indices]), 
         col = incr_df$col[segment_indices], lwd = 1.2)

arrows(x0 = xpos[arrow_indices], x1 = xpos[arrow_indices], 
       y1 = log10(incr_df$inc_LRI_NSP_upper[arrow_indices]), 
       y0 = ymin, 
       code = 1,
       length = 0.04, col = incr_df$col[arrow_indices], lwd = 1.2)

mtext("J. NSP", at = xmin - 0.45*(xmax-xmin), side = 3, line = 0.5, adj = 0, cex = 0.5, font = 1)

div_x <- c(1.5, 4.5)
segments(x0 = div_x, y0 = ymin+0.01, 
         x1 = div_x, y1 = ymax, 
         col = "grey70", lwd = 1)



#### 17. (Fig.K) #####

par(mar = c(9, 4, 2, 1))  
par(mgp = c(3, 0.3, 0))     
par(tck = -0.025)
par(lwd = 1)

xmin <- 0.5; xmax <- length(xpos) + 0.5
ymin <- log10(0.001) ; ymax <- log10(1)

plot(log10(inc_LRI_Cephalosporins_med) ~ xpos, 
     data =  incr_df,
     xlab = "", 
     ylab = "",
     pch = 16, col = rgb(0.8, 0.8, 0.8, 0),
     xlim = c(xmin,xmax),
     ylim = c(ymin,ymax),
     lwd = "1",
     axes=F,ann=F,
     main = "")



axis(1, at = seq(xmin-0.1, xmax, by = 3), labels = FALSE, las = 2, pos = ymin, tck = 0)
axis(1, at = xpos, labels = FALSE, pos = ymin)
text(x = xpos, y = ymin - 0.06 * (ymax - ymin), labels = incr_df$Category, srt = 45, adj = 1, xpd = TRUE, cex = 0.6)

yticks <- log10(c(0.001, 0.01, 0.1, 1))   
ylabels <- c(0.001, 0.01, 0.1, 1)         
minor_labels <- as.vector(sapply(c(0.001, 0.01, 0.1), function(v) v * 2:9))

axis(2,  at = yticks, labels = FALSE, las = 2, pos = xmin)
axis(2, at = log10(minor_labels), labels = FALSE, las = 2, pos = xmin, tck = -0.015, lwd = 0.7)
text(y = yticks, x = xmin-0.04*(xmax-xmin),  
     labels = ylabels, srt = 0, adj = 1, xpd = TRUE, cex = 0.65)

segments(x0 = xpos, x1 = xpos, y0 = log10(incr_df$inc_LRI_Cephalosporins_lower), y1 = log10(incr_df$inc_LRI_Cephalosporins_upper), col = incr_df$col, lwd = 1.2)
points(xpos, log10(incr_df$inc_LRI_Cephalosporins_med), pch = 16, col = incr_df$col, cex = 0.8)

mtext("K. Cephalosporins", at = xmin - 0.45*(xmax-xmin), side = 3, line = 0.5, adj = 0, cex = 0.5, font = 1)

div_x <- c(1.5, 4.5)
segments(x0 = div_x, y0 = ymin+0.01, 
         x1 = div_x, y1 = ymax, 
         col = "grey70", lwd = 1)



#### 18. (Fig.L) #####

par(mar = c(9, 4, 2, 1))  
par(mgp = c(3, 0.3, 0))     
par(tck = -0.025)
par(lwd = 1)

xmin <- 0.5; xmax <- length(xpos) + 0.5
ymin <- log10(0.001) ; ymax <- log10(1)

plot(log10(inc_LRI_Aminoglycosides_med) ~ xpos, 
     data =  incr_df,
     xlab = "", 
     ylab = "",
     pch = 16, col = rgb(0.8, 0.8, 0.8, 0),
     xlim = c(xmin,xmax),
     ylim = c(ymin,ymax),
     lwd = "1",
     axes=F,ann=F,
     main = "")


axis(1, at = seq(xmin-0.1, xmax, by = 3), labels = FALSE, las = 2, pos = ymin, tck = 0)
axis(1, at = xpos, labels = FALSE, pos = ymin)
text(x = xpos, y = ymin - 0.06 * (ymax - ymin), labels = incr_df$Category, srt = 45, adj = 1, xpd = TRUE, cex = 0.6)

yticks <- log10(c(0.001, 0.01, 0.1, 1))   
ylabels <- c(0.001, 0.01, 0.1, 1)         
minor_labels <- as.vector(sapply(c(0.001, 0.01, 0.1), function(v) v * 2:9))

axis(2,  at = yticks, labels = FALSE, las = 2, pos = xmin)
axis(2, at = log10(minor_labels), labels = FALSE, las = 2, pos = xmin, tck = -0.015, lwd = 0.7)
text(y = yticks, x = xmin-0.04*(xmax-xmin),  
     labels = ylabels, srt = 0, adj = 1, xpd = TRUE, cex = 0.65)

segments(x0 = xpos, x1 = xpos, y0 = log10(incr_df$inc_LRI_Aminoglycosides_lower), y1 = log10(incr_df$inc_LRI_Aminoglycosides_upper), col = incr_df$col, lwd = 1.2)
points(xpos, log10(incr_df$inc_LRI_Aminoglycosides_med), pch = 16, col = incr_df$col, cex = 0.8)

mtext("L. Aminoglycosides", at = xmin - 0.45*(xmax-xmin), side = 3, line = 0.5, adj = 0, cex = 0.5, font = 1)

div_x <- c(1.5, 4.5)
segments(x0 = div_x, y0 = ymin+0.01, 
         x1 = div_x, y1 = ymax, 
         col = "grey70", lwd = 1)

dev.off()


