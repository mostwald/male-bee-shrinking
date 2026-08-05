# Bee body size change in California, 1900-2024
# Analysis script
# August 2026

#### Setup ####

library(tidyverse)
library(mgcv)
library(emmeans)
library(car)
library(nnet)
library(zoo)
library(sf)
library(terra)
library(tigris)

#### Load raw data ####

df <- read.csv("grid_filtered_CA_size_data.csv")

df$sex            <- as.factor(df$sex)
df$genus          <- as.factor(df$genus)
df$grid_id        <- as.factor(df$grid_id)
df$scientificName <- as.factor(df$scientificName)


#### extract specimen-level precipitation and temperature ####

#Requires raw files from PRISM, HISDAC, HYDE datasets
ppt_files <- list.files(
  "PRISM_ppt", #change file name if needed.
  pattern = ".*ppt.*\\.bil$", recursive = TRUE, full.names = TRUE)

tmean_files <- list.files(
  "PRISM_tmean",
  pattern = "\\.bil$", recursive = TRUE, full.names = TRUE)

cat("File counts - ppt:", length(ppt_files), "| tmean:", length(tmean_files), "\n")

# Precipitation: specimen point, annual sum of 12 monthly rasters, CRS-corrected
extract_specimen_ppt <- function(df, ppt_files) {
  df$ppt_annual <- NA_real_
  df_valid <- df %>% filter(!is.na(developmentYear))
  for (yr in unique(df_valid$developmentYear)) {
    rows <- df_valid %>% filter(developmentYear == yr)
    coords_vect <- vect(rows, geom = c("decimalLongitude", "decimalLatitude"), crs = "EPSG:4326")
    yr_files <- ppt_files[grepl(paste0(yr, "\\d{2}"), ppt_files) & grepl("\\.bil$", ppt_files)]
    if (length(yr_files) > 0) {
      r_stack <- rast(yr_files)
      coords_proj <- project(coords_vect, crs(r_stack))
      vals <- extract(r_stack, coords_proj)
      annual_val <- rowSums(vals[, -1, drop = FALSE], na.rm = TRUE)
      df$ppt_annual[df$developmentYear == yr & !is.na(df$developmentYear)] <- annual_val
    }
  }
  df
}

# Temperature: specimen point, annual mean of 12 monthly rasters, CRS-corrected
extract_specimen_tmean <- function(df, tmean_files) {
  df$tmean_annual <- NA_real_
  df_valid <- df %>% filter(!is.na(developmentYear))
  for (yr in unique(df_valid$developmentYear)) {
    rows <- df_valid %>% filter(developmentYear == yr)
    coords_vect <- vect(rows, geom = c("decimalLongitude", "decimalLatitude"), crs = "EPSG:4326")
    yr_files <- tmean_files[grepl(paste0(yr, "[0-1][0-9]_bil\\.bil$"), tmean_files)]
    if (length(yr_files) > 0) {
      r_stack <- rast(yr_files)
      coords_proj <- project(coords_vect, crs(r_stack))
      vals <- extract(r_stack, coords_proj)
      annual_val <- rowMeans(vals[, -1, drop = FALSE], na.rm = TRUE)
      df$tmean_annual[df$developmentYear == yr & !is.na(df$developmentYear)] <- annual_val
    }
  }
  df
}

df <- extract_specimen_ppt(df, ppt_files)

df <- extract_specimen_tmean(df, tmean_files)

cat("ppt_annual range:  ", range(df$ppt_annual, na.rm = TRUE), "\n")
cat("tmean_annual range:", range(df$tmean_annual, na.rm = TRUE), "\n")
cat("NAs - ppt:", sum(is.na(df$ppt_annual)), "| tmean:", sum(is.na(df$tmean_annual)), "\n")



#### Document scaling parameters ####

scaling_params <- df %>%
  mutate(log_ppt = log(ppt_annual + 1)) %>%
  summarise(
    year_mean = mean(developmentYear, na.rm = TRUE), year_sd = sd(developmentYear, na.rm = TRUE),
    tmean_mean = mean(tmean_annual, na.rm = TRUE),    tmean_sd = sd(tmean_annual, na.rm = TRUE),
    tmax_mean = mean(tmax_annual, na.rm = TRUE),      tmax_sd = sd(tmax_annual, na.rm = TRUE),
    ppt_mean = mean(log_ppt, na.rm = TRUE),           ppt_sd = sd(log_ppt, na.rm = TRUE),
    urban_mean = mean(urban_fraction, na.rm = TRUE),  urban_sd = sd(urban_fraction, na.rm = TRUE),
    cropland_mean = mean(cropland_fraction, na.rm = TRUE), cropland_sd = sd(cropland_fraction, na.rm = TRUE)
  )
print(scaling_params)

#### Scale predictors ####

df <- df %>%
  mutate(
    log_ppt = log(ppt_annual + 1),
    developmentYear_scaled = as.numeric(scale(developmentYear)),
    tmean_scaled     = as.numeric(scale(tmean_annual)),
    tmax_scaled      = as.numeric(scale(tmax_annual)),
    log_ppt_scaled   = as.numeric(scale(log_ppt)),
    urban_scaled     = as.numeric(scale(urban_fraction)),
    cropland_scaled  = as.numeric(scale(cropland_fraction)),
    log_itd          = log(median_itd)
  )


#### Within-cell standardisation (grid-cell-levelclimate baselines) ####

best_grid <- st_read("best_grid_150km_may2026.gpkg")
analysis_grid_ids <- unique(as.integer(as.character(df$grid_id)))

r_check_prism   <- rast(ppt_files[1])
grid_vect_prism <- project(vect(best_grid), crs(r_check_prism))
grid_vect_prism_filtered <- grid_vect_prism[grid_vect_prism$grid_id %in% analysis_grid_ids, ]

extract_grid_annual_ppt <- function(files, grid_vect, years = 1900:2024, agg_fun = mean) {
  map_dfr(years, function(yr) {
    yr_files <- files[grepl(paste0(yr, "\\d{2}"), files) & grepl("\\.bil$", files)]
    if (length(yr_files) == 0) return(NULL)
    r <- rast(yr_files); r_agg <- app(r, fun = agg_fun)
    cell_vals <- extract(r_agg, grid_vect, fun = mean, na.rm = TRUE, ID = FALSE)
    tibble(grid_id = grid_vect$grid_id, year = yr, value = cell_vals[[1]])
  })
}

cat("\nExtracting grid-cell precipitation...\n")
grid_climate_ppt <- extract_grid_annual_ppt(ppt_files, grid_vect_prism_filtered, agg_fun = sum) %>%
  rename(ppt_annual = value)
write.csv(grid_climate_ppt, "grid_precip_full_record.csv", row.names = FALSE)

cat("Extracting grid-cell temperature...\n")
grid_climate_tmean <- extract_grid_annual_ppt(tmean_files, grid_vect_prism_filtered, agg_fun = mean) %>%
  rename(tmean_annual = value)
write.csv(grid_climate_tmean, "grid_tmean_full_record.csv", row.names = FALSE)

hyde_cropland_files <- list.files(
  "HYDE_cropland_data",
  pattern = "cropland\\d{4}AD\\.asc$", recursive = TRUE, full.names = TRUE)
r_check_hyde <- rast(hyde_cropland_files[1])
grid_vect_hyde <- project(vect(best_grid), crs(r_check_hyde))
grid_vect_hyde_filtered <- grid_vect_hyde[grid_vect_hyde$grid_id %in% analysis_grid_ids, ]

extract_grid_cropland <- function(files, grid_vect) {
  map_dfr(files, function(f) {
    yr <- as.integer(str_extract(basename(f), "\\d{4}"))
    r <- rast(f)
    cell_vals <- extract(r, grid_vect, fun = mean, na.rm = TRUE, ID = FALSE)
    tibble(grid_id = grid_vect$grid_id, year = yr, cropland = cell_vals[[1]])
  })
}
cat("Extracting grid-cell cropland...\n")
grid_climate_cropland <- extract_grid_cropland(hyde_cropland_files, grid_vect_hyde_filtered)
write.csv(grid_climate_cropland, "grid_cropland_full_record.csv", row.names = FALSE)

hisdac_files <- list.files(
  "HISDAC_urbanization_data",
  pattern = "\\.tif$", full.names = TRUE)
r_check_hisdac <- rast(hisdac_files[1])
grid_vect_hisdac <- project(vect(best_grid), crs(r_check_hisdac))
grid_vect_hisdac_filtered <- grid_vect_hisdac[grid_vect_hisdac$grid_id %in% analysis_grid_ids, ]

extract_grid_urban <- function(files, grid_vect) {
  map_dfr(files, function(f) {
    yr <- as.integer(str_extract(basename(f), "^\\d{4}"))
    r <- rast(f)
    cell_vals <- extract(r, grid_vect, fun = mean, na.rm = TRUE, ID = FALSE)
    tibble(grid_id = grid_vect$grid_id, year = yr, urban_bupr = cell_vals[[1]])
  })
}
cat("Extracting grid-cell urbanisation...\n")
grid_climate_urban_5yr <- extract_grid_urban(hisdac_files, grid_vect_hisdac_filtered)
grid_climate_urban <- grid_climate_urban_5yr %>%
  group_by(grid_id) %>% complete(year = 1900:2024) %>% arrange(grid_id, year) %>%
  mutate(urban_bupr = na.approx(urban_bupr, x = year, na.rm = FALSE)) %>%
  fill(urban_bupr, .direction = "downup") %>% ungroup()
write.csv(grid_climate_urban, "grid_urban_full_record.csv", row.names = FALSE)

# Baselines
ppt_baseline <- grid_climate_ppt %>% group_by(grid_id) %>%
  summarise(ppt_clim_mean = mean(log(ppt_annual + 1), na.rm = TRUE))
tmean_baseline <- grid_climate_tmean %>% group_by(grid_id) %>%
  summarise(tmean_clim_mean = mean(tmean_annual, na.rm = TRUE))
cropland_baseline <- grid_climate_cropland %>% group_by(grid_id) %>%
  summarise(cropland_clim_mean = mean(cropland, na.rm = TRUE))
urban_baseline <- grid_climate_urban %>% group_by(grid_id) %>%
  summarise(urban_clim_mean = mean(urban_bupr, na.rm = TRUE))

df <- df %>% select(-any_of(c("ppt_clim_mean", "tmean_clim_mean", "cropland_clim_mean", "urban_clim_mean")))

df <- df %>%
  mutate(grid_id_int = as.integer(as.character(grid_id))) %>%
  left_join(ppt_baseline, by = c("grid_id_int" = "grid_id")) %>%
  left_join(tmean_baseline, by = c("grid_id_int" = "grid_id")) %>%
  left_join(cropland_baseline, by = c("grid_id_int" = "grid_id")) %>%
  left_join(urban_baseline, by = c("grid_id_int" = "grid_id")) %>%
  select(-grid_id_int) %>%
  mutate(
    tmean_within = tmean_annual - tmean_clim_mean,
    ppt_within = log_ppt - ppt_clim_mean,
    cropland_within = cropland_fraction - cropland_clim_mean,
    urban_within = urban_fraction - urban_clim_mean
  ) %>%
  mutate(
    tmean_within_scaled = as.numeric(scale(tmean_within)),
    ppt_within_scaled = as.numeric(scale(ppt_within)),
    cropland_within_scaled = as.numeric(scale(cropland_within)),
    urban_within_scaled = as.numeric(scale(urban_within))
  )

df <- df %>% mutate(ppt_within_check = log_ppt - ppt_clim_mean)
stopifnot(isTRUE(all.equal(df$ppt_within, df$ppt_within_check)))
df <- df %>% select(-ppt_within_check)

df$resid_itd <- residuals(lm(log(median_itd) ~ genus, data = df))

grid_trends <- grid_climate_ppt %>% group_by(grid_id) %>%
  group_modify(~ {
    m <- lm(ppt_annual ~ year, data = .x)
    tibble(slope = coef(m)["year"], p_value = summary(m)$coefficients["year", "Pr(>|t|)"], n_years = nrow(.x))
  }) %>% ungroup()
print(grid_trends, n = Inf)

#### Random effects selection ####

m_re0 <- bam(log_itd ~ genus + sex + s(developmentYear_scaled), data = df, method = "fREML")
m_re1 <- bam(log_itd ~ genus + sex + s(developmentYear_scaled) + s(scientificName, bs = "re"), data = df, method = "fREML")
m_re2 <- bam(log_itd ~ genus + sex + s(developmentYear_scaled) + s(grid_id, bs = "re"), data = df, method = "fREML")
m_re3 <- bam(log_itd ~ genus + sex + s(developmentYear_scaled) + s(scientificName, bs = "re") + s(grid_id, bs = "re"), data = df, method = "fREML")
AIC(m_re0, m_re1, m_re2, m_re3)

#### Question 1: Has body size changed over time, and does this differ between sexes? ####

m_t1 <- bam(log_itd ~ genus + sex + developmentYear_scaled + s(scientificName, bs = "re") + s(grid_id, bs = "re"), data = df, method = "fREML")
m_t2 <- bam(log_itd ~ genus + sex + genus:developmentYear_scaled + s(scientificName, bs = "re") + s(grid_id, bs = "re"), data = df, method = "fREML")
m_t3 <- bam(log_itd ~ genus + sex + sex:developmentYear_scaled + s(scientificName, bs = "re") + s(grid_id, bs = "re"), data = df, method = "fREML")
m_t4 <- bam(log_itd ~ genus + sex + genus:developmentYear_scaled + sex:developmentYear_scaled + s(scientificName, bs = "re") + s(grid_id, bs = "re"), data = df, method = "fREML")
AIC(m_t1, m_t2, m_t3, m_t4)

best_temporal_model <- m_t4
summary(best_temporal_model); summary(best_temporal_model)$r.sq; anova(best_temporal_model)

year_sd <- scaling_params$year_sd
year_span <- max(df$developmentYear) - min(df$developmentYear)

em_year_sex <- emtrends(best_temporal_model, specs = ~ sex, var = "developmentYear_scaled")
summary(em_year_sex, infer = TRUE); pairs(em_year_sex)

sex_trend_df <- as.data.frame(em_year_sex) %>%
  mutate(slope_per_year = developmentYear_scaled.trend / year_sd, se_per_year = SE / year_sd,
         total_percent_change = (exp(slope_per_year * year_span) - 1) * 100,
         total_lower = (exp((slope_per_year - 1.96 * se_per_year) * year_span) - 1) * 100,
         total_upper = (exp((slope_per_year + 1.96 * se_per_year) * year_span) - 1) * 100)
print(sex_trend_df)

sex_coef <- coef(best_temporal_model)["sexmale"]
cat("Males are", round(abs((exp(sex_coef) - 1) * 100), 1), "% smaller than females\n")

em_year_genus_sex <- emtrends(best_temporal_model, specs = ~ genus + sex, var = "developmentYear_scaled")
trend_df_sex <- as.data.frame(summary(em_year_genus_sex, infer = TRUE)) %>%
  mutate(slope_per_year = developmentYear_scaled.trend / year_sd, se_per_year = SE / year_sd,
         total_percent_change = (exp(slope_per_year * year_span) - 1) * 100,
         total_lower = (exp((slope_per_year - 1.96 * se_per_year) * year_span) - 1) * 100,
         total_upper = (exp((slope_per_year + 1.96 * se_per_year) * year_span) - 1) * 100)

#### Question 2: What environmental factors modulate the male declines? ####

df_male <- df %>% filter(sex == "male")

m_mechanism <- bam(
  log_itd ~ genus +
    developmentYear_scaled * tmean_within_scaled +
    developmentYear_scaled * ppt_within_scaled +
    developmentYear_scaled * cropland_within_scaled +
    developmentYear_scaled * urban_within_scaled +
    s(scientificName, bs = "re") + s(grid_id, bs = "re"),
  data = df_male, method = "fREML"
)
summary(m_mechanism); anova(m_mechanism)

# VIF check
X_mech <- model.matrix(m_mechanism)
check_vif <- function(target_col, mat) {
  other_cols <- setdiff(colnames(mat), target_col)
  r_squared <- summary(lm(mat[, target_col] ~ mat[, other_cols]))$r.squared
  1 / (1 - r_squared)
}
mech_vif_terms <- c("developmentYear_scaled:tmean_within_scaled", "developmentYear_scaled:ppt_within_scaled",
                     "developmentYear_scaled:cropland_within_scaled", "developmentYear_scaled:urban_within_scaled")
sapply(mech_vif_terms, check_vif, mat = X_mech)

# Marginal predictions at -1 SD / mean / +1 SD, precipitation and temperature
emtrends_ppt <- emtrends(m_mechanism, specs = ~ ppt_within_scaled, var = "developmentYear_scaled",
                          at = list(ppt_within_scaled = c(-1, 0, 1)))
emtrends_ppt_df <- as.data.frame(summary(emtrends_ppt, infer = TRUE)) %>%
  mutate(slope_peryear = developmentYear_scaled.trend / scaling_params$year_sd,
         pct_per_decade = (exp(slope_peryear * 10) - 1) * 100,
         lower_pct = (exp((lower.CL / scaling_params$year_sd) * 10) - 1) * 100,
         upper_pct = (exp((upper.CL / scaling_params$year_sd) * 10) - 1) * 100)
print(emtrends_ppt_df)

emtrends_temp <- emtrends(m_mechanism, specs = ~ tmean_within_scaled, var = "developmentYear_scaled",
                            at = list(tmean_within_scaled = c(-1, 0, 1)))
emtrends_temp_df <- as.data.frame(summary(emtrends_temp, infer = TRUE)) %>%
  mutate(slope_peryear = developmentYear_scaled.trend / scaling_params$year_sd,
         pct_per_decade = (exp(slope_peryear * 10) - 1) * 100,
         lower_pct = (exp((lower.CL / scaling_params$year_sd) * 10) - 1) * 100,
         upper_pct = (exp((upper.CL / scaling_params$year_sd) * 10) - 1) * 100)
print(emtrends_temp_df)

#### Question 3: How do spatial environmental gradients predict sex-specific size variation? ####

m_w1 <- bam(log_itd ~ genus + sex + tmean_scaled + log_ppt_scaled + urban_scaled + cropland_scaled +
              s(developmentYear_scaled) + s(scientificName, bs = "re") + s(grid_id, bs = "re"), data = df, method = "fREML")
m_w2 <- bam(log_itd ~ genus + sex + genus:tmean_scaled + genus:log_ppt_scaled + genus:urban_scaled + genus:cropland_scaled +
              s(developmentYear_scaled) + s(scientificName, bs = "re") + s(grid_id, bs = "re"), data = df, method = "fREML")
m_w3 <- bam(log_itd ~ genus + sex + tmean_scaled + log_ppt_scaled + urban_scaled + cropland_scaled +
              sex:tmean_scaled + sex:log_ppt_scaled + sex:urban_scaled + sex:cropland_scaled +
              s(developmentYear_scaled) + s(scientificName, bs = "re") + s(grid_id, bs = "re"), data = df, method = "fREML")
m_w4 <- bam(log_itd ~ genus + sex + genus:tmean_scaled + genus:log_ppt_scaled + genus:urban_scaled + genus:cropland_scaled +
              sex:tmean_scaled + sex:log_ppt_scaled + sex:urban_scaled + sex:cropland_scaled +
              s(developmentYear_scaled) + s(scientificName, bs = "re") + s(grid_id, bs = "re"), data = df, method = "fREML")
AIC(m_w1, m_w2, m_w3, m_w4)

# VIF check on the full (m_w4) model before simplifying
X_driver <- model.matrix(m_w4)
driver_vif_terms <- c("sexmale:tmean_scaled", "sexmale:log_ppt_scaled", "sexmale:cropland_scaled", "sexmale:urban_scaled")
sapply(driver_vif_terms, check_vif, mat = X_driver)
# temperature VIF ~4.5, precipitation ~5.0, cropland ~3.4, urbanisation ~2.8

# Temperature's sex-interaction is non-significant (P=0.567) in m_w4; refit without it
m_driver_no_temp <- bam(log_itd ~ genus + sex + genus:tmean_scaled + genus:log_ppt_scaled + genus:urban_scaled + genus:cropland_scaled +
              sex:log_ppt_scaled + sex:urban_scaled + sex:cropland_scaled +
              s(developmentYear_scaled) + s(scientificName, bs = "re") + s(grid_id, bs = "re"), data = df, method = "fREML")
AIC(m_w4, m_driver_no_temp)
# ΔAIC ~1.5 slightly favors exclusion; precipitation's coefficient and significance are essentially unchanged - adopt as final model.

best_driver_model <- m_driver_no_temp
summary(best_driver_model); summary(best_driver_model)$r.sq; summary(best_driver_model)$dev.expl
anova(best_driver_model)

plot(fitted(best_driver_model), residuals(best_driver_model), main = "Driver model: residuals vs fitted")
qqPlot(residuals(best_driver_model))
plot(best_driver_model, select = 1, shade = TRUE, xlab = "Development year (scaled)",
     ylab = "Partial effect on log(ITD)", main = "Driver model: residual temporal smooth")
abline(h = 0, lty = 2, col = "red")

m_baseline <- bam(log_itd ~ genus + sex + s(developmentYear_scaled) + s(scientificName, bs = "re") + s(grid_id, bs = "re"), data = df, method = "fREML")
summary(m_baseline)$r.sq; summary(best_driver_model)$r.sq

back_transform <- function(em_object, var_name) {
  as.data.frame(em_object) %>% rename(slope = matches("\\.trend$"), SE = SE) %>%
    mutate(predictor = var_name, percent_per_sd = (exp(slope) - 1) * 100,
           lower_percent = (exp(slope - 1.96 * SE) - 1) * 100, upper_percent = (exp(slope + 1.96 * SE) - 1) * 100)
}

# Genus-level slopes (supplementary)
em_temp     <- emtrends(best_driver_model, specs = ~ genus, var = "tmean_scaled")
em_ppt      <- emtrends(best_driver_model, specs = ~ genus, var = "log_ppt_scaled")
em_urban    <- emtrends(best_driver_model, specs = ~ genus, var = "urban_scaled")
em_cropland <- emtrends(best_driver_model, specs = ~ genus, var = "cropland_scaled")
summary(em_temp, infer = TRUE); summary(em_ppt, infer = TRUE); summary(em_urban, infer = TRUE); summary(em_cropland, infer = TRUE)

# Sex-level slopes: precipitation, cropland, urbanisation ONLY (final model)
em_ppt_sex      <- emtrends(best_driver_model, specs = ~ sex, var = "log_ppt_scaled")
em_urban_sex    <- emtrends(best_driver_model, specs = ~ sex, var = "urban_scaled")
em_cropland_sex <- emtrends(best_driver_model, specs = ~ sex, var = "cropland_scaled")
summary(em_ppt_sex, infer = TRUE); summary(em_urban_sex, infer = TRUE); summary(em_cropland_sex, infer = TRUE)
pairs(em_ppt_sex); pairs(em_urban_sex); pairs(em_cropland_sex)

# Temperature sex-level slopes: from the FULL model (m_w4) for figure only
em_temp_sex_display <- emtrends(m_w4, specs = ~ sex, var = "tmean_scaled")
summary(em_temp_sex_display, infer = TRUE)
bt_temp_sex_display <- back_transform(em_temp_sex_display, "Temperature")

bt_ppt_sex      <- back_transform(em_ppt_sex,      "Precipitation")
bt_urban_sex    <- back_transform(em_urban_sex,    "Urbanisation")
bt_cropland_sex <- back_transform(em_cropland_sex, "Cropland")
sex_effects <- bind_rows(bt_ppt_sex, bt_urban_sex, bt_cropland_sex)
sex_effects_with_temp <- bind_rows(bt_temp_sex_display, sex_effects)  # for Fig 4 display only

bt_temp     <- back_transform(em_temp,     "Temperature")
bt_ppt      <- back_transform(em_ppt,      "Precipitation")
bt_urban    <- back_transform(em_urban,    "Urbanisation")
bt_cropland <- back_transform(em_cropland, "Cropland")
genus_effects <- bind_rows(bt_temp, bt_ppt, bt_urban, bt_cropland)

write.csv(genus_effects %>% select(predictor, genus, percent_per_sd, lower_percent, upper_percent),
          "genus_effects_driver_model.csv", row.names = FALSE)
write.csv(sex_effects %>% select(predictor, sex, percent_per_sd, lower_percent, upper_percent),
          "sex_effects_driver_model.csv", row.names = FALSE)

# Joint test: does sex modulate environmental response overall?
m_reduced_v2 <- bam(log_itd ~ genus + sex + genus:tmean_scaled + genus:log_ppt_scaled + genus:urban_scaled + genus:cropland_scaled +
                   s(developmentYear_scaled) + s(scientificName, bs = "re") + s(grid_id, bs = "re"), data = df, method = "fREML")
AIC(m_reduced_v2, best_driver_model)
anova(m_reduced_v2, best_driver_model, test = "Chisq")

#### Sensitivity analyses ####

composition_tests <- map_dfr(levels(df$genus), function(g) {
  genus_df <- df %>% filter(genus == g) %>% mutate(scientificName = droplevels(scientificName))
  if (nlevels(genus_df$scientificName) < 2) return(tibble(genus = g, LR = NA, p = NA, n_species = 1))
  m_null <- multinom(scientificName ~ sex, data = genus_df, trace = FALSE)
  m_year <- multinom(scientificName ~ sex + developmentYear_scaled, data = genus_df, trace = FALSE)
  test <- anova(m_null, m_year)
  tibble(genus = g, LR = test$`LR stat.`[2], p = test$`Pr(Chi)`[2], n_species = nlevels(genus_df$scientificName))
})
print(composition_tests)

#sensitivity analysis with species as a fixed effect
m_temporal_fixed <- gam(log_itd ~ scientificName + sex * developmentYear + s(grid_id, bs = "re"), data = df, method = "REML")
anova(m_temporal_fixed)

#### Local environmental trends (male specimens) ####

m_ppt_trend <- lm(ppt_within ~ developmentYear, data = df_male)
summary(m_ppt_trend)
cat("Local precipitation:", (exp(coef(m_ppt_trend)["developmentYear"] * 10) - 1) * 100,
    "% per decade, R2 =", summary(m_ppt_trend)$r.squared, "\n")

m_tmean_trend <- lm(tmean_within ~ developmentYear, data = df_male)
summary(m_tmean_trend)

m_cropland_trend <- lm(cropland_within ~ developmentYear, data = df_male)
summary(m_cropland_trend)

m_urban_trend <- lm(urban_within ~ developmentYear, data = df_male)
summary(m_urban_trend)


write.csv(df, "grid_filtered_CA_size_data_FINAL.csv", row.names = FALSE)

