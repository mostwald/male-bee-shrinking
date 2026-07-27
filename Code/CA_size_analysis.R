# Bee body size change in California, 1900-2024
# Analysis script
# May 2026

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

#### Load and prepare data ####
df <- read.csv("grid_filtered_CA_size_data.csv")

df$sex            <- as.factor(df$sex)
df$genus          <- as.factor(df$genus)
df$grid_id        <- as.factor(df$grid_id)
df$scientificName <- as.factor(df$scientificName)

#### Document scaling parameters ####

scaling_params <- df %>%
  mutate(log_ppt = log(ppt_annual + 1)) %>%
  summarise(
    year_mean     = mean(developmentYear,    na.rm = TRUE),
    year_sd       = sd(developmentYear,      na.rm = TRUE),
    tmean_mean    = mean(tmean_annual,        na.rm = TRUE),
    tmean_sd      = sd(tmean_annual,          na.rm = TRUE),
    tmax_mean     = mean(tmax_annual,         na.rm = TRUE),
    tmax_sd       = sd(tmax_annual,           na.rm = TRUE),
    ppt_mean      = mean(log_ppt,             na.rm = TRUE),
    ppt_sd        = sd(log_ppt,               na.rm = TRUE),
    urban_mean    = mean(urban_fraction,      na.rm = TRUE),
    urban_sd      = sd(urban_fraction,        na.rm = TRUE),
    cropland_mean = mean(cropland_fraction,   na.rm = TRUE),
    cropland_sd   = sd(cropland_fraction,     na.rm = TRUE)
  )
print(scaling_params)

#### Scale predictors ####

df <- df %>%
  mutate(
    log_ppt                = log(ppt_annual + 1),
    developmentYear_scaled = as.numeric(scale(developmentYear)),
    tmean_scaled           = as.numeric(scale(tmean_annual)),
    tmax_scaled            = as.numeric(scale(tmax_annual)),
    log_ppt_scaled         = as.numeric(scale(log_ppt)),
    urban_scaled           = as.numeric(scale(urban_fraction)),
    cropland_scaled        = as.numeric(scale(cropland_fraction)),
    log_itd                = log(median_itd)
  )

#### Within-cell standardisation ####

# For each grid cell, subtract the long-term cell mean from each environmental variable, then scale to unit variance globally. I
df <- df %>%
  group_by(grid_id) %>%
  mutate(
    tmean_within    = tmean_annual      - mean(tmean_annual,      na.rm = TRUE),
    tmax_within     = tmax_annual       - mean(tmax_annual,       na.rm = TRUE),
    ppt_within      = log_ppt           - mean(log_ppt,           na.rm = TRUE),
    urban_within    = urban_fraction    - mean(urban_fraction,    na.rm = TRUE),
    cropland_within = cropland_fraction - mean(cropland_fraction, na.rm = TRUE)
  ) %>%
  ungroup() %>%
  mutate(
    tmean_within_scaled    = as.numeric(scale(tmean_within)),
    tmax_within_scaled     = as.numeric(scale(tmax_within)),
    ppt_within_scaled      = as.numeric(scale(ppt_within)),
    urban_within_scaled    = as.numeric(scale(urban_within)),
    cropland_within_scaled = as.numeric(scale(cropland_within))
  )

df$resid_itd <- residuals(lm(log(median_itd) ~ genus, data = df))

#### Random effects selection ####

m_re0 <- bam(log_itd ~ genus + sex +
               s(developmentYear_scaled),
             data = df, method = "fREML")

m_re1 <- bam(log_itd ~ genus + sex +
               s(developmentYear_scaled) +
               s(scientificName, bs = "re"),
             data = df, method = "fREML")

m_re2 <- bam(log_itd ~ genus + sex +
               s(developmentYear_scaled) +
               s(grid_id, bs = "re"),
             data = df, method = "fREML")

m_re3 <- bam(log_itd ~ genus + sex +
               s(developmentYear_scaled) +
               s(scientificName, bs = "re") +
               s(grid_id, bs = "re"),
             data = df, method = "fREML")

AIC(m_re0, m_re1, m_re2, m_re3)
# m_re3 has lowest AIC — both random effects justified, carried forward to all models

#### Question 1: Has body size changed over time, and does this differ between sexes? ####

m_t1 <- bam(log_itd ~ genus + sex +
              developmentYear_scaled +
              s(scientificName, bs = "re") +
              s(grid_id, bs = "re"),
            data = df, method = "fREML")

m_t2 <- bam(log_itd ~ genus + sex +
              genus:developmentYear_scaled +
              s(scientificName, bs = "re") +
              s(grid_id, bs = "re"),
            data = df, method = "fREML")

m_t3 <- bam(log_itd ~ genus + sex +
              sex:developmentYear_scaled +
              s(scientificName, bs = "re") +
              s(grid_id, bs = "re"),
            data = df, method = "fREML")

m_t4 <- bam(log_itd ~ genus + sex +
              genus:developmentYear_scaled +
              sex:developmentYear_scaled +
              s(scientificName, bs = "re") +
              s(grid_id, bs = "re"),
            data = df, method = "fREML")

AIC(m_t1, m_t2, m_t3, m_t4)
# m_t4 has lowest AIC — genus x year and sex x year interactions both justified

best_temporal_model <- m_t4
summary(best_temporal_model)
summary(best_temporal_model)$r.sq
anova(best_temporal_model)

plot(fitted(best_temporal_model), residuals(best_temporal_model),
     main = "Temporal model: residuals vs fitted")
qqPlot(residuals(best_temporal_model))

# Per-sex temporal trends

year_sd   <- scaling_params$year_sd
year_span <- max(df$developmentYear) - min(df$developmentYear)

em_year_sex <- emtrends(best_temporal_model,
                        specs = ~ sex,
                        var = "developmentYear_scaled")
summary(em_year_sex, infer = TRUE)
pairs(em_year_sex)

sex_trend_df <- as.data.frame(em_year_sex) %>%
  mutate(
    slope_per_year       = developmentYear_scaled.trend / year_sd,
    se_per_year          = SE / year_sd,
    total_percent_change = (exp(slope_per_year * year_span) - 1) * 100,
    total_lower          = (exp((slope_per_year - 1.96 * se_per_year) * year_span) - 1) * 100,
    total_upper          = (exp((slope_per_year + 1.96 * se_per_year) * year_span) - 1) * 100
  )
print(sex_trend_df)
# Males: -4.95% (95% CI: -8.81 to -0.92%), p = 0.017
# Females: +1.44% (95% CI: -2.13 to +5.16%), p = 0.434
# Sex divergence: p = 0.002

#Overall sex size difference
sex_coef     <- coef(best_temporal_model)["sexmale"]
percent_diff <- (exp(sex_coef) - 1) * 100
cat("Males are", round(abs(percent_diff), 1), "% smaller than females\n")
# Males are 11.4% smaller than females

# Per-genus x sex temporal trends (supplementary)
em_year_genus_sex <- emtrends(best_temporal_model,
                              specs = ~ genus + sex,
                              var = "developmentYear_scaled")
summary(em_year_genus_sex, infer = TRUE)

trend_df_sex <- as.data.frame(summary(em_year_genus_sex,
                                      infer = TRUE)) %>%
  mutate(
    slope_per_year       = developmentYear_scaled.trend / year_sd,
    se_per_year          = SE / year_sd,
    total_percent_change = (exp(slope_per_year * year_span) - 1) * 100,
    total_lower          = (exp((slope_per_year - 1.96 * se_per_year) * year_span) - 1) * 100,
    total_upper          = (exp((slope_per_year + 1.96 * se_per_year) * year_span) - 1) * 100
  )

#### Question 2: Do environmental conditions modulate temporal trends in body size? ####  
# Fit for males only — females showed no significant temporal trend.
# Within-cell standardised predictors isolate temporal signal.

m_mechanism <- bam(
  log_itd ~ genus +
    developmentYear_scaled * tmean_within_scaled +
    developmentYear_scaled * ppt_within_scaled +
    developmentYear_scaled * cropland_within_scaled +
    developmentYear_scaled * urban_within_scaled +
    s(scientificName, bs = "re") +
    s(grid_id, bs = "re"),
  data = df %>% filter(sex == "male") %>% droplevels(),
  method = "fREML"
)

summary(m_mechanism)
anova(m_mechanism)

#### Question 3: How do environmental conditions predict body size variation across spatial and temporal gradients? ####
# Raw scaled predictors to capture both spatial and temporal variation. Temporal smooth included as nuisance term.

m_w1 <- bam(log_itd ~ genus + sex +
              tmean_scaled + log_ppt_scaled +
              urban_scaled + cropland_scaled +
              s(developmentYear_scaled) +
              s(scientificName, bs = "re") +
              s(grid_id, bs = "re"),
            data = df, method = "fREML")

m_w2 <- bam(log_itd ~ genus + sex +
              genus:tmean_scaled + genus:log_ppt_scaled +
              genus:urban_scaled + genus:cropland_scaled +
              s(developmentYear_scaled) +
              s(scientificName, bs = "re") +
              s(grid_id, bs = "re"),
            data = df, method = "fREML")

m_w3 <- bam(log_itd ~ genus + sex +
              tmean_scaled + log_ppt_scaled +
              urban_scaled + cropland_scaled +
              sex:tmean_scaled + sex:log_ppt_scaled +
              sex:urban_scaled + sex:cropland_scaled +
              s(developmentYear_scaled) +
              s(scientificName, bs = "re") +
              s(grid_id, bs = "re"),
            data = df, method = "fREML")

m_w4 <- bam(log_itd ~ genus + sex +
              genus:tmean_scaled + genus:log_ppt_scaled +
              genus:urban_scaled + genus:cropland_scaled +
              sex:tmean_scaled + sex:log_ppt_scaled +
              sex:urban_scaled + sex:cropland_scaled +
              s(developmentYear_scaled) +
              s(scientificName, bs = "re") +
              s(grid_id, bs = "re"),
            data = df, method = "fREML")

AIC(m_w1, m_w2, m_w3, m_w4)

best_driver_model <- m_w4
summary(best_driver_model)
summary(best_driver_model)$r.sq
summary(best_driver_model)$dev.expl
anova(best_driver_model)

# Diagnostics

plot(fitted(best_driver_model), residuals(best_driver_model),
     main = "Driver model: residuals vs fitted")
qqPlot(residuals(best_driver_model))
plot(best_driver_model, select = 1, shade = TRUE,
     xlab = "Development year (scaled)",
     ylab = "Partial effect on log(ITD)",
     main = "Driver model: residual temporal smooth")
abline(h = 0, lty = 2, col = "red")

back_transform <- function(em_object, var_name) {
  as.data.frame(em_object) %>%
    rename(slope = matches("\\.trend$"), SE = SE) %>%
    mutate(
      predictor      = var_name,
      percent_per_sd = (exp(slope) - 1) * 100,
      lower_percent  = (exp(slope - 1.96 * SE) - 1) * 100,
      upper_percent  = (exp(slope + 1.96 * SE) - 1) * 100
    )
}

# Per-genus slopes(supplementary)
em_temp     <- emtrends(best_driver_model, specs = ~ genus,
                        var = "tmean_scaled")
em_ppt      <- emtrends(best_driver_model, specs = ~ genus,
                        var = "log_ppt_scaled")
em_urban    <- emtrends(best_driver_model, specs = ~ genus,
                        var = "urban_scaled")
em_cropland <- emtrends(best_driver_model, specs = ~ genus,
                        var = "cropland_scaled")

summary(em_temp,     infer = TRUE)
summary(em_ppt,      infer = TRUE)
summary(em_urban,    infer = TRUE)
summary(em_cropland, infer = TRUE)

# Per-sex slopes
em_temp_sex     <- emtrends(best_driver_model, specs = ~ sex,
                            var = "tmean_scaled")
em_ppt_sex      <- emtrends(best_driver_model, specs = ~ sex,
                            var = "log_ppt_scaled")
em_urban_sex    <- emtrends(best_driver_model, specs = ~ sex,
                            var = "urban_scaled")
em_cropland_sex <- emtrends(best_driver_model, specs = ~ sex,
                            var = "cropland_scaled")

summary(em_temp_sex,     infer = TRUE)
summary(em_ppt_sex,      infer = TRUE)
summary(em_urban_sex,    infer = TRUE)
summary(em_cropland_sex, infer = TRUE)

pairs(em_temp_sex)
pairs(em_ppt_sex)
pairs(em_urban_sex)
pairs(em_cropland_sex)

# Back-transform
bt_temp     <- back_transform(em_temp,     "Temperature")
bt_ppt      <- back_transform(em_ppt,      "Precipitation")
bt_urban    <- back_transform(em_urban,    "Urbanisation")
bt_cropland <- back_transform(em_cropland, "Cropland")
genus_effects <- bind_rows(bt_temp, bt_ppt, bt_urban, bt_cropland)

bt_temp_sex     <- back_transform(em_temp_sex,     "Temperature")
bt_ppt_sex      <- back_transform(em_ppt_sex,      "Precipitation")
bt_urban_sex    <- back_transform(em_urban_sex,    "Urbanisation")
bt_cropland_sex <- back_transform(em_cropland_sex, "Cropland")
sex_effects <- bind_rows(bt_temp_sex, bt_ppt_sex,
                         bt_urban_sex, bt_cropland_sex)


#### Sensitivity analyses (supplementary) ####

# Test whether sampled species composition within each genus changed significantly over time.
composition_tests <- map_dfr(levels(df$genus), function(g) {
  genus_df <- df %>%
    filter(genus == g) %>%
    mutate(scientificName = droplevels(scientificName))
  
  if (nlevels(genus_df$scientificName) < 2) {
    return(tibble(genus = g, LR = NA, p = NA, n_species = 1))
  }
  
  m_null <- multinom(scientificName ~ sex,
                     data = genus_df, trace = FALSE)
  m_year <- multinom(scientificName ~ sex + developmentYear_scaled,
                     data = genus_df, trace = FALSE)
  
  test <- anova(m_null, m_year)
  
  tibble(
    genus     = g,
    LR        = test$`LR stat.`[2],
    p         = test$`Pr(Chi)`[2],
    n_species = nlevels(genus_df$scientificName)
  )
})
print(composition_tests)
# All genera show significant compositional change over time

# Is it genuine shrinkage, or just species turnover?
m_temporal_fixed <- gam(
  log_itd ~ scientificName + sex * developmentYear + s(grid_id, bs = "re"),
  data = df,
  method = "REML"
)
anova(m_temporal_fixed)




#### California environmental trends ####

# Load environmental data

env_ca <- read.csv("california_environmental_trends.csv") %>%
  mutate(variable = factor(variable,
                           levels = c("Temperature (°C)",
                                      "Precipitation (mm)",
                                      "Urbanisation intensity",
                                      "Cropland intensity")))

ppt_ca <- env_ca %>%
  filter(variable == "Precipitation (mm)") %>%
  select(year, value)

# Pixel-wise precipitation variability rasters
ppt_var_trend <- rast("ppt_var_trend_pixelwise.tif")
ppt_abs_resid <- rast("ppt_abs_resid_stack.tif")

# Statewide trend tests
trend_tests <- env_ca %>%
  group_by(variable) %>%
  group_modify(~ {
    m <- lm(value ~ year, data = .x)
    tibble(
      slope             = coef(m)["year"],
      se                = summary(m)$coefficients["year", "Std. Error"],
      t_value           = summary(m)$coefficients["year", "t value"],
      p_value           = summary(m)$coefficients["year", "Pr(>|t|)"],
      r_squared         = summary(m)$r.squared,
      change_per_decade = coef(m)["year"] * 10
    )
  }) %>%
  ungroup()
print(trend_tests)

env_ca %>%
  group_by(variable) %>%
  summarise(
    year_start   = min(year),
    year_end     = max(year),
    value_start  = value[which.min(year)],
    value_end    = value[which.max(year)],
    total_change = value_end - value_start,
    pct_change   = (value_end - value_start) / value_start * 100
  )

#Pixel-wise precipitation variability
ppt_annual <- ppt_ca %>%
  mutate(abs_dev = abs(value - mean(value, na.rm = TRUE)))

# Summarise pixel-wise variability trend
var_trend_vals <- values(ppt_var_trend, na.rm = TRUE)
t_result       <- t.test(var_trend_vals)

cat("Mean pixel-wise variability trend:",
    round(mean(var_trend_vals), 5), "mm/year\n")
cat("% pixels positive:",
    round(mean(var_trend_vals > 0) * 100, 1), "%\n")
print(t_result)

