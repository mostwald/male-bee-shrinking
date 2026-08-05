# Bee body size change in California, 1900-2024
# Figures script - FINAL
# August 2026

# Run CA_size_analysis_FINAL.R first to generate all model objects.


library(tidyverse)
library(mgcv)
library(emmeans)
library(sf)
library(terra)
library(tigris)
library(rnaturalearth)
library(wesanderson)
library(patchwork)
library(ggridges)
library(ggtext)
library(ggh4x)
library(broom)

#### Colour palettes ####

genus_colours <- c(
  "Andrena" = "#E76F51", "Anthophora" = "#F28E2B", "Bombus" = "#90BE6D",
  "Hesperapis" = "#FCCA46", "Lasioglossum" = "#335F70", "Osmia" = "#995D81", "Xylocopa" = "#2A9D8F"
)

sex_colours <- c("female" = "#E76F51", "male" = "#335F70")
sex_labels  <- c("female" = "Female", "male" = "Male")

#### Shared spatial context ####

california <- states(cb = TRUE, resolution = "500k") %>% filter(NAME == "California")
western_states <- states(cb = TRUE, resolution = "500k") %>%
  filter(NAME %in% c("Oregon", "Nevada", "Arizona", "Idaho", "Utah"))
mexico <- ne_countries(country = "mexico", returnclass = "sf", scale = "medium")

best_grid  <- st_read("best_grid_150km_may2026.gpkg")
final_grid <- best_grid %>% filter(grid_id %in% levels(df$grid_id)) %>% st_transform(4326)
california_vect <- states(cb = TRUE, resolution = "500k") %>%
  filter(NAME == "California") %>%
  vect()

#### Figure 1: Map and sampling  ####

set.seed(42)
df_map <- df %>%
  filter(!is.na(decimalLatitude) & !is.na(decimalLongitude)) %>%
  mutate(lat_jitter = jitter(decimalLatitude, amount = 0.05), lon_jitter = jitter(decimalLongitude, amount = 0.05))

ggplot() +
  geom_sf(data = mexico, fill = "grey88", colour = "grey60", linewidth = 0.3) +
  geom_sf(data = western_states, fill = "grey88", colour = "grey60", linewidth = 0.3) +
  geom_sf(data = california, fill = "white", colour = "grey40", linewidth = 0.4) +
  geom_point(data = df_map, aes(x = lon_jitter, y = lat_jitter, colour = genus), size = 0.65
             , alpha = 0.35) +
  scale_colour_manual(values = genus_colours, name = "Genus") +
  coord_sf(xlim = c(-125.5, -113), ylim = c(31.5, 42.5), expand = FALSE, clip = "on") +
  guides(colour = guide_legend(override.aes = list(size = 3, alpha = 0.8))) +
  theme_void(base_size = 12) +
  theme(panel.background = element_rect(fill = "white", colour = NA), legend.position = "right",
        legend.title = element_text(face = "bold", size = 11), legend.text = element_text(face = "italic", size = 10),
        panel.border = element_rect(colour = "grey40", fill = NA, linewidth = 0.5),
        plot.background = element_rect(fill = "#d6eaf8", colour = NA))

ggplot() +
  geom_sf(data = mexico, fill = "grey88", colour = "grey60", linewidth = 0.3) +
  geom_sf(data = western_states, fill = "grey88", colour = "grey60", linewidth = 0.3) +
  geom_sf(data = california, fill = "white", colour = "grey40", linewidth = 0.4) +
  geom_sf(data = final_grid, fill = NA, colour = "grey30", linewidth = 0.4, linetype = "dashed") +
  geom_point(data = df_map, aes(x = lon_jitter, y = lat_jitter, colour = genus), size = 0.6, alpha = 0.4) +
  scale_colour_manual(values = genus_colours, name = "Genus") +
  coord_sf(xlim = c(-125.5, -113), ylim = c(31.5, 42.5), expand = FALSE, clip = "on") +
  guides(colour = guide_legend(override.aes = list(size = 3, alpha = 1))) +
  theme_void(base_size = 12) +
  theme(panel.background = element_rect(fill = "#d6eaf8", colour = NA), legend.position = "right",
        legend.title = element_text(face = "bold", size = 11), legend.text = element_text(face = "italic", size = 10),
        panel.border = element_rect(colour = "grey40", fill = NA, linewidth = 0.5),
        plot.background = element_rect(fill = "white", colour = NA))

genus_order <- c("Xylocopa", "Anthophora", "Bombus", "Osmia", "Andrena", "Lasioglossum", "Hesperapis")

df %>%
  mutate(genus = factor(genus, levels = rev(genus_order))) %>%
  ggplot(aes(x = developmentYear, y = genus, fill = genus, alpha = sex)) +
  geom_density_ridges(scale = 1.7, colour = "white", linewidth = 0.3, bandwidth = 5, position = "identity") +
  scale_fill_manual(values = genus_colours, guide = "none") +
  scale_alpha_manual(values = c("female" = 0.7, "male" = 0.4), labels = c("Female", "Male"), name = "Sex") +
  scale_x_continuous(breaks = seq(1900, 2020, by = 40)) +
  scale_y_discrete(position = "right") +
  labs(x = "Development year", y = NULL) +
  theme_classic(base_size = 12) +
  theme(axis.text.y = element_text(face = "italic", size = 11), legend.position = "right", axis.text.x = element_text(angle=40, hjust=1),
        plot.background = element_rect(fill = "white", colour = NA))

#### Figure 2: Body size over time ####

ggplot(df, aes(x = developmentYear, y = resid_itd, colour = sex, fill = sex)) +
  geom_point(size = 0.8, alpha = 0.3, position = position_jitter(width = 0.3, height = 0)) +
  geom_smooth(method = "lm", se = TRUE, linewidth = 0.7, alpha = 0.2) +
  scale_colour_manual(values = sex_colours, labels = sex_labels, name = "Sex") +
  scale_fill_manual(values = sex_colours, labels = sex_labels, name = "Sex") +
  scale_x_continuous(breaks = seq(1900, 2025, by = 25)) +
  scale_y_continuous(limits = c(-0.75, 0.75)) +
  labs(x = "Development year", y = "Relative body size\n(residual log ITD)") +
  theme_classic(base_size = 12) +
  theme(legend.position = "right", plot.background = element_rect(fill = "white", colour = NA))

ggplot(df, aes(x = sex, y = median_itd, fill = sex, colour = sex)) +
  geom_violin(alpha = 0.8, trim = FALSE) +
  scale_fill_manual(values = sex_colours, labels = sex_labels, guide = "none") +
  scale_colour_manual(values = sex_colours, labels = sex_labels, guide = "none") +
  scale_x_discrete(labels = sex_labels) +
  labs(x = NULL, y = "Body size (ITD, mm)") +
  theme_classic(base_size = 12) +
  theme(plot.background = element_rect(fill = "white", colour = NA))

year_start <- min(df$developmentYear); year_end <- max(df$developmentYear)
years_seq <- seq(year_start, year_end, length.out = 100)

pred_from_emtrends <- sex_trend_df %>%
  select(sex, slope_per_year, se_per_year) %>%
  crossing(years = years_seq) %>%
  mutate(fit_pct = (exp(slope_per_year * (years - year_start)) - 1) * 100,
         upper_pct = (exp((slope_per_year + 1.96 * se_per_year) * (years - year_start)) - 1) * 100,
         lower_pct = (exp((slope_per_year - 1.96 * se_per_year) * (years - year_start)) - 1) * 100)

ggplot(pred_from_emtrends, aes(x = years, y = fit_pct, colour = sex, fill = sex)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40", linewidth = 0.5) +
  geom_ribbon(aes(ymin = lower_pct, ymax = upper_pct), alpha = 0.2, colour = NA) +
  geom_line(linewidth = 1) +
  scale_colour_manual(values = sex_colours, labels = sex_labels, name = "Sex") +
  scale_fill_manual(values = sex_colours, labels = sex_labels, name = "Sex") +
  scale_x_continuous(breaks = seq(1900, 2025, by = 25)) +
  scale_y_continuous(labels = function(x) paste0(x, "%"), limits = c(-12, 8)) +
  labs(x = "Development year", y = "Predicted % change in\nbody size relative to 1900") +
  theme_classic(base_size = 12) +
  theme(legend.position = "right", plot.background = element_rect(fill = "white", colour = NA))

#### Figure 3: Local precipitation and temperature drive the temporal decline in male body size   ####

get_marginal_df <- function(model, predictor_name, seq_range = c(-2, 2)) {
  pred_seq <- seq(seq_range[1], seq_range[2], length.out = 100)
  b_year <- coef(model)["developmentYear_scaled"]
  b_int <- coef(model)[paste0("developmentYear_scaled:", predictor_name)]
  vcov_m <- vcov(model)
  se_year <- sqrt(vcov_m["developmentYear_scaled", "developmentYear_scaled"])
  se_int <- sqrt(vcov_m[paste0("developmentYear_scaled:", predictor_name), paste0("developmentYear_scaled:", predictor_name)])
  cov_yi <- vcov_m["developmentYear_scaled", paste0("developmentYear_scaled:", predictor_name)]
  se_slope <- sqrt(se_year^2 + pred_seq^2 * se_int^2 + 2 * pred_seq * cov_yi)
  slope <- (b_year + b_int * pred_seq) / scaling_params$year_sd
  slope_upper <- (b_year + b_int * pred_seq + 1.96 * se_slope) / scaling_params$year_sd
  slope_lower <- (b_year + b_int * pred_seq - 1.96 * se_slope) / scaling_params$year_sd
  tibble(predictor = predictor_name, x = pred_seq, slope = (exp(slope * 10) - 1) * 100,
         upper = (exp(slope_upper * 10) - 1) * 100, lower = (exp(slope_lower * 10) - 1) * 100)
}

all_marginal <- bind_rows(
  get_marginal_df(m_mechanism, "ppt_within_scaled"),
  get_marginal_df(m_mechanism, "tmean_within_scaled"),
  get_marginal_df(m_mechanism, "cropland_within_scaled"),
  get_marginal_df(m_mechanism, "urban_within_scaled")
) %>%
  mutate(predictor = factor(predictor,
           levels = c("ppt_within_scaled", "tmean_within_scaled", "cropland_within_scaled", "urban_within_scaled"),
           labels = c("Local precip.", "Local temp.", "Local cropland", "Local urban.")),
         sig = predictor %in% c("Local precip.", "Local temp."))  # BOTH significant now

local_env_df <- bind_rows(
  df_male %>% transmute(variable = "Local precip.", year = developmentYear, value = ppt_within),
  df_male %>% transmute(variable = "Local temp.", year = developmentYear, value = tmean_within),
  df_male %>% transmute(variable = "Local cropland", year = developmentYear, value = cropland_within),
  df_male %>% transmute(variable = "Local urban.", year = developmentYear, value = urban_within)
) %>%
  mutate(variable = factor(variable, levels = c("Local precip.", "Local temp.", "Local cropland", "Local urban.")))

make_local_trend_plot <- function(data, var_label, sig, show_x = FALSE, y_zoom = NULL) {
  m <- lm(value ~ year, data = data)
  pred_df <- tibble(year = seq(min(data$year), max(data$year), length.out = 100))
  preds <- predict(m, newdata = pred_df, se.fit = TRUE)
  pred_df <- pred_df %>% mutate(fit = preds$fit, lower = preds$fit - 1.96 * preds$se.fit, upper = preds$fit + 1.96 * preds$se.fit)
  p <- ggplot() +
    geom_point(data = data, aes(x = year, y = value), size = 0.8, alpha = 0.3, colour = "grey50") +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40", linewidth = 0.4) +
    geom_ribbon(data = pred_df, aes(x = year, ymin = lower, ymax = upper), fill = ifelse(sig, "#2A9D8F", "grey80"), alpha = 0.5) +
    geom_line(data = pred_df, aes(x = year, y = fit), colour = ifelse(sig, "#2A9D8F", "grey60"), linewidth = 1) +
    labs(title = var_label, x = if (show_x) "Development year" else NULL, y = NULL) +
    theme_classic(base_size = 11) +
    theme(plot.title = element_text(face = "bold", size = 10, hjust = 0.5),
          axis.text.x = if (show_x) element_text(angle = 45, hjust = 1) else element_blank(),
          axis.ticks.x = if (show_x) element_line() else element_blank(), plot.margin = margin(1, 2, 1, 2))
  if (!is.null(y_zoom)) p <- p + coord_cartesian(ylim = y_zoom)
  p
}

make_mech_plot <- function(data, var_label, sig, show_x = FALSE, y_limits = NULL) {
  p <- ggplot(data, aes(x = x, y = slope)) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40", linewidth = 0.4) +
    geom_ribbon(aes(ymin = lower, ymax = upper), fill = ifelse(sig, "#2A9D8F", "grey80"), alpha = 0.3) +
    geom_line(colour = ifelse(sig, "#2A9D8F", "grey60"), linewidth = 1) +
    scale_x_continuous(breaks = c(-1, 0, 1), limits = c(-1.6, 1.6), labels = c("Low", "Mean", "High")) +
    labs(title = var_label, x = if (show_x) "Deviation from cell mean" else NULL, y = NULL) +
    theme_classic(base_size = 11) +
    theme(plot.title = element_text(face = "bold", size = 10, hjust = 0.5),
          axis.text.x = if (show_x) element_text(angle = 45, hjust = 1) else element_blank(),
          axis.ticks.x = if (show_x) element_line() else element_blank(), plot.margin = margin(1, 2, 1, 2))
  if (!is.null(y_limits)) p <- p + scale_y_continuous(limits = y_limits)
  p
}

p_trend_precip <- make_local_trend_plot(local_env_df %>% filter(variable == "Local precip."), "Local precip.", sig = TRUE)
p_mech_precip  <- make_mech_plot(all_marginal %>% filter(predictor == "Local precip."), "Local precip.", sig = TRUE, y_limits = c(-1.3, 1.0))

p_trend_temp <- make_local_trend_plot(local_env_df %>% filter(variable == "Local temp."), "Local temp.", sig = TRUE)
p_mech_temp  <- make_mech_plot(all_marginal %>% filter(predictor == "Local temp."), "Local temp.", sig = TRUE, y_limits = c(-1.3, 1.0))

p_trend_crop <- make_local_trend_plot(local_env_df %>% filter(variable == "Local cropland"), "Local cropland", sig = TRUE)
p_mech_crop  <- make_mech_plot(all_marginal %>% filter(predictor == "Local cropland"), "Local cropland", sig = FALSE, y_limits = c(-1.3, 1.0))

p_trend_urban <- make_local_trend_plot(local_env_df %>% filter(variable == "Local urban."), "Local urban.", sig = TRUE, show_x = TRUE)
p_mech_urban  <- make_mech_plot(all_marginal %>% filter(predictor == "Local urban."), "Local urban.", sig = FALSE, show_x = TRUE, y_limits = c(-1.3, 1.0))

right_col <- p_trend_precip / p_trend_temp / p_trend_crop / p_trend_urban
left_col  <- p_mech_precip  / p_mech_temp  / p_mech_crop  / p_mech_urban
fig3a <- left_col | right_col
fig3a + plot_annotation(title = "Local environmental trends and effects on male body size",
                         theme = theme(plot.title = element_text(face = "bold", size = 13)))

# forest plot of interaction coefficients
interaction_terms <- c("developmentYear_scaled:ppt_within_scaled", "developmentYear_scaled:tmean_within_scaled", "developmentYear_scaled:cropland_within_scaled", "developmentYear_scaled:urban_within_scaled")

forest_df <- tidy(m_mechanism, parametric = TRUE, conf.int = TRUE) %>%
  filter(term %in% interaction_terms) %>%
  mutate(predictor = case_when(
           term == "developmentYear_scaled:ppt_within_scaled" ~ "Local precip.",
           term == "developmentYear_scaled:tmean_within_scaled" ~ "Local temp.",
           term == "developmentYear_scaled:cropland_within_scaled" ~ "Local cropland",
           term == "developmentYear_scaled:urban_within_scaled" ~ "Local urban."),
         predictor = fct_relevel(predictor, "Local precip.", "Local temp.", "Local cropland", "Local urban."),
         sig = p.value < 0.05)

fig3b <- ggplot(forest_df, aes(x = predictor, y = estimate)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50", linewidth = 0.4) +
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high, colour = sig), width = 0.15, linewidth = 1) +
  geom_point(aes(colour = sig), size = 2.5) +
  scale_colour_manual(values = c("TRUE" = "#2A9D8F", "FALSE" = "grey60"), guide = "none") +
  labs(x = NULL, y = "Year \u00d7 environment interaction estimate") +
  theme_classic(base_size = 12) +
  theme(axis.text.x = element_text(size = 13), plot.background = element_rect(fill = "white", colour = NA), plot.margin = margin(10, 15, 10, 10))
fig3b

#grid-cell precipitation trend map
specimen_counts <- df %>% mutate(grid_id = as.integer(as.character(grid_id))) %>% count(grid_id, name = "n_specimens")
grid_trends_sf <- best_grid %>% filter(grid_id %in% grid_trends$grid_id) %>%
  left_join(grid_trends, by = "grid_id") %>% left_join(specimen_counts, by = "grid_id")
california_aligned <- st_transform(california, st_crs(grid_trends_sf))
grid_trends_clipped <- st_intersection(grid_trends_sf, california_aligned)

male_points <- df %>% filter(sex == "male") %>% filter(!is.na(decimalLatitude), !is.na(decimalLongitude)) %>%
  st_as_sf(coords = c("decimalLongitude", "decimalLatitude"), crs = 4326) %>% st_transform(st_crs(grid_trends_clipped))

ggplot() +
  geom_sf(data = mexico, fill = "grey90", colour = "grey65", linewidth = 0.3) +
  geom_sf(data = western_states, fill = "grey90", colour = "grey65", linewidth = 0.3) +
  geom_sf(data = california, fill = "white", colour = NA) +
  geom_sf(data = grid_trends_clipped, aes(fill = slope), colour = "grey40", linewidth = 0.3, alpha = 0.75) +
  geom_sf(data = male_points, size = 1, alpha = 0.3, colour = "grey10") +
  geom_sf(data = california, fill = NA, colour = "grey15", linewidth = 0.5) +
  scale_fill_gradient2(low = "#E76F51", mid = "#FCCA46", high = "#2A9D8F", midpoint = 0, name = "Local precip.\ntrend (mm/yr)") +
  coord_sf(xlim = c(-125.5, -113), ylim = c(31.5, 42.5), expand = FALSE, clip = "on") +
  theme_void(base_size = 12) +
  theme(legend.position = "right", panel.background = element_rect(fill = "#d6eaf8", colour = NA),
        panel.border = element_rect(colour = "grey40", fill = NA, linewidth = 0.5), plot.background = element_rect(fill = "white", colour = NA))


#### Figure 4: Male size varies across environmental gradients      ####
# Uses the SIMPLIFIED driver model (temperature's sex-interaction excluded). Temperature shown in the dot plot for completeness only

# Panel a: sex-specific marginal effect curves (precip, cropland, urban only)
get_marginal_by_sex <- function(em_object, predictor_label, seq_range = c(-2, 2)) {
  em_df <- as.data.frame(summary(em_object, infer = TRUE))
  slope_col <- names(em_df)[grepl("trend$", names(em_df))]
  pred_seq <- seq(seq_range[1], seq_range[2], length.out = 100)
  map_dfr(seq_len(nrow(em_df)), function(i) {
    b <- em_df[[slope_col]][i]; se <- em_df$SE[i]
    tibble(sex = em_df$sex[i], x = pred_seq, slope = (exp(b * pred_seq) - 1) * 100,
           upper = (exp((b + 1.96 * se) * pred_seq) - 1) * 100, lower = (exp((b - 1.96 * se) * pred_seq) - 1) * 100)
  }) %>% mutate(predictor = predictor_label)
}

sex_marginal_all <- bind_rows(
  get_marginal_by_sex(em_ppt_sex,      "Precipitation"),
  get_marginal_by_sex(em_cropland_sex, "Cropland"),
  get_marginal_by_sex(em_urban_sex,    "Urbanisation")
) %>% mutate(predictor = factor(predictor, levels = c("Precipitation", "Cropland", "Urbanisation")))

ggplot(sex_marginal_all, aes(x = x, y = slope, colour = sex, fill = sex)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40", linewidth = 0.4) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2, colour = NA) +
  geom_line(linewidth = 1) +
  facet_wrap(~predictor, nrow = 1) +
  scale_colour_manual(values = sex_colours, labels = sex_labels, name = "Sex") +
  scale_fill_manual(values = sex_colours, labels = sex_labels, name = "Sex") +
  labs(x = "Deviation from mean (SD)", y = "Predicted % change\nin body size") +
  theme_classic(base_size = 11) +
  theme(strip.background = element_blank(), strip.text = element_text(face = "bold", size = 10),
        legend.position = "right", plot.background = element_rect(fill = "white", colour = NA))


# Panel b: sex-effects dot plot, precip/temp/cropland/urban. Temperature included for completeness (not part of final model).

sex_effects_with_temp %>%
  mutate(predictor = factor(predictor, levels = c("Urbanisation", "Cropland", "Temperature", "Precipitation")),
         sex = factor(sex, levels = c("female", "male"))) %>%
  ggplot(aes(x = percent_per_sd, y = predictor, colour = sex)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey40", linewidth = 0.5) +
  geom_errorbarh(aes(xmin = lower_percent, xmax = upper_percent), height = 0.2, linewidth = 0.8, position = position_dodge(width = 0.5)) +
  geom_point(size = 3, position = position_dodge(width = 0.5)) +
  scale_colour_manual(values = sex_colours, labels = sex_labels, name = "Sex") +
  scale_x_continuous(labels = function(x) paste0(x, "%")) +
  labs(x = "% change in body size\nper SD of predictor", y = NULL) +
  theme_classic(base_size = 12) +
  theme(legend.position = "right", panel.grid.major.y = element_blank(), axis.line.x = element_line(colour = "grey40"),
        plot.background = element_rect(fill = "white", colour = NA), axis.text.y = element_text(angle = 45, size = 13))

# Panel c: environmental maps

most_recent_ppt <- ppt_files[grepl("2024", ppt_files)]
r_ppt_recent <- rast(most_recent_ppt) %>% app(fun = sum)
r_ppt_recent_ca <- crop(r_ppt_recent, project(california_vect, crs(r_ppt_recent))) %>% mask(project(california_vect, crs(r_ppt_recent))) %>% project("EPSG:4326")
ppt_recent_df <- as.data.frame(r_ppt_recent_ca, xy = TRUE) %>%
  rename(value = 3) %>%
  filter(!is.na(value)) %>%
  mutate(log_value = log(value + 1))

most_recent_tmean <- tmean_files[grepl("2024", tmean_files)]
r_tmean_recent <- rast(most_recent_tmean) %>% app(fun = mean)
r_tmean_recent_ca <- crop(r_tmean_recent, project(california_vect, crs(r_tmean_recent))) %>% mask(project(california_vect, crs(r_tmean_recent))) %>% project("EPSG:4326")
tmean_recent_df <- as.data.frame(r_tmean_recent_ca, xy = TRUE) %>% rename(value = 3) %>% filter(!is.na(value))

most_recent_hyde <- hyde_cropland_files[grepl("2023", hyde_cropland_files)]
ca_hyde <- project(california_vect, crs(rast(most_recent_hyde)))
cropland_recent <- rast(most_recent_hyde) %>% crop(ca_hyde) %>% mask(ca_hyde) %>% project("EPSG:4326")
cropland_recent_df <- as.data.frame(cropland_recent, xy = TRUE) %>% rename(value = 3) %>% filter(!is.na(value))

most_recent_hisdac <- hisdac_files[grepl("2020", hisdac_files)]
california_hisdac <- project(california_vect, crs(rast(most_recent_hisdac)))
urban_recent <- rast(most_recent_hisdac) %>% crop(california_hisdac) %>% mask(california_hisdac) %>% project("EPSG:4326")
urban_recent_df <- as.data.frame(urban_recent, xy = TRUE) %>% rename(value = 3) %>% filter(!is.na(value)) %>% mutate(value_log = log1p(value))

map_theme <- function() {
  theme_void(base_size = 10) +
    theme(plot.title = element_text(face = "bold", size = 9, hjust = 0.5), legend.position = "right",
          legend.key.width = unit(0.3, "cm"), legend.title = element_text(size = 8), legend.text = element_text(size = 7),
          panel.background = element_rect(fill = "#d6eaf8", colour = NA),
          panel.border = element_rect(colour = "grey40", fill = NA, linewidth = 0.4),
          plot.background = element_rect(fill = "white", colour = NA), plot.margin = margin(1, 1, 1, 1))
}
map_background <- function() {
  list(geom_sf(data = mexico, fill = "grey88", colour = "grey60", linewidth = 0.3),
       geom_sf(data = western_states, fill = "grey88", colour = "grey60", linewidth = 0.3))
}

env_gradient <- function(name, trans = "identity", stops = c(0, 0.15, 0.55, 1)) {
  scale_fill_gradientn(
    colours = c("#F5EFE6", "#DD4D40", "#FCCA46", "#2A9D8F"),
    values  = stops, name = name, trans = trans, na.value = "grey90"
  )
}

p_map_temp <- ggplot() + map_background() +
  geom_raster(data = tmean_recent_df, aes(x = x, y = y, fill = value), alpha = 0.85) +
  geom_sf(data = california, fill = NA, colour = "grey30", linewidth = 0.4) +
  env_gradient("Temp.\n(\u00b0C)") +
  coord_sf(xlim = c(-125.5, -113), ylim = c(31.5, 42.5), expand = FALSE, clip = "on") +
  ggtitle("Temperature (2024)") + map_theme()

p_map_ppt <- ggplot() + map_background() +
  geom_raster(data = ppt_recent_df, aes(x = x, y = y, fill = log_value), alpha = 0.85) +
  geom_sf(data = california, fill = NA, colour = "grey30", linewidth = 0.4) +
  env_gradient("Precip.\n(log mm)") +
  coord_sf(xlim = c(-125.5, -113), ylim = c(31.5, 42.5), expand = FALSE, clip = "on") +
  ggtitle("Precipitation (2024)") + map_theme()

p_map_crop <- ggplot() + map_background() +
  geom_raster(data = cropland_recent_df, aes(x = x, y = y, fill = value), alpha = 0.85) +
  geom_sf(data = california, fill = NA, colour = "grey30", linewidth = 0.4) +
  env_gradient("Cropland\nintensity") +
  coord_sf(xlim = c(-125.5, -113), ylim = c(31.5, 42.5), expand = FALSE, clip = "on") +
  ggtitle("Precipitation (2024)") + map_theme()

p_map_crop <- ggplot() + map_background() +
  geom_raster(data = cropland_recent_df, aes(x = x, y = y, fill = value)) +
  geom_sf(data = california, fill = NA, colour = "grey30", linewidth = 0.4) +
  env_gradient("Cropland\nintensity") +
  ggtitle("Cropland intensity (2023)") + map_theme()

p_map_urban <- ggplot() + map_background() +
  geom_raster(data = urban_recent_df, aes(x = x, y = y, fill = value_log)) +
  geom_sf(data = california, fill = NA, colour = "grey30", linewidth = 0.4) +
  env_gradient("Urban.\n(log)") +
  coord_sf(xlim = c(-125.5, -113), ylim = c(31.5, 42.5), expand = FALSE, clip = "on") +
  ggtitle("Urbanisation intensity (2020)") + map_theme()

fig4c <- (p_map_ppt | p_map_temp) / (p_map_crop | p_map_urban)
fig4c





#### Supplementary figures ####

heatmap_df <- genus_effects %>%
  mutate(predictor = factor(predictor, levels = c("Temperature", "Precipitation", "Urbanisation", "Cropland")),
         genus = factor(genus, levels = c("Xylocopa", "Bombus", "Andrena", "Anthophora", "Osmia", "Lasioglossum", "Hesperapis")))

#### Rebuild genus_effects with p-values attached, for significance asterisks ####

em_temp_df     <- as.data.frame(summary(em_temp,     infer = TRUE)) %>% mutate(predictor = "Temperature")
em_ppt_df      <- as.data.frame(summary(em_ppt,      infer = TRUE)) %>% mutate(predictor = "Precipitation")
em_urban_df    <- as.data.frame(summary(em_urban,    infer = TRUE)) %>% mutate(predictor = "Urbanisation")
em_cropland_df <- as.data.frame(summary(em_cropland, infer = TRUE)) %>% mutate(predictor = "Cropland")

# Standardise column names before combining (each has a differently-named .trend column)
em_temp_df     <- em_temp_df     %>% rename(estimate = tmean_scaled.trend)
em_ppt_df      <- em_ppt_df      %>% rename(estimate = log_ppt_scaled.trend)
em_urban_df    <- em_urban_df    %>% rename(estimate = urban_scaled.trend)
em_cropland_df <- em_cropland_df %>% rename(estimate = cropland_scaled.trend)

genus_pvals <- bind_rows(em_temp_df, em_ppt_df, em_urban_df, em_cropland_df) %>%
  select(genus, predictor, p.value)

# Join p-values onto genus_effects (already has percent_per_sd, from back_transform())
heatmap_df <- genus_effects %>%
  left_join(genus_pvals, by = c("genus", "predictor")) %>%
  mutate(
    sig = case_when(p.value < 0.001 ~ "***", p.value < 0.01 ~ "**", p.value < 0.05 ~ "*", TRUE ~ ""),
    predictor = factor(predictor, levels = c("Temperature", "Precipitation", "Urbanisation", "Cropland")),
    genus = factor(genus, levels = c("Xylocopa", "Bombus", "Andrena", "Anthophora", "Osmia", "Lasioglossum", "Hesperapis"))
  )

ggplot(heatmap_df, aes(x = predictor, y = genus, fill = percent_per_sd)) +
  geom_tile(colour = "white", linewidth = 0.5) +
  geom_text(aes(label = sig), size = 5, vjust = 0.75, colour = "black") +
  scale_fill_gradient2(low = "#E76F51", mid = "#FCCA46", high = "#2A9D8F", midpoint = 0,
                       name = "% change\nper SD", limits = c(-5, 5), oob = scales::squish) +
  labs(x = NULL, y = NULL) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(size = 11, angle = 45, hjust = 1, colour = "grey20"),
    axis.text.y = element_text(size = 11, face = "italic"),
    panel.grid = element_blank(),
    legend.position = "right",
    plot.background = element_rect(fill = "white", colour = NA)
  )


trend_df_sex %>%
  ggplot(aes(x = reorder(genus, total_percent_change), y = total_percent_change, colour = sex)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40", linewidth = 0.5) +
  geom_errorbar(aes(ymin = total_lower, ymax = total_upper), width = 0.3, linewidth = 0.9, position = position_dodge(width = 0.5)) +
  geom_point(size = 3, position = position_dodge(width = 0.5)) +
  scale_colour_manual(values = sex_colours, labels = sex_labels, name = "Sex") +
  scale_y_continuous(labels = function(x) paste0(x, "%")) +
  labs(x = NULL, y = paste0("Change in body size (", min(df$developmentYear), "-", max(df$developmentYear), ")")) +
  coord_flip() +
  theme_classic(base_size = 12) +
  theme(strip.background = element_blank(), legend.position = "right", panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank(), axis.text.y = element_text(face = "italic", size = 11),
        axis.line.x = element_line(colour = "grey40"), plot.background = element_rect(fill = "white", colour = NA))

