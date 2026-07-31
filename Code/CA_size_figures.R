# Bee body size change in California, 1900-2024
# Figures script
# May 2026

# RunCA_size_analysis.R first to generate all model objects

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
library(patchwork)

#### Colour palettes ####

genus_colours <- c(
  "Andrena"      = "#E76F51",
  "Anthophora"   = "#F28E2B",
  "Bombus"       = "#90BE6D",
  "Hesperapis"   = "#FCCA46",
  "Lasioglossum" = "#335F70",
  "Osmia"        = "#995D81",
  "Xylocopa"     = "#2A9D8F"
)


sex_colours <- c("female" = "#E76F51", "male" = "#335F70")
sex_labels  <- c("female" = "Female", "male" = "Male")

#### Figure 1: Map and sampling ####

#### Figure 1a: California specimen map ####

california <- states(cb = TRUE, resolution = "500k") %>%
  filter(NAME == "California")

western_states <- states(cb = TRUE, resolution = "500k") %>%
  filter(NAME %in% c("Oregon", "Nevada", "Arizona", "Idaho", "Utah"))

mexico <- ne_countries(country = "mexico", returnclass = "sf",
                       scale = "medium")

best_grid <- st_read("best_grid_150km_may2026.gpkg")
final_grid <- best_grid %>%
  filter(grid_id %in% levels(df$grid_id)) %>%
  st_transform(4326)

set.seed(42)
df_map <- df %>%
  filter(!is.na(decimalLatitude) & !is.na(decimalLongitude)) %>%
  mutate(
    lat_jitter = jitter(decimalLatitude,  amount = 0.05),
    lon_jitter = jitter(decimalLongitude, amount = 0.05)
  )

ggplot() +
  geom_sf(data = mexico,
          fill = "grey88", colour = "grey60", linewidth = 0.3) +
  geom_sf(data = western_states,
          fill = "grey88", colour = "grey60", linewidth = 0.3) +
  geom_sf(data = california,
          fill = "white", colour = "grey40", linewidth = 0.4) +
  geom_point(data = df_map,
             aes(x = lon_jitter, y = lat_jitter, colour = genus),
             size = 1, alpha = 0.5) +
  scale_colour_manual(values = genus_colours, name = "Genus") +
  coord_sf(xlim = c(-125.5, -113), ylim = c(31.5, 42.5),
           expand = FALSE, clip = "on") +
  guides(colour = guide_legend(override.aes = list(size = 3,
                                                   alpha = 1))) +
  theme_void(base_size = 12) +
  theme(
    panel.background = element_rect(fill = "#d6eaf8", colour = NA),
    legend.position  = "right",
    legend.title     = element_text(face = "bold", size = 11),
    legend.text      = element_text(face = "italic", size = 10),
    panel.border     = element_rect(colour = "grey40", fill = NA,
                                    linewidth = 0.5),
    plot.background  = element_rect(fill = "white", colour = NA)
  )


#sup figure with grid cells
ggplot() +
  geom_sf(data = mexico,
          fill = "grey88", colour = "grey60", linewidth = 0.3) +
  geom_sf(data = western_states,
          fill = "grey88", colour = "grey60", linewidth = 0.3) +
  geom_sf(data = california,
          fill = "white", colour = "grey40", linewidth = 0.4) +
  geom_sf(data = final_grid,
          fill = NA, colour = "grey30", linewidth = 0.4,
          linetype = "dashed") +
  geom_point(data = df_map,
             aes(x = lon_jitter, y = lat_jitter, colour = genus),
             size = 1, alpha = 0.5) +
  scale_colour_manual(values = genus_colours, name = "Genus") +
  coord_sf(xlim = c(-125.5, -113), ylim = c(31.5, 42.5),
           expand = FALSE, clip = "on") +
  guides(colour = guide_legend(override.aes = list(size = 3,
                                                   alpha = 1))) +
  theme_void(base_size = 12) +
  theme(
    panel.background = element_rect(fill = "#d6eaf8", colour = NA),
    legend.position  = "right",
    legend.title     = element_text(face = "bold", size = 11),
    legend.text      = element_text(face = "italic", size = 10),
    panel.border     = element_rect(colour = "grey40", fill = NA,
                                    linewidth = 0.5),
    plot.background  = element_rect(fill = "white", colour = NA))
    
    
#ridge plot

# reorder
genus_order <- c("Xylocopa", "Anthophora", "Bombus", "Osmia",
                 "Andrena", "Lasioglossum", "Hesperapis")

df %>%
  mutate(genus = factor(genus, levels = rev(genus_order))) %>%
  ggplot(aes(x = developmentYear,
             y = genus,
             fill = genus,
             alpha = sex)) +
  geom_density_ridges(scale = 1.7,
                      colour = "white", linewidth = 0.3,
                      bandwidth = 5,
                      position = "identity") +
  scale_fill_manual(values = genus_colours, guide = "none") +
  scale_alpha_manual(values = c("female" = 0.9, "male" = 0.4),
                     labels = c("Female", "Male"),
                     name = "Sex") +
  scale_x_continuous(breaks = seq(1900, 2020, by = 40)) +
  labs(x = "Development year", y = NULL) +
  theme_classic(base_size = 12) +
  theme(
    axis.text.y     = element_text(face = "italic", size = 11),
    legend.position = "right",
    plot.background = element_rect(fill = "white", colour = NA)
  )

#### Figure 2: Body size over time ####

ggplot(df, aes(x = developmentYear, y = resid_itd,
               colour = sex, fill = sex)) +
  geom_point(size = 0.8, alpha = 0.3,
             position = position_jitter(width = 0.3, height = 0)) +
  geom_smooth(method = "lm", se = TRUE, linewidth = 0.7, alpha = 0.2) +
  scale_colour_manual(values = sex_colours, labels = sex_labels,
                      name = "Sex") +
  scale_fill_manual(values = sex_colours, labels = sex_labels,
                    name = "Sex") +
  scale_x_continuous(breaks = seq(1900, 2025, by = 25)) +
  scale_y_continuous(limits = c(-0.75, 0.75)) +
  labs(x = "Development year",
       y = "Relative body size\n(residual log ITD)") +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "right",
    plot.background = element_rect(fill = "white", colour = NA)
  )

ggplot(df, aes(x = sex, y = median_itd,
               fill = sex, colour = sex)) +
  geom_violin(alpha = 0.8,
              trim = FALSE) +
  scale_fill_manual(values = sex_colours, labels = sex_labels,
                    guide = "none") +
  scale_colour_manual(values = sex_colours, labels = sex_labels,
                      guide = "none") +
  scale_x_discrete(labels = sex_labels) +
  labs(x = NULL, y = "Body size (ITD, mm)") +
  theme_classic(base_size = 12) +
  theme(
    plot.background = element_rect(fill = "white", colour = NA)
  )

# predicted % change

# Calculate trends 

year_start <- min(df$developmentYear)
year_end   <- max(df$developmentYear)
years_seq  <- seq(year_start, year_end, length.out = 100)

pred_from_emtrends <- sex_trend_df %>%
  select(sex, slope_per_year, se_per_year) %>%
  crossing(years = years_seq) %>%
  mutate(
    fit_pct   = (exp(slope_per_year * (years - year_start)) - 1) * 100,
    upper_pct = (exp((slope_per_year + 1.96 * se_per_year) *
                       (years - year_start)) - 1) * 100,
    lower_pct = (exp((slope_per_year - 1.96 * se_per_year) *
                       (years - year_start)) - 1) * 100
  )

ggplot(pred_from_emtrends,
       aes(x = years, y = fit_pct,
           colour = sex, fill = sex)) +
  geom_hline(yintercept = 0, linetype = "dashed",
             colour = "grey40", linewidth = 0.5) +
  geom_ribbon(aes(ymin = lower_pct, ymax = upper_pct),
              alpha = 0.2, colour = NA) +
  geom_line(linewidth = 1) +
  #annotate("text", x = 2024, y = -4.3,
           #label = "-4.95% (P", hjust = 1, size = 3.5,
           #fontface = "bold",
           #colour = sex_colours["male"]) +
  #annotate("text", x = 2024, y = 2.1,
           #label = "+1.44%", hjust = 1, size = 3.5,
           #fontface = "bold",
           #colour = sex_colours["female"]) +
  scale_colour_manual(values = sex_colours, labels = sex_labels,
                      name = "Sex") +
  scale_fill_manual(values = sex_colours, labels = sex_labels,
                    name = "Sex") +
  scale_x_continuous(breaks = seq(1900, 2025, by = 25)) +
  scale_y_continuous(labels = function(x) paste0(x, "%"),
                     limits = c(-12, 8)) +
  labs(x = "Development year",
       y = "Predicted % change in\nbody size relative to 1900") +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "right",
    plot.background = element_rect(fill = "white", colour = NA)
  )


#### Figure 3 ####


# Marginal effects for all four predictors
get_marginal_df <- function(model, predictor_name,
                             seq_range = c(-2, 2)) {
  pred_seq <- seq(seq_range[1], seq_range[2], length.out = 100)
  b_year   <- coef(model)["developmentYear_scaled"]
  b_int    <- coef(model)[paste0("developmentYear_scaled:",
                                  predictor_name)]
  vcov_m   <- vcov(model)
  se_year  <- sqrt(vcov_m["developmentYear_scaled",
                           "developmentYear_scaled"])
  se_int   <- sqrt(vcov_m[paste0("developmentYear_scaled:",
                                  predictor_name),
                           paste0("developmentYear_scaled:",
                                  predictor_name)])
  cov_yi   <- vcov_m["developmentYear_scaled",
                      paste0("developmentYear_scaled:",
                             predictor_name)]
  se_slope    <- sqrt(se_year^2 + pred_seq^2 * se_int^2 +
                        2 * pred_seq * cov_yi)
  slope       <- (b_year + b_int * pred_seq) / scaling_params$year_sd
  slope_upper <- (b_year + b_int * pred_seq + 1.96 * se_slope) /
    scaling_params$year_sd
  slope_lower <- (b_year + b_int * pred_seq - 1.96 * se_slope) /
    scaling_params$year_sd
  tibble(
    predictor = predictor_name,
    x         = pred_seq,
    slope     = (exp(slope * 10) - 1) * 100,
    upper     = (exp(slope_upper * 10) - 1) * 100,
    lower     = (exp(slope_lower * 10) - 1) * 100
  )
}

all_marginal <- bind_rows(
  get_marginal_df(m_mechanism, "ppt_within_scaled"),
  get_marginal_df(m_mechanism, "tmean_within_scaled"),
  get_marginal_df(m_mechanism, "cropland_within_scaled"),
  get_marginal_df(m_mechanism, "urban_within_scaled")
) %>%
  mutate(
    predictor = factor(predictor,
                       levels = c("ppt_within_scaled",
                                  "tmean_within_scaled",
                                  "cropland_within_scaled",
                                  "urban_within_scaled"),
                       labels = c("Annual precip.",
                                  "Mean annual temp.",
                                  "Cropland intensity",
                                  "Urban. intensity")),
    sig = predictor == "Annual precip."
  )

p_labels <- tibble(
  predictor = factor(
    c("Annual precip.", "Mean annual temp.",
      "Cropland intensity", "Urban. intensity"),
    levels = c("Annual precip.", "Mean annual temp.",
               "Cropland intensity", "Urban. intensity")
  ),
  label    = c("***P* = 0.010**", "*P* = 0.135",
               "*P* = 0.515",     "*P* = 0.115"),
  colour   = c("#2A9D8F", "grey30", "grey30", "grey30"),
  x        = -1.5,
  slope    =  0.9
)



# environmental trends
env_ca <- read.csv("california_environmental_trends.csv")

env_marginal <- env_ca %>%
  mutate(
    variable = case_when(
      variable == "Temperature (°C)"       ~ "Mean annual temp.",
      variable == "Precipitation (mm)"     ~ "Annual precip.",
      variable == "Urbanisation intensity" ~ "Urban. intensity",
      variable == "Cropland intensity"     ~ "Cropland intensity",
      TRUE ~ as.character(variable)
    ),
    variable = factor(variable,
                      levels = c("Annual precip.",
                                 "Mean annual temp.",
                                 "Cropland intensity",
                                 "Urban. intensity"))
  )

# Decadal means for cropland only
cropland_decadal <- env_ca %>%
  filter(variable == "Cropland intensity") %>%
  mutate(decade = floor(year / 10) * 10) %>%
  group_by(decade) %>%
  summarise(value = mean(value, na.rm = TRUE),
            .groups = "drop") %>%
  rename(year = decade) %>%
  mutate(variable = factor("Cropland intensity",
                           levels = levels(env_marginal$variable)))

# Combined — annual for all except cropland
env_plot <- env_marginal %>%
  filter(variable != "Cropland intensity") %>%
  bind_rows(cropland_decadal)


#precip variability
ppt_annual <- ppt_ca %>%
  select(year, value) %>%
  mutate(abs_dev = abs(value - mean(value, na.rm = TRUE)))

ggplot(ppt_annual, aes(x = year, y = abs_dev)) +
  geom_point(size = 1.2, colour = "grey50", alpha = 0.6) +
  geom_smooth(method = "lm", se = TRUE,
              colour = "#2A9D8F", fill = "#2A9D8F",
              linewidth = 1, alpha = 0.2)  +
  labs(x = "Year",
       y = "Absolute dev. from\nmean precip. (mm)") +
  theme_classic(base_size = 12) +
  ylim(0,47)+
  scale_x_continuous(breaks = c(1900, 1925, 1950, 1975, 2000, 2025))+
  theme(
    plot.background = element_rect(fill = "white", colour = NA)
  )

#pixel wise precip variability plot
ppt_var_df <- as.data.frame(ppt_var_trend, xy = TRUE) %>%
  rename(slope = 3) %>%
  filter(!is.na(slope))

# California outline
california_sf <- states(cb = TRUE, resolution = "500k") %>%
  filter(NAME == "California") %>%
  st_transform(4326)

# Plot
ggplot() +
  geom_sf(data = mexico,
          fill = "grey88", colour = "grey60", linewidth = 0.3) +
  geom_sf(data = western_states,
          fill = "grey88", colour = "grey60", linewidth = 0.3) +
  geom_raster(data = ppt_var_df,
              aes(x = x, y = y, fill = slope)) +
  geom_sf(data = california,
          fill = NA, colour = "grey30", linewidth = 0.4) +
   scale_fill_gradient2(
     low      = "#2A9D8F",
     mid      = "#FCCA46",
     high     = "#E76F51",
     midpoint = 0,
     name = "Change in\nvariability\n(mm/year)",
     oob      = scales::squish,
     limits   = c(-0.05, 0.05)
   ) +
  coord_sf(xlim = c(-125.5, -113),
           ylim = c(31.5, 42.5),
           expand = FALSE, clip = "on") +
  theme_void(base_size = 12) +
  theme(
    plot.title       = element_text(face = "bold", size = 11,
                                    hjust = 0.5),
    legend.position  = "right",
    panel.background = element_rect(fill = "#d6eaf8", colour = NA),
    panel.border     = element_rect(colour = "grey40", fill = NA,
                                    linewidth = 0.5),
    plot.background  = element_rect(fill = "white", colour = NA)
  )


## 4x2 environmental trend plots
# Helper: trend plot for one variable
make_trend_plot <- function(data, var_label, sig, show_x = FALSE,
                            y_limits = NULL, y_expand = 0.6) {
  p <- ggplot(data, aes(x = year, y = value)) +
    geom_point(size = 0.8, alpha = 0.4, colour = "grey50") +
    geom_smooth(method = "loess", se = TRUE, span = 0.75,
                colour = ifelse(sig, "#2A9D8F", "grey60"),
                fill   = ifelse(sig, "#2A9D8F", "grey80"),
                linewidth = 1) +
    labs(title = var_label, x = if (show_x) "Year" else NULL, y = NULL) +
    theme_classic(base_size = 11) +
    theme(
      plot.title   = element_text(face = "bold", size = 10, hjust = 0.5),
      axis.text.x  = if (show_x) element_text(angle = 45, hjust = 1) else element_blank(),
      axis.ticks.x = if (show_x) element_line() else element_blank(),
      plot.margin  = margin(1, 2, 1, 2)
    )
  
  if (!is.null(y_limits)) {
    p <- p + scale_y_continuous(limits = y_limits,
                                expand = expansion(mult = 0.1))
  } else {
    p <- p + scale_y_continuous(expand = expansion(mult = y_expand))
  }
  p
}

# Helper: mechanism effect plot for one variable
make_mech_plot <- function(data, var_label, sig, show_x = FALSE,
                           y_limits = NULL) {
  p <- ggplot(data, aes(x = x, y = slope)) +
    geom_hline(yintercept = 0, linetype = "dashed",
               colour = "grey40", linewidth = 0.4) +
    geom_ribbon(aes(ymin = lower, ymax = upper),
                fill = ifelse(sig, "#2A9D8F", "grey80"), alpha = 0.3) +
    geom_line(colour = ifelse(sig, "#2A9D8F", "grey60"), linewidth = 1) +
    scale_x_continuous(breaks = c(-1, 0, 1),
                       limits = c(-1.6, 1.6),
                       labels = c("Low", "Mean", "High")) +
    labs(title = var_label,
         x = if (show_x) "Deviation from cell mean" else NULL, y = NULL) +
    theme_classic(base_size = 11) +
    theme(
      plot.title   = element_text(face = "bold", size = 10, hjust = 0.5),
      axis.text.x  = if (show_x) element_text(angle = 45, hjust = 1) else element_blank(),
      axis.ticks.x = if (show_x) element_line() else element_blank(),
      plot.margin  = margin(1, 2, 1, 2)
    )
  
  if (!is.null(y_limits)) {
    p <- p + scale_y_continuous(limits = y_limits)
  }
  p
}

# --- Build the 8 panels ---

p_trend_precip <- make_trend_plot(env_plot %>% filter(variable == "Annual precip."),
                                  "Annual precip.", sig = FALSE)
p_mech_precip  <- make_mech_plot(all_marginal %>% filter(predictor == "Annual precip."),
                                 "Annual precip.", sig = TRUE, y_limits = c(-1.3, 1.0))

p_trend_temp <- make_trend_plot(env_plot %>% filter(variable == "Mean annual temp."),
                                "Mean annual temp.", sig = TRUE)
p_mech_temp  <- make_mech_plot(all_marginal %>% filter(predictor == "Mean annual temp."),
                               "Mean annual temp.", sig = FALSE, y_limits = c(-1.3, 1.0))

p_trend_crop <- make_trend_plot(env_plot %>% filter(variable == "Cropland intensity"),
                                "Cropland intensity", sig = TRUE,
                                y_limits = c(5, 12))
p_mech_crop  <- make_mech_plot(all_marginal %>% filter(predictor == "Cropland intensity"),
                               "Cropland intensity", sig = FALSE, y_limits = c(-1.3, 1.0))

p_trend_urban <- make_trend_plot(env_plot %>% filter(variable == "Urban. intensity"),
                                 "Urban. intensity", sig = TRUE, show_x = TRUE)
p_mech_urban  <- make_mech_plot(all_marginal %>% filter(predictor == "Urban. intensity"),
                                "Urban. intensity", sig = FALSE, show_x = TRUE, y_limits = c(-1.3, 1.0))

# Assemble into 4x2 grid 
right_col  <- p_trend_precip / p_trend_temp / p_trend_crop / p_trend_urban
left_col <- p_mech_precip  / p_mech_temp  / p_mech_crop  / p_mech_urban

p_grid <- left_col | right_col

p_grid + plot_annotation(
  title = "Environmental trends and effects on male body size",
  theme = theme(plot.title = element_text(face = "bold", size = 13))
)


#### Figure 4: Spatial drivers of body size ####

# Sex effects dot plot

p_sex_effects <- sex_effects %>%
  mutate(
    predictor = factor(predictor,
                       levels = c("Temperature", "Precipitation",
                                  "Urbanisation", "Cropland")),
    sex = factor(sex, levels = c("female", "male"))
  ) %>%
  ggplot(aes(x = percent_per_sd, y = predictor,
             colour = sex)) +
  geom_vline(xintercept = 0, linetype = "dashed",
             colour = "grey40", linewidth = 0.5) +
  geom_errorbarh(aes(xmin = lower_percent, xmax = upper_percent),
                 height = 0.2, linewidth = 0.8,
                 position = position_dodge(width = 0.5)) +
  geom_point(size = 3,
             position = position_dodge(width = 0.5)) +
  scale_colour_manual(values = sex_colours, labels = sex_labels,
                      name = "Sex") +
  scale_x_continuous(labels = function(x) paste0(x, "%")) +
  labs(x = "% change in body size\nper SD of predictor",
       y = NULL) +
  theme_classic(base_size = 12) +
  theme(
    legend.position    = "right",
    panel.grid.major.y = element_blank(),
    axis.line.x        = element_line(colour = "grey40"),
    plot.background    = element_rect(fill = "white", colour = NA),
    axis.text.y = element_text(angle=45)
  )


# Cropland and urbanisation maps 
hyde_cropland_files <- list.files(
  "HYDE_cropland_data",
  pattern = "cropland\\d{4}AD\\.asc$",
  recursive = TRUE, full.names = TRUE)

hisdac_files <- list.files(
  "HISDAC_urbanization_data",
  pattern = "\\.tif$", full.names = TRUE)

california_vect <- states(cb = TRUE, resolution = "500k") %>%
  filter(NAME == "California") %>%
  vect()

# Cropland: most recent year (2023) 

most_recent_hyde <- hyde_cropland_files[grepl("2023",
                                              hyde_cropland_files)]
cat("Cropland file:", basename(most_recent_hyde), "\n")

ca_hyde         <- project(california_vect,
                           crs(rast(most_recent_hyde)))
cropland_recent <- rast(most_recent_hyde) %>%
  crop(ca_hyde) %>%
  mask(ca_hyde) %>%
  project("EPSG:4326")

cropland_recent_df <- as.data.frame(cropland_recent, xy = TRUE) %>%
  rename(value = 3) %>%
  filter(!is.na(value))

# Urbanisation: most recent year (2020)

most_recent_hisdac <- hisdac_files[grepl("2020", hisdac_files)]
cat("Urbanisation file:", basename(most_recent_hisdac), "\n")

california_hisdac <- project(california_vect,
                             crs(rast(most_recent_hisdac)))
urban_recent <- rast(most_recent_hisdac) %>%
  crop(california_hisdac) %>%
  mask(california_hisdac) %>%
  project("EPSG:4326")

urban_recent_df <- as.data.frame(urban_recent, xy = TRUE) %>%
  rename(value = 3) %>%
  filter(!is.na(value))

# Cropland map

cropland_max <- max(cropland_recent_df$value, na.rm = TRUE)

ggplot() +
  geom_sf(data = mexico,
          fill = "grey88", colour = "grey60", linewidth = 0.3) +
  geom_sf(data = western_states,
          fill = "grey88", colour = "grey60", linewidth = 0.3) +
  geom_raster(data = cropland_recent_df,
              aes(x = x, y = y, fill = value)) +
  geom_sf(data = california,
          fill = NA, colour = "grey30", linewidth = 0.4) +
  scale_fill_gradient2(
    low      = "#335F70",
    mid      = "#E76F51",
    high     = "#FCCA46",
    midpoint = cropland_max / 2,
    name     = "Cropland\nintensity",
    limits   = c(0, cropland_max),
    na.value = "white"
  ) +
  coord_sf(xlim = c(-125.5, -113), ylim = c(31.5, 42.5),
           expand = FALSE, clip = "on") +
  ggtitle("Cropland intensity (2023)") +
  theme_void(base_size = 12) +
  theme(
    plot.title       = element_text(face = "bold", size = 11,
                                    hjust = 0.5),
    legend.position  = "right",
    panel.background = element_rect(fill = "#d6eaf8", colour = NA),
    panel.border     = element_rect(colour = "grey40", fill = NA,
                                    linewidth = 0.5),
    plot.background  = element_rect(fill = "white", colour = NA)
  )


#  Urbanisation map 
urban_recent_df <- urban_recent_df %>%
  mutate(value_log = log1p(value))  # log(x + 1) to handle zeros

# Get the actual max value
urban_max <- quantile(urban_recent_df$value_log, 0.99, na.rm = TRUE)
cat("Truncating at:", urban_max, "\n")

ggplot() +
  geom_sf(data = mexico,
          fill = "grey88", colour = "grey60", linewidth = 0.3) +
  geom_sf(data = western_states,
          fill = "grey88", colour = "grey60", linewidth = 0.3) +
  geom_raster(data = urban_recent_df,
              aes(x = x, y = y, fill = value_log)) +
  geom_sf(data = california,
          fill = NA, colour = "grey30", linewidth = 0.4) +
  scale_fill_gradient2(
    low      = "#335F70",
    mid      = "#E76F51",
    high     = "#FCCA46",
    midpoint = urban_max / 2,
    name     = "Urbanisation\nintensity\n(log scale)",
    limits   = c(0, urban_max),
    oob      = scales::squish
  ) +
  coord_sf(xlim = c(-125.5, -113), ylim = c(31.5, 42.5),
           expand = FALSE, clip = "on") +
  ggtitle("Urbanisation intensity (2020)") +
  theme_void(base_size = 12) +
  theme(
    plot.title       = element_text(face = "bold", size = 11,
                                    hjust = 0.5),
    legend.position  = "right",
    panel.background = element_rect(fill = "#d6eaf8", colour = NA),
    panel.border     = element_rect(colour = "grey40", fill = NA,
                                    linewidth = 0.5),
    plot.background  = element_rect(fill = "white", colour = NA)
  )



#### Supplementary figures ####


#Genus effects heatmap

heatmap_df <- genus_effects %>%
  mutate(
    sig = case_when(
      genus == "Bombus"     & predictor == "Temperature"   ~ "*",
      genus == "Hesperapis" & predictor == "Precipitation" ~ "*",
      genus == "Bombus"     & predictor == "Urbanisation"  ~ "*",
      genus == "Hesperapis" & predictor == "Cropland"      ~ "*",
      TRUE ~ ""
    ),
    predictor = factor(predictor,
                       levels = c("Temperature", "Precipitation",
                                  "Urbanisation", "Cropland")),
    genus = factor(genus,
                   levels = c("Xylocopa", "Bombus", "Andrena",
                              "Anthophora", "Osmia",
                              "Lasioglossum", "Hesperapis"))
  )

ggplot(heatmap_df, aes(x = predictor, y = genus,
                       fill = percent_per_sd)) +
  geom_tile(colour = "white", linewidth = 0.5) +
  geom_text(aes(label = sig), size = 5, vjust = 0.75,
            colour = "black") +
  scale_fill_gradient2(
    low = "#3A9AB2", mid = "white", high = "#EF5703",
    midpoint = 0, name = "% change\nper SD",
    limits = c(-3.5, 3.5), oob = scales::squish
  ) +
  labs(x = NULL, y = NULL) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x     = element_text(size = 11, angle = 45,
                                   hjust = 1, colour = "grey20"),
    axis.text.y     = element_text(size = 11, face = "italic"),
    panel.grid      = element_blank(),
    legend.position = "right",
    plot.background = element_rect(fill = "white", colour = NA)
  )

#Genus-level temporal trends 

trend_df_sex %>%
  ggplot(aes(x = reorder(genus, total_percent_change),
             y = total_percent_change, colour = sex)) +
  geom_hline(yintercept = 0, linetype = "dashed",
             colour = "grey40", linewidth = 0.5) +
  geom_errorbar(aes(ymin = total_lower, ymax = total_upper),
                width = 0.3, linewidth = 0.9,
                position = position_dodge(width = 0.5)) +
  geom_point(size = 3, position = position_dodge(width = 0.5)) +
  scale_colour_manual(values = sex_colours, labels = sex_labels,
                      name = "Sex") +
  scale_y_continuous(labels = function(x) paste0(x, "%")) +
  labs(x = NULL,
       y = paste0("Change in body size (",
                  min(df$developmentYear), "–",
                  max(df$developmentYear), ")")) +
  coord_flip() +
  theme_classic(base_size = 12) +
  theme(
    strip.background   = element_blank(),
    legend.position    = "right",
    panel.grid.minor   = element_blank(),
    panel.grid.major.y = element_blank(),
    axis.text.y        = element_text(face = "italic", size = 11),
    axis.line.x        = element_line(colour = "grey40"),
    plot.background    = element_rect(fill = "white", colour = NA)
  )
