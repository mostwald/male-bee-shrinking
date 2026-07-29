# Bee body size change in California, 1900-2024
# Data extraction script
# May 2026

# This script extracts environmental data from raw gridded datasets and saves pre-processed outputs used in the main analysis. 
# You don't need to run this script to reproduce the main analysis and figures. Pre-processed data are available with source repository.

# To run this script from scratch you will need to download:
#   PRISM climate data (4km gridded, 1900-2024):
#   https://prism.oregonstate.edu/
#   Variables: tmean, tmax, ppt (monthly)

#   HISDAC-US urbanisation data:
#   https://dataverse.harvard.edu/dataverse/hisdac-us
#   Uhl et al. (2021) doi:10.7910/DVN/TMPGN0

#   HYDE 3.3 cropland data:
#   https://www.pbl.nl/en/image/links/hyde
#   Klein Goldewijk et al. (2017) doi:10.5194/essd-9-927-2017

#### Setup ####

library(tidyverse)
library(terra)
library(sf)
library(tigris)
library(prism)

#### File paths ####
ppt_files <- list.files(
  "PRISM_ppt", # update path here
  pattern = ".*ppt.*\\.bil$",
  recursive = TRUE, full.names = TRUE)

tmean_files <- list.files(
  "PRISM_tmean",
  pattern = ".*tmean.*\\.bil$",
  recursive = TRUE, full.names = TRUE)

tmax_files <- list.files(
  "PRISM_tmax",
  pattern = ".*tmax.*\\.bil$",
  recursive = TRUE, full.names = TRUE)

hisdac_files <- list.files(
  "HISDAC_urbanization_data",
  pattern = "\\.tif$", full.names = TRUE)

hyde_cropland_files <- list.files(
  "HYDE_cropland_data",
  pattern = "cropland\\d{4}AD\\.asc$",
  recursive = TRUE, full.names = TRUE)

#### California boundary ####

california_vect <- states(cb = TRUE, resolution = "500k") %>%
  filter(NAME == "California") %>%
  vect()


#### Extract statewide annual means ####

extract_annual_mean <- function(files, years = 1900:2024) {
  map_dfr(years, function(yr) {
    yr_files <- files[grepl(paste0(yr, "\\d{2}"), files) &
                        grepl("\\.bil$", files)]
    if (length(yr_files) == 0) {
      message("No files found for year: ", yr)
      return(NULL)
    }
    if (length(yr_files) < 12) {
      message("Only ", length(yr_files), " months for year: ", yr)
    }
    r      <- rast(yr_files)
    r_mean <- mean(r)
    r_ca   <- mask(crop(r_mean, california_vect), california_vect)
    tibble(year = yr,
           value = global(r_ca, "mean", na.rm = TRUE)$mean)
  })
}

tmean_ca <- extract_annual_mean(tmean_files) %>%
  mutate(variable = "Temperature (°C)")

ppt_ca <- extract_annual_mean(ppt_files) %>%
  mutate(variable = "Precipitation (mm)")

urban_ca <- map_dfr(hisdac_files, function(f) {
  yr <- as.integer(gsub("(\\d{4})_BUPR\\.tif", "\\1", basename(f)))
  if (yr < 1900 | yr > 2024) return(NULL)
  r    <- rast(f)
  r_ca <- mask(crop(r, project(california_vect, crs(r))),
               project(california_vect, crs(r)))
  tibble(year  = yr,
         value = global(r_ca, "mean", na.rm = TRUE)$mean)
}) %>% mutate(variable = "Urbanisation intensity")

cropland_ca <- map_dfr(hyde_cropland_files, function(f) {
  yr <- as.integer(gsub(".*cropland(\\d{4})AD\\.asc", "\\1",
                         basename(f)))
  if (yr < 1900 | yr > 2023) return(NULL)
  r    <- rast(f)
  r_ca <- mask(crop(r, project(california_vect, crs(r))),
               project(california_vect, crs(r)))
  tibble(year  = yr,
         value = global(r_ca, "mean", na.rm = TRUE)$mean)
}) %>% mutate(variable = "Cropland intensity")

env_ca <- bind_rows(tmean_ca, ppt_ca, urban_ca, cropland_ca) %>%
  mutate(variable = factor(variable,
                           levels = c("Temperature (°C)",
                                      "Precipitation (mm)",
                                      "Urbanisation intensity",
                                      "Cropland intensity")))

write.csv(env_ca, "california_environmental_trends.csv",
          row.names = FALSE)
