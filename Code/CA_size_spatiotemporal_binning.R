# Bee body size change in California, 1900-2024
# Spatiotemporal binning script (after Guralnick et al 2020)
# May 2026

# Assigns specimens to optimal 150x150km grid cells and filters by genus-level temporal coverage criteria.

# Grid cell parameters:
# Cell size: 150km
# Min specimens per genus per cell: 10
# Min decades per genus per cell: 3
# Random offsets tested: 10,000

#### Setup ####

library(sf)
library(tidyverse)
library(terra)
library(prism)


df <- read.csv("before_filtering_CA_size_data.csv")

# Convert to sf and project to California Albers
sf_df <- df %>%
  filter(!is.na(decimalLongitude) & !is.na(decimalLatitude)) %>%
  st_as_sf(coords = c("decimalLongitude", "decimalLatitude"),
           crs = 4326) %>%
  st_transform(3310)  # California Albers Equal Area (metres)

sf_df$decade <- floor(sf_df$developmentYear / 10) * 10

cat("Specimens with valid coordinates:", nrow(sf_df), "\n")
cat("Specimens dropped (missing coords):",
    nrow(df) - nrow(sf_df), "\n")

# Parameters 
cellsize <- 150000   # 150km in metres
n_iter   <- 10000    # random grid offsets to test
bbox     <- st_bbox(sf_df)

# Helper: generate random grid offset
random_grid <- function(sf_obj, cellsize, bbox) {
  dx     <- runif(1, 0, cellsize)
  dy     <- runif(1, 0, cellsize)
  offset <- c(bbox["xmin"] + dx, bbox["ymin"] + dy)
  g      <- st_make_grid(sf_obj,
                         cellsize = cellsize,
                         square   = TRUE,
                         offset   = offset)
  st_sf(grid_id = seq_along(g), geometry = g)
}

# Helper: genus-level filtering 
apply_genus_filter <- function(joined,
                               min_specimens = 10,
                               min_decades   = 3) {
  joined %>%
    st_drop_geometry() %>%
    group_by(grid_id, genus, decade) %>%
    summarise(n = n(), .groups = "drop") %>%
    group_by(grid_id, genus) %>%
    summarise(
      total     = sum(n),
      n_decades = n_distinct(decade),
      .groups   = "drop"
    ) %>%
    filter(total >= min_specimens & n_decades >= min_decades)
}

# Optimisation loop 

best_count   <- -Inf
best_pts     <- NULL
best_grid    <- NULL
best_summary <- NULL

set.seed(42)

for (i in seq_len(n_iter)) {
  
  grid_sf          <- random_grid(sf_df, cellsize, bbox)
  joined           <- st_join(sf_df, grid_sf,
                              join = st_within, left = FALSE)
  genus_cell_valid <- apply_genus_filter(joined)
  
  n_retained <- joined %>%
    st_drop_geometry() %>%
    semi_join(genus_cell_valid, by = c("grid_id", "genus")) %>%
    nrow()
  
  if (n_retained > best_count) {
    best_count   <- n_retained
    best_grid    <- grid_sf
    best_pts     <- joined
    best_summary <- genus_cell_valid
    message("Iteration ", i, ": new best — ", n_retained,
            " specimens across ",
            n_distinct(genus_cell_valid$grid_id), " cells.")
  }
  
  if (i %% 1000 == 0) {
    message("Completed ", i, " iterations. Best = ",
            best_count, " specimens retained.")
  }
}

message("Finished ", n_iter, " iterations.")
message("Best grid retained ", best_count, " specimens.")

# Apply genus-level filter to best grid 
joined_best <- st_join(sf_df, best_grid,
                       join = st_within, left = FALSE)

genus_cell_valid_final <- apply_genus_filter(joined_best)

df_grid_filtered <- joined_best %>%
  semi_join(genus_cell_valid_final, by = c("grid_id", "genus"))

#### Save outputs ####

# Grid geometry
st_write(best_grid, "best_grid_150km_may2026.gpkg",
         delete_layer = TRUE)

# Filtered points with geometry
st_write(df_grid_filtered, "best_points_150km_may2026.gpkg",
         delete_layer = TRUE)

# Transform to WGS84 and extract decimal degree coordinates
df_grid_filtered_wgs84 <- st_transform(df_grid_filtered, 4326) %>%
  mutate(
    decimalLongitude = st_coordinates(geometry)[, 1],
    decimalLatitude  = st_coordinates(geometry)[, 2]
  )

# Drop geometry and save final analysis dataset
df_final <- st_drop_geometry(df_grid_filtered_wgs84)

write.csv(df_final,
          "grid_filtered_CA_size_data.csv",
          row.names = FALSE)

