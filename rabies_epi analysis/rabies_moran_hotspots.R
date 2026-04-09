#!/usr/bin/env Rscript

# =========================================================
# rabies_moran_hotspots.R
# County-level spatial clustering of raccoon rabies cases
# in Alabama using Global and Local Moran's I
# Uses labels_coords.csv to define study counties
# =========================================================

suppressPackageStartupMessages({
  library(sf)
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(spdep)
  library(stringr)
  library(tibble)
  library(tigris)
})

# -----------------------------
# User-defined paths
# -----------------------------
root_dir <- "/scratch/Raccoon_andrea/gabe_trial"
work_dir <- file.path(root_dir, "scripts/PopulationGenetics/rabies_epi")

rabies_csv <- "/scratch/Raccoon_andrea/gabe_trial/scripts/PopulationGenetics/AL_rabies_2009_2025.csv"
labels_csv <- "/scratch/Raccoon_andrea/gabe_trial/scripts/PopulationGenetics/labels_coords.csv"

out_dir <- file.path(work_dir, "outputs/moran_hotspots")
plot_dir <- file.path(out_dir, "plots")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

# optional cache directory for tigris
tigris_cache_dir <- file.path(root_dir, "data/tigris_cache")
dir.create(tigris_cache_dir, recursive = TRUE, showWarnings = FALSE)

options(tigris_use_cache = TRUE)
options(tigris_class = "sf")
Sys.setenv(TIGRIS_CACHE_DIR = tigris_cache_dir)

# -----------------------------
# Helper functions
# -----------------------------
clean_county <- function(x) {
  x %>%
    toupper() %>%
    str_replace_all("\\.", "") %>%
    str_replace_all("-", "") %>%
    str_replace_all("\\s+", "") %>%
    str_replace_all("COUNTY", "") %>%
    trimws()
}

make_lisa_class <- function(x, lag_x, pval, alpha = 0.05) {
  out <- rep("Not significant", length(x))
  sig <- pval <= alpha

  out[sig & x > 0 & lag_x > 0] <- "High-High"
  out[sig & x < 0 & lag_x < 0] <- "Low-Low"
  out[sig & x > 0 & lag_x < 0] <- "High-Low"
  out[sig & x < 0 & lag_x > 0] <- "Low-High"

  factor(
    out,
    levels = c("High-High", "Low-Low", "High-Low", "Low-High", "Not significant")
  )
}

# -----------------------------
# Read Alabama counties from tigris
# -----------------------------
cat("Downloading/reading Alabama counties from TIGER...\n")
counties_sf <- tigris::counties(state = "AL", cb = TRUE, year = 2020) %>%
  mutate(
    county = NAME,
    county_std = clean_county(NAME)
  )

# -----------------------------
# Read rabies totals
# -----------------------------
cat("Reading rabies totals...\n")
rabies <- read_csv(rabies_csv, show_col_types = FALSE)

if (!all(c("county", "total") %in% names(rabies))) {
  stop("rabies_counts_2009_2025.csv must contain columns named 'county' and 'total'.")
}

rabies <- rabies %>%
  mutate(
    county_std = clean_county(county),
    total = as.numeric(total)
  )

# -----------------------------
# Read labels_coords and define study counties
# -----------------------------
cat("Reading labels_coords.csv...\n")
labels <- read_csv(labels_csv, show_col_types = FALSE)

if (!"county" %in% names(labels)) {
  stop("labels_coords.csv must contain a column named 'county'.")
}

labels <- labels %>%
  mutate(county_std = clean_county(county)) %>%
  distinct(county_std)

# -----------------------------
# Join study counties with rabies data
# -----------------------------
counties_use <- counties_sf %>%
  filter(county_std %in% labels$county_std) %>%
  left_join(rabies %>% select(county_std, total), by = "county_std") %>%
  mutate(total = ifelse(is.na(total), 0, total))

if (nrow(counties_use) < 5) {
  stop("Too few counties after joining labels_coords and rabies data.")
}

cat("Number of study counties:", nrow(counties_use), "\n")

# -----------------------------
# Project geometry and create county points
# -----------------------------
counties_use <- st_make_valid(counties_use)
counties_proj <- st_transform(counties_use, 5070)  # NAD83 / Conus Albers

county_pts <- st_point_on_surface(counties_proj)
coords <- st_coordinates(county_pts)

# -----------------------------
# Build k-nearest neighbor spatial weights
# -----------------------------
k_val <- 5
if (nrow(counties_proj) <= k_val) {
  k_val <- max(1, nrow(counties_proj) - 1)
}

cat("Building k-nearest-neighbor network with k =", k_val, "\n")

knn_obj <- knearneigh(coords, k = k_val)
nb_obj <- knn2nb(knn_obj)
lw <- nb2listw(nb_obj, style = "W", zero.policy = TRUE)

# -----------------------------
# Variables for analysis
# -----------------------------
x_raw <- counties_proj$total
x_log <- log1p(x_raw)

# -----------------------------
# Global Moran's I
# -----------------------------
cat("Running Global Moran's I...\n")

global_raw <- moran.mc(x_raw, lw, nsim = 1000, zero.policy = TRUE)
global_log <- moran.mc(x_log, lw, nsim = 1000, zero.policy = TRUE)

global_results <- tibble(
  variable = c("raw_total", "log1p_total"),
  moran_I = c(unname(global_raw$statistic), unname(global_log$statistic)),
  p_value = c(global_raw$p.value, global_log$p.value),
  permutations = 1000,
  n_counties = nrow(counties_proj),
  k_neighbors = k_val
)

write_csv(global_results, file.path(out_dir, "global_moransI_results.csv"))

# -----------------------------
# Local Moran's I helper
# -----------------------------
run_local_moran <- function(x, lw, counties_sf, var_name) {
  x_centered <- as.numeric(scale(x, center = TRUE, scale = FALSE))
  lag_x <- lag.listw(lw, x_centered, zero.policy = TRUE)

  local_res <- tryCatch(
    {
      localmoran_perm(x_centered, lw, nsim = 1000, zero.policy = TRUE)
    },
    error = function(e) {
      message("localmoran_perm unavailable; using localmoran().")
      localmoran(x_centered, lw, zero.policy = TRUE)
    }
  )

  local_df <- as.data.frame(local_res)
  names(local_df) <- make.names(names(local_df))

  p_col <- names(local_df)[grepl("^Pr", names(local_df))]
  if (length(p_col) == 0) {
    stop("Could not identify Local Moran p-value column.")
  }

  tibble(
    county = counties_sf$county,
    county_std = counties_sf$county_std,
    total = counties_sf$total,
    value_used = x,
    centered_value = x_centered,
    lag_centered_value = lag_x,
    local_I = local_df$Ii,
    z_score = if ("Z.Ii" %in% names(local_df)) local_df$Z.Ii else NA_real_,
    p_value = local_df[[p_col[1]]],
    lisa_class = make_lisa_class(x_centered, lag_x, local_df[[p_col[1]]]),
    variable = var_name
  )
}

# -----------------------------
# Local Moran's I
# -----------------------------
cat("Running Local Moran's I...\n")

local_raw_df <- run_local_moran(x_raw, lw, counties_proj, "raw_total")
local_log_df <- run_local_moran(x_log, lw, counties_proj, "log1p_total")

local_results <- bind_rows(local_raw_df, local_log_df)
write_csv(local_results, file.path(out_dir, "local_moransI_results.csv"))

# join results back to sf for mapping
map_raw <- counties_proj %>%
  left_join(
    local_raw_df %>% select(county_std, local_I, p_value, lisa_class),
    by = "county_std"
  )

map_log <- counties_proj %>%
  left_join(
    local_log_df %>% select(county_std, local_I, p_value, lisa_class),
    by = "county_std"
  )

  # full Alabama counties for background
counties_all_proj <- st_transform(counties_sf, 5070)

# state outline
state_outline <- st_union(counties_all_proj) |> st_sf()

# -----------------------------
# Plot 1: rabies case totals
# -----------------------------
p_cases <- ggplot() +
  # background: all Alabama counties
  geom_sf(data = counties_all_proj, fill = "#f7f4f2", color = "white", linewidth = 0.25) +
  # study counties with rabies totals
  geom_sf(data = counties_proj, aes(fill = total), color = "white", linewidth = 0.35) +
  # state outline
  geom_sf(data = state_outline, fill = NA, color = "black", linewidth = 0.8) +
  scale_fill_gradient(
    low = "#fee8c8",
    high = "#8c2d04",
    name = "Total rabies cases\n(2009–2025)"
  ) +
  theme_minimal(base_size = 12) +
  labs(
    title = "Reported raccoon rabies cases across Alabama study counties",
    subtitle = "County-level cumulative totals, 2009–2025"
  ) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

ggsave(
  file.path(plot_dir, "rabies_case_totals_map.png"),
  p_cases,
  width = 7,
  height = 7,
  dpi = 300
)

# -----------------------------
# Plot 2: Local Moran hotspot map (raw)
# -----------------------------
p_lisa_raw <- ggplot() +
  # background: all Alabama counties
  geom_sf(data = counties_all_proj, fill = "#f7f4f2", color = "white", linewidth = 0.25) +
  # study counties with LISA classes
  geom_sf(data = map_raw, aes(fill = lisa_class), color = "white", linewidth = 0.35) +
  # state outline
  geom_sf(data = state_outline, fill = NA, color = "black", linewidth = 0.8) +
  scale_fill_manual(
    values = c(
      "High-High" = "#b2182b",
      "Low-Low" = "#2166ac",
      "High-Low" = "#ef8a62",
      "Low-High" = "#67a9cf",
      "Not significant" = "grey85"
    ),
    drop = FALSE,
    name = "Local Moran's I\n(raw totals)"
  ) +
  theme_minimal(base_size = 12) +
  labs(
    title = "Local spatial clustering of reported raccoon rabies cases",
    subtitle = "Hotspot map based on county-level raw totals"
  ) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

ggsave(
  file.path(plot_dir, "local_moransI_hotspots_raw.png"),
  p_lisa_raw,
  width = 7,
  height = 7,
  dpi = 300
)

# -----------------------------
# Plot 3: Local Moran hotspot map (log)
# -----------------------------
p_lisa_log <- ggplot() +
  # background: all Alabama counties
  geom_sf(data = counties_all_proj, fill = "grey92", color = "white", linewidth = 0.25) +
  # study counties with LISA classes
  geom_sf(data = map_log, aes(fill = lisa_class), color = "white", linewidth = 0.35) +
  # state outline
  geom_sf(data = state_outline, fill = NA, color = "black", linewidth = 0.8) +
  scale_fill_manual(
    values = c(
      "High-High" = "#b2182b",
      "Low-Low" = "#2166ac",
      "High-Low" = "#ef8a62",
      "Low-High" = "#67a9cf",
      "Not significant" = "grey85"
    ),
    drop = FALSE,
    name = "Local Moran's I\n(log totals)"
  ) +
  theme_minimal(base_size = 12) +
  labs(
    title = "Local spatial clustering of reported raccoon rabies cases",
    subtitle = "Hotspot map based on log-transformed county totals"
  ) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

ggsave(
  file.path(plot_dir, "local_moransI_hotspots_log.png"),
  p_lisa_log,
  width = 7,
  height = 7,
  dpi = 300
)

# -----------------------------
# Save geospatial outputs
# -----------------------------
st_write(
  map_raw,
  dsn = file.path(out_dir, "rabies_local_moransI_raw.gpkg"),
  delete_dsn = TRUE,
  quiet = TRUE
)

st_write(
  map_log,
  dsn = file.path(out_dir, "rabies_local_moransI_log.gpkg"),
  delete_dsn = TRUE,
  quiet = TRUE
)

# -----------------------------
# Console summary
# -----------------------------
cat("\nDone.\n")
cat("\nGlobal Moran's I results:\n")
print(global_results)
print(local_results)

cat("\nOutputs written to:\n")
cat(out_dir, "\n")