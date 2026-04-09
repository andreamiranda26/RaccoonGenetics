#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(sf)
  library(tigris)
  library(ggspatial)
  library(patchwork)
  library(rnaturalearth)
  library(rnaturalearthdata)
})

options(tigris_use_cache = TRUE)

ROOT_DIR <- "/scratch/Raccoon_andrea/gabe_trial"
IN_CSV   <- file.path(ROOT_DIR, "results/rabies_epi/county_rabies_orv_distance.csv")
SAMP_CSV <- file.path(ROOT_DIR, "results/PCA_Moran_county/county_sample_sizes.csv")
OUT_DIR  <- file.path(ROOT_DIR, "results/rabies_epi")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

OUT_PNG  <- file.path(OUT_DIR, "map_rabies_heat_orv_outline_samples_inset.png")

clean_county <- function(x) {
  x <- tolower(trimws(x))
  x <- gsub("[[:punct:]]", "", x)
  x <- gsub("\\s+", "", x)
  x
}

# ---- manually define ORV counties here ----
orv_counties <- c(
  "Blount", "Cherokee", "Cullman", "DeKalb",
  "Etowah", "Jackson", "Jefferson", "Marshall",
  "Morgan", "St. Clair", "Shelby",
  "Baldwin"
)
orv_keys <- clean_county(orv_counties)

if (!file.exists(IN_CSV)) {
  stop("Input CSV not found: ", IN_CSV,
       "\nRun your ORV/rabies script first so it writes county_rabies_orv_distance.csv")
}

if (!file.exists(SAMP_CSV)) {
  stop("Sample-size CSV not found: ", SAMP_CSV,
       "\nMake sure county_sample_sizes.csv exists from your PCA script.")
}

message("Reading rabies data: ", IN_CSV)
dat <- read_csv(IN_CSV, show_col_types = FALSE) |> as.data.frame()
names(dat) <- tolower(names(dat))

message("Reading sample sizes: ", SAMP_CSV)
samp <- read_csv(SAMP_CSV, show_col_types = FALSE) |> as.data.frame()
names(samp) <- tolower(names(samp))

# ---- flexible column detection for rabies file ----
if (!("county_key" %in% names(dat))) {
  county_like <- grep("^county($|\\.)|county_name|name$", names(dat), value = TRUE)
  if (length(county_like) == 0) {
    stop(
      "Could not find county or county_key columns in: ", IN_CSV,
      "\nColumns found: ", paste(names(dat), collapse = ", "),
      "\nFix: ensure the rabies script writes either 'county' or 'county_key'."
    )
  }
  dat$county_key <- clean_county(dat[[county_like[1]]])
}

if (!("total_cases" %in% names(dat))) {
  alt_cases <- c("cases", "total", "total_rabies", "rabies_total", "rabies_totals")
  hit <- intersect(alt_cases, names(dat))
  if (length(hit) > 0) {
    dat$total_cases <- as.numeric(dat[[hit[1]]])
  } else {
    stop(
      "Could not find 'total_cases' in: ", IN_CSV,
      "\nColumns found: ", paste(names(dat), collapse = ", "),
      "\nFix: make sure your rabies script outputs a numeric 'total_cases' column."
    )
  }
}

dat2 <- dat |>
  transmute(
    county_key = clean_county(county_key),
    total_cases = as.numeric(total_cases)
  )

# ---- clean sample-size file ----
if (!all(c("county", "n") %in% names(samp))) {
  stop(
    "Sample-size file must contain columns named 'county' and 'N'.",
    "\nColumns found: ", paste(names(samp), collapse = ", ")
  )
}

samp2 <- samp |>
  transmute(
    county_key = clean_county(county),
    sample_n = as.integer(n)
  )

message("Loading Alabama county polygons (tigris)...")
al <- tigris::counties(state = "AL", year = 2020, class = "sf") |>
  st_transform(4326) |>
  mutate(county_key = clean_county(NAME))

al2 <- al |>
  left_join(dat2, by = "county_key") |>
  left_join(samp2, by = "county_key") |>
  mutate(orv_any = ifelse(county_key %in% orv_keys, 1, 0))

# ---- label points ----
label_pts <- st_point_on_surface(al2) |>
  filter(!is.na(sample_n) & sample_n > 0)

# ---- quick check of outlined ORV counties ----
message("ORV counties being outlined:")
print(al2 |>
  filter(orv_any == 1) |>
  st_drop_geometry() |>
  select(NAME) |>
  arrange(NAME))

# ---- load river centerlines ----
message("Loading river centerlines...")
rivers <- rnaturalearth::ne_download(
  scale = 10,
  type = "rivers_lake_centerlines",
  category = "physical",
  returnclass = "sf"
) |>
  st_transform(4326)

# crop to Alabama extent plus small buffer
al_bbox <- st_bbox(al2)
bbox_poly <- st_as_sfc(al_bbox) |>
  st_buffer(0.5)

rivers_al <- st_intersection(rivers, bbox_poly)

# keep Alabama-Coosa river system
target_patterns <- c("alabama", "coosa", "tallapoosa")
river_keep <- grepl(
  paste(target_patterns, collapse = "|"),
  tolower(rivers_al$name)
)

river_lines <- rivers_al[river_keep, ]

message("River features retained:")
if (nrow(river_lines) > 0) {
  print(unique(river_lines$name))
} else {
  warning("No river features matched the requested names.")
}

# ---- main Alabama map ----
p_main <- ggplot() +
  geom_sf(data = al, fill = "gray95", color = "#8B8878", linewidth = 0.2) +
  geom_sf(
    data = al2 |> filter(!is.na(total_cases)),
    aes(fill = total_cases),
    color = "gray95",
    linewidth = 0.2
  ) +
  geom_sf(
    data = river_lines,
    color = "royalblue3",
    linewidth = 1.2,
    lineend = "round"
  ) +
  geom_sf(
    data = al2 |> filter(orv_any == 1),
    fill = NA,
    color = "black",
    linewidth = 0.9
  ) +
  geom_sf_text(
    data = label_pts,
    aes(label = sample_n),
    size = 2.8,
    color = "black",
    fontface = "bold"
  ) +
  scale_fill_gradient(
    low  = "#FF8C00",
    high = "#8B4500",
    name = "Total rabies cases\n(2009–2025)",
    trans = "sqrt"
  ) +
  annotation_scale(
    location = "bl",
    width_hint = 0.25,
    line_width = 0.6,
    text_cex = 0.8
  ) +
  annotation_north_arrow(
    location = "br",
    which_north = "true",
    pad_x = unit(0.2, "in"),
    pad_y = unit(0.2, "in"),
    style = north_arrow_fancy_orienteering
  ) +
  coord_sf() +
  theme_minimal(base_size = 12) +
  theme(
    panel.background = element_rect(fill = "white", color = NA),
    plot.background  = element_rect(fill = "white", color = NA),
    panel.grid       = element_blank(),
    axis.text        = element_blank(),
    axis.title       = element_blank(),
    axis.ticks       = element_blank(),
    legend.title     = element_text(size = 10),
    legend.text      = element_text(size = 9)
  )

# ---- USA inset map ----
message("Loading USA states for inset...")
usa <- tigris::states(cb = TRUE, year = 2020, class = "sf") |>
  st_transform(4326)

alabama_state <- usa |>
  filter(STUSPS == "AL")

p_inset <- ggplot() +
  geom_sf(data = usa, fill = "gray90", color = "white", linewidth = 0.15) +
  geom_sf(data = alabama_state, fill = "#CC5500", color = "black", linewidth = 0.3) +
  coord_sf(xlim = c(-125, -66), ylim = c(24, 50), expand = FALSE) +
  theme_void() +
  theme(
    panel.background = element_rect(fill = "white", color = "black", linewidth = 0.3),
    plot.background  = element_rect(fill = "white", color = NA)
  )

# ---- combine with inset ----
p_final <- p_main +
  inset_element(
    p_inset,
    left = 0.78,
    right = 1.00,
    bottom = 0.62,
    top = 0.92,
    align_to = "full"
  )

ggsave(OUT_PNG, p_final, width = 8.5, height = 6.5, dpi = 300)

message("Saved: ", OUT_PNG)
message("DONE.")