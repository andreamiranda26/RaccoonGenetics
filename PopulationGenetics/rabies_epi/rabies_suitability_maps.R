#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(readr)
})

ROOT_DIR <- "/scratch/Raccoon_andrea/gabe_trial"

LABELS_FILE <- file.path(ROOT_DIR, "scripts/PopulationGenetics/labels_with_landcover2.csv")
RABIES_FILE <- file.path(ROOT_DIR, "scripts/PopulationGenetics/AL_rabies_2009_2025.csv")

OUT_DIR <- file.path(ROOT_DIR, "results/rabies_epi")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

OUT_JOINED_CSV <- file.path(OUT_DIR, "county_rabies_landcover_joined.csv")
OUT_STATS_TXT  <- file.path(OUT_DIR, "rabies_suitability_stats.txt")
FIG_MAP_SUIT_PNG    <- file.path(OUT_DIR, "map_movement_suitability_by_county.png")
FIG_SCAT_NB_PNG     <- file.path(OUT_DIR, "rabies_vs_suitability_nb.png")
FIG_COMBINED_PNG    <- file.path(OUT_DIR, "figure_suitability_map_plus_scatter.png")


# -------------------------
# helpers
# -------------------------
clean_county <- function(x) {
  x <- tolower(trimws(x))
  x <- gsub("[[:punct:]]", "", x)  # remove punctuation: "st. clair" -> "st clair"
  x <- gsub("\\s+", "", x)         # remove spaces: "st clair" -> "stclair"
  x
}

# -------------------------
# load data
# -------------------------
message("Reading labels: ", LABELS_FILE)
labs <- readr::read_csv(LABELS_FILE, show_col_types = FALSE) |> as.data.frame()
names(labs) <- tolower(names(labs))

message("Reading rabies totals: ", RABIES_FILE)
rabies <- readr::read_csv(RABIES_FILE, show_col_types = FALSE) |> as.data.frame()
names(rabies) <- tolower(names(rabies))

# Ensure county col exists
if (!("county" %in% names(labs))) stop("labels file missing column: county")
if (!("county" %in% names(rabies))) stop("rabies file missing column: county")

# rabies total col: assume the non-county column is the total
rabies_count_col <- setdiff(names(rabies), "county")
if (length(rabies_count_col) != 1) {
  stop("Rabies file should have exactly 2 columns (county + totals). Found: ",
       paste(names(rabies), collapse = ", "))
}
rabies <- rabies |> rename(total_cases = all_of(rabies_count_col))
rabies$total_cases <- as.numeric(rabies$total_cases)

# Keep ONE row per county from labels (centroid + landcover)
req_cols <- c(
  "county","lon","lat",
  "dev_open","dev_low","dev_med","dev_high",
  "water",
  "forest_decid","forest_ever","forest_mix",
  "wet_woody","wet_emerg",
  "agriculture"
)
miss <- setdiff(req_cols, names(labs))
if (length(miss) > 0) {
  stop("labels_with_landcover.csv missing expected columns: ", paste(miss, collapse=", "),
       "\n(If your landcover column names differ, tell me what they are and I’ll adjust.)")
}

county_tbl <- labs |>
  transmute(
    county = trimws(county),
    county_key = clean_county(county),
    lon = as.numeric(lon),
    lat = as.numeric(lat),

    # NLCD class-level %
    dev_open   = as.numeric(dev_open),
    dev_low    = as.numeric(dev_low),
    dev_med    = as.numeric(dev_med),
    dev_high   = as.numeric(dev_high),
    water      = as.numeric(water),

    forest_decid = as.numeric(forest_decid),
    forest_ever  = as.numeric(forest_ever),
    forest_mix   = as.numeric(forest_mix),

    wet_woody  = as.numeric(wet_woody),
    wet_emerg  = as.numeric(wet_emerg),

    agriculture = as.numeric(agriculture)
  ) |>
  distinct(county_key, .keep_all = TRUE) |>
  filter(!is.na(lon), !is.na(lat))

rabies2 <- rabies |>
  transmute(
    county = trimws(county),
    county_key = clean_county(county),
    total_cases = as.numeric(total_cases)
  )

dat <- county_tbl |>
  left_join(rabies2, by = "county_key", suffix = c("", "_rabies")) |>
  mutate(total_cases = ifelse(is.na(total_cases), 0, total_cases))

# -------------------------
# movement suitability index (simple + interpretable)
# -------------------------
# Rationale: raccoons are generalists; movement likely facilitated in mixed mosaics.
# This index emphasizes cover + edge/mosaic landscapes:
# forest + wetland + agriculture + developed (lightly) with mild penalty on developed to avoid "purely urban" dominance.
dat <- dat |>
  mutate(
    # combine related NLCD classes
    forest_total = forest_decid + forest_ever + forest_mix,
    wet_total    = wet_woody + wet_emerg,

    # high suitability: forest + wetlands + low/mod developed (edge/suburban)
    dev_low_med  = dev_open + dev_low + dev_med,

    # low suitability: high-intensity developed + open water
    low_barrier  = dev_high + water,

    # Option B scoring (simple, easy to explain)
    # high classes weighted more than moderate (agriculture), low classes given minimal/zero contribution
    suitability_raw = 2*(forest_total + wet_total + dev_low_med) + 1*(agriculture) + 0*(low_barrier),

    suitability_z = as.numeric(scale(suitability_raw))
  )

# write joined table
write.csv(dat, OUT_JOINED_CSV, row.names = FALSE)
message("Wrote joined county table: ", OUT_JOINED_CSV)

# -------------------------
# stats: correlation + simple models
# -------------------------
sink(OUT_STATS_TXT)
cat("Rabies totals (2009–2025) vs movement suitability (county-level)\n")
cat("Date:", format(Sys.time()), "\n\n")
cat("Counties in analysis:", nrow(dat), "\n\n")

cat("Spearman correlation (total_cases ~ suitability_raw):\n")
print(cor.test(dat$total_cases, dat$suitability_raw, method = "spearman"))
cat("\n")

cat("Pearson correlation (total_cases ~ suitability_raw):\n")
print(cor.test(dat$total_cases, dat$suitability_raw, method = "pearson"))
cat("\n")

cat("Poisson GLM (counts ~ suitability_raw):\n")
pois <- glm(total_cases ~ suitability_raw, data = dat, family = poisson())
print(summary(pois))
cat("\n")

# If MASS is available, fit negative binomial (often better for count overdispersion)
if (requireNamespace("MASS", quietly = TRUE)) {
  cat("Negative binomial GLM (counts ~ suitability_raw):\n")
  nb <- MASS::glm.nb(total_cases ~ suitability_raw, data = dat)
  print(summary(nb))
  cat("\n")
} else {
  cat("MASS not installed; skipped negative binomial model.\n\n")
}
sink()
message("Wrote stats: ", OUT_STATS_TXT)

# -------------------------
# Figure: scatter (rabies vs suitability)
# -------------------------
# -------------------------
# Figure: scatter with Negative Binomial fit (matches reported model)
# -------------------------
if (!requireNamespace("MASS", quietly = TRUE)) {
  stop("MASS package not installed; cannot fit negative binomial model.")
}

nb <- MASS::glm.nb(total_cases ~ suitability_raw, data = dat)

# prediction grid across observed range
grid <- data.frame(
  suitability_raw = seq(min(dat$suitability_raw, na.rm = TRUE),
                        max(dat$suitability_raw, na.rm = TRUE),
                        length.out = 200)
)

# predict on link scale, then back-transform to counts
pred <- predict(nb, newdata = grid, type = "link", se.fit = TRUE)
grid$fit_link <- pred$fit
grid$se_link  <- pred$se.fit

grid$fit <- exp(grid$fit_link)
grid$lwr <- exp(grid$fit_link - 1.96 * grid$se_link)
grid$upr <- exp(grid$fit_link + 1.96 * grid$se_link)

p_scatter_nb <- ggplot(dat, aes(x = suitability_raw, y = total_cases)) +
  geom_point(size = 3, alpha = 0.85, color = "#9A649A") +
  geom_ribbon(data = grid, aes(x = suitability_raw, ymin = lwr, ymax = upr),
              inherit.aes = FALSE, alpha = 0.25, fill = "#D9D9D9") +
 geom_line(data = grid, aes(x = suitability_raw, y = fit),
          inherit.aes = FALSE, linewidth = 1.1, color = "black") +
  theme_classic() +
  labs(
    x = "Movement suitability index (county)",
    y = "Total rabies cases (2009–2025)"
  )

ggsave(FIG_SCAT_NB_PNG, p_scatter_nb, width = 7.5, height = 5.2, dpi = 300)
message("Saved: ", FIG_SCAT_NB_PNG)

# -------------------------
# Figure A: Alabama county choropleth of movement suitability
# -------------------------
did_map <- FALSE
p_map_suit <- NULL

if (requireNamespace("sf", quietly = TRUE) && requireNamespace("tigris", quietly = TRUE)) {
  suppressPackageStartupMessages({
    library(sf)
    library(tigris)
  })
  options(tigris_use_cache = TRUE)

  message("Building Alabama county polygons via tigris...")
  al <- tigris::counties(state = "AL", year = 2020, class = "sf") |>
    st_transform(4326)

  al$key <- clean_county(al$NAME)

 al2 <- al |>
  left_join(dat |> select(county_key, total_cases, suitability_raw), by = c("key" = "county_key"))

p_map_suit <- ggplot(al2) +
  geom_sf(aes(fill = suitability_raw), color = "white", linewidth = 0.25) +
  geom_sf(fill = NA, color = "grey55", linewidth = 0.35) +
  theme_void() +
  labs(
    fill = "Movement\nsuitability",
    title = "County-level movement suitability"
  ) +
  scale_fill_gradient(
    low = "#FE8A03",
    high = "#4058C8",
    na.value = "grey85"
  )

  ggsave(FIG_MAP_SUIT_PNG, p_map_suit, width = 6.5, height = 6.2, dpi = 300)
  message("Saved: ", FIG_MAP_SUIT_PNG)
  did_map <- TRUE
} else {
  message("sf/tigris not available; skipping map.")
}

# -------------------------
# Combine map + scatter into one figure
# -------------------------
if (!requireNamespace("patchwork", quietly = TRUE)) {
  install.packages("patchwork")
}
library(patchwork)

if (!is.null(p_map_suit)) {
  p_combined <- p_map_suit + p_scatter_nb +
    plot_annotation(tag_levels = "A")

  ggsave(FIG_COMBINED_PNG, p_combined, width = 12.5, height = 6.2, dpi = 300)
  message("Saved combined figure: ", FIG_COMBINED_PNG)
}

message("DONE.")