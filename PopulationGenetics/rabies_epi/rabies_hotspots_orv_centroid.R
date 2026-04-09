#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(terra)
})

ROOT_DIR <- "/scratch/Raccoon_andrea/gabe_trial"

LABELS <- file.path(ROOT_DIR, "scripts/PopulationGenetics/labels_with_landcover2.csv")
RABIES <- file.path(ROOT_DIR, "scripts/PopulationGenetics/AL_rabies_2009_2025.csv")

OUT_DIR <- file.path(ROOT_DIR, "results/rabies_epi")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

OUT_TXT    <- file.path(OUT_DIR, "orv_distance_models.txt")
OUT_CSV    <- file.path(OUT_DIR, "county_rabies_orv_distance.csv")
FIG_BUBBLE <- file.path(OUT_DIR, "rabies_centroid_bubble_orv.png")

# NEW figures
FIG_SCAT_DIST     <- file.path(OUT_DIR, "rabies_vs_dist_to_orv_scatter_NB.png")
FIG_SCAT_SUIT_ORV <- file.path(OUT_DIR, "rabies_vs_suitability_by_orv_scatter_NB.png")

# Toggle: plot y-axis as raw total cases, or log(1 + cases)
USE_LOG_Y <- FALSE

clean_county <- function(x) {
  x <- tolower(trimws(x))
  x <- gsub("[[:punct:]]", "", x)
  x <- gsub("\\s+", "", x)
  x
}

message("Reading labels: ", LABELS)
labs <- read_csv(LABELS, show_col_types = FALSE) |> as.data.frame()
names(labs) <- tolower(names(labs))

# expected columns in labels file
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
  stop(
    "labels_with_landcover2.csv is missing expected columns: ",
    paste(miss, collapse = ", "),
    "\nIf your column names differ, paste names(labs) and I’ll adjust."
  )
}

county_tbl <- labs |>
  transmute(
    county = trimws(county),
    county_key = clean_county(county),
    lon = as.numeric(lon),
    lat = as.numeric(lat),

    dev_open = as.numeric(dev_open),
    dev_low  = as.numeric(dev_low),
    dev_med  = as.numeric(dev_med),
    dev_high = as.numeric(dev_high),
    water    = as.numeric(water),

    forest_decid = as.numeric(forest_decid),
    forest_ever  = as.numeric(forest_ever),
    forest_mix   = as.numeric(forest_mix),

    wet_woody = as.numeric(wet_woody),
    wet_emerg = as.numeric(wet_emerg),

    agriculture = as.numeric(agriculture)
  ) |>
  distinct(county_key, .keep_all = TRUE) |>
  filter(!is.na(lon), !is.na(lat))

message("Reading rabies totals: ", RABIES)
rab <- read_csv(RABIES, show_col_types = FALSE) |> as.data.frame()
names(rab) <- tolower(names(rab))

if (!("county" %in% names(rab))) stop("Rabies file must include a 'county' column.")

rabies_count_col <- setdiff(names(rab), "county")
if (length(rabies_count_col) != 1) {
  stop(
    "Rabies file should have exactly 2 columns (county + totals). Found: ",
    paste(names(rab), collapse = ", ")
  )
}

rab <- rab |>
  rename(total_cases = all_of(rabies_count_col)) |>
  transmute(
    county = trimws(county),
    county_key = clean_county(county),
    total_cases = as.numeric(total_cases)
  )

dat <- county_tbl |>
  left_join(rab, by = "county_key") |>
  mutate(total_cases = ifelse(is.na(total_cases), 0, total_cases))

# ----------------------------
# Suitability index (Option B)
# high=2, medium=1, low=0 (then rescale 0–1)
# ----------------------------
dat <- dat |>
  mutate(
    forest_total = forest_decid + forest_ever + forest_mix,
    wet_total    = wet_woody + wet_emerg,
    dev_lowish   = dev_open + dev_low + dev_med,

    suit_raw = 2*(forest_total + wet_total + dev_lowish) +
               1*(agriculture) +
               0*(dev_high + water),

    suitability = (suit_raw - min(suit_raw, na.rm = TRUE)) /
                  (max(suit_raw, na.rm = TRUE) - min(suit_raw, na.rm = TRUE)),

    log_cases = log1p(total_cases)
  )

# ----------------------------
# ORV ZONE CODING (county-level)
# ----------------------------
orv_ne <- c(
 "Blount", "Cherokee", "Cullman", "DeKalb",
  "Etowah", "Jackson", "Jefferson", "Marshall",
  "Morgan", "St. Clair", "Shelby"
)

# SW coastal bait zone (per your methods)
orv_sw <- c("Baldwin")

orv_ne_key <- clean_county(orv_ne)
orv_sw_key <- clean_county(orv_sw)

dat <- dat |>
  mutate(
    orv_ne  = as.integer(county_key %in% orv_ne_key),
    orv_sw  = as.integer(county_key %in% orv_sw_key),
    orv_any = as.integer(orv_ne == 1 | orv_sw == 1)
  )

# ----------------------------
# DISTANCE TO ORV (centroid-based; km)
# ----------------------------
pts <- terra::vect(dat, geom = c("lon","lat"), crs = "EPSG:4326")

min_dist_km <- function(all_pts, target_idx) {
  d <- terra::distance(all_pts, all_pts[target_idx]) # meters, matrix (n x k)
  if (is.null(dim(d))) d <- matrix(d, ncol = 1)
  apply(d, 1, min, na.rm = TRUE) / 1000
}

idx_ne <- which(dat$orv_ne == 1)
idx_sw <- which(dat$orv_sw == 1)

dat$dist_to_orv_ne_km  <- if (length(idx_ne) > 0) min_dist_km(pts, idx_ne) else NA_real_
dat$dist_to_orv_sw_km  <- if (length(idx_sw) > 0) min_dist_km(pts, idx_sw) else NA_real_
dat$dist_to_orv_any_km <- pmin(dat$dist_to_orv_ne_km, dat$dist_to_orv_sw_km, na.rm = TRUE)

write.csv(dat, OUT_CSV, row.names = FALSE)
message("Wrote: ", OUT_CSV)

# ----------------------------
# MODELS (Poisson + NB)
# ----------------------------
pois <- glm(total_cases ~ dist_to_orv_any_km + suitability, data = dat, family = poisson())

have_mass <- requireNamespace("MASS", quietly = TRUE)
if (have_mass) {
  nb <- MASS::glm.nb(total_cases ~ dist_to_orv_any_km + suitability, data = dat)
} else {
  nb <- NULL
}

# write stats output
sink(OUT_TXT)
cat("County-level ORV + rabies (2009–2025 totals)\n")
cat("Counties in analysis:", nrow(dat), "\n\n")

cat("Compare cases inside vs outside ORV_ANY (Wilcoxon):\n")
print(wilcox.test(total_cases ~ orv_any, data = dat))
cat("\n")

cat("Spearman: cases vs distance-to-ORV_ANY (km):\n")
print(cor.test(dat$total_cases, dat$dist_to_orv_any_km, method = "spearman"))
cat("\n")

cat("Poisson GLM: cases ~ distance_to_ORV_ANY + suitability\n")
print(summary(pois))
cat("\n")

if (!is.null(nb)) {
  cat("NegBin GLM: cases ~ distance_to_ORV_ANY + suitability\n")
  print(summary(nb))
  cat("\n")
} else {
  cat("MASS not installed; skipped NegBin.\n\n")
}

cat("NOTE: ORV zones are often placed where rabies risk is high, so association is descriptive,\n")
cat("not causal. Interpret as overlap with enzootic area/front.\n")
sink()
message("Wrote: ", OUT_TXT)

# ----------------------------
# FIGURE 1: centroid bubble map (shape = ORV in/out; size = cases)
# ----------------------------
p_bubble <- ggplot(dat, aes(lon, lat)) +
  geom_point(aes(size = total_cases, shape = factor(orv_any)), alpha = 0.75) +
  theme_classic() +
  labs(
    x = "Longitude", y = "Latitude",
    size = "Total cases\n(2009–2025)",
    shape = "In ORV zone\n(county-level)",
    title = "Rabies totals by study-county centroids",
    subtitle = "Point size = total cases; shape indicates counties in ORV zone lists"
  )

ggsave(FIG_BUBBLE, p_bubble, width = 7.5, height = 5.2, dpi = 300)
message("Saved: ", FIG_BUBBLE)

# ----------------------------
# Helper: build NB prediction curve + CI for plotting
# ----------------------------
make_pred_curve <- function(model, xname, xseq, hold_vals) {

  # build a data frame with one row per x value
  newdat <- data.frame(
    dist_to_orv_any_km = rep(hold_vals$dist_to_orv_any_km, length(xseq)),
    suitability        = rep(hold_vals$suitability,        length(xseq))
  )

  # set the x variable you’re sweeping over
  newdat[[xname]] <- xseq

  # predict on link scale so CI is correct
  pr <- predict(model, newdata = newdat, type = "link", se.fit = TRUE)

  newdat$fit_link <- pr$fit
  newdat$se_link  <- pr$se.fit

  # back-transform to mean counts
  newdat$fit <- exp(newdat$fit_link)
  newdat$lwr <- exp(newdat$fit_link - 1.96 * newdat$se_link)
  newdat$upr <- exp(newdat$fit_link + 1.96 * newdat$se_link)

  newdat
}

# choose model for plotting: NB preferred, Poisson fallback
plot_model <- if (!is.null(nb)) nb else pois
plot_model_name <- if (!is.null(nb)) "Negative binomial" else "Poisson"

# y variable for points
yvar <- if (USE_LOG_Y) "log_cases" else "total_cases"
ylab <- if (USE_LOG_Y) "log(1 + total rabies cases) (2009–2025)" else "Total rabies cases (2009–2025)"

# hold values at medians for conditional effect plots
hold <- list(
  dist_to_orv_any_km = median(dat$dist_to_orv_any_km, na.rm = TRUE),
  suitability        = median(dat$suitability, na.rm = TRUE)
)

# ----------------------------
# FIGURE 2: rabies vs distance to ORV (NB curve holding suitability constant)
# ----------------------------
xseq_dist <- seq(min(dat$dist_to_orv_any_km, na.rm = TRUE),
                 max(dat$dist_to_orv_any_km, na.rm = TRUE),
                 length.out = 200)

pred_dist <- make_pred_curve(plot_model, "dist_to_orv_any_km", xseq_dist, hold)

# transform curve if using log y
if (USE_LOG_Y) {
  pred_dist <- pred_dist |>
    mutate(fit = log1p(fit), lwr = log1p(lwr), upr = log1p(upr))
}

p_dist <- ggplot(dat, aes(x = dist_to_orv_any_km, y = .data[[yvar]], shape = factor(orv_any))) +
  geom_point(size = 3, alpha = 0.85) +
  geom_ribbon(data = pred_dist, aes(x = dist_to_orv_any_km, ymin = lwr, ymax = upr),
              inherit.aes = FALSE, alpha = 0.20) +
  geom_line(data = pred_dist, aes(x = dist_to_orv_any_km, y = fit),
            inherit.aes = FALSE, linewidth = 1) +
  theme_classic() +
  labs(
    x = "Distance to nearest ORV county centroid (km)",
    y = ylab,
    shape = "In ORV zone\n(county-level)",
    title = "Rabies burden vs. proximity to ORV counties",
    subtitle = paste0(plot_model_name, " fit (holding suitability at median)")
  )

ggsave(FIG_SCAT_DIST, p_dist, width = 7.5, height = 5.2, dpi = 300)
message("Saved: ", FIG_SCAT_DIST)

# ----------------------------
# FIGURE 3: rabies vs suitability (NB curve holding distance constant), faceted by ORV
# ----------------------------
xseq_suit <- seq(min(dat$suitability, na.rm = TRUE),
                 max(dat$suitability, na.rm = TRUE),
                 length.out = 200)

pred_suit <- make_pred_curve(plot_model, "suitability", xseq_suit, hold)

if (USE_LOG_Y) {
  pred_suit <- pred_suit |>
    mutate(fit = log1p(fit), lwr = log1p(lwr), upr = log1p(upr))
}

p_suit_orv <- ggplot(dat, aes(x = suitability, y = .data[[yvar]], color = factor(orv_any))) +
  
  # points
  geom_point(size = 3, alpha = 0.6) +
  
  # confidence ribbon
  geom_ribbon(
    data = pred_suit,
    aes(x = suitability, ymin = lwr, ymax = upr),
    inherit.aes = FALSE,
    fill = "grey60",
    alpha = 0.25
  ) +
  
  # fitted line
  geom_line(
    data = pred_suit,
    aes(x = suitability, y = fit),
    inherit.aes = FALSE,
    linewidth = 1,
    color = "black"
  ) +
  
  theme_classic() +
  
  scale_color_manual(
    values = c("0" = "#FC8A03", "1" = "#4058C8"),
    labels = c("0" = "Outside ORV", "1" = "Inside ORV"),
    name = "ORV zone"
  ) +
  
  labs(
    x = "Movement suitability index",
    y = ylab,
    title = "Rabies burden vs. movement suitability",
    subtitle = paste0(plot_model_name, " fit (holding distance-to-ORV at median)")
  )

ggsave(FIG_SCAT_SUIT_ORV, p_suit_orv, width = 8.5, height = 5.2, dpi = 300)
message("Saved: ", FIG_SCAT_SUIT_ORV)

message("DONE.")