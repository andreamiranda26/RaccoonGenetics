#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(maps)
  library(dplyr)
  library(readr)
})

ROOT_DIR <- "/scratch/Raccoon_andrea/gabe_trial"

LABELS_FILE <- file.path(ROOT_DIR, "scripts/PopulationGenetics/labels_with_landcover2.csv")
RABIES_FILE <- file.path(ROOT_DIR, "scripts/PopulationGenetics/AL_rabies_2009_2025.csv")

OUT_DIR <- file.path(ROOT_DIR, "results/spatial_epi")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

OUT_PNG <- file.path(OUT_DIR, "map_AL_rabies_ORV_points_baseMaps.png")

# ---- read data ----
labels <- read_csv(LABELS_FILE, show_col_types = FALSE) %>%
  mutate(
    county = trimws(county),
    lon = as.numeric(lon),
    lat = as.numeric(lat)
  )

rabies <- read_csv(RABIES_FILE, show_col_types = FALSE) %>%
  mutate(county = trimws(county))

# Expected rabies columns: county, total
if (!all(c("county","total") %in% names(rabies))) {
  stop("Rabies CSV must contain columns: county, total")
}

county_tbl <- labels %>%
  select(county, lon, lat) %>%
  distinct()

dat <- county_tbl %>%
  left_join(rabies, by="county") %>%
  mutate(total = ifelse(is.na(total), 0, total))

# ORV counties list you provided (county-level membership)
orv_any <- c(
  # NE / North-central zone (county list)
  "Blount","Cherokee","Cullman","DeKalb",
  "Etowah","Jackson","Jefferson","Marshall","Morgan","St. Clair","Shelby",
  # SW coastal zone (place list -> approximate counties; keep as best-effort)
  "Baldwin"
)

dat <- dat %>%
  mutate(
    orv_any = ifelse(county %in% orv_any, 1, 0)
  )

# ---- point sizes (sqrt scaling for readability) ----
max_cases <- max(dat$total, na.rm=TRUE)
cex_vals <- if (max_cases == 0) rep(1, nrow(dat)) else (sqrt(dat$total / max_cases) * 4 + 0.8)

# ---- plot ----
png(OUT_PNG, width=2000, height=1200, res=200)

# expand right margin (last value controls right side)
par(mar = c(2.5, 2.5, 2.5, 14))  # bigger right margin

maps::map("county", region="alabama", fill=TRUE, col="white", border="gray80", lwd=0.8)
maps::map("state", region="alabama", add=TRUE, col="gray30", lwd=1.8)

points(dat$lon, dat$lat,
       pch=21,
       bg=ifelse(dat$orv_any==1, "tomato3", "steelblue4"),
       col="white",
       cex=cex_vals)

title("Alabama raccoon rabies totals (2009–2025) by study counties")

# ---- LEGENDS (robust placement outside the map) ----
par(xpd = NA)  # allow drawing in the figure region (incl. margins)

usr <- par("usr")
x_leg <- usr[2] + 0.06 * diff(usr[1:2])
y_top <- usr[4]

# ORV legend (top)
legend(
  x = x_leg, y = y_top,
  legend = c("Not in ORV county list", "In ORV county list"),
  pt.bg  = c("steelblue4", "tomato3"),
  pch    = 21,
  pt.cex = 1.5,
  pt.lwd = 0.5,
  bty    = "n",
  title  = "ORV status",
  y.intersp = 1.1
)

# Size legend (give big bubbles more vertical breathing room)
size_vals <- c(0, 20, 40, 60)
size_cex  <- (sqrt(size_vals / max_cases) * 4 + 0.8)

legend(
  x = x_leg,
  y = y_top - 0.25 * diff(usr[3:4]),
  legend = as.character(size_vals),
  pt.bg  = "gray60",
  pch    = 21,
  pt.cex = size_cex,
  bty    = "n",
  title  = "Total cases\n(2009–2025)",
  y.intersp = 1.6,   # <-- increase spacing between rows (fixes overlap)
  x.intersp = 1.0,
  cex = 0.95         # optional: slightly smaller text so it breathes
)

dev.off()