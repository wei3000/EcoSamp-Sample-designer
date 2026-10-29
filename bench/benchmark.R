################################################################################
# EcoSamp benchmark
#
# Times each stage of ecosamp() and records peak memory, on the Vale do Candiru
# test data. Every scenario runs in a FRESH R process (via callr) so peak memory
# and timings are not polluted by earlier runs.
#
# It loads the package from your local git clone with pkgload::load_all(), so
# whatever is checked out (main = Hans's original, dev = my changes) is what
# gets benchmarked. Each row in the results CSV records the git branch + commit.
#
# Lives in the repo at bench/benchmark.R. Test data and the elevation map stay
# OUTSIDE the repo (in ROOT) and are never committed.
#
# Usage (RStudio): set ROOT below, then Source the file (Ctrl+Shift+S).
# Results (git-ignored): bench/results/bench_summary.csv  (one row per run)
#                        bench/results/bench_stages.csv   (one row per stage)
#                        bench/results/points_<scenario>_<commit>_rep<k>.rds
################################################################################

## ---- Settings ---------------------------------------------------------------
ROOT <- Sys.getenv("ECOSAMP_ROOT",
                   "C:/Computer/Imperial/Year 2/Work/UROP/EcoSamp")
REPO     <- file.path(ROOT, "EcoSamp-Sample-designer")      # git clone
PKG_DIR  <- file.path(REPO, "ecosamp")
DATA_DIR <- file.path(ROOT, "TestData_VdC")                   # not in repo
ELEV_TIF <- file.path(ROOT, "BRA_alt", "elevation", "BRA_elv_msk.tif")  # not in repo
OUT_DIR  <- file.path(REPO, "bench", "results")               # git-ignored

SEED <- 1
REPS <- 3           # repeats per scenario (raise to 3 for final numbers)
PROFILE <- FALSE    # TRUE = also write an Rprof hotspot summary per scenario

# Hans's trial-script parameters (with max_distance -> max_distance_sampling fixed)
base_args <- list(
  n_total = 15, point_number_min = 2, point_dist_min = 500,
  map_res = 5, sd_crs = "EPSG:32723",
  edge_dist_min = 100,
  road_dist_min = 50, road_dist_max = 500,
  water_dist_min = 100,
  slope_max = 30,
  max_distance_sampling = FALSE, area_weighting = TRUE,
  plot_results = FALSE
)

# Comment scenarios in/out as needed. Each entry overrides base_args.
scenarios <- list(
  trial         = list(),                                 # Hans's exact setup
  trial_maxdist = list(max_distance_sampling = TRUE),     # max-distance mode
  n30           = list(n_total = 30),                     # loop scaling (points)
  res2.5        = list(map_res = 2.5)                     # raster scaling (4x cells)
)
# Run only some scenarios, e.g. ONLY <- c("trial"); NULL = run all
ONLY <- NULL
if (nzchar(Sys.getenv("ECOSAMP_ONLY"))) ONLY <- strsplit(Sys.getenv("ECOSAMP_ONLY"), ",")[[1]]
if (!is.null(ONLY)) scenarios <- scenarios[ONLY]

## ---- One run (executed inside a fresh R process) ----------------------------
run_one <- function(pkg_dir, data_dir, elev_tif, args, seed, profile, prof_file) {
  t_start <- Sys.time()
  # The package NAMESPACE imports nothing, so ecosamp() relies on these being
  # attached (crs(), ext(), values(), n(), %>% ...). Mirror the trial script.
  suppressPackageStartupMessages({
    library(methods); library(sf); library(terra); library(raster)
    library(dplyr);   library(sp); library(fasterize)
  })
  pkgload::load_all(pkg_dir, quiet = TRUE, export_all = FALSE)

  peak_rss_mb <- function() {           # process peak memory, incl. C++ (terra/GDAL)
    if (.Platform$OS.type == "windows") {
      mi <- ps::ps_memory_info(ps::ps_handle())
      return(unname(mi[["peak_wset"]]) / 1024^2)
    }
    s <- readLines("/proc/self/status")
    as.numeric(gsub("\\D", "", grep("^VmHWM", s, value = TRUE))) / 1024
  }

  ## Inputs, prepared exactly as in ecosamp_test_trial.R
  t0 <- Sys.time()
  lc <- sf::read_sf(dsn = data_dir, layer = "vale_candiru_v3_05012023")
  lc <- sf::st_zm(lc, drop = TRUE, what = "ZM")
  lc$Treatment <- ifelse(grepl("0", lc$Classe), "Control", "FullPlantation")
  lc$Habitat   <- substr(lc$Classe, start = 3, stop = nchar(lc$Classe))
  map_habitat <- lc %>% dplyr::select(Habitat, geometry)
  map_treatmt <- lc %>% dplyr::select(Treatment, geometry)
  map_roads   <- sf::read_sf(dsn = data_dir, layer = "valecandiru_roads")
  map_waters  <- sf::read_sf(dsn = data_dir, layer = "waters_valecandiru2")
  map_elevation <- terra::rast(elev_tif)   # same object geodata::elevation_30s returns

  # Name-matched index (same numbers as the trial script, but robust to row order)
  hab_idx <- c("Infraestrutura" = 0, "APP" = 5.1, "Pasto sujo alaga" = 1.2,
               "Estrada" = 0, "Rede el\u00e9trica" = 0, "Floresta" = 5,
               "Pasto sujo desmatamento" = 1.1, "Pasto limpo desmatamento" = 2.1,
               "Pasto limpo alaga" = 2.2, "Pasto limpo" = 2, "Pasto sujo" = 1,
               "Juquira velha" = 4, "Juquira nova" = 3)
  u <- unique(map_habitat$Habitat)
  if (!all(u %in% names(hab_idx))) stop("Unmapped habitat(s): ",
                                        paste(setdiff(u, names(hab_idx)), collapse = ", "))
  index_habitat <- data.frame(Field = u, Index = unname(hab_idx[u]))
  index_treatmt <- data.frame(Field = c("Control", "FullPlantation"), Index = c(1, 2))
  t_inputs <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

  ## Stage timing: ecosamp() print()s a status line at the start of each stage.
  ## Run a copy of ecosamp() whose print() also timestamps those lines.
  ## (Non-invasive: the package code is not modified.)
  stamps <- data.frame(label = character(), t = numeric())
  f <- ecosamp
  shim <- new.env(parent = environment(f))
  shim$print <- function(x, ...) {
    if (is.character(x) && length(x) == 1)
      stamps[nrow(stamps) + 1, ] <<- list(x, as.numeric(Sys.time()))
    base::print(x, ...)
  }
  environment(f) <- shim

  call_args <- c(args, list(map_habitat = map_habitat, index_habitat = index_habitat,
                            map_treatmt = map_treatmt, index_treatmt = index_treatmt,
                            map_roads = map_roads, map_waters = map_waters,
                            map_elevation = map_elevation))

  gc(reset = TRUE)
  set.seed(seed)
  if (profile) Rprof(prof_file, interval = 0.02)
  t_call <- as.numeric(Sys.time())
  err <- NULL
  out <- tryCatch(do.call(f, call_args), error = function(e) { err <<- conditionMessage(e); NULL })
  t_end <- as.numeric(Sys.time())
  if (profile) Rprof(NULL)
  gcm <- gc()

  ## Collapse print() timestamps into named stages
  stage_of <- function(lbl) {
    if (grepl("^Current CRS", lbl))                   return("1_setup_rasterize")
    if (grepl("^CRS unit", lbl))                      return(NA_character_)
    if (grepl("Processing habitat", lbl))             return("2_reproject_inputs")
    if (grepl("calculating slope", lbl))              return("2b_slope")
    if (grepl("^Optional input", lbl))                return("2_reproject_inputs")
    if (grepl("distance to roads", lbl))              return("3_dist_roads")
    if (grepl("distance to waters", lbl))             return("4_dist_waters")
    if (grepl("habitat edge", lbl))                   return("5_dist_edge")
    if (grepl("Generating sample", lbl))              return("6_sampling_loop")
    if (grepl("^Complete", lbl))                      return("7_output")
    if (grepl("^Plotting", lbl))                      return("8_plot")
    NA_character_
  }
  # Everything before the first print (rasterizing habitat/treatment) is stage 1;
  # the gap between "Optional input"/edge stage and "Generating" includes masking.
  st <- stamps
  st$stage <- vapply(st$label, stage_of, character(1))
  st <- st[!is.na(st$stage), ]
  st$stage[1] <- "1_setup_rasterize"
  times <- data.frame(stage = st$stage,
                      secs = diff(c(t_call, st$t[-1], t_end)))
  times <- stats::aggregate(secs ~ stage, times, sum)

  ## Constraint checks on the output (validates any future rewrite too)
  checks <- list()
  if (!is.null(out)) {
    xy <- as.matrix(out[, c("x", "y")])
    d  <- stats::dist(xy)
    tol <- args$map_res                     # allow one pixel of discretisation
    checks$n_points          <- nrow(out)
    checks$min_pair_dist_m   <- if (length(d)) round(min(d), 1) else NA
    checks$ok_point_dist     <- !length(d) || min(d) >= args$point_dist_min - tol
    if ("Dist_road"  %in% names(out)) checks$ok_road  <- all(out$Dist_road >= args$road_dist_min - tol &
                                                             out$Dist_road <= args$road_dist_max + tol)
    if ("Dist_water" %in% names(out)) checks$ok_water <- all(out$Dist_water >= args$water_dist_min - tol)
    if ("Dist_edge"  %in% names(out)) checks$ok_edge  <- all(out$Dist_edge >= args$edge_dist_min - tol)
    if ("Slope"      %in% names(out)) checks$ok_slope <- all(out$Slope <= args$slope_max + 1e-6, na.rm = TRUE)
    bins <- table(paste(out$Habitat, out$Treatment, sep = " | "))
    checks$n_bins_sampled    <- length(bins)
    checks$min_points_in_bin <- min(bins)
    checks$bins              <- paste(sprintf("%s=%d", names(bins), bins), collapse = "; ")
  }

  hot <- NULL
  if (profile) {
    s <- summaryRprof(prof_file)$by.total
    s <- s[!grepl("^\"(do.call|f|tryCatch|tryCatchList|tryCatchOne|doTryCatch|run_one)\"$",
                  rownames(s)), ]
    hot <- utils::head(s, 25)
  }

  list(times = times, inputs_secs = t_inputs,
       total_secs = t_end - t_call,
       peak_r_heap_mb = sum(gcm[, ncol(gcm)]),   # "max used" (Mb), R objects only
       peak_rss_mb = peak_rss_mb(),             # whole process, incl. terra/GDAL
       checks = checks, error = err, points = out, hot = hot)
}

## ---- Driver ------------------------------------------------------------------
for (p in c("callr", "pkgload", "ps")) if (!requireNamespace(p, quietly = TRUE))
  stop("Please install.packages('", p, "')")
stopifnot(dir.exists(PKG_DIR), dir.exists(DATA_DIR), file.exists(ELEV_TIF))
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

git <- function(...) tryCatch(system2("git", c("-C", shQuote(PKG_DIR), ...),
                                      stdout = TRUE, stderr = FALSE)[1],
                              error = function(e) NA, warning = function(w) NA)
commit <- git("rev-parse", "--short", "HEAD")
branch <- git("rev-parse", "--abbrev-ref", "HEAD")
dirty  <- length(suppressWarnings(system2("git", c("-C", shQuote(PKG_DIR), "status",
                                  "--porcelain", "--", "R"), stdout = TRUE))) > 0
tag <- paste0(commit, if (isTRUE(dirty)) "-dirty")
cat(sprintf("Benchmarking %s @ %s (%s)\n\n", PKG_DIR, tag, branch))

rows <- list(); summary_rows <- list()
for (sc in names(scenarios)) for (rep in seq_len(REPS)) {
  args <- utils::modifyList(base_args, scenarios[[sc]])
  cat(sprintf("== %s (rep %d/%d) ... ", sc, rep, REPS)); flush.console()
  prof_file <- file.path(OUT_DIR, sprintf("rprof_%s_%s.out", sc, tag))
  res <- callr::r(run_one,
                  args = list(PKG_DIR, DATA_DIR, ELEV_TIF, args, SEED + rep - 1,
                              PROFILE, prof_file),
                  show = FALSE)
  if (!is.null(res$error)) {
    cat("ERROR:", res$error, "\n")
  } else {
    cat(sprintf("%.1f s, peak RSS %.0f MB, %d points\n",
                res$total_secs, res$peak_rss_mb, res$checks$n_points))
  }
  meta <- data.frame(time = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
                     branch = branch, commit = tag, scenario = sc, rep = rep,
                     seed = SEED + rep - 1, os = R.version$platform)
  rows[[length(rows) + 1]] <- cbind(meta, rbind(res$times,
      data.frame(stage = "TOTAL", secs = res$total_secs)))
  ch <- res$checks
  summary_rows[[length(summary_rows) + 1]] <- cbind(meta, data.frame(
      total_secs = round(res$total_secs, 2),
      peak_rss_mb = round(res$peak_rss_mb), peak_r_heap_mb = round(res$peak_r_heap_mb),
      n_points = if (is.null(ch$n_points)) NA else ch$n_points,
      min_pair_dist_m = if (is.null(ch$min_pair_dist_m)) NA else ch$min_pair_dist_m,
      constraints_ok = if (is.null(res$error)) all(unlist(ch[grepl("^ok_", names(ch))])) else NA,
      min_points_in_bin = if (is.null(ch$min_points_in_bin)) NA else ch$min_points_in_bin,
      error = if (is.null(res$error)) "" else res$error))
  if (!is.null(res$points))
    saveRDS(res$points, file.path(OUT_DIR, sprintf("points_%s_%s_rep%d.rds", sc, tag, rep)))
  if (!is.null(res$hot)) {
    cat("   top hotspots (total secs incl. callees):\n")
    print(utils::head(res$hot[, c("total.time", "total.pct")], 12))
    utils::write.csv(res$hot, file.path(OUT_DIR, sprintf("hotspots_%s_%s.csv", sc, tag)))
  }
}

stages  <- do.call(rbind, rows)
summ    <- do.call(rbind, summary_rows)
append_csv <- function(df, f) utils::write.table(df, f, sep = ",", row.names = FALSE,
                                                 col.names = !file.exists(f), append = file.exists(f))
append_csv(stages, file.path(OUT_DIR, "bench_stages.csv"))
append_csv(summ,   file.path(OUT_DIR, "bench_summary.csv"))

cat("\n---- Summary ----\n")
print(summ[, c("scenario", "rep", "total_secs", "peak_rss_mb", "n_points",
               "min_pair_dist_m", "constraints_ok", "min_points_in_bin")], row.names = FALSE)
cat("\n---- Seconds per stage ----\n")
wide <- stats::reshape(stats::aggregate(secs ~ scenario + stage, stages, mean),
                       idvar = "stage", timevar = "scenario", direction = "wide")
names(wide) <- sub("^secs\\.", "", names(wide))
wide[, -1] <- round(wide[, -1], 2)
print(wide[order(wide$stage), ], row.names = FALSE)
cat("\nWritten to", OUT_DIR, "\n")
