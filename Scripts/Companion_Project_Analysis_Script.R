# ========================================================================
# COMPANION PROJECT - BLOCK-ADJUSTED RCBD REANALYSIS
# Sites fixed; blocks nested in sites; plots are experimental units.
# Plant/plug observations are nested within plots for secondary analyses.
# Logit GLMMs, beta-binomial seed fill, DHARMa diagnostics and sensitivity.
# DO NOT cite results until this script is executed on the final Excel data,
# model fits/convergence checked, and model-specific outputs inspected.
# ========================================================================

required_packages <- c(
  "tidyverse",
  "readxl",
  "glmmTMB",
  "emmeans",
  "performance",
  "broom",
  "broom.mixed",
  "writexl",
  "scales",
  "patchwork"
)

invisible(
  lapply(required_packages, library, character.only = TRUE)
)

options(
  contrasts = c("contr.sum", "contr.poly"),
  na.action = "na.exclude",
  dplyr.summarise.inform = FALSE
)

set.seed(20261120)

# Example:
# Rscript Companion_Project_Analysis_RCBD_Updated.R "Companion Project(20261120-154847).xlsx"
# Output is anchored to the folder containing this R script, NOT to getwd().
# To explicitly choose another output location, set COMPANION_OUTPUT_DIR.
# DHARMa is strongly recommended; missing/failed diagnostics are explicitly
# flagged but will not suppress the manuscript figures.
# Check output/session_info.txt before stating the exact R version used.

# Manual clean_names() replacement (janitor unavailable offline)
clean_names <- function(df) {
  nm <- names(df)
  nm <- gsub("[^A-Za-z0-9]+", "_", nm)
  nm <- gsub("([a-z])([A-Z])", "\\1_\\2", nm)
  nm <- tolower(nm)
  nm <- gsub("_+", "_", nm)
  nm <- gsub("^_|_$", "", nm)
  names(df) <- nm
  df
}

# Locate this script consistently when executed with Rscript or sourced in RStudio.
# If neither location can be detected, use getwd() with an explicit message.
script_directory <- function() {
  script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(script_arg)) {
    script_path <- sub("^--file=", "", script_arg[1])
    return(dirname(normalizePath(script_path, mustWork = TRUE)))
  }
  source_files <- lapply(sys.frames(), function(fr) fr$ofile)
  source_files <- Filter(function(x) !is.null(x) && length(x) == 1 &&
                           nzchar(x) && file.exists(x), source_files)
  if (length(source_files)) {
    return(dirname(normalizePath(tail(source_files, 1)[[1]], mustWork = TRUE)))
  }
  if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
    editor_path <- tryCatch(rstudioapi::getActiveDocumentContext()$path,
                            error = function(e) "")
    if (length(editor_path) == 1 && nzchar(editor_path) && file.exists(editor_path)) {
      return(dirname(normalizePath(editor_path, mustWork = TRUE)))
    }
  }
  warning("Could not determine the script location. Using working directory: ", getwd())
  getwd()
}

project_dir <- script_directory()
output_override <- Sys.getenv("COMPANION_OUTPUT_DIR", unset = "")
output_dir <- if (nzchar(output_override)) {
  path.expand(output_override)
} else {
  file.path(project_dir, "output")
}
# Fully qualify output paths so they cannot change with setwd().
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
output_dir <- normalizePath(output_dir, winslash = "/", mustWork = TRUE)
figures_dir <- file.path(output_dir, "figures")
tables_dir <- file.path(output_dir, "tables")
models_dir <- file.path(output_dir, "models")
for (d in c(figures_dir, tables_dir, models_dir)) {
  dir.create(d, showWarnings = FALSE, recursive = TRUE)
  if (!dir.exists(d)) stop("Unable to create output directory: ", d)
}
message("Working directory: ", normalizePath(getwd(), winslash = "/"))
message("Script/project directory: ", project_dir)
message("ALL FIGURES WILL BE SAVED TO: ", normalizePath(figures_dir, winslash = "/"))

# Test folder write access before running lengthy statistical models.
probe_file <- file.path(figures_dir, ".figure_write_test")
if (!file.create(probe_file)) {
  stop("Cannot write to figure output folder: ", figures_dir)
}
unlink(probe_file)

save_plot <- function(filename, plot, width = 8, height = 6, dpi = 300) {
  # basename avoids accidentally creating output/figures/output/figures/.
  target <- file.path(figures_dir, basename(filename))
  ggplot2::ggsave(filename = target, plot = plot, width = width,
                  height = height, dpi = dpi, bg = "white", device = "png")
  if (!file.exists(target) || is.na(file.info(target)$size) ||
      file.info(target)$size <= 0) {
    stop("Figure could not be verified after save: ", target)
  }
  message("FIGURE SAVED: ", normalizePath(target, winslash = "/"))
  invisible(target)
}

write_table <- function(x, filename) {
  readr::write_csv(x, file.path(tables_dir, filename), na = "")
}


format_p <- function(p) {
  if (is.na(p)) return("p = NA")
  if (p < 0.001) return("p < 0.001")
  paste0("p = ", formatC(p, format = "f", digits = 3))
}

# Build compact significance-letter displays from Tukey-adjusted pairwise
# comparisons. This avoids requiring multcompView and works for the three
# treatment levels used in this experiment. Treatments sharing a letter are
# not significantly different at alpha = 0.05.
make_cld_three <- function(pair_obj, emm_obj,
                           treatment_levels = c("N", "F", "P"),
                           alpha = 0.05) {
  pair_df <- as.data.frame(pair_obj)
  emm_df <- as.data.frame(emm_obj)
  
  if (!"site" %in% names(pair_df) || !"site" %in% names(emm_df)) {
    stop("Pairwise comparisons and EMMs must both contain a 'site' column.")
  }
  
  estimate_col <- intersect(c("emmean", "prob", "response", "rate"), names(emm_df))[1]
  if (is.na(estimate_col)) {
    stop("Could not identify the estimated-mean column in emmeans output.")
  }
  
  map_dfr(levels(dat$site), function(site_name) {
    pd <- pair_df |> filter(as.character(site) == site_name)
    ed <- emm_df |> filter(as.character(site) == site_name)
    
    sig <- matrix(
      FALSE,
      nrow = length(treatment_levels),
      ncol = length(treatment_levels),
      dimnames = list(treatment_levels, treatment_levels)
    )
    
    if (nrow(pd) > 0) {
      for (i in seq_len(nrow(pd))) {
        parts <- strsplit(as.character(pd$contrast[i]), "\\s[-/]\\s")[[1]]
        if (length(parts) == 2 && all(parts %in% treatment_levels)) {
          sig[parts[1], parts[2]] <- !is.na(pd$p.value[i]) && pd$p.value[i] < alpha
          sig[parts[2], parts[1]] <- sig[parts[1], parts[2]]
        }
      }
    }
    
    all_subsets <- unlist(
      lapply(seq_along(treatment_levels), function(k) {
        combn(treatment_levels, k, simplify = FALSE)
      }),
      recursive = FALSE
    )
    
    valid_subsets <- Filter(function(z) {
      if (length(z) == 1) return(TRUE)
      all(combn(z, 2, function(pair) !sig[pair[1], pair[2]]))
    }, all_subsets)
    
    maximal_subsets <- Filter(function(z) {
      !any(vapply(valid_subsets, function(other) {
        length(other) > length(z) && all(z %in% other)
      }, logical(1)))
    }, valid_subsets)
    
    # Order letters so the clique containing the largest estimated mean is "a".
    est_lookup <- setNames(ed[[estimate_col]], as.character(ed$treatment))
    clique_score <- vapply(maximal_subsets, function(z) {
      max(est_lookup[z], na.rm = TRUE)
    }, numeric(1))
    maximal_subsets <- maximal_subsets[order(clique_score, decreasing = TRUE)]
    
    letter_names <- letters[seq_along(maximal_subsets)]
    group_letters <- setNames(rep("", length(treatment_levels)), treatment_levels)
    
    for (j in seq_along(maximal_subsets)) {
      for (trt in maximal_subsets[[j]]) {
        group_letters[trt] <- paste0(group_letters[trt], letter_names[j])
      }
    }
    
    tibble(
      site = factor(site_name, levels = levels(dat$site)),
      treatment = factor(treatment_levels, levels = treatment_levels),
      significance_group = unname(group_letters[treatment_levels])
    )
  })
}

args <- commandArgs(trailingOnly = TRUE)
file_from_env <- Sys.getenv("COMPANION_DATA_FILE", unset = "")
if (length(args)) {
  # An explicitly supplied file must exist: never fall back to another version.
  data_file <- args[1]
  if (!file.exists(data_file) && file.exists(file.path(project_dir, data_file))) {
    data_file <- file.path(project_dir, data_file)
  }
  if (!file.exists(data_file)) stop("Supplied workbook not found: ", data_file)
} else if (nzchar(file_from_env)) {
  data_file <- file_from_env
  if (!file.exists(data_file) && file.exists(file.path(project_dir, data_file))) {
    data_file <- file.path(project_dir, data_file)
  }
  if (!file.exists(data_file)) stop("COMPANION_DATA_FILE not found: ", data_file)
} else {
  candidates <- c(
    "Companion Project(20261120-154847).xlsx",  # verified final workbook
    "Data/Companion Project(20261120-154847).xlsx",
    "Data/Companion Project.xlsx"  # legacy path, use only if final workbook absent
  )
  found_files <- candidates[file.exists(candidates)]
  if (!length(found_files)) {
    found_files <- file.path(project_dir, candidates)
    found_files <- found_files[file.exists(found_files)]
  }
  if (!length(found_files)) {
    stop("Final verified workbook not found. Pass its path as the first Rscript argument. ",
         "Do not silently substitute an older workbook.")
  }
  data_file <- found_files[1]
}
message("Reading: ", normalizePath(data_file))
raw <- read_excel(data_file, sheet = "Complete Data") |> clean_names()

needed_columns <- c("site", "block", "plot", "treatment", "plant", "mortality",
                    "produced_seed", "filled", "empty", "weight_g",
                    "filled_proportion", "seed_total")
missing_columns <- setdiff(needed_columns, names(raw))
if (length(missing_columns)) {
  stop("Workbook is missing required columns: ", paste(missing_columns, collapse = ", "))
}

# Rename field sites for publication-facing outputs.
# Raw data labels are retained only at import; all downstream tables, models,
# and figures use Panhandle and Peninsula.
raw <- raw |>
  mutate(
    site = recode(
      as.character(site),
      "Joe" = "Panhandle",
      "Terry" = "Peninsula",
      .default = as.character(site)
    )
  )

dat <- raw |>
  transmute(
    site = factor(site, levels = c("Panhandle", "Peninsula")),
    block = factor(as.character(block)),
    plot = as.integer(plot),
    treatment = factor(treatment, levels = c("N", "F", "P")),
    plant = as.integer(plant),
    seed_total = as.numeric(seed_total),
    empty = as.numeric(empty),
    filled = as.numeric(filled),
    mortality = as.integer(mortality),
    produced_seed = as.integer(produced_seed),
    filled_proportion = as.numeric(filled_proportion),
    weight_g = as.numeric(weight_g),
    plot_id = interaction(site, plot, drop = TRUE, sep = "_"),
    site_block = interaction(site, block, drop = TRUE, sep = "_"),
    survived = 1L - mortality,
    set_seed = produced_seed,
    positive_weight_g = weight_g,
    # Structural zero only when the plant did not produce seed. If a plant is
    # recorded as producing seed but its weight is missing, retain NA rather
    # than incorrectly converting the missing measurement to 0 g.
    yield_g = case_when(
      produced_seed == 0 ~ 0,
      produced_seed == 1 & !is.na(weight_g) ~ weight_g,
      TRUE ~ NA_real_
    )
  )

unique_plants <- dat |> distinct(site, plot, plant) |> nrow()

validation <- tibble(
  check = c(
    "Total plant records", "Unique plant IDs", "Sites", "Plots",
    "Treatments", "Plants per plot", "Mortality coded 0/1",
    "Produced seed coded 0/1", "Dead plants did not produce seed",
    "Non-seed producers lack weights", "Seed producers missing weights"
  ),
  result = as.character(c(
    nrow(dat), unique_plants, n_distinct(dat$site), n_distinct(dat$plot_id),
    n_distinct(dat$treatment), min(count(dat, plot_id)$n),
    all(dat$mortality %in% c(0, 1)),
    all(dat$produced_seed %in% c(0, 1)),
    all(dat$produced_seed[dat$mortality == 1] == 0),
    all(is.na(dat$weight_g[dat$produced_seed == 0])),
    sum(dat$produced_seed == 1 & is.na(dat$weight_g))
  ))
)
print(validation)
write_table(validation, "00_validation.csv")

data_issues <- dat |>
  filter(
    (produced_seed == 1 & is.na(weight_g)) |
      (mortality == 1 & produced_seed == 1) |
      (produced_seed == 0 & !is.na(weight_g))
  ) |>
  select(site, block, plot, treatment, plant, mortality, produced_seed,
         seed_total, empty, filled, filled_proportion, weight_g)

write_table(data_issues, "00_data_issues.csv")
if (nrow(data_issues) > 0) {
  warning(
    "Data consistency issues detected. Review output/tables/00_data_issues.csv. ",
    "Plots containing a seed-producing plant with missing seed weight are ",
    "excluded from the primary seed-yield model until the missing weight is resolved."
  )
}

# ---- STRICT RANDOMIZED COMPLETE BLOCK VALIDATION ----
# There must be six blocks per site, each with exactly one N, F and P plot,
# and exactly 10 planting-position records per plot. Plot IDs can repeat
# between sites but NEVER across blocks within a site.
design_check <- dat |>
  distinct(site, block, site_block, plot, plot_id, treatment) |>
  count(site, block, treatment, name = "n_plots") |>
  tidyr::complete(site, block, treatment, fill = list(n_plots = 0L)) |>
  arrange(site, block, treatment)

plot_check <- dat |>
  group_by(site, plot_id) |>
  summarise(n_records = n(), n_blocks = n_distinct(site_block),
            n_treatments = n_distinct(treatment), .groups = "drop")

block_check <- dat |>
  distinct(site, site_block, plot_id, treatment) |>
  group_by(site, site_block) |>
  summarise(n_plots = n(), n_treatments = n_distinct(treatment),
            n_n = sum(treatment == "N"), n_f = sum(treatment == "F"),
            n_p = sum(treatment == "P"), .groups = "drop")

write_table(design_check, "00_design_check.csv")
write_table(plot_check, "00_plot_consistency_check.csv")
write_table(block_check, "00_block_completeness.csv")
print(block_check)

assertion <- c(
  nrow(dat) == 360,
  n_distinct(dat$site) == 2,
  n_distinct(dat$plot_id) == 36,
  n_distinct(dat$site_block) == 12,
  nrow(block_check) == 12,
  all(block_check$n_plots == 3L & block_check$n_treatments == 3L &
        block_check$n_n == 1L & block_check$n_f == 1L & block_check$n_p == 1L),
  all(plot_check$n_records == 10L & plot_check$n_blocks == 1L &
        plot_check$n_treatments == 1L),
  nrow(dat |> distinct(site, plot, plant)) == 360L,
  !anyNA(dat$site) && !anyNA(dat$block) && !anyNA(dat$treatment),
  !anyNA(dat$mortality) && !anyNA(dat$produced_seed),
  all(dat$mortality %in% c(0L, 1L)),
  all(dat$produced_seed %in% c(0L, 1L)),
  all(dat$produced_seed[dat$mortality == 1L] == 0L),
  sum(dat$produced_seed == 1L & is.na(dat$weight_g)) == 0L
)
if (!all(assertion)) {
  stop("Final data / RCBD validation FAILED. Inspect output/tables/00_block_completeness.csv ",
       "and output/tables/00_plot_consistency_check.csv. Do not run inferential models.")
}
message("PASS: 360 plug records, 36 plots, 12 complete RCBD blocks, ",
        "N/F/P each represented once per block.")

# Track the exact validated input used to reproduce these results.
writeLines(c(paste("Input:", normalizePath(data_file)),
             paste("MD5:", unname(tools::md5sum(data_file))),
             paste("Analysis timestamp:", Sys.time())),
           file.path(output_dir, "input_manifest.txt"))

plant_treatment_summary <- dat |>
  group_by(treatment) |>
  summarise(
    planted = n(),
    n_survived = sum(survived),
    n_mortality = sum(mortality),
    n_seed_producers = sum(produced_seed),
    n_survived_no_seed = sum(survived == 1 & produced_seed == 0),
    survival_rate = mean(survived),
    mortality_rate = mean(mortality),
    produced_seed_rate = mean(produced_seed),
    total_measured_weight_g = sum(weight_g, na.rm = TRUE),
    n_seed_producers_missing_weight = sum(produced_seed == 1 & is.na(weight_g)),
    mean_weight_per_planted_g = mean(yield_g, na.rm = TRUE),
    median_weight_per_planted_g = median(yield_g, na.rm = TRUE),
    mean_weight_producers_g = mean(weight_g[produced_seed == 1], na.rm = TRUE),
    median_weight_producers_g = median(weight_g[produced_seed == 1], na.rm = TRUE),
    sd_weight_producers_g = sd(weight_g[produced_seed == 1], na.rm = TRUE),
    mean_filled_proportion = mean(filled_proportion[produced_seed == 1], na.rm = TRUE),
    sd_filled_proportion = sd(filled_proportion[produced_seed == 1], na.rm = TRUE)
  )

print(plant_treatment_summary)
write_table(plant_treatment_summary, "01_treatment_descriptive_statistics.csv")

site_treatment_summary <- dat |>
  group_by(site, treatment) |>
  summarise(
    planted = n(),
    n_survived = sum(survived),
    n_mortality = sum(mortality),
    n_seed_producers = sum(produced_seed),
    n_survived_no_seed = sum(survived == 1 & produced_seed == 0),
    mortality_rate = mean(mortality),
    survival_rate = mean(survived),
    produced_seed_rate = mean(produced_seed),
    total_measured_weight_g = sum(weight_g, na.rm = TRUE),
    n_seed_producers_missing_weight = sum(produced_seed == 1 & is.na(weight_g)),
    mean_weight_per_planted_g = mean(yield_g, na.rm = TRUE),
    mean_weight_producers_g = mean(weight_g[produced_seed == 1], na.rm = TRUE),
    mean_filled_proportion = mean(filled_proportion[produced_seed == 1], na.rm = TRUE)
  )

print(site_treatment_summary)
write_table(site_treatment_summary, "02_site_treatment_descriptive_statistics.csv")

plot_dat <- dat |>
  group_by(site, block, site_block, plot, plot_id, treatment) |>
  summarise(
    planted = n(),
    n_survived = sum(survived),
    n_mortality = sum(mortality),
    n_seed_producers = sum(produced_seed),
    n_survived_no_seed = sum(survived == 1 & produced_seed == 0),
    mortality_rate = mean(mortality),
    survival_rate = mean(survived),
    produced_seed_rate = mean(produced_seed),
    weight_complete = !any(produced_seed == 1 & is.na(weight_g)),
    total_weight_g = if_else(
      weight_complete,
      sum(if_else(produced_seed == 0, 0, weight_g), na.rm = TRUE),
      NA_real_
    ),
    mean_weight_per_planted_g = if_else(
      weight_complete,
      total_weight_g / planted,
      NA_real_
    ),
    mean_weight_producers_g = if_else(
      sum(produced_seed == 1 & !is.na(weight_g)) > 0,
      mean(weight_g[produced_seed == 1], na.rm = TRUE),
      NA_real_
    ),
    total_filled = sum(filled, na.rm = TRUE),
    total_empty = sum(empty, na.rm = TRUE),
    total_seed_examined = total_filled + total_empty,
    pooled_filled_proportion = if_else(
      total_seed_examined > 0, total_filled / total_seed_examined, NA_real_
    )
  )

write_table(plot_dat, "03_plot_level_dataset.csv")

plot_summary <- plot_dat |>
  group_by(site, treatment) |>
  summarise(
    n_plots = n(),
    n_plots_weight_complete = sum(weight_complete),
    mean_total_weight = mean(total_weight_g, na.rm = TRUE),
    sd_total_weight = sd(total_weight_g, na.rm = TRUE),
    se_total_weight = sd_total_weight / sqrt(n_plots_weight_complete),
    mean_mortality = mean(mortality_rate),
    sd_mortality = sd(mortality_rate),
    se_mortality = sd_mortality / sqrt(n_plots),
    mean_produced_seed = mean(produced_seed_rate),
    sd_produced_seed = sd(produced_seed_rate),
    se_produced_seed = sd_produced_seed / sqrt(n_plots),
    mean_filled = mean(pooled_filled_proportion, na.rm = TRUE),
    sd_filled = sd(pooled_filled_proportion, na.rm = TRUE),
    n_plots_with_seed = sum(!is.na(pooled_filled_proportion)),
    se_filled = sd_filled / sqrt(n_plots_with_seed)
  )
write_table(plot_summary, "04_plot_level_summary_statistics.csv")

# ---- COMMON STATISTICAL FUNCTIONS ----
# Nested models are compared with identical random-effects terms, data and
# likelihood family. LRT p-values are asymptotic and should be read alongside
# DHARMa diagnostics, convergence checks and sensitivity analyses.

save_lrt <- function(reduced, full, filename) {
  res <- anova(reduced, full)
  print(res)
  write_table(as.data.frame(res), filename)
  res
}

model_fit_status <- function(fit, label) {
  tibble(model = label,
         converged = isTRUE(fit$fit$convergence == 0),
         positive_definite_hessian = isTRUE(fit$sdr$pdHess),
         AIC = AIC(fit), logLik = as.numeric(logLik(fit)))
}

# Pearson dispersion on the initial fixed-effect GLMs detects extra-binomial
# variance before accounting for plot clustering. This is *not* a goodness-
# of-fit test of the fitted GLMM, for which DHARMa is used separately.
initial_binomial_dispersion <- function(response) {
  f <- reformulate(c("treatment * site", "site:block"), response = response)
  g <- glm(f, data = dat, family = binomial(link = "logit"))
  chisq <- sum(residuals(g, type = "pearson")^2)
  df <- df.residual(g)
  tibble(response = response, initial_model = "fixed-block binomial GLM",
         pearson_chisq = chisq, residual_df = df,
         pearson_dispersion = chisq / df,
         approximate_upper_tail_p = pchisq(chisq, df, lower.tail = FALSE))
}

# DHARMa simulation-based residual assessment of *fitted* mixed models.
# A small p-value indicates potential model inadequacy, not proof thereof.
run_dharma <- function(fit, label, nsim = 1000L) {
  placeholder <- function(reason) {
    warning("DHARMa diagnostics incomplete for ", label, ": ", reason,
            ". Figures will still be saved; do NOT finalize model inference.")
    result <- tibble(model = label, check = "diagnostics_incomplete",
                     statistic = NA_real_, p_value = NA_real_, method = reason)
    write_table(result, paste0("diagnostics_", label, ".csv"))
    result
  }
  if (!requireNamespace("DHARMa", quietly = TRUE)) {
    return(placeholder("DHARMa package not installed"))
  }
  tryCatch({
    sim <- DHARMa::simulateResiduals(fittedModel = fit, n = nsim, refit = FALSE)
    diag_file <- file.path(figures_dir, paste0("DHARMa_", label, ".png"))
    grDevices::png(diag_file, width = 2000, height = 1500, res = 220)
    tryCatch(plot(sim), finally = grDevices::dev.off())
    message("DIAGNOSTIC FIGURE SAVED: ", diag_file)
    tests <- list(
      dispersion = DHARMa::testDispersion(sim, plot = FALSE),
      uniformity = DHARMa::testUniformity(sim, plot = FALSE),
      outliers = DHARMa::testOutliers(sim, plot = FALSE),
      zero_inflation = DHARMa::testZeroInflation(sim, plot = FALSE)
    )
    out <- purrr::imap_dfr(tests, function(z, nm) {
      tibble(model = label, check = nm,
             statistic = if (length(z$statistic)) as.numeric(z$statistic[[1]]) else NA_real_,
             p_value = z$p.value, method = as.character(z$method))
    })
    write_table(out, paste0("diagnostics_", label, ".csv"))
    out
  }, error = function(e) placeholder(conditionMessage(e)))
}

# Hierarchical random effects: site_block (six blocks in each site),
# and plot_id (36 experimental plots), with 10 plug outcomes per plot.
# A block variance near zero is a boundary estimate; do not suppress warnings.
rcbd_random <- "(1 | site_block) + (1 | plot_id)"

# ---- MORTALITY (PLUG LEVEL) ----
mortality_model_site_only <- glmmTMB(
  mortality ~ site + (1 | site_block) + (1 | plot_id),
  family = binomial(link = "logit"), data = dat)
mortality_model_additive <- glmmTMB(
  mortality ~ treatment + site + (1 | site_block) + (1 | plot_id),
  family = binomial(link = "logit"), data = dat)
mortality_model_interaction <- glmmTMB(
  mortality ~ treatment * site + (1 | site_block) + (1 | plot_id),
  family = binomial(link = "logit"), data = dat)

mortality_treatment_test <- save_lrt(mortality_model_site_only, mortality_model_additive,
                                     "08b_mortality_overall_treatment_test.csv")
mortality_interaction_test <- save_lrt(mortality_model_additive, mortality_model_interaction,
                                       "05a_mortality_interaction_test.csv")
mortality_emm_site <- emmeans(mortality_model_interaction, ~ treatment | site, type = "response")
mortality_pairs_site <- pairs(mortality_emm_site, adjust = "tukey")
mortality_emm_overall <- emmeans(mortality_model_additive, ~ treatment, type = "response")
mortality_pairs_overall <- pairs(mortality_emm_overall, adjust = "tukey")
write_table(as.data.frame(mortality_emm_site), "05_mortality_emmeans_by_site.csv")
write_table(as.data.frame(mortality_pairs_site), "06_mortality_pairwise_by_site.csv")
write_table(as.data.frame(mortality_emm_overall), "07_mortality_emmeans_overall.csv")
write_table(as.data.frame(mortality_pairs_overall), "08_mortality_pairwise_overall.csv")

# ---- SEED-PRODUCTION SUCCESS (PLUG LEVEL; DISTINCT FROM MORTALITY) ----
produced_seed_model_site_only <- glmmTMB(
  produced_seed ~ site + (1 | site_block) + (1 | plot_id),
  family = binomial(link = "logit"), data = dat)
produced_seed_model_additive <- glmmTMB(
  produced_seed ~ treatment + site + (1 | site_block) + (1 | plot_id),
  family = binomial(link = "logit"), data = dat)
produced_seed_model_interaction <- glmmTMB(
  produced_seed ~ treatment * site + (1 | site_block) + (1 | plot_id),
  family = binomial(link = "logit"), data = dat)
produced_seed_treatment_test <- save_lrt(produced_seed_model_site_only,
                                         produced_seed_model_additive,
                                         "09e_produced_seed_overall_treatment_test.csv")
produced_seed_interaction_test <- save_lrt(produced_seed_model_additive,
                                           produced_seed_model_interaction,
                                           "09f_produced_seed_interaction_test.csv")
produced_seed_emm_site <- emmeans(produced_seed_model_interaction, ~ treatment | site,
                                  type = "response")
produced_seed_pairs_site <- pairs(produced_seed_emm_site, adjust = "tukey")
produced_seed_emm_overall <- emmeans(produced_seed_model_additive, ~ treatment,
                                     type = "response")
produced_seed_pairs_overall <- pairs(produced_seed_emm_overall, adjust = "tukey")
write_table(as.data.frame(produced_seed_emm_site), "09a_produced_seed_emmeans_by_site.csv")
write_table(as.data.frame(produced_seed_pairs_site), "09b_produced_seed_pairwise_by_site.csv")
write_table(as.data.frame(produced_seed_emm_overall), "09c_produced_seed_emmeans_overall.csv")
write_table(as.data.frame(produced_seed_pairs_overall), "09d_produced_seed_pairwise_overall.csv")

# ---- FILLED-SEED PROPORTION (ONLY PLUGS THAT PRODUCED SEED) ----
seed_dat <- dat |>
  filter(produced_seed == 1L, !is.na(filled), !is.na(empty), filled + empty > 0) |>
  droplevels()
if (!nrow(seed_dat)) stop("No eligible seed-fill observations available.")
if (any(seed_dat$filled < 0 | seed_dat$empty < 0)) {
  stop("Filled and empty seed counts must be nonnegative.")
}
quality_model_additive <- glmmTMB(
  cbind(filled, empty) ~ treatment + site + (1 | site_block) + (1 | plot_id),
  family = betabinomial(link = "logit"), data = seed_dat)
quality_model_interaction <- glmmTMB(
  cbind(filled, empty) ~ treatment * site + (1 | site_block) + (1 | plot_id),
  family = betabinomial(link = "logit"), data = seed_dat)
quality_interaction_test <- save_lrt(quality_model_additive, quality_model_interaction,
                                     "10a_filled_seed_interaction_test.csv")
quality_emm_site <- emmeans(quality_model_interaction, ~ treatment | site,
                            type = "response")
quality_pairs_site <- pairs(quality_emm_site, adjust = "tukey")
quality_emm_overall <- emmeans(quality_model_additive, ~ treatment,
                               type = "response")
quality_pairs_overall <- pairs(quality_emm_overall, adjust = "tukey")
write_table(as.data.frame(quality_emm_site), "10_filled_seed_emmeans_by_site.csv")
write_table(as.data.frame(quality_pairs_site), "11_filled_seed_pairwise_by_site.csv")
write_table(as.data.frame(quality_emm_overall), "12_filled_seed_emmeans_overall.csv")
write_table(as.data.frame(quality_pairs_overall), "13_filled_seed_pairwise_overall.csv")

# ---- MODEL DIAGNOSTICS, DISPERSION, AND FIT STATUS ----
initial_dispersion <- bind_rows(
  initial_binomial_dispersion("mortality"),
  initial_binomial_dispersion("produced_seed"))
write_table(initial_dispersion, "09_initial_binomial_overdispersion.csv")

mortality_diagnostics <- run_dharma(mortality_model_interaction, "mortality")
produced_seed_diagnostics <- run_dharma(produced_seed_model_interaction, "seed_production")
quality_diagnostics <- run_dharma(quality_model_interaction, "filled_seed")
secondary_model_status <- bind_rows(
  model_fit_status(mortality_model_site_only, "mortality_site_only"),
  model_fit_status(mortality_model_additive, "mortality_additive"),
  model_fit_status(mortality_model_interaction, "mortality_interaction"),
  model_fit_status(produced_seed_model_site_only, "produced_seed_site_only"),
  model_fit_status(produced_seed_model_additive, "produced_seed_additive"),
  model_fit_status(produced_seed_model_interaction, "produced_seed_interaction"),
  model_fit_status(quality_model_additive, "filled_seed_additive"),
  model_fit_status(quality_model_interaction, "filled_seed_interaction"))
write_table(secondary_model_status, "09_model_fit_status.csv")
residual_alerts <- bind_rows(mortality_diagnostics,
                             produced_seed_diagnostics,
                             quality_diagnostics) |>
  mutate(needs_review = check == "diagnostics_incomplete" |
           (!is.na(p_value) & p_value < 0.05))
write_table(residual_alerts, "09_residual_alerts.csv")
if (any(residual_alerts$needs_review)) {
  warning("DHARMa diagnostics are missing or detected at least one potential model violation. ",
          "Inspect output/tables/09_residual_alerts.csv before reporting results.")
}
if (!all(secondary_model_status$converged & secondary_model_status$positive_definite_hessian)) {
  warning("At least one mixed model has convergence or Hessian concerns; inspect 09_model_fit_status.csv.")
}

# ---- BETA-BINOMIAL SENSITIVITY ANALYSIS FOR BINARY OUTCOMES ----
# Plot-aggregated successes/trials with *fixed* blocks nested in site offers a
# contrast to the random-block, plot-random-intercept binomial GLMM.
# Fixed-block beta-binomial models cannot be interpreted as exact matches
# to the main binomial models, but can flag fragile likelihood-ratio p-values.
fit_betabinomial_sensitivity <- function(success, failure, label) {
  response <- paste0("cbind(", success, ", ", failure, ")")
  formula_add <- as.formula(paste(response, "~ treatment + site + site:block"))
  formula_int <- as.formula(paste(response, "~ treatment * site + site:block"))
  tryCatch({
    m0 <- glmmTMB(formula_add, family = betabinomial(link = "logit"), data = plot_dat)
    m1 <- glmmTMB(formula_int, family = betabinomial(link = "logit"), data = plot_dat)
    res <- save_lrt(m0, m1, paste0("sensitivity_", label, "_LRT.csv"))
    write_table(bind_rows(model_fit_status(m0, paste0(label, "_additive")),
                          model_fit_status(m1, paste0(label, "_interaction"))),
                paste0("sensitivity_", label, "_fit_status.csv"))
    run_dharma(m1, paste0("sensitivity_", label))
    invisible(list(additive = m0, interaction = m1, test = res))
  }, error = function(e) {
    warning("Beta-binomial sensitivity fit failed for ", label, ": ", conditionMessage(e))
    NULL
  })
}
plot_dat <- plot_dat |>
  mutate(n_alive = planted - n_mortality,
         n_no_seed = planted - n_seed_producers)
mortality_sensitivity <- fit_betabinomial_sensitivity("n_mortality", "n_alive", "mortality")
produced_seed_sensitivity <- fit_betabinomial_sensitivity("n_seed_producers", "n_no_seed",
                                                       "seed_production")

# ---- PRIMARY SEED-MASS MODEL (FIXED-BLOCK RCBD, PLOT AS UNIT) ----
# Distinguish true zeros (non-producing plugs) from missing weights.
weight_model_dat <- plot_dat |> filter(weight_complete, !is.na(total_weight_g))
if (nrow(weight_model_dat) != 36L) {
  stop("Primary model must include all 36 plots; missing yield must be resolved.")
}
weight_model_additive <- lm(log1p(total_weight_g) ~ treatment + site + site:block,
                            data = weight_model_dat)
weight_model_interaction <- lm(log1p(total_weight_g) ~ treatment * site + site:block,
                               data = weight_model_dat)
weight_interaction_test <- anova(weight_model_additive, weight_model_interaction)
weight_anova <- anova(weight_model_interaction)
print(weight_interaction_test)
print(weight_anova)
write_table(as.data.frame(weight_interaction_test), "15a_weight_interaction_test.csv")
write_table(as.data.frame(weight_anova), "15b_weight_anova.csv")
weight_emm_site <- emmeans(weight_model_interaction, ~ treatment | site)
weight_pairs_site <- pairs(weight_emm_site, adjust = "tukey")
weight_emm_overall <- emmeans(weight_model_additive, ~ treatment)
weight_pairs_overall <- pairs(weight_emm_overall, adjust = "tukey")
write_table(as.data.frame(weight_emm_site), "15_weight_emmeans_by_site_log_scale.csv")
write_table(as.data.frame(weight_pairs_site), "16_weight_pairwise_by_site.csv")
write_table(as.data.frame(weight_emm_overall), "17_weight_emmeans_overall_log_scale.csv")
write_table(as.data.frame(weight_pairs_overall), "18_weight_pairwise_overall.csv")
write_table(broom::tidy(weight_model_interaction), "19_weight_model_coefficients.csv")

weight_diag_file <- file.path(figures_dir, "Supplement_Weight_Model_Diagnostics.png")
grDevices::png(weight_diag_file, width = 2000, height = 1600, res = 250)
tryCatch({
  par(mfrow = c(2, 2))
  plot(weight_model_interaction)
}, finally = grDevices::dev.off())
message("DIAGNOSTIC FIGURE SAVED: ", weight_diag_file)
weight_shapiro <- shapiro.test(residuals(weight_model_interaction))
weight_diagnostics <- tibble(test = "Shapiro-Wilk residual normality",
                              statistic = unname(weight_shapiro$statistic),
                              p_value = weight_shapiro$p.value)
write_table(weight_diagnostics, "20_weight_model_diagnostics.csv")

# ---- CONDITIONAL SEED MASS AMONG SEED-PRODUCING PLUGS ----
positive_dat <- dat |>
  filter(produced_seed == 1L, !is.na(weight_g), weight_g > 0) |>
  droplevels()
positive_weight_model_additive <- glmmTMB(
  log(weight_g) ~ treatment + site + site:block + (1 | plot_id),
  family = gaussian(), data = positive_dat)
positive_weight_model_interaction <- glmmTMB(
  log(weight_g) ~ treatment * site + site:block + (1 | plot_id),
  family = gaussian(), data = positive_dat)
positive_weight_interaction_test <- save_lrt(
  positive_weight_model_additive, positive_weight_model_interaction,
  "21a_positive_weight_interaction_test.csv")
positive_weight_emm <- emmeans(positive_weight_model_interaction, ~ treatment | site)
positive_weight_pairs <- pairs(positive_weight_emm, adjust = "tukey")
write_table(as.data.frame(positive_weight_emm), "21_positive_weight_emmeans.csv")
write_table(as.data.frame(positive_weight_pairs), "22_positive_weight_pairwise.csv")
positive_weight_diagnostics <- run_dharma(positive_weight_model_interaction,
                                          "positive_seed_mass")
positive_weight_fit_status <- bind_rows(
  model_fit_status(positive_weight_model_additive, "positive_mass_additive"),
  model_fit_status(positive_weight_model_interaction, "positive_mass_interaction"))
write_table(positive_weight_fit_status, "21b_positive_weight_fit_status.csv")
if (!all(positive_weight_fit_status$converged & positive_weight_fit_status$positive_definite_hessian)) {
  warning("Conditional positive-mass model needs fit review.")
}

# ---- SIGNIFICANCE LETTERS FOR FIGURES ----
# Letters are based on Tukey-adjusted within-site pairwise comparisons.
# Treatments sharing a letter are not significantly different (alpha = 0.05).
weight_letters <- make_cld_three(weight_pairs_site, weight_emm_site)
mortality_letters <- make_cld_three(mortality_pairs_site, mortality_emm_site)
produced_seed_letters <- make_cld_three(produced_seed_pairs_site, produced_seed_emm_site)
quality_letters <- make_cld_three(quality_pairs_site, quality_emm_site)
positive_weight_letters <- make_cld_three(positive_weight_pairs, positive_weight_emm)

write_table(weight_letters, "25_weight_significance_letters.csv")
write_table(mortality_letters, "26_mortality_significance_letters.csv")
write_table(produced_seed_letters, "26b_produced_seed_significance_letters.csv")
write_table(quality_letters, "27_filled_seed_significance_letters.csv")
write_table(positive_weight_letters, "28_positive_weight_significance_letters.csv")

weight_interaction_p <- weight_interaction_test$`Pr(>F)`[2]
mortality_interaction_p <- mortality_interaction_test$`Pr(>Chisq)`[2]
produced_seed_interaction_p <- produced_seed_interaction_test$`Pr(>Chisq)`[2]
quality_interaction_p <- quality_interaction_test$`Pr(>Chisq)`[2]


# ---- FIGURE 1: Total seed weight ----
weight_letter_pos <- weight_model_dat |>
  group_by(site, treatment) |>
  summarise(y = max(total_weight_g, na.rm = TRUE)) |>
  left_join(weight_letters, by = c("site", "treatment")) |>
  group_by(site) |>
  mutate(
    panel_max = max(weight_model_dat$total_weight_g[weight_model_dat$site == first(site)], na.rm = TRUE),
    y = y + pmax(panel_max * 0.08, 2)
  ) |>
  ungroup()

fig_weight <- ggplot(weight_model_dat, aes(x = treatment, y = total_weight_g, fill = treatment)) +
  geom_boxplot(width = 0.55, alpha = 0.65, outlier.shape = NA) +
  geom_jitter(aes(shape = site), width = 0.10, size = 2.7, alpha = 0.8) +
  geom_text(
    data = weight_letter_pos,
    aes(x = treatment, y = y, label = significance_group),
    inherit.aes = FALSE,
    fontface = "bold",
    size = 5
  ) +
  facet_wrap(~ site, nrow = 1, scales = "free_y") +
  scale_y_continuous(expand = expansion(mult = c(0.03, 0.16))) +
  labs(
    title = NULL,
    subtitle = paste0(
      "Different letters indicate Tukey-adjusted differences within site; treatment x site ",
      format_p(weight_interaction_p)
    ),
    x = "Treatment",
    y = "Total harvested seed weight per plot (g)",
    shape = "Site"
  ) +
  theme_classic(base_size = 13) +
  theme(
    legend.position = "none",
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = 16),
    axis.title = element_text(face = "bold", size = 16),
    axis.text = element_text(size = 14)
  ) +
  guides(fill = "none")
print(fig_weight)
save_plot("Figure_1_Total_Seed_Weight.png", fig_weight, width = 9, height = 5.5)

# ---- FIGURE 2: Mortality ----
mortality_plot_dat <- plot_summary |>
  left_join(mortality_letters, by = c("site", "treatment")) |>
  mutate(label_y = pmin(mean_mortality + se_mortality + 0.08, 1.08))

fig_mortality <- ggplot(mortality_plot_dat, aes(x = treatment, y = mean_mortality, fill = site)) +
  geom_col(position = position_dodge(width = 0.75), width = 0.65) +
  geom_errorbar(
    aes(
      ymin = pmax(mean_mortality - se_mortality, 0),
      ymax = pmin(mean_mortality + se_mortality, 1)
    ),
    position = position_dodge(width = 0.75), width = 0.15
  ) +
  geom_point(
    data = plot_dat,
    aes(x = treatment, y = mortality_rate, group = site, fill = site),
    position = position_jitterdodge(jitter.width = 0.08, dodge.width = 0.75),
    inherit.aes = FALSE, alpha = 0.6, size = 1.8, shape = 21, color = "black"
  ) +
  geom_text(
    aes(y = label_y, label = significance_group, group = site),
    position = position_dodge(width = 0.75),
    fontface = "bold",
    size = 5,
    vjust = 0
  ) +
  scale_y_continuous(
    labels = percent_format(accuracy = 1),
    limits = c(0, 1.12),
    breaks = seq(0, 1, 0.2),
    expand = expansion(mult = c(0, 0.02))
  ) +
  labs(
    title = "Liatris plug mortality",
    subtitle = paste0(
      "Bars are plot-level means +/- SE; different letters indicate Tukey-adjusted differences within site; treatment x site ",
      format_p(mortality_interaction_p)
    ),
    x = "Treatment",
    y = "Mortality (%)",
    fill = "Site"
  ) +
  theme_classic(base_size = 13) +
  theme(
    legend.position = "top",
    axis.title = element_text(face = "bold", size = 16),
    axis.text = element_text(size = 14),
    legend.title = element_text(face = "bold", size = 15),
    legend.text = element_text(size = 14)
  )
print(fig_mortality)
save_plot("Figure_2_Mortality.png", fig_mortality, width = 8.5, height = 6)

# ---- SUPPLEMENT FIGURE: Probability of producing seed ----
produced_seed_plot_dat <- plot_summary |>
  left_join(produced_seed_letters, by = c("site", "treatment")) |>
  mutate(label_y = pmin(mean_produced_seed + se_produced_seed + 0.08, 1.08))

fig_produced_seed <- ggplot(
  produced_seed_plot_dat,
  aes(x = treatment, y = mean_produced_seed, fill = site)
) +
  geom_col(position = position_dodge(width = 0.75), width = 0.65) +
  geom_errorbar(
    aes(
      ymin = pmax(mean_produced_seed - se_produced_seed, 0),
      ymax = pmin(mean_produced_seed + se_produced_seed, 1)
    ),
    position = position_dodge(width = 0.75), width = 0.15
  ) +
  geom_point(
    data = plot_dat,
    aes(x = treatment, y = produced_seed_rate, group = site, fill = site),
    position = position_jitterdodge(jitter.width = 0.08, dodge.width = 0.75),
    inherit.aes = FALSE, alpha = 0.6, size = 1.8, shape = 21, color = "black"
  ) +
  geom_text(
    aes(y = label_y, label = significance_group, group = site),
    position = position_dodge(width = 0.75),
    fontface = "bold", size = 5, vjust = 0
  ) +
  scale_y_continuous(
    labels = percent_format(accuracy = 1), limits = c(0, 1.12),
    breaks = seq(0, 1, 0.2), expand = expansion(mult = c(0, 0.02))
  ) +
  labs(
    title = "Liatris seed-production success",
    subtitle = paste0(
      "Bars are plot-level means +/- SE; letters compare treatments within site; treatment x site ",
      format_p(produced_seed_interaction_p)
    ),
    x = "Treatment", y = "Plugs producing seed (%)", fill = "Site"
  ) +
  theme_classic(base_size = 13) +
  theme(legend.position = "top", axis.title = element_text(face = "bold"))

print(fig_produced_seed)
save_plot("Supplement_Figure_Seed_Production_Success.png", fig_produced_seed, width = 8.5, height = 6)

# ---- FIGURE 3: Filled-seed proportion ----
quality_plot_dat <- plot_summary |>
  left_join(quality_letters, by = c("site", "treatment")) |>
  mutate(label_y = pmin(mean_filled + se_filled + 0.08, 1.08))

fig_filled <- ggplot(quality_plot_dat, aes(x = treatment, y = mean_filled, fill = site)) +
  geom_col(position = position_dodge(width = 0.75), width = 0.65) +
  geom_errorbar(
    aes(
      ymin = pmax(mean_filled - se_filled, 0),
      ymax = pmin(mean_filled + se_filled, 1)
    ),
    position = position_dodge(width = 0.75), width = 0.15
  ) +
  geom_point(
    data = plot_dat |> filter(!is.na(pooled_filled_proportion)),
    aes(x = treatment, y = pooled_filled_proportion, group = site, fill = site),
    position = position_jitterdodge(jitter.width = 0.08, dodge.width = 0.75),
    inherit.aes = FALSE, alpha = 0.6, size = 1.8, shape = 21, color = "black"
  ) +
  geom_text(
    aes(y = label_y, label = significance_group, group = site),
    position = position_dodge(width = 0.75),
    fontface = "bold",
    size = 5,
    vjust = 0
  ) +
  scale_y_continuous(
    labels = percent_format(accuracy = 1),
    limits = c(0, 1.12),
    breaks = seq(0, 1, 0.2),
    expand = expansion(mult = c(0, 0.02))
  ) +
  labs(
    title = "Filled-seed proportion",
    subtitle = paste0(
      "Bars are plot-level means +/- SE; different letters indicate Tukey-adjusted differences within site; treatment x site ",
      format_p(quality_interaction_p)
    ),
    x = "Treatment",
    y = "Filled seeds (%)",
    fill = "Site"
  ) +
  theme_classic(base_size = 13) +
  theme(
    legend.position = "top",
    axis.title = element_text(face = "bold", size = 16),
    axis.text = element_text(size = 14),
    legend.title = element_text(face = "bold", size = 15),
    legend.text = element_text(size = 14)
  )
print(fig_filled)
save_plot("Figure_3_Filled_Seed_Proportion.png", fig_filled, width = 8.5, height = 6)

# ---- SUPPLEMENT FIGURE 1: Weight among reproductive plants ----
positive_weight_letter_pos <- positive_dat |>
  group_by(site, treatment) |>
  summarise(y = max(weight_g, na.rm = TRUE)) |>
  left_join(positive_weight_letters, by = c("site", "treatment")) |>
  group_by(site) |>
  mutate(
    panel_max = max(positive_dat$weight_g[positive_dat$site == first(site)], na.rm = TRUE),
    y = y + pmax(panel_max * 0.08, 0.5)
  ) |>
  ungroup()

fig_positive_weight <- positive_dat |>
  ggplot(aes(x = treatment, y = weight_g)) +
  geom_boxplot(width = 0.55, outlier.shape = NA) +
  geom_jitter(aes(shape = site), width = 0.12, alpha = 0.65, size = 2) +
  geom_text(
    data = positive_weight_letter_pos,
    aes(x = treatment, y = y, label = significance_group),
    inherit.aes = FALSE,
    fontface = "bold",
    size = 5
  ) +
  facet_wrap(~ site, nrow = 1, scales = "free_y") +
  scale_y_continuous(expand = expansion(mult = c(0.03, 0.16))) +
  labs(
    title = "Seed mass among reproductive Liatris plugs",
    subtitle = "Different letters indicate Tukey-adjusted differences within site; nonreproductive plugs excluded",
    x = "Treatment",
    y = "Harvested seed mass per plug (g)",
    shape = "Site"
  ) +
  theme_classic(base_size = 13) +
  theme(legend.position = "none", axis.title = element_text(face = "bold"),
        strip.background = element_blank(), strip.text = element_text(face = "bold"))
print(fig_positive_weight)
save_plot("Supplement_Figure_1_Weight_Among_Producers.png", fig_positive_weight, width = 9, height = 5.5)

# ---- SUPPLEMENT FIGURE 2: Site x treatment interaction ----
weight_interaction_summary <- weight_model_dat |>
  group_by(site, treatment) |>
  summarise(mean_weight = mean(total_weight_g), se_weight = sd(total_weight_g) / sqrt(n())) |>
  left_join(weight_letters, by = c("site", "treatment"))

fig_interaction <- ggplot(
  weight_interaction_summary,
  aes(x = treatment, y = mean_weight, group = site, shape = site)
) +
  geom_line(aes(linetype = site), linewidth = 0.8) +
  geom_point(size = 3.2) +
  geom_errorbar(
    aes(ymin = pmax(mean_weight - se_weight, 0), ymax = mean_weight + se_weight),
    width = 0.08
  ) +
  geom_text(
    aes(y = mean_weight + se_weight + 4, label = significance_group),
    fontface = "bold",
    size = 4.5,
    show.legend = FALSE
  ) +
  scale_y_continuous(expand = expansion(mult = c(0.03, 0.15))) +
  labs(
    title = "Site x treatment response in seed production",
    subtitle = paste0(
      "Mean total harvested seed weight per plot +/- SE; letters compare treatments within site; treatment x site ",
      format_p(weight_interaction_p)
    ),
    x = "Treatment",
    y = "Seed weight per plot (g)",
    shape = "Site",
    linetype = "Site"
  ) +
  theme_classic(base_size = 13) +
  theme(legend.position = "top", axis.title = element_text(face = "bold"))
print(fig_interaction)
save_plot("Supplement_Figure_2_Site_Treatment_Interaction.png", fig_interaction, width = 8.5, height = 6)

# ---- QUANTITY vs FILL CORRELATION ----
cor_dat <- seed_dat |>
  filter(!is.na(weight_g), !is.na(filled_proportion))
weight_quality_cor <- cor.test(cor_dat$weight_g, cor_dat$filled_proportion,
                              method = "spearman", exact = FALSE)
correlation_results <- tibble(
  statistic = unname(weight_quality_cor$statistic),
  rho = unname(weight_quality_cor$estimate),
  p_value = weight_quality_cor$p.value,
  n = nrow(cor_dat)
)
print(correlation_results)
write_table(correlation_results, "23_weight_vs_filled_seed_correlation.csv")

fig_quality_weight <- cor_dat |>
  ggplot(aes(x = weight_g, y = filled_proportion, shape = site)) +
  geom_point(alpha = 0.65, size = 2.2) +
  geom_smooth(method = "lm", se = TRUE) +
  facet_wrap(~ treatment) +
  scale_y_continuous(labels = percent_format()) +
  labs(
    title = "Relationship between seed quantity and seed filling",
    subtitle = "Seed-producing Liatris plugs only",
    x = "Harvested seed mass per plug (g)", y = "Filled-seed proportion", shape = "Site"
  ) +
  theme_classic(base_size = 13) +
  theme(legend.position = "top", strip.background = element_blank(), strip.text = element_text(face = "bold"))
print(fig_quality_weight)
save_plot("Supplement_Figure_3_Weight_vs_Filled_Seed.png", fig_quality_weight, width = 9, height = 5.5)

# ---- ESTIMATED FILLED-SEED MASS (DESCRIPTIVE PROXY, NOT VIABILITY) ----
dat <- dat |>
  mutate(estimated_filled_weight_g = case_when(
    produced_seed == 0 ~ 0,
    !is.na(weight_g) & !is.na(filled_proportion) ~ weight_g * filled_proportion,
    TRUE ~ NA_real_
  ))

estimated_filled_summary <- dat |>
  group_by(site, treatment) |>
  summarise(
    estimated_filled_weight_g = if_else(any(is.na(estimated_filled_weight_g)),
                                        NA_real_, sum(estimated_filled_weight_g)),
    estimated_filled_weight_per_planted_g = if_else(
      any(is.na(estimated_filled_weight_g)), NA_real_,
      mean(estimated_filled_weight_g)),
    n_missing_filled_weight = sum(is.na(estimated_filled_weight_g))
  )
write_table(estimated_filled_summary, "24_estimated_filled_seed_mass.csv")


# ---- COMBINED FIGURE ----
combined_figure <- (fig_weight + fig_mortality + fig_filled + plot_layout(guides = "collect")) &
  theme(legend.position = "bottom")
save_plot("Figure_4_Combined_Performance.png", combined_figure, width = 15, height = 6)

# Confirm all manuscript/supplementary figure exports, and create a searchable
# manifest BEFORE any later Excel/model exports can fail.
expected_figures <- c(
  "Figure_1_Total_Seed_Weight.png",
  "Figure_2_Mortality.png",
  "Figure_3_Filled_Seed_Proportion.png",
  "Figure_4_Combined_Performance.png",
  "Supplement_Figure_Seed_Production_Success.png",
  "Supplement_Figure_1_Weight_Among_Producers.png",
  "Supplement_Figure_2_Site_Treatment_Interaction.png",
  "Supplement_Figure_3_Weight_vs_Filled_Seed.png"
)
figure_manifest <- tibble(
  figure = expected_figures,
  absolute_path = file.path(figures_dir, expected_figures),
  exists = file.exists(file.path(figures_dir, expected_figures)),
  size_bytes = as.numeric(file.info(file.path(figures_dir, expected_figures))$size)
)
write_table(figure_manifest, "99_figure_manifest.csv")
if (!all(figure_manifest$exists & !is.na(figure_manifest$size_bytes) &
         figure_manifest$size_bytes > 0)) {
  stop("Some figures failed to save. See: ", file.path(tables_dir, "99_figure_manifest.csv"))
}
message("SUCCESS: All ", nrow(figure_manifest),
        " manuscript/supplementary figures saved to: ", figures_dir)

# ---- SAVE MODELS ----
saveRDS(mortality_model_interaction, file.path(models_dir, "mortality_model.rds"))
saveRDS(produced_seed_model_interaction, file.path(models_dir, "produced_seed_model.rds"))
saveRDS(quality_model_interaction, file.path(models_dir, "filled_seed_model.rds"))
saveRDS(weight_model_interaction, file.path(models_dir, "plot_weight_model.rds"))
saveRDS(positive_weight_model_interaction, file.path(models_dir, "positive_weight_model.rds"))

# ---- CONSOLIDATED EXCEL ----
results_excel <- list(
  Validation = validation,
  Data_Issues = data_issues,
  Design = design_check,
  Block_Completeness = block_check,
  Plot_Validation = plot_check,
  Initial_Overdispersion = initial_dispersion,
  GLMM_Fit_Status = secondary_model_status,
  Positive_Mass_Fit_Status = positive_weight_fit_status,
  DHARMa_Alerts = residual_alerts,
  Mortality_DHARMa = mortality_diagnostics,
  Seed_Success_DHARMa = produced_seed_diagnostics,
  Seed_Fill_DHARMa = quality_diagnostics,
  Positive_Mass_DHARMa = positive_weight_diagnostics,
  Mortality_Interaction_LRT = as.data.frame(mortality_interaction_test),
  Seed_Success_Interaction_LRT = as.data.frame(produced_seed_interaction_test),
  Seed_Fill_Interaction_LRT = as.data.frame(quality_interaction_test),
  Seed_Mass_Interaction_F = as.data.frame(weight_interaction_test),
  Treatment_Summary = plant_treatment_summary,
  Site_Treatment_Summary = site_treatment_summary,
  Plot_Data = plot_dat,
  Plot_Summary = plot_summary,
  Mortality_EMM_Site = as.data.frame(mortality_emm_site),
  Mortality_Pairs_Site = as.data.frame(mortality_pairs_site),
  Mortality_Diagnostics = mortality_diagnostics,
  Mortality_Treatment_Test = as.data.frame(mortality_treatment_test),
  Produced_Seed_EMM_Site = as.data.frame(produced_seed_emm_site),
  Produced_Seed_Pairs_Site = as.data.frame(produced_seed_pairs_site),
  Filled_EMM_Site = as.data.frame(quality_emm_site),
  Filled_Pairs_Site = as.data.frame(quality_pairs_site),
  Filled_Diagnostics = quality_diagnostics,
  Weight_EMM_Site = as.data.frame(weight_emm_site),
  Weight_Pairs_Site = as.data.frame(weight_pairs_site),
  Positive_Weight_EMM = as.data.frame(positive_weight_emm),
  Positive_Weight_Pairs = as.data.frame(positive_weight_pairs),
  Weight_Quality_Correlation = correlation_results,
  Estimated_Filled_Weight = estimated_filled_summary,
  Weight_Significance_Letters = weight_letters,
  Mortality_Significance_Letters = mortality_letters,
  Seed_Success_Letters = produced_seed_letters,
  Filled_Significance_Letters = quality_letters,
  Positive_Mass_Letters = positive_weight_letters
)
write_xlsx(results_excel, file.path(output_dir, "Companion_Project_Final_Analysis.xlsx"))

capture.output(sessionInfo(), file = file.path(output_dir, "session_info.txt"))

cat("\n\nANALYSIS COMPLETE\n")
cat("Plug records analyzed: ", nrow(dat), "\n")
cat("Plots analyzed: ", nrow(plot_dat), "\n")
cat("Surviving plugs: ", sum(dat$survived), "\n")
cat("Mortality: ", sum(dat$mortality), "\n")
cat("Plugs recorded as producing seed: ", sum(dat$produced_seed), "\n")
cat("Surviving plugs that did not produce seed: ", sum(dat$survived == 1 & dat$produced_seed == 0), "\n")
cat("Seed-producing plugs missing seed mass: ", sum(dat$produced_seed == 1 & is.na(dat$weight_g)), "\n")
# Snapshot the fitted random-effects variance components for interpretation.
capture.output(VarCorr(mortality_model_interaction),
               file = file.path(models_dir, "mortality_random_variances.txt"))
capture.output(VarCorr(produced_seed_model_interaction),
               file = file.path(models_dir, "produced_seed_random_variances.txt"))
capture.output(VarCorr(quality_model_interaction),
               file = file.path(models_dir, "filled_seed_random_variances.txt"))
capture.output(VarCorr(positive_weight_model_interaction),
               file = file.path(models_dir, "conditional_seed_mass_random_variances.txt"))
cat("\nRead output/tables/09_model_fit_status.csv and diagnostics_* files ",
    "before reporting p-values.\n")

message("FIGURES DIRECTORY: ", normalizePath(figures_dir, winslash = "/"))
