# ========================================================================
# COMPANION PROJECT - UPDATED FINAL ANALYSIS (MORTALITY + SEED PRODUCTION)
# (Adapted: janitor::clean_names() and DHARMa replaced with equivalents
#  because those two packages are not reachable in this sandboxed
#  environment - no CRAN network access. Everything else is unchanged
#  from the user's original script.)
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

set.seed(20260814)

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

dir.create("output", showWarnings = FALSE)
dir.create("output/figures", showWarnings = FALSE, recursive = TRUE)
dir.create("output/tables", showWarnings = FALSE, recursive = TRUE)
dir.create("output/models", showWarnings = FALSE, recursive = TRUE)

save_plot <- function(filename, plot, width = 8, height = 6, dpi = 300) {
  ggsave(
    filename = file.path("output", "figures", filename),
    plot = plot, width = width, height = height, dpi = dpi, bg = "white"
  )
}

write_table <- function(x, filename) {
  write_csv(x, file.path("output", "tables", filename), na = "")
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

data_file <- "Data/Companion Project.xlsx"
message("Reading: ", data_file)

raw <- read_excel(data_file, sheet = "Complete Data") |> clean_names()

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
  select(site, plot, treatment, plant, mortality, produced_seed,
         seed_total, empty, filled, filled_proportion, weight_g)

write_table(data_issues, "00_data_issues.csv")
if (nrow(data_issues) > 0) {
  warning(
    "Data consistency issues detected. Review output/tables/00_data_issues.csv. ",
    "Plots containing a seed-producing plant with missing seed weight are ",
    "excluded from the primary seed-yield model until the missing weight is resolved."
  )
}

design_check <- dat |> distinct(site, plot, treatment) |>
  count(site, treatment, name = "n_plots")
print(design_check)
write_table(design_check, "00_design_check.csv")

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
    mean_weight_per_planted_g = mean(replace_na(yield_g, 0), na.rm = TRUE),
    median_weight_per_planted_g = median(replace_na(yield_g, 0), na.rm = TRUE),
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
    mean_weight_per_planted_g = mean(replace_na(yield_g, 0), na.rm = TRUE),
    mean_weight_producers_g = mean(weight_g[produced_seed == 1], na.rm = TRUE),
    mean_filled_proportion = mean(filled_proportion[produced_seed == 1], na.rm = TRUE)
  )

print(site_treatment_summary)
write_table(site_treatment_summary, "02_site_treatment_descriptive_statistics.csv")

plot_dat <- dat |>
  group_by(site, plot, plot_id, treatment) |>
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

# ---- MORTALITY MODEL ----
# Mortality is now analyzed directly from the Mortality field. It is no longer
# used as a proxy for whether a surviving plant produced seed.
mortality_model_site_only <- glmmTMB(mortality ~ site + (1 | plot_id), family = binomial, data = dat)
mortality_model_additive <- glmmTMB(mortality ~ treatment + site + (1 | plot_id), family = binomial, data = dat)
mortality_model_interaction <- glmmTMB(mortality ~ treatment * site + (1 | plot_id), family = binomial, data = dat)

# Overall treatment test (additive model vs site-only model) and interaction test.
mortality_treatment_test <- anova(mortality_model_site_only, mortality_model_additive)
mortality_interaction_test <- anova(mortality_model_additive, mortality_model_interaction)
print(mortality_treatment_test)
print(mortality_interaction_test)

mortality_emm_site <- emmeans(mortality_model_interaction, ~ treatment | site, type = "response")
print(mortality_emm_site)
mortality_pairs_site <- pairs(mortality_emm_site, adjust = "tukey")
print(mortality_pairs_site)
mortality_emm_overall <- emmeans(mortality_model_additive, ~ treatment, type = "response")
mortality_pairs_overall <- pairs(mortality_emm_overall, adjust = "tukey")

write_table(as.data.frame(mortality_emm_site), "05_mortality_emmeans_by_site.csv")
write_table(as.data.frame(mortality_pairs_site), "06_mortality_pairwise_by_site.csv")
write_table(as.data.frame(mortality_emm_overall), "07_mortality_emmeans_overall.csv")
write_table(as.data.frame(mortality_pairs_overall), "08_mortality_pairwise_overall.csv")
write_table(as.data.frame(mortality_treatment_test), "08b_mortality_overall_treatment_test.csv")

# ---- MORTALITY MODEL DIAGNOSTICS (performance package substitute for DHARMa) ----
tryCatch({
  png("output/figures/Supplement_Mortality_Model_Diagnostics.png", width = 2000, height = 1600, res = 250)
  print(check_model(mortality_model_interaction))
  dev.off()
}, error = function(e) message("Diagnostics plot skipped: ", conditionMessage(e)))

mortality_diagnostics <- tibble(
  test = c("Dispersion ratio", "Pearson chi-sq", "df"),
  value = {
    r <- residuals(mortality_model_interaction, type = "pearson")
    disp <- sum(r^2) / df.residual(mortality_model_interaction)
    c(disp, sum(r^2), df.residual(mortality_model_interaction))
  }
)
write_table(mortality_diagnostics, "09_mortality_model_diagnostics.csv")

# ---- SEED-PRODUCTION SUCCESS MODEL ----
# This is distinct from mortality: a surviving plant may still fail to produce seed.
produced_seed_model_additive <- glmmTMB(
  produced_seed ~ treatment + site + (1 | plot_id),
  family = binomial, data = dat
)
produced_seed_model_interaction <- glmmTMB(
  produced_seed ~ treatment * site + (1 | plot_id),
  family = binomial, data = dat
)
produced_seed_interaction_test <- anova(
  produced_seed_model_additive, produced_seed_model_interaction
)
print(produced_seed_interaction_test)

produced_seed_emm_site <- emmeans(
  produced_seed_model_interaction, ~ treatment | site, type = "response"
)
produced_seed_pairs_site <- pairs(produced_seed_emm_site, adjust = "tukey")
produced_seed_emm_overall <- emmeans(
  produced_seed_model_interaction, ~ treatment, type = "response"
)
produced_seed_pairs_overall <- pairs(produced_seed_emm_overall, adjust = "tukey")
print(produced_seed_emm_site)
print(produced_seed_pairs_site)

write_table(as.data.frame(produced_seed_emm_site), "09a_produced_seed_emmeans_by_site.csv")
write_table(as.data.frame(produced_seed_pairs_site), "09b_produced_seed_pairwise_by_site.csv")
write_table(as.data.frame(produced_seed_emm_overall), "09c_produced_seed_emmeans_overall.csv")
write_table(as.data.frame(produced_seed_pairs_overall), "09d_produced_seed_pairwise_overall.csv")

# ---- FILLED-SEED QUALITY MODEL ----
seed_dat <- dat |> filter(produced_seed == 1, !is.na(filled), !is.na(empty))

quality_model_additive <- glmmTMB(cbind(filled, empty) ~ treatment + site + (1 | plot_id),
                                  family = betabinomial(link = "logit"), data = seed_dat)
quality_model_interaction <- glmmTMB(cbind(filled, empty) ~ treatment * site + (1 | plot_id),
                                     family = betabinomial(link = "logit"), data = seed_dat)
quality_interaction_test <- anova(quality_model_additive, quality_model_interaction)
print(quality_interaction_test)

quality_emm_site <- emmeans(quality_model_interaction, ~ treatment | site, type = "response")
quality_pairs_site <- pairs(quality_emm_site, adjust = "tukey")
print(quality_emm_site)
print(quality_pairs_site)
quality_emm_overall <- emmeans(quality_model_interaction, ~ treatment, type = "response")
quality_pairs_overall <- pairs(quality_emm_overall, adjust = "tukey")

write_table(as.data.frame(quality_emm_site), "10_filled_seed_emmeans_by_site.csv")
write_table(as.data.frame(quality_pairs_site), "11_filled_seed_pairwise_by_site.csv")
write_table(as.data.frame(quality_emm_overall), "12_filled_seed_emmeans_overall.csv")
write_table(as.data.frame(quality_pairs_overall), "13_filled_seed_pairwise_overall.csv")

tryCatch({
  png("output/figures/Supplement_Filled_Seed_Model_Diagnostics.png", width = 2000, height = 1600, res = 250)
  print(check_model(quality_model_interaction))
  dev.off()
}, error = function(e) message("Diagnostics plot skipped: ", conditionMessage(e)))

quality_diagnostics <- tibble(
  test = c("Dispersion ratio", "Pearson chi-sq", "df"),
  value = {
    r <- residuals(quality_model_interaction, type = "pearson")
    disp <- sum(r^2) / df.residual(quality_model_interaction)
    c(disp, sum(r^2), df.residual(quality_model_interaction))
  }
)
write_table(quality_diagnostics, "14_filled_seed_model_diagnostics.csv")

# ---- PRIMARY SEED-YIELD MODEL ----
# Exclude plots with a seed-producing plant whose seed weight is missing;
# do not treat an unmeasured yield as a biological zero.
weight_model_dat <- plot_dat |> filter(weight_complete, !is.na(total_weight_g))

weight_model_additive <- lm(log1p(total_weight_g) ~ treatment + site, data = weight_model_dat)
weight_model_interaction <- lm(log1p(total_weight_g) ~ treatment * site, data = weight_model_dat)
weight_interaction_test <- anova(weight_model_additive, weight_model_interaction)
print(weight_interaction_test)
weight_anova <- anova(weight_model_interaction)
print(weight_anova)

weight_emm_site <- emmeans(weight_model_interaction, ~ treatment | site)
weight_pairs_site <- pairs(weight_emm_site, adjust = "tukey")
weight_emm_overall <- emmeans(weight_model_interaction, ~ treatment)
weight_pairs_overall <- pairs(weight_emm_overall, adjust = "tukey")
print(weight_emm_site)
print(weight_pairs_site)

write_table(as.data.frame(weight_emm_site), "15_weight_emmeans_by_site_log_scale.csv")
write_table(as.data.frame(weight_pairs_site), "16_weight_pairwise_by_site.csv")
write_table(as.data.frame(weight_emm_overall), "17_weight_emmeans_overall_log_scale.csv")
write_table(as.data.frame(weight_pairs_overall), "18_weight_pairwise_overall.csv")
write_table(broom::tidy(weight_model_interaction), "19_weight_model_coefficients.csv")

png("output/figures/Supplement_Weight_Model_Diagnostics.png", width = 2000, height = 1600, res = 250)
par(mfrow = c(2, 2))
plot(weight_model_interaction)
dev.off()

weight_shapiro <- shapiro.test(residuals(weight_model_interaction))
weight_diagnostics <- tibble(
  test = "Shapiro-Wilk residual normality",
  statistic = unname(weight_shapiro$statistic),
  p_value = weight_shapiro$p.value
)
write_table(weight_diagnostics, "20_weight_model_diagnostics.csv")

# ---- CONDITIONAL WEIGHT ANALYSIS ----
positive_dat <- dat |> filter(produced_seed == 1, !is.na(weight_g), weight_g > 0)
positive_weight_model_additive <- glmmTMB(log(weight_g) ~ treatment + site + (1 | plot_id),
                                          family = gaussian(), data = positive_dat)
positive_weight_model_interaction <- glmmTMB(log(weight_g) ~ treatment * site + (1 | plot_id),
                                             family = gaussian(), data = positive_dat)
positive_weight_interaction_test <- anova(positive_weight_model_additive, positive_weight_model_interaction)
print(positive_weight_interaction_test)

positive_weight_emm <- emmeans(positive_weight_model_interaction, ~ treatment | site)
positive_weight_pairs <- pairs(positive_weight_emm, adjust = "tukey")
print(positive_weight_emm)
print(positive_weight_pairs)
write_table(as.data.frame(positive_weight_emm), "21_positive_weight_emmeans.csv")
write_table(as.data.frame(positive_weight_pairs), "22_positive_weight_pairwise.csv")


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
    title = "Liatris mortality",
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
    x = "Treatment", y = "Plants producing seed (%)", fill = "Site"
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
    title = "Seed weight among reproductive Liatris plants",
    subtitle = "Different letters indicate Tukey-adjusted differences within site; nonreproductive plants excluded",
    x = "Treatment",
    y = "Harvested seed weight per plant (g)",
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

# ---- QUANTITY vs QUALITY CORRELATION ----
weight_quality_cor <- cor.test(seed_dat$weight_g, seed_dat$filled_proportion, method = "spearman", exact = FALSE)
correlation_results <- tibble(
  statistic = unname(weight_quality_cor$statistic),
  rho = unname(weight_quality_cor$estimate),
  p_value = weight_quality_cor$p.value,
  n = nrow(seed_dat)
)
print(correlation_results)
write_table(correlation_results, "23_weight_vs_filled_seed_correlation.csv")

fig_quality_weight <- seed_dat |>
  ggplot(aes(x = weight_g, y = filled_proportion, shape = site)) +
  geom_point(alpha = 0.65, size = 2.2) +
  geom_smooth(method = "lm", se = TRUE) +
  facet_wrap(~ treatment) +
  scale_y_continuous(labels = percent_format()) +
  labs(
    title = "Relationship between seed quantity and seed filling",
    subtitle = "Individual plants that produced harvestable seed",
    x = "Harvested seed weight per plant (g)", y = "Filled-seed proportion", shape = "Site"
  ) +
  theme_classic(base_size = 13) +
  theme(legend.position = "top", strip.background = element_blank(), strip.text = element_text(face = "bold"))
print(fig_quality_weight)
save_plot("Supplement_Figure_3_Weight_vs_Filled_Seed.png", fig_quality_weight, width = 9, height = 5.5)

# ---- ESTIMATED FILLED-SEED MASS ----
dat <- dat |>
  mutate(estimated_filled_weight_g = case_when(
    produced_seed == 0 ~ 0,
    !is.na(weight_g) & !is.na(filled_proportion) ~ weight_g * filled_proportion,
    TRUE ~ NA_real_
  ))

estimated_filled_summary <- dat |>
  group_by(site, treatment) |>
  summarise(
    estimated_filled_weight_g = sum(estimated_filled_weight_g, na.rm = TRUE),
    estimated_filled_weight_per_planted_g = mean(replace_na(estimated_filled_weight_g, 0))
  )
write_table(estimated_filled_summary, "24_estimated_filled_seed_mass.csv")


# ---- COMBINED FIGURE ----
combined_figure <- (fig_weight + fig_mortality + fig_filled + plot_layout(guides = "collect")) &
  theme(legend.position = "bottom")
save_plot("Figure_4_Combined_Performance.png", combined_figure, width = 15, height = 6)

# ---- SAVE MODELS ----
saveRDS(mortality_model_interaction, "output/models/mortality_model.rds")
saveRDS(produced_seed_model_interaction, "output/models/produced_seed_model.rds")
saveRDS(quality_model_interaction, "output/models/filled_seed_model.rds")
saveRDS(weight_model_interaction, "output/models/plot_weight_model.rds")
saveRDS(positive_weight_model_interaction, "output/models/positive_weight_model.rds")

# ---- CONSOLIDATED EXCEL ----
results_excel <- list(
  Validation = validation,
  Data_Issues = data_issues,
  Design = design_check,
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
  Produced_Seed_Significance_Letters = produced_seed_letters,
  Filled_Significance_Letters = quality_letters,
  Positive_Weight_Significance_Letters = positive_weight_letters
)
write_xlsx(results_excel, "output/Companion_Project_Final_Analysis.xlsx")

capture.output(sessionInfo(), file = "output/session_info.txt")

cat("\n\nANALYSIS COMPLETE\n")
cat("Plants analyzed: ", nrow(dat), "\n")
cat("Plots analyzed: ", nrow(plot_dat), "\n")
cat("Surviving plants: ", sum(dat$survived), "\n")
cat("Mortality: ", sum(dat$mortality), "\n")
cat("Plants recorded as producing seed: ", sum(dat$produced_seed), "\n")
cat("Survivors that did not produce seed: ", sum(dat$survived == 1 & dat$produced_seed == 0), "\n")
cat("Seed producers missing seed weight: ", sum(dat$produced_seed == 1 & is.na(dat$weight_g)), "\n")