suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})
source("R/utils.R")
ensure_directories()

model_order <- c(
  "Training mean",
  "Feature elastic net",
  "Feature random forest",
  "Hybrid elastic net",
  "Identity ridge"
)
model_colors <- c(
  "Training mean" = "#777777",
  "Feature elastic net" = "#D55E00",
  "Feature random forest" = "#3B6FB6",
  "Hybrid elastic net" = "#8A5A44",
  "Identity ridge" = "#007C83"
)

theme_chapter <- function(base_size = 8) {
  theme_minimal(base_size = base_size, base_family = "Helvetica") +
    theme(
      text = element_text(color = "#000000"),
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_blank(),
      axis.title = element_text(color = "#000000", size = 8),
      axis.title.x = element_text(margin = margin(t = 7)),
      axis.title.y = element_text(margin = margin(r = 5)),
      axis.text = element_text(color = "#000000", size = 7.5),
      legend.title = element_text(color = "#000000", size = 8),
      legend.text = element_text(color = "#000000", size = 7.5),
      plot.margin = margin(7, 6, 7, 6),
      plot.tag = element_text(face = "bold", color = "#000000", size = 10.5),
      plot.title = element_blank()
    )
}

figure_tables <- getOption("section3.figure_tables", "outputs/tables")
performance_ci <- fread(file.path(figure_tables, "model_performance_pair_bootstrap_ci.csv"))[metric == "rmse"]
performance_ci[, model := factor(model, levels = model_order)]

performance_plot <- ggplot(
  performance_ci,
  aes(x = estimate, y = model, color = model)
) +
  geom_errorbar(
    aes(xmin = lower_95, xmax = upper_95),
    orientation = "y",
    width = 0.13,
    linewidth = 0.45
  ) +
  geom_point(size = 2.25) +
  geom_text(
    aes(label = sprintf("%.2f", estimate)),
    color = "#000000",
    hjust = 0.5,
    nudge_y = 0.24,
    size = 2.6
  ) +
  scale_color_manual(values = model_colors, guide = "none") +
  scale_x_continuous(limits = c(14.8, 24.3), breaks = seq(15, 23, 2)) +
  scale_y_discrete(expand = expansion(add = c(0.52, 0.78))) +
  labs(x = "Leave-pair-out RMSE", y = NULL, tag = "A") +
  theme_chapter()

identity_pairs <- fread(file.path(figure_tables, "chapter_pair_points.csv"))
assert_true(nrow(identity_pairs) == 60L, "Expected 60 top-decile pair candidates")

candidate_plot <- ggplot(
  identity_pairs,
  aes(x = pair_prediction, y = selected_observed_loewe, color = selected_css)
) +
  geom_abline(intercept = 0, slope = 1, color = "#777777", linetype = "dashed", linewidth = 0.45) +
  geom_point(size = 2.15, alpha = 0.88) +
  scale_color_viridis_c(option = "C", direction = -1, name = "Observed CSS") +
  labs(
    x = "Predicted pair score",
    y = "Observed mean\nLoewe score",
    tag = "B"
  ) +
  theme_chapter() +
  theme(legend.position = "right")

ggsave(
  "outputs/figures/model_performance_rmse.png",
  performance_plot + labs(tag = NULL),
  width = 6.5,
  height = 4.0,
  dpi = 320,
  bg = "white"
)
ggsave(
  "outputs/figures/model_performance_rmse.pdf",
  performance_plot + labs(tag = NULL),
  width = 6.5,
  height = 4.0,
  device = "pdf",
  useDingbats = FALSE,
  bg = "white"
)
ggsave(
  "outputs/figures/candidate_pair_confirmation.png",
  candidate_plot + labs(tag = NULL),
  width = 6.5,
  height = 4.6,
  dpi = 320,
  bg = "white"
)
ggsave(
  "outputs/figures/candidate_pair_confirmation.pdf",
  candidate_plot + labs(tag = NULL),
  width = 6.5,
  height = 4.6,
  device = "pdf",
  useDingbats = FALSE,
  bg = "white"
)

save_combined <- function(path, type = c("png", "pdf")) {
  type <- match.arg(type)
  if (type == "png") {
    png(path, width = 6.5, height = 2.65, units = "in", res = 320, bg = "white")
  } else {
    pdf(path, width = 6.5, height = 2.65, bg = "white", useDingbats = FALSE)
  }
  grid::grid.newpage()
  layout <- grid::grid.layout(1L, 2L, widths = grid::unit(c(1.06, 1.14), "null"))
  grid::pushViewport(grid::viewport(layout = layout))
  print(performance_plot, vp = grid::viewport(layout.pos.row = 1L, layout.pos.col = 1L))
  print(candidate_plot, vp = grid::viewport(layout.pos.row = 1L, layout.pos.col = 2L))
  dev.off()
}

save_combined("outputs/figures/worked_example_results.png", "png")
save_combined("outputs/figures/worked_example_results.pdf", "pdf")

if (!isTRUE(getOption("section3.figure_only", FALSE))) {
ranking <- fread("outputs/tables/triplet_top_decile_overall.csv")
ranking[, model := factor(model, levels = model_order[model_order %in% ranking$model])]
baseline_rate <- unique(round(ranking$baseline_positive_rate, 12))
assert_true(length(baseline_rate) == 1L, "Baseline positive rate is inconsistent")

ranking_plot <- ggplot(ranking, aes(x = positive_rate, y = model, color = model)) +
  geom_vline(xintercept = baseline_rate, color = "#777777", linetype = "dashed", linewidth = 0.6) +
  geom_segment(aes(x = baseline_rate, xend = positive_rate, yend = model), linewidth = 0.8) +
  geom_point(size = 3.0) +
  geom_text(
    aes(label = sprintf("%.1f%%", 100 * positive_rate)),
    color = "#202020",
    hjust = 0,
    nudge_x = 0.018,
    size = 3.2
  ) +
  scale_color_manual(values = model_colors, guide = "none") +
  scale_x_continuous(
    limits = c(0.30, 1.00),
    breaks = seq(0.4, 1.0, 0.1),
    labels = function(x) sprintf("%d%%", round(100 * x))
  ) +
  labs(x = "Observed Loewe score > 10 among top-decile predictions", y = NULL) +
  theme_chapter()
ggsave(
  "outputs/figures/top_decile_ranking.png",
  ranking_plot,
  width = 6.5,
  height = 3.5,
  dpi = 320,
  bg = "white"
)
ggsave(
  "outputs/figures/top_decile_ranking.pdf",
  ranking_plot,
  width = 6.5,
  height = 3.5,
  device = "pdf",
  useDingbats = FALSE,
  bg = "white"
)

outcomes <- readRDS("data/processed/oneil_analysis_outcomes.rds")
endpoint_plot <- ggplot(outcomes, aes(x = loewe, y = drugcomb_loewe)) +
  geom_hline(yintercept = 0, color = "#B0B0B0", linewidth = 0.35) +
  geom_vline(xintercept = 0, color = "#B0B0B0", linewidth = 0.35) +
  geom_point(color = "#2C6E8F", alpha = 0.13, size = 0.55) +
  geom_smooth(method = "lm", se = FALSE, color = "#D55E00", linewidth = 0.75) +
  labs(x = "Published DeepSynergy Loewe score", y = "DrugComb v1.5 recalculated Loewe score") +
  theme_chapter()
ggsave(
  "outputs/figures/loewe_endpoint_reconciliation.png",
  endpoint_plot,
  width = 6.0,
  height = 4.5,
  dpi = 320,
  bg = "white"
)
ggsave(
  "outputs/figures/loewe_endpoint_reconciliation.pdf",
  endpoint_plot,
  width = 6.0,
  height = 4.5,
  device = "pdf",
  useDingbats = FALSE,
  bg = "white"
)

predictions <- readRDS("data/processed/out_of_fold_predictions.rds")
calibration_models <- c(
  pred_identity_ridge = "Identity ridge",
  pred_feature_elastic_net = "Feature elastic net",
  pred_hybrid_elastic_net = "Hybrid elastic net",
  pred_feature_forest = "Feature random forest"
)
calibration <- rbindlist(lapply(names(calibration_models), function(model_column) {
  rbindlist(lapply(0:4, function(fold) {
    fold_data <- predictions[outer_fold == fold]
    fold_data[, prediction_decile := cut(
      get(model_column),
      breaks = quantile(get(model_column), probs = seq(0, 1, 0.1), na.rm = TRUE),
      include.lowest = TRUE,
      labels = FALSE
    )]
    fold_data[, .(
      mean_predicted = mean(get(model_column)),
      mean_observed = mean(observed_loewe),
      n = .N
    ), by = prediction_decile][, `:=`(
      outer_fold = fold,
      model = calibration_models[model_column]
    )]
  }))
}))
fwrite(calibration, "outputs/tables/calibration_deciles.csv")

calibration_plot <- ggplot(
  calibration,
  aes(x = mean_predicted, y = mean_observed, group = outer_fold, color = model)
) +
  geom_abline(intercept = 0, slope = 1, color = "#888888", linetype = "dashed", linewidth = 0.45) +
  geom_line(alpha = 0.38, linewidth = 0.55) +
  geom_point(alpha = 0.6, size = 1.1) +
  facet_wrap(~model, ncol = 2) +
  scale_color_manual(values = model_colors, guide = "none") +
  labs(x = "Mean predicted Loewe score by fold-specific decile", y = "Mean observed Loewe score") +
  theme_chapter(base_size = 9.5)
ggsave(
  "outputs/figures/calibration_deciles.png",
  calibration_plot,
  width = 7.5,
  height = 6.2,
  dpi = 320,
  bg = "white"
)

figure_inventory <- data.table(
  file = c(
    "worked_example_results.png",
    "model_performance_rmse.png",
    "candidate_pair_confirmation.png",
    "top_decile_ranking.png",
    "loewe_endpoint_reconciliation.png",
    "calibration_deciles.png"
  ),
  role = c(
    "Primary two-panel chapter candidate",
    "Standalone model comparison",
    "Standalone retrospective candidate confirmation",
    "Ranking enrichment diagnostic",
    "Endpoint-version diagnostic",
    "Calibration diagnostic"
  ),
  manuscript_status = c(
    "Candidate",
    "Supporting candidate",
    "Supporting candidate",
    "Analysis only unless space permits",
    "Analysis only unless discussed",
    "Analysis quality control"
  )
)
fwrite(figure_inventory, "outputs/tables/figure_inventory.csv")

writeLines(
  c(
    "Generated the primary two-panel results figure at 6.5 by 2.65 inches and five standalone diagnostics.",
    "The primary figure uses moderately enlarged, consistent pure-black text; Panel A estimates are placed above thinner confidence intervals.",
    "Confidence-interval details and the three-cell-line selection rule are reserved for the manuscript caption.",
    "Both 320 dpi PNG and vector PDF versions were saved for chapter-candidate figures."
  ),
  "logs/06_make_figures.log"
)
}
