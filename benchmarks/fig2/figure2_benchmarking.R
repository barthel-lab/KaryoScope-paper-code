#!/usr/bin/env Rscript
# ===========================================================================
# Figure 2 -- Benchmarking KaryoScope against RepeatMasker and CenSat
#
# Panel A: Per-class concordance with RepeatMasker (stacked bars)
# Panel B: Per-satellite concordance with CenSat (stacked bars)
# Panel C: RepeatMasker confusion matrix with hierarchy brackets
# Panel D: CenSat confusion matrix with hierarchy brackets
#
# Nature Biotechnology figure specs:
#   Single column: 88 mm  |  1.5 column: 120 mm  |  Full page: 180 mm
#   Font: 5-7 pt labels, max 8 pt; sans-serif preferred
#   Resolution: 300 DPI minimum (we use 600 for print)
#   Colour-blind friendly palette
# ===========================================================================

suppressPackageStartupMessages({
  .libPaths(c(Sys.getenv("R_LIBS_USER"), .libPaths()))
  library(tidyverse)
  library(patchwork)
  library(scales)
  library(optparse)
})

# -- Command-line arguments ---------------------------------------------------

option_list <- list(
  make_option("--rm_hap1", type = "character",
              help = "RepeatMasker contrast TSV for haplotype 1"),
  make_option("--rm_hap2", type = "character",
              help = "RepeatMasker contrast TSV for haplotype 2"),
  make_option("--cs_hap1", type = "character",
              help = "CenSat contrast TSV for haplotype 1"),
  make_option("--cs_hap2", type = "character",
              help = "CenSat contrast TSV for haplotype 2"),
  make_option("--hierarchy", type = "character",
              help = "KaryoScope hierarchy file (feature_set, child, parent)"),
  make_option("--outdir", type = "character", default = "main",
              help = "Output directory for figures [default: %default]")
)

opt <- parse_args(OptionParser(option_list = option_list))

required <- c("rm_hap1", "rm_hap2", "cs_hap1", "cs_hap2", "hierarchy")
missing <- required[!required %in% names(opt) | sapply(opt[required], is.null)]
if (length(missing) > 0) {
  stop("Missing required arguments: --", paste(missing, collapse = ", --"))
}

RM_CONTRAST_HAP1 <- opt[["rm_hap1"]]
RM_CONTRAST_HAP2 <- opt[["rm_hap2"]]
CS_CONTRAST_HAP1 <- opt[["cs_hap1"]]
CS_CONTRAST_HAP2 <- opt[["cs_hap2"]]
HIERARCHY_FILE   <- opt[["hierarchy"]]
OUT_DIR          <- opt[["outdir"]]

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)


# -- Nature Biotech theme ----------------------------------------------------

theme_natbiotech <- function(base_size = 7) {
  theme_classic(base_size = base_size, base_family = "sans") +
    theme(
      axis.line = element_line(linewidth = 0.3, colour = "black"),
      axis.ticks = element_line(linewidth = 0.3, colour = "black"),
      axis.text = element_text(size = rel(0.85), colour = "black"),
      axis.title = element_text(size = rel(1), face = "plain"),
      legend.title = element_text(size = rel(0.9), face = "bold"),
      legend.text = element_text(size = rel(0.85)),
      legend.key.size = unit(3, "mm"),
      legend.background = element_blank(),
      legend.position = "right",
      strip.background = element_blank(),
      strip.text = element_text(size = rel(1), face = "bold"),
      panel.grid = element_blank(),
      panel.spacing = unit(1, "mm"),
      plot.title = element_text(size = rel(1.1), face = "bold", hjust = 0),
      plot.tag = element_text(size = rel(1.3), face = "bold"),
      plot.margin = margin(2, 4, 2, 2, unit = "mm")
    )
}


# -- Colour palettes (colour-blind friendly) ---------------------------------

rm_compat_colours <- c(
  "Equivalent"  = "#4DAF4A",
  "Compatible"  = "#A6D854",
  "Novel"       = "#984EA3",
  "Mismatch"    = "#E41A1C"
)

censat_compat_colours <- c(
  "Equivalent"         = "#4DAF4A",
  "Compatible"         = "#A6D854",
  "Novel"              = "#984EA3",
  "Wrong centromeric"  = "#FF7F00",
  "Wrong arm"          = "#FFFF33",
  "Arm as centromeric" = "#E41A1C",
  "Centromeric as arm" = "#377EB8",
  "Other mismatch"     = "#999999"
)


# -- Parse hierarchy file ----------------------------------------------------

parse_hierarchy <- function(filepath) {
  df <- read_tsv(filepath, col_types = "ccc", show_col_types = FALSE)
  result <- list()
  for (fs in unique(df$feature_set)) {
    sub <- df %>% filter(feature_set == fs)
    pm <- setNames(sub$parent, sub$child)
    result[[fs]] <- pm
  }
  result
}


# -- Helper functions --------------------------------------------------------

read_contrast <- function(filepath, haplotype) {
  read_tsv(filepath, col_types = "cccdc", show_col_types = FALSE) %>%
    rename(bp = length) %>%
    mutate(haplotype = haplotype)
}

classify_rm_compat <- function(compatibility) {
  case_when(
    compatibility == "equivalent" ~ "Equivalent",
    compatibility == "compatible" ~ "Compatible",
    compatibility == "novel"      ~ "Novel",
    TRUE                          ~ "Mismatch"
  )
}

classify_cs_compat <- function(compatibility) {
  case_when(
    compatibility == "equivalent"         ~ "Equivalent",
    compatibility == "compatible"         ~ "Compatible",
    compatibility == "novel"              ~ "Novel",
    compatibility == "wrong_centromeric"  ~ "Wrong centromeric",
    compatibility == "wrong_arm"          ~ "Wrong arm",
    compatibility == "arm_as_centromeric" ~ "Arm as centromeric",
    compatibility == "centromeric_as_arm" ~ "Centromeric as arm",
    TRUE                                  ~ "Other mismatch"
  )
}


# -- Label orders (hierarchy DFS traversal) ----------------------------------

# Prediction axis: full DFS including internal nodes
RM_PREDICTION_ORDER <- c(
  "repeat", "Interspersed_Repeat",
  "RNA", "rRNA", "scRNA", "snRNA", "srpRNA", "tRNA",
  "Transposable_Element", "Class_I_Retrotransposition",
  "LINE", "LINE-dependent_Retroposon", "Retroposon", "SINE", "LTR",
  "Class_II_DNA_Transposition", "DNA", "RC",
  "Unknown",
  "tandem_repeat",
  "nonrepeat",
  "novel"
)

# Truth axis: leaf nodes only, DFS order
RM_TRUTH_ORDER <- c(
  "rRNA", "scRNA", "snRNA", "srpRNA", "tRNA",
  "LINE", "Retroposon", "SINE", "LTR",
  "DNA", "RC",
  "Unknown",
  "tandem_repeat",
  "nonrepeat"
)

# Panel A: truth order by genome abundance (for stacked bars)
RM_TRUTH_ORDER_ABUNDANCE <- c(
  "LINE", "SINE", "LTR", "DNA", "tandem_repeat", "RC", "Retroposon",
  "rRNA", "snRNA", "scRNA", "srpRNA", "tRNA",
  "Unknown",
  "nonrepeat"
)

# CenSat: leaf nodes
CENSAT_LEAF_ORDER <- c(
  "active_hor", "dhor", "hor", "mixedAlpha",
  "mon",
  "bSat", "cenSat", "ct", "gSat",
  "HSat1A", "HSat1B", "HSat2", "HSat3",
  "rDNA",
  "p_arm", "q_arm"
)

# CenSat: full DFS including internal nodes (for prediction axis)
CENSAT_PREDICTION_ORDER <- c(
  "centromeric", "aSat", "alpha_hor",
  "active_hor", "dhor", "hor", "mixedAlpha",
  "mon",
  "bSat", "cenSat", "ct", "gSat",
  "HSat", "HSat1", "HSat1A", "HSat1B", "HSat2", "HSat3",
  "rDNA",
  "arm", "p_arm", "q_arm",
  "novel"
)


# -- Bracket annotation builder ----------------------------------------------

compute_brackets <- function(visible_labels, parent_map, min_span = 2) {
  label_pos <- setNames(seq_along(visible_labels), visible_labels)

  get_depth <- function(node) {
    d <- 0
    current <- node
    while (current %in% names(parent_map) && parent_map[[current]] != "categorized") {
      current <- parent_map[[current]]
      d <- d + 1
    }
    d + 1
  }

  is_descendant_or_self <- function(node, ancestor) {
    if (node == ancestor) return(TRUE)
    current <- node
    while (current %in% names(parent_map)) {
      current <- parent_map[[current]]
      if (current == ancestor) return(TRUE)
    }
    FALSE
  }

  internal_nodes <- unique(unname(parent_map))
  internal_nodes <- setdiff(internal_nodes, "categorized")

  brackets <- list()
  for (inode in internal_nodes) {
    desc <- visible_labels[sapply(visible_labels, function(l) {
      is_descendant_or_self(l, inode)
    })]
    if (length(desc) < min_span) next
    positions <- label_pos[desc]
    depth <- get_depth(inode)
    brackets <- append(brackets, list(tibble(
      node = inode, xmin = min(positions), xmax = max(positions), depth = depth
    )))
  }

  if (length(brackets) == 0) {
    return(tibble(node = character(), xmin = numeric(), xmax = numeric(),
                  depth = integer(), y = numeric()))
  }

  bind_rows(brackets) %>%
    mutate(y = max(depth) - depth + 1)
}

build_bracket_plot <- function(bracket_data, n_labels, label_size = 1.1) {
  if (nrow(bracket_data) == 0) {
    return(ggplot() + theme_void())
  }

  max_y <- max(bracket_data$y)

  ggplot(bracket_data) +
    geom_segment(aes(x = xmin - 0.4, xend = xmax + 0.4, y = y, yend = y),
                 linewidth = 0.2, color = "grey40") +
    geom_segment(aes(x = xmin - 0.4, xend = xmin - 0.4, y = y, yend = y + 0.15),
                 linewidth = 0.2, color = "grey40") +
    geom_segment(aes(x = xmax + 0.4, xend = xmax + 0.4, y = y, yend = y + 0.15),
                 linewidth = 0.2, color = "grey40") +
    geom_text(aes(x = (xmin + xmax) / 2, y = y - 0.25, label = node),
              size = label_size, color = "grey30") +
    scale_x_continuous(limits = c(0.5, n_labels + 0.5), expand = c(0, 0)) +
    scale_y_continuous(limits = c(0.3, max_y + 0.5), expand = c(0, 0)) +
    labs(x = "KaryoScope prediction") +
    coord_cartesian(clip = "off") +
    theme_void() +
    theme(
      axis.title.x = element_text(size = 5, family = "sans", margin = margin(t = 2)),
      plot.margin = margin(0, 2, 2, 2, unit = "mm")
    )
}


# -- Load data ---------------------------------------------------------------

message("Loading data...")
hierarchy <- parse_hierarchy(HIERARCHY_FILE)
repeat_parents <- hierarchy[["repeat"]]
# Add tandem_repeat as a child of 'repeat' for bracket annotations
# (merged from Satellite, Simple_repeat, Low_complexity, Noninterspersed)
repeat_parents[["tandem_repeat"]] <- "repeat"
region_parents <- hierarchy[["region"]]

rm_data <- bind_rows(
  read_contrast(RM_CONTRAST_HAP1, "hap1"),
  read_contrast(RM_CONTRAST_HAP2, "hap2")
)
cs_data <- bind_rows(
  read_contrast(CS_CONTRAST_HAP1, "hap1"),
  read_contrast(CS_CONTRAST_HAP2, "hap2")
)

message(sprintf("  RepeatMasker: %s rows, %.1f Gb total",
                format(nrow(rm_data), big.mark = ","), sum(rm_data$bp) / 1e9))
message(sprintf("  CenSat:       %s rows, %.1f Gb total",
                format(nrow(cs_data), big.mark = ","), sum(cs_data$bp) / 1e9))


# ============================================================================
# Panel A -- RepeatMasker concordance (stacked bars)
# ============================================================================

build_panel_A <- function(rm_data, stage_filter = "smoothed") {
  plot_data <- rm_data %>%
    filter(stage == stage_filter) %>%
    mutate(compat_display = classify_rm_compat(compatibility)) %>%
    group_by(haplotype, truth, compat_display) %>%
    summarise(bp = sum(bp), .groups = "drop") %>%
    group_by(haplotype, truth) %>%
    mutate(total_bp = sum(bp), pct = bp / total_bp) %>%
    ungroup() %>%
    filter(truth %in% RM_TRUTH_ORDER_ABUNDANCE) %>%
    mutate(
      truth = factor(truth, levels = rev(RM_TRUTH_ORDER_ABUNDANCE)),
      compat_display = factor(compat_display, levels = names(rm_compat_colours))
    )

  ggplot(plot_data, aes(x = pct, y = truth, fill = compat_display)) +
    geom_col(position = "stack", width = 0.7) +
    facet_wrap(~ haplotype, ncol = 2, labeller = labeller(
      haplotype = c("hap1" = "Haplotype 1 (paternal)",
                     "hap2" = "Haplotype 2 (maternal)")
    )) +
    scale_x_continuous(labels = percent_format(), expand = expansion(mult = c(0, 0.02))) +
    scale_fill_manual(values = rm_compat_colours, name = "Concordance") +
    labs(x = "Fraction of annotated bases", y = NULL) +
    theme_natbiotech() +
    theme(legend.position = "bottom", axis.text.y = element_text(size = rel(0.9)))
}


# ============================================================================
# Panel B -- CenSat concordance (stacked bars)
# ============================================================================

build_panel_B <- function(cs_data, stage_filter = "smoothed") {
  plot_data <- cs_data %>%
    filter(stage == stage_filter) %>%
    mutate(compat_display = classify_cs_compat(compatibility)) %>%
    group_by(haplotype, truth, compat_display) %>%
    summarise(bp = sum(bp), .groups = "drop") %>%
    group_by(haplotype, truth) %>%
    mutate(total_bp = sum(bp), pct = bp / total_bp) %>%
    ungroup() %>%
    filter(truth %in% CENSAT_LEAF_ORDER) %>%
    mutate(
      truth = factor(truth, levels = rev(CENSAT_LEAF_ORDER)),
      compat_display = factor(compat_display, levels = names(censat_compat_colours))
    )

  ggplot(plot_data, aes(x = pct, y = truth, fill = compat_display)) +
    geom_col(position = "stack", width = 0.7) +
    facet_wrap(~ haplotype, ncol = 2, labeller = labeller(
      haplotype = c("hap1" = "Haplotype 1 (paternal)",
                     "hap2" = "Haplotype 2 (maternal)")
    )) +
    scale_x_continuous(labels = percent_format(), expand = expansion(mult = c(0, 0.02))) +
    scale_fill_manual(values = censat_compat_colours, name = "Compatibility", drop = FALSE) +
    labs(x = "Fraction of annotated bases", y = NULL) +
    theme_natbiotech() +
    theme(legend.position = "bottom", axis.text.y = element_text(size = rel(0.9))) +
    guides(fill = guide_legend(nrow = 2))
}


# ============================================================================
# Panel C / D -- Confusion matrix heatmap with hierarchy brackets
# ============================================================================

build_confusion_matrix <- function(data, truth_order, prediction_order,
                                   parent_map = NULL,
                                   stage_filter = "smoothed", title = "") {
  agg_data <- data %>%
    filter(stage == stage_filter) %>%
    group_by(truth, prediction) %>%
    summarise(bp = sum(bp), .groups = "drop")

  truth_labels <- intersect(truth_order, unique(agg_data$truth))
  pred_labels <- intersect(prediction_order, unique(agg_data$prediction))

  plot_data <- agg_data %>%
    filter(truth %in% truth_labels, prediction %in% pred_labels) %>%
    group_by(truth) %>%
    mutate(pct = bp / sum(bp)) %>%
    ungroup() %>%
    mutate(
      x_pos = as.numeric(factor(prediction, levels = pred_labels)),
      truth = factor(truth, levels = rev(truth_labels))
    )

  n_pred <- length(pred_labels)

  p_heat <- ggplot(plot_data, aes(x = x_pos, y = truth, fill = pct)) +
    geom_tile(colour = "white", linewidth = 0.3) +
    geom_text(
      aes(label = ifelse(pct >= 0.005, percent(pct, accuracy = 1), "")),
      size = 1.3, colour = "black"
    ) +
    scale_fill_gradient(
      low = "white", high = "#2166AC",
      labels = percent_format(),
      name = "Fraction\nof truth",
      limits = c(0, 1)
    ) +
    scale_x_continuous(
      breaks = seq_along(pred_labels),
      labels = pred_labels,
      limits = c(0.5, n_pred + 0.5),
      expand = c(0, 0)
    ) +
    labs(x = NULL, y = "Truth annotation", title = title) +
    theme_natbiotech() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = rel(0.7)),
      axis.text.y = element_text(size = rel(0.8)),
      axis.line.x = element_blank(),
      panel.border = element_rect(fill = NA, colour = "grey70", linewidth = 0.3),
      plot.title = element_text(size = rel(1), hjust = 0.5),
      plot.margin = margin(2, 2, 0, 2, unit = "mm")
    )

  if (!is.null(parent_map)) {
    bracket_data <- compute_brackets(pred_labels, parent_map)
    p_brackets <- build_bracket_plot(bracket_data, n_pred)

    combined <- p_heat / p_brackets +
      plot_layout(heights = c(5, 2))

    return(combined)
  }

  p_heat
}


# ============================================================================
# Assemble Figure 2
# ============================================================================

message("Building panels...")

pA <- build_panel_A(rm_data, "smoothed")
pB <- build_panel_B(cs_data, "smoothed")
pC <- build_confusion_matrix(rm_data, RM_TRUTH_ORDER, RM_PREDICTION_ORDER,
                              parent_map = repeat_parents,
                              stage_filter = "smoothed",
                              title = "RepeatMasker")
pD <- build_confusion_matrix(cs_data, CENSAT_LEAF_ORDER, CENSAT_PREDICTION_ORDER,
                              parent_map = region_parents,
                              stage_filter = "smoothed",
                              title = "CenSat")

figure2 <- (
  (pA / pB) | (wrap_elements(pC) / wrap_elements(pD))
) +
  plot_annotation(
    tag_levels = "A",
    theme = theme(plot.tag = element_text(size = 10, face = "bold", family = "sans"))
  ) +
  plot_layout(widths = c(1, 1.3))


# -- Save outputs ------------------------------------------------------------

message("Saving figures...")

ggsave(file.path(OUT_DIR, "figure2_benchmarking.pdf"),
       figure2, width = 210, height = 210, units = "mm", device = cairo_pdf)

ggsave(file.path(OUT_DIR, "figure2_benchmarking.png"),
       figure2, width = 210, height = 210, units = "mm", dpi = 600)

ggsave(file.path(OUT_DIR, "figure2A_repeatmasker_concordance.pdf"),
       pA, width = 120, height = 70, units = "mm", device = cairo_pdf)

ggsave(file.path(OUT_DIR, "figure2B_censat_concordance.pdf"),
       pB, width = 120, height = 80, units = "mm", device = cairo_pdf)

ggsave(file.path(OUT_DIR, "figure2C_confusion_repeatmasker.pdf"),
       pC, width = 150, height = 130, units = "mm", device = cairo_pdf)

ggsave(file.path(OUT_DIR, "figure2D_confusion_censat.pdf"),
       pD, width = 130, height = 130, units = "mm", device = cairo_pdf)

message("Done! Figures saved to: ", OUT_DIR)


# ============================================================================
# Supplementary: presmoothed vs smoothed comparison
# ============================================================================

build_smoothing_comparison <- function(data, classify_fn, compat_colours,
                                       panel_label = "") {
  plot_data <- data %>%
    mutate(compat_display = classify_fn(compatibility)) %>%
    group_by(haplotype, stage, compat_display) %>%
    summarise(bp = sum(bp), .groups = "drop") %>%
    mutate(
      Mb = bp / 1e6,
      stage = factor(stage, levels = c("presmoothed", "smoothed"),
                     labels = c("Pre-smoothed", "Smoothed"))
    )

  ggplot(plot_data, aes(x = stage, y = Mb, fill = compat_display)) +
    geom_col(position = "stack", width = 0.6) +
    facet_wrap(~ haplotype) +
    scale_fill_manual(values = compat_colours, name = "Concordance") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
    labs(x = NULL, y = "Megabases", title = panel_label) +
    theme_natbiotech() +
    theme(legend.position = "bottom")
}

p_smoothing_rm <- build_smoothing_comparison(
  rm_data, classify_rm_compat, rm_compat_colours, "RepeatMasker"
)
p_smoothing_cs <- build_smoothing_comparison(
  cs_data, classify_cs_compat, censat_compat_colours, "CenSat"
)

p_smoothing <- p_smoothing_rm / p_smoothing_cs +
  plot_annotation(tag_levels = "A")

ggsave(file.path(OUT_DIR, "supplementary_smoothing_comparison.pdf"),
       p_smoothing, width = 140, height = 120, units = "mm", device = cairo_pdf)

message("Supplementary smoothing figure saved.")
