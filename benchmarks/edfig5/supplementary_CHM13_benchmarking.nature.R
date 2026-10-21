#!/usr/bin/env Rscript
# ===========================================================================
# Supplementary Figure -- CHM13 self-benchmark
#
# Validates that the KaryoScope database accurately reproduces the CHM13
# annotations it was built from.  Single haplotype (CHM13 is haploid).
#
# Panel A: RepeatMasker concordance (stacked bars)
# Panel B: CenSat concordance (stacked bars)
# Panel C: RepeatMasker confusion matrix
# Panel D: CenSat confusion matrix
#
# Nature version (Extended Data Fig. 5) of KaryoScope-archive
# workflow/scripts/supplementary_CHM13_benchmarking_old.R, the version that drew the published
# figure, with the cluster's Nature size constants and the Quartz PDF device.
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
  make_option("--rm", type = "character",
              help = "RepeatMasker contrast TSV"),
  make_option("--cs", type = "character",
              help = "CenSat contrast TSV"),
  make_option("--hierarchy", type = "character",
              help = "KaryoScope hierarchy file"),
  make_option("--outdir", type = "character", default = "main",
              help = "Output directory [default: %default]")
)

opt <- parse_args(OptionParser(option_list = option_list))

required <- c("rm", "cs", "hierarchy")
missing <- required[!required %in% names(opt) | sapply(opt[required], is.null)]
if (length(missing) > 0) {
  stop("Missing required arguments: --", paste(missing, collapse = ", --"))
}

OUT_DIR <- opt[["outdir"]]
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)


# -- Theme and palettes -------------------------------------------------------

theme_natbiotech <- function(base_size = 7) {
  theme_classic(base_size = base_size, base_family = "Arial") +
    theme(
      axis.line = element_line(linewidth = 0.3, colour = "black"),
      axis.ticks = element_line(linewidth = 0.3, colour = "black"),
      axis.text = element_text(size = rel(0.85), colour = "black"),
      axis.title = element_text(size = rel(1), face = "plain"),
      legend.title = element_text(size = rel(0.9), face = "bold"),
      legend.text = element_text(size = rel(0.85)),
      legend.key.size = unit(3, "mm"),
      legend.background = element_blank(),
      panel.grid = element_blank(),
      plot.title = element_text(size = rel(1), face = "bold", hjust = 0),
      plot.tag = element_text(size = rel(8 / 7), face = "bold"),
      plot.margin = margin(2, 4, 2, 2, unit = "mm")
    )
}

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


# -- Label orders -------------------------------------------------------------

RM_TRUTH_ORDER <- c(
  "LINE", "SINE", "LTR", "DNA", "tandem_repeat", "RC", "Retroposon",
  "rRNA", "snRNA", "scRNA", "srpRNA", "tRNA",
  "Unknown", "nonrepeat"
)

RM_PREDICTION_ORDER <- c(
  "repeat", "Interspersed_Repeat",
  "RNA", "rRNA", "scRNA", "snRNA", "srpRNA", "tRNA",
  "Transposable_Element", "Class_I_Retrotransposition",
  "LINE", "LINE-dependent_Retroposon", "Retroposon", "SINE", "LTR",
  "Class_II_DNA_Transposition", "DNA", "RC",
  "Unknown", "tandem_repeat",
  "nonrepeat", "novel"
)

CENSAT_LEAF_ORDER <- c(
  "active_hor", "dhor", "hor", "mixedAlpha", "mon",
  "bSat", "cenSat", "ct", "gSat",
  "HSat1A", "HSat1B", "HSat2", "HSat3",
  "rDNA", "p_arm", "q_arm"
)

CENSAT_PREDICTION_ORDER <- c(
  "centromeric", "aSat", "alpha_hor",
  "active_hor", "dhor", "hor", "mixedAlpha", "mon",
  "bSat", "cenSat", "ct", "gSat",
  "HSat", "HSat1", "HSat1A", "HSat1B", "HSat2", "HSat3",
  "rDNA", "arm", "p_arm", "q_arm", "novel"
)


# -- Load data ---------------------------------------------------------------

message("Loading data...")

read_contrast <- function(filepath) {
  read_tsv(filepath, col_types = "cccdc", show_col_types = FALSE) %>%
    rename(bp = length)
}

rm_data <- read_contrast(opt[["rm"]])
cs_data <- read_contrast(opt[["cs"]])

message(sprintf("  RepeatMasker: %s rows, %.1f Gb",
                format(nrow(rm_data), big.mark = ","), sum(rm_data$bp) / 1e9))
message(sprintf("  CenSat:       %s rows, %.1f Gb",
                format(nrow(cs_data), big.mark = ","), sum(cs_data$bp) / 1e9))


# -- Panel A: RM stacked bars ------------------------------------------------

build_panel_A <- function(data, stage_filter = "smoothed") {
  plot_data <- data %>%
    filter(stage == stage_filter) %>%
    mutate(compat_display = classify_rm_compat(compatibility)) %>%
    group_by(truth, compat_display) %>%
    summarise(bp = sum(bp), .groups = "drop") %>%
    group_by(truth) %>%
    mutate(pct = bp / sum(bp)) %>%
    ungroup() %>%
    filter(truth %in% RM_TRUTH_ORDER) %>%
    mutate(
      truth = factor(truth, levels = rev(RM_TRUTH_ORDER)),
      compat_display = factor(compat_display, levels = names(rm_compat_colours))
    )

  ggplot(plot_data, aes(x = pct, y = truth, fill = compat_display)) +
    geom_col(position = "stack", width = 0.7) +
    scale_x_continuous(labels = percent_format(), expand = expansion(mult = c(0, 0.02))) +
    scale_fill_manual(values = rm_compat_colours, name = "Concordance") +
    labs(x = "Fraction of annotated bases", y = NULL, title = "RepeatMasker") +
    theme_natbiotech() +
    theme(legend.position = "bottom")
}


# -- Panel B: CenSat stacked bars --------------------------------------------

build_panel_B <- function(data, stage_filter = "smoothed") {
  plot_data <- data %>%
    filter(stage == stage_filter) %>%
    mutate(compat_display = classify_cs_compat(compatibility)) %>%
    group_by(truth, compat_display) %>%
    summarise(bp = sum(bp), .groups = "drop") %>%
    group_by(truth) %>%
    mutate(pct = bp / sum(bp)) %>%
    ungroup() %>%
    filter(truth %in% CENSAT_LEAF_ORDER) %>%
    mutate(
      truth = factor(truth, levels = rev(CENSAT_LEAF_ORDER)),
      compat_display = factor(compat_display, levels = names(censat_compat_colours))
    )

  ggplot(plot_data, aes(x = pct, y = truth, fill = compat_display)) +
    geom_col(position = "stack", width = 0.7, show.legend = TRUE) +
    scale_x_continuous(labels = percent_format(), expand = expansion(mult = c(0, 0.02))) +
    scale_fill_manual(values = censat_compat_colours, name = "Compatibility", drop = FALSE) +
    labs(x = "Fraction of annotated bases", y = NULL, title = "CenSat") +
    theme_natbiotech() +
    theme(legend.position = "bottom", legend.title.position = "top") +
    guides(fill = guide_legend(nrow = 3))
}


# -- Panel C/D: Confusion matrices -------------------------------------------

build_confusion_matrix <- function(data, truth_order, prediction_order,
                                   stage_filter = "smoothed", title = "") {
  plot_data <- data %>%
    filter(stage == stage_filter) %>%
    group_by(truth, prediction) %>%
    summarise(bp = sum(bp), .groups = "drop")

  truth_labels <- intersect(truth_order, unique(plot_data$truth))
  pred_labels <- intersect(prediction_order, unique(plot_data$prediction))

  plot_data <- plot_data %>%
    filter(truth %in% truth_labels, prediction %in% pred_labels) %>%
    group_by(truth) %>%
    mutate(pct = bp / sum(bp)) %>%
    ungroup() %>%
    mutate(
      truth = factor(truth, levels = rev(truth_labels)),
      prediction = factor(prediction, levels = pred_labels)
    )

  ggplot(plot_data, aes(x = prediction, y = truth, fill = pct)) +
    geom_tile(colour = "white", linewidth = 0.3) +
    geom_text(
      aes(label = ifelse(pct >= 0.005, percent(pct, accuracy = 1, suffix = ""), "")),
      size = 1.8, colour = "black"
    ) +
    scale_fill_gradient(low = "white", high = "#2166AC",
                        labels = percent_format(suffix = ""), name = "% of truth",
                        limits = c(0, 1)) +
    labs(x = "KaryoScope prediction", y = "Truth annotation", title = title) +
    theme_natbiotech() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = rel(0.85)),
      axis.text.y = element_text(size = rel(0.85)),
      panel.border = element_rect(fill = NA, colour = "grey70", linewidth = 0.3),
      legend.position = "bottom",
      legend.key.width = unit(8, "mm"), legend.key.height = unit(2, "mm"),
      plot.title = element_text(size = rel(1), hjust = 0.5)
    )
}


# -- Assemble -----------------------------------------------------------------

message("Building panels...")

pA <- build_panel_A(rm_data, "smoothed")
pB <- build_panel_B(cs_data, "smoothed")
pC <- build_confusion_matrix(rm_data, RM_TRUTH_ORDER, RM_PREDICTION_ORDER,
                              stage_filter = "smoothed", title = "RepeatMasker")
pD <- build_confusion_matrix(cs_data, CENSAT_LEAF_ORDER, CENSAT_PREDICTION_ORDER,
                              stage_filter = "smoothed", title = "CenSat")

figure <- (
  (pA / pB) | (pC / pD)
) +
  plot_annotation(
    tag_levels = "a",
    theme = theme(
      plot.tag = element_text(size = 8, face = "bold", family = "Arial")
    )
  ) +
  plot_layout(widths = c(1, 1.3))


# -- Save ---------------------------------------------------------------------

FIG_HEIGHT_MM <- as.numeric(Sys.getenv("FIG_HEIGHT_MM", "160.25"))

message("Saving figures...")

# cairo_pdf rounds glyph advances to whole points, which spaces 5-7 pt labels unevenly ("S INE").
# The macOS Quartz device places glyphs at the font's true widths and embeds Arial.
pdf_device <- if (capabilities("aqua")) {
  function(filename, ...) grDevices::quartz(file = filename, type = "pdf", ...)
} else {
  warning("Quartz unavailable: falling back to cairo_pdf, whose glyph spacing is uneven at small sizes")
  cairo_pdf
}

ggsave(file.path(OUT_DIR, "supplementary_CHM13_benchmarking.pdf"),
       figure, width = 183, height = FIG_HEIGHT_MM, units = "mm", device = pdf_device)

ggsave(file.path(OUT_DIR, "supplementary_CHM13_benchmarking.png"),
       figure, width = 183, height = FIG_HEIGHT_MM, units = "mm", dpi = 600)

message("Done! Figures saved to: ", OUT_DIR)
