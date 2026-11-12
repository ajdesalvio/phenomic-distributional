csv<-inputs[['contrasts']]
d <- read.csv(csv, check.names = FALSE)
needed <- c("VI", "Coefficient", "Reference", "Probability", "Estimate", "RawEstimate",
            "CR2_SE", "LowerPointwise95", "UpperPointwise95", "LowerFamily15_95",
            "UpperFamily15_95", "PointwiseCritical", "Family15Critical", "Status")
stopifnot(all(needed %in% names(d)), nrow(d) == 1500L,
          all(d$VI == "ExG2"), all(d$Reference == "5016_TM1"),
          all(d$Status == "APPROXIMATE_CR2_PLANT_WILD_BOOTSTRAP_T_FIXED_WEIGHTS_AND_SMOOTHERS"))
stopifnot(all(grepl("^Pedigree50(0[1-9]|1[0-5])$", d$Coefficient)))
d$Genotype <- as.integer(sub("^Pedigree", "", d$Coefficient))
d$Quantile_ID <- as.integer(round(100 * d$Probability))
stopifnot(setequal(unique(d$Genotype), 5001:5015),
          all(abs(d$Probability - d$Quantile_ID / 100) < 1e-12),
          !anyDuplicated(d[c("Genotype", "Quantile_ID")]))
for (g in 5001:5015) stopifnot(identical(sort(d$Quantile_ID[d$Genotype == g]), 1:100))
num <- c("Probability", "Estimate", "RawEstimate", "CR2_SE", "LowerPointwise95",
         "UpperPointwise95", "LowerFamily15_95", "UpperFamily15_95", "PointwiseCritical", "Family15Critical")
stopifnot(all(vapply(d[num], function(x) is.numeric(x) && all(is.finite(x)), logical(1))),
          all(d$CR2_SE > 0), all(d$PointwiseCritical > 0), all(d$Family15Critical > 0))
tol <- 1e-12
stopifnot(all(d$LowerFamily15_95 <= d$LowerPointwise95 + tol),
          all(d$LowerPointwise95 <= d$Estimate + tol),
          all(d$Estimate <= d$UpperPointwise95 + tol),
          all(d$UpperPointwise95 <= d$UpperFamily15_95 + tol))
formula_errors <- c(
  pointwise_lower = max(abs(d$LowerPointwise95 - (d$Estimate - d$PointwiseCritical * d$CR2_SE))),
  pointwise_upper = max(abs(d$UpperPointwise95 - (d$Estimate + d$PointwiseCritical * d$CR2_SE))),
  family_lower = max(abs(d$LowerFamily15_95 - (d$Estimate - d$Family15Critical * d$CR2_SE))),
  family_upper = max(abs(d$UpperFamily15_95 - (d$Estimate + d$Family15Critical * d$CR2_SE)))
)
stopifnot(max(formula_errors) < tol, length(unique(d$Family15Critical)) == 1L)
d <- d[order(d$Genotype, d$Quantile_ID), ]
row.names(d) <- NULL

# Exact saved bands and smoothed estimates: no new smoothing, interpolation,
# clipping, rescaling, multiple-testing calculation, or inferred intercept CI.
# Common scales permit direct comparison of contrast magnitudes across panels.
y_limits <- c(-0.04, 0.06)
stopifnot(min(d$LowerFamily15_95) >= y_limits[1], max(d$UpperFamily15_95) <= y_limits[2])
make_panel <- function(g) {
  z <- d[d$Genotype == g, ]
  ggplot2::ggplot(z, ggplot2::aes(x = Probability)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = LowerFamily15_95, ymax = UpperFamily15_95), fill = "#D6D6D6") +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = LowerPointwise95, ymax = UpperPointwise95), fill = "#8F8F8F") +
    ggplot2::geom_hline(yintercept = 0, colour = "#444444", linetype = "dashed", linewidth = 0.45) +
    ggplot2::geom_line(ggplot2::aes(y = Estimate), colour = "black", linewidth = 0.9) +
    ggplot2::scale_x_continuous(limits = c(0, 1), breaks = c(0, 0.25, 0.5, 0.75, 1),
                              labels = c("0", ".25", ".50", ".75", "1"), expand = ggplot2::expansion(mult = 0.02)) +
    ggplot2::scale_y_continuous(limits = y_limits, breaks = seq(-0.04, 0.06, by = 0.02),
                              labels = function(x) sprintf("%.2f", x), expand = ggplot2::expansion(mult = 0.025)) +
    ggplot2::labs(title = as.character(g), x = "Quantile level", y = bquote(beta[.(g)](p))) +
    ggplot2::theme_classic(base_size = 10, base_family = "sans") +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, face = "bold", size = 14, margin = ggplot2::margin(b = 3)),
                   axis.title = ggplot2::element_text(size = 10),
                   axis.title.x = ggplot2::element_text(margin = ggplot2::margin(t = 3)),
                   axis.title.y = ggplot2::element_text(margin = ggplot2::margin(r = 2)),
                   axis.text = ggplot2::element_text(size = 8, colour = "#333333"),
                   axis.line = ggplot2::element_line(linewidth = 0.4),
                   axis.ticks = ggplot2::element_line(linewidth = 0.4),
                   axis.ticks.length = grid::unit(0.06, "in"),
                   plot.margin = ggplot2::margin(5, 7, 4, 4))
}
legend_panel <- ggplot2::ggplot() +
  ggplot2::annotate("segment", x = 0.10, xend = 0.28, y = 0.69, yend = 0.69, linewidth = 0.9, colour = "black") +
  ggplot2::annotate("rect", xmin = 0.10, xmax = 0.28, ymin = 0.465, ymax = 0.555, fill = "#8F8F8F") +
  ggplot2::annotate("rect", xmin = 0.10, xmax = 0.28, ymin = 0.285, ymax = 0.375, fill = "#D6D6D6") +
  ggplot2::annotate("text", x = 0.34, y = c(0.69, 0.51, 0.33),
                    label = c("Estimate", "Pointwise 95%", "Joint 95%"), hjust = 0,
                    size = 3.5, family = "sans") +
  ggplot2::coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE, clip = "off") +
  ggplot2::theme_void() + ggplot2::theme(plot.margin = ggplot2::margin(5, 7, 4, 4))
panels <- c(lapply(5001:5015, make_panel), list(legend_panel))
figure <- gridExtra::arrangeGrob(grobs = panels, nrow = 4, ncol = 4)

if (!dir.create(out, recursive = TRUE, showWarnings = FALSE)) stop("Cannot create output: ", out)
out <- normalizePath(out, winslash = "/", mustWork = TRUE)
stem <- file.path(out, "Figure3")
ggplot2::ggsave(paste0(stem, ".png"), figure, width = 10, height = 7, units = "in", dpi = 600, bg = "white")
ggplot2::ggsave(paste0(stem, ".jpg"), figure, width = 10, height = 7, units = "in", dpi = 600, bg = "white", quality = 100)
ggplot2::ggsave(paste0(stem, ".pdf"), figure, width = 10, height = 7, units = "in", device = grDevices::cairo_pdf, bg = "white")
write.csv(d, file.path(out, "Figure3_Plotted_Data_V1.csv"), row.names = FALSE)
