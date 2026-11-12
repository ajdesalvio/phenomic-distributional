recalculate_within_date_h2 <-
function (components) 
{
    if (!requireNamespace("data.table", quietly = TRUE)) {
        stop("The data.table package is required; no packages are installed automatically.")
    }
    d <- data.table::as.data.table(data.table::copy(components))
    required <- c("grp", "vcov", "Vegetation.Index")
    if (!all(required %in% names(d)) || nrow(d) == 0L) {
        stop("Input must contain nonempty full variance components with grp, vcov, and Vegetation.Index.")
    }
    date_columns <- intersect(c("DAT", "DAP"), names(d))
    if (length(date_columns) != 1L) {
        stop("Exactly one date column (DAT or DAP) is required; pooled-date fits are unsupported.")
    }
    allowed_groups <- c("Pedigree", "Genotype", "Range", "Row", "Pass", "Replicate", "Residual")
    if (anyNA(d$grp) || any(!d$grp %in% allowed_groups)) {
        stop("Only separate within-date intercept components are supported; pooled or interaction groups are rejected.")
    }
    genotype_group <- intersect(c("Pedigree", "Genotype"), unique(d$grp))
    if (length(genotype_group) != 1L || !"Residual" %in% d$grp) {
        stop("Exactly one genotype label and full residual components are required.")
    }
    if (!is.numeric(d$vcov)) 
        stop("vcov must be numeric, unrounded variance components.")
    key_columns <- intersect(c("Crop", "crop", "Env", "Environment", "Vegetation.Index", "DAT", "DAP", 
        "Probability"), names(d))
    if (anyNA(d[, ..key_columns])) 
        stop("Fit keys must not be missing.")
    date_values <- d[[date_columns]]
    if (!is.numeric(date_values) || any(!is.finite(date_values))) {
        stop("Date keys must be finite numeric values.")
    }
    if ("Probability" %in% names(d) && (!is.numeric(d$Probability) || any(!is.finite(d$Probability)) || 
        any(d$Probability < 0 | d$Probability > 1))) {
        stop("Probability must be finite and in [0, 1].")
    }
    added_columns <- c("H2_Individual_Adjusted", "H2_Recalculation_Status", "H2_Definition")
    if (any(c(added_columns, ".audit_row_id") %in% names(d))) {
        stop("Input already contains recalculation columns; use the unchanged original input.")
    }
    if (anyDuplicated(d[, c(key_columns, "grp"), with = FALSE])) {
        stop("Duplicate fit/component keys found; resolve provenance before recalculation.")
    }
    counts <- d[, .(n_genotype = sum(grp == genotype_group), n_residual = sum(grp == "Residual")), by = key_columns]
    if (any(counts$n_genotype != 1L | counts$n_residual != 1L)) {
        stop("Every fit must have exactly one genotype component and one residual component.")
    }
    audit <- d[, {
        vg <- vcov[grp == genotype_group]
        ve <- vcov[grp == "Residual"]
        status <- if (any(!is.finite(vcov) | vcov < 0)) {
            "invalid_variance_component"
        }
        else if (vg + ve <= 0) {
            "zero_genotype_plus_residual"
        }
        else {
            "computed_from_saved_components"
        }
        value <- if (status == "computed_from_saved_components") 
            vg/(vg + ve)
        else NA_real_
        list(H2_Individual_Adjusted = value, H2_Recalculation_Status = status)
    }, by = key_columns]
    d[, `:=`(.audit_row_id, .I)]
    out <- merge(d, audit, by = key_columns, all.x = TRUE, sort = FALSE)
    data.table::setorder(out, .audit_row_id)
    out[, `:=`(.audit_row_id, NULL)]
    out[, `:=`(H2_Definition, "Vg/(Vg+Ve); individual unit within date; conditional on nuisance design terms")]
    data.table::setcolorder(out, c(names(components), added_columns))
    out[]
}
