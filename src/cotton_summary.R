cotton_summary <- function(d) {
q <- d[, .(MedianH2 = median(H2Individual), Fits = .N), by = .(Vegetation.Index, Probability)]
q[, `:=`(Rank, rank(-MedianH2, ties.method = "min")), by = Vegetation.Index]
pooled <- d[, .(MedianH2 = median(H2Individual), Fits = .N), by = Probability]
pooled[, `:=`(Rank, rank(-MedianH2, ties.method = "min"))]
five <- d[DAT %in% c(57, 86, 96, 113, 130), .(MedianH2 = median(H2Individual), Fits = .N), by = Probability]
five[, `:=`(Rank, rank(-MedianH2, ties.method = "min"))]
vi <- d[, .(MedianH2 = median(H2Individual), MinH2 = min(H2Individual), MaxH2 = max(H2Individual), Fits = .N, 
    Dates = uniqueN(DAT), SingularFits = sum(Singular), ZeroGenotypeFits = sum(Vg == 0), ZeroBlockFits = sum(Vb == 
        0)), by = Vegetation.Index]
main <- merge(data.table(Vegetation.Index = c("PSRI", "CVI", "ExG2"), Formula = c("(R - G)/RE", "NIR * R / G^2", 
    "(2*G - R - B)/(R + G + B)"), Reference = c("NV5 Geospatial Solutions (2025)", "Soft.Farm (2025)", 
    "Formula stated explicitly")), vi[, .(Vegetation.Index, MedianH2)], by = "Vegetation.Index", sort = FALSE)
main <- main[match(c("PSRI", "CVI", "ExG2"), Vegetation.Index)]
objects <- list(Table2_Cotton_Focal_VIs_V1 = main, Cotton_Focal_VI_Summary_V1 = vi, Cotton_Focal_Quantile_Summary_V1 = q, 
    Cotton_Combined_Quantile_Summary_V1 = pooled, Cotton_Five_Date_Quantile_Summary_V1 = five, Cotton_Focal_Date_Summary_V1 = d[, 
        .(Fits = .N, MedianH2 = median(H2Individual), SingularFits = sum(Singular)), by = .(Vegetation.Index, 
            DAT)])
objects
}
