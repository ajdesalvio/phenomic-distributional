summarize_maize <- function(all,envs) {
d <- all[!Vegetation.Index %in% c("GCC","NDI")]
stopifnot(data.table::uniqueN(d$Vegetation.Index)==31L)
q <- rbindlist(lapply(c(envs, "Combined"), function(e) {
    z <- if (e == "Combined") 
        d
    else d[Env == e]
    z[, .(MedianH2 = median(H2Individual), Fits = .N), by = Probability][, `:=`(Env, e)][]
}))
q[, `:=`(Rank, rank(-MedianH2, ties.method = "min")), by = Env]
summary <- rbindlist(lapply(c(envs, "Combined"), function(e) {
    z <- if (e == "Combined") 
        d
    else d[Env == e]
    f <- q[Env == e & abs(Probability - 0.5) < 1e-10]
    data.table(Environment = e, N_Dates = uniqueN(paste(z$Env, z$DAP)), N_VIs = uniqueN(z$Vegetation.Index), 
        N_Fits = nrow(z), Median_H2_All_Quantiles = median(z$H2Individual), H2_Q50 = f$MedianH2, Rank_Q50 = f$Rank)
}))
eq <- rbindlist(lapply(list(c("ExG2", "GCC"), c("NGRDI", "NDI")), function(pair) {
    a <- merge(all[Vegetation.Index == pair[1]], all[Vegetation.Index == pair[2]], by = c("Env", "DAP", 
        "Probability"))
    data.table(Canonical = pair[1], Redundant = pair[2], MatchedFits = nrow(a), MedianAbsoluteDifference = median(abs(a$H2Individual.x - 
        a$H2Individual.y)), MaximumAbsoluteDifference = max(abs(a$H2Individual.x - a$H2Individual.y)), 
        N_DifferenceAbove1e_6 = sum(abs(a$H2Individual.x - a$H2Individual.y) > 1e-06))
}))
list(Table_S2_Maize_Individual_Unit_Heritability_V1=summary,
Maize_Canonical_31_Quantile_Summary_V1=q,
Maize_Canonical_31_VI_Summary_V1=d[,.(Fits=.N,MedianH2=median(H2Individual),MedianGenotypeFraction=median(PercentUnrounded)),by=.(Env,Vegetation.Index)],
Maize_Canonical_31_Date_Summary_V1=d[,.(Fits=.N,MedianH2=median(H2Individual)),by=.(Env,DAP)],
Maize_Affine_Saved_H2_Comparison_V1=eq)
}
