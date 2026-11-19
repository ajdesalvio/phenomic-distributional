# Maize heritability

Install `data.table`, `digest`, and `jsonlite`. If you keep the saved variance components and flight times outside their default locations, edit their paths in `config/maize_heritability.json`. Relative paths start from that file’s directory.

From the package root:

```sh
Rscript scripts/maize/recalculate_heritability.R config/maize_heritability.json output/maize_h2
```

Use a new output directory. The script computes `Vg / (Vg + Ve)` directly from the unrounded variance components saved for each date. It puts the individual-unit estimates in `H2Individual` and keeps the original `Heritability` column.

Outputs include Table S2, quantile/VI/date summaries, Figure 4 component data, complete HAM curves, and the flight-time calculations. The 11 manuscript environments and excluded dates are explicit in the configuration. The summaries use all valid fitted quantiles. The heatmaps use only curves with estimates for all 100 quantiles.

The summaries use 31 distinct representations from the 33 original labels. They keep ExG2 for the ExG2/GCC pair and NGRDI for the NGRDI/NDI pair. Original labels remain in the saved-fit export. Formula and equivalence metadata are in `metadata/`.

Combined summaries pool fit-level estimates across environments. The genotype percentage in the stacked variance plots uses total variance, including nuisance components, as its denominator. Heritability uses only `Vg + Ve`. The audit verifies the calculations from saved components. Assessing convergence and model adequacy requires the underlying fits. See [verification](../../provenance/maize_ham_cleanup.md).
