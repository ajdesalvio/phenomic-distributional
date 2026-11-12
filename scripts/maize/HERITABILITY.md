# Maize heritability

Install `data.table`, `digest`, and `jsonlite`. Set the saved component and flight-time paths in `config/maize_heritability.json`; relative paths start from that file.

From the package root:

```sh
Rscript scripts/maize/recalculate_heritability.R config/maize_heritability.json output/maize_h2
```

Use a new output directory. This recalculates `Vg / (Vg + Ve)` from unrounded saved within-date variance components. It fits no models and uses no replication divisor. The historical `Heritability` column is preserved; updated calculations use `H2Individual` explicitly.

Outputs include Table S2, quantile/VI/date summaries, Figure 4 component data, complete HAM curves, and the flight-time calculations. The 11 manuscript environments and excluded dates are explicit in the configuration. Valid fitted quantiles remain in summaries even when another quantile is missing; HAM requires complete 100-quantile curves.

Canonical summaries retain ExG2 instead of GCC and NGRDI instead of NDI: 31 representations from 33 labels. Original labels remain in the saved-fit export. Formula and equivalence metadata are in `metadata/`.

Combined summaries pool fit-level estimates across environments. Genotype percentage in a stacked variance plot includes nuisance variance in its denominator and is distinct from heritability. Saved components do not establish historical convergence or model adequacy. See [verification](../../provenance/maize_ham_cleanup.md).
