# Heritability heatmaps

Install the Python packages in `environment/requirements-ham.txt`. Set the cotton and maize image directories in `config/ham.json`. The required eleven filenames and hashes are listed there; the small curve and date inputs are bundled.

From the package root:

```sh
python scripts/figures/build_ham.py config/ham.json output/ham
Rscript scripts/figures/render_ham.R output/ham
```

The R renderer requires `png`, `jsonlite`, and `digest` (`grid` ships with R). Use a new output directory. The commands create the eleven RGB/HAM pairs, numeric arrays, curve/summary tables, and the original Figure 5 layout as 600-dpi PNG/JPEG plus PDF.

Each image uses its own valid pixels to calculate right-continuous empirical ranks, then interpolates the matching date/VI heritability curve. Masks, endpoint handling, band order, rotations, and the common jet color scale are unchanged. These maps display quantile level heritability by pixel rank.

The bundled curves reproduce the manuscript's five cotton PSRI and six maize ExGR maps. Their source hashes and selection counts are in `provenance/ham_input_excerpts.json`. To use newly calculated curves, explicitly update the configured paths and hashes after checking them. See [verification](../../provenance/maize_ham_cleanup.md).
