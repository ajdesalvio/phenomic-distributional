# Heritability heatmaps

Install the Python packages in `environment/requirements-ham.txt`. The repository includes the small curve and date tables. Get the eleven images from the [data release](../../data/DOWNLOADS.md). If you store the images elsewhere, update the cotton and maize directories in `config/ham.json`, which also lists the filenames and checksums.

Run these commands from the repository root with a new output directory.

```sh
python scripts/figures/build_ham.py config/ham.json output/ham
Rscript scripts/figures/render_ham.R output/ham
```

The R renderer uses `png`, `jsonlite`, `digest`, and R’s built-in `grid`. The workflow writes eleven RGB and heatmap pairs, numeric arrays, and tables with curves and summaries. It exports Figure 5 as 600 dpi PNG and JPEG files and as a PDF.

For each image, the script calculates right-continuous empirical ranks from that image’s valid pixels. It assigns colors by interpolating the heritability curve for the corresponding date and vegetation index at those ranks. The workflow preserves the masks, endpoint handling, band order, rotations, and common jet color scale.

The supplied curves reproduce the five cotton PSRI and six maize ExGR maps. The source checksums and row counts are in `provenance/ham_input_excerpts.json`. If you use newly calculated curves, check their values before updating the paths and checksums. See the [validation notes](../../provenance/maize_ham_cleanup.md).
