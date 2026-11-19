# Manuscript analysis code

Use the saved analysis inputs to reproduce the results in the maize and cotton manuscript. The workflow guides also explain how to fit new models. Refitting ExGR models can change the BLUEs. See the [optimizer sensitivity audit](provenance/optimizer_sensitivity.md).

1. Install the software listed in the [setup guide](environment/README.md).
2. Choose a workflow below and check which files it needs. The [processed-data guide](data/DOWNLOADS.md) explains how to extract the release once the journal DOI is available. The [input list](data/configured_inputs.csv) gives the expected paths and checksums. See the [G2F collections](data/G2F_COLLECTIONS.md) for source imagery.
3. Run commands from the repository root. Quote paths that contain spaces. Choose a new output location for each run unless you are resuming it as described in the guide. If you keep inputs elsewhere, edit the relevant file in `config/`. Relative input paths are resolved from that file’s directory.

| Task | Guide |
|---|---|
| Field join, frozen BLUE recovery, raw support, yield BLUEs | [Preparation](scripts/preprocess/README.md) |
| Optional fresh ExGR fits | [Refitting](scripts/preprocess/REFIT.md) |
| Within-environment nested CV and aggregation | [Nested CV](scripts/maize/README.md) |
| Across-environment predictions and optional FPCA or LOEO refits | [Across-environment prediction](scripts/maize/LOEO.md) |
| Export cotton results or fit new models | [Cotton](scripts/cotton/README.md) |
| Maize heritability from saved components and Table S2 | [Heritability](scripts/maize/HERITABILITY.md) |
| Figure 5 maps | [HAM](scripts/figures/HAM.md) |
| Remaining figures, correlations, and tables | [Figures and tables](scripts/figures/README.md) |

Run the entry points in `scripts/`. Their shared functions live in `src/`, and `config/` holds the input paths and settings. Find curated tables in `metadata/`, selected results in `results/reference/`, and validation records in `provenance/`. The examples write to `output/`.

Export and rendering commands use saved results. Nested CV, yield generation, and commands marked as fitting steps run models. To check the preparation inputs, run the optional validator below.

```sh
Rscript scripts/preprocess/validate_inputs.R config/published_inputs.json output/input_check.json
```

Editable Figures 1–2 are in [figures/source](figures/source/README.md).

Historical scripts and figures are preserved separately in [legacy](legacy/README.md).

The [MIT license](LICENSE) covers original code. [CC BY 4.0](LICENSE-DOCUMENTATION.md) covers original documentation and figures. See the [licensing guide](LICENSING.md) for data and third-party terms.

For the original custom Python functions, see the [Python setup and usage instructions](legacy/pre_revision/README.md#installation), including shapefile creation, image cropping, soil masking, and vegetation index extraction.
