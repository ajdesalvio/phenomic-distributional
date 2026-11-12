# Manuscript analysis code

Reproduce the revised manuscript from its frozen inputs, or run the separately documented model fits. Fresh ExGR fitting can produce different BLUEs as noted in [optimizer sensitivity](provenance/optimizer_sensitivity.md).

1. Follow [software setup](environment/README.md).
2. Choose a workflow below and obtain its inputs. [Configured inputs](data/configured_inputs.csv) lists paths and hashes across all configurations; [G2F imagery DOIs](data/G2F_COLLECTIONS.md) are supplied; [processed-data download and extraction instructions](data/DOWNLOADS.md) describe the prepared release (journal DOI pending). Edit paths in the relevant `config/*.json`. Relative paths start from that configuration file.
3. Run commands from this package's root. Quote paths containing spaces. Use new output locations unless a guide explicitly describes checkpoint resume.

| Task | Guide |
|---|---|
| Field join, frozen BLUE recovery, raw support, yield BLUEs | [Preparation](scripts/preprocess/README.md) |
| Optional fresh ExGR fits | [Refitting](scripts/preprocess/REFIT.md) |
| Within-environment nested CV and aggregation | [Nested CV](scripts/maize/README.md) |
| Saved across-environment predictions; optional FPCA/LOEO refits | [Across-environment prediction](scripts/maize/LOEO.md) |
| Cotton saved result export and optional model fits | [Cotton](scripts/cotton/README.md) |
| Maize saved component heritability and Table S2 | [Heritability](scripts/maize/HERITABILITY.md) |
| Figure 5 maps | [HAM](scripts/figures/HAM.md) |
| Remaining figures, correlations, and tables | [Figures and tables](scripts/figures/README.md) |

`scripts/` contains runnable commands; `src/` contains their shared functions. `config/` defines inputs and settings, `metadata/` holds curated tables, and `results/reference/` contains selected small reference results. `provenance/` records verification and limitations. Generated files belong in `output/`.

Saved result exports and renderers do not refit models. Nested CV, yield generation, and the explicitly marked fitting commands do. The preparation-only validator is optional when using other workflows:

```sh
Rscript scripts/preprocess/validate_inputs.R config/published_inputs.json output/input_check.json
```

Editable Figures 1–2 are in [figures/source](figures/source/README.md).

Historical scripts and figures are preserved separately in [legacy](legacy/README.md).

Original code is licensed under [MIT](LICENSE); original documentation and figures under [CC BY 4.0](LICENSE-DOCUMENTATION.md). See [license scope](LICENSING.md) for data and third-party exceptions.

For the original custom Python functions, see the [Python setup and usage instructions](legacy/pre_revision/README.md#installation), including shapefile creation, image cropping, soil masking, and vegetation index extraction.
