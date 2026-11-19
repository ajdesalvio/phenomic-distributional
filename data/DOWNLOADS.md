# Processed-data downloads

**Journal repository DOI pending.** The archives are ready for deposit. Download links will appear here once the DOI is available. For the original maize imagery, see the [G2F collections](G2F_COLLECTIONS.md).

## Choose archives

| Archive | Contents | Download size |
| --- | --- | ---: |
| `01_figures_tables_heritability.zip` | Figure/table inputs, maize heritability components, and HAM images | 192 MB |
| `02_maize_prediction.zip` | BLUEs, raw support, yield, saved LOEO models, and FPCA projections | 4.45 GB |
| `03_cotton_models.zip` | Cotton inputs, field map, and saved heritability/functional results | 573 MB |
| `04_raw_support_lineage.zip` | Raw ExGR quantiles, maize field map, and support provenance | 440 MB |
| `05_optional_reconstruction.zip` | Phenotype inputs, BLUE workspace, and original projection parts | 2.87 GB |
| `06_nested_cv_saved_evidence.zip` | Saved nested CV predictions, coefficients, fold assignments, and tuning results | 337 MB |

Downloading all six archives takes **8.86 GB**, and the extracted files take **31.78 GB**. For figures and Table S3, use archive 01. For cotton Table 2, use 03. Archive 02 supplies nested CV and the correlation calculations, and the preparation validator uses both 02 and 04. Use 05 and 06 for optional reconstruction or inspection.

## Download and extract

1. Download the chosen ZIPs plus `FILE_MANIFEST.csv`, `ARCHIVE_CHECKSUMS.txt`, and `verify_release.py` from the deposit into one download folder. The deposit also includes `README_RELEASE.txt` and `ARCHIVE_SUMMARY.csv`.
2. Check the ZIP’s checksum with Python’s standard library.

   ```sh
   python "/path/to/downloads/verify_release.py" --archives --group 01_figures_tables_heritability.zip
   ```

3. Open a terminal in the **repository root**, beside `README.md`, `config/`, and `scripts/`. Extract each ZIP there.

   ```sh
   python -m zipfile -e "/path/to/downloads/01_figures_tables_heritability.zip" .
   ```

   Substitute the filename of each archive you downloaded. Each ZIP already contains `data/`, so extracting at the root puts inputs where the configurations expect them. Use a fresh checkout or back up any locally modified inputs first.

4. Verify the extracted files.

   ```sh
   python "/path/to/downloads/verify_release.py" --installed . --group 01_figures_tables_heritability.zip
   ```

For `--archives`, omit `--group` to check all six downloads. For `--installed`, omit it once all six archives have been extracted. The manifest lists each file’s destination and checksum. Python’s ZIP reader supports these ZIP64 archives.

After extraction, choose a [workflow](../README.md). Archive 06 contains the saved nested CV outputs for inspection and for rebuilding feature summaries. To resume a fit, use checkpoints from your own run. We verified the release against commit `9a1d239bc5d530a123a81f368ed86151a326ccbe`. The deposit will provide the data license and citation.
