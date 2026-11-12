# Processed-data downloads

**Journal repository DOI/download link: pending.** The release archives are prepared but are not yet publicly available. This page will link the deposit when its DOI is assigned. Original maize imagery is listed separately in [G2F collections](G2F_COLLECTIONS.md).

## Choose archives

| Archive | Contents | Download size |
| --- | --- | ---: |
| `01_figures_tables_heritability.zip` | Figure/table inputs, maize heritability components, and HAM images | 192 MB |
| `02_maize_prediction.zip` | BLUEs, raw support, yield, saved LOEO models, and FPCA projections | 4.45 GB |
| `03_cotton_models.zip` | Cotton inputs, field map, and saved heritability/functional results | 573 MB |
| `04_raw_support_lineage.zip` | Raw ExGR quantiles, maize field map, and support provenance | 440 MB |
| `05_optional_reconstruction.zip` | Phenotype inputs, BLUE workspace, and original projection parts | 2.87 GB |
| `06_nested_cv_saved_evidence.zip` | Saved V15 predictions, coefficients, folds, and tuning evidence | 337 MB |

All six total **8.86 GB downloaded / 31.78 GB extracted**. For figure rendering and S3, start with 01; cotton Table 2 uses 03. Nested CV and correlation calculations use 02. Preparation validation needs 02 and 04. Archives 05 and 06 support optional reconstruction and inspection.

## Download and extract

1. Download the chosen ZIPs plus `FILE_MANIFEST.csv`, `ARCHIVE_CHECKSUMS.txt`, and `verify_release.py` from the deposit into one download folder. The deposit also includes `README_RELEASE.txt` and `ARCHIVE_SUMMARY.csv`.
2. Verify a download (Python standard library only):

   ```sh
   python "/path/to/downloads/verify_release.py" --archives --group 01_figures_tables_heritability.zip
   ```

3. Open a terminal in the **code repository root**, beside `README.md`, `config/`, and `scripts/`. Extract each selected archive there:

   ```sh
   python -m zipfile -e "/path/to/downloads/01_figures_tables_heritability.zip" .
   ```

   Repeat with each chosen ZIP filename. Each archive already contains `data/`; extracting into `data/` would incorrectly create `data/data/`. Use a fresh checkout or preserve any locally modified input files before extraction. No configuration path changes are needed with this layout.

4. Verify the extracted group:

   ```sh
   python "/path/to/downloads/verify_release.py" --installed . --group 01_figures_tables_heritability.zip
   ```

Omit `--group` only when all six archives have been downloaded or extracted, respectively. The manifest lists every file's destination and checksum. These ZIP64 archives are supported by Python's ZIP reader.

Follow the [workflow guides](../README.md) after extraction. Saved V15 evidence in 06 is not a resumable checkpoint set for the cleaned code; it supports inspection and feature-summary reconstruction. The release was verified against code commit `9a1d239bc5d530a123a81f368ed86151a326ccbe`. Data license and final citation details will be supplied with the deposit.
