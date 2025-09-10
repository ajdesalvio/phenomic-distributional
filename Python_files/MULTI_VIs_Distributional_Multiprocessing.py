import rasterio
import numpy as np
import os
import pandas as pd
import concurrent.futures
import multiprocessing.reduction as reduction
import cloudpickle
import csv

reduction.ForkingPickler.dumps = cloudpickle.dumps

def get_multiband_file_paths(folder):
    """
    Instead of loading images immediately, return file paths and base filenames.
    """
    valid_extensions = ('.tif', '.tiff')
    file_paths = []
    img_files = []
    for filename in os.listdir(folder):
        if filename.lower().endswith(valid_extensions):
            file_paths.append(os.path.join(folder, filename))
            img_files.append(os.path.splitext(filename)[0])
    return file_paths, img_files

def safe_mean_median(array):
    valid = ~np.isnan(array) & ~np.isinf(array)
    if np.any(valid):
        array = array[valid]
        return np.mean(array), np.median(array)
    else:
        return np.nan, np.nan


base_vi_functions = {
    "BI": lambda scale, bands: np.sqrt(((bands['R']/scale)**2 + (bands['G']/scale)**2 + (bands['B']/scale)**2) / 3),
    "GLI": lambda scale, bands: (2*(bands['G']/scale) - (bands['R']/scale) - (bands['B']/scale)) / (2*(bands['G']/scale) + (bands['R']/scale) + (bands['B']/scale)),
    "NGRDI": lambda scale, bands: ((bands['G']/scale) - (bands['R']/scale)) / ((bands['G']/scale) + (bands['R']/scale)),
    "VARI": lambda scale, bands: ((bands['G']/scale) - (bands['R']/scale)) / ((bands['G']/scale) + (bands['R']/scale) - (bands['B']/scale)),
    "BGI": lambda scale, bands: ((bands['B']/scale) / (bands['G']/scale)),
    "PSRI": lambda scale, bands: ((bands['R']/scale) - (bands['G']/scale)) / (bands['RE']/scale),
    "NDVI": lambda scale, bands: ((bands['NIR']/scale) - (bands['R']/scale)) / ((bands['NIR']/scale) + (bands['R']/scale)),
    "GNDVI": lambda scale, bands: ((bands['NIR']/scale) - (bands['G']/scale)) / ((bands['NIR']/scale) + (bands['G']/scale)),
    "RVI": lambda scale, bands: (bands['NIR']/scale) / (bands['R']/scale),
    "NDRE": lambda scale, bands: ((bands['NIR']/scale) - (bands['RE']/scale)) / ((bands['NIR']/scale) + (bands['RE']/scale)),
    "TVI": lambda scale, bands: 0.5 * (120 * ((bands['NIR']/scale) - (bands['G']/scale)) - 200 * ((bands['R']/scale) - (bands['G']/scale))),
    "CVI": lambda scale, bands: ((bands['NIR']/scale) * (bands['R']/scale)) / ((bands['G']/scale)**2),
    "EVI": lambda scale, bands: 2.5 * ((bands['NIR']/scale) - (bands['R']/scale)) / ((bands['NIR']/scale) + 6*(bands['R']/scale) - 7.5*(bands['B']/scale) + 1),
    "CIG": lambda scale, bands: ((bands['NIR']/scale) / (bands['G']/scale)) - 1,
    "CIRE": lambda scale, bands: ((bands['NIR']/scale) / (bands['RE']/scale)) - 1,
    "DVI": lambda scale, bands: (bands['NIR']/scale) - (bands['RE']/scale),
    "BCC": lambda scale, bands: ((bands['B']/scale) / ((bands['R']/scale) + (bands['G']/scale) + (bands['B']/scale))),
    "CIVE": lambda scale, bands: ((0.441*(bands['R']/scale)) - (0.811*(bands['G']/scale)) + (0.385*(bands['B']/scale)) + 18.78745),
    "COM1": lambda scale, bands: (((2*(bands['G']/scale)) - (bands['R']/scale) - (bands['B']/scale)) +
                           ((0.441*(bands['R']/scale)) - (0.811*(bands['G']/scale)) + (0.385*(bands['B']/scale)) + 18.78745) +
                           ((3*(bands['G']/scale)) - (2.4*(bands['R']/scale)) - (bands['B']/scale)) +
                           ((bands['G']/scale) / ((bands['R']/scale)**0.667 * (bands['B']/scale)**0.334))),
    "COM2": lambda scale, bands: ((0.36 * ((2*(bands['G']/scale)) - (bands['R']/scale) - (bands['B']/scale))) +
                           (0.47 * ((0.441*(bands['R']/scale)) - (0.811*(bands['G']/scale)) + (0.385*(bands['B']/scale)) + 18.78745)) +
                           (0.17 * ((bands['G']/scale) / ((bands['R']/scale)**0.667 * (bands['B']/scale)**0.334)))),
    "EBI": lambda scale, bands: (((bands['B']/scale) - (bands['G']/scale)) / ((bands['B']/scale) - (bands['R']/scale))),
    "EGI": lambda scale, bands: (((bands['G']/scale) - (bands['R']/scale)) / ((bands['R']/scale) - (bands['B']/scale))),
    "ExG_Old": lambda scale, bands: ((2*(bands['G']/scale)) - (bands['R']/scale) - (bands['B']/scale)),
    "ExG_TP": lambda scale, bands: 2 * (bands['G']/scale / (bands['R']/scale + bands['G']/scale + bands['B']/scale)) - 
                            (bands['R']/scale / (bands['R']/scale + bands['G']/scale + bands['B']/scale)) - 
                            (bands['B']/scale / (bands['R']/scale + bands['G']/scale + bands['B']/scale)),
    "ExG2": lambda scale, bands: (((2*(bands['G']/scale)) - (bands['R']/scale) - (bands['B']/scale)) / 
                           ((bands['R']/scale) + (bands['G']/scale) + (bands['B']/scale))),
    "ExGR": lambda scale, bands: ((3*(bands['G']/scale)) - (2.4*(bands['R']/scale)) - (bands['B']/scale)),
    "EXR": lambda scale, bands: ((1.4*(bands['R']/scale)) - (bands['G']/scale)),
    "GminusB": lambda scale, bands: ((bands['G']/scale) - (bands['B']/scale)),
    "GminusR": lambda scale, bands: ((bands['G']/scale) - (bands['R']/scale)),
    "GdivB": lambda scale, bands: ((bands['G']/scale) / (bands['B']/scale)),
    "GdivR": lambda scale, bands: ((bands['G']/scale) / (bands['R']/scale)),
    "GCC": lambda scale, bands: ((bands['G']/scale) / ((bands['R']/scale) + (bands['G']/scale) + (bands['B']/scale))),
    "MExG": lambda scale, bands: ((1.262*(bands['G']/scale)) - (0.884*(bands['R']/scale)) - (0.311*(bands['B']/scale))),
    "MGVRI": lambda scale, bands: (((bands['G']/scale)**2 - (bands['R']/scale)**2) / ((bands['G']/scale)**2 + (bands['R']/scale)**2)),
    "NDI": lambda scale, bands: 128 * ((((bands['G']/scale) - (bands['R']/scale)) / ((bands['G']/scale) + (bands['R']/scale))) + 1),
    "NDBRI": lambda scale, bands: ((bands['R']/scale) - (bands['B']/scale)) / ((bands['R']/scale) + (bands['B']/scale)),
    "NGBDI": lambda scale, bands: ((bands['G']/scale) - (bands['B']/scale)) / ((bands['G']/scale) + (bands['B']/scale)),
    "RminusB": lambda scale, bands: ((bands['R']/scale) - (bands['B']/scale)),
    "RdivB": lambda scale, bands: ((bands['R']/scale) / (bands['B']/scale)),
    "RCC": lambda scale, bands: ((bands['R']/scale) / ((bands['R']/scale) + (bands['G']/scale) + (bands['B']/scale))),
    "MRCCbyAlper": lambda scale, bands: ((bands['R']/scale)**3) / ((bands['R']/scale) + (bands['G']/scale) + (bands['B']/scale)),
    "RGBVI": lambda scale, bands: ((((bands['G']/scale)**2) - ((bands['R']/scale)*(bands['B']/scale))) / (((bands['G']/scale)**2) + ((bands['R']/scale)*(bands['B']/scale)))),
    "TGI": lambda scale, bands: ((bands['G']/scale) - (0.39*(bands['R']/scale) - 0.69*(bands['B']/scale))),
    "VEG": lambda scale, bands: ((bands['G']/scale) / ((bands['R']/scale)**0.667 * (bands['B']/scale)**0.334)),
    "NRMBI": lambda scale, bands: (((bands['R']/scale) - (bands['B']/scale)) / (bands['G']/scale)),
    "GRVI": lambda scale, bands: ((bands['NIR']/scale) / (bands['G']/scale)),
    "GDVI": lambda scale, bands: ((bands['NIR']/scale) - (bands['G']/scale)),
    "GWDRVI": lambda scale, bands: (((0.12*(bands['NIR']/scale)) - (bands['G']/scale)) / ((0.12*(bands['NIR']/scale)) + (bands['G']/scale))),
    "MSR_G": lambda scale, bands: ((((bands['NIR']/scale) / (bands['G']/scale)) - 1) / np.sqrt(((bands['NIR']/scale) / (bands['G']/scale)) + 1)),
    "GSAVI": lambda scale, bands: 1.5 * (((bands['NIR']/scale) - (bands['G']/scale)) / ((bands['NIR']/scale) + (bands['G']/scale) + 0.5)),
    "MSAVI": lambda scale, bands: 0.5 * (2*(bands['NIR']/scale) + 1 - np.sqrt((2*(bands['NIR']/scale) + 1)**2 - 8*((bands['NIR']/scale) - (bands['R']/scale)))),
    "MSR": lambda scale, bands: (((bands['NIR']/scale) / (bands['R']/scale)) - 1) / np.sqrt(((bands['NIR']/scale) / (bands['R']/scale)) + 1),
    "TNDVI": lambda scale, bands: np.sqrt(((bands['NIR']/scale) - (bands['R']/scale)) / ((bands['NIR']/scale) + (bands['R']/scale)) + 0.5),
    "REPR": lambda scale, bands: ((bands['R']/scale) + (bands['NIR']/scale)) / 2,
    "NLI": lambda scale, bands: (((bands['NIR']/scale)**2 - (bands['R']/scale)) / ((bands['NIR']/scale)**2 + (bands['R']/scale))),
    "MNLI": lambda scale, bands: 1.5 * (((bands['NIR']/scale)**2 - (bands['R']/scale)) / ((bands['NIR']/scale)**2 + (bands['R']/scale) + 0.5)),
    "NDVIxRVI": lambda scale, bands: (((bands['NIR']/scale)**2 - (bands['R']/scale)) / ((bands['NIR']/scale) + (bands['R']/scale)**2)),
    "SAVIxSR": lambda scale, bands: (((bands['NIR']/scale)**2 - (bands['R']/scale)) / (((bands['NIR']/scale) + (bands['R']/scale) + 0.5) * (bands['R']/scale))),
    "RERVI": lambda scale, bands: ((bands['NIR']/scale) / (bands['RE']/scale)),
    "RERDVI": lambda scale, bands: ((bands['NIR']/scale) - (bands['RE']/scale)) / np.sqrt((bands['NIR']/scale) + (bands['RE']/scale)),
    "REWDRVI": lambda scale, bands: ((0.12*(bands['NIR']/scale)) - (bands['RE']/scale)) / ((0.12*(bands['NIR']/scale)) + (bands['RE']/scale)),
    "RESAVI": lambda scale, bands: (1.5 * ((bands['NIR']/scale) - (bands['RE']/scale))) / ((bands['NIR']/scale) + (bands['RE']/scale) + 0.5),
    "REOSAVI": lambda scale, bands: (1 + 0.16) * ((bands['NIR']/scale) - (bands['RE']/scale)) / ((bands['NIR']/scale) + (bands['RE']/scale) + 0.16),
    "MRESAVI": lambda scale, bands: 0.5 * (2*(bands['NIR']/scale) + 1 - np.sqrt((2*(bands['NIR']/scale) + 1)**2 - 8*((bands['NIR']/scale) - (bands['RE']/scale)))),
    "MSRRE": lambda scale, bands: (((bands['NIR']/scale) / (bands['RE']/scale)) - 1) / np.sqrt(((bands['NIR']/scale) / (bands['RE']/scale)) + 1),
    "REVIopt": lambda scale, bands: 100 * (np.log10((bands['NIR']/scale)) - np.log10((bands['RE']/scale))),
    "GRDVI": lambda scale, bands: ((bands['NIR']/scale) - (bands['G']/scale)) / np.sqrt((bands['NIR']/scale) + (bands['G']/scale)),
    "GOSAVI": lambda scale, bands: (1 + 0.16) * ((bands['NIR']/scale) - (bands['G']/scale)) / ((bands['NIR']/scale) + (bands['G']/scale) + 0.16),
    "MGSAVI": lambda scale, bands: 0.5 * (2*(bands['NIR']/scale) + 1 - np.sqrt((2*(bands['NIR']/scale) + 1)**2 - 8*((bands['NIR']/scale) - (bands['G']/scale)))),
    "MSRG": lambda scale, bands: (((bands['NIR']/scale) / (bands['G']/scale)) - 1) / np.sqrt(((bands['NIR']/scale) / (bands['G']/scale)) + 1),
    "RENDVI": lambda scale, bands: ((bands['RE']/scale) - (bands['R']/scale)) / ((bands['RE']/scale) + (bands['R']/scale)),
    "RESR": lambda scale, bands: ((bands['RE']/scale) / (bands['R']/scale)),
    "MREDVI": lambda scale, bands: ((bands['RE']/scale) - (bands['R']/scale)),
    "MSRGR": lambda scale, bands: np.sqrt((bands['G']/scale) / (bands['R']/scale)),
    "TNDGR": lambda scale, bands: np.sqrt(((bands['G']/scale) - (bands['R']/scale)) / ((bands['G']/scale) + (bands['R']/scale)) + 0.5),
    "MTCI": lambda scale, bands: ((bands['NIR']/scale) - (bands['RE']/scale)) / ((bands['RE']/scale) - (bands['R']/scale)),
    "DATT": lambda scale, bands: ((bands['NIR']/scale) - (bands['RE']/scale)) / ((bands['NIR']/scale) - (bands['R']/scale)),
    "mNDVI1": lambda scale, bands: ((bands['NIR']/scale) - (bands['R']/scale) + 2*(bands['G']/scale)) / ((bands['NIR']/scale) + (bands['R']/scale) - 2*(bands['G']/scale)),
    "MCARI": lambda scale, bands: (((bands['RE']/scale) - (bands['R']/scale)) - 0.2*((bands['RE']/scale) - (bands['G']/scale))) * ((bands['RE']/scale) / (bands['R']/scale)),
    "MCARI1": lambda scale, bands: 1.2 * (2.5*((bands['NIR']/scale) - (bands['R']/scale)) - 1.3*((bands['NIR']/scale) - (bands['G']/scale))),
    "MCARI2": lambda scale, bands: 1.5 * (2.5*((bands['NIR']/scale) - (bands['R']/scale)) - 1.3*((bands['NIR']/scale) - (bands['G']/scale))) / np.sqrt((2*(bands['NIR']/scale) + 1)**2 - (6*(bands['NIR']/scale) - 5*np.sqrt((bands['R']/scale))) - 0.5),
    "TCARI": lambda scale, bands: 3 * (((bands['RE']/scale) - (bands['R']/scale)) - 0.2*((bands['RE']/scale) - (bands['G']/scale)) * ((bands['RE']/scale) / (bands['R']/scale))),
    "TCI": lambda scale, bands: 1.2*((bands['RE']/scale) - (bands['G']/scale)) - 1.5*((bands['R']/scale) - (bands['G']/scale)) * np.sqrt((bands['RE']/scale) / (bands['R']/scale)),
    "NRBDI": lambda scale, bands: ((bands['RE']/scale) - (bands['B']/scale)) / ((bands['RE']/scale) + (bands['B']/scale)),
    "NNBDI": lambda scale, bands: ((bands['NIR']/scale) - (bands['B']/scale)) / ((bands['NIR']/scale) + (bands['B']/scale)),
    "NNIRI": lambda scale, bands: (bands['NIR']/scale) / ((bands['NIR']/scale) + (bands['RE']/scale) + (bands['R']/scale)),
    "NREI": lambda scale, bands: (bands['RE']/scale) / ((bands['NIR']/scale) + (bands['RE']/scale) + (bands['R']/scale)),
    "NRI": lambda scale, bands: (bands['R']/scale) / ((bands['NIR']/scale) + (bands['RE']/scale) + (bands['R']/scale)),
    "MNDI": lambda scale, bands: ((bands['NIR']/scale) - (bands['RE']/scale)) / ((bands['NIR']/scale) + (bands['RE']/scale) - 2*(bands['R']/scale)),
    "MEVI": lambda scale, bands: 2.5 * ((bands['NIR']/scale) - (bands['RE']/scale)) / ((bands['NIR']/scale) + 6*(bands['RE']/scale) - 7.5*(bands['G']/scale) + 1),
    "MNDRE2": lambda scale, bands: ((bands['NIR']/scale) - (bands['RE']/scale) + 2*(bands['R']/scale)) / ((bands['NIR']/scale) + (bands['RE']/scale) - 2*(bands['R']/scale)),
    "SAVI": lambda scale, bands: 1.5 * (((bands['NIR']/scale) - (bands['R']/scale)) / ((bands['NIR']/scale) + (bands['R']/scale) + 0.5)),
    "IPVI": lambda scale, bands: ((((bands['NIR']/scale) - (bands['R']/scale)) / ((bands['NIR']/scale) + (bands['R']/scale))) + 1) / 2,
    "REGDVI": lambda scale, bands: (bands['RE']/scale) - (bands['G']/scale),
    "NGI": lambda scale, bands: (bands['G']/scale) / ((bands['NIR']/scale) + (bands['RE']/scale) + (bands['G']/scale)),
    "GNREI": lambda scale, bands: (bands['RE']/scale) / ((bands['NIR']/scale) + (bands['RE']/scale) + (bands['G']/scale)),
    "NNIR": lambda scale, bands: (bands['NIR']/scale) / ((bands['NIR']/scale) + (bands['RE']/scale) + (bands['G']/scale)),
    "REGRVI": lambda scale, bands: (bands['RE']/scale) / (bands['G']/scale)
}

def process_single_multiband_image(args):
    """
    Processes a single multispectral image file (loaded on the fly) and returns its VI stats
    and (optionally) distribution rows.
    """
    file_path, fname, final_vi_dict, band_map, distributions, probability_steps, scale_param = args
    # Load image on the fly with rasterio
    with rasterio.open(file_path) as src:
        img = src.read()  # shape: (bands, rows, cols)
    bands, rows, cols = img.shape
    total_pixels = rows * cols
    black_pixel_mask = np.all(img == 0, axis=0)
    black_pixels = np.sum(black_pixel_mask)
    non_black_pixels = total_pixels - black_pixels

    row_data = {
        'file_name': fname,
        'total_pixels': total_pixels,
        'black_pixels': black_pixels,
        'non_black_pixels': non_black_pixels
    }
    distribution_rows = []

    if non_black_pixels > 0:
        # Build dictionary of band arrays for non-black pixels
        band_arrays = {}
        for name, idx in band_map.items():
            band_2d = img[idx, :, :].astype(float)
            band_arrays[name] = band_2d[~black_pixel_mask]

        for vi_name, vi_func in final_vi_dict.items():
            vi_vals = vi_func(scale_param, band_arrays)
            if isinstance(vi_vals, (int, float)):
                vi_vals = np.full(non_black_pixels, vi_vals, dtype=float)
            vi_vals = np.array(vi_vals, dtype=float)
            vi_vals[np.isinf(vi_vals)] = np.nan
            vi_vals[np.isnan(vi_vals)] = np.nan
            mean_val, median_val = safe_mean_median(vi_vals)
            row_data[f'{vi_name}_mean'] = mean_val
            row_data[f'{vi_name}_median'] = median_val

            if distributions:
                valid_vals = vi_vals[~np.isnan(vi_vals)]
                if valid_vals.size == 0:
                    quantiles = np.full(probability_steps.shape, np.nan)
                else:
                    quantiles = np.quantile(valid_vals, probability_steps)
                for prob, quant in zip(probability_steps, quantiles):
                    distribution_rows.append({
                        'file_name': fname,
                        'vi_name': vi_name,
                        'probability': prob,
                        'quantile': quant
                    })
    else:
        for vi_name in final_vi_dict.keys():
            row_data[f'{vi_name}_mean'] = np.nan
            row_data[f'{vi_name}_median'] = np.nan

    return row_data, distribution_rows

def compute_multi_vegetation_indices_parallel(
    image_folder: str,
    output_folder: str,
    output_csv_name: str,
    vi_functions: dict = None,
    replace_vi_dictionary: bool = False,
    band_map: dict = None,
    distributions: bool = False,
    output_csv_name_distribution: str = None,
    max_workers: int = None,
    batch_size: int = 100,  # Process tasks in batches to reduce memory pressure
    scale_param: float = 65535
) -> pd.DataFrame:
    """
    Computes vegetation indices from multispectral images in image_folder using parallel processing,
    and streams results directly to CSV files in output_folder with name output_csv_name.
    Additionally, if distributions is True, computes quantile distributions for each VI using
    probability steps from 0.01 to 1.00 and exports a separate CSV if output_csv_name_distribution is provided.

    Parameters
    ----------
    image_folder : str [Required]
        Path to the folder containing input image files.
    output_folder : str [Required]
        Path to the directory where the output CSV files will be saved.
    output_csv_name : str [Required]
        Filename for the main vegetation index statistics CSV (mean and median values only; no distributional data).
    vi_functions : dict [Optional]
        Custom dictionary of vegetation index functions to extend or override the defaults.
    replace_vi_dictionary : bool [Optional]
        If True, the provided vi_functions will entirely replace the default functions. If False and VIs are
        provided in vi_functions, then these user-defined VIs will be appended to the default list.
    band_map : dict [Optional]
        A dictionary that specifies the order in which bands are saved in the raster, e.g.,
        my_band_map = {'B':0, 'G':1, 'R':2, 'RE':3, 'NIR':4}. If a band map isn't supplied,
        it defaults to {'B':0, 'G':1, 'R':2}
    distributions : bool [Optional]
        If True, computes and streams quantile distribution data for each vegetation index.
    output_csv_name_distribution : str [Optional]
        Filename for the CSV containing distribution statistics (used if distributions is True).
    max_workers : int [Optional]
        Maximum number of worker processes to use for parallel processing.
    batch_size : int [Optional]
        Number of image tasks processed per batch to minimize memory overhead.
    scale_param : float [Optional]
        Scaling factor used to normalize pixel values (e.g., 255 for RGB images).

    Returns
    -------
    pd.DataFrame
        DataFrame containing the main vegetation index statistics reloaded from the output CSV.
        (Note: The results are primarily written to disk to maintain memory efficiency.)
    """
    if replace_vi_dictionary:
        final_vi_dict = vi_functions if vi_functions is not None else base_vi_functions
    else:
        final_vi_dict = {**base_vi_functions, **(vi_functions or {})}

    if band_map is None:
        band_map = {'B':0, 'G':1, 'R':2}

    file_paths, img_files = get_multiband_file_paths(image_folder)
    probability_steps = np.linspace(0.01, 1.00, 100) if distributions else None

    # Prepare CSV file for streaming main VI stats
    os.makedirs(output_folder, exist_ok=True)
    out_path = os.path.join(output_folder, output_csv_name)
    main_fieldnames = ["file_name", "total_pixels", "black_pixels", "non_black_pixels"]
    for vi in final_vi_dict.keys():
        main_fieldnames.extend([f'{vi}_mean', f'{vi}_median'])
    main_file = open(out_path, "w", newline="")
    main_csv_writer = csv.DictWriter(main_file, fieldnames=main_fieldnames)
    main_csv_writer.writeheader()

    # Prepare CSV file for streaming distribution rows if requested
    distribution_file = None
    distribution_csv_writer = None
    if distributions and output_csv_name_distribution is not None:
        out_path_dist = os.path.join(output_folder, output_csv_name_distribution)
        distribution_file = open(out_path_dist, "w", newline="")
        dist_fieldnames = ["file_name", "vi_name", "probability", "quantile"]
        distribution_csv_writer = csv.DictWriter(distribution_file, fieldnames=dist_fieldnames)
        distribution_csv_writer.writeheader()

    tasks = [
        (file_path, fname, final_vi_dict, band_map, distributions, probability_steps, scale_param)
        for file_path, fname in zip(file_paths, img_files)
    ]

    # Process tasks in batches to reduce the number of futures and recycle workers
    for i in range(0, len(tasks), batch_size):
        batch_tasks = tasks[i:i+batch_size]
        with concurrent.futures.ProcessPoolExecutor(max_workers=max_workers) as executor:
            # Using executor.map minimizes the accumulation of Future objects.
            for row_data, dist_rows in executor.map(process_single_multiband_image, batch_tasks):
                main_csv_writer.writerow(row_data)
                main_file.flush()
                if distribution_csv_writer is not None:
                    for row in dist_rows:
                        distribution_csv_writer.writerow(row)
                    distribution_file.flush()

    main_file.close()
    print(f"Saved VI stats to: {out_path}")
    if distribution_file is not None:
        distribution_file.close()
        print(f"Saved VI distribution stats to: {out_path_dist}")

    # Optionally, load the main CSV file into a DataFrame for further processing
    df = pd.read_csv(out_path)
    return df

#### Example usage ####
image_folder = 'C:/Users/aaron.desalvio/Documents/TESTING/Cotton/CS24_Masked_Images_ExG/Small_Subset'
output_folder = 'C:/Users/aaron.desalvio/Documents/TESTING/Cotton/CS24_Masked_Images_ExG'
output_csv_name = 'CS24-STEL-MULTI-VIs-Subset1.csv'
my_band_map = {'B':0, 'G':1, 'R':2, 'RE':3, 'NIR':4}
output_csv_name_distribution = 'CS24-STEL-MULTI-VIs-Distributions-Subset1.csv'

safe_num_workers = os.process_cpu_count() - 12

compute_multi_vegetation_indices_parallel(
    image_folder=image_folder,
    output_folder=output_folder,
    output_csv_name=output_csv_name,
    vi_functions=None,
    replace_vi_dictionary=False,
    band_map=my_band_map,
    distributions=True,
    output_csv_name_distribution=output_csv_name_distribution,
    max_workers=safe_num_workers,
    batch_size=100,
    scale_param=65535
)
