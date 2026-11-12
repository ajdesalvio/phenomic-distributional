import os
import rasterio
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import PowerNorm

VI_FUNCTIONS = {
    # BI
    "BI": lambda bands, scale: np.sqrt(((bands['R']/scale)**2 + (bands['G']/scale)**2 + (bands['B']/scale)**2) / 3),
    "GLI": lambda bands, scale: (2*(bands['G']/scale) - (bands['R']/scale) - (bands['B']/scale)) / (2*(bands['G']/scale) + (bands['R']/scale) + (bands['B']/scale)),
    "NGRDI": lambda bands, scale: ((bands['G']/scale) - (bands['R']/scale)) / ((bands['G']/scale) + (bands['R']/scale)),
    "VARI": lambda bands, scale: ((bands['G']/scale) - (bands['R']/scale)) / ((bands['G']/scale) + (bands['R']/scale) - (bands['B']/scale)),
    "BGI": lambda bands, scale: ((bands['B']/scale) / (bands['G']/scale)),
    "PSRI": lambda bands, scale: ((bands['R']/scale) - (bands['G']/scale)) / (bands['RE']/scale),
    "NDVI": lambda bands, scale: ((bands['NIR']/scale) - (bands['R']/scale)) / ((bands['NIR']/scale) + (bands['R']/scale)),
    "GNDVI": lambda bands, scale: ((bands['NIR']/scale) - (bands['G']/scale)) / ((bands['NIR']/scale) + (bands['G']/scale)),
    "RVI": lambda bands, scale: (bands['NIR']/scale) / (bands['R']/scale),
    "NDRE": lambda bands, scale: ((bands['NIR']/scale) - (bands['RE']/scale)) / ((bands['NIR']/scale) + (bands['RE']/scale)),
    "TVI": lambda bands, scale: 0.5 * (120 * ((bands['NIR']/scale) - (bands['G']/scale)) - 200 * ((bands['R']/scale) - (bands['G']/scale))),
    "CVI": lambda bands, scale: ((bands['NIR']/scale) * (bands['R']/scale)) / ((bands['G']/scale)**2),
    "EVI": lambda bands, scale: 2.5 * ((bands['NIR']/scale) - (bands['R']/scale)) / ((bands['NIR']/scale) + 6*(bands['R']/scale) - 7.5*(bands['B']/scale) + 1),
    "CIG": lambda bands, scale: ((bands['NIR']/scale) / (bands['G']/scale)) - 1,
    "CIRE": lambda bands, scale: ((bands['NIR']/scale) / (bands['RE']/scale)) - 1,
    "DVI": lambda bands, scale: (bands['NIR']/scale) - (bands['RE']/scale),
    "BCC": lambda bands, scale: ((bands['B']/scale) / ((bands['R']/scale) + (bands['G']/scale) + (bands['B']/scale))),
    "CIVE": lambda bands, scale: ((0.441*(bands['R']/scale)) - (0.811*(bands['G']/scale)) + (0.385*(bands['B']/scale)) + 18.78745),
    "COM1": lambda bands, scale: (((2*(bands['G']/scale)) - (bands['R']/scale) - (bands['B']/scale)) +
                           ((0.441*(bands['R']/scale)) - (0.811*(bands['G']/scale)) + (0.385*(bands['B']/scale)) + 18.78745) +
                           ((3*(bands['G']/scale)) - (2.4*(bands['R']/scale)) - (bands['B']/scale)) +
                           ((bands['G']/scale) / ((bands['R']/scale)**0.667 * (bands['B']/scale)**0.334))),
    "COM2": lambda bands, scale: ((0.36 * ((2*(bands['G']/scale)) - (bands['R']/scale) - (bands['B']/scale))) +
                           (0.47 * ((0.441*(bands['R']/scale)) - (0.811*(bands['G']/scale)) + (0.385*(bands['B']/scale)) + 18.78745)) +
                           (0.17 * ((bands['G']/scale) / ((bands['R']/scale)**0.667 * (bands['B']/scale)**0.334)))),
    "EBI": lambda bands, scale: (((bands['B']/scale) - (bands['G']/scale)) / ((bands['B']/scale) - (bands['R']/scale))),
    "EGI": lambda bands, scale: (((bands['G']/scale) - (bands['R']/scale)) / ((bands['R']/scale) - (bands['B']/scale))),
    "ExG_Old": lambda bands, scale: ((2*(bands['G']/scale)) - (bands['R']/scale) - (bands['B']/scale)),
    "ExG_TP": lambda bands, scale: 2 * (bands['G']/scale / (bands['R']/scale + bands['G']/scale + bands['B']/scale)) - 
                            (bands['R']/scale / (bands['R']/scale + bands['G']/scale + bands['B']/scale)) - 
                            (bands['B']/scale / (bands['R']/scale + bands['G']/scale + bands['B']/scale)),
    "ExG2": lambda bands, scale: (((2*(bands['G']/scale)) - (bands['R']/scale) - (bands['B']/scale)) / 
                           ((bands['R']/scale) + (bands['G']/scale) + (bands['B']/scale))),
    "ExGR": lambda bands, scale: ((3*(bands['G']/scale)) - (2.4*(bands['R']/scale)) - (bands['B']/scale)),
    "EXR": lambda bands, scale: ((1.4*(bands['R']/scale)) - (bands['G']/scale)),
    "GminusB": lambda bands, scale: ((bands['G']/scale) - (bands['B']/scale)),
    "GminusR": lambda bands, scale: ((bands['G']/scale) - (bands['R']/scale)),
    "GdivB": lambda bands, scale: ((bands['G']/scale) / (bands['B']/scale)),
    "GdivR": lambda bands, scale: ((bands['G']/scale) / (bands['R']/scale)),
    "GCC": lambda bands, scale: ((bands['G']/scale) / ((bands['R']/scale) + (bands['G']/scale) + (bands['B']/scale))),
    "MExG": lambda bands, scale: ((1.262*(bands['G']/scale)) - (0.884*(bands['R']/scale)) - (0.311*(bands['B']/scale))),
    "MGVRI": lambda bands, scale: (((bands['G']/scale)**2 - (bands['R']/scale)**2) / ((bands['G']/scale)**2 + (bands['R']/scale)**2)),
    "NDI": lambda bands, scale: 128 * ((((bands['G']/scale) - (bands['R']/scale)) / ((bands['G']/scale) + (bands['R']/scale))) + 1),
    "NDBRI": lambda bands, scale: ((bands['R']/scale) - (bands['B']/scale)) / ((bands['R']/scale) + (bands['B']/scale)),
    "NGBDI": lambda bands, scale: ((bands['G']/scale) - (bands['B']/scale)) / ((bands['G']/scale) + (bands['B']/scale)),
    "RminusB": lambda bands, scale: ((bands['R']/scale) - (bands['B']/scale)),
    "RdivB": lambda bands, scale: ((bands['R']/scale) / (bands['B']/scale)),
    "RCC": lambda bands, scale: ((bands['R']/scale) / ((bands['R']/scale) + (bands['G']/scale) + (bands['B']/scale))),
    "MRCCbyAlper": lambda bands, scale: ((bands['R']/scale)**3 / ((bands['R']/scale) + (bands['G']/scale) + (bands['B']/scale))),
    "RGBVI": lambda bands, scale: ((((bands['G']/scale)**2) - ((bands['R']/scale)*(bands['B']/scale))) / (((bands['G']/scale)**2) + ((bands['R']/scale)*(bands['B']/scale)))),
    "TGI": lambda bands, scale: ((bands['G']/scale) - (0.39*(bands['R']/scale) - 0.69*(bands['B']/scale))),
    "VEG": lambda bands, scale: ((bands['G']/scale) / ((bands['R']/scale)**0.667 * (bands['B']/scale)**0.334)),
    "NRMBI": lambda bands, scale: (((bands['R']/scale) - (bands['B']/scale)) / (bands['G']/scale)),
    "GRVI": lambda bands, scale: ((bands['NIR']/scale) / (bands['G']/scale)),
    "GDVI": lambda bands, scale: ((bands['NIR']/scale) - (bands['G']/scale)),
    "GWDRVI": lambda bands, scale: (((0.12*(bands['NIR']/scale)) - (bands['G']/scale)) / ((0.12*(bands['NIR']/scale)) + (bands['G']/scale))),
    "MSR_G": lambda bands, scale: ((((bands['NIR']/scale) / (bands['G']/scale)) - 1) / np.sqrt(((bands['NIR']/scale) / (bands['G']/scale)) + 1)),
    "GSAVI": lambda bands, scale: 1.5 * (((bands['NIR']/scale) - (bands['G']/scale)) / ((bands['NIR']/scale) + (bands['G']/scale) + 0.5)),
    "MSAVI": lambda bands, scale: 0.5 * (2*(bands['NIR']/scale) + 1 - np.sqrt((2*(bands['NIR']/scale) + 1)**2 - 8*((bands['NIR']/scale) - (bands['R']/scale)))),
    "MSR": lambda bands, scale: (((bands['NIR']/scale) / (bands['R']/scale)) - 1) / np.sqrt(((bands['NIR']/scale) / (bands['R']/scale)) + 1),
    "TNDVI": lambda bands, scale: np.sqrt(((bands['NIR']/scale) - (bands['R']/scale)) / ((bands['NIR']/scale) + (bands['R']/scale)) + 0.5),
    "REPR": lambda bands, scale: ((bands['R']/scale) + (bands['NIR']/scale)) / 2,
    "NLI": lambda bands, scale: (((bands['NIR']/scale)**2 - (bands['R']/scale)) / ((bands['NIR']/scale)**2 + (bands['R']/scale))),
    "MNLI": lambda bands, scale: 1.5 * (((bands['NIR']/scale)**2 - (bands['R']/scale)) / ((bands['NIR']/scale)**2 + (bands['R']/scale) + 0.5)),
    "NDVIxRVI": lambda bands, scale: (((bands['NIR']/scale)**2 - (bands['R']/scale)) / ((bands['NIR']/scale) + (bands['R']/scale)**2)),
    "SAVIxSR": lambda bands, scale: (((bands['NIR']/scale)**2 - (bands['R']/scale)) / (((bands['NIR']/scale) + (bands['R']/scale) + 0.5) * (bands['R']/scale))),
    "RERVI": lambda bands, scale: ((bands['NIR']/scale) / (bands['RE']/scale)),
    "RERDVI": lambda bands, scale: ((bands['NIR']/scale) - (bands['RE']/scale)) / np.sqrt((bands['NIR']/scale) + (bands['RE']/scale)),
    "REWDRVI": lambda bands, scale: ((0.12*(bands['NIR']/scale)) - (bands['RE']/scale)) / ((0.12*(bands['NIR']/scale)) + (bands['RE']/scale)),
    "RESAVI": lambda bands, scale: (1.5 * ((bands['NIR']/scale) - (bands['RE']/scale))) / ((bands['NIR']/scale) + (bands['RE']/scale) + 0.5),
    "REOSAVI": lambda bands, scale: (1 + 0.16) * ((bands['NIR']/scale) - (bands['RE']/scale)) / ((bands['NIR']/scale) + (bands['RE']/scale) + 0.16),
    "MRESAVI": lambda bands, scale: 0.5 * (2*(bands['NIR']/scale) + 1 - np.sqrt((2*(bands['NIR']/scale) + 1)**2 - 8*((bands['NIR']/scale) - (bands['RE']/scale)))),
    "MSRRE": lambda bands, scale: (((bands['NIR']/scale) / (bands['RE']/scale)) - 1) / np.sqrt(((bands['NIR']/scale) / (bands['RE']/scale)) + 1),
    "REVIopt": lambda bands, scale: 100 * (np.log10((bands['NIR']/scale)) - np.log10((bands['RE']/scale))),
    "GRDVI": lambda bands, scale: ((bands['NIR']/scale) - (bands['G']/scale)) / np.sqrt((bands['NIR']/scale) + (bands['G']/scale)),
    "GOSAVI": lambda bands, scale: (1 + 0.16) * ((bands['NIR']/scale) - (bands['G']/scale)) / ((bands['NIR']/scale) + (bands['G']/scale) + 0.16),
    "MGSAVI": lambda bands, scale: 0.5 * (2*(bands['NIR']/scale) + 1 - np.sqrt((2*(bands['NIR']/scale) + 1)**2 - 8*((bands['NIR']/scale) - (bands['G']/scale)))),
    "MSRG": lambda bands, scale: (((bands['NIR']/scale) / (bands['G']/scale)) - 1) / np.sqrt(((bands['NIR']/scale) / (bands['G']/scale)) + 1),
    "RENDVI": lambda bands, scale: ((bands['RE']/scale) - (bands['R']/scale)) / ((bands['RE']/scale) + (bands['R']/scale)),
    "RESR": lambda bands, scale: ((bands['RE']/scale) / (bands['R']/scale)),
    "MREDVI": lambda bands, scale: ((bands['RE']/scale) - (bands['R']/scale)),
    "MSRGR": lambda bands, scale: np.sqrt((bands['G']/scale) / (bands['R']/scale)),
    "TNDGR": lambda bands, scale: np.sqrt(((bands['G']/scale) - (bands['R']/scale)) / ((bands['G']/scale) + (bands['R']/scale)) + 0.5),
    "MTCI": lambda bands, scale: ((bands['NIR']/scale) - (bands['RE']/scale)) / ((bands['RE']/scale) - (bands['R']/scale)),
    "DATT": lambda bands, scale: ((bands['NIR']/scale) - (bands['RE']/scale)) / ((bands['NIR']/scale) - (bands['R']/scale)),
    "mNDVI1": lambda bands, scale: ((bands['NIR']/scale) - (bands['R']/scale) + 2*(bands['G']/scale)) / ((bands['NIR']/scale) + (bands['R']/scale) - 2*(bands['G']/scale)),
    "MCARI": lambda bands, scale: (((bands['RE']/scale) - (bands['R']/scale)) - 0.2*((bands['RE']/scale) - (bands['G']/scale))) * ((bands['RE']/scale) / (bands['R']/scale)),
    "MCARI1": lambda bands, scale: 1.2 * (2.5*((bands['NIR']/scale) - (bands['R']/scale)) - 1.3*((bands['NIR']/scale) - (bands['G']/scale))),
    "MCARI2": lambda bands, scale: 1.5 * (2.5*((bands['NIR']/scale) - (bands['R']/scale)) - 1.3*((bands['NIR']/scale) - (bands['G']/scale))) / np.sqrt((2*(bands['NIR']/scale) + 1)**2 - (6*(bands['NIR']/scale) - 5*np.sqrt((bands['R']/scale))) - 0.5),
    "TCARI": lambda bands, scale: 3 * (((bands['RE']/scale) - (bands['R']/scale)) - 0.2*((bands['RE']/scale) - (bands['G']/scale)) * ((bands['RE']/scale) / (bands['R']/scale))),
    "TCI": lambda bands, scale: 1.2*((bands['RE']/scale) - (bands['G']/scale)) - 1.5*((bands['R']/scale) - (bands['G']/scale)) * np.sqrt((bands['RE']/scale) / (bands['R']/scale)),
    "NRBDI": lambda bands, scale: ((bands['RE']/scale) - (bands['B']/scale)) / ((bands['RE']/scale) + (bands['B']/scale)),
    "NNBDI": lambda bands, scale: ((bands['NIR']/scale) - (bands['B']/scale)) / ((bands['NIR']/scale) + (bands['B']/scale)),
    "NNIRI": lambda bands, scale: (bands['NIR']/scale) / ((bands['NIR']/scale) + (bands['RE']/scale) + (bands['R']/scale)),
    "NREI": lambda bands, scale: (bands['RE']/scale) / ((bands['NIR']/scale) + (bands['RE']/scale) + (bands['R']/scale)),
    "NRI": lambda bands, scale: (bands['R']/scale) / ((bands['NIR']/scale) + (bands['RE']/scale) + (bands['R']/scale)),
    "MNDI": lambda bands, scale: ((bands['NIR']/scale) - (bands['RE']/scale)) / ((bands['NIR']/scale) + (bands['RE']/scale) - 2*(bands['R']/scale)),
    "MEVI": lambda bands, scale: 2.5 * ((bands['NIR']/scale) - (bands['RE']/scale)) / ((bands['NIR']/scale) + 6*(bands['RE']/scale) - 7.5*(bands['G']/scale) + 1),
    "MNDRE2": lambda bands, scale: ((bands['NIR']/scale) - (bands['RE']/scale) + 2*(bands['R']/scale)) / ((bands['NIR']/scale) + (bands['RE']/scale) - 2*(bands['R']/scale)),
    "SAVI": lambda bands, scale: 1.5 * (((bands['NIR']/scale) - (bands['R']/scale)) / ((bands['NIR']/scale) + (bands['R']/scale) + 0.5)),
    "IPVI": lambda bands, scale: ((((bands['NIR']/scale) - (bands['R']/scale)) / ((bands['NIR']/scale) + (bands['R']/scale))) + 1) / 2,
    "REGDVI": lambda bands, scale: (bands['RE']/scale) - (bands['G']/scale),
    "NGI": lambda bands, scale: (bands['G']/scale) / ((bands['NIR']/scale) + (bands['RE']/scale) + (bands['G']/scale)),
    "GNREI": lambda bands, scale: (bands['RE']/scale) / ((bands['NIR']/scale) + (bands['RE']/scale) + (bands['G']/scale)),
    "NNIR": lambda bands, scale: (bands['NIR']/scale) / ((bands['NIR']/scale) + (bands['RE']/scale) + (bands['G']/scale)),
    "REGRVI": lambda bands, scale: (bands['RE']/scale) / (bands['G']/scale)
}


# Function to compute a global ECDF for a specific flight_date
def compute_global_ecdf(main_folder, flight_date, vi_func, band_map, scale, allowed_files=None):
    """
    Now only uses files that both contain flight_date AND are in allowed_files (if provided).
    Quietly returns (empty_array, 0) if nothing matches.
    """
    all_fnames = os.listdir(main_folder)
    # select only those names that match our two criteria
    matched = [
        fn for fn in all_fnames
        if fn.lower().endswith('.tif')
           and flight_date in fn
           and (allowed_files is None or fn in allowed_files)
    ]
    if not matched:
        print(f"  ‣ Warning: no env‐filtered images found for date {flight_date}")
        return np.array([]), 0

    global_values = []
    for fn in matched:
        fp = os.path.join(main_folder, fn)
        with rasterio.open(fp) as src:
            img = src.read()
        black_mask = np.all(img == 0, axis=0)
        nb = ~black_mask
        band_arrays = {k: img[idx].astype(float) for k, idx in band_map.items()}
        bn = {k: arr[nb] for k, arr in band_arrays.items()}
        vi_vals = vi_func(bn)
        if np.isscalar(vi_vals):
            vi_vals = np.full(nb.sum(), vi_vals, dtype=float)
        else:
            vi_vals = np.array(vi_vals, dtype=float)
        global_values.append(vi_vals)

    all_vals = np.concatenate(global_values)
    return np.sort(all_vals), len(all_vals)

# Visualization function that uses the global ECDF
def visualize_heritability_global(image_path, heritability_df, dat, vegetation_index, vi_func, band_map, 
                                  global_sorted, total_count, gamma, scale=65535):
    """
    For an individual image, compute its VI, assign a global quantile (based on the global_sorted array)
    to each non-black pixel, then interpolate a heritability value for each pixel.
    Displays side-by-side the RGB image and the heritability overlay.
    
    Parameters:
      image_path (str): Path to the focal image.
      heritability_df (DataFrame): Contains columns 'DAP', 'Vegetation.Index', 'Probability', and 'Heritability'.
      dat (numeric): The DAP used to filter heritability_df.
      vegetation_index (str): The VI used.
      vi_func (callable): Function to compute VI.
      band_map (dict): Mapping for band indices.
      global_sorted (np.array): Sorted global VI values (from all images of the flight date).
      total_count (int): Total number of VI samples in the global distribution.
      gamma (int): Gamma > 1 emphasizes variation in the HIGH heritability range; gamma = 1 is linear, gamma < 1 emphasizes
          variation in the LOW heritability range
      scale (numeric): Scale factor for display.
    
    Returns:
      fig, vi_image, herit_image
    """
    # Load the focal image
    with rasterio.open(image_path) as src:
        img = src.read()
    bands, rows, cols = img.shape

    # Create an RGB composite (assumes multispectral with R=2, G=1, B=0)
    rgb_img = np.stack([img[2, :, :], img[1, :, :], img[0, :, :]], axis=-1).astype(np.float32)
    rgb_img = np.clip(rgb_img / scale, 0, 1)
    
    # Mask non-black pixels
    black_mask = np.all(img == 0, axis=0)
    non_black_mask = ~black_mask
    
    # Build band arrays
    band_arrays = {k: img[idx, :, :].astype(float) for k, idx in band_map.items()}
    vi_image = np.full((rows, cols), np.nan)
    bands_non_black = {k: arr[non_black_mask] for k, arr in band_arrays.items()}
    vi_vals = vi_func(bands_non_black)
    if np.isscalar(vi_vals):
        vi_vals = np.full(np.sum(non_black_mask), vi_vals, dtype=float)
    else:
        vi_vals = np.array(vi_vals, dtype=float)
    vi_image[non_black_mask] = vi_vals

    # For each non-black pixel in the focal image, obtain its rank in the global_sorted array
    local_values = vi_image[non_black_mask]
    pixel_ranks = np.searchsorted(global_sorted, local_values, side='right')
    global_probabilities = pixel_ranks / float(total_count)
    
    # Subset the heritability DataFrame
    df_sub = heritability_df[
        (heritability_df['DAP'] == dat) & (heritability_df['Vegetation.Index'] == vegetation_index)
    ]
    if df_sub.empty:
        raise ValueError(f"No heritability data available for DAP = {dat} and VI = {vegetation_index}")
    df_sub = df_sub.sort_values('Probability')
    prob_vals = df_sub['Probability'].values
    herit_vals = df_sub['Heritability'].values

    # Interpolate heritability values
    herit_interpolated = np.interp(global_probabilities, prob_vals, herit_vals)
    herit_image = np.full((rows, cols), np.nan)
    herit_image[non_black_mask] = herit_interpolated
    herit_masked = np.ma.masked_invalid(herit_image)
    
    extent = (0, cols, rows, 0)
    fig, axes = plt.subplots(1, 2, figsize=(18,8), gridspec_kw={'width_ratios':[1.0,1.1]}, sharex=True, sharey=True)
    
    axes[0].imshow(rgb_img, extent=extent, interpolation='none')
    axes[0].set_title('RGB Image')
    axes[0].axis('off')
    axes[0].set_aspect('equal')
    
    # Right subplot: RGB image with heritability overlay
    axes[1].imshow(rgb_img, interpolation='none', extent=extent)
    cmap_heat = plt.cm.jet.copy()
    cmap_heat.set_bad(color='black')
    norm = PowerNorm(gamma=gamma, vmin=0, vmax=1) # Changing gamma to experiment with showing subtle differences in variation for high heritability values
    im = axes[1].imshow(herit_masked, cmap=cmap_heat, alpha=1, # Setting alpha=1 makes it completely opaque and covers the RGB image 
                          interpolation='none', norm=norm, extent=extent)
    axes[1].set_title(f'Heritability Overlay (DAP = {dat}, VI = {vegetation_index})')
    axes[1].axis('off')
    axes[1].set_aspect('equal')
    fig.colorbar(im, ax=axes[1], label='Heritability', fraction=0.046, pad=0.04)
    
    plt.tight_layout()
    return fig, vi_image, herit_image


# Example usage:
if __name__ == '__main__':
    
    # Path segment used for loading and saving files
    data_path = 'C:/Users/aaron.desalvio/Downloads/All_Distributional_Files'
    
    # Directory with all G2F images
    main_image_dir = 'E:/2020_2021_G2F_Distributional/2020_2021_G2F_TIFs_Renamed_Hue_Masked_Images/2020_2021_G2F_TIFs_Renamed_Hue_Masked_Images'
    
    # Path to save heatmaps
    export_path = (f'{data_path}')
    
    # VI / scale setup (same as before)
    scale = 255
    selected_vi = 'ExGR'
    vi_func = lambda bands: VI_FUNCTIONS[selected_vi](bands, scale)
    
    # Read the DAP-to-date conversion CSV
    conversion_csv_path = (f'{data_path}/Flight_Date_DAP_Conversion_G2F.csv')
    dat_conversion = pd.read_csv(conversion_csv_path)
    
    # Environment, plot, and DAP list
    env = 'TXH1.2020'
    plant_id = 'CS20-G2F1-031'
    selected_dat_values = [19,28,38,47,50,56,61,68,76,79,82,88,96,102,106,111,117]
    #selected_dat_values = [28]
    
    # Load data frame that will allow us to filter only TIFs that belong to a specific Env
    file_name_csv = f'{data_path}/File_Name_TIFs_JPGs.csv'
    fn_df = pd.read_csv(file_name_csv)
    # Extract just those filenames for our env
    allowed_files = fn_df.loc[fn_df['Env']==env, 'File.Name.TIF.Unified.Masked'].tolist()
    
    # Heritability lookup table
    heritability_df = pd.read_csv(f'{data_path}/VarComp_per_Quantile_Within_DAP_RGB_VIs_G2F_PedigreeOnly.csv')
    
    # Band map for multispectral
    multispec_band_map = {'R':0, 'G':1, 'B':2}
    
    # Loop
    for selected_dat in selected_dat_values:
        # 1) look up flight_date as a string, filtering by both DAP and Env
        flight_date = str(
            dat_conversion.loc[
                (dat_conversion['DAP'] == selected_dat) &
                (dat_conversion['Env'] == env),
                'Flight_Date_YYYYMMDD'
            ].values[0]
        )
        
        # 2) compute global ECDF over all images whose names contain that flight_date
        global_sorted, total_count = compute_global_ecdf(
            main_image_dir,
            flight_date,
            vi_func,
            multispec_band_map,
            scale,
            allowed_files=allowed_files
        )
        
        # 3) build the focal image path for this plant/date
        focal_image_path = os.path.join(
            main_image_dir,
            f"{flight_date}-{plant_id}_masked.tif"
        )
        
        # 4) generate the figure using the global ECDF
        fig, vi_image, herit_image = visualize_heritability_global(
            image_path=focal_image_path,
            heritability_df=heritability_df,
            dat=selected_dat,
            vegetation_index=selected_vi,
            vi_func=vi_func,
            band_map=multispec_band_map,
            global_sorted=global_sorted,
            total_count=total_count,
            gamma=1.0,
            scale=scale
        )
        
        # 5) save with dynamic filename
        out_filename = os.path.join(
            export_path,
            f"Heritability_Heatmap_{plant_id}_{env}_DAP{selected_dat}_{selected_vi}.jpg"
        )
        fig.savefig(out_filename, dpi=600, bbox_inches='tight')
        plt.close(fig)
