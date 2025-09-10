import rasterio
import numpy as np
import os
import cv2


def load_multiband_images_from_folder(folder):
    valid_extensions = ('.tif', '.tiff')
    images = []
    img_files = []

    for filename in os.listdir(folder):
        if filename.lower().endswith(valid_extensions):
            path = os.path.join(folder, filename)
            with rasterio.open(path) as src:
                data = src.read()  # Rasterio reads as (bands, rows, cols)
            images.append(data)
            img_files.append(os.path.splitext(filename)[0])

    return images, img_files

##### ExG Soil Removal Function #####
def rem_soil_index_multi(img, 
                   threshold=0.1, 
                   index_name='ExG', 
                   comparison='less',
                   custom_indices=None,
                   band_map=None):
    """
    Compute a vegetation index for each pixel and remove pixels (set to np.nan) 
    based on a threshold. Also removes all-black pixels.
    
    Parameters
    ----------
    img : np.ndarray
        Input image in shape (bands, rows, cols).
    threshold : float
        Threshold used for filtering.
    index_name : str
        Indicates which index function to use from 'custom_indices'.
    comparison : str
        'less' or 'greater' - direction of threshold check.
    custom_indices : dict
        A dictionary of lambdas {index_name: lambda R,G,B,...: ...}.
        If None, defaults to a dictionary containing only 'ExG'.
    band_map : dict
        Dictionary mapping channel names (e.g. 'B','G','R') to integer indices in 'img'. 
        For example: {'B':0, 'G':1, 'R':2, 'RE':3, 'NIR':4}.
    
    Returns
    -------
    np.ndarray: A 2D array of the computed vegetation index (rows, cols) 
        with problematic pixels set to np.nan.
    """
    # If no custom dictionary is provided, use a default dictionary with ExG
    if custom_indices is None:
        custom_indices = {
            'ExG': lambda R, G, B: 2*(G/(R+G+B)) - (R/(R+G+B)) - (B/(R+G+B))
        }

    # If the user doesn't pass a band_map, assume a simple B=0, G=1, R=2
    if band_map is None:
        band_map = {'B':0, 'G':1, 'R':2}

    if index_name not in custom_indices:
        raise ValueError(f"Index '{index_name}' not found in custom_indices. "
                         f"Available keys: {list(custom_indices.keys())}")

    # Extract individual channels based on band_map
    B = img[band_map['B'], :, :].astype(np.float32)
    G = img[band_map['G'], :, :].astype(np.float32)
    R = img[band_map['R'], :, :].astype(np.float32)

    # Calculate the vegetation index using the selected lambda
    remsoil_index = custom_indices[index_name](R, G, B)

    # Identify all-black pixels in these three bands
    black_pixel_mask = (R == 0) & (G == 0) & (B == 0)

    # Apply threshold in the chosen direction
    if comparison == 'less':
        remove_mask = remsoil_index < threshold
    elif comparison == 'greater':
        remove_mask = remsoil_index > threshold
    else:
        raise ValueError("comparison must be 'less' or 'greater'.")

    # Mark black pixels OR threshold-failing pixels as np.nan
    remsoil_index[black_pixel_mask | remove_mask] = np.nan

    return remsoil_index

def create_binary_mask_multi(remsoil_index):
    """
    Create a binary mask from the vegetation index.
    Valid pixels are 1, and NA (numpy.nan) values are 0.

    Parameters
    ----------
    remsoil_index : np.ndarray
        Vegetation index with possible NA values with shape (rows, cols).

    Returns
    -------
    np.ndarray
        Binary mask (rows, cols) with 1 for valid, 0 for invalid.
    """
    binary_mask = np.zeros_like(remsoil_index, dtype=np.uint8)
    valid_mask = ~np.isnan(remsoil_index)
    binary_mask[valid_mask] = 1
    return binary_mask

def process_and_save_masks_multi(
                           input_images: str,
                           masks_out: str,
                           masked_images_out: str,
                           threshold: float = 0.1,
                           index_name: str = 'ExG',
                           comparison: str = 'less',
                           custom_indices: dict = None,
                           output_type: str = '.tif',
                           band_map: dict = None):
    """
    1) Loads images from 'input_images' folder via rasterio.
    2) Computes vegetation index for each pixel within each image using rem_soil_index_multi.
    3) Creates a binary mask and saves it to 'masks_out' folder.
    4) Applies the mask to the image and saves the masked result to 'masked_images_out'.

    Parameters
    ----------
    input_images : str [Required]
        Path to folder containing images (e.g., .tif, .tiff).
    masks_out : str [Required]
        Path to folder where binary masks will be saved.
    masked_images_out : str [Required]
        Path to folder where masked images will be saved.
    threshold : float [Optional]
        Threshold above or below which a pixel is marked for filtering. Default is 0.1 (ExG).
    index_name : str [Optional]
        Name of the index to use from custom_indices. Default 'ExG'.
    comparison : str [Optional]
        'less' or 'greater'. Defaults to 'less' according to the ExG index.
    custom_indices : dict [Optional]
        Optional dict of lambdas {index_name: lambda R,G,B: ...}. If None, will use
        a default dict with only ExG.
    output_type : str [Optional]
        File extension for output images (e.g., '.jpg', '.tif'). Defaults to '.tif'.
    band_map : dict [Optional]
        Dictionary mapping channel names to integer band indices.
        e.g., {'B':0, 'G':1, 'R':2, 'RE':3, 'NIR':4}.
        
    Returns
    -------
    None
    """

    # Ensure output directories exist
    os.makedirs(masks_out, exist_ok=True)
    os.makedirs(masked_images_out, exist_ok=True)

    # Load images with rasterio-based loader
    images, img_files = load_multiband_images_from_folder(input_images)

    # Process each image
    for i, img in enumerate(images):
        # Calculate index
        remsoil = rem_soil_index_multi(
            img, 
            threshold=threshold, 
            index_name=index_name, 
            comparison=comparison,
            custom_indices=custom_indices,
            band_map=band_map
        )

        # Create binary mask
        binary_mask = create_binary_mask_multi(remsoil)

        # Save the mask to a file
        mask_output_path = os.path.join(masks_out, f"{img_files[i]}_mask{output_type}")
        
        # Save the mask as a typical 8-bit grayscale image using OpenCV:
        cv2.imwrite(mask_output_path, (binary_mask * 255).astype(np.uint8))

        # Apply mask to each band of the original image
        # img shape: (bands, rows, cols), or CHW format (channels, height, width)
        masked_img = img.copy().astype(np.float32)
        # broadcast multiply: for each band, multiply by the 0/1 mask
        # make sure shapes are compatible for broadcasting
        # masked_img shape is (bands, rows, cols), binary_mask is (rows, cols)
        masked_img *= binary_mask[np.newaxis, :, :]
        
        # Save the image in either .jpg/.png/.jpeg OR in .tif format
        # The type the user enters will determine how the image is handled before
        # it is exported.

        if masked_img.shape[0] >= 3 and output_type.lower() in ['.jpg', '.png', '.jpeg']:
            # Reorder to (rows, cols, channels) and cast to 8-bit
            transposed_img = np.transpose(masked_img[:3], (1, 2, 0))  # shape (rows, cols, 3)
            # Scale or clip if necessary
            masked_img_8u = np.clip(transposed_img, 0, 255).astype(np.uint8)
            
            masked_img_output_path = os.path.join(
                masked_images_out, f"{img_files[i]}_masked{output_type}"
            )
            cv2.imwrite(masked_img_output_path, masked_img_8u)
        
        else:
            # Save multi-band data as a GeoTIFF
            masked_img_output_path = os.path.join(
                masked_images_out, f"{img_files[i]}_masked.tif"
            )
            # Simple, multi-band TIF without georeferencing:
            bands, rows, cols = masked_img.shape
            # Writes a generic 32-bit float tiff
            with rasterio.open(
                masked_img_output_path,
                'w',
                driver='GTiff',
                height=rows,
                width=cols,
                count=bands,
                dtype=masked_img.dtype
            ) as dst:
                dst.write(masked_img)
        
        print(f"Processed and saved mask for {img_files[i]} at {mask_output_path}")
        

#### Example usage ####
input_images = 'C:/Users/aaron.desalvio/Documents/TESTING/Python_Phenomics/MULTI/Cropped'
masks_out = 'C:/Users/aaron.desalvio/Documents/TESTING/Python_Phenomics/MULTI/Masks'
masked_images_out = 'C:/Users/aaron.desalvio/Documents/TESTING/Python_Phenomics/MULTI/Masked_Images'
my_band_map = {'B':0, 'G':1, 'R':2, 'RE':3, 'NIR':4}

# Use the default masking vegetation index (ExG)
process_and_save_masks_multi(
                           input_images=input_images,
                           masks_out=masks_out,
                           masked_images_out=masked_images_out,
                           threshold=0.1,
                           index_name='ExG',
                           comparison='less',
                           custom_indices=None,
                           output_type='.tif',
                           band_map=my_band_map)

# Use a custom masking vegetation index (HUE)
my_index = {'HUE': lambda R,G,B: np.arctan(2*(B-G-R)/30.5*(G-R))}
process_and_save_masks_multi(
                           input_images=input_images,
                           masks_out=masks_out,
                           masked_images_out=masked_images_out,
                           threshold=0.7,
                           index_name='HUE',
                           comparison='greater',
                           custom_indices=my_index,
                           output_type='.tif',
                           band_map=my_band_map)







