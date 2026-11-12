from pathlib import Path
import csv
import hashlib
import re
import numpy as np


VI_FUNCTIONS = {'PSRI': lambda bands, scale: (bands['R'] / scale - bands['G'] / scale) / (bands['RE'] / scale), 'ExGR': lambda bands, scale: 3 * (bands['G'] / scale) - 2.4 * (bands['R'] / scale) - bands['B'] / scale}

def sha256_file(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()

def read_rows(path):
    with Path(path).open(newline='', encoding='utf-8-sig') as handle:
        return list(csv.DictReader(handle))

def lookup_curve(rows, day, vi, day_column, value_column, env=None):
    """Return one unique finite probability/H2 curve; isolate environment first."""
    if not rows:
        raise ValueError('Heritability table is empty.')
    required = {day_column, 'Vegetation.Index', 'Probability', value_column}
    missing = required - rows[0].keys()
    if missing:
        raise ValueError(f'Missing heritability columns: {sorted(missing)}')
    if env is not None:
        if 'Env' not in rows[0]:
            raise ValueError('An environment was supplied but the H2 table has no Env column.')
        rows = [row for row in rows if row['Env'] == env]
    elif 'Env' in rows[0] and len({row['Env'] for row in rows}) > 1:
        raise ValueError('Multiple environments are present: supply --env.')
    rows = [row for row in rows if float(row[day_column]) == float(day)
            and row['Vegetation.Index'] == vi]
    if rows and 'grp' in rows[0]:
        rows = [row for row in rows if row['grp'] in ('Pedigree', 'Genotype')]
    if not rows:
        raise ValueError(f'No genotype H2 curve for {env=}, {day=}, {vi=}.')
    pairs = sorted((float(row['Probability']), float(row[value_column])) for row in rows)
    values = np.asarray(pairs, dtype=float)
    if not np.isfinite(values).all():
        raise ValueError('H2 curve contains missing or nonfinite values; audit the affected fits.')
    p, h = values[:, 0], values[:, 1]
    if len(p) < 2 or np.any(np.diff(p) <= 0):
        raise ValueError('Expected at least two UNIQUE probability levels for one environment/day/VI.')
    if np.any((p < 0) | (p > 1)) or np.any((h < 0) | (h > 1)):
        raise ValueError('Probability and H2 values must be between zero and one; no clipping is applied.')
    return p, h

def valid_vi_image(img, band_map, vi_func, raster_masks=None):
    """Retain nonblack, nodata-valid, finite band/VI pixels consistently."""
    img = np.asarray(img)
    if img.ndim != 3:
        raise ValueError('Image must have dimensions (bands, rows, columns).')
    if not {'R', 'G', 'B'}.issubset(band_map):
        raise ValueError('band_map must define R, G, and B for display.')
    indices = list(band_map.values())
    if any(not isinstance(i, int) or i < 0 or i >= img.shape[0] for i in indices):
        raise ValueError('Band index is outside the image.')
    selected = img[indices].astype(float)
    valid = np.any(selected != 0, axis=0) & np.all(np.isfinite(selected), axis=0)
    if raster_masks is not None:
        masks = np.asarray(raster_masks)
        if masks.shape != img.shape:
            raise ValueError('Raster mask shape does not match image shape.')
        valid &= np.all(masks[indices] > 0, axis=0)
    candidates = int(valid.sum())
    out = np.full(img.shape[1:], np.nan, dtype=float)
    if candidates:
        bands = {name: img[index].astype(float)[valid] for name, index in band_map.items()}
        with np.errstate(divide='ignore', invalid='ignore', over='ignore'):
            vi_values = np.asarray(vi_func(bands), dtype=float)
        if vi_values.ndim == 0:
            vi_values = np.full(candidates, vi_values.item())
        if vi_values.shape != (candidates,):
            raise ValueError('VI function did not return one value per candidate pixel.')
        out[valid] = vi_values
    valid &= np.isfinite(out)
    out[~valid] = np.nan
    if not valid.any():
        raise ValueError('Image has no valid finite VI pixels.')
    return out, valid, {'candidate_pixels': candidates, 'valid_vi_pixels': int(valid.sum()),
                        'nonfinite_vi_pixels_excluded': candidates - int(valid.sum())}

def rgb_composite(img, band_map, display_scale):
    if not np.isfinite(display_scale) or display_scale <= 0:
        raise ValueError('display_scale must be finite and positive.')
    rgb = np.stack([img[band_map[name]] for name in ('R', 'G', 'B')], axis=-1).astype(float)
    return np.clip(np.nan_to_num(rgb / display_scale, nan=0, posinf=1, neginf=0), 0, 1)

def map_quantile_heritability(vi_image, probability, heritability, reference_values=None):
    """Use right-continuous eCDF F(x)=#{v<=x}/n, then interpolate H2(F(x)).

The recorded curve can start at p=.01. Lower ranks use the nearest endpoint,
as in np.interp; endpoint clipping counts are included in metadata.
"""
    valid = np.isfinite(vi_image)
    focal = np.asarray(vi_image)[valid]
    reference = focal if reference_values is None else np.asarray(reference_values, dtype=float).ravel()
    reference = reference[np.isfinite(reference)]
    if not focal.size or not reference.size:
        raise ValueError('Focal and reference distributions must contain finite values.')
    p = np.asarray(probability, dtype=float)
    h = np.asarray(heritability, dtype=float)
    if p.ndim != 1 or p.shape != h.shape or len(p) < 2 or not np.isfinite(p).all() or not np.isfinite(h).all():
        raise ValueError('A finite matching probability/H2 curve is required.')
    if np.any(np.diff(p) <= 0) or np.any((p < 0) | (p > 1)) or np.any((h < 0) | (h > 1)):
        raise ValueError('Curve probabilities must be unique/increasing and probability/H2 within [0,1].')
    ranks = np.searchsorted(np.sort(reference), focal, side='right') / float(reference.size)
    output = np.full(np.asarray(vi_image).shape, np.nan, dtype=float)
    output[valid] = np.interp(ranks, p, h)
    info = {'reference_finite_pixel_count': int(reference.size),
            'pixels_below_curve_min_probability': int((ranks < p[0]).sum()),
            'pixels_above_curve_max_probability': int((ranks > p[-1]).sum()),
            'ecdf_ties': 'right: fraction of finite valid pixels <= focal value',
            'interpolation': 'linear; endpoint hold outside tabulated probability range'}
    return output, ranks, info

def flight_for_day(rows, day, day_column, env):
    if not rows or day_column not in rows[0] or 'Flight_Date_YYYYMMDD' not in rows[0]:
        raise ValueError('Date conversion table lacks required columns.')
    if env is not None:
        if 'Env' not in rows[0]:
            raise ValueError('Date conversion table must contain Env for maize.')
        rows = [row for row in rows if row['Env'] == env]
    elif 'Env' in rows[0] and len({row['Env'] for row in rows}) > 1:
        raise ValueError('Date conversion table has multiple environments; supply --env.')
    flights = {row['Flight_Date_YYYYMMDD'].removesuffix('.0') for row in rows
               if float(row[day_column]) == float(day)}
    if len(flights) != 1:
        raise ValueError(f'Expected exactly one calendar date for {env=}, {day=}; found {flights}.')
    flight = flights.pop()
    if not re.fullmatch(r'\d{8}', flight):
        raise ValueError(f'Unexpected calendar date: {flight}')
    return flight
