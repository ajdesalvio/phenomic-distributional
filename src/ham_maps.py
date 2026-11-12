"""Build the eleven Figure 5 maps from explicit inputs; no model fitting."""
import argparse
import csv
import json
from pathlib import Path
import sys
import warnings

import numpy as np
from PIL import Image

sys.path.insert(0, str(Path(__file__).resolve().parent))
import ham
HELPER = Path(ham.__file__)

def strict_curve(rows, day, vi, day_column, env):
    p, h = ham.lookup_curve(rows, day, vi, day_column, 'H2Individual', env)
    q = np.rint(100 * p).astype(int)
    if len(q) != 100 or not np.array_equal(q, np.arange(1, 101)) or not np.allclose(p, q / 100, atol=1e-12, rtol=0):
        raise ValueError(f'Incomplete/noncanonical 100-quantile curve: {env}/{day}/{vi}')
    return p, h


def jet_lut():
    # Classic jet's channel control points; 256 colors, common linear [0,1] scale.
    controls = [([0,.35,.66,.89,1],[0,0,1,1,.5]),
                ([0,.125,.375,.64,.91,1],[0,0,1,1,0,0]),
                ([0,.11,.34,.65,1],[.5,1,1,0,0])]
    return np.round(255*np.stack([np.interp(np.linspace(0,1,256),x,y) for x,y in controls],axis=1)).astype(np.uint8)


def write_csv(path, rows):
    with path.open('w', newline='', encoding='utf-8') as fh:
        writer=csv.DictWriter(fh,fieldnames=list(rows[0])); writer.writeheader(); writer.writerows(rows)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output-dir',type=Path,required=True)
    parser.add_argument('--cotton-images',type=Path,required=True)
    parser.add_argument('--maize-images',type=Path,required=True)
    parser.add_argument('--cotton-h2',type=Path,required=True)
    parser.add_argument('--maize-h2',type=Path,required=True)
    parser.add_argument('--cotton-dates',type=Path,required=True)
    parser.add_argument('--maize-dates',type=Path,required=True)
    parser.add_argument('--image-manifest',type=Path,required=True)
    args=parser.parse_args()
    import rasterio
    if args.output_dir.exists(): raise FileExistsError('Use a new output directory; prior exports are never overwritten.')
    cotton_conversion=args.cotton_dates
    maize_conversion=args.maize_dates
    filename_manifest=args.image_manifest
    configs=[dict(species='cotton',vi='PSRI',env=None,day_column='DAT',days=[57,86,96,113,130],plot='337_03',
                  images=args.cotton_images,h2=args.cotton_h2,conversion=cotton_conversion,
                  band_map={'B':0,'G':1,'R':2,'RE':3,'NIR':4},vi_scale=65535,display_scale=3276,
                  template='{date}_337_03_masked.tif',rotation=0),
             dict(species='maize',vi='ExGR',env='TXH1.2020',day_column='DAP',days=[28,47,68,88,96,111],plot='CS20-G2F1-031',
                  images=args.maize_images,h2=args.maize_h2,conversion=maize_conversion,
                  band_map={'R':0,'G':1,'B':2},vi_scale=255,display_scale=255,
                  template='{date}-CS20-G2F1-031_masked.tif',rotation=1)]
    input_paths={HELPER,Path(__file__).resolve(),args.cotton_h2,args.maize_h2,cotton_conversion,maize_conversion,filename_manifest}
    memberships={}
    for r in ham.read_rows(filename_manifest): memberships.setdefault(r['File.Name.TIF.Unified.Masked'],set()).add(r['Env'])
    jobs=[]
    for c in configs:
        # Retain only requested keys; do not load the broad maize table into memory.
        with c['h2'].open(newline='',encoding='utf-8-sig') as fh:
            curves=[r for r in csv.DictReader(fh) if r['Vegetation.Index']==c['vi'] and float(r[c['day_column']]) in c['days'] and (c['env'] is None or r['Env']==c['env'])]
        conversions=ham.read_rows(c['conversion'])
        for day in c['days']:
            p,h=strict_curve(curves,day,c['vi'],c['day_column'],c['env'])
            selected=[r for r in curves if float(r[c['day_column']])==day]
            if c['species']=='cotton':
                if any(r['Converged']!='TRUE' or r['Warnings'] or (r['ConvergenceMessages'] and not (r['Singular']=='TRUE' and r['ConvergenceMessages']=="boundary (singular) fit: see help('isSingular')")) for r in selected):
                    raise ValueError('Cotton curve contains a failed or unresolved-warning fit.')
                if any(r['Formula']!='response ~ 1 + (1 | Pedigree) + (1 | Block)' for r in selected):
                    raise ValueError('Cotton curve is not the selected block-only model.')
                vg=np.array([float(r['Vg']) for r in selected]); ve=np.array([float(r['Ve']) for r in selected])
                hh=np.array([float(r['H2Individual']) for r in selected])
                if not np.allclose(hh,vg/(vg+ve),atol=1e-12,rtol=0): raise ValueError('Incorrect cotton H2 arithmetic.')
            flight=ham.flight_for_day(conversions,day,c['day_column'],c['env'])
            path=c['images']/c['template'].format(date=flight)
            if not path.is_file(): raise FileNotFoundError(path)
            if c['env'] and memberships.get(path.name)!= {c['env']}:
                raise ValueError(f'Ambiguous/wrong environment membership: {path}')
            input_paths.add(path); jobs.append((c,day,flight,path,p,h))
    before={str(p):ham.sha256_file(p) for p in sorted(input_paths,key=str)}
    args.output_dir.mkdir(parents=True)
    assets=args.output_dir/'Map_Assets'; assets.mkdir()
    lut=jet_lut(); Image.fromarray(np.repeat(lut[::-1,None,:],20,axis=1)).save(assets/'Colorbar_Jet_V1.png')
    records=[]; curve_rows=[]; metadata=[]
    for c,day,flight,path,p,h in jobs:
        with warnings.catch_warnings():
            warnings.filterwarnings('ignore',category=rasterio.errors.NotGeoreferencedWarning)
            with rasterio.open(path) as src:
                img=src.read(); masks=src.read_masks()
                raster_info=dict(band_descriptions=src.descriptions,color_interpretations=[x.name for x in src.colorinterp],nodata=src.nodata,dtypes=src.dtypes)
        if img.shape[0]!=len(c['band_map']): raise ValueError('Unexpected band count.')
        if c['species']=='maize' and raster_info['color_interpretations']!=['red','green','blue']:
            raise ValueError('Maize TIFF RGB channel metadata mismatch.')
        vi,valid,counts=ham.valid_vi_image(img,c['band_map'],lambda bands:ham.VI_FUNCTIONS[c['vi']](bands,c['vi_scale']),masks)
        mapped,ranks,rank_info=ham.map_quantile_heritability(vi,p,h)
        rgb=ham.rgb_composite(img,c['band_map'],c['display_scale']); rgb[~valid]=0
        heat=np.zeros((*vi.shape,3),dtype=np.uint8)
        heat[valid]=lut[np.minimum(255,(256*mapped[valid]).astype(int))]
        rgb8=np.round(255*rgb).astype(np.uint8)
        stem=f"{c['species']}_{c['day_column']}{day}_{c['vi']}_V1"
        # Rotation for display only; NPZ numeric arrays remain in native TIFF orientation.
        Image.fromarray(np.rot90(rgb8,c['rotation'])).save(assets/f'{stem}_RGB.png')
        Image.fromarray(np.rot90(heat,c['rotation'])).save(assets/f'{stem}_HAM.png')
        rank_array=np.full(vi.shape,np.nan); rank_array[valid]=ranks
        np.savez_compressed(assets/f'{stem}_Numeric.npz',vi=vi,rank=rank_array,mapped_h2=mapped,valid_mask=valid,probability=p,h2_curve=h)
        record=dict(Species=c['species'],Env=c['env'] or 'CS24_Cotton',VI=c['vi'],DayColumn=c['day_column'],Day=day,FlightDate=flight,PlotID=c['plot'],
                    ValidPixels=int(valid.sum()),TotalPixels=int(valid.size),NonfiniteVIExcluded=counts['nonfinite_vi_pixels_excluded'],
                    UniqueVIValues=int(np.unique(vi[valid]).size),BelowProbability001=int((ranks<.01).sum()),
                    MedianTraitH2=float(h[49]),MinimumCurveH2=float(h.min()),MaximumCurveH2=float(h.max()),
                    MeanMappedH2=float(mapped[valid].mean()),MinimumMappedH2=float(mapped[valid].min()),MaximumMappedH2=float(mapped[valid].max()),
                    RGBFile=f'Map_Assets/{stem}_RGB.png',HAMFile=f'Map_Assets/{stem}_HAM.png',NumericFile=f'Map_Assets/{stem}_Numeric.npz')
        records.append(record)
        curve_rows.extend(dict(Species=c['species'],Env=c['env'] or 'CS24_Cotton',VI=c['vi'],DayColumn=c['day_column'],Day=day,Probability=float(pp),H2Individual=float(hh)) for pp,hh in zip(p,h))
        metadata.append(dict(**record,SourceImage=str(path),SourceH2=str(c['h2']),BandMap=c['band_map'],VIScale=c['vi_scale'],DisplayScale=c['display_scale'],DisplayRotationCCW=90*c['rotation'],RasterMetadata=raster_info,PixelCounts=counts,RankInfo=rank_info))
    after={p:ham.sha256_file(p) for p in before}
    if before!=after: raise RuntimeError('Inputs changed during map generation.')
    write_csv(args.output_dir/'Figure5_Map_Summary_V1.csv',records)
    write_csv(args.output_dir/'Figure5_H2_Curves_V1.csv',curve_rows)
    provenance=dict(Status='COMPLETE_MAP_DATA',MapCount=len(records),CurveRows=len(curve_rows),RankReference='finite valid pixels of each focal image only',
                    H2Definition='Vg/(Vg+Ve), conditional on the fitted nuisance design terms; no replication divisor',
                    Interpretation='Quantile-level H2 displayed by pixel rank, not pixel-specific heritability or individual genetic merit.',
                    Grid='exactly p=.01,.02,...,1.00; linear interpolation; nearest endpoint below .01',
                    ColorScale='classic jet 256 colors, linear 0 to 1; invalid pixels black',
                    CottonBandEvidence='MULTI_VIs_Distributional_Multiprocessing.py lines238 and320: B,G,R,RE,NIR. TIFF band descriptions are absent.',
                    SourcesSHA256=before,InputHashesUnchanged=True,NumpyVersion=np.__version__,RasterioVersion=rasterio.__version__,Maps=metadata)
    provenance['OutputsSHA256']={str(p.relative_to(args.output_dir)):ham.sha256_file(p) for p in args.output_dir.rglob('*') if p.is_file()}
    (args.output_dir/'Figure5_Map_Provenance_V1.json').write_text(json.dumps(provenance,indent=2),encoding='utf-8')
    print(json.dumps(dict(Status=provenance['Status'],MapCount=len(records),Output=str(args.output_dir))))


if __name__=='__main__': main()
