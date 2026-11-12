import pandas as pd
import plotly.graph_objects as go
from plotly.colors import sample_colorscale
import os
import numpy as np
from scipy.interpolate import UnivariateSpline


# Define the data path and file name
data_path = 'path/to/python/files/'
file_name = 'CS24-STEL-MULTI-RTK-Distributions-NDVI-FieldCoords.csv'
file_path = data_path + file_name

# Path for saving output files
output_path = 'path/for/saving/many/image/files'
os.makedirs(output_path, exist_ok=True)

# Read the CSV data into a pandas DataFrame
df = pd.read_csv(file_path)

# Sort the DataFrame by 'DAT'
df = df.sort_values(by='DAT')

# Define filtering criteria
ped_num = 5016
DATs_exclude = [46, 71, 72, 99, 100]  # equivalent to c(46, 71:72, 99:100) in R

# Filter the DataFrame: select rows where Pedigree equals ped_num,
# DAT is not in the exclusion list, and Replicate equals 7
df_subset = df[(df['Pedigree'] == ped_num) &
               (~df['DAT'].isin(DATs_exclude)) &
               (df['Replicate'] == 7)]

# Define a custom colorscale for NDVI quantile
colorscale = [
    [0.0, 'brown'],
    [0.35, 'sandybrown'],
    [0.7, 'lightgreen'],
    [0.85, 'green'],
    [1.0, 'darkgreen']
]

fig = go.Figure()
# 3-D quantile ribbons, one trace per DAT
show_scale = True
for dat_val, group in df_subset.groupby('DAT'):
    group = group.sort_values(by='probability')
    fig.add_trace(
        go.Scatter3d(
            x=group['probability'],
            y=group['DAT'],
            z=group['quantile'],
            mode='lines',
            line=dict(
                color=group['quantile'],
                colorscale=colorscale,
                cmin=0,
                cmax=1,
                width=4,
                showscale=show_scale,
                colorbar=dict(
                    thickness=20,
                    len=0.55,
                    x=1.1,
                    y=0.5,
                    title="NDVI Value",
                    tickfont=dict(size=16)
                )
            ),
            showlegend=False,
            name=str(dat_val)
        )
    )
    show_scale = False

# Black line through the medians (p = 0.5)
median_df = (
    df_subset[np.isclose(df_subset['probability'], 0.5)]
    .sort_values('DAT')
)

dat_vals = median_df['DAT'].values
ndvi_vals = median_df['quantile'].values
spline = UnivariateSpline(dat_vals, ndvi_vals, s=0.5 * len(dat_vals))
dat_dense = np.linspace(dat_vals.min(), dat_vals.max(), 200)
ndvi_smooth = spline(dat_dense)

fig.add_trace(go.Scatter3d(
    x=[0.5] * len(dat_dense),
    y=dat_dense,
    z=ndvi_smooth, # smoothed NDVI
    mode='lines',
    line=dict(color='black', width=6),
    name='Smoothed median NDVI',
    showlegend=True
))

# Layout tweaks
fig.update_layout(
    scene=dict(
        xaxis=dict(
            title=dict(text='Quantile Level', font=dict(size=18, family='Arial, bold')),
            tickfont=dict(size=14)
        ),
        yaxis=dict(
            title=dict(text='DAT', font=dict(size=18, family='Arial, bold')),
            tickfont=dict(size=14)
        ),
        zaxis=dict(
            title=dict(text='NDVI Value', font=dict(size=18, family='Arial, bold')),
            tickfont=dict(size=14)
        ),
        camera=dict(eye=dict(x=2, y=-1.5, z=0.45))
    )
)

# Show / save
fig.show()
fig.write_image(
    os.path.join(output_path, "Example.jpg"),
    width=1200,
    height=1200,
    scale=5
)

# Define number of frames for the full rotation
n_frames = 600

# Loop to capture frames at different camera angles
for i in range(n_frames):
    angle = i * 360 / n_frames  # angle in degrees for a full rotation
    # Compute new camera eye position for a smooth rotation around the z-axis
    eye_x = 2 * np.cos(np.deg2rad(angle))
    eye_y = 2 * np.sin(np.deg2rad(angle))
    
    # Update camera position; adjust z as needed for your view
    fig.update_layout(scene=dict(camera=dict(eye=dict(x=eye_x, y=eye_y, z=1.45))))
    
    # Save the current frame as a PNG file
    frame_filename = os.path.join(output_path, f"frame_{i:03d}.png")
    fig.write_image(frame_filename, width=1200, height=800, scale=4)
    

# 2D plot with all DATs shown in the same XY plane

# Use RGB strings (no named colors) so sample_colorscale can interpolate
colorscale_rgb = [
    [0.00, 'rgb(165,42,42)'],   # brown
    [0.35, 'rgb(244,164,96)'],  # sandybrown
    [0.70, 'rgb(144,238,144)'], # lightgreen
    [0.85, 'rgb(0,128,0)'],     # green
    [1.00, 'rgb(0,100,0)']      # darkgreen
]

# Robust color normalization to your data range
vmin = float(df_subset['quantile'].min())
vmax = float(df_subset['quantile'].max())
rng  = vmax - vmin if vmax > vmin else 1.0  # avoid div-by-zero

def norm01(v):
    return float(np.clip((v - vmin) / rng, 0.0, 1.0))

fig2 = go.Figure()

# Thin if needed for performance (higher -> fewer segments)
step = 1

for dat_val, g in df_subset.groupby('DAT'):
    g = g.sort_values('probability').iloc[::step]
    p = g['probability'].to_numpy()
    q = g['quantile'].to_numpy()  # NDVI value

    # Draw short segments with NDVI-based color
    for i in range(len(p) - 1):
        ndvi_mid = 0.5 * (q[i] + q[i+1])
        col = sample_colorscale(colorscale_rgb, norm01(ndvi_mid))[0]  # returns 'rgb(r,g,b)'

        fig2.add_trace(
            go.Scatter(
                x=[p[i], p[i+1]],
                y=[q[i], q[i+1]],
                mode='lines',
                line=dict(color=col, width=3),
                hovertemplate=f"DAT: {dat_val}<br>p: %{{x:.2f}}<br>NDVI: %{{y:.3f}}<extra></extra>",
                showlegend=False
            )
        )

# Add a proper NDVI colorbar via an invisible marker trace
fig2.add_trace(
    go.Scatter(
        x=[None, None],
        y=[None, None],
        mode='markers',
        marker=dict(
            colorscale=colorscale_rgb,
            cmin=vmin, cmax=vmax,
            color=[vmin, vmax],  # span the scale
            size=0,
            showscale=True,
            colorbar=dict(
                title="NDVI Value",
                thickness=20,
                len=0.85,
                y=0.5,
                tickfont=dict(size=15)
            ),
        ),
        hoverinfo='skip',
        showlegend=False
    )
)

# Optional guide at p = 0.5
fig2.add_vline(x=0.5, line_width=1, line_dash='dot', line_color='black')

fig2.update_layout(
    xaxis=dict(
        title=dict(text='Quantile Level', font=dict(size=19, family='Arial, bold')),
        tickfont=dict(size=15),
        range=[0, 1]
    ),
    yaxis=dict(
        title=dict(text='NDVI Value', font=dict(size=19, family='Arial, bold')),
        tickfont=dict(size=15),
        range=[vmin, vmax]  # or [0, 1] if NDVI is normalized
    ),
    template='simple_white',
    margin=dict(l=70, r=60, t=40, b=60),
)

fig2.show()

# Save
fig2.write_image(
    os.path.join(output_path, "Quantiles_2D_AllDATs_coloredByNDVI_V4.jpg"),
    width=600, height=450, scale=8
)
