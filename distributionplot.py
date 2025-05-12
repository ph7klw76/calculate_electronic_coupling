import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
import matplotlib.gridspec as gridspec
from sklearn.cluster import KMeans
import glob, os

# ─────────────────────────────────────────────────────────────
# 1. USER SETTINGS  (edit these paths / parameters)
# ─────────────────────────────────────────────────────────────
folder_path = "frame"   # <‑‑ the folder with many .gro files
atom_filter = {"O1","C10","C5","C6","C1","C2","C3","C4",
               "O2","C7","C13","C12","C11","C8","C9","C14"}
bins_2d     = 80          # histogram resolution
k_clusters  = 3           # K‑means clusters
out_txt     = "frame/aggregated_nearest_neighbours_with_clusters.txt"
# ─────────────────────────────────────────────────────────────

def centroid(xyz): return xyz.mean(axis=0)

def plane_normal(xyz):
    centred = xyz - centroid(xyz)
    w, v = np.linalg.eigh(centred.T @ centred)
    n = v[:,0]                                # eigen‑vector with smallest eigen‑value
    return n / np.linalg.norm(n)

def parse_filtered_gro(path, wanted_atoms):
    """Return (residue_ids, list-of‑xyz‑arrays, box-vector) for a single .gro"""
    with open(path) as fh:
        fh.readline()                        # title
        nat = int(fh.readline())
        coords, res_map = [], {}
        for _ in range(nat):
            ln = fh.readline()
            resid = int(ln[0:5])
            aname = ln[10:15].strip()
            if aname not in wanted_atoms:     # skip anything not in filter
                continue
            xyz = tuple(float(ln[s:e]) for s,e in ((20,28),(28,36),(36,44)))
            res_map.setdefault(resid, []).append(xyz)
        box = np.fromstring(fh.readline(), sep=" ")
    # keep only residues that had ≥1 filtered atom
    resid_list = sorted(res_map)
    xyz_list   = [np.asarray(res_map[r]) for r in resid_list]
    return resid_list, xyz_list, box

# ── pass 1: collect every distance / angle from every frame ──
all_d, all_a, frame_field = [], [], []     # aggregated pools
pair_rows = []                             # to save full neighbour list

gro_files = sorted(glob.glob(os.path.join(folder_path, "*.gro")))
if not gro_files:
    raise RuntimeError(f"No .gro files found in: {folder_path}")

for f_idx, grof in enumerate(gro_files):
    res_id, xyz_sets, box = parse_filtered_gro(grof, atom_filter)

    centres = np.array([centroid(p)      for p in xyz_sets])
    normals = np.array([plane_normal(p)  for p in xyz_sets])

    # minimum‑image distances
    dv   = centres[:,None,:] - centres[None,:,:]
    dv  -= box * np.round(dv/box)
    dist = np.linalg.norm(dv, axis= -1)
    np.fill_diagonal(dist, np.inf)
    nidx = np.argmin(dist, axis=1)
    nd   = dist[np.arange(len(res_id)), nidx]

    ang  = np.degrees(np.arccos(
             np.clip(np.abs(np.einsum("ij,ij->i", normals, normals[nidx])), 0, 1)))

    # book‑keeping
    all_d.extend(nd)
    all_a.extend(ang)
    frame_field.extend([f_idx]*len(nd))
    for i,(j,d,a) in enumerate(zip([res_id[k] for k in nidx], nd, ang)):
        pair_rows.append((f_idx, res_id[i], j, d, a))

all_d = np.asarray(all_d)
all_a = np.asarray(all_a)
data  = np.column_stack((all_d, all_a))

# ── clustering on aggregated data ────────────────────────────
zs   = (data - data.mean(0)) / data.std(0)
km   = KMeans(n_clusters=k_clusters, n_init=10, random_state=0).fit(zs)
labs = km.labels_
cent = km.cluster_centers_ * data.std(0) + data.mean(0)

# append cluster to records & write txt
with open(out_txt, "w") as fh:
    fh.write("Frame\tMol_i\tMol_j_nearest\tDistance_nm\tAngle_deg\tCluster\n")
    for (fr,mi,mj,dg,ag),cl in zip(pair_rows, labs):
        fh.write(f"{fr}\t{mi}\t{mj}\t{dg:.6f}\t{ag:.3f}\t{cl}\n")
print(f"Saved neighbour list → {out_txt}")

# ── 2‑D histogram + marginals, à‑la Figure ───────────────────
# Assuming necessary imports (numpy as np, matplotlib.pyplot as plt, etc.)
# Assuming `all_d`, `all_a`, `data`, `bins_2d` are defined

# Define bins based on data range (or desired range)
# Make sure ybins starts exactly where you want the plot to start (e.g., 0)
xbins = np.linspace(all_d.min(), all_d.max(), bins_2d + 1) # +1 for edges
ybins = np.linspace(0, 90, bins_2d + 1)                   # +1 for edges
H, xe, ye = np.histogram2d(all_d, all_a, bins=[xbins, ybins])

# ellipse stats
μ   = data.mean(0)
Σ   = np.cov(data.T)
w, v = np.linalg.eigh(Σ); idx = w.argsort()[::-1]
w, v = w[idx], v[:,idx]

def add_sigma(ax, s, **kw):
    width, height = 2*s*np.sqrt(w)
    ang = np.degrees(np.arctan2(*v[:,0][::-1]))
    # Ensure ellipse center μ is within plot limits if needed, though usually okay
    ax.add_patch(Ellipse(μ, width, height, angle=ang, fill=False, **kw)) # Use 'angle' keyword

fig = plt.figure(figsize=(8,8), constrained_layout=True)
gs  = gridspec.GridSpec(2,2, width_ratios=[4,1.3], height_ratios=[1.3,4],
                        wspace=0.05, hspace=0.05) # Note: constrained_layout often handles spacing well

axH = fig.add_subplot(gs[1,0])
axX = fig.add_subplot(gs[0,0], sharex=axH)
axY = fig.add_subplot(gs[1,1], sharey=axH)

# Plot the 2D histogram
# Use shading='flat' or 'gouraud' if pcolormesh gives warnings about dimensions
# For histogram data where xe/ye are bin edges, 'flat' is appropriate.
# H.T is needed because histogram2d returns (xbins, ybins) but pcolormesh expects (y, x) shapes for C
pc  = axH.pcolormesh(xe, ye, H.T, cmap="viridis", shading='flat')
cbar = fig.colorbar(pc, ax=axH, # Associate with main axes
                    orientation='horizontal', # Horizontal layout
                    location='bottom',       # Position below the axes
                    shrink=0.6,             # Fraction of original width
                    aspect=30,              # Ratio of long to short dimension
                    pad=0.18)               # Padding below axH (increase if it overlaps xlabel)

# --- FIX: Set axis limits to match the data extent ---
axH.set_xlim(xe.min(), xe.max())
axH.set_ylim(ye.min(), ye.max())
cbar.set_label("Frequency")
# --- End Fix ---


# Add ellipses and text
for s,ls in zip((1,2,3),("--",":","-.")):
    add_sigma(axH, s, color="white", lw=1, ls=ls)
# Check if μ is within the calculated limits before plotting text
if xe.min() <= μ[0] <= xe.max() and ye.min() <= μ[1] <= ye.max():
    axH.text(*μ, "μ", color="white", ha="center", va="center", fontsize=9)

# Set labels for the main plot
axH.set_xlabel("Nearest‑neighbour distance (nm)")
axH.set_ylabel("Angle between plane normals (°)")
axY.set_xlabel("Frequency")
# Plot marginal histograms
axX.hist(all_d, bins=xbins, color="steelblue")


axY.hist(all_a, bins=ybins, orientation="horizontal", color="firebrick")
axX.set_ylabel("Frequency")
# Optional: Remove ticks from shared axes if they appear due to sharing before axis('off')
plt.setp(axX.get_xticklabels(), visible=False)
plt.setp(axY.get_yticklabels(), visible=False)


plt.show()
fig.savefig("2Dplot.tiff", dpi=600, bbox_inches='tight')
# quick cluster summary
summary = pd.DataFrame({
    "Cluster": range(k_clusters),
    "Centre_dist_nm": cent[:,0],
    "Centre_angle_deg": cent[:,1],
    "Members": np.bincount(labs)
})
print(summary)
