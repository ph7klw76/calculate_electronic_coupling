import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.cluster import KMeans

# ---------- re‑use parsing from previous cell ----------
file_path = "output_whole.gro"
atom_filter = {"O1","C10","C5","C6","C1","C2","C3","C4","O2","C7",
               "C13","C12","C11","C8","C9","C14"}

def parse_filtered(filepath, atom_names):
    coords_by_res = {}
    with open(filepath) as f:
        f.readline()
        natoms = int(f.readline())
        for _ in range(natoms):
            line = f.readline()
            resid = int(line[0:5])
            atom_name = line[10:15].strip()
            if atom_name not in atom_names:
                continue
            x = float(line[20:28])
            y = float(line[28:36])
            z = float(line[36:44])
            coords_by_res.setdefault(resid, []).append((x, y, z))
        box = np.array(list(map(float, f.readline().strip().split())))
    residue_ids = sorted(coords_by_res.keys())
    coords = [np.array(coords_by_res[r]) for r in residue_ids]
    return residue_ids, coords, box

def centroid(points):
    return points.mean(axis=0)

def plane_normal(points):
    c = centroid(points)
    centered = points - c
    cov = centered.T @ centered
    _, vecs = np.linalg.eigh(cov)
    n = vecs[:, 0]
    return n / np.linalg.norm(n)

# Parse file
res_ids, residues_coords, box = parse_filtered(file_path, atom_filter)
nres = len(res_ids)
centers = np.array([centroid(p) for p in residues_coords])
normals = np.array([plane_normal(p) for p in residues_coords])

# Distance matrix under PBC
diff = centers[:, None, :] - centers[None, :, :]
diff -= box * np.round(diff / box)
dist_mat = np.sqrt((diff ** 2).sum(axis=-1))
np.fill_diagonal(dist_mat, np.inf)

nearest_idx = np.argmin(dist_mat, axis=1)
nearest_dist = dist_mat[np.arange(nres), nearest_idx]

# Angles
dot = np.abs(np.einsum("ij,ij->i", normals, normals[nearest_idx]))
angles = np.degrees(np.arccos(np.clip(dot, 0, 1)))

# ---------- Clustering ----------
data = np.vstack((nearest_dist, angles)).T
# Standardize
data_std = (data - data.mean(axis=0)) / data.std(axis=0)
# Choose k=3 (compact elbow in exploratory tests)
kmeans = KMeans(n_clusters=3, random_state=0, n_init=10).fit(data_std)
labels = kmeans.labels_
centroids_std = kmeans.cluster_centers_
centroids = centroids_std * data.std(axis=0) + data.mean(axis=0)

# ---------- Save with clusters ----------
out_path = "nearest_neighbours_with_clusters.txt"
with open(out_path, "w") as f:
    f.write("Mol_i\tMol_j_nearest\tDistance_nm\tAngle_deg\tCluster\n")
    for i,(j,d,a,c) in enumerate(zip([res_ids[idx] for idx in nearest_idx], nearest_dist, angles, labels)):
        f.write(f"{res_ids[i]}\t{j}\t{d:.6f}\t{a:.3f}\t{c}\n")

# ---------- Probability contour ----------
bins_d = np.linspace(nearest_dist.min(), nearest_dist.max(), 60)
bins_a = np.linspace(0, 90, 60)
H, xedges, yedges = np.histogram2d(nearest_dist, angles, bins=[bins_d, bins_a], density=True)

X, Y = np.meshgrid((xedges[:-1] + xedges[1:]) / 2, (yedges[:-1] + yedges[1:]) / 2)

plt.figure(figsize=(7,6))
plt.contourf(X, Y, H.T, levels=20)
plt.scatter(centroids[:,0], centroids[:,1], s=80, marker='x')
for idx,(cx,cy) in enumerate(centroids):
    plt.text(cx, cy+2, f"C{idx}", ha='center')
plt.xlabel("Nearest‑neighbour distance (nm)")
plt.ylabel("Angle between plane normals (°)")
plt.title("Joint probability density with K‑means cluster centroids")
plt.colorbar(label="Probability density")
plt.tight_layout()
plt.show()

# ---------- Summary table ----------
summary_df = pd.DataFrame({
    "Cluster": range(3),
    "Centre_dist_nm": centroids[:,0],
    "Centre_angle_deg": centroids[:,1],
    "Members": np.bincount(labels)
})

import ace_tools_open as tools; tools.display_dataframe_to_user(name="Cluster summary", dataframe=summary_df)
