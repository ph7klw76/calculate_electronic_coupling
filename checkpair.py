import numpy as np
import glob
import os
import csv

def read_xyz_coordinates(filename):
    """
    Reads an XYZ file and returns a NumPy array of shape (N, 3)
    containing the x, y, z coordinates of each atom.
    Skips first two header lines.
    """
    coords = []
    with open(filename, 'r') as f:
        # Skip atom count and comment
        next(f)
        next(f)
        for line in f:
            parts = line.split()
            if len(parts) < 4:
                continue
            x, y, z = map(float, parts[1:4])
            coords.append([x, y, z])
    return np.array(coords)

def compute_centroid(coords):
    """Return the mean position of an (nÃ—3) array of coordinates."""
    return coords.mean(axis=0)

def process_folder(folder_path, output_csv="distances.csv"):
    """
    Loops over all .xyz files in folder_path, computes centroid-distance,
    and writes results to a CSV.
    """
    xyz_files = glob.glob(os.path.join(folder_path, "*.xyz"))
    results = []

    for xyz in xyz_files:
        coords = read_xyz_coordinates(xyz)
        if coords.shape[0] < 176:
            print(f"Skipping {os.path.basename(xyz)}: only {coords.shape[0]} atoms")
            continue

        first_88 = coords[:88]
        last_88  = coords[-88:]
        c1 = compute_centroid(first_88)
        c2 = compute_centroid(last_88)
        dist = np.linalg.norm(c1 - c2)

        results.append((os.path.basename(xyz), dist))

    # Write out CSV
    with open(output_csv, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(["filename", "distance"])
        writer.writerows(results)

    print(f"Processed {len(results)} files. Results saved to {output_csv}")

if __name__ == "__main__":
    # Change this to your target directory:
    folder = "xyz_folder"
    process_folder(folder)
