# Note of Flow of calcuation

After MD is completed. You neeed to remove booundaries of your final production run using command

```bash
echo -e "0\nq" | gmx_mpi trjconv -f npt.gro -s npt.tpr -pbc mol -o output_whole.gro
sleep 10
echo 0 | gmx_mpi trjconv -f npt.trr -s npt.tpr -pbc mol -sep -o frame.gro
```
where npt.gro is your final result and you have at least 100 tracked structure in the .trr file

run extractdistancebetween2molecules.py. Before running that file, you need to change the atom you need to used such as


```python
# ---------- re‑use parsing from previous cell ----------
file_path = "output_whole.gro"
atom_filter = {"O1","C10","C5","C6","C1","C2","C3","C4","O2","C7",
               "C13","C12","C11","C8","C9","C14"}
```

the atom_filters define the plane or structure you want to investigate.

If it is a plane, then the normal between two plane is obtained . 

# Plane Normal Calculation via PCA

In that snippet you’re feeding an $N \times 3$ array `xyz` of Cartesian coordinates (one row per atom).  
The helper functions do two things:

## Step-by-step Explanation

| Step               | Math                                                                                   | Meaning                                                           |
|--------------------|----------------------------------------------------------------------------------------|-------------------------------------------------------------------|
| `centroid(xyz)`    | $$ \mathbf{r}_{\text{cm}} = \frac{1}{N} \sum_{i=1}^{N} \mathbf{r}_i $$                 | Geometric centre of the selected atoms                            |
| `plane_normal(xyz)`| 1. Centre the coordinates: $$ \tilde{\mathbf{r}}_i = \mathbf{r}_i - \mathbf{r}_{\text{cm}} $$ <br> 2. Build the $3 \times 3$ covariance matrix: $$ C = \tilde{R}^\top \tilde{R} $$ <br> 3. Diagonalise $C$: $$ C \mathbf{v}_k = \lambda_k \mathbf{v}_k $$ <br> 4. Choose the eigen‑vector $\mathbf{v}_{\min}$ with the smallest eigen‑value | Principal-component analysis (PCA) of the point cloud |

Because the smallest eigen‑value corresponds to the axis with the **least variance**,  
its eigen‑vector $\mathbf{v}_{\min}$ is orthogonal to the best‑fit plane that minimizes  
the sum of squared distances to all points (i.e., a **total‑least‑squares fit**).


```python
n = v[:, 0]                 # eigen‑vector with lowest eigen‑value
return n / np.linalg.norm(n)  # unit‑length normal
```

So plane_normal returns a unit vector that points perpendicular to the plane that best fits the chosen atoms (e.g., the peptide/ligand ring defined by N2, C30, …, C38).

Applications of the Plane Normal
You can use this normal vector $\hat{n}$ to:

Compute the tilt of the ring relative to another surface

Find dihedral angles between two such planes

Cluster similar orientations (used later in your k-means step)

In short:
“Plane normal” here is just a direction vector $\hat{n}$ at right angles to the molecular plane defined by the filtered atoms,
obtained through a PCA-style least-squares fit.
