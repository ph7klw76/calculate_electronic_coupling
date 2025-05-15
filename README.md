# Note of Flow of calcuation

### 1. After MD is completed. You neeed to remove booundaries of your final production run using command

```bash
echo -e "0\nq" | gmx_mpi trjconv -f npt.gro -s npt.tpr -pbc mol -o output_whole.gro
sleep 10
echo 0 | gmx_mpi trjconv -f npt.trr -s npt.tpr -pbc mol -sep -o frame.gro
```
where npt.gro is your final result and you have at least 100 tracked structure in the .trr file
Put all you frame files into anuther folder called frame for better storage.

### 2. copy extractdistancebetween2molecules.py in the same folder of the data. Before running that file, you need to change the atom you need to used such as


```python
# ---------- re‑use parsing from previous cell ----------
file_path = "output_whole.gro"
atom_filter = {"O1","C10","C5","C6","C1","C2","C3","C4","O2","C7",
               "C13","C12","C11","C8","C9","C14"}
```

the atom_filters define the plane or structure you want to investigate.

If it is a plane, then the normal between two plane is obtained . 

In that snippet you’re feeding an $N \times 3$ array `xyz` of Cartesian coordinates (one row per atom).  
The helper functions do two things:


![image](https://github.com/user-attachments/assets/155063c7-296d-4f0e-a704-deefa2f35334)

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

### 3. run the ./extract_pairs2.sh output_whole.gro index.ndx nearest_molecule_distances.txt 0.40 0.40 0062100000

the nearest_molecule_distances.txt is the file you obatined from extractdistancebetween2molecules.py and 0.4 is the target distance followed by tolerance and the last value is the tuned omega value using wb97x functional for the molecule following gaussian 09 convention 
This will create all the caussian file.

### 4. Check the molecular pair structure by running checkpair.py. this will check any outlinier results which because of matching boundary results in being far away of each other.

this can be done by

 ```python
        if coords.shape[0] < 176: # change this , to be total number of atoms
            print(f"Skipping {os.path.basename(xyz)}: only {coords.shape[0]} atoms")
            continue

        first_88 = coords[:88]  # change this to be the number of atoms of first pair
        last_88  = coords[-88:] # change this to be the number of atoms of second pair
```

### 5. then run the delete.sh file
   
### 6. To submit for calculation in that folder run submit.sh note there will be a lot of files.you might want to break them into a few folder

### 7. In the frame folder you can run distributionplot.py. susing the same modified atom_filter in step 2. This will create a 2-D contour plot such as

   ![image](https://github.com/user-attachments/assets/8d4fdd7a-3335-4c20-b0cf-19d0ddcb0b34)





