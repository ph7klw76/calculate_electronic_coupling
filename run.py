import shutil
import os
from pathlib import Path
import ECOUPLING as EC
import glob

# Define the current directory
current_directory = os.getcwd()

# Read the 'nearest_neighbor_distances.txt' file and process line by line
with open('nearest_neighbours_with_clusters-acceptor.txt', 'r') as f1:
    next(f1) 
    for line in f1:
        # Split line into components and retrieve necessary values
        line_parts = line.split()
        if len(line_parts) < 5:
            continue  # Skip malformed lines (if any)
        
        line0, line1 = line_parts[0], line_parts[1]
        
        # Set up source folder and files to copy
        source_folder = f'{line0}_{line1}'
        file_to_copy = ['FILE.53', 'FILE.54', f'fort.7-{line0}', f'fort.7-{line1}']
        
        # Copy files from the source folder to the current directory
        for file_name in file_to_copy:
            source_file_path = Path(source_folder) / file_name
            try:
                shutil.copy(source_file_path, current_directory)
                print(f"Copied {source_file_path} to {current_directory}")
            except Exception as e:
                print(f"Error copying {source_file_path}: {e}")

        # Create or overwrite 'inFile.in'
        with open('inFile.in', 'w') as f2:
            f2.write('FILE.53\n')
            f2.write('FILE.54\n')
            f2.write(f'fort.7-{line0}\n')
            f2.write(f'fort.7-{line1}\n')

        # Call ECOUPLING method
        try:
            EC.coupling(['a', 'inFile.in', 'paramFile.txt', f'output{line0}_{line1}'])
            print(f"Coupling executed for {line0}_{line1}")
        except Exception as e:
            print(f"Error executing ECOUPLING for {line0}_{line1}: {e}")

# Remove all files starting with 'fort'
files_to_remove = glob.glob('fort*')
for file in files_to_remove:
    try:
        os.remove(file)
        print(f"Removed: {file}")
    except Exception as e:
        print(f"Error removing {file}: {e}")
