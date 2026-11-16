from concurrent.futures import ProcessPoolExecutor, as_completed
import pandas as pd
from pathlib import Path
import warnings
import numpy as np
from tqdm import tqdm
import argparse
import os
import subprocess

import biotite.structure as struc
import biotite.structure.io as strucio
import biotite.structure.io.pdbx as pdbx
import biotite.structure.io.pdb as pdb

def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-i",
        "--input",
        type=Path,
        required=True,
        help="Path to the folder with af3 output containing confidence metrics",
    )
    parser.add_argument(
            "-o",
            "--output",
            type=Path,
            required=True,
            help="Path to the output folder which will contain the converted pdb files (folder can be created with this script)",
        )
    return parser.parse_args()


def cif_to_pdb(path_model, output_dir):
    cif_file = pdbx.CIFFile.read(path_model)
    model = pdbx.get_structure(cif_file)[0]

    pdb_file = pdb.PDBFile()
    pdb.set_structure(pdb_file, model)

    base_name = os.path.splitext(os.path.basename(path_model))[0] + ".pdb"
    out_pdb_path = os.path.join(output_dir, base_name)
    pdb_file.write(out_pdb_path)

    # pdb_file.write(f"{os.path.splitext(path_model)[0]}.pdb")
    # return f"{os.path.splitext(path_model)[0]}.pdb"

def main():
    args = parse_args()
    input_path = args.input
    outdir = args.output
    os.makedirs(outdir, exist_ok = True)
    #Iterate over the files in the input path
    
    for m in os.listdir(input_path):
        if m.endswith(".cif"):
            cif_to_pdb(os.path.join(input_path, m), outdir)

if __name__ == "__main__":
    main()
