import argparse

import numpy as np
import MDAnalysis as mda
from MDAnalysis.analysis.distances import dist

def get_cms_diff(u: mda.Universe):
    target = u.segments[0]
    binder = u.segments[1]
    
    return np.linalg.norm(target.atoms.center_of_mass() - binder.atoms.center_of_mass())

if __name__ == "__main__":
   
    parsero = argparse.ArgumentParser()
    parsero.add_argument("in_pdb", type = str, help = "Input PDB file")

    args = parsero.parse_args()

    u = mda.Universe(args.in_pdb)
    print(get_cms_diff(u) / 10)
