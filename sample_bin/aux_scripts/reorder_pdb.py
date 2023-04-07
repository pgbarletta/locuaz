from pathlib import Path
import string

import MDAnalysis as mda

input_path = Path("")
# Change this:
reordering_indices = [2, 0, 1, 3]

u = mda.Universe(str(input_path / "input.pdb"))
seg_unis = []
for i, char in zip(reordering_indices, list(string.ascii_uppercase)):
    s = u.segments[i]
    v = mda.Universe.empty(len(s.atoms))
    v.atoms = s.atoms
    if s.segid != '':
        v.atoms.segments.segids = v.atoms.chainIDs = char
    seg_unis.append(v)
out_uni = mda.Merge(*[ uni.atoms for uni in seg_unis])
out_uni.dimensions = u.dimensions
out_uni.atoms.write(str(input_path / "fixed.pdb"))
