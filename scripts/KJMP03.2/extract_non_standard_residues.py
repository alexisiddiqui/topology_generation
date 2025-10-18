#


"""
Pymol operations
Remove:
# dodgy termini
sel resi 301+302+303+304
# hetatoms
sel resn DMS+MPO+GOL+PGE


"""

non_standard_resnames = ["YFT", "NMR", "NMO", "1RJ"]

cleaned_mol2_path = "raw/KJMP03.2/exclude_ter_from_mol2/KJMP03.2castruc_TER_removed.mol2"
