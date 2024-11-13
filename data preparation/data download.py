from Bio.PDB import PDBList

pdbl = PDBList()
entries = pdbl.get_all_entries()
print(entries)
print(len(entries))