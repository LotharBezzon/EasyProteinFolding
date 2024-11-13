import numpy as np
from Bio.PDB import PDBParser, PDBList, parse_pdb_header
from Bio.PDB.Polypeptide import is_aa

pdbl = PDBList()
pdbl.retrieve_pdb_file("1FAT", file_format='pdb', pdir=f'data\\')

# Parse the PDB file
parser = PDBParser()
structure = parser.get_structure('1FAT', f'data\pdb1fat.ent')
h = parse_pdb_header('data\pdb1fat.ent')
missing_res = [(res['ssseq'], res['res_name'], res['chain']) for res in h['missing_residues']]

# Sort missing_res by the first index (sequence position)
missing_res = sorted(missing_res, key=lambda x: x[0])

chain = next(structure.get_chains())
sequence_list = []
for residue in chain:
    if is_aa(residue, standard=True):
        sequence_list.append(residue.resname)

# Insert missing residues at the specified positions
for res in missing_res:
    if res[2] == chain.id:
        position, res_name, chain_id = res
        sequence_list.insert(position - 1, res_name) 

sequence = ''.join(sequence_list)
print(sequence)

with open('data\pdb1fat.ent', 'r') as file:
    helices = []
    strands = []
    for line in file:
        # Search for helices of type ' 1' (alpha), belonging to chain 'A'
        if line[:5]=='HELIX' and line[19]=='A' and line[31]=='A' and line[38:40]==' 1':
            initSeqNum, endSeqNum = int(line[21:25]), int(line[33:37])
            helices.append((initSeqNum, endSeqNum))
            print(line, initSeqNum, endSeqNum)
        # Search for strands (alpha) belonging to chain 'A'
        if line[:5]=='SHEET' and line[21]=='A' and line[32]=='A':
            initSeqNum, endSeqNum = int(line[22:26]), int(line[33:37])
            # must exclude duplicate from beta barrels and biforcations
            if (initSeqNum, endSeqNum) not in strands:
                strands.append((initSeqNum, endSeqNum))
            print(line, initSeqNum, endSeqNum)

annotated_sequence = [sequence, helices, strands]

def print_sequence(seq, helices, strands):
    colored_seq = list(seq)
    
    # Apply red background for helices
    for start, end in helices:
        for i in range((start - 1) * 3, end * 3):
            colored_seq[i] = f"\033[41m{colored_seq[i]}\033[0m"  # 41 is the code for red background
    
    # Apply blue background for strands
    for start, end in strands:
        for i in range((start - 1) * 3, end * 3):
            colored_seq[i] = f"\033[44m{colored_seq[i]}\033[0m"  # 44 is the code for blue background
    
    print("".join(colored_seq))

print_sequence(sequence, helices, strands)

