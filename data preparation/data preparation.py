import warnings
from Bio import BiopythonWarning
from Bio.PDB.PDBExceptions import PDBConstructionWarning

# Ignore PDBConstructionWarning
warnings.simplefilter('ignore', PDBConstructionWarning)

import numpy as np
from Bio.PDB import PDBParser, PDBList, parse_pdb_header
from Bio.PDB.Polypeptide import is_aa

class ProteinAnnotator:
    def __init__(self, pdb_dir='data\\'):
        self.pdb_dir = pdb_dir
        self.pdbl = PDBList()
        self.parser = PDBParser()

    def retrieve_and_parse_pdb(self, protein_name):
        self.pdbl.retrieve_pdb_file(protein_name, file_format='pdb', pdir=self.pdb_dir)
        structure = self.parser.get_structure(protein_name, f'{self.pdb_dir}pdb{protein_name.lower()}.ent')
        header = parse_pdb_header(f'{self.pdb_dir}pdb{protein_name.lower()}.ent')
        return structure, header

    def get_missing_residues(self, header):
        missing_res = [(res['ssseq'], res['res_name'], res['chain']) for res in header['missing_residues']]
        return sorted(missing_res, key=lambda x: x[0])

    def get_sequence(self, structure, missing_res):
        chain = next(structure.get_chains())
        sequence_list = []
        for residue in chain:
            if is_aa(residue, standard=True):
                sequence_list.append(residue.resname)
        
        for res in missing_res:
            if res[2] == chain.id:
                position, res_name, chain_id = res
                sequence_list.insert(position - 1, res_name)
        
        return ''.join(sequence_list)

    def get_helices_and_strands(self, protein_name):
        helices = []
        strands = []
        with open(f'{self.pdb_dir}pdb{protein_name.lower()}.ent', 'r') as file:
            for line in file:
                if line[:5] == 'HELIX' and line[19] == 'A' and line[31] == 'A' and line[38:40] == ' 1':
                    initSeqNum, endSeqNum = int(line[21:25]), int(line[33:37])
                    helices.append((initSeqNum, endSeqNum))
                if line[:5] == 'SHEET' and line[21] == 'A' and line[32] == 'A':
                    initSeqNum, endSeqNum = int(line[22:26]), int(line[33:37])
                    if (initSeqNum, endSeqNum) not in strands:
                        strands.append((initSeqNum, endSeqNum))
        return helices, strands

    def get_annotated_sequence(self, protein_name):
        structure, header = self.retrieve_and_parse_pdb(protein_name)
        missing_res = self.get_missing_residues(header)
        sequence = self.get_sequence(structure, missing_res)
        helices, strands = self.get_helices_and_strands(protein_name)
        return [sequence, helices, strands]

    def print_sequence(self, seq, helices, strands):
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

# Example usage:
annotator = ProteinAnnotator()
annotated_sequence = annotator.get_annotated_sequence('1FAT')
print(annotated_sequence)
annotator.print_sequence(*annotated_sequence)