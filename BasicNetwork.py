import numpy as np
from Bio.PDB import PDBParser, PDBList, parse_pdb_header
from Bio.PDB.Polypeptide import is_aa
from DataPreparation.DataPreparation import ProteinAnnotator


annotator = ProteinAnnotator()
annotated_sequence = annotator.get_annotated_sequence('1FAT')
annotator.print_sequence(*annotated_sequence)