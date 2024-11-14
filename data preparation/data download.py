from Bio.PDB import PDBList
import pickle
import os
import random

class data_download:
    def __init__(self, pdb_dir='data\\'):
        self.pdb_dir = pdb_dir
        self.pdbl = PDBList()
        self.entries = self.get_entries()
    
    def get_entries(self):
        if os.path.exists('entries.pkl'):
            # Load entries from the file if it exists
            with open('entries.pkl', 'rb') as f:
                entries = pickle.load(f)
        else:
            # Get all entries and save them to the file if it doesn't exist
            entries = self.pdbl.get_all_entries()
            with open('entries.pkl', 'wb') as f:
                pickle.dump(entries, f)
        return entries

    # Download n random PDB files (if not already downloaded)
    def download_random_pdb_files(self, n=10):
        random_entries = random.sample(self.entries, n)
        self.pdbl.download_pdb_files(random_entries, pdir=self.pdb_dir, file_format='pdb')

    def get_protein_names(self):
        file_names = os.listdir(self.pdb_dir)
        return [name[3:7].upper() for name in file_names]

data = data_download()
data.download_random_pdb_files(n=100)
names = data.get_protein_names()   
print(names)
print(len(names))