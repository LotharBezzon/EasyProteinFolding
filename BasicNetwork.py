import numpy as np
import torch
import torch.nn as nn
from DataPreparation.DataPreparation import ProteinAnnotator

# Define a mapping from three-letter amino acid codes to integers
amino_acid_to_index = {
    'ALA': 1, 'CYS': 2, 'ASP': 3, 'GLU': 4, 'PHE': 5, 'GLY': 6, 'HIS': 7, 
    'ILE': 8, 'LYS': 9, 'LEU': 10, 'MET': 11, 'ASN': 12, 'PRO': 13, 'GLN': 14, 
    'ARG': 15, 'SER': 16, 'THR': 17, 'VAL': 18, 'TRP': 19, 'TYR': 20
}

def split_sequence(sequence):
    return [sequence[i:i+3] for i in range(0, len(sequence), 3)]

def tokenize_sequence(sequence):
    return [amino_acid_to_index[aa] for aa in split_sequence(sequence)]

def positional_encoding(max_len, d_model):
    pos = np.arange(max_len)[:, np.newaxis]
    i = np.arange(d_model)[np.newaxis, :]
    angle_rates = 1 / np.power(10000, (2 * (i // 2)) / np.float32(d_model))
    angle_rads = pos * angle_rates

    # Apply sin to even indices in the array; 2i
    angle_rads[:, 0::2] = np.sin(angle_rads[:, 0::2])

    # Apply cos to odd indices in the array; 2i+1
    angle_rads[:, 1::2] = np.cos(angle_rads[:, 1::2])

    pos_encoding = angle_rads[np.newaxis, ...]

    return torch.tensor(pos_encoding, dtype=torch.float32)

class TransformerInputLayer(nn.Module):
    def __init__(self, vocab_size, d_model, max_len):
        super(TransformerInputLayer, self).__init__()
        self.embedding = nn.Embedding(vocab_size, d_model)
        self.pos_encoding = positional_encoding(max_len, d_model)
        self.layer_norm = nn.LayerNorm(d_model)

    def forward(self, x):
        seq_len = x.size(1)
        x = self.embedding(x)  # (batch_size, input_seq_len, d_model)
        x *= torch.sqrt(torch.tensor(x.size(-1), dtype=torch.float32))
        x += self.pos_encoding[:, :seq_len, :]
        x = self.layer_norm(x)
        return x

# Example usage:
annotator = ProteinAnnotator()
sequence, helices, strands = annotator.get_annotated_sequence('1FAT')

# Tokenize the sequence
tokenized_sequence = tokenize_sequence(sequence)

# Define the embedding layer
embedding_dim = 128  # Dimension of the embedding vectors
vocab_size = len(amino_acid_to_index) + 1  # Vocabulary size (+1 for padding)
max_len = len(tokenized_sequence)  # Maximum length of the sequence

# Convert tokenized sequence to tensor
tokenized_sequence_tensor = torch.tensor([tokenized_sequence], dtype=torch.long)

# Create the transformer input layer
transformer_input_layer = TransformerInputLayer(vocab_size, embedding_dim, max_len)

# Get the embedded sequence with positional encoding and layer normalization
embedded_sequence = transformer_input_layer(tokenized_sequence_tensor)
print("Embedded Sequence with Positional Encoding and Layer Normalization:", embedded_sequence)