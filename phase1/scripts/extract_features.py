from Bio import SeqIO
import pandas as pd
from collections import Counter
from itertools import product

AMINO_ACIDS = 'ACDEFGHIKLMNPQRSTVWY'

def get_kmers(seq, k=3):
    return [seq[i:i+k] for i in range(len(seq) - k + 1)]

def compute_kmer_counts(seq, k=3):
    kmers = get_kmers(seq, k)
    return dict(Counter(kmers))

def compute_aa_composition(seq):
    aa_counts = Counter(seq)
    total = len(seq)
    return {f'AA_{aa}': aa_counts.get(aa, 0) / total for aa in AMINO_ACIDS}

#set limit to 1000 to limit for this machine, can increase
def extract_features_from_fasta(fasta_path, k=3, limit=1000):
    feature_list = []
    records = SeqIO.parse(fasta_path, "fasta")
    for i, record in enumerate(records):
        if i == limit:
            break
        seq = str(record.seq)
        features = Counter([seq[j:j + k] for j in range(len(seq) - k + 1)])
        features = dict(features)
        features['ProteinID'] = record.id  # Add the protein ID
        feature_list.append(features)
    return pd.DataFrame(feature_list).fillna(0)



if __name__ == "__main__":
    fasta_path = "phase1/data/uniprot_proteins.fasta"
    features_df = extract_features_from_fasta(fasta_path, k=3)
    features_df.to_csv("phase1/data/sequence_features.csv", index=False)
    print("Extracted sequence features.")
