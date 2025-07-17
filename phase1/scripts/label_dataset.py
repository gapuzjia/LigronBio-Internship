import pandas as pd

def label_proteins(sequence_feature_file, known_degrons_file):
    features_df = pd.read_csv(sequence_feature_file)
    known_df = pd.read_csv(known_degrons_file)

    # Set of protein IDs that have at least one known degron
    known_protein_ids = set(known_df['ProteinID'])

    # Assign labels
    features_df['Label'] = features_df['ProteinID'].apply(lambda pid: 1 if pid in known_protein_ids else 0)

    return features_df

if __name__ == "__main__":
    labeled_df = label_proteins(
        "phase1/data/sequence_features.csv",
        "phase1/data/known_degrons_clean.csv"
    )
    labeled_df.to_csv("phase1/data/features/labeled_sequence_features.csv", index=False)
    print("Labeled dataset saved.")
