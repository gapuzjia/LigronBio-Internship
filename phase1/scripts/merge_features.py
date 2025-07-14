import pandas as pd

def load_mutation_data(mutation_file, degfi_file):
    mut_df = pd.read_csv(mutation_file)
    degfi_df = pd.read_csv(degfi_file)


    # Clean and reduce to necessary fields
    mut_df = mut_df[['Entry', 'DEGRON']]
    mut_df = mut_df.rename(columns={'DEGRON': 'mut_in_degron'})

    mut_df['mut_in_degron'] = mut_df['mut_in_degron'].notnull().astype(int)
    mut_df = mut_df.groupby('Entry', as_index=False).max()

    degfi_df = degfi_df[['Entry', 'DegFI_score']]
    degfi_df = degfi_df.groupby('Entry', as_index=False).mean()

    return mut_df.rename(columns={'Entry': 'ProteinID'}), degfi_df.rename(columns={'Entry': 'ProteinID'})

def merge_features(base_file, mutation_df, degfi_df):
    base_df = pd.read_csv(base_file)

    merged = base_df.merge(mutation_df, on='ProteinID', how='left')
    merged = merged.merge(degfi_df, on='ProteinID', how='left')

    merged['mut_in_degron'] = merged['mut_in_degron'].fillna(0)
    merged['DegFI_score'] = merged['DegFI_score'].fillna(0)

    return merged

if __name__ == "__main__":
    base_file = "phase1/data/features/labeled_sequence_features.csv"
    mutation_file = "phase1/data/all_missense_mutation_related_degron_part.csv"
    degfi_file = "phase1/data/all_missense_mutation_related_degron_DegFI_score_part.csv"


    mutation_df, degfi_df = load_mutation_data(mutation_file, degfi_file)
    merged_df = merge_features(base_file, mutation_df, degfi_df)

    merged_df.to_csv("phase1/data/FEATURE_MATRIX.csv", index=False)
    merged_df.to_csv("FEATURE_MATRIX.csv", index=False)
    print("Final feature matrix saved.")
