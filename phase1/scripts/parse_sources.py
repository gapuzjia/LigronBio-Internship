import pandas as pd
from Bio import SeqIO

# === Original parsing functions (unchanged) ===

def parse_predicted():
    df = pd.read_csv("phase1/data/predicted_degrons_properties.csv")
    return pd.DataFrame({
        'ProteinID': df['Entry'],
        'Start': df['START'],
        'End': df['END'],
        'Motif': df['DEGRON'],
        'Hit': df['Hit'],
        'Prob_DEGRON': df['Prob_DEGRON'],
        'Prob_RANDOM': df['Prob_RANDOM'],
        'Predicted_Class': df['Predicted_Class'],
        'SolventAccessibility': df['ASA_SCORE'],
        'Conservation': df['CONS_SCORE'],
        'Rigidity': df['RIG_SCORE'],
        'Disorder': df['ANCHOR_SCORE'],
        'PTMs': df['ptms_flanking'],
        'Nearby_Lysines': df['nflanking_ub_lysines'],
        'Source': 'predicted'
    })

def parse_known():
    df = pd.read_csv("phase1/data/Known degron_instance.txt", sep='\t')
    df['Label'] = 1
    return pd.DataFrame({
        'ProteinID': df['Entry'],
        'Start': df['START'],
        'End': df['END'],
        'Motif': df['DEGRON'],
        'Label': df['Label'],
        'Source': 'known'
    })

def parse_motif_classes():
    df = pd.read_csv("phase1/data/Degron_class.txt", sep='\t')
    return df[['ELMIdentifier', 'Motif', 'Regex']]

# === New: Load UniProt FASTA ===
def load_fasta(fasta_path):
    seqs = {}
    for record in SeqIO.parse(fasta_path, "fasta"):
        uid = record.id.split("|")[1] if "|" in record.id else record.id
        seqs[uid] = str(record.seq)
    return seqs

# === New: Append degron sequence and save as separate file ===
def save_known_with_sequences(fasta_dict, input_path, output_path):
    df = pd.read_csv(input_path, sep='\t')
    degron_seqs = []
    for _, row in df.iterrows():
        uid = row['Entry']
        start, end = int(row['START']), int(row['END'])
        full_seq = fasta_dict.get(uid)
        if full_seq and end <= len(full_seq):
            degron_seqs.append(full_seq[start - 1:end])
        else:
            degron_seqs.append("N/A")
    df['Degron_Seq'] = degron_seqs
    df.to_csv(output_path, sep='\t', index=False)
    print(f"Saved: {output_path}")

# === Main script logic ===
if __name__ == "__main__":
    pred = parse_predicted()
    known = parse_known()
    motifs = parse_motif_classes()

    pred.to_csv("phase1/data/predicted_degrons_clean.csv", index=False)
    known.to_csv("phase1/data/known_degrons_clean.csv", index=False)
    motifs.to_csv("phase1/data/degron_motif_classes.csv", index=False)

    # Save known degrons with sequence
    fasta_dict = load_fasta("phase1/data/uniprot_proteins.fasta")
    save_known_with_sequences(
        fasta_dict,
        "phase1/data/Known degron_instance.txt",
        "phase1/data/known_degrons_with_seq.csv"
    )

    print("Parsed all degron source files.")
