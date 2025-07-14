import pandas as pd

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

if __name__ == "__main__":
    pred = parse_predicted()
    known = parse_known()
    motifs = parse_motif_classes()

    pred.to_csv("phase1/data/predicted_degrons_clean.csv", index=False)
    known.to_csv("phase1/data/known_degrons_clean.csv", index=False)
    motifs.to_csv("phase1/data/degron_motif_classes.csv", index=False)

    print("Parsed all degron source files.")
