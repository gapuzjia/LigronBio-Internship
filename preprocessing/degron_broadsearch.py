import pandas as pd
import requests
import time

#load data
df = pd.read_csv('degrons.csv')
df = df.dropna(subset=['Degron', 'UniProtID'])
df['Degron'] = df['Degron'].astype(str).str.strip()
df['UniProtID'] = df['UniProtID'].astype(str).str.strip()

#create list of known degrons
unique_degrons = df['Degron'].unique()

#create list of known IDs
unique_uniprot_ids = df['UniProtID'].unique()

#fetch sequences for each UniProtID
protein_seqs = {}
for uid in unique_uniprot_ids:
    url = f"https://rest.uniprot.org/uniprotkb/{uid}.fasta"
    try:
        response = requests.get(url)
        if response.ok:
            lines = response.text.splitlines()
            seq = ''.join(lines[1:])
            protein_seqs[uid] = seq
        else:
            protein_seqs[uid] = None
    except Exception as e:
        print(f"Error fetching {uid}: {e}")
        protein_seqs[uid] = None
    time.sleep(0.5)

#search each degron across degron databases
results = []
for degron in unique_degrons:
    matched_ids = []
    for uid, seq in protein_seqs.items():
        if seq and degron in seq:
            matched_ids.append(uid)
    results.append({
        'Degron': degron,
        'Match_Count': len(matched_ids),
        'Matched_UniProtIDs': ', '.join(matched_ids)
    })

#write results
df_result = pd.DataFrame(results)
df_result.to_csv('degron_matches_from_uniprot.csv', index=False)

print("Saved matches to degron_matches_from_uniprot.csv")
