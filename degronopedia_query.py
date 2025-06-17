import pandas as pd
import requests
import time

# Load your CSV
df = pd.read_csv('degrons.csv')

# Ensure UniProtID is treated as string
df['UniProtID'] = df['UniProtID'].astype(str)

#Get non-empty UniProt IDs
proteins = df['UniProtID'].dropna()
proteins = [p for p in proteins if p.strip() != '' and p != 'nan']

#API endpoint
base_url = "https://rest.uniprot.org/uniprotkb/search"

results = []

for protein in proteins:
    query = f"accession:{protein} AND organism_id:9606"
    params = {
        'query': query,
        'format': 'tsv',
        'fields': 'accession,id,protein_name,gene_names,organism_name'
    }

    response = requests.get(base_url, params=params)
    time.sleep(0.5)  # Rate limit buffer

    if response.ok and response.text.strip():
        lines = response.text.strip().split('\n')
        if len(lines) > 1:
            fields = lines[1].split('\t')  # Skip header
            results.append({
                'Query': protein,
                'UniProt_ID': fields[0],
                'Entry_Name': fields[1],
                'Protein_Name': fields[2],
                'Gene_Names': fields[3],
                'Organism': fields[4]
            })
        else:
            results.append({
                'Query': protein,
                'UniProt_ID': protein,
                'Entry_Name': None,
                'Protein_Name': None,
                'Gene_Names': None,
                'Organism': None
            })
    else:
        results.append({
            'Query': protein,
            'UniProt_ID': protein,
            'Entry_Name': None,
            'Protein_Name': None,
            'Gene_Names': None,
            'Organism': None
        })

#Save to CSV
df_results = pd.DataFrame(results)
df_results.to_csv('protein_to_uniprot.csv', index=False)

print("Saved UniProt ID metadata to protein_to_uniprot.csv")
