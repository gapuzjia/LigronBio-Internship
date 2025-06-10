import pandas as pd
import requests
import time

# Load protein names from your CSV
df = pd.read_csv('UniProtIDProteins.csv')
proteins = df['Protein'].tolist()

# UniProt API endpoint
base_url = "https://rest.uniprot.org/uniprotkb/search"

results = []

for protein in proteins:
    query = f"{protein}+AND+organism_id:9606"
    params = {
        'query': query,
        'format': 'tsv',
        'fields': 'accession,id,protein_name,gene_names,organism_name'
    }

    response = requests.get(base_url, params=params)

    # Throttle to avoid getting rate-limited
    time.sleep(0.5)

    if response.ok and response.text.strip():
        lines = response.text.strip().split('\n')
        if len(lines) > 1:
            fields = lines[1].split('\t')  # skip header
            results.append({
                'Query': protein,
                'UniProt_ID': fields[0],
                'Entry_Name': fields[1],
                'Protein_Name': fields[2],
                'Gene_Names': fields[3],
                'Organism': fields[4]
            })
    else:
        results.append({'Query': protein, 'UniProt_ID': None})

# Save to CSV
df_results = pd.DataFrame(results)
df_results.to_csv('protein_to_uniprot.csv', index=False)

print("Saved UniProt ID mappings to protein_to_uniprot.csv")
