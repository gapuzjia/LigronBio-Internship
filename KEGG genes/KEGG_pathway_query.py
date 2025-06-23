import requests
import time
import pandas as pd

# Load KEGG gene list from the .txt file
with open("kegg_genes_compiled.txt", "r") as f:
    lines = f.readlines()

# Extract gene IDs like hsa:351
gene_ids = [line.strip().split('\t')[1] for line in lines]

results = []

# Fetch amino acid sequence for each gene
for i, gene_id in enumerate(gene_ids):
    print(f"[{i+1}/{len(gene_ids)}] Fetching {gene_id}...")
    url = f"http://rest.kegg.jp/get/{gene_id}/aaseq"
    response = requests.get(url)

    if response.ok and response.text.startswith('>'):
        seq = ''.join(response.text.split('\n')[1:])
        results.append({'Gene_ID': gene_id, 'Protein_Sequence': seq})
    else:
        print(f"⚠️ Failed to fetch {gene_id}")
        results.append({'Gene_ID': gene_id, 'Protein_Sequence': None})

    time.sleep(1)

# Save to CSV
df = pd.DataFrame(results)
df.to_csv("kegg_genes_proteins.csv", index=False)
print("Done! Saved to kegg_genes_proteins.csv")
