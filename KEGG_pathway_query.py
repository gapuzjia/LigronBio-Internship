import requests
from time import sleep
from Bio import SeqIO
from io import StringIO

#set KEGG pathway
pathway_id = "hsa05010"

#extract gene entries linked to the pathway
print(f"🔍 Fetching gene entries for pathway {pathway_id}...")
link_url = f"http://rest.kegg.jp/link/genes/{pathway_id}"
response = requests.get(link_url)

if not response.ok:
    raise Exception(f"Failed to get gene list from KEGG. Status: {response.status_code}")

gene_entries = [line.split('\t')[1].strip() for line in response.text.strip().split('\n')]
print(f"Found {len(gene_entries)} genes.")

#extract protein sequences (aaseq) for each gene
fasta_records = []

for i, gene_id in enumerate(gene_entries):
    print(f"[{i+1}/{len(gene_entries)}] Fetching sequence for {gene_id}")
    aaseq_url = f"http://rest.kegg.jp/get/{gene_id}/aaseq"
    r = requests.get(aaseq_url)

    if r.ok and r.text.startswith('>'):
        try:
            record = SeqIO.read(StringIO(r.text), "fasta")
            fasta_records.append(record)
        except Exception as e:
            print(f"Could not parse sequence for {gene_id}: {e}")
    else:
        print(f"No sequence found for {gene_id}")

    #KEGG rate limits
    sleep(1)  

#save into file
output_file = f"{pathway_id}_proteins.fasta"
with open(output_file, "w") as handle:
    SeqIO.write(fasta_records, handle, "fasta")

print(f"\nSaved {len(fasta_records)} protein sequences to {output_file}")
