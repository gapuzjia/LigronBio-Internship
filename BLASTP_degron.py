import pandas as pd
from Bio.Blast import NCBIWWW, NCBIXML
from Bio import Entrez, SeqIO
import time

# Set your email for NCBI access (required)
Entrez.email = "your_email@example.com"  # Replace with your real email

# Load degrons
df = pd.read_csv('degrons.csv')
df['Degron'] = df['Degron'].astype(str).str.strip()
degrons = [d for d in df['Degron'].unique() if d and d.lower() != 'nan']

results = []

for i, degron in enumerate(degrons):
    print(f"🔍 {i+1}/{len(degrons)} | BLASTing: {degron}")

    # Wrap the full block in try-except
    try:
        # Step 1: Run BLASTP
        blast_result = NCBIWWW.qblast("blastp", "nr", degron, format_type="XML", hitlist_size=1)
        blast_record = NCBIXML.read(blast_result)

        if blast_record.alignments:
            alignment = blast_record.alignments[0]
            accession = alignment.accession  # safer than trying to find "gi|"

            # Step 2: Fetch full protein GenBank record
            handle = Entrez.efetch(db="protein", id=accession, rettype="gb", retmode="text")
            record = SeqIO.read(handle, "genbank")

            # Step 3: Store data
            results.append({
                "Degron": degron,
                "Matched_Protein_ID": record.id,
                "Description": record.description,
                "Full_Sequence": str(record.seq)
            })

            handle.close()
        else:
            results.append({
                "Degron": degron,
                "Matched_Protein_ID": None,
                "Description": "No hit",
                "Full_Sequence": None
            })

    except Exception as e:
        print(f"⚠️ Error processing {degron}: {e}")
        results.append({
            "Degron": degron,
            "Matched_Protein_ID": "Error",
            "Description": str(e),
            "Full_Sequence": None
        })

    # Delay to avoid NCBI rate-limiting
    time.sleep(15)

# Save results
output_df = pd.DataFrame(results)
output_df.to_csv('degron_blast_origin_results.csv', index=False)
print("✅ Results saved to degron_blast_origin_results.csv")
