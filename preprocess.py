import pandas as pd

# Read the text file
with open('G-loop.txt', 'r') as file:
    lines = file.read().strip().splitlines()

# Parse into list of dictionaries
entries = []
for i in range(0, len(lines), 2):
    if lines[i].startswith(">") and i+1 < len(lines):
        name = lines[i][1:].strip()
        seq = lines[i+1].strip()
        entries.append({'Protein': name, 'Degron': seq})

# Convert to DataFrame
df = pd.DataFrame(entries)

# Save to CSV
df.to_csv('degrons.csv', index=False)

df_proteins = df[['Protein']]
df.to_csv('UniProtID Degrons', index=False)