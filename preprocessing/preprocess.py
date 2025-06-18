import pandas as pd
import requests
import time

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

raw_df = pd.read_csv('rawdata.csv')
UniProtIDs = raw_df['UniProtID'].tolist()

df['UniProtID'] = UniProtIDs

# Save to CSV
df.to_csv('degrons.csv', index=False)


