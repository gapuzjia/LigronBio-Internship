import math
import matplotlib.pyplot as plt
from collections import defaultdict
from Bio import SeqIO
import csv
import os
import re

#load UniProt sequences
def load_fasta(path):
    seqs = {}
    for record in SeqIO.parse(path, "fasta"):
        uid = record.id.split("|")[1] if "|" in record.id else record.id
        seqs[uid] = str(record.seq)
    return seqs

#extract degron sequences by motif
def extract_degrons(txt_file, seq_dict):
    motif_seqs = defaultdict(list)
    with open(txt_file, newline='') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            uid = row['Entry']
            start, end = int(row['START']), int(row['END'])
            motif = row['ID_MANUSCRIPT']
            if uid in seq_dict and end <= len(seq_dict[uid]):
                subseq = seq_dict[uid][start-1:end]
                motif_seqs[motif].append(subseq)
    return motif_seqs

#filter groups with aligned sequences
def filter_aligned(motif_seqs):
    return {
        motif: seqs
        for motif, seqs in motif_seqs.items()
        if len(set(len(s) for s in seqs)) == 1
    }

#position frequency matrix (PFM)
#count number of occurences of an amino acid in a specific position
def compute_pfm(seqs):
    pfm = [defaultdict(int) for _ in range(len(seqs[0]))]
    for seq in seqs:
        for i, aa in enumerate(seq):
            pfm[i][aa] += 1
    return pfm

#position probability matrix (PPM)
#convert to probabilities 
def compute_ppm(pfm):
    return [
        {aa: count / sum(pos.values()) for aa, count in pos.items()}
        for pos in pfm
    ]

#compute info content (conservation)
def compute_ic(ppm):
    return [
        math.log2(20) + sum(p * math.log2(p) for p in pos.values() if p > 0)
        for pos in ppm
    ]

#plot sequence Llogo
def plot_logo(ppm, ic, title):
    fig, ax = plt.subplots(figsize=(len(ppm) * 0.5, 4))
    y_offset = [0] * len(ppm)

    for i, pos in enumerate(ppm):
        sorted_aa = sorted(pos.items(), key=lambda x: x[1] * ic[i])
        for aa, prob in sorted_aa:
            height = prob * ic[i]
            ax.text(i, y_offset[i], aa,
                    fontsize=10 + height * 20,
                    ha='center', va='bottom')
            y_offset[i] += height

    ax.set_xticks(range(len(ppm)))
    ax.set_xticklabels(range(1, len(ppm)+1))
    ax.set_ylim(0, max(y_offset) + 1)
    ax.set_ylabel("Info (bits)")
    ax.set_title(f"Sequence Logo: {title}")
    plt.tight_layout()

    #create output folder if needed
    output_dir = os.path.join("..", "seq_logos")
    os.makedirs(output_dir, exist_ok=True)

    #make filename safe
    safe_title = re.sub(r'[^\w\-_.]', '_', title)
    filename = os.path.join(output_dir, f"{safe_title}.svg")

    fig.savefig(filename, format='svg')
    print(f"Saved: {filename}")
    plt.close(fig)


#main
if __name__ == "__main__":
    fasta_path = "phase1/data/uniprot_proteins.fasta"
    degron_file = "phase1/data/Known degron_instance.txt"

    seq_dict = load_fasta(fasta_path)
    degrons = extract_degrons(degron_file, seq_dict)
    aligned = filter_aligned(degrons)

    for motif, seqs in aligned.items():
        print(f"Generating logo for: {motif} ({len(seqs)} sequences)")
        pfm = compute_pfm(seqs)
        ppm = compute_ppm(pfm)
        ic = compute_ic(ppm)
        plot_logo(ppm, ic, motif)
