import pandas as pd
import re
import math
from collections import Counter

pileup = pd.read_csv('data-work/CLIP-let7g-gene.pileup', sep='\t', names=['chrom', 'pos', '_ref', 'count', 'basereads', 'quals'])

# Remove pileup metadata tokens from the base-read string.
toremove = re.compile('[<>$*#^]')
pileup['matches'] = pileup['basereads'].apply(lambda x: toremove.sub('', x))

# Convert pileup read characters to canonical bases/symbols.
def pileup_base_counts(row):
    raw = row['matches']
    ref_base = row['_ref'].upper()
    counts = Counter()
    for base in raw:
        if base in '.,':
            counts[ref_base] += 1
        else:
            base = base.upper()
            if base in 'ACGTN':
                counts[base] += 1
    return counts

# Compute Shannon entropy for a distribution of base counts.
def shannon_entropy(counts):
    total = sum(counts.values())
    if total == 0:
        return 0.0
    entropy = 0.0
    for count in counts.values():
        if count == 0:
            continue
        p = count / total
        entropy -= p * math.log2(p)
    return entropy

pileup['base_counts'] = pileup.apply(pileup_base_counts, axis=1)
pileup['total_reads'] = pileup['base_counts'].apply(lambda c: sum(c.values()))
pileup['shannon_entropy'] = pileup['base_counts'].apply(shannon_entropy)

print(pileup[['chrom', 'pos', 'base_counts', 'total_reads', 'shannon_entropy']].head())

# Write bedGraph output: chrom, start, end, entropy.
output_file = 'results/CLIP-let7g-gene.shannon.bedgraph'
with open(output_file, 'w') as out:
    out.write('track type=bedGraph name="let7g_shannon_entropy" description="Shannon entropy at each position"\n')
    for _, row in pileup.iterrows():
        start = int(row['pos']) - 1
        end = int(row['pos'])
        out.write(f"{row['chrom']}\t{start}\t{end}\t{row['shannon_entropy']:.6f}\n")

print(f'BedGraph written to {output_file}')
