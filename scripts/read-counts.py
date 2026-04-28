import pandas as pd
from matplotlib import pyplot as plt
import numpy as np

cnts = pd.read_csv('data-work/read-counts.txt', sep='\t', comment='#', index_col=0)

cnts['clip_enrichment'] = cnts['data-work/CLIP-35L33G.bam'] / cnts['data-work/RNA-control.bam']
cnts['rden_change'] = (cnts['data-work/RPF-siLin28a.bam'] / cnts['data-work/RNA-siLin28a.bam']) / (cnts['data-work/RPF-siLuc.bam'] / cnts['data-work/RNA-siLuc.bam'])

print(cnts.head())

fig, ax = plt.subplots(1, 1, figsize=(5, 5))
ax.scatter(np.log2(cnts['clip_enrichment']),
           np.log2(cnts['rden_change']))
ax.set_xlabel('log2(CLIP enrichment)')
ax.set_ylabel('log2(ribosome occupancy change)')
ax.set_title('CLIP enrichment vs ribosome occupancy change')

plt.tight_layout()
plt.savefig('results/clip_ribosome_plot.png', dpi=150)
print("Figure saved to results/clip_ribosome_plot.png")


