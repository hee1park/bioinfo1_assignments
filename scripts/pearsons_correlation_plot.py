import pandas as pd
from matplotlib import pyplot as plt
import numpy as np
from scipy import stats

# Load read counts data
cnts = pd.read_csv('data-work/read-counts.txt', sep='\t', comment='#', index_col=0)
print(f"Total genes: {len(cnts)}") 


# RNA-seq filtering : exclude transcripts with low read counts in RNA-seq (<30 raw reads)
rna_cols = ['data-work/RNA-control.bam', 'data-work/RNA-siLin28a.bam', 'data-work/RNA-siLuc.bam']
rna_filter = (cnts[rna_cols] >= 30).all(axis=1)

# RNA footprint filtering : exclude transcripts with low ribosome footprints (<80 raw footprint tags in siLuc)
footprint_filter = (cnts['data-work/RPF-siLuc.bam'] >= 80)

cnts = cnts[rna_filter & footprint_filter]

# total counts
CLIP_total = cnts['data-work/CLIP-35L33G.bam'].sum()
RNA_CTRL_total = cnts['data-work/RNA-control.bam'].sum()
RPF_KD_total = cnts['data-work/RPF-siLin28a.bam'].sum()
RNA_KD_total = cnts['data-work/RNA-siLin28a.bam'].sum()
RPF_CTRL_total = cnts['data-work/RPF-siLuc.bam'].sum()
RNA_CTRL_total = cnts['data-work/RNA-siLuc.bam'].sum()
# normalize
cnts["CLIP_norm"] = cnts['data-work/CLIP-35L33G.bam'] / CLIP_total
cnts["RNA_CTRL_norm"] = cnts['data-work/RNA-control.bam'] / RNA_CTRL_total
cnts["RPF_KD_norm"] = cnts['data-work/RPF-siLin28a.bam'] / RPF_KD_total
cnts["RNA_KD_norm"] = cnts['data-work/RNA-siLin28a.bam'] / RNA_KD_total
cnts["RPF_CTRL_norm"] = cnts['data-work/RPF-siLuc.bam'] / RPF_CTRL_total
cnts["RNA_CTRL_norm"] = cnts['data-work/RNA-siLuc.bam'] / RNA_CTRL_total

# Calculate metrics
cnts['clip_enrichment'] = (cnts["CLIP_norm"]) / (cnts["RNA_CTRL_norm"])
cnts['rden_change'] = (cnts["RPF_KD_norm"] / cnts["RNA_KD_norm"]) / (cnts["RPF_CTRL_norm"] / cnts["RNA_CTRL_norm"])

#exclude genes with zero or negative values in either metric before log transformation
valid_data = cnts[(cnts['clip_enrichment'] > 0) & (cnts['rden_change'] > 0)].copy() 
valid_data = valid_data.replace([np.inf, -np.inf], np.nan).dropna()
print(f"Valid genes for correlation: {len(valid_data)}")

# Calculate Pearson correlation
log2_clip = np.log2(valid_data['clip_enrichment'])
log2_rden = np.log2(valid_data['rden_change'])
# log2_rden_centered = valid_data['log2_rden'] - valid_data['log2_rden'].mean()
# r_value, p_value = stats.pearsonr(log2_clip, log2_rden_centered)
r_value, p_value = stats.pearsonr(log2_clip, log2_rden)
print(f"Pearson correlation: r = {r_value:.4f}, p-value = {p_value:.2e}")

# Create plot
fig, ax = plt.subplots(1, 1, figsize=(6, 6))

ax.scatter(log2_clip, log2_rden, alpha=0.5, s=10, c='steelblue', edgecolors='none')

# Add horizontal and vertical reference lines at 0
ax.axhline(y=0, color='gray', linestyle=':', linewidth=0.8, alpha=0.7)
ax.axvline(x=0, color='gray', linestyle=':', linewidth=0.8, alpha=0.7)

# Labels and title
ax.set_xlabel('LIN28A CLIP enrichment (log$_2$)', fontsize=12)
ax.set_ylabel('Ribosome density change upon $\\it{LIN28a}$ knockdown (log$_2$)', fontsize=12)
ax.set_title('CLIP and ribosome footprinting upon $\\it{LIN28a}$ knockdown', fontsize=13)

# Add correlation stats as text
ax.legend(loc='lower right', labels=[f"r = {r_value:.4f}"], fontsize=12, frameon=False)
ax.set_xlim(left=max(log2_clip.min() - 0.5, -6))
ax.set_ylim(bottom=log2_rden.min() - 0.5, top=log2_rden.max() + 0.5)

plt.tight_layout()
plt.savefig('results/pearsons_correlation_plot.png', dpi=150)
print("Figure saved to results/pearsons_correlation_plot.png")