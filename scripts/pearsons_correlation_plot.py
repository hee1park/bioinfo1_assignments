import pandas as pd
from matplotlib import pyplot as plt
import numpy as np
from scipy import stats

# Load read counts data
cnts = pd.read_csv('data-work/read-counts.txt', sep='\t', comment='#', index_col=0)

# Define BAM columns
bam_cols = ['data-work/CLIP-35L33G.bam', 'data-work/RNA-control.bam',
            'data-work/RNA-siLin28a.bam', 'data-work/RNA-siLuc.bam',
            'data-work/RPF-siLin28a.bam', 'data-work/RPF-siLuc.bam']

# Exclude genes with 0 expression in all BAM files
cnts['total_expr'] = cnts[bam_cols].sum(axis=1)
cnts = cnts[cnts['total_expr'] > 0]
print(f"Genes after filtering zero expression: {len(cnts)}")

# Calculate metrics
cnts['clip_enrichment'] = cnts['data-work/CLIP-35L33G.bam'] / cnts['data-work/RNA-control.bam']
cnts['rden_change'] = (cnts['data-work/RPF-siLin28a.bam'] / cnts['data-work/RNA-siLin28a.bam']) / (cnts['data-work/RPF-siLuc.bam'] / cnts['data-work/RNA-siLuc.bam'])

# Filter out zeros and infinities for correlation analysis
valid_data = cnts[(cnts['clip_enrichment'] > 0) & (cnts['rden_change'] > 0)].copy()
valid_data = valid_data.replace([np.inf, -np.inf], np.nan).dropna()

print(f"Total genes: {len(cnts)}")
print(f"Valid genes for correlation: {len(valid_data)}")

# Calculate Pearson correlation
log2_clip = np.log2(valid_data['clip_enrichment'])
log2_rden = np.log2(valid_data['rden_change'])

r_value, p_value = stats.pearsonr(log2_clip, log2_rden)

print(f"Pearson correlation: r = {r_value:.4f}, p-value = {p_value:.2e}")

# Create plot
fig, ax = plt.subplots(1, 1, figsize=(6, 6))

ax.scatter(log2_clip, log2_rden, alpha=0.5, s=20, c='steelblue', edgecolors='none')

# Add regression line
slope, intercept = np.polyfit(log2_clip, log2_rden, 1)
x_line = np.linspace(log2_clip.min(), log2_clip.max(), 100)
y_line = slope * x_line + intercept
ax.plot(x_line, y_line, color='red', linewidth=2, linestyle='--', label=f'Linear fit (r={r_value:.3f})')

# Add horizontal and vertical reference lines at 0
ax.axhline(y=0, color='gray', linestyle=':', linewidth=0.8, alpha=0.7)
ax.axvline(x=0, color='gray', linestyle=':', linewidth=0.8, alpha=0.7)

# Labels and title
ax.set_xlabel('log2(CLIP enrichment)', fontsize=12)
ax.set_ylabel('log2(ribosome occupancy change)', fontsize=12)
ax.set_title('LIN28A interaction vs ribosome occupancy change\n(Pearson correlation)', fontsize=13)

# Add correlation stats as text
textstr = f'Pearson r = {r_value:.4f}\np-value = {p_value:.2e}\nn = {len(valid_data)}'
props = dict(boxstyle='round', facecolor='wheat', alpha=0.8)
ax.text(0.05, 0.95, textstr, transform=ax.transAxes, fontsize=10,
        verticalalignment='top', bbox=props)

ax.legend(loc='lower right')
ax.set_xlim(left=max(log2_clip.min() - 0.5, -5))
ax.set_ylim(bottom=min(log2_rden.min() - 0.5, -5), top=max(log2_rden.max() + 0.5, 5))

plt.tight_layout()
plt.savefig('results/pearsons_correlation_plot.png', dpi=150)
print("Figure saved to results/pearsons_correlation_plot.png")