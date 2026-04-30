import pandas as pd
from matplotlib import pyplot as plt
import numpy as np
from scipy import stats
import ssl

# Load read counts data
cnts = pd.read_csv('data-work/read-counts.txt', sep='\t', comment='#', index_col=0)
print(f"Total genes: {len(cnts)}") 

# Load mouse localization data
ssl._create_default_https_context = ssl._create_unverified_context
mouselocal = pd.read_csv('https://hyeshik.qbio.io/binfo/mouselocalization-20210507.txt', sep='\t')
print(f"Localization data: {len(mouselocal)} entries")
print(f"Localization columns: {mouselocal.columns.tolist()}") 


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

# Merge with localization data
valid_data = valid_data.reset_index()
# Strip version number from Geneid for matching
valid_data['gene_id_stripped'] = valid_data['Geneid'].str.split('.').str[0]
valid_data = valid_data.merge(mouselocal, left_on='gene_id_stripped', right_on='gene_id', how='left')
print(f"Genes with localization info: {valid_data['type'].notna().sum()}")

# Calculate Pearson correlation
log2_clip = np.log2(valid_data['clip_enrichment'])
log2_rden = np.log2(valid_data['rden_change'])
r_value, p_value = stats.pearsonr(log2_clip, log2_rden)
print(f"Pearson correlation: r = {r_value:.4f}, p-value = {p_value:.2e}")

# Create plot
fig, ax = plt.subplots(1, 1, figsize=(8, 6))

# Color by localization if available
localizations = valid_data['type'].dropna().unique()
colors = plt.cm.tab20(np.linspace(0, 1, len(localizations)))
color_map = {loc: colors[i] for i, loc in enumerate(localizations)}

# Plot with colors based on localization
for loc in localizations:
    mask = valid_data['type'] == loc
    ax.scatter(log2_clip[mask], log2_rden[mask], alpha=0.4, s=15, 
               c=[color_map[loc]], edgecolors='none', label=loc)

# Plot genes without localization info in gray
no_loc = valid_data['type'].isna()
# ax.scatter(log2_clip[no_loc], log2_rden[no_loc], alpha=0, s=10, 
#            c='lightgray', edgecolors='none', label='Unknown')

# Add horizontal and vertical reference lines at 0
ax.axhline(y=0, color='gray', linestyle=':', linewidth=0.8, alpha=0.7)
ax.axvline(x=0, color='gray', linestyle=':', linewidth=0.8, alpha=0.7)

# Labels and title
ax.set_xlabel('LIN28A CLIP enrichment (log$_2$)', fontsize=12)
ax.set_ylabel('Ribosome density change upon $\\it{LIN28a}$ knockdown (log$_2$)', fontsize=12)
ax.set_title('CLIP and ribosome footprinting upon $\\it{LIN28a}$ knockdown', fontsize=13)

# Add correlation stats as text
textstr = f'Pearson r = {r_value:.4f}\np-value = {p_value:.2e}\nn = {len(valid_data)}'
props = dict(boxstyle='round', facecolor='wheat', alpha=0.8)
ax.text(0.05, 0.95, textstr, transform=ax.transAxes, fontsize=10,
        verticalalignment='top', bbox=props)

# Limit legend to top 10 localizations to avoid clutter
handles, labels = ax.get_legend_handles_labels()
loc_counts = valid_data['type'].value_counts()
top_locs = loc_counts.head(10).index.tolist()
legend_handles = []
legend_labels = []
for h, l in zip(handles, labels):
    if l in top_locs or l == 'Unknown':
        legend_handles.append(h)
        legend_labels.append(l)
ax.legend(legend_handles, legend_labels, loc='center left', bbox_to_anchor=(1, 0.5), 
          fontsize=8, title='Localization')

ax.set_xlim(left=max(log2_clip.min() - 0.5, -6))
ax.set_ylim(bottom=log2_rden.min() - 0.5, top=log2_rden.max() + 0.5)

plt.tight_layout()
plt.savefig('results/protein_localization_plot_without_unknown.png', dpi=150)
print("Figure saved to results/protein_localization_plot_without_unknown.png")