import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path

# Check input argument

input_file = "results/activity_ratios_set12.csv"
# Read CSV file
df = pd.read_csv(input_file)

output_file = Path(input_file).stem + "_scatter.png"

# Create scatter plot
plt.figure(figsize=(10, 5))

plt.scatter(
    range(len(df)),
    df["log2(exp/nc)"],
)

# Add horizontal reference line at 0 and -1
plt.axhline(0, linestyle="-", color="black", linewidth=0.5)
plt.axhline(-1, linestyle="--", color="gray")

# Label points with strain names (optional)
for i, strain in enumerate(df["Strain"]):
    plt.text(
        i,
        df["log2(exp/nc)"].iloc[i],
        strain,
        fontsize=5,
        rotation=90,
        ha='center',
        va='bottom'
    )

# Axis labels and title
plt.xlabel("Isolates")
plt.ylabel("log2(exp/nc)")
plt.title("Scatter Plot of log2(exp/nc)")

# Improve layout
plt.tight_layout()

# Save figure
plt.savefig(output_file, dpi=300)
