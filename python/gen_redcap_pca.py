import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# Load your newly generated integrated matrix
df = pd.read_csv("pca/pca_with_metadata.csv")

# Set up the plotting canvas
plt.figure(figsize=(10, 8), dpi=300)

# Define what trait to color the points by (e.g., Ventilatory Support or Sex)
# Replace this column string name with any of your 199 matched clinical features
color_trait = "Did or does the patient receive ventilatory support?"

# Create the scatter plot
sns.scatterplot(
    data=df,
    x="PC1",
    y="PC2",
    hue=color_trait,
    palette="viridis",
    alpha=0.7,
    edgecolor=None,
    s=25
)

# Style the chart labels using your explained variance outputs
plt.title(f"Transcriptomic PCA Alignment colored by:\n{color_trait}", fontsize=14, weight='bold')
plt.xlabel("PC1 (37.62% Var Explained)", fontsize=12)
plt.ylabel("PC2 (10.38% Var Explained)", fontsize=12)
plt.legend(title=color_trait, bbox_to_anchor=(1.05, 1), loc='upper left')
plt.tight_layout()

# Save the visualization to disk
output_plot = "pca/pca_transcript_plot.png"
plt.savefig(output_plot, bbox_inches='tight')
print(f"PCA Scatter Plot successfully saved to: {output_plot}")
