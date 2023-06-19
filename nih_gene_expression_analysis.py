"""
NIH Data Analysis Program

Author: Kevin Mastascusa
Date: 6/19/2023

This program allows you to analyze gene expression data from a GEO dataset.
Provide the GEO dataset ID as input to the 'analyze_gene_expression' function for analysis and visualization.

Example Usage:
    Enter the GEO dataset ID: GSE2553

Output:
- Heatmap Visualization:
   - The heatmap represents gene expression values for different samples and genes.
   - Each cell corresponds to the expression level of a specific gene in a particular sample.
   - Cooler colors (e.g., blue) indicate lower expression, while warmer colors (e.g., red) indicate higher expression.

TODO:
- Implement additional analysis and visualization functionalities.
- Customize the plot appearance and enhance the analysis based on specific requirements.
"""

import GEOparse
import matplotlib.pyplot as plt
import seaborn as sns


def analyze_gene_expression(gse_id: str):
    """
    Analyze gene expression data from a GEO dataset.

    :param gse_id: The GEO dataset ID.
    """
    # Fetch the GEO dataset using the provided ID
    gse = GEOparse.get_GEO(geo=gse_id)
    print(f"Analyzing gene expression data from GEO dataset {gse_id}...")

    # Perform data analysis and visualization
    # Example: Plot gene expression heatmap
    expression_values = gse.pivot_samples('VALUE')
    sns.heatmap(expression_values, cmap='coolwarm')
    plt.xlabel("Samples")
    plt.ylabel("Genes")
    plt.title("Gene Expression Heatmap")
    plt.show()


if __name__ == "__main__":
    # Prompt the user to enter the GEO dataset ID for analysis
    dataset_id = input("Enter the GEO dataset ID: ")
    analyze_gene_expression(dataset_id)
