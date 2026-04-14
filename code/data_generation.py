import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

def generate_gene_data(n_genes=50, seed=42):
    np.random.seed(seed)

    gene_data = pd.DataFrame({
        "Gene": [f"G{i}" for i in range(1, n_genes+1)],
        "Control": np.random.poisson(lam=50, size=n_genes),
        "Treatment": np.random.poisson(lam=60, size=n_genes)
    })
     # 📊 Title added here
    plt.title("Gene Expression: Control vs Treatment")
    plt.xlabel("Gene Index")
    plt.ylabel("Expression Count")
    plt.plot(gene_data["Control"], label="Control")
    plt.plot(gene_data["Treatment"], label="Treatment")
    plt.legend()
    plt.show()
    return gene_data
