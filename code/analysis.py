import numpy as np
import matplotlib.pyplot as plt
def compute_log2fc(df):
    df["log2FC"] = np.log2((df["Treatment"] + 1) / (df["Control"] + 1))
    plt.bar(df["Gene"], df["log2FC"])
    plt.xticks(rotation=90)
    # 📊 Add title here
    plt.title("Log2 Fold Change (Treatment vs Control)")

    plt.xlabel("Genes")
    plt.ylabel("log2FC")
    plt.show()
    return df

def select_genes(df, threshold=0.5):
    return df[abs(df["log2FC"]) > threshold]
