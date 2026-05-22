import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

from Bio import SeqIO, Phylo
from Bio.Seq import Seq
from Bio.Align import MultipleSeqAlignment, PairwiseAligner
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor

# =========================
# LOAD CSV + LOG2FC
# =========================
gene_data = pd.read_csv("gene_data.csv")

def compute_log2fc(df):
    df["log2FC"] = np.log2((df["Treatment"] + 1) / (df["Control"] + 1))

    plt.figure()
    plt.bar(df["Gene"], df["log2FC"])
    plt.xticks(rotation=90)
    plt.title("Log2 Fold Change (Treatment vs Control)")
    plt.xlabel("Genes")
    plt.ylabel("log2FC")
    plt.tight_layout()
    plt.savefig("log2fc.png")
    plt.show()

    return df

def select_genes(df, threshold=0.5):
    return df[abs(df["log2FC"]) > threshold]

gene_data = compute_log2fc(gene_data)
selected_genes = select_genes(gene_data)

print("\nSelected genes:")
print(selected_genes.head())

# =========================
# FASTA INPUT
# =========================
fasta_file = "selected_genes.fasta"
motif = "ATG"

sequences = list(SeqIO.parse(fasta_file, "fasta"))

print([len(r.seq) for r in sequences])
data = []

# =========================
# SEQUENCE ANALYSIS
# =========================
for record in sequences:
    seq_id = record.id
    seq = str(record.seq).upper()

    length = len(seq)

    gc = (seq.count("G") + seq.count("C")) / length * 100
    rna = seq.replace("T", "U")
    protein = str(Seq(seq).translate(to_stop=True))
    motif_count = seq.count(motif)

    data.append([seq_id, length, round(gc, 2), motif_count, rna, protein])

df = pd.DataFrame(
    data,
    columns=["ID", "Length", "GC%", "Motif_Count", "RNA", "Protein"]
)

print("\n🧾 RESULTS:")
print(df)

# =========================
# N50
# =========================
def calculate_n50(lengths):
    lengths = sorted(lengths, reverse=True)
    total = sum(lengths)
    cum = 0
    for l in lengths:
        cum += l
        if cum >= total / 2:
            return l

n50 = calculate_n50(df["Length"])
print("\n📊 N50:", n50)

df.to_csv("results.csv", index=False)

# =========================
# 1. PAIRWISE ALIGNMENT + PNG
# =========================
if len(sequences) >= 2:
    aligner = PairwiseAligner()
    alignment = aligner.align(sequences[0].seq, sequences[1].seq)[0]

    print("\n🔬 Pairwise Alignment:")
    print(alignment)

    # SAVE TEXT
    with open("pairwise_alignment.txt", "w") as f:
        f.write(str(alignment))

    # 🔥 SAVE ALIGNMENT AS IMAGE
    plt.figure(figsize=(10,4))
    plt.text(0.01, 0.5, str(alignment), fontsize=10, family="monospace")
    plt.axis("off")
    plt.title("Pairwise Alignment")
    plt.tight_layout()
    plt.savefig("pairwise_alignment.png")
    plt.show()

# =========================
# 2. DISTANCE MATRIX + NUMBERS
# =========================
alignment_all = MultipleSeqAlignment(sequences)

calculator = DistanceCalculator("identity")
dist_matrix = calculator.get_distance(alignment_all)

print("\n📊 Distance Matrix:")
print(dist_matrix)

dist_df = pd.DataFrame(
    np.array(dist_matrix),
    index=dist_matrix.names,
    columns=dist_matrix.names
)

dist_df.to_csv("distance_matrix.csv")

# 🔥 HEATMAP WITH NUMBERS
plt.figure(figsize=(7,6))
plt.imshow(dist_df.values)

plt.title("Distance Matrix Heatmap")

plt.xticks(range(len(dist_df.columns)), dist_df.columns, rotation=90)
plt.yticks(range(len(dist_df.index)), dist_df.index)

# Add numbers inside cells
for i in range(len(dist_df)):
    for j in range(len(dist_df)):
        plt.text(j, i, f"{dist_df.iloc[i, j]:.2f}",
                 ha="center", va="center", color="black", fontsize=8)

plt.colorbar()
plt.tight_layout()
plt.savefig("distance_matrix.png")
plt.show()

# =========================
# 3. PHYLOGENETIC TREE
# =========================
constructor = DistanceTreeConstructor()
tree = constructor.nj(dist_matrix)

print("\n🌳 PHYLOGENETIC TREE:")
Phylo.draw_ascii(tree)

plt.figure(figsize=(8,6))
Phylo.draw(tree, do_show=False)
plt.savefig("phylogenetic_tree.png")
plt.show()

# =========================
# 4. PLOTS
# =========================

# Length distribution
plt.figure()
plt.hist(df["Length"], bins=10)
plt.title("Sequence Length Distribution")
plt.savefig("length_distribution.png")
plt.show()

# GC trend
plt.figure()
plt.plot(df["ID"], df["GC%"], marker="o")
plt.xticks(rotation=45)
plt.title("GC Content Trend")
plt.savefig("gc_trend.png")
plt.show()

# Length vs GC
plt.figure()
plt.scatter(df["Length"], df["GC%"])
plt.title("Length vs GC%")
plt.savefig("length_vs_gc.png")
plt.show()

# Top sequences
top = df.sort_values(by="Length", ascending=False).head(5)

plt.figure()
plt.bar(top["ID"], top["Length"])
plt.title("Top 5 Sequences")
plt.xticks(rotation=45)
plt.savefig("top_sequences.png")
plt.show()

# =========================
# FINAL MESSAGE
# =========================
print("\n✅ ANALYSIS COMPLETE")
print("Generated files:")
print("- results.csv")
print("- pairwise_alignment.png")
print("- distance_matrix.png (with values)")
print("- phylogenetic_tree.png")
print("- log2fc.png")
print("- gc_trend.png")
print("- length_distribution.png")
print("- length_vs_gc.png")
print("- top_sequences.png")
