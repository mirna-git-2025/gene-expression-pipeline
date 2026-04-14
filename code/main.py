from data_generation import generate_gene_data
from analysis import compute_log2fc, select_genes
from sequence_utils import generate_sequences, gc_content
from alignment import pairwise_alignment
from phylo_tree import build_tree, print_tree
from utils import print_gc_contents
import matplotlib.pyplot as plt

# 1. Generate Data
gene_data = generate_gene_data()

# 2. Compute log2FC
gene_data = compute_log2fc(gene_data)

# 3. Select genes
selected_genes = select_genes(gene_data, threshold=0.5)

print("Selected genes:")
print(selected_genes.head())

# 4. Generate sequences
sequences = generate_sequences(selected_genes["Gene"])

# 5. GC Content
print_gc_contents(sequences, gc_content)

gc_values = [gc_content(r) for r in sequences]
gene_names = [r.id for r in sequences]

plt.bar(gene_names, gc_values)
plt.title("GC Content per Gene")
plt.xticks(rotation=90)
plt.show()

# 6. Pairwise alignment
if len(sequences) >= 2:
    alignment = pairwise_alignment(sequences[0].seq, sequences[1].seq)
    print("\nPairwise Alignment:")
    print(alignment)
    print("\n=== Score ===")
    print(alignment.score)

# 7. Phylogenetic tree
if len(sequences) >= 3:
    tree = build_tree(sequences)
    print("\nPhylogenetic Tree:")
    print_tree(tree)

# 8. Save data
gene_data.to_csv("gene_expression_no_stats.csv", index=False)
