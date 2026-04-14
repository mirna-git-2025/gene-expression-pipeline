from Bio.Align import MultipleSeqAlignment
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
from Bio import Phylo

def build_tree(sequences, n=5):
    alignment = MultipleSeqAlignment(sequences[:n])
    calculator = DistanceCalculator('identity')
    dist_matrix = calculator.get_distance(alignment)
     # 📊 PRINT DISTANCE MATRIX
    print("=== Distance Matrix ===")
    print(dist_matrix)
    constructor = DistanceTreeConstructor()
    tree = constructor.nj(dist_matrix)

    return tree

def print_tree(tree):
    Phylo.draw_ascii(tree)
