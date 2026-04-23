import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

from Bio import SeqIO
from Bio.Align import MultipleSeqAlignment
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
from Bio import Phylo

# ==============================
# 1. LOAD GROWTH DATA
# ==============================

df = pd.read_csv("growth_data.csv")
df.columns = [c.strip() for c in df.columns]

# fix column name if needed
if "Time (hours)" in df.columns:
    df.rename(columns={"Time (hours)": "Time"}, inplace=True)

print("\n=== GROWTH DATA ===")
print(df)

# ==============================
# 2. GROWTH RATE (NUMPY)
# ==============================

time = df["Time"].to_numpy(float)
od = df["OD600"].to_numpy(float)

growth_rate = np.diff(od) / np.diff(time)

print("\n=== GROWTH RATES ===")
print(growth_rate)

# ==============================
# 3. LOAD SEQUENCES
# ==============================

records = list(SeqIO.parse("sequence.fasta", "fasta"))

print("\n=== RAW SEQUENCES ===")
for r in records:
    print(r.id, r.seq)

# ==============================
# 4. ALIGNMENT PREPARATION (FIX LENGTHS)
# ==============================

max_len = max(len(r.seq) for r in records)

for r in records:
    r.seq = r.seq + ("-" * (max_len - len(r.seq)))

print("\n=== PADDED SEQUENCES (WITH GAPS) ===")
for r in records:
    print(r.id, r.seq)

# ==============================
# 5. CREATE MULTIPLE SEQUENCE ALIGNMENT
# ==============================

alignment = MultipleSeqAlignment(records)

# ==============================
# 6. DISTANCE MATRIX (CORRECT METHOD)
# ==============================

calculator = DistanceCalculator("identity")
dm = calculator.get_distance(alignment)

print("\n=== DISTANCE MATRIX ===")
print(dm)

# ==============================
# 7. PHYLOGENETIC TREE
# ==============================

constructor = DistanceTreeConstructor()
tree = constructor.nj(dm)

print("\n=== PHYLOGENETIC TREE (ASCII) ===")
Phylo.draw_ascii(tree)

# Optional graphical tree
Phylo.draw(tree)
