import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

from Bio import SeqIO
from Bio.Align import MultipleSeqAlignment
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
from Bio import Phylo

# ==============================
# 1. LOAD METADATA (CSV)
# ==============================

df = pd.read_csv("data.csv")

# Clean ALL column names
df.columns = df.columns.str.strip().str.lower()

print("\n=== METADATA ===")
print(df)

# ==============================
# 2. COMPUTE CONCENTRATION (g/cm^3)
# ==============================

df["concentration"] = df["mass"] / df["volume"]

print("\n=== CONCENTRATION (g/cm^3) ===")
print(df[["species", "concentration"]])

# ==============================
# 3. LOAD FASTA SEQUENCES
# ==============================

records = list(SeqIO.parse("projectsequences.fasta", "fasta"))

print("\n=== RAW SEQUENCES ===")
for r in records:
    print(r.id, r.seq)



# ==============================
# 4. MULTIPLE SEQUENCE ALIGNMENT
# ==============================

alignment = MultipleSeqAlignment(records)

# ==============================
#5. DISTANCE MATRIX (IDENTITY)
# ==============================

calculator = DistanceCalculator("identity")
dm = calculator.get_distance(alignment)

print("\n=== DISTANCE MATRIX ===")
print(dm)

# ==============================
# 6. PHYLOGENETIC TREE (NJ)
# ==============================

constructor = DistanceTreeConstructor()
tree = constructor.nj(dm)

print("\n=== PHYLOGENETIC TREE (ASCII) ===")
Phylo.draw_ascii(tree)

# ==============================
# 7. GRAPHICAL TREE
# ==============================

Phylo.draw(tree)

# ==============================
# 8. BAR CHART (CONCENTRATION)
# ==============================

plt.figure(figsize=(10, 5))
plt.bar(df["species"], df["concentration"], color="skyblue")

plt.xticks(rotation=45, ha="right")
plt.ylabel("Concentration (g/cm³)")
plt.title("Concentration per Organism")

plt.tight_layout()
plt.show()

# ==============================
# 9. SCATTER PLOT (MASS vs VOLUME)
# ==============================

plt.figure(figsize=(6, 5))
plt.scatter(df["volume"], df["mass"], color="green")

for i, species in enumerate(df["species"]):
    plt.text(df["volume"][i], df["mass"][i], species)

plt.xlabel("Volume (cm³)")
plt.ylabel("Mass (g)")
plt.title("Mass vs Volume")

plt.grid(True)
plt.show()

# ==============================
# 10. PIE CHART (CONCENTRATION DISTRIBUTION)
# ==============================

plt.figure(figsize=(6, 6))
plt.pie(df["concentration"], labels=df["species"], autopct="%1.1f%%")

plt.title("Relative Concentration Distribution")
plt.show()
