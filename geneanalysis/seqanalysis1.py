from Bio import SeqIO
from Bio.Seq import Seq
import pandas as pd
import matplotlib.pyplot as plt

# =========================
# INPUT FILE
# =========================
fasta_file = "seq1.fasta"
motif = "ATG"

data = []

# =========================
# READ FASTA
# =========================
for record in SeqIO.parse(fasta_file, "fasta"):
    seq_id = record.id
    seq = str(record.seq).upper()

    length = len(seq)
   

    # GC content
    gc = (seq.count("G") + seq.count("C")) / length * 100

    # Protein translation
    protein = str(Seq(seq).translate(to_stop=True))

    # Motif count
    motif_count = seq.count(motif)

    data.append([seq_id, length, round(gc, 2), motif_count, protein])

# =========================
# DATAFRAME
# =========================
df = pd.DataFrame(data, columns=["ID", "Length", "GC%", "Motif_Count", "Protein"])

print("\n🧾 Results:")
print(df)

# =========================
# N50 FUNCTION
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

# =========================
# SAVE RESULTS
# =========================
df.to_csv("results.csv", index=False)

# =========================================
# 📊 1. LENGTH DISTRIBUTION
# =========================================
plt.figure(figsize=(7,4))
plt.hist(df["Length"], bins=10, color="skyblue", edgecolor="black")
plt.title("Sequence Length Distribution")
plt.xlabel("Length")
plt.ylabel("Count")
plt.tight_layout()
plt.show()

# =========================================
# 📊 2. GC% TREND LINE (BEST CHOICE)
# =========================================
plt.figure(figsize=(10,5))

plt.plot(
    df["ID"],
    df["GC%"],
    marker="o",
    linestyle="-",
    color="seagreen"
)

plt.title("GC Content Trend Across Sequences")
plt.xlabel("Sequence ID")
plt.ylabel("GC %")

plt.xticks(rotation=45, ha="right")
plt.ylim(0, 100)
plt.grid(alpha=0.3)

plt.tight_layout()
plt.show()

# =========================================
# 📊 3. LENGTH vs GC
# =========================================
plt.figure(figsize=(6,4))
plt.scatter(df["Length"], df["GC%"], color="purple")
plt.title("Length vs GC%")
plt.xlabel("Length")
plt.ylabel("GC %")
plt.grid(alpha=0.3)
plt.tight_layout()
plt.show()

# =========================================
# 📊 4. TOP SEQUENCES BY LENGTH
# =========================================
top = df.sort_values(by="Length", ascending=False).head(5)

plt.figure(figsize=(7,4))
plt.bar(top["ID"], top["Length"], color="orange")
plt.title("Top 5 Long Sequences")
plt.xlabel("Sequence ID")
plt.ylabel("Length")
plt.xticks(rotation=45)
plt.tight_layout()
plt.show()

print("\n✅ Analysis complete. File saved: results.csv")
