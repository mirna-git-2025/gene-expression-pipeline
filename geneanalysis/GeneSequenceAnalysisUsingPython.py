# =========================================
# 🧬 COMPLETE BIOINFORMATICS PIPELINE (WITH OUTPUT DISPLAY)
# =========================================

import pandas as pd
import numpy as np
from Bio import Entrez, SeqIO, Phylo
from Bio.Seq import Seq
from Bio.Align import MultipleSeqAlignment
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor

# ---------- CONFIGURATION ----------
#Entrez.email = "mirna.mokhtar@liu.edu.lb"
Entrez.email = "Aya.alassaad@st.ul.edu.lb"
Entrez.tool = "PhD_Bioinformatics_Pipeline"

# ---------- GLOBAL STORAGE ----------
downloaded_files = []
sequence_lengths = {}
sequence_gc = {}
protein_sequences = {}
motif_counts = {}

# ---------- NCBI SEARCH ----------
def search_ncbi(gene, organism, retmax=5):
    print("\n🔎 Searching NCBI...")
    handle = Entrez.esearch(
        db="nucleotide",
        term=f"{gene}[Gene] AND {organism}[Organism]",
        retmax=retmax
    )
    record = Entrez.read(handle)
    handle.close()

    ids = record["IdList"]

    if not ids:
        print("❌ No results found.")
    else:
        print(f"✅ Found {len(ids)} sequences:")
        print(ids)

    return ids

# ---------- FETCH FASTA ----------
def fetch_fasta(ids, output_file):
    print("\n📥 Fetching FASTA sequences...")
    handle = Entrez.efetch(
        db="nucleotide",
        id=",".join(ids),
        rettype="fasta",
        retmode="text"
    )
    data = handle.read()
    handle.close()

    with open(output_file, "w") as f:
        f.write(data)

    downloaded_files.append(output_file)
    print(f"✅ FASTA saved: {output_file}")

# ---------- PARSE FASTA + METRICS ----------
def parse_fasta(file):
    print("\n📊 Parsing FASTA & Computing Metrics...")

    for record in SeqIO.parse(file, "fasta"):
        seq_id = record.id
        seq = str(record.seq)
        length = len(seq)
        gc = (seq.count("G") + seq.count("C")) / length * 100 if length > 0 else 0

        sequence_lengths[seq_id] = length
        sequence_gc[seq_id] = round(gc, 2)

    df = pd.DataFrame({
        "Length": sequence_lengths,
        "GC%": sequence_gc
    })

    print("\n🧾 Sequence Metrics:")
    print(df)

# ---------- N50 CALCULATION ----------
def calculate_n50():
    lengths = sorted(sequence_lengths.values(), reverse=True)
    total = sum(lengths)
    cumulative = 0

    for l in lengths:
        cumulative += l
        if cumulative >= total / 2:
            return l
    return 0

# ---------- PROTEIN TRANSLATION ----------
def translate_sequences(file):
    print("\n🧬 Translating DNA → Protein...")

    for record in SeqIO.parse(file, "fasta"):
        seq = record.seq
        trimmed_seq = seq[:len(seq)//3*3]
        protein = str(trimmed_seq.translate(to_stop=True))
        protein_sequences[record.id] = protein

    df = pd.DataFrame.from_dict(protein_sequences, orient="index", columns=["Protein"])
    print("\n🧾 Protein Sequences:")
    print(df.head())

# ---------- MOTIF ANALYSIS ----------
def find_motif(file, motif="ATG"):
    print(f"\n🔍 Searching for motif: {motif}")

    for record in SeqIO.parse(file, "fasta"):
        motif_counts[record.id] = str(record.seq).count(motif)

    df = pd.DataFrame.from_dict(motif_counts, orient="index", columns=["Motif Count"])
    print("\n🧾 Motif Counts:")
    print(df)

# ---------- MULTIPLE SEQUENCE ALIGNMENT ----------
def build_alignment(file):
    print("\n🧬 Building Alignment (equal length padding)...")

    records = list(SeqIO.parse(file, "fasta"))
    max_len = max(len(r.seq) for r in records)

    for r in records:
        r.seq = r.seq + "-" * (max_len - len(r.seq))

    alignment = MultipleSeqAlignment(records)

    print("✅ Alignment Preview:")
    for rec in alignment[:2]:
        print(rec.seq[:60], "...")

    return alignment

# ---------- PHYLOGENETIC TREE ----------
def build_tree(alignment):
    print("\n🌳 Building Phylogenetic Tree...")

    calculator = DistanceCalculator('identity')
    dist_matrix = calculator.get_distance(alignment)

    print("\n📏 Distance Matrix:")
    print(dist_matrix)

    constructor = DistanceTreeConstructor()
    tree = constructor.nj(dist_matrix)

    print("\n🌳 Phylogenetic Tree (ASCII):")
    Phylo.draw_ascii(tree)

    return tree

# ---------- SEQUENCE STATISTICS ----------
def sequence_statistics():
    print("\n📊 Sequence Statistics:")

    lengths = list(sequence_lengths.values())

    print("Mean length:", np.mean(lengths))
    print("Max length:", np.max(lengths))
    print("Min length:", np.min(lengths))
    print("N50:", calculate_n50())

# ---------- EXPORT RESULTS ----------
def export_results():
    print("\n💾 Exporting Results...")

    pd.DataFrame.from_dict(sequence_lengths, orient="index", columns=["Length"]).to_csv("lengths.csv")
    pd.DataFrame.from_dict(sequence_gc, orient="index", columns=["GC%"]).to_csv("gc_content.csv")
    pd.DataFrame.from_dict(motif_counts, orient="index", columns=["Motif_Count"]).to_csv("motif_counts.csv")
    pd.DataFrame.from_dict(protein_sequences, orient="index", columns=["Protein"]).to_csv("proteins.csv")

    print("✅ Files generated:")
    print(" - lengths.csv")
    print(" - gc_content.csv")
    print(" - motif_counts.csv")
    print(" - proteins.csv")

# ---------- MAIN PIPELINE ----------
def main():
    gene = "BRCA1"
    organism = "Homo sapiens"
    output_file = "sequences.fasta"

    ids = search_ncbi(gene, organism)
    if not ids:
        return

    fetch_fasta(ids, output_file)
    parse_fasta(output_file)
    translate_sequences(output_file)
    find_motif(output_file, motif="ATG")

    alignment = build_alignment(output_file)
    tree = build_tree(alignment)

    sequence_statistics()
    export_results()

    print("\n📂 Downloaded Files:", downloaded_files)

# ---------- RUN ----------
if __name__ == "__main__":
    main()
