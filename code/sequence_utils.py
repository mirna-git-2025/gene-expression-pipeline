import numpy as np
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def generate_sequences(genes, length=50):
    sequences = []
    for gene in genes:
        seq = ''.join(np.random.choice(list('ATGC'), length))
        sequences.append(SeqRecord(Seq(seq), id=gene))
    return sequences

def gc_content(seq_record):
    seq = seq_record.seq
    return (seq.count("G") + seq.count("C")) / len(seq) * 100
