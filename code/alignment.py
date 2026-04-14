from Bio.Align import PairwiseAligner

def pairwise_alignment(seq1, seq2):
    aligner = PairwiseAligner()
    alignment = aligner.align(seq1, seq2)[0]
    return alignment
