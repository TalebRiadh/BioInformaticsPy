# DNA Toolkit file
import collections
from structures import *


# Check the sequence to make sure it is a DNA string
def validateSeq(dna_seq):
    tmpseq = dna_seq.upper()
    for nuc in tmpseq:
        if nuc not in Nucleotides:
            return False
    return tmpseq


def countNucFrequency(seq):
    tmpFreqDict = {"A": 0, "C": 0, "G": 0, "T": 0}
    for nuc in seq:
        tmpFreqDict[nuc] += 1
    return tmpFreqDict
    # return dict(collections.Counter(seq))


def transcription(seq):
    # DNA -> RNA Transcription
    return seq.replace("T", "U")


def reverse_complement(seq):
   # return ''.join([DNA_ReverseComplement[nuc] for nuc in seq])[::-1]

   mapping = str.maketrans('ATCG', 'TAGC')
   return seq.translate(mapping)[::1]


def gc_content(seq):
    return round((seq.count('C') + seq.count('G') / len(seq) * 100))


def gc_content_subsec(seq, k=20):
    """GC Content in a DNA/RNA sub-sequence length k, k =20 by default"""
    res = []
    for i in range(0, len(seq) - k + 1, k):
        subseq = seq[i:i + k]
        res.append(gc_content(subseq))
    return res


def translate_seq(seq, init_pos=0):
    """Translate a DNA sequence into an aminoacid sequence"""
    return [DNA_Codons[seq[pos:pos + 3]] for pos in range(init_pos, len(seq) -2, 3)]


def codon_usage(seq, aminoacid):
    """Provides the frequency of each codon encoding a given aminoacid in a DNA  sequece"""
    tmpList = []
    for i in range(0, len(seq) -2, 3):
        if DNA_Codons[seq[i:i +3]] == aminoacid:
            tmpList.append(seq[i:i + 3])
    freqDict = dict(collections.Counter(tmpList))
    totalWight = sum(freqDict.values())
    for seq in freqDict:
        freqDict[seq] = round(freqDict[seq] / totalWight, 2)
        return freqDict


def gen_reading_frames(seq):
    """Generate the six reading frames of a DNA sequence, including the reverse complement"""
    frames = []
    frames.append(translate_seq(seq, 0))
    frames.append(translate_seq(seq, 1))
    frames.append(translate_seq(seq, 2))
    frames.append(translate_seq(reverse_complement(seq), 0))
    frames.append(translate_seq(reverse_complement(seq), 1))
    frames.append(translate_seq(reverse_complement(seq), 2))
    return frames


def proteins_from_rf(aa_seq):
    """Compute all possible proteins in an aminoacid seq and return a list of possible proteins"""
    current_prot = []
    proteins = []
    for aa in aa_seq:
        if aa == "_":
            # STOP accumulating amino acids if _ - STOP was found
            if current_prot:
                for p in current_prot:
                    proteins.append(p)
                current_prot = []
        else:
            # START accumulating amino acids if M - START was found
            if aa == "M":
                current_prot.append("")
            for i in range(len(current_prot)):
                current_prot[i] += aa
    return proteins


def all_proteins_from_orfs(seq, startReadPos=0, endReadPos=0, ordered=False):
    """Compute all possible proteins for all open reading frames"""
    if endReadPos > startReadPos:
        rfs = gen_reading_frames(seq[startRead: endRead])
    else:
        rfs = gen_reading_frames(seq)
    res = []
    for rf in rfs:
        prots = proteins_from_rf(rf)
        for p in prots:
            res.append(p)
    if ordered:
        return sorted(res, key=len, reverse=True)
    return res