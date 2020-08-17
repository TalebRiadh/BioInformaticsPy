import collections

from bio_structs import DNA_Codons, DNA_Nucleotides
import random


# DNA sequence class, Default value: ATCG, DNA, No label
class bio_seq:
    def __init__(self, seq="ATCG", seq_type="DNA", label="No_Label"):

        self.seq = seq.upper()
        self.label = label
        self.seq_type = seq_type
        self.is_valid = self.__validate()
        assert self.is_valid, f"Provided data does not seem to be a correct {self.seq_type}"

    # Check the sequence to make sure it is a valid DNA string
    def __validate(self):
        return set(DNA_Nucleotides).issuperset(self.seq)

    # Return sequence type
    def get_seq_biotype(self):
        return self.seq_type

    # Returns 4 strings, full sequence information
    def get_seq_info(self):
        return f"[label]: {self.label}\n[Sequence]: {self.seq}\n[Biotype]: {self.seq_type}"

    # Generate a random DNA sequence, provided the length
    def generate_rnd_seq(self, length=10, seq_type="DNA"):
        seq = ''.join([random.choice(DNA_Nucleotides)
                       for x in  range(length)])
        self.__init__(seq, seq_type, "Randomly generated sequence")

    # Count nucleotides in a given sequence, return a dictionary
    def nucleotide_frequency(self):
        return dict(collections.Counter(self.seq))

    # DNA -> RNA Transcription
    # return a copier of the sequence
    def transcription(self):
        return self.seq.replace("T", "U")

    # return ''.join([DNA_ReverseComplement[nuc] for nuc in seq])[::-1]
    def reverse_complement(self):
        mapping = str.maketrans('ATCG', 'TAGC')
        return self.seq.translate(mapping)[::1]


    def gc_content(self):
        return round((self.seq.count('C') + self.seq.count('G') / len(self.seq) * 100))

    # GC Content in a DNA/RNA sub_sequence length k, k=20 by default
    def gc_content_subsec(self, k=20):
        res = []
        for i in range(0, len(self.seq) - k + 1, k):
            subseq = self.seq[i:i + k]
            res.append(
                round((subseq.count('C') + subseq.count('G') / len(subseq) * 100))
            )
        return res

    # Translate a DNA sequence into a amino_acid sequence
    def translate_seq(self, init_pos=0):
        return [DNA_Codons[self.seq[pos:pos + 3]] for pos in range(init_pos, len(self.seq) - 2, 3)]

    # Provides the frequency of each codon encoding a given amino_acid in a DNA sequence
    def codon_usage(self, aminoacid):
        tmpList = []
        for i in range(0, len(self.seq) - 2, 3):
            if DNA_Codons[self.seq[i:i + 3]] == aminoacid:
                tmpList.append(self.seq[i:i + 3])
        freqDict = dict(collections.Counter(tmpList))
        totalWight = sum(freqDict.values())
        for seq in freqDict:
            freqDict[seq] = round(freqDict[seq] / totalWight, 2)
            return freqDict