import collections

from bio_structs import *
import random
class bio_seq:
    """DNA sequence class, Default value: ATCH, DNA, No label"""
    def __init__(self, seq="ATCG", seq_type="DNA", label="No_Label"):

        self.seq = seq.upper()
        self.label = label
        self.seq_type = seq_type
        self.is_valid = self.__validate()
        assert self.is_valid, f"Provided data does not seem to be a correct {self.seq_type}"


    def __validate(self):
        """check the sequence to make sure it is a valid DNA string"""
        return set(DNA_Nucleotides).issuperset(self.seq)


    def get_seq_biotype(self):
        """Returns sequence type"""
        return self.seq_type


    def get_seq_info(self):
        """Returns 4 strings , Full sequence information"""
        return f"[label]: {self.label}\n[Sequence]: {self.seq}\n[Biotype]: {self.seq_type}"


    def generate_rnd_seq(self, length=10, seq_type="DNA"):
        """Generate a random DNA sequence, provided the length"""
        seq = ''.join([random.choice(DNA_Nucleotides)
                       for x in  range(length)])
        self.__init__(seq, seq_type, "Randomly generated sequence")


    def countNucFrequency(self):
        return dict(collections.Counter(self.seq))


    def transcription(self):
        # DNA -> RNA Transcription
        rna = self.seq.replace("T", "U")
        return rna


    def reverse_complement(self):
        # return ''.join([DNA_ReverseComplement[nuc] for nuc in seq])[::-1]
        tmp_seq = self.seq
        mapping = str.maketrans('ATCG', 'TAGC')
        return tmp_seq.translate(mapping)[::1]


    def gc_content(self):
        return round((self.seq.count('C') + self.seq.count('G') / len(self.seq) * 100))


    def gc_content_subsec(self, k=20):
        """GC Content in a DNA/RNA sub_sequence length k, k=20 by default"""
        res = []
        for i in range(0, len(self.seq) - k + 1, k):
            subseq = self.seq[i:i + k]
            res.append(self.gc_content(subseq))
        return res


    def translate_seq(self, init_pos=0):
        """Translate a DNA sequence into an aminoacid sequence"""
        return [DNA_Codons[self.seq[pos:pos + 3]] for pos in range(init_pos, len(self.seq) - 2, 3)]
