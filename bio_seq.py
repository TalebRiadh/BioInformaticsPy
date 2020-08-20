import collections

from bio_structs import DNA_Codons, RNA_Codons, NUCLEOTIDE_BASE
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
        return set(NUCLEOTIDE_BASE[self.seq_type]).issuperset(self.seq)

    # Return sequence type
    def get_seq_biotype(self):
        return self.seq_type

    # Returns 4 strings, full sequence information
    def get_seq_info(self):
        return f"[label]: {self.label}\n[Sequence]: {self.seq}\n[Biotype]: {self.seq_type}"

    # Generate a random DNA sequence, provided the length
    def generate_rnd_seq(self, length=10, seq_type="DNA"):
        seq = ''.join([random.choice(NUCLEOTIDE_BASE[seq_type])
                       for x in range(length)])
        self.__init__(seq, seq_type, "Randomly generated sequence")

    # Count nucleotides in a given sequence, return a dictionary
    def nucleotide_frequency(self):
        return dict(collections.Counter(self.seq))

    # DNA -> RNA Transcription
    # return a copier of the sequence
    def transcription(self):
        if self.seq_type == "DNA":
            return self.seq.replace("T", "U")
        return "Not a DNA sequence"
    # return ''.join([DNA_ReverseComplement[nuc] for nuc in seq])[::-1]
    def reverse_complement(self):
        if self.seq_type == "DNA":
            mapping = str.maketrans("ATCG", "TAGC")
        else:
            mapping = str.maketrans("AUCG", "UAGC")
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
        if self.seq_type == "DNA":
            return [DNA_Codons[self.seq[pos:pos + 3]] for pos in range(init_pos, len(self.seq) - 2, 3)]
        elif self.seq_type == "RNA":
            return [RNA_Codons[self.seq[pos:pos + 3]] for pos in range(init_pos, len(self.seq) - 2, 3)]

    # Provides the frequency of each codon encoding a given amino_acid in a DNA sequence
    def codon_usage(self, aminoacid):
        tmpList = []
        if self.seq_type == "DNA":
            for i in range(0, len(self.seq) - 2, 3):
                if DNA_Codons[self.seq[i:i + 3]] == aminoacid:
                    tmpList.append(self.seq[i:i + 3])
        elif self.seq_type == "RNA":
            for i in range(0, len(self.seq) - 2, 3):
                if RNA_Codons[self.seq[i:i + 3]] == aminoacid:
                    tmpList.append(self.seq[i:i + 3])
        freqDict = dict(collections.Counter(tmpList))
        totalWight = sum(freqDict.values())
        for seq in freqDict:
            freqDict[seq] = round(freqDict[seq] / totalWight, 2)
            return freqDict

    # Generate the six reading frames of a DNA sequence, including reverse the reverse complement
    def gen_reading_frames(self):
        frames = []
        frames.append(self.translate_seq(0))
        frames.append(self.translate_seq(1))
        frames.append(self.translate_seq(2))
        tmp_seq = bio_seq(self.reverse_complement(), self.seq_type)
        frames.append(tmp_seq.translate_seq(0))
        frames.append(tmp_seq.translate_seq(1))
        frames.append(tmp_seq.translate_seq(2))
        del tmp_seq
        return frames

    # Compute al possible proteins on a amino_acid sequence and return a list of possible proteins
    def proteins_from_rf(self, aa_seq):
        current_prot = []
        proteins = []
        for aa in aa_seq:
            if aa == "_":
                # STOP accumulating amino acids if _
                if current_prot:
                    for p in current_prot:
                        proteins.append(p)
                    current_prot = []
            else:
                # START accumulating amino acids if M
                if aa == "M":
                    current_prot.append("")
                for i in range(len(current_prot)):
                    current_prot[i] += aa
        return proteins

    # Compute all possible proteins for all open reading frames
    def all_proteins_from_orfs(self, startReadPos=0, endReadPos=0, ordered=False):
        if endReadPos > startReadPos:
            tmp_seq = bio_seq(self.seq[startReadPos: endReadPos], self.seq_type)
            rfs = tmp_seq.gen_reading_frames()
        else:
            rfs = self.gen_reading_frames()
        res = []
        for rf in rfs:
            prots = self.proteins_from_rf(rf)
            for p in prots:
                res.append(p)
        if ordered:
            return sorted(res, key=len, reverse=True)
        return res