# DNA Toolset/Code testing file
from bio_seq import *

test_dna = bio_seq()
test_dna.generate_rnd_seq(40,"DNA")
print(test_dna.get_seq_info())
print(test_dna.transcription())
print(test_dna.reverse_complement())
print(test_dna.translate_seq())
print(test_dna.get_seq_info())
print(test_dna.codon_usage('L'))

