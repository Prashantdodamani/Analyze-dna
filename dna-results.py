Python 3.12.3 (tags/v3.12.3:f6650f9, Apr  9 2024, 14:05:25) [MSC v.1938 64 bit (AMD64)] on win32
Type "help", "copyright", "credits" or "license()" for more information.
>>> from dnabox import *
>>> dna=dna_seq_generator(100)
>>> analyze_dna_sequence(dna)
{'GC_percent': 49.0, 'complement_dna': 'tagaggtggtagccatactggcgtgctagcatcacctagagttttatgccctcacaagggtcctgctactggtccaccgcctaaaagtatgggtaaatta', 'primer': 'tagaggtggtagccatactggcgt', 'primer_GC_percent': 54.166666666666664, 'primer_annealing_temperature': 74, 'mRNA': 'uagaggugguagccauacuggcgugcuagcaucaccuagaguuuuaugcccucacaaggguccugcuacugguccaccgccuaaaaguauggguaaauua', 'aminoacid-chain': 'X Arg Trp X Pro Tyr Trp X X Ser Thr Thr X Ser Phe Met Pro Ser Gln Gly Ser Cys Tyr Trp Ser Thr Ala X Lys Tyr Gly X Ile X'}
