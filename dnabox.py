def gc_percent(dna):
    gc_count=dna.count('g')+dna.count('c')
    gc_percent=(gc_count/len(dna))*100.0
    return gc_percent

def stop_codon(dna):
    stop_codon_found=False
    stop_codons=['tag','taa','tga']
    for i in range(0,len(dna),3):
        codon=dna[i:1+3].lower()
        if codon in stop_codons:
            stop_codon_fount=True
            break
        return stop_codon_found

def complement(dna):
    basecomplement={'A':'T','T':'A','C':'G','G':'C','a':'t','t':'a','c':'g','g':'c'}
    letters=list(dna)
    letters=[basecomplement[base] for base in letters]
    return''.join(letters)

def primer_generator(dna):
    "this function generates a complementary primer of length 24 bases for a given dna sequence"
    complement={'a':'t','A':'T','t':'a','T':'A','c':'g','C':'G','g':'c','G':'C'}
    primer=""
    for base in dna[:24]:
        complement_base=complement[base]
        primer+=complement_base
            
    return primer

    
import random
def dna_seq_generator(length):
    alphabet= "atcg"
    dna=''.join(random.choice(alphabet) for i in range (length))
    return dna

def annealing_temperature(primer):
    n_A = primer.lower().count('a')
    n_T = primer.lower().count('t')
    n_G = primer.lower().count('g')
    n_C = primer.lower().count('c')
    annealing_temperature = (2 * (n_A + n_T) + 4 * (n_G + n_C))
    return annealing_temperature


from Bio.Blast import NCBIWWW, NCBIXML
def find_species(sequence):
    result_handle = NCBIWWW.qblast("blastn", "nt", sequence)
    blast_record = NCBIXML.read(result_handle)
    species_list = [alignment.title.split('|')[-1] for alignment in blast_record.alignments]
    result_handle.close()
    return species_list


def find_orfs(sequence):
    orfs = []
    for i in range(len(sequence) - 2):
        if sequence[i:i+3] == "ATG":
            for j in range(i+3, len(sequence), 3):
                codon = sequence[j:j+3]
                if codon in ["TAA", "TAG", "TGA"]:
                    orfs.append(sequence[i:j+3])
                    break
    return orfs
def longest_orf_length(sequence):
    orfs = find_orfs(sequence)
    if not orfs:
        return 0
    
    return max(len(orf) for orf in orfs)
    
def aminoacid_chain(mRNA):
    aminoacidcomplement = {
        'ggc': 'Gly', 'gcc': 'Ala', 'gca': 'Ala', 'gcg': 'Ala',
        'cgc': 'Arg', 'cga': 'Arg', 'cgg': 'Arg', 'aga': 'Arg', 'agg': 'Arg',
        'aau': 'Asn', 'aac': 'Asn',
        'gau': 'Asp', 'gac': 'Asp',
        'ugu': 'Cys', 'ugc': 'Cys',
        'caa': 'Gln', 'cag': 'Gln',
        'gaa': 'Glu', 'gag': 'Glu',
        'ggu': 'Gly', 'gga': 'Gly', 'ggg': 'Gly',
        'cau': 'His', 'cac': 'His',
        'auu': 'Ile', 'auc': 'Ile', 'aua': 'Ile',
        'uua': 'Leu', 'uug': 'Leu', 'cuu': 'Leu', 'cuc': 'Leu', 'cua': 'Leu', 'cug': 'Leu',
        'aaa': 'Lys', 'aag': 'Lys',
        'aug': 'Met',
        'uuu': 'Phe', 'uuc': 'Phe',
        'ccu': 'Pro', 'ccc': 'Pro', 'cca': 'Pro', 'ccg': 'Pro',
        'ucu': 'Ser', 'ucc': 'Ser', 'uca': 'Ser', 'ucg': 'Ser', 'agu': 'Ser', 'agc': 'Ser',
        'auc': 'Thr', 'acc': 'Thr', 'aca': 'Thr', 'acg': 'Thr',
        'ugg': 'Trp',
        'uau': 'Tyr', 'uac': 'Tyr',
        'guu': 'Val', 'guc': 'Val', 'gua': 'Val', 'gug': 'Val'
    }

    codons = [mRNA[i:i + 3] for i in range(0, len(mRNA), 3)]
    aminoacids = [aminoacidcomplement.get(codon.lower(), 'X') for codon in codons]
    return ' '.join(aminoacids)


def mRNA_sequence(dna):
    basecomplement = {
        'a': 'u', 't': 'a', 'c': 'g', 'g': 'c'
    }
    letters = list(dna)
    letters = [basecomplement[base] for base in letters]
    return ''.join(letters)


from dnabox import gc_percent, complement, primer_generator, annealing_temperature

def analyze_dna_sequence(dna):
    result = {}
    result['GC_percent'] = gc_percent(dna)
    result['complement_dna'] = complement(dna)
    primer = primer_generator(dna)
    result['primer'] = primer
    result['primer_GC_percent'] = gc_percent(primer)
    primer_annealing_temperature = annealing_temperature(primer)
    result['primer_annealing_temperature'] = primer_annealing_temperature
    mRNA = mRNA_sequence(dna)
    result['mRNA']= mRNA
    result['aminoacid-chain'] = aminoacid_chain(mRNA)

    return result




 








