# NB! index key template: <smth0>_<smth1>_<e|m|n><replicate_number>

indexList = {"trip4_1_m1":"CAAGATAA", "trip4_1_m2":"GGACAACG"}
all_indexes = ['AGTCGCCG', 'TAAACATC', 'ACAATTCG', 'TACTTGTC', 'GGTATGTT', 'GAGGGACC', 'TAGCTCTA', 'TAATTGCG', 'CAAGATAA', 'GGACAACG', 'AGCGAGCT', 'CTGCACGT']

# promotor index
pmi = ['AGCTC', 'ACGTA', 'CTGCT', 'AGTCA', 'TCAA', 'TTGAG', 'TCGCT']
pmi_rev = ['GAGCT', 'TACGT', 'AGCAG', 'TGACT', 'TTGA', 'CTCAA', 'AGCGA']
pmiSubst = 0

indexError = 2
barcodeError = 2
barcodeLength = 18
pmiLength = 5
mutationLength = 8
readsValue = 0
mutationProbability = 0.80
separator = " "
probability = 0.5
# Max distance between paired-end reads into bowtie alignment
p_e_distance = 5000
# Expected minimal genome length got by the sequence parsing
expected_min_genome_len = 50

# Path to Rscript
rscript = '/usr/bin/Rscript'

# Path to paired_sequence_match.py
psmPy = '/usr/local/bin/paired_sequence_match.py'
bwIndex = '/home/anton/data/DAM/indexes/dmel-r6.19/dmel-r6.19'
rfplIndex = '/home/anton/data/DAM/indexes/rfpl/rfpl'

const_1, const_2, const_3, const_4 = "gtcacaagggccggccacaactc", "ctc", "gtacgtcacaatatgattatctttctagggttaa", "gatc"
const_1Error, const_2Error, const_3Error = 3, 1, 2

# Normal regular expression variables
regExpIndex = {x:"^({}){{s<={}}}({}){{s<={}}}".format(x, str(indexError), const_1.upper(), str(const_1Error)) for x in indexList.values()}

# expression for promotor index BcMut
regExpBcMut_fwd = "^(?P<barcode>.*)([ATGC]{{4}})(?P<pIndex>{}){{s<={}}}({}){{s<={}}}(?P<genome>{}.*)".format("|".join(pmi_rev), str(pmiSubst), const_2.upper(), str(const_2Error), const_4.upper())
regExpBcMut_rev = "^({}){{s<={}}}(?P<genome>.*)".format(const_3.upper(), str(const_3Error))

# expression for promotor index Barcode
regExpBc = "^(?P<barcode>.*)([ATGC]{{4}})(?P<pIndex>{}){{s<={}}}({}){{s<={}}}".format("|".join(pmi_rev), str(pmiSubst), const_2.upper(), str(const_2Error))

# expression for promotor index Filter
fltExp = ["TCTAGACCCTCCGGATTACTTGTACAGCTCGTCCATGCCGAGAGTGATCCCGGCGGCGGTCACGA", "AAGAGCTACCAACTCTTTTTCCGAAGGTAACTGGCTTCAGCAGAGCGCAGATACCAAATACTGTC", "CTTTTTTTCTGCGCGTAATCTGCTGCTTGCAAACAAAAAAACCACCGCTACCAGCGGTGGTTTGT", "TTCTTGAGATCCTTTTTTTCTGCGCGTAATCTGCTGCTTGCAAACAAAAAAACCACCGCTACCAG"]
fltExp_Error = 5

reFilter = ".*({})([ATGC]{{3}}GATC)({}){{e<={}}}.*".format("|".join(pmi_rev).upper(), "|".join(fltExp).upper(), str(fltExp_Error))

