# A1-4 mapping index
indexList = {"m18-1-1":"TTCGGAGT", "m18-1-2":"ACTCATTT", "m18-2-1":"GGGATCCG", "m18-2-2":"TCAAGCAA"}

# A10-13 normalization index
# indexList = {"n18-1-1":"AGTCGCCG", "n18-1-2":"TAAACATC", "n18-2-1":"ACAATTCG", "n18-2-2":"TACTTGTC"}

#A20-23 expression index
# indexList = {"e18-1-1":"GGTATGTT", "e18-1-2":"GAGGGACC", "e18-2-1":"TAGCTCTA", "e18-2-2":"TAATTGCG"}

# A5-6
# indexList = ["CAAGATAA", "GGACAACG"]

# A14-16
# indexList = {"norm":"GTACCGTT", "expr_rep1": "GCCACATA", "expr_rep2": "CCTATGGT"}

# promotor index
pmi = ['AGCTC', 'ACGTA', 'CTGCT', 'AGTCA', 'TCAAA', 'TTGAG', 'TCGCT']
pmi_rev = ['GAGCT', 'TACGT', 'AGCAG', 'TGACT', 'TTTGA', 'CTCAA', 'AGCGA']
pmiSubst = 1

indexError = 0
barcodeError = 2
barcodeLength = 18
# barcodeLength = 20
pmiLength = 5
mutationLength = 8
readsValue = 2
mutationProbability = 0.95
separator = " "
probability = 0.5

# Path to Rscript
rscript = '/usr/bin/Rscript'

# Path to paired_sequence_match.py
psmPy = '/usr/local/bin/paired_sequence_match.py'
bwIndex = '/home/anton/data/DAM/indexes/dmel-r6.19/dmel-r6.19'
rfplIndex = '/home/anton/data/DAM/indexes/rfpl/rfpl'

# BC21-mut21, BC15-mut15 - A1-4
# const_1, const_2, const_3 = "GACACTCGAGGATCGAG", "GAGTTGTGGCCGGCCCTTGTGACTGGGAAAACCCTGGCGTAAATAAAATACGAAATGACTAGTcatgcgtc", "gcatgattatctttaacgtacgtcACAAT"
# const_1Error, const_2Error, const_3Error = 2, 5, 3

# Lib deltaC, wt - A14-16
# const_1, const_2, const_3 = "CGCCAGGGTTTTCCCAGTCACAAGGGCCGGCCACAACTC", "CTCGATCCTCGAGTGTCACCTAAATCGTATGCGGCCGCGAATTCTTACTTGTACAGCTCGTCCATGCCGAGAGTGATCCCGGCG", ""
# const_1Error, const_2Error, const_3Error = 3, 5, 0

# Lib cDNA 29-36 A5-6
# const_1, const_2, const_3 = "CGCCAGGGTTTTCCCAGTCACAAGGGCCGGCCACAACTC", "CTCGATCCTCGAGTGTCACCTAAATCGTATGCGGCCGCGAATTCTTACTTGTACAGCTCGTCCATGCCGAGAGTGATCCCGGCGGC", ""
# const_1Error, const_2Error, const_3Error = 3, 5, 0

# Lib cDNA 29-36 A14-16
# const_1, const_2, const_3 = "cgccagggttttcccagtcacaagggccggccacaactc", "ctcgatcctcgagtgtcacctaaatcgtatgcggccgcgaattcttacttgtacagctcgtccatgccgagagtgatcccggcggc", ""
# const_1Error, const_2Error, const_3Error = 3, 5, 0

# Genome Lib Normalization&Expression A10-13, A20-23
# const_1, const_2, const_3 = "cgccagggttttcccagtcacaagggccggccacaactc", "ctcgatctctagaccctccggattacttgtacagctcgtccatgccgagagtgatcccggcggcggtcacgaactcc", ""
# const_1Error, const_2Error, const_3Error = 3, 5, 0

# Genome Lib Mapping A1-4
const_1, const_2, const_3 = "gtcacaagggccggccacaactc", "ctcgatc", "ttaaccctagaaagataatc"
const_1Error, const_2Error, const_3Error = 3, 1, 2

# Normal regular expression variables
regExpIndex = {x:"^({}){{s<={}}}({}){{s<={}}}".format(x, str(indexError), const_1.upper(), str(const_1Error)) for x in indexList.values()}

# casual expression
# regExpBcMut = "^(?P<barcode>.*)({}){{s<={}}}(?P<mutation>.*)({}){{s<={}}}.*$".format(const_2.upper(), str(const_2Error), const_3.upper(), str(const_3Error))

# expression for promotor index BcMut
# regExpBcMut = "^(?P<barcode>.*)([ATGC]{{4}})(?P<pIndex>{}){{s<={}}}({}){{s<={}}}(?P<genome>.*)({}|$){{s<={}}}".format("|".join(pmi_rev), str(pmiSubst), const_2.upper(), str(const_2Error), const_3.upper(), str(const_3Error))
regExpBcMut = "^(?P<barcode>.*)([ATGC]{{4}})(?P<pIndex>{}){{s<={}}}({}){{s<={}}}(((?P<genome1>.*)({}){{s<={}}}.*)|((?P<genome2>.*)$))".format("|".join(pmi_rev), str(pmiSubst), const_2.upper(), str(const_2Error), const_3.upper(), str(const_3Error))

# casual expression
# regExpBc = "^(?P<barcode>.*)({}){{s<={}}}".format(const_2.upper(), str(const_2Error))

# expression for promotor index Barcode
regExpBc = "^(?P<barcode>.*)([ATGC]{{4}})(?P<pIndex>{}){{s<={}}}({}){{s<={}}}".format("|".join(pmi_rev), str(pmiSubst), const_2.upper(), str(const_2Error))

# expression for promotor index Filter
fltExp = ["TCTAGACCCTCCGGATTACTTGTACAGCTCGTCCATGCCGAGAGTGATCCCGGCGGCGGTCACGA", "AAGAGCTACCAACTCTTTTTCCGAAGGTAACTGGCTTCAGCAGAGCGCAGATACCAAATACTGTC", "CTTTTTTTCTGCGCGTAATCTGCTGCTTGCAAACAAAAAAACCACCGCTACCAGCGGTGGTTTGT", "TTCTTGAGATCCTTTTTTTCTGCGCGTAATCTGCTGCTTGCAAACAAAAAAACCACCGCTACCAG"]
fltExp_Error = 5
# eGFP = "TCTAGACCCTCCGGATTACTTGTACAGCTCGTCCATGCCGAGAGTGATCCCGGCGGCGGTCACGA"
# dpn1 = "AAGAGCTACCAACTCTTTTTCCGAAGGTAACTGGCTTCAGCAGAGCGCAGATACCAAATACTGTC"
# dpn2 = "CTTTTTTTCTGCGCGTAATCTGCTGCTTGCAAACAAAAAAACCACCGCTACCAGCGGTGGTTTGT"
# dpn3 = "TTCTTGAGATCCTTTTTTTCTGCGCGTAATCTGCTGCTTGCAAACAAAAAAACCACCGCTACCAG"

# re_eGFP = ".*({}).*{}.*".format("|".join(pmi_rev), eGFP.upper())
reFilter = ".*({})([ATGC]{{3}}GATC)({}){{e<={}}}.*".format("|".join(pmi_rev).upper(), "|".join(fltExp).upper(), str(fltExp_Error))


###########################################################
# args = p.parse_args(["--input", "/home/anton/backup/input/trip/RUN_2017-11-27/sample_S1_L001_R1_001.fastq.gz", "--mode", "genome", "--output", "/home/anton/backup/input/trip/RUN_2017-11-27/results", "--experiment_label", "Genome_Norm_A10-13", "--make_barcode_library", "--reverse_barcode"])