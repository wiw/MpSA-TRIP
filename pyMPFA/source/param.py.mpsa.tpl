# NB! index key template: <smth0>_<smth1>_<e|m|n><replicate_number>
indexList = {"Lib29_36_n1": "GTACCGTT", "Lib29_36_n2": "ATGCGTCA"}

indexError = 0
barcodeError = 2
barcodeLength = 18
mutationLength = 8
readsValue = 2
mutationProbability = 0.8
separator = " "
probability = 0.5

# Path to Rscript
rscript = '/usr/bin/Rscript'

# Path to paired_sequence_match.py
psmPy = '/usr/local/bin/paired_sequence_match.py'

const_1, const_2, const_3 = "cgccagggttttcccagtcacaagggccggccacaactc", "ctcgatcctcgagtgtcacctaaatcgtatgcggccgcgaattcttacttgtacagctcgtccatgccgagagtgatcccggcggc", ""
const_1Error, const_2Error, const_3Error = 3, 5, 0

# Normal regular expression variables
regExpIndex = {x:"^({}){{s<={}}}({}){{s<={}}}".format(x, str(indexError), const_1.upper(), str(const_1Error)) for x in indexList.values()}
regExpBcMut = "^(?P<barcode>.*)({}){{s<={}}}(?P<mutation>.*)({}){{s<={}}}.*$".format(const_2.upper(), str(const_2Error), const_3.upper(), str(const_3Error))
regExpBc = "^(?P<barcode>.*)({}){{s<={}}}".format(const_2.upper(), str(const_2Error))
