# HOW TO MAKE CONFIG FILE FOR TRIP/MpSA PROJECT

At first, *config file* contains **six** main blocks: *project*, *logging*, *core*, *input\_file*, *output\_dir* and *content*.

Second, the file must be saved with the ***json*** extension. The recommended name is ***config***.

## Section *project*:

You can define two values. 
**TRIP** - if you use *"Thousands of Reporters Integrated in Parallel"* or **MpSA** - if you are using *"Multiplex Sequence Analysis"*

#### *Example*

    {
        "project": "TRIP"
    }

## Section *logging*:

Contain single parameter - path to logging system config file. Default value is *"logging.json"*. This file is provided by default **AS IS** and can not be changed.

#### *Example*

    {
        "logging": "logging.json"
    }

## Section *core*:

Contain **five** dictionaries who defined main programm parameters.
### *arguments*
> contains required parameters

1. *experiment_name* - Some unique label for your experiment run
2. *paired* - This is **true** if you have forward**&**reverse reads after Illumina Mi/Hi-Seq or another sequence system, else it's **false**
3. *barcode_reversed* - If, when you were creating a library, you built in barcodes in reversed mode that set this option to **true**, else - **false**

#### *Example*

    {
        "core": {
            "arguments": {
                "experiment_name": "some_experiment_name",
                "paired": false,
                "barcode_reversed": true
            }
        }
    }
***
### *errors*
> contains some error reads values

1. *index_error* - error count in *index* sequence. The recommended value is **10%** of the length of the sequence.
2. *barcode_error* - error count in *barcode* sequence to filter out potential mutants. The recommended value is **10%** of the length of the sequence.
3. *sub\_index\_error* - ***Optional***, error count in *sub\_index* sequence. Use only if you have somewhere *sub\_index* sequence in your data. The recommended value is **10%** of the length of the sequence.

#### *Example*

    {
        "core": {
            "errors": {
                "index_error": 0,
                "sub_index_error": 1,
                "barcode_error": 2
            }
        }
    }
***
### *length*
> contains parameters of sequence lengths

1. *barcode_length* - the length of the barcode. The program looks for barcodes which are *barcode_length* ± 10%. 
2. *mutation_length* - ***Optional***, but ***required*** if you are using project **MpSA**. Determines the length of a mutation
3. *sub\_index\_length* - ***Optional***, but ***required*** if you specify *sub\_index\_error* value. Determines the length of a sub\_index

#### *Example*

    {
        "core": {
            "length": {
                "sub_index_length": 5,
                "barcode_length": 18
            }
        }
    }
***
### *thresholds*
> contains parameters that determine the restrictions applied to the output of the results

1. *minimum_reads* -  the minimum number of reads for considering a barcode genuine. An arbitrary recommendation is 5 reads. But if your sequencing depth is low than you might want to reduce it to 3 or 4.
2. *sequence_reliability* - the threshold below which the studied sequence (mutated region or region of the genome) is considered unreliable for a specific bar code

#### *Example*

    {
        "core": {
            "thresholds": {
                "minimum_reads": 3,
                "sequence_reliability": 0.8
            }
        }
    }
***
### *external*
> link to external libraries and additional files

1. *R* - specify the absolute path to you R installation
2. *bowtie* - **Required** for project *TRIP*, specify the absolute path to you *bowtie* installation
3. *bowtie_indexes* - **Required** for project *TRIP*, specify the absolute path to location of you *bowtie indexes*
4. *bowtie_filter_indexes* - ***Optional***, it is used only if you need to filter out a foreign genome, for example plasmids. Specify the absolute path to location of you *bowtie filter indexes*

#### *Example*

    {
        "core": {
            "external": {
                "R": "/usr/bin/Rscript",
                "bowtie": "/usr/bin/bowtie",
                "bowtie_indexes": "/path/to/your/indexes/folder/dmel-r6.19/dmel-r6.19",
                "bowtie_filter_indexes": "/path/to/your/indexes/folder/rfpl/rfpl"
            }
        }
    }

## Section *input_file*:

Contains a dictionary that indicates the location of the source data. The formats used are *.fastq*, .*fastq.gz*, *.picle*. If you have one file that contains the entire data set, use the *single* option. If there are many, then *multiple*.

### **ATTENTION!!!**
Use **ONLY ONE** option: *single* or *multiple*

### *single*
> contains a link to a file with a full set of data

* *path* - the absolute path to file

#### *Example*

    {
        "input_file": {
            "single": {
                "path": "/path/to/your/fastq/fastq.gz/or/pickle.file"
            }
        }
    }

**or you can use *multiple* option**

### *multiple*
> contains a links to a files with a separate set of data. If you have multiple replicates for each data set or for individual data, you can specify the path to each replicate, separated by a comma, according to the example.

* *fmap* - the absolute path to file/files that contains mapping data
* *fnorm* - the absolute path to file/files that contains normalization data
* *fexpr* - the absolute path to file/files that contains expression data

#### *Example*

    {
        "input_file": {
            "multiple": {
                "fmap": ["/path/to/your/fastq/fastq.gz/or/pickle.file", "/path/to/your/fastq/fastq.gz/or/pickle.file"],
                "fnorm": ["/path/to/your/fastq/fastq.gz/or/pickle.file", "/path/to/your/fastq/fastq.gz/or/pickle.file"],
                "fexpr": ["/path/to/your/fastq/fastq.gz/or/pickle.file", "/path/to/your/fastq/fastq.gz/or/pickle.file"]
            }
        }
    }

## Section *output_dir*:

Contain path to directory result

#### *Example*

    {
        "output_dir": "path/to/output/directory"
    }

## Section *content*:

Section in which you can determine the structure of your data. Data can consist of sets of mapping, normalization, expression or a combination thereof. The settings for each data set are determined by a dictionary that has a permanent structure.

### *mapping* or *normalization* or *expression*

* *replicates* - define replicates count for each data set who you use

#### *Example*

    {
        "content": {
            "mapping": {
                "replicates": 2
            }
        }
    }

* *index* - defined by an arbitrary index name and its sequence. It is important that the number of indices correspond to the number of replicas, even if the indices have the same sequences.

#### *Example*

    {
        "content": {
            "mapping": {
                "index": {
                    "m1": "TTCGGAGT",
                    "m2": "ACTCATTT" 
                }
            }
        }
    }

* *sub_index* - ***Optional***, but ***required*** if you specify *sub\_index\_error* value. Defines sequences of sub-indexes and their aliases.

#### *Example*

    {
        "content": {
            "mapping": {
                "sub_index": {
                    "PCNA": "AGTCA",
                    "Hsp70": "ACGTA",
                    "MtnA": "CTGCT",
                    "HexA": "AGCTC",
                    "Tbp": "TTGAG",
                    "Pyk": "TCAAA",
                    "Promoterless": "TCGCT"
                }
            }
        }
    }

* *constant_1* - *constant_**N*** - defines sequences of constant parts of a read. They can be any number instead ***"N"***, but it is better to adhere to a reasonable number of constant parts - no more than 4-5 sequences per read. Each constant part has an alias, sequence and an error level, which is recommended to be set as 10% of the length of the sequence

#### *Example*

    {
        "content": {
            "mapping": {
                "constant_1": ["gtcacaagggccggccacaactc", 3],
                "constant_2": ["ctcgatc", 1],
                "constant_3": ["ttaaccctagaaagataatc", 2],
            }
        }
    }

* *read_structure* - here we determine the order of the elements of which the reader consists.

To do this, use the following standard keywords and notations:
* **index** - index sequence
* **constant\_N** - constant part of read, where **"N"** is number of part
* **sub\_index** - ***optional*** defined sub-index in read
* **barcode** - barcode location in read
* **sequence** - the investigated sequence, *"mutation"* in **"MpSA"** project or *"genome"* in **"TRIP"**
* **N** - a fixed number of characters from dictionary - *"ATGC"*, not related to any of the above categories. 

#### *Example*

    {
        "content": {
            "mapping": {
                "read_structure": ["index", "constant_1", "barcode", 7, "sub_index", "constant_2", "sequence", "constant_3"]
            }
        }
    }
***
***

## And in the end, an example of a fully assembled config

    {
        "project": "TRIP",
        "logging": "logging.json",
        "core": {
            "arguments": {
                "experiment_name": "some_experiment_name",
                "paired": false,
                "barcode_reversed": true
            },
            "errors": {
                "index_error": 0,
                "sub_index_error": 1,
                "barcode_error": 2
            },
            "length": {
                "sub_index_length": 5,
                "barcode_length": 18,
            },
            "thresholds": {
                "minimum_reads": 3,
                "sequence_reliability": 0.8
            },
            "external": {
                "R": "/usr/bin/Rscript",
                "bowtie": "/usr/bin/bowtie",
                "bowtie_indexes": "/path/to/your/indexes/folder/dmel-r6.19/dmel-r6.19",
                "bowtie_filter_indexes": "/path/to/your/indexes/folder/rfpl/rfpl"
            }
        },
        "input_file": {
            "single": {
                "path": "path_to_fastq/fastq.gz/pickle_source"
                }
        },
        "output_dir": "path_to_output_directory",
        "content":{
            "mapping": {
                "replicates": 2,
                "index": {
                    "m1": "TTCGGAGT",
                    "m2": "ACTCATTT"
                },
                "sub_index": {
                    "AGTCA": "PCNA",
                    "ACGTA": "Hsp70",
                    "CTGCT": "MtnA",
                    "AGCTC": "HexA",
                    "TTGAG": "Tbp",
                    "TCAAA": "Pyk",
                    "TCGCT": "Promoterless"
                },
                "constant_1": ["gtcacaagggccggccacaactc", 3],
                "constant_2": ["ctcgatc", 1],
                "constant_3": ["ttaaccctagaaagataatc", 2],
                "read_structure": ["index", "constant_1", "barcode", 7, "sub_index", "constant_2", "sequence", "constant_3"]
            },
            "normalization": {
                "replicates": 2,
                "index": {
                    "n1": "AGTCGCCG",
                    "n2": "TAAACATC"
                },
                "sub_index": {
                    "PCNA": "AGTCA",
                    "Hsp70": "ACGTA",
                    "MtnA": "CTGCT",
                    "HexA": "AGCTC",
                    "Tbp": "TTGAG",
                    "Pyk": "TCAAA",
                    "Promoterless": "TCGCT"
                },
                "constant_1": ["gtcacaagggccggccacaactc", 3],
                "constant_2": ["ctcgatc", 1],
                "constant_3": ["ttaaccctagaaagataatc", 2],
                "read_structure": ["index", "constant_1", "barcode", 7, "sub_index", "constant_2"]
            },
            "expression": {
                "replicates": 2,
                "index": {
                    "e1": "GGTATGTT",
                    "e2": "GAGGGACC"
                },
                "sub_index": {
                    "AGTCA": "PCNA",
                    "ACGTA": "Hsp70",
                    "CTGCT": "MtnA",
                    "AGCTC": "HexA",
                    "TTGAG": "Tbp",
                    "TCAAA": "Pyk",
                    "TCGCT": "Promoterless"
                },
                "constant_1": ["gtcacaagggccggccacaactc", 3],
                "constant_2": ["ctcgatc", 1],
                "constant_3": ["ttaaccctagaaagataatc", 2],
                "read_structure": ["index", "constant_1", "barcode", 7, "sub_index", "constant_2"]
            }
        }
    }