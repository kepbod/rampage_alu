## A schematic flow shows the pipeline

## Prerequisites

### Softwares

* [fseq](http://fureylab.web.unc.edu/software/fseq/)

### Python libraries 

* [docopt](http://docopt.org/)
* [seqlib](https://github.com/kepbod/seqlib)
* [pysam](https://pysam.readthedocs.io/en/latest/index.html)
* [pybedtools](https://daler.github.io/pybedtools/)
* [numpy](http://www.numpy.org/)
* [joblib](https://joblib.readthedocs.io/en/latest/)

## Usage

### Step 1:

Fetch proper read pairs and remove PCR dunplicates.

```
Usage: rm_pcr.py [options] <rampage>...

Options:
    -h --help                      Show help message.
    --version                      Show version.
    -p THREAD --thread=THREAD      Threads. [default: 5]
    -o OUTPUT --output=OUTPUT      Output directory. [default: rampage_peak]
    --min=MIN                      Minimum read counts. [default: 1]
```

* Inputs: BAM files of RAMPAGE (`<rampage>...`)
* Output: A output folder containing relevant files (`-o OUTPUT`)
    * `rampage_plus_5end.bed`: BED file of the 5' end of plus strand read pairs
    * `rampage_plus_3read.bed`: BED file of the 3' end of plus strand read pairs
    * `rampage_minus_5end.bed`: BED file of the 5' end of minus strand read pairs
    * `rampage_minus_3read.bed`: BED file of the 3' end of minus strand read pairs
    * `rampage_link.bed`: BED file linking the 5' and 3' ends of read pairs

Example: 

```
rm_pcr.py -o rampage_peak rampage_rep1.bam rampage_rep2.bam
```

### Step 2:

Call peaks using 5' end of RAMPAGE read pairs

```
Usage: call_peak.py [options] <rampagedir>

Options:
    -h --help                      Show help message.
    -v --version                   Show version.
    -l LENGTH                      Feature length for F-seq. [default: 30]
    --wig                          Create Wig files.
    -p PERCENT                     Retained percent of reads in resized peaks.
                                   [default: 0.95]
```

* Input: the output folder created by `rm_pcr.py`
* Output: `rampage_peaks.txt` under the input folder

Format of `rampage_peaks.txt`:

| Field       | Description                   |
| :---------: | :---------------------------- |
| Chrom       | Chromosome                    |
| Start       | Start of peak region          |
| End         | End of peak region            |
| Name        | peak                          |
| Score       | 0                             |
| Strand      | Strand of peak                |
| Peak        | peak site                     |
| Height      | Height of peak site           |
| Peak reads  | Reads of (peak site ± 2 bp)   |
| Total       | Total reads of peak region    |
| Start_Fseq  | Start of F-seq peak region    |
| End_Fseq    | End of F-seq peak region      |
| RPM         | RPM of peak region            |

Example: 

```
call_peak.py rampage_peak
```

### Step 3:

Calculate entropy for RAMPAGE peaks

```
Usage: entropy.py [options] <rampagedir>

Options:
    -h --help                      Show help message.
    --version                      Show version.
    -p THREAD --thread=THREAD      Threads. [default: 5]
```

* Input: the output folder created by `rm_pcr.py`
* Output: `rampage_entropy.txt` under the input folder

Format of `rampage_entropy.txt`:

The first thirteen columns of `rampage_entropy.txt` are the same as `rampage_peaks.txt`.

The additional two columns are listed below
| Field       | Description                   |
| :---------: | :---------------------------- |
| Entropy     | Entropy of RAMPAGE peak       |
| 3' end      | 3' end of read pairs in peak  |

Example: 

```
entropy.py rampage_peak
```

### Step 4:

Annotate expressed Alu elements

```
Usage: annotate_alu.py [options] -f ref (-a alu | -r rep) <rampagedir>

Options:
    -h --help                      Show help message.
    --version                      Show version.
    -f ref                         Gene annotations.
    -t type                        File type of gene annotations.
                                   [default: ref]
    --promoter region              Promoter region. [default: 250]
    -a alu                         Alu annotations (BED format).
    -r rep                         Repeatmasker annotations (RMSK format).
    --extend length                Alu extended length. [default: 50]
    --entropy entropy              Entropy cutoff. [default: 2.5]
    --span span                    Span cutoff. [default: 1000]
    --coverage coverage            Coverage cutoff. [default: 0.5]
    -o out                         Output file. [default: alu_peak.txt]
```

* Input: 
    * gene annotation file (`-f ref`)
    * Alu annotation file (`-a alu` or `-r rep`)
    * the output folder created by `rm_pcr.py`
* Output: expressed Alu file (`-o out`) 

Format of expressed Alu file:

The first fifteen columns are the same as `rampage_entropy.txt`.

The additional six columns are listed below
| Field       | Description                   |
| :---------: | :---------------------------- |
| Chrom       | Chromosome of Alu             |
| Start       | Start of Alu                  |
| End         | End of Alu                    |
| Name        | Name of Alu                   |
| Score       | 0                             |
| Strand      | Strand of Alu                 |

Example: 

```
annotate_alu.py -f ref.txt -a alu.bed -o alu_peak.txt rampage_peak
```

## License
Copyright (C) 2018-2019 Xiao-Ou Zhang. See the [LICENSE](https://github.com/kepbod/rampage_alu/blob/master/LICENSE) file for license rights and limitations (MIT).