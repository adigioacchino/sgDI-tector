# sgDI-tector: a tool for sgRNA discovery and quantification from NGS data

sgDI-tector is a python 3 script for viral subgenomic RNA (sgRNA) detection
from next-generation sequencing (NGS) data.
It can be used also to extract information about
the leader sequence and transcriptional regulatory sequence (TRS) from NGS data.
sgDI-tector works by first running [DI-tector](https://rnajournal.cshlp.org/content/24/10/1285)
and then analyzing its output files.
Notice that DI-tector can be run independently, and sgDI-tector will
not run it again if its standard output files are found in the
output folder (unless the option --Force_DItector is provided).
More details and examples of usage are given in the [paper describing the
package](https://rnajournal.cshlp.org/content/28/3/277).

### Input files
Minimum input files are the virus sequence (in fasta format)
and the NGS data (in fastq.gz or fastq format). Notice that the fastq.gz file is
used only when DI-tector is run. Optionally, a file with reference ORFs
translated by sgRNA expression can be given as input
(through the -r or --sgRNA_Ref option), in fasta (CDS) format.
If this file is provided, sgDI-tector attempts aligning
all ORFs corresponding to the sgRNAs found against the reference, and the
output file will be annotated accordingly.
Examples of input files can be found within this repository in the folder
"input_examples".

### Output files
sgDI-tector outputs the following files:
- DI-tector output files (if DI-tector is run), see
[DI-tector documentation](https://rnajournal.cshlp.org/content/24/10/1285);
- `leader_sequence.txt`: the leader sequence (from the 5' end of viral genome)
is given.
- `sgRNA_junctions.csv`: for each family of sgRNA, the BP (breakpoint)
and RI (reinitiation) sites are given, together with the number of counts
of the sgRNA and the position of the first ATG codon in the sgRNA produced.
Moreovoer, a name identifying the ORF is given, and this name is obtained
after comparison with the reference subgenomic ORFs if the
option -r or --sgRNA_Ref is used.
- `detected_TRS_motifs.csv`: for each family of sgRNA, the longest subsequence
in the 20 nts around the RI site having an identical overlapping subsequence
in the last part of the leader sequence is given, together with some info
about the sgRNA family it comes from. Notice that subsequences shortest
than a threshold (fixed by the algorithm so that less than 5% of all possible
20-nt subsequences would give a hit) are not considered.
The output folder can be chosen by using the -o or --Output_Directory
option, if no folder is provided the files will be created in the working
directory.

### Options
The full list of option can be obtained by running
`python3 sgDI-tector.py -h`. Notice that all the options whose description
starts with "[DI-tector]" only influence DI-tector behaviour, so they are
useless if DI-tector is not run by sgDI-tector.

## Repository structure:
The sgDI-tector script file is `sgDI-tector.py`.

The folder "input_examples"
contains a SARS-CoV-2 reference sequence (obtained from NCBI) fasta file,
`NC_045512.2.fasta`, and a small example of input fastq data named `SARS2_subsampling.fastq.gz`.
Moreover, this folder contains an example of file that can be given as
input with the -r or --sgRNA_Ref option is `NC_045512.2_sgRNA.fasta`: it
is obtained from the same NCBI record, and it contains only the ORF
known to be expressed thorugh sgRNAs. sgRNA families found by sgDI-tector
which can be aligned well to an entry of this file will have as name the
description of the entry (that is, the string after the symbol '>').

## Dependencies:
sgDI-tector can be run on Linux or Windows. In the second case, all the
bioinformatic tools (and DI-tector) are run through [WSL](https://ubuntu.com/wsl),
which must be installed. All tools discussed below must be installed
within WSL, and in PATH.

sgDI-tector needs several other tools/libraries to be run. First, it
needs DI-tector (version 0.6) if the output files of DI-tector are not
already available. The DI-tector script, must be in the running
directory of sgDI-tector and named `DI-tector_06.py`.
DI-tector can be obtained from [this website](http://www.di-tector.cyame.eu/).

The other tools needed are common bioinformatic softwares:
- [bwa](http://www.htslib.org/) and
- [samtools](http://www.htslib.org/).

Notice that bwa and samtools are only needed to run DI-tector, so if
DI-tector output data are already available there is no need for these
softwares.

Some python3 libraries are needed to run this script, in particular
`numpy`, `pandas` and `biopython`.

## Example:
The command to use sgDI-tector tool on the provided example files is
`python3 sgDI-tector.py -r ./input_examples/NC_045512.2_sgRNA.fasta ./input_examples/NC_045512.2.fasta ./input_examples/SARS2_subsampling.fastq.gz`.
Notice that for this to work, the DI-tector script named `DI-tector_06.py` must be in the running directory of sgDI-tector (see section "Dependencies" above).

## How to cite this software:
If you  use this software for an academic publication, please cite the
[paper](https://rnajournal.cshlp.org/content/28/3/277) that describes it:

> Di Gioacchino, Andrea, et al. "sgDI-tector: defective interfering viral
  genome bioinformatics for detection of coronavirus subgenomic RNAs."
  RNA 28.3 (2022): 277-289.

## Contacts:
For comments or question, feel free to [contact me](mailto:andrea.dgioacchino@gmail.com).
