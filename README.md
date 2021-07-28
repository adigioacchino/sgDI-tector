# sgDI-tector: a tool for sgRNA discovery and quantification from NGS data

sgDI-tector is a python 3 script for viral subgenomic RNA (sgRNA) detection
from next-generation sequencing (NGS) data. 
It can be used also to extract informations about
the leader sequence and transcriptional regulatory sequence (TRS) from NGS data.
sgDI-tector works by first running [DI-tector](https://rnajournal.cshlp.org/content/24/10/1285) 
and then analyzing its output files. 
Notice that DI-tector can be run independently, and sgDI-tector will
not run it again if its standard output files are found in the
output folder (unless the option --Force_DItector is provided).

### Input files
Minimum input files are the virus sequence (in fasta format)
and the NGS data (in fastq.gz format). Notice that the fastq.gz file is
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
- `TRS_logo.pdf`: the logo of the 20 nt preceding and following each RI
position is given, after an alignment step.
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
`NC_045512.2.fasta`. Moreover, an example of file that can be given as
input with the -r or --sgRNA_Ref option is `NC_045512.2_sgRNA.fasta`: it
is obtained from the same NCBI record, and it contains only the ORF 
known to be expressed thorugh sgRNAs. sgRNA families found by sgDI-tector
which can be aligned well to an entry of this file will have as name the
description of the entry (that is, the string after the symbol '>').

The folder "utilities" contains the package `RNA_sequence_logo.py`, that
is needed to print the TRS sequence logo.


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
- [bwa](http://www.htslib.org/),
- [samtools](http://www.htslib.org/) and
- [mafft](https://mafft.cbrc.jp/alignment/software/).  

Notice that bwa and samtools are only needed to run DI-tector, so if
DI-tector output data are already available there is no need for these
softwares.

Some python3 libraries are needed to run this script, in particular
`numpy`, `pandas` and `biopython`. Moreover, the package `RNA_sequence_logo.py`
is needed to print the TRS logo. This last package must be in the running
folder of sgDI-tector, or in a subfolder called "utilities".

## Acknowledgements:
The utility python package used for printing RNA logos (`RNA_sequence_logo.py`) 
has been written by Jerome Tubiana, see the licence note in the file.

## Contacts:
For comments or question, feel free to [contact me](mailto:andrea.dgioacchino@gmail.com).
