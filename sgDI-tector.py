"""
 Copyright 2021 - by Andrea Di Gioacchino (andrea.dgioacchino@gmail.com)
     All rights reserved
     
     Permission is granted for anyone to copy, use, or modify this
     software for any uncommercial purposes, provided this copyright 
     notice is retained, and note is made of any changes that have 
     been made. This software is distributed without any warranty, 
     express or implied. In no event shall the author or contributors be 
     liable for any damage arising out of the use of this software.
     
     The publication of research using this software, modified or not, must include 
     appropriate citations to: 
"""

#=================================

import sys, os, subprocess, argparse, shlex
from shutil import which
from pathlib import Path

import bisect
import pandas as pd
import numpy as np

from Bio import SeqIO
from Bio.Seq import Seq
from Bio import pairwise2

sys.path.append('./utilities/')
import RNA_sequence_logo


#=================================

def align_junction_LeftRight(viral_seq, bp_pos, ri_pos, align_to="L"):
    """If align_to="L", put all the ambiguous nucleotides in the
    'body' part of the junction defined by bp_pos and ri_pos,
    that is the junction point is moved as close to the 5' 
    end (of viral_seq) as possible. If align_to="R", do the opposite."""    
    py_ri = ri_pos-1 # in this way viral_seq[py_ri] is the first nucleotide after the junction obtained by DI-tector
    py_bp = bp_pos-1 # in this way viral_seq[py_bp] is the last nucleotide before the junction obtained by DI-tector    
    assert (align_to == "L" or align_to == "R"), "Plese enter R or L to align as right as possible or as left as possible."
    new_bp_pos = py_bp
    new_ri_pos = py_ri
    try_next_alignement = True
    while try_next_alignement:
        if align_to == "L":
            if vir_seq[new_bp_pos] == vir_seq[new_ri_pos-1]:
                new_bp_pos -= 1
                new_ri_pos -= 1
            else:
                try_next_alignement = False
        elif align_to == "R":
            if vir_seq[new_bp_pos+1] == vir_seq[new_ri_pos]:
                new_bp_pos += 1
                new_ri_pos += 1
            else:
                try_next_alignement = False     
    new_bp_pos += 1 # in this way I am using a fixed convention
    new_ri_pos += 1 # in this way I am using a fixed convention
    return new_bp_pos, new_ri_pos


def findall(p, s):
    """Yields all the positions of
    the pattern p in the string s."""
    i = s.find(p)
    while i != -1:
        yield i
        i = s.find(p, i+1)


def find_top_BP_intervals(df_bpcount, interval_size=20, plot_flag=False,
                         key_bp="BP_Pos", key_count="count"):
    """
    Return the starting position of the interval of length interval_size of the reference sequence
    such that the maximum number of BP sites are in the interval. When plot_flag=True, returns all
    the number of BP sites in each interval (sliding windows).
    """
    max_start_pos = max(df_bpcount[key_bp])
    max_nbp = 0
    max_nbp_pos = [1]
    if plot_flag:
        plot_nbp = []
    for pos in range(1, max_start_pos):
        p_start = pos
        p_end = pos + interval_size
        t_nbp = df_bpcount[(df_bpcount[key_bp] > p_start) & (df_bpcount[key_bp] < p_end)][key_count].sum()        
        if plot_flag:
            plot_nbp.append(t_nbp)
        if t_nbp > max_nbp:
            max_nbp_pos = [pos]
            max_nbp = t_nbp
        elif t_nbp == max_nbp:
            max_nbp_pos.append(pos)
    if plot_flag:
        return (max_nbp_pos, plot_nbp)
    else:
        return max_nbp_pos


def cluster_subgenomic_families(df_start, vir_seq, key_bp="BP_pos", key_ri="RI_pos", key_count="count"):
    """
    Cluster together families which have the same first ATG after the RI_Pos. 
    For each cluster/family, the counts are summed together and the BP_Pos and RI_Pos shown
    are those of the member with the largest number of counts before the clustering.
    """
    # first step: find all ATG positions (the +1 is to have starting pos as 1 instead of 0)
    atg_startpos = [i+1 for i in findall('ATG', vir_seq)]
    # for each entry in df_start, find the first atg after that and add this info as a column to the df
    ri_poss = df_start[key_ri]
    first_atg_df = []
    for ri in ri_poss:
        if ri > atg_startpos[-1]:
            first_atg_df.append(-1)
        else:
            first_atg_df.append(atg_startpos[bisect.bisect_left(atg_startpos, ri)])
    t_df = df_start.copy()
    t_df["atg_pos"] = first_atg_df
    # create final dataframe, by keeping for each atg_pos only the one with largest number of counts
    fam_bp_poss = []
    fam_ri_poss = []
    fam_counts = []
    unique_atg_pos = list(set(first_atg_df))
    for atgp in unique_atg_pos:
        t_df2 = t_df[t_df["atg_pos"] == atgp].reset_index(drop=True)
        t_df2 = t_df2.sort_values(by=key_count, ascending=False)
        fam_counts.append(t_df2[key_count].sum())
        fam_bp_poss.append(t_df2[key_bp].to_list()[0])
        fam_ri_poss.append(t_df2[key_ri].to_list()[0])
    family_df = pd.DataFrame({'bp_pos' : fam_bp_poss, 'ri_pos' : fam_ri_poss,
                              'counts' : fam_counts, 'atg_pos' : unique_atg_pos})
    return family_df.sort_values(by='counts', ascending=False).reset_index(drop=True)


def seq2num(seq):
    """Transform a RNA sequence or list of RNA sequences
    (seq) into a list of numbers."""
    aa = ['A', 'C', 'G', 'U', '-']
    aadict = {aa[k]:k for k in range(len(aa))}
    if type(seq) == str:
        return np.array([aadict[x] for x in seq], dtype=np.int16)[np.newaxis,:]
    elif type(seq) ==list:
        return np.array([[aadict[x] for x in seq_] for seq_ in seq], dtype=np.int16)

    
def average(MSA2num, q):
    """Compute the average value in each column of the MSA
    MSA2num provided, q being the size of the alphabet."""
    B = MSA2num.shape[0]
    N = MSA2num.shape[1]
    out = np.zeros((N, q), dtype=np.float32)
    for i in range(q):
        out.T[i] = np.sum(MSA2num == i, axis=0)
    out /= B
    return out

#=================================

if __name__ =='__main__':
    # parser for command-line options
    parser = argparse.ArgumentParser()
    parser.add_argument("Virus_Ref", help="Virus genome reference sequence in FASTA format.")
    parser.add_argument("Input_Data", help="File containing single reads in FASTQ format.")
    parser.add_argument("-g", "--Host_Ref", help="[DI-tector] Host genome reference sequence in FASTA format (optional).")
    parser.add_argument("-r", "--sgRNA_Ref", help="sgRNAs of reference viral sequence, in FASTA format (optional).")

    parser.add_argument("-s", "--Min_Segment", help="[DI-tector] Minimum segment length. Default is 15.", type=int, default=15)
    parser.add_argument("-m", "--Min_MAPQ", help="[DI-tector] Skip alignments with MAPQ smaller than INT. Default is 25.", type=int, default=25)
    parser.add_argument("-n", "--Nb_Reads", help="[DI-tector] Show only DVGs with counts reads greater or equal to INT. Default is 2.))", type=int, default=2)
    parser.add_argument("-o", "--Output_Directory", help="[DI-tector] Directory that will contains all output files.")
    parser.add_argument("-t", "--Tag", help="[DI-tector] Tag name that will be appended before each DI-tector output file. Default is 'DI-tector'.", default="DI-tector")
    parser.add_argument("-d", "--DVG_sequences", action='store_true', help="[DI-tector] Generate multi-fasta file with DVG sequences. Default is (OFF).")
    parser.add_argument("-l", "--InDel_Length", help="[DI-tector] Skip alignments with size of InDels smaller or equal to INT. Default size is 1.", type=int, default=1)
    parser.add_argument("-f", "--Fasta", action='store_true', help="[DI-tector] Use '-f' if data is in FASTA format fasta. Default is FASTQ.")
    parser.add_argument("-p", "--Polarity", help="[DI-tector] [0] Positive strand genome / [1] Negative strand genome. Default is 0.", type=int, default=0)
    parser.add_argument("-q", "--No_Quantification", action='store_true', help="[DI-tector] Inactive percentage quantification. Needs bedtools. Default is (ON).")
    parser.add_argument("-k", "--Keep_files", action="store_true", help="[DI-tector] Keep intermediate files (i.e. alignment etc...). Default is (OFF).")
    parser.add_argument("-x", "--Nb_threads", help="Number of threads. Default is 1.", type=int, default=1)
    parser.add_argument("--Force_DItector", action="store_true", help="Run DI-tector even if DI-tector output files are found in the output directory. The DI-tector output files will be overwritten. Default is (OFF).")
    args = parser.parse_args()
    
    # windows accepts both '/' and '\' for specifying paths, while linux accepts only '/'.
    # so let's convert possible windows paths in linux ones
    if sys.platform == 'win32':
        args.Virus_Ref = os.path.normpath(args.Virus_Ref).replace("\\",'/')
        args.Input_Data = os.path.normpath(args.Input_Data).replace("\\",'/')
        if args.Host_Ref:
            args.Host_Ref = os.path.normpath(args.Host_Ref).replace("\\",'/')
        if args.sgRNA_Ref:
            args.sgRNA_Ref = os.path.normpath(args.sgRNA_Ref).replace("\\",'/')
        if args.Output_Directory:
            args.Output_Directory = os.path.normpath(args.Output_Directory).replace("\\",'/')
    
    # prepare options for calling DI-tector
    ditect_argv = sys.argv.copy()[1:]
    try:
        t_index = ditect_argv.index("-r")
        del ditect_argv[t_index: t_index + 2]
    except ValueError:
        pass
    try:
        t_index = ditect_argv.index("--sgRNA_Ref")
        del ditect_argv[t_index: t_index + 2]
    except ValueError:
        pass
    try:
        t_index = ditect_argv.index("--Force_DItector")
        del ditect_argv[t_index]
    except ValueError:
        pass
    ditect_options = " ".join(ditect_argv).replace("\\","/")
    
    # sgDI-tector starts!
    # sgDI-tector parameters
    leader_window = args.Min_Segment + 5 # that is, default: 20.
    threshold_sgRNAnames = 0.9 # threshold for assigning sgRNA name
    junct_window = 20 # junction window to consider for alignment and logo production
    gap_thresh = 0.5 # fraction of rows containing gaps necessary to exclude a column in the logo
    
    print("###########################################################")
    print("This is sgDI-tector.")
    if sys.platform == 'win32':
        print("Operating system Windows detected.")
        print("sgDI-tector uses WSL (https://ubuntu.com/wsl) to call bwa,",
              "DI-tector and MAFFT, please check that WSL is installed",
              "and all necessary software (bwa, samtools, mafft)",
              "is in PATH within WSL. DI-tector (version 0.6) must be",
              "in the same folder of sgDI-tector, and named",
              "DI-tector_06.py.")
    elif sys.platform == 'linux':
        print("Please check that",
              "all necessary software (bwa, samtools, mafft)",
              "is in PATH. DI-tector (version 0.6) must be",
              "in the same folder of sgDI-tector, and named",
              "DI-tector_06.py.")
    else:
        sys.exit("Operating system not supported. Sorry!")
    print("###########################################################")
    

    # check if all required softwares are installed
    software_needed = ["bwa", "samtools", "mafft"]
    if sys.platform == 'win32':
        # required
        if not which("wsl"):
            print("WARNING: WSL is needed and seems not to be installed, likely to become a fatal error...")
        for soft in software_needed:
            proc = subprocess.Popen("wsl", stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)    
            wh = proc.communicate(input=str.encode("which " + soft))[0]
            if len(wh)==0:
                print("WARNING: " + soft + " is needed and seems not to be installed, likely to become a fatal error...")
    elif sys.platform == 'linux':
        for soft in software_needed:
            if not which(soft):
                print("WARNING: " + soft + " is needed and seems not to be installed, likely to become a fatal error...")

    # check if DI-tector must be run
    if args.Output_Directory:
        out_dir = args.Output_Directory
        if not out_dir.endswith("/"):
            out_path = out_dir+'/'
    else:
        out_path = "./"
    di_tect_outfiles = ['DI-tector_output_sorted.txt']
    if args.Force_DItector:
        run_ditector = True
    elif not Path(out_path).exists():
        run_ditector = True
    elif set(di_tect_outfiles) <= set(os.listdir(out_path)):
        run_ditector = False
    else:
        run_ditector = True
    
    if run_ditector:
        # print an error if DI-tector_06.py is not found
        assert Path("DI-tector_06.py").exists(), "ERROR: DI-tector_06.py needs to be in the same folder as sgDI-tector, exiting..."
        
        # create viral sequence index
        print("Indexing viral sequence...", end=' ', flush=True)
        t_command = "bwa index " + args.Virus_Ref
        if sys.platform == 'win32':
            proc = subprocess.Popen("wsl", stdin=subprocess.PIPE) # With Popen seems that .bashrc is loaded...
            proc.communicate(input=str.encode(t_command))
        elif sys.platform == 'linux':
            subprocess.run(shlex.split(t_command))
        assert Path(args.Virus_Ref+".bwt").exists(), "ERROR: viral sequence indexing not succesfull, exiting..."
        print("Done!")
        
        print("Running DI-tector...")
        print("This is DI-tector:")
        t_command = "python3 DI-tector_06.py " + ditect_options
        print(t_command)
        if sys.platform == 'win32':
            proc = subprocess.Popen("wsl", stdin=subprocess.PIPE)
            proc.communicate(input=str.encode(t_command))
        elif sys.platform == 'linux':
            subprocess.run(shlex.split(t_command))
        print("Done!")
    else:
        print("DI-tector output files found, not running DI-tector again (to change this, use the --Force_DItector option).")
    
    print("Loading DVG data...", end=' ', flush=True)
    ditect_outfile_outsor = out_path + 'DI-tector_output_sorted.txt'
    # extract deletions (fwd and rev strand), keep only bp, ri, sequence, group by them and add count column
    t_df = pd.read_csv(ditect_outfile_outsor, sep='\t')
    if len(t_df) == 0: # exit if DI-tector didn't find any DVG sequence
        print("No DVG data found from DI-tector output file! Exiting...")
        sys.exit(0)
    del_df = t_df[(t_df["DVG's type"] == 'Deletion DVG (Fwd. strand)') | (t_df["DVG's type"] == 'Deletion DVG (Rev. strand)')]
    del_small_df = del_df[["BP_Pos", "RI_Pos", "SEQ_FL_ori"]].copy() # here I am loosing info about fwd/rev strand deletions (but the sequences could be different...)
    del_small_df["count"] = 1
    del_small_df = del_small_df.groupby(by=["BP_Pos", "RI_Pos", "SEQ_FL_ori"], as_index=False).sum()
    # correct for ambiguous BP and RI pos (move as close to 3' as possible)
    t_seqs = []
    for record in SeqIO.parse(args.Virus_Ref, "fasta"):
        t_seqs.append(str(record.seq))
    assert len(t_seqs) == 1, "Please provide the viral sequence in FASTA format, as a single entry."
    vir_seq = t_seqs[0]
    rabps = []
    raris = []
    for i in range(len(del_small_df)):
        bp = del_small_df["BP_Pos"][i]
        ri = del_small_df["RI_Pos"][i]
        t_bp, t_ri = align_junction_LeftRight(vir_seq, bp, ri, "R")
        rabps.append(t_bp)
        raris.append(t_ri)
    del_small_df["BP_Pos"] = rabps
    del_small_df["RI_Pos"] = raris
    print("Done!")
    
    print("Searching for leader sequence...", end=' ', flush=True)
    # find end of leader sequence
    max_nbp_pos_list = find_top_BP_intervals(del_small_df, leader_window)
    # check that a single TRS-zone is found
    l_extra = 0
    for i in range(len(max_nbp_pos_list)-1):
        pos1 = max_nbp_pos_list[i]
        pos2 = max_nbp_pos_list[i+1]
        if abs(pos1 - pos2) > leader_window:
            print("WARNING: two separate TRS could be present here! I am considering only the first one...")
        else:
            l_extra += 1
    max_nbp_pos = max_nbp_pos_list[0]
    leader_window += l_extra # if several positions close together have the same number of junctions
    # print file with leader sequence
    with open(out_path + 'leader_sequence.txt','w') as f:
        f.write(vir_seq[:max_nbp_pos+leader_window])
    print("Done!")
        
    print("Filtering and clustering sequences...", end=' ', flush=True)
    # filter out everything with BP_Pos not in the interval, then cluster for ATG position
    del_small_filtered_df = del_small_df[(del_small_df["BP_Pos"] >= max_nbp_pos) & (del_small_df["BP_Pos"] < max_nbp_pos+leader_window)].reset_index(drop=True)
    del_small_filtered_families_df = cluster_subgenomic_families(del_small_filtered_df, vir_seq, key_bp="BP_Pos", key_ri="RI_Pos")
    print("Done!")
    
    print("Annotating ORF names...", end=' ', flush=True)
    # now annotate ORFs: take the coding sequences from each ATG to the first end codon, then align against each of the
    # reference cds (nucleotide) in the file, and put a threshold for the match.
    # collect all full translated sgRNAs
    hit_translated_sgRNAs = []
    for atgp in del_small_filtered_families_df["atg_pos"]:
        t_s = vir_seq[atgp-1:]
        t_s = t_s[:int((len(t_s) // 3) * 3)]
        t_cdna = Seq(t_s)
        hit_translated_sgRNAs.append(t_cdna.translate(to_stop=True))
    # add names to unknown families - all families which have no name will be called with the position of the ORF
    atg_startpos = [i+1 for i in findall('ATG', vir_seq)]
    startpos_to_exaname = ['U'+hex(x)[2:].rjust(4, '0') for x in atg_startpos]
    nams = [startpos_to_exaname[atg_startpos.index(sp)] for sp in del_small_filtered_families_df["atg_pos"]]
    del_small_filtered_families_df["ORF"] = nams
    # give names to known families - the highest score in the alignment will be called ORF x, the second ORF x.1 and so on,
    # for all the scores above the threshold
    if args.sgRNA_Ref:
        # collect all full translated standard sgORFs
        subgen_ref_nams = []
        subgen_ref_trans_seqs = []
        for record in SeqIO.parse(args.sgRNA_Ref, "fasta"):
            subgen_ref_nams.append(str(record.description))
            t_cdna = Seq(str(record.seq))
            subgen_ref_trans_seqs.append(t_cdna.translate(to_stop=True))    
        for st_nam, st_prot in zip(subgen_ref_nams, subgen_ref_trans_seqs):
            t_scores = []
            for i, t_prot in enumerate(hit_translated_sgRNAs):
                tt_score = pairwise2.align.globalxx(st_prot, t_prot, score_only=True) / max(len(st_prot), len(t_prot))
                if tt_score >= threshold_sgRNAnames:
                    t_scores.append([i, tt_score])
            t_scores = np.array(t_scores)
            if len(t_scores)==0:
                continue
            t_scores = t_scores[t_scores[:, 1].argsort()[::-1]]
            for i, pos in enumerate(t_scores[:, 0]):
                if i == 0:
                    del_small_filtered_families_df.at[int(pos), "ORF"] = st_nam
                else:
                    del_small_filtered_families_df.at[int(pos), "ORF"] = st_nam + "-" + str(i)
    else:
        print("(no reference sgRNAs provided, so sgDI-tector will not annotate junctions with standard names)", end=' ', flush=True)
    # print file with sgRNA dataframe
    del_small_filtered_families_df.to_csv(out_path + "sgRNA_junctions.csv", index=False)
    print("Done!")
    
    print("Alignment of sgRNA junctions...")
    print("This is MAFFT:")
    # last part: align and print the logo of the sequences
    # write temp file in fasta format with leader_window
    junc_seq = [vir_seq[rip-junct_window:rip]+vir_seq[rip:rip+junct_window] for rip in del_small_filtered_df["RI_Pos"]]
    with open('temp.fasta','w') as f:
        for i, s in enumerate(junc_seq):
            f.write('> junction' + str(i) + '\n')
            f.write(s + '\n')
    # align all seqs with MAFFT fftnsi routine (https://mafft.cbrc.jp/alignment/software/manual/manual.html)
    t_command = "einsi --thread " + str(args.Nb_threads) + " temp.fasta > temp_ali.fasta"
    if sys.platform == 'win32':
        proc = subprocess.Popen("wsl", stdin=subprocess.PIPE) # This trick make sure that .bashrc is loaded!
        proc.communicate(input=str.encode(t_command))
    elif sys.platform == 'linux':
        with open("temp_ali.fasta", 'w') as mafft_outtempfile: 
            subprocess.run(shlex.split(t_command)[:-2], stdout=mafft_outtempfile, text=True)
    assert Path("temp_ali.fasta").exists(), "ERROR: mafft alignment not succesfull, exiting..."
    # load MSA, trim columns with fraction of gaps > gap_thresh
    MSA_junctions = [str(record.seq).upper() for record in SeqIO.parse("temp_ali.fasta", "fasta")]
    RNA_MSA_junctions = [s.replace("T","U") for s in MSA_junctions]
    RNA_MSA_junc_2num = seq2num(RNA_MSA_junctions)
    trimmed_RNA_MSA_junc2num = RNA_MSA_junc_2num[:, np.sum(RNA_MSA_junc_2num == 4, axis=0) <= len(RNA_MSA_junc_2num) * gap_thresh]
    # sequence logo - not weighted
    mu = average(trimmed_RNA_MSA_junc2num, q=5)
    fig = RNA_sequence_logo.Sequence_logo(mu, ticks_every=5, figsize=(12,5), ticks_labels_size=18, title_size=22, show=False)
    # remove temporary files
    try:
        os.remove("temp.fasta")
    except OSError:
        pass
    try:
        os.remove("temp_ali.fasta")
    except OSError:
        pass
    # save figure with logo
    fig.savefig(out_path + 'TRS_logo.pdf', format='pdf', bbox_inches='tight')
    print("Done!")
    
    print("All done!")
