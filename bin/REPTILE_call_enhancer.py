#!/usr/bin/env python

######################################################################################
## REPTILE
## - Regulatory Element Prediction based on TIssue-specific Local Epigenetic marks
## by Yupeng He (yupeng.he.bioinfo at gmail)
## Salk Institute for Biological Studies
######################################################################################

def print_error(error_message=""):
    sys.stderr.write(error_message)
    sys.exit()

def get_peaks_from_block(block_bins,peak_dist):
    peak_inds = []
    scores = [ tmp_bin[4] for tmp_bin in block_bins ]
    max_score = 0
    while True:
        peak_ind = argmax(scores)
        max_score = scores[peak_ind]
        if max_score == -1:
            break
        peak_inds.append(peak_ind)
        ## Set scores of bins adjacent to peak to be -1
        scores[peak_ind] = -1
        for ind in range(len(block_bins)):
            if scores[ind] == -1:
                continue
            elif not ((int(block_bins[ind][1]) - int(block_bins[peak_ind][2])) > peak_dist or\
                      (int(block_bins[peak_ind][1]) - int(block_bins[ind][2])) > peak_dist):
                scores[ind] = -1
    return(peak_inds)
    
######################################################################################
#Usage information
import sys
import argparse

description ="""Used to call enhancers based on REPTILE enhancer confidence scores 
of DMRs and sliding windows. Please email Yupeng He (yupeng.he.bioinfo at gmail) for 
feedbacks, questions or bugs.
"""
parser = argparse.ArgumentParser(description=description,
                                 formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument("query_region_file",
                    help="bed file of query regions with enhancer scores such as the\n"
                    + "*.R.bed files from REPTILE_compute_score.R.\n"
                    + "format:\n"
                    + "  <chromosome><\\t><start><\\t><end>\\t<region id>\n"
                    + "example:\n"
                    + "  chr1 3172000 3173000 reg_0\n"
                    + "  ...\n"
                    + "  chr9 124412000 124413000 reg_302031\n"
)
parser.add_argument("-d","--DMR-file",
                    type=str,
                    default=None,
                    help="bed file of DMRs with enhancer scores such as the\n"
                    + "*.DMR.bed file from REPTILE_compute_score.R.\n"
                    + "format:\n"
                    + "  <chromosome><\\t><start><\\t><end>\\t<DMR id>\n"
                    + "example:\n"
                    + "  chr1 5143000 5143500 dmr_0\n"
                    + "  ...\n"
                    + "  chr19 4412008 4412403 dmr_2052\n"
)
parser.add_argument("-o", "--output-file",
                    type=str,
                    default=None,
                    help="Name of output BED file of predicted enhancers and their scores.\n")

parser.add_argument("-p","--cutoff",
                    type=float,
                    default=0.5,
                    help="Minimum enhancer confidence score for a sliding window\n"
                    + "to be consider as a peak\n"
                    + "default: 0.5\n"
)
parser.add_argument("-s","--peak-spacing",
                    type=int,
                    default=100,
                    help="Minimum distance for two sliding windows\n"
                    + "to be called as separate peaks\n"
)


try:
    options = parser.parse_args()
except:
    print("\nAdd option --help or -h to get help information\n")
    sys.exit()


from numpy import argmax
import pandas as pd
import math
import os
import subprocess
import shlex

peak_dist = options.peak_spacing

tmp_files = []


if options.DMR_file is not None:
    # Filter DMRs by enhancer score
    DMR = pd.read_table(options.DMR_file,
                        names = ["chr","start","end","id","score"])
    DMR = DMR.ix[DMR["score"]>=options.cutoff]
    DMR.to_csv(options.DMR_file+".filtered",sep="\t",header=False,index=False)
    tmp_files.append(options.DMR_file+".filtered")
    del DMR

    # Sort filtered DMRs
    cmd = " ".join(["/usr/bin/env",
                    "bedtools",
                    "sort",
                    "-i",
                    options.DMR_file+".filtered"])
    ohandle = open(options.DMR_file+".filtered.sorted",'w')
    tmp_files.append(options.DMR_file+".filtered.sorted")
    subprocess.check_call(shlex.split(cmd),stdout=ohandle)
    ohandle.close()

    # Merge DMRs that pass the filter
    cmd = " ".join(["/usr/bin/env",
                    "bedtools",
                    "merge",
                    "-c 4,5",
                    "-o collapse,max",
                    "-i",options.DMR_file+".filtered.sorted"])
    ohandle = open(options.DMR_file+".filtered.sorted.merged",'w')
    tmp_files.append(options.DMR_file+".filtered.sorted.merged")
    subprocess.check_call(shlex.split(cmd),stdout=ohandle)
    ohandle.close()

# Filter query regions by score
fhandle = open(options.query_region_file,'r')
ohandle = open(options.query_region_file+".filtered",'w')
tmp_files.append(options.query_region_file+".filtered")
for line in fhandle:
    line = line.rstrip()
    (chrom,start,end,reg_id,score) = line.split("\t")
    score = float(score)
    if score >= options.cutoff:
        ohandle.write(line+"\n")
fhandle.close()
ohandle.close()

# Sort
cmd = " ".join(["/usr/bin/env",
                "bedtools",
                "sort",
                "-i",
                options.query_region_file+".filtered"])
ohandle = open(options.query_region_file+".filtered.sorted",'w')
tmp_files.append(options.query_region_file+".filtered.sorted")
subprocess.check_call(shlex.split(cmd),stdout=ohandle)
ohandle.close()

# Peak calling
fhandle = open(options.query_region_file+".filtered.sorted",'r')
ohandle = open(options.query_region_file+".peaks",'w')
tmp_files.append(options.query_region_file+".peaks")
cur_chrom=""
block_bins=[]
block_bin_end = 0

for line in fhandle:
    line = line.rstrip()
    (chrom,start,end,reg_id,score) = line.split("\t")
    start = int(start)
    end = int(end)
    score = float(score)
    
    if chrom != cur_chrom:
        if len(block_bins) > 0:
            peak_inds = get_peaks_from_block(block_bins,peak_dist)
            for ind in peak_inds:
                ohandle.write("\t".join(map(str,block_bins[ind]))+"\n")
        block_bins = []
        block_bin_end = 0
        cur_chrom = chrom
        
    if len(block_bins) > 0:
        if start - block_bin_end < peak_dist:
            block_bins.append([chrom,start,end,reg_id,score])
            block_bin_end = end
        else:
            peak_inds = get_peaks_from_block(block_bins,peak_dist)
            for ind in peak_inds:
                ohandle.write("\t".join(map(str,block_bins[ind]))+"\n")
            block_bins = [[chrom,start,end,reg_id,score]]
            block_bin_end = end
    else:
        block_bins = [[chrom,start,end,reg_id,score]]
        block_bin_end = end

if len(block_bins) > 0: ## Boundary condition
    peak_inds = get_peaks_from_block(block_bins,peak_dist)
    for ind in peak_inds:
        ohandle.write("\t".join(map(str,block_bins[ind]))+"\n")

fhandle.close()
ohandle.close()


if options.DMR_file is not None:
    # Get peaks that do not overlap any DMRs
    ohandle = open(options.query_region_file+".peaks.DMR",'w')
    tmp_files.append(options.query_region_file+".peaks.DMR")
    cmd = " ".join(["/usr/bin/env",
                    "bedtools",
                    "intersect",
                    "-a",options.query_region_file+".peaks",
                    "-b",options.DMR_file+".filtered.sorted.merged",
                    "-v"])
    subprocess.check_call(shlex.split(cmd),stdout=ohandle)
    
    # Final results
    fhandle = open(options.DMR_file+".filtered.sorted.merged",'r')
    ohandle.write(fhandle.read())
    ohandle.close()

    # Sort
    cmd = " ".join(["/usr/bin/env",
                    "bedtools",
                    "sort",
                    "-i",
                    options.query_region_file+".peaks.DMR"])
    if options.output_file is not None:
        ohandle = open(options.output_file,'w')
        subprocess.check_call(shlex.split(cmd),stdout=ohandle)
        ohandle.close()
    else:
        subprocess.check_call(shlex.split(cmd))
else:    
    cmd = " ".join(["/usr/bin/env",
                    "bedtools",
                    "sort",
                    "-i",
                    options.query_region_file+".peaks"])
    if options.output_file is not None:
        ohandle = open(options.output_file,'w')
        subprocess.check_call(shlex.split(cmd),stdout=ohandle)
        ohandle.close()
    else:
        subprocess.check_call(shlex.split(cmd))

        
# Remove temporary files
subprocess.check_call(["rm"]+tmp_files)
