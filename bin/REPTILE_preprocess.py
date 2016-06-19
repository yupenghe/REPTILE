#!/usr/bin/env python

######################################################################################
## REPTILE
## - Regulatory Element Prediction based on TIssue-specific Local Epigenetic marks
## by Yupeng He (yupeng.he.bioinfo at gmail)
## Salk Institute for Biological Studies
######################################################################################
def print_error(error_message=""):
    sys.stderr.write(error_message)
    sys.exit(1)

def get_epimark_profile(query_bed_file,epimark_id,bw_file,output_filename):
    intervals = pd.read_table(query_bed_file,
                              names=["chr","start","end","id"])

    cmd = " ".join(["/usr/bin/env",
                    "bigWigAverageOverBed",
                    bw_file,
                    query_bed_file,
                    "stdout"])

    ## The code regarding subprocess.Popen is based from oarevalo's answer on:
    ## http://stackoverflow.com/questions/16198546/get-exit-code-and-stderr-from-subprocess-call
    pipes = subprocess.Popen(shlex.split(cmd),
                             stdout=subprocess.PIPE,
                             stderr=subprocess.PIPE,
                             universal_newlines=True)
    std_out, std_err = pipes.communicate()
    if pipes.returncode != 0:
        print_error(std_err)
    score = pd.read_table(StringIO(std_out),
                          names=["id","size","covered","sum","mean0",mark + "_" + sample]
    )
     
    #score = pd.read_table(
    #    StringIO(subprocess.check_output(shlex.split(cmd),
    #                                     universal_newlines=True)),
    #    names=["id","size","covered","sum","mean0",mark + "_" + sample]
    #)
        
    #intervals = intervals["id"]
    score = pd.merge(intervals,score,on="id")
    del intervals
    #score = score[["id",mark + "_" + sample]]
    score = score[mark + "_" + sample]
    score.to_csv(output_filename,sep="\t",header=False,index=False)
    del score

def merge_split_epimark_files(query_bed_file,epimark_file_prefix,
                              epimark_id_list,header,output_filename):
    outhandle = open(output_filename,'w')
    outhandle.write(header)
    num_line = 0
    inhandle_list = []
    for epimark_id in epimark_id_list:
        inhandle_list.append(open(epimark_file_prefix+str(epimark_id)+".txt",'r'))
    fhandle = open(query_bed_file,'r')
    for line in fhandle:
        num_line += 1
        line = line.rstrip()
        for i in range(len(inhandle_list)):
            line = line + "\t" + inhandle_list[i].readline().rstrip()    
        outhandle.write(line+"\n")
    for i in range(len(inhandle_list)):
        inhandle_list[i].close()
    outhandle.close()
    return(num_line)
    
######################################################################################
#Usage information
import argparse
import sys

description ="""Preprocess data into input format for REPTILE. Please email 
Yupeng He (yupeng.he.bioinfo at gmail) for feedback, question and bug.
"""
parser = argparse.ArgumentParser(description=description,
                                 formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument("data_info_file",
                    help="tab-separated file providing information about samples,\n"
                    + "epigenetic marks and path to corresponding bigwig files.\n"
                    + "No duplicated mark name is allowed for each sample and\n"
                    + "no duplicated sample name is allowed for each mark.\n"
                    + "format:\n"
                    + "  <sample name><\\t><mark name><\\t><path to bigwig file>\n"
                    + "example (header is required):\n"
                    + "  sample mark bw_file\n"
                    + "  heart H3K4me1 bigwig/heart_H3K4me1.bw\n"
                    + "  ...\n"
                    + "  brain H3K27ac bigwig/brain_H3K27ac.bw\n"
)
parser.add_argument("query_region_file",
                    help="bed file of regions used for training and/or prediction,\n"
                    + "format:\n"
                    + "  <chromosome><\\t><start><\\t><end>\\t<region id>\n"
                    + "example:\n"
                    + "  chr1 3172000 3173000 reg_0\n"
                    + "  ...\n"
                    + "  chr9 124412000 124413000 reg_302031\n"
)
parser.add_argument("output_prefix",
                    help="prefix of output files\n"
                    +"example: heart_enhancer\n"
                    +"  two output files will be:\n"
                    +"  - heart_enhancer.DMR_with_epimark.tsv\n"
                    +"  - heart_enhancer.region_with_epimark.tsv\n"
)
parser.add_argument("-d","--DMR-file",
                    type=str,
                    default=None,
                    help="bed file specifying differentially methylated regions (DMRs)\n"
                    + "DMRs are used as high-resolution enhancer candidates to \n"
                    + "increase the resolution of training and prediction.\n"
                    + "File can also contain high-resolution candidate loci from other\n"
                    + "assays such DNase-seq and ATAC-seq, or loci from motif analysis.\n"
                    + "format:\n"
                    + "  <chromosome><\\t><start><\\t><end>\\t<region id>\n"
                    + "example:\n"
                    + "  chr1 3172266 3172488 dmr_0\n"
                    + "  ...\n"
                    + "  chr19 61316546 61316778 dmr_513260\n"
)
parser.add_argument("-f","--DMR-overlap-fraction",
                    type=float,
                    default=1.0,
                    help="Minimum fraction of DMR in overlapped interval required to\n"
                    + "consider the DMR to be overlapping with region.\n"
                    + "default: 1.0"
)

parser.add_argument("-e","--extending-DMR",
                    type=int,
                    default=0,
                    help="Based to extend for DMRs to include flanking regions in \n"
                    + "calculation of the signal for epigenetic marks. It is ususally\n"
                    + "helpful because DMRs are often nucleosome free and informative\n"
                    + "informative epigenetic marks locate in flanking regions.\n"
                    + "Start of DMRs will be subtracted by EXTENDING_DMR and the end\n"
                    + "will be increase by EXTENDING_DMR in epigenetic mark signal\n"
                    + "calculation. EXTENDING_DMR can be negative if you like\n"
                    + "only consider center regions of DMRs. Be careful about choosing\n"
                    + "value in that situation. The length of some DMRs may because negative\n"
                    + "and if so, this program will stop and pump out an error\n"
                    + "This option does not affect overlap between DMRs and query regions.\n"
                    + "default: 0"
)

parser.add_argument("-p","--num-process",
                    type=int,
                    default=1,
                    help="Number of threads/processes"
)

parser.add_argument("-n","--num-output-file",
                    type=int,
                    default=1,
                    help="If greater than 1, split output files equally. \n"
)

parser.add_argument("-g","--for-genome-wide-prediction",
                    action='store_true',
                    help="This option is used for generating the input files used for genome-wide\n"
                    + "prediction. Please don't use this option when trying to generate inputs for\n"
                    + "model training or prediction on given query regions.\n"
                    + "  1. When this option is used, this program will not try to overlap DMRs\n"
                    + "with query regions. Instead, it will just output the values epigenetic\n"
                    + "marks of query regions and DMRs.\n"
                    + "  2. When this option is not used, DMRs will be overlapped with query\n"
                    + "regions and the values of epigenetic marks will be calculated for the DMR\n"
                    + "in each overlap instance. If a DMR overlaps with multiple query regions,\n"
                    + "its values of epigenetic marks will be calculated multiple times. In the\n"
                    + "output file(s), the DMR id will be merged with the id of overlapping query\n"
                    + "region to specify the overlapping relationship. For example:\n"
                    + "  dmr_23:query_reg_33\n"
                    + "  dmr_23:query_reg_34\n"
                    + "  dmr_23:query_reg_35\n"
                    + "  dmr_24:query_reg_36\n"
                    + "  dmr_24:query_reg_37\n"
                    + "The overlapping information is required for training. The DMRs with no\n"
                    + "overlapping query regions will be ignored.\n"
)


try:
    options = parser.parse_args()
except:
    print("\nAdd option --help or -h to get help information\n")
    sys.exit(1)

#Check whether fraction (-f) is within range
if options.DMR_overlap_fraction > 1.0 or options.DMR_overlap_fraction < 0.0:
    print("Error: Invalid value for -f: " + str(options.DMR_overlap_fraction))
    print("       Should be between 0.0 and 1.0")

######################################################################################
#Import modules
    
import multiprocessing
import subprocess
import shlex
import os
import pandas as pd
try:
    from StringIO import StringIO
except:
    from io import StringIO
    
######################################################################################
#Read sample information
try:
    data_info = pd.read_table(options.data_info_file)
except:
    print_error("Error! Cannot open file \"" + options.data_info_file + "\"!\n")

#Check if there are duplicated sample and mark combination
if sum(data_info.duplicated(["sample","mark"])) > 0:
    print(data_info.ix[data_info.duplicated(["sample","mark"])])
    print_error("Error! Found duplicated sample and mark combination!\n")

#Check if every bigwig file exists
for (sample,mark,bw_file) in data_info.itertuples(index=False):
    if not (os.path.isfile(bw_file)):
        print_error("Error! bigwig file \"" + bw_file + "\" not exist!\n")

######################################################################################
if options.DMR_file is not None:
    if not options.for_genome_wide_prediction:
        cmd = " ".join(["/usr/bin/env",
                        "bedtools",
                        "intersect",
                        "-a",options.DMR_file,
                        "-b",options.query_region_file,
                        "-wa -wb",
                        "-f",str(options.DMR_overlap_fraction)])

        #Get DMRs that overlap regions
        ## The code regarding subprocess.Popen is based from oarevalo's answer on:
        ## http://stackoverflow.com/questions/16198546/get-exit-code-and-stderr-from-subprocess-call
        pipes = subprocess.Popen(shlex.split(cmd),
                                 stdout=subprocess.PIPE,
                                 stderr=subprocess.PIPE,
                                 universal_newlines=True)
        std_out, std_err = pipes.communicate()
        if pipes.returncode != 0:
            print_error(std_err)

        DMR = pd.read_table(StringIO(std_out),
                            names = ["dmr_chr","dmr_start","dmr_end","dmr_id",
                                     "reg_chr","reg_start","reg_end","reg_id"])

        DMR["combined_id"] = DMR["dmr_id"]+":"+DMR["reg_id"]
        DMR = DMR[["dmr_chr","dmr_start","dmr_end","combined_id"]]
    else:
        DMR = pd.read_table(options.DMR_file,
                            names = ["dmr_chr","dmr_start","dmr_end","dmr_id"])

    DMR["dmr_start"] = DMR["dmr_start"] - options.extending_DMR
    DMR.loc[DMR["dmr_start"] < 0,"dmr_start"] = 0
    DMR["dmr_end"] = DMR["dmr_end"] + options.extending_DMR
    DMR.to_csv(options.output_prefix+".DMR.bed",sep="\t",header=False,index=False)
    del DMR
######################################################################################
#
#print(min(options.num_process,data_info.shape[0]))
pool = multiprocessing.Pool(min(options.num_process,data_info.shape[0]))

for (sample,mark,bw_file) in data_info.itertuples(index=False):
    epimark_id = mark + "_" + sample
    pool.apply_async(
        get_epimark_profile,
        (options.query_region_file,
         epimark_id,
         bw_file,
         options.output_prefix+".region."+str(epimark_id)+".txt"),
    )
    if options.DMR_file is not None:
        pool.apply_async(
            get_epimark_profile,
            (options.output_prefix+".DMR.bed",
             epimark_id,
             bw_file,
             options.output_prefix+".DMR."+str(epimark_id)+".txt"),
        )


pool.close()
pool.join()
######################################################################################
# Merging
## get header
epimark_id_list = [mark+"_"+sample for (sample,mark,bw_file) in data_info.itertuples(index=False)]
header = "\t".join(["chr","start","end","id"]+epimark_id_list) + "\n"

## begin merging
pool = multiprocessing.Pool(min(options.num_process,2))

num_region = pool.apply_async(
    merge_split_epimark_files,
    (options.query_region_file,
     options.output_prefix+".region.",
     epimark_id_list,
     header,
     options.output_prefix+".region_with_epimark.tsv"),
)

if options.DMR_file is not None:
    num_DMR = pool.apply_async(
        merge_split_epimark_files,
        (options.output_prefix+".DMR.bed",
         options.output_prefix+".DMR.",
         epimark_id_list,
         header,
         options.output_prefix+".DMR_with_epimark.tsv"),
    )

pool.close()
pool.join()
num_region = num_region.get(timeout=1)

if options.DMR_file is not None:
    num_DMR = num_DMR.get(timeout=1)
    subprocess.check_call(["rm"]+
                          [options.output_prefix+".DMR.bed"]+
                          [options.output_prefix+".DMR."+str(epimark_id)+".txt" for epimark_id in epimark_id_list]+
                          [options.output_prefix+".region."+str(epimark_id)+".txt" for epimark_id in epimark_id_list])
else:
    subprocess.check_call(["rm"]+
                          [options.output_prefix+".region."+str(epimark_id)+".txt" for epimark_id in epimark_id_list])

######################################################################################
# Split files
if options.num_output_file > 1:
    out_line_cutoff = int(float(num_region)/float(options.num_output_file))+1
    out_line_count = 0
    out_file_id = 0
    fhandle = open(options.output_prefix+".region_with_epimark.tsv",'r')
    fhandle.readline()
    outhandle = open(options.output_prefix+".region_with_epimark."+str(out_file_id)+".tsv",'w')
    outhandle.write(header)
    regid2fileid = {}
    for line in fhandle:
        if out_line_count >= out_line_cutoff:
            outhandle.close()
            out_line_count = 0
            out_file_id += 1
            outhandle = open(options.output_prefix+".region_with_epimark."+str(out_file_id)+".tsv",'w')
            outhandle.write(header)
        outhandle.write(line)
        out_line_count += 1
        regid2fileid[line.split("\t")[3]] = out_file_id
        
    outhandle.close()
    fhandle.close()        

    ## DMR
    if options.DMR_file is not None:
        if not options.for_genome_wide_prediction:
            fhandle = open(options.output_prefix+".DMR_with_epimark.tsv",'r')
            fhandle.readline() ## Skip header
            out_line_count = 0
            outhandles = []
            for out_file_id in range(options.num_output_file):
                outhandles.append(open(options.output_prefix+".DMR_with_epimark."+str(out_file_id)+".tsv",'w'))
                outhandles[-1].write(header)
        
            for line in fhandle:
                region_id = line.split("\t")[3].split(":")[1]
                outhandles[regid2fileid[region_id]].write(line)
            fhandle.close()
    
            for i in range(len(outhandles)):
                outhandles[i].close()

        else:
            ## For genome-wide enhancer prediction
            ## DMR (without link to query regions)
            out_line_cutoff = int(float(num_DMR)/float(options.num_output_file))+1
            out_line_count = 0
            out_file_id = 0
            fhandle = open(options.output_prefix+".DMR_with_epimark.tsv",'r')
            fhandle.readline()
            outhandle = open(options.output_prefix+".DMR_with_epimark."+str(out_file_id)+".tsv",'w')
            outhandle.write(header)
            regid2fileid = {}
            for line in fhandle:
                if out_line_count >= out_line_cutoff:
                    outhandle.close()
                    out_line_count = 0
                    out_file_id += 1
                    outhandle = open(options.output_prefix+".DMR_with_epimark."+str(out_file_id)+".tsv",'w')
                    outhandle.write(header)
                outhandle.write(line)
                out_line_count += 1
            
            outhandle.close()
            fhandle.close()
