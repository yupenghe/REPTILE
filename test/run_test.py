#!/usr/bin/env python

######################################################################################
## REPTILE
## - Regulatory Element Prediction based on TIssue-specific Local Epigenetic marks
## by Yupeng He (yupeng.he.bioinfo at gmail)
## Salk Institute for Biological Studies
######################################################################################

import subprocess
import shlex
import sys

def print_error(error_message = ""):
    sys.stderr.write(error_message)
    sys.exit()

def test_command_list(cmds,cmd_names):
    for ind in range(len(cmds)):
        pipes = subprocess.Popen(shlex.split(cmds[ind]),
                                 stdout = subprocess.PIPE,
                                 stderr = subprocess.PIPE,
                                 universal_newlines = True)
        std_out, std_err = pipes.communicate()
        if pipes.returncode != 0:
            print_error(
                " - " + 
                cmd_names[ind] +
                ": Failed!\n\n" +
                "Command:\n" +
                cmds[ind] +
                "\n\n" +
                "Error message:\n" +
                std_err)
        else:
            print(" - " + cmd_names[ind] + ": Pass!")    

######################################################################################
# Create folder for temoprary file
subprocess.check_call(shlex.split("mkdir -p tmp/"))
print("Test start!")

# Test 0
print("\nTesting REPTILE requirements")
## Tests bedtools and bigWigAverageOverBed
test_bed_file = "data/DMR_CG_mESC_E11_5_ext150.bed"
test_bw_file = "data/bw/mESC_Meth.bw"
cmds = [
    # bedtools without parameters
    " ".join(["/usr/bin/env",
              "bedtools"]),
    # bedtools intersect
    " ".join(["/usr/bin/env",
              "bedtools",
              "intersect",
              "-a",test_bed_file,
              "-b",test_bed_file,
              "-wa -wb"]),
    # bedtools merge
    " ".join(["/usr/bin/env",
              "bedtools",
              "merge",
              "-c 4",
              "-o collapse,max",
              "-i",test_bed_file]),
    # bedtools sort
    " ".join(["/usr/bin/env",
              "bedtools",
              "sort",
              "-i",
              test_bed_file]),
    # bigWigAverageOverBed
    " ".join(["/usr/bin/env",
              "bigWigAverageOverBed",
              test_bw_file,
              test_bed_file,
              "stdout"])

]
cmd_names = ["run bedtools",
             "bedtools interesct",
             "bedtools merge",
             "bedtools sort",
             "bigWigAverageOverBed"]

test_command_list(cmds,cmd_names)

######################################################################################        
# Test 1.5
## numpy and pandas
try:
    import numpy
    print(" - numpy import: Pass!")    
except:
    print_error(
        " - " + 
        "numpy import" +
        ": Failed!\n\n" +
        "Error message:\n" +
        str(sys.exc_info()[1]) +
        "\n"
    )
    
try:
    import pandas
    print(" - pandas import: Pass!")    
except:
    print_error(
        " - " + 
        "pandas import" +
        ": Failed!\n\n" +
        "Error message:\n" +
        str(sys.exc_info()[1]) +
        "\n"
    )



######################################################################################        
# Test 1
## REPTILE_preprocess.py
print("\nTesting script for preprocessing (REPTILE_preprocess.py)")
dmr_file = "data/DMR_CG_mESC_E11_5_ext150.bed"
data_info_file = "data/data_info_mESC.tsv"
test_region = "data/test_data/test_region.bed"
output_prefix = "tmp/dbug_preproces"
## Preprocessing training data
cmds = [
    " ".join(["../bin/REPTILE_preprocess.py", data_info_file,
              test_region, "-d", dmr_file, output_prefix]),
    " ".join(["../bin/REPTILE_preprocess.py", data_info_file,
              test_region, "-d", dmr_file, output_prefix, "-g"])
]
cmd_names = ["REPTILE_preprocess.py",
             "REPTILE_preprocess.py (-g)"]

test_command_list(cmds,cmd_names)

######################################################################################        
# Test 2
## REPTILE_train.R
print("\nTesting script for training enhancer model (REPTILE_train.R)")
data_info_file = "data/test_data/data_info_mESC_E11_5.tsv"
train_region = "data/test_data/test_region.region_with_epimark.tsv"
train_DMR = "data/test_data/test_region.DMR_with_epimark.tsv"
train_label = "data/test_data/test_region_label.tsv"
output_prefix = "tmp/dbug_enhancer_model"
## Preprocessing training data
cmds = [
    " ".join(["../bin/REPTILE_train.R",
              "-i", data_info_file,
              "-a", train_region,
              "-d", train_DMR,
              "-l", train_label,
              "-s", "E11_5_FB",
              "-o", output_prefix])
]
cmd_names = ["REPTILE_train.py"]
test_command_list(cmds,cmd_names)

######################################################################################        
# Test 3
## REPTILE_compute_score.R
print("\nTesting script for computing enhancer scroe (REPTILE_compute_score.R)")
data_info_file = "data/test_data/data_info_mESC_E11_5.tsv"
model_file = "tmp/dbug_enhancer_model.reptile"
test_region = "data/test_data/test_region.region_with_epimark.tsv"
test_DMR = "data/test_data/test_region.DMR_with_epimark.tsv"
test_label = "data/test_data/test_region_label.tsv"
output_prefix = "tmp/dbg_test_region_score"
## Preprocessing training data
cmds = [
    " ".join(["../bin/REPTILE_compute_score.R",
              "-i", data_info_file,
              "-a", test_region,
              "-d", test_DMR,
              "-m", model_file,
              "-s", "E11_5_HT",
              "-o", output_prefix,
              "-w"]),
    " ".join(["../bin/REPTILE_compute_score.R",
              "-i", data_info_file,
              "-a", test_region,
              "-d", test_DMR,
              "-m", model_file,
              "-s", "E11_5_HT",
              "-o", output_prefix])

]
cmd_names = ["REPTILE_compute_score.R (-w)",
             "REPTILE_compute_score.R"]
test_command_list(cmds,cmd_names)

######################################################################################        
# Test 4
## REPTILE_call_enhancer.py
print("\nTesting script for calling enhancers (REPTILE_call_enhancer.py)")
region_score = "tmp/dbg_test_region_score.R.bed"
DMR_score = "tmp/dbg_test_region_score.DMR.bed"
output_file = "tmp/dbg_test_region_score"
## Preprocessing training data
cmds = [
    " ".join(["../bin/REPTILE_call_enhancer.py",              
              region_score,
              "-d", DMR_score,
              "-o",
              output_file])
]
cmd_names = ["REPTILE_call_enhancer.R"]
test_command_list(cmds,cmd_names)

######################################################################################        
# Test 5
## REPTILE_evaluate_prediction.py
print("\nTesting script for evaluating predictions (REPTILE_evaluate_prediction.R)")
region_score = "tmp/dbg_test_region_score.D.bed"
test_label = "data/test_data/test_region_label.tsv"
## Preprocessing training data
cmds = [
    " ".join(["../bin/REPTILE_evaluate_prediction.R",
              "-s", "E11_5_HT",
              "-p", region_score,
              "-l", test_label])
]
cmd_names = ["REPTILE_evaluate_prediction.R"]
test_command_list(cmds,cmd_names)

######################################################################################
# Remove temporary files
subprocess.check_call(shlex.split("rm -r tmp/"))
print("\nDone!")
