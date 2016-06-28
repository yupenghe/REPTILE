# REPTILE
### - Regulatory Element Prediction based on TIssue-specific Local Epigenetic marks

https://github.com/yupenghe/REPTILE

REPTILE is a tool to identify the precise location of enhancers by integrating histone modification data 
and base-resolution DNA methylation profiles.

Please contact [Yupeng He](mailto:yupeng.he.bioinfo@gmail.com) for for feedbacks, questions or bugs.

Please cite "Yupeng He, David U. Gorkin, Joseph R. Nery, Rosa Castanon, Ah Young Lee, Yin Shen, Axel Visel, Len A. Pennacchio, Bing Ren, Joseph R. Ecker, REPTILE: Regulatory Element Prediction based on TIssue-specific Local Epigenetic marks, In preparation".

## Overview
The document includes an example of using REPTILE to predict enhancers in mouse tissues. 
This example is based on the exact dataset used in the REPTILE manuscript. First, REPTILE
will be run to learn an enhancer model based on the data in mouse embryonic stem cells (mESCs).
In the training dataset, EP300 binding sites are used as representative active enhancers and
the promoters regions along with random genomic intervals are used as inative instances. In the
next step, the enhancer model learned in mESCs will be applied to predict the enhancer activity of
545 genomic elements in the heart tissue from mouse embryo at E11.5 developmental stages (E11_5_HT).
These elements have been experimentally validated using transgenic reporter assay. We will compare
the predictions with experimental results to evaluate the accuracy of REPTILE. Lastly, we will use
REPTILE to generate enhancer predictions across the genome.

The dataset used in the example includes DNA methylation (Meth) and 6 histone modifications (H3K4me1, H3K4me2, H3K4me3, H3K27ac, H3K27me3 and H3K36me3) of in total nine samples. These samples are mESCs and 8 mouse tissues from embryo at E11.5 developmental stage:
* E11_5_HT: E11.5 heart
* E11_5_LM: E11.5 limb
* E11_5_FB: E11.5 forebrain
* E11_5_MB: E11.5 midbrain
* E11_5_HB: E11.5 hindbrain
* E11_5_NT: E11.5 neural tube
* E11_5_CF: E11.5 embryonic facial prominence
* E11_5_LV: E11.5 liver

In addition, the dataset contains 40,000 genomic regions for training REPTILE and 545 genmoic regions
to evaluation its prediction accurary. Differetially methylated regions (DMRs), which are needed to improve the resolution of REPTILE predictions, were called by comparing the DNA methylation profile of all 9 samples.


## Requirements
It is recommended to run this example in a computer cluster (server) because of the runtime and the requirement for memory.
Finishing the entire example requires minimum 6 Gb memory (single CPU) and 60 Gb space on the hard drive. Although
6 Gb memory is enough, to achieve the memory requirement requires only one CPU to use and the runtime will be >20 hours.
Roughly, the peak of memory usage is the product of number of processors (`8` in the code below) and 6 Gb. 
Each processor will take around 6 Gb memory. If the total memory requirement cannot be met, it is better
to reduce the number of processors to use by change the value passed to the `-p` option when running `REPTILE_compute_score.R`
and `REPTILE_preprocess.py`. 

If the server meets the requirement of 48 Gb memory (6 Gb x 8 CPU to use) and 60 Gb space, you can run the example by simply
doing:
```bash
cd example/
sh run_example.sh > log 2> err
```
Otherwise, you will need to change the value assigned to `num_procs` in `run_example.sh` to reduce the need of memory before running the script. 

Below is the step by step guide of the example.

## Download example data
Create a folder to run the example.
```bash
mkdir REPTILE_example/
cd REPTILE_example/
mkdir -p tmp/ results/
```

Download and unzip the data needed for running the example.
```bash
wget neomorph.salk.edu/yupeng/share/REPTILE_example_data.tar
tar xf REPTILE_example_data.tar
rm REPTILE_example_data.tar
```

## Description of example data
In the newly created `REPTILE_example/` folder, you should see the below files/folders.
* `data/bw/` is folder contains the bigWig files of all epigenetic marks of all nine samples
* `data/data_info_mESC_E11_5.tsv` is the data info file indicating the samples and marks involved as well as the corresponding bigWig files.
* `data/DMR_CG_mESC_E11_5_ext150.bed` contains genomic coordinates of all DMRs.
* `data/mm10_chrLen.tsv` has the length of each chromosome in mouse mm10 reference.
* `data/training_data/` and `data/test_data/` contain the data (files) needed for training and accuracy evaluation respectively.
* `data/training_data/mESC_region_for_train.bed` and `data/test_data/vista_enhancer_mm10_mm.bed` contain the genomic coordinates of regions for training and regions for evaluating REPTILE accuracy.
* `data/training_data/mESC_region_for_train_label.tsv` and `data/test_data/vista_enhancer_state.tsv` are the label files for regions contained in the two files above. 


## Model training
```bash
REPTILE_preprocess.py \
	data/data_info_mESC_E11_5.tsv \
	data/training_data/mESC_region_for_train.bed \
	tmp/training_region \
	-d data/DMR_CG_mESC_E11_5_ext150.bed \
	-p 8
```

```bash
REPTILE_train.R \
	-i data/data_info_mESC_E11_5.tsv \
	-a tmp/training_region.region_with_epimark.tsv \
	-d tmp/training_region.DMR_with_epimark.tsv \
	-l data/training_data/mESC_region_for_train_label.tsv \
	-s mESC \
	-o tmp/REPTILE_model
```

## Predict the enhancer activity of VISTA elements

```bash
REPTILE_preprocess.py \
	data/data_info_mESC_E11_5.tsv \
	data/test_data/vista_enhancer_mm10_mm.bed \
	tmp/test_region \
	-d data/DMR_CG_mESC_E11_5_ext150.bed
```

```bash
## Generating enhancer scores
REPTILE_compute_score.R \
	-i data/data_info_mESC_E11_5.tsv \
	-m tmp/REPTILE_model.reptile \
	-d tmp/test_region.DMR_with_epimark.tsv \
	-a tmp/test_region.region_with_epimark.tsv \
	-s E11_5_HT \
	-o results/E11_5_HT_pred \
	-w
		
## Evaluate the prediction results
echo -n "E11_5_HT "
echo \
`REPTILE_evaluate_prediction.R \
	-p results/E11_5_HT_pred.D.bed \
	-s E11_5_HT \
	-l data/test_data/vista_enhancer_state.tsv`
done
```

`results/results/E11_5_HT_pred.D.bed`


## Generate putative enhancers across genome

```bash
## Generate sliding windows across mouse genome
bedtools makewindows -w 2000 -s 100 -g data/mm10_chrLen.tsv |awk '{print $_"\tbin_"i++}' > tmp/mm10_w2kb_s100bp.bed

## Preprocessing
REPTILE_preprocess.py \
	data/data_info_mESC_E11_5.tsv \
	tmp/mm10_w2kb_s100bp.bed \
	tmp/mm10_w2kb_s100bp \
	-d data/DMR_CG_mESC_E11_5_ext150.bed \
	-n 12 \
	-p 8 \
	-g
```

```bash
## Generating enhancer scores
REPTILE_compute_score.R \
	-i data/data_info_mESC_E11_5.tsv \
	-m tmp/REPTILE_model.reptile \
	-d tmp/mm10_w2kb_s100bp.DMR_with_epimark \
	-a tmp/mm10_w2kb_s100bp.region_with_epimark \
	-s E11_5_HT \
	-o tmp/E11_5_HT_pred \
	-p 8 \
	-n 12
    
## Call putative enhancers based on the scores
REPTILE_call_enhancer.py \
	tmp/E11_5_HT_pred.R.bed \
	-d tmp/E11_5_HT_pred.DMR.bed \
	-p 0.5 \
	-o results/enhancer_E11_5_HT.bed
done
```

`results/enhancer_E11_5_HT.bed`


## More about this example
`E11_5_HT` can be 
