# REPTILE
### - Regulatory Element Prediction based on TIssue-specific Local Epigenetic marks

[https://github.com/yupenghe/REPTILE]

REPTILE is a tool to identify the precise location of enhancers by integrating histone modification data 
and base-resolution DNA methylation profiles.

Please contact [Yupeng He](mailto:yupeng.he.bioinfo@gmail.com) for for feedbacks, questions or bugs.

## Data information file (DATA_INFO)
## Query region file and DMR file
## Label file
## Other related formats
### bed format
### bigWig file



## Preprocessing
### REPTILE_preprocess.py
```
usage: REPTILE_preprocess.py [-h] [-d DMR_FILE] [-f DMR_OVERLAP_FRACTION]
                             [-e EXTENDING_DMR] [-p NUM_PROCESS]
                             [-n NUM_OUTPUT_FILE] [-g]
                             data_info_file query_region_file output_prefix

Preprocess data into input format for REPTILE. Please email 
Yupeng He (yupeng.he.bioinfo at gmail) for feedback, question and bug.

positional arguments:
  data_info_file        tab-separated file providing information about samples,
                        epigenetic marks and path to corresponding bigwig files.
                        No duplicated mark name is allowed for each sample and
                        no duplicated sample name is allowed for each mark.
                        format:
                          <sample name><\t><mark name><\t><path to bigwig file>
                        example (header is required):
                          sample mark bw_file
                          heart H3K4me1 bigwig/heart_H3K4me1.bw
                          ...
                          brain H3K27ac bigwig/brain_H3K27ac.bw
  query_region_file     bed file of regions used for training and/or prediction,
                        format:
                          <chromosome><\t><start><\t><end>\t<region id>
                        example:
                          chr1 3172000 3173000 reg_0
                          ...
                          chr9 124412000 124413000 reg_302031
  output_prefix         prefix of output files
                        example: heart_enhancer
                          two output files will be:
                          - heart_enhancer.DMR_with_epimark.tsv
                          - heart_enhancer.region_with_epimark.tsv

optional arguments:
  -h, --help            show this help message and exit
  -d DMR_FILE, --DMR-file DMR_FILE
                        bed file specifying differentially methylated regions (DMRs)
                        DMRs are used as high-resolution enhancer candidates to 
                        increase the resolution of training and prediction.
                        File can also contain high-resolution candidate loci from other
                        assays such DNase-seq and ATAC-seq, or loci from motif analysis.
                        format:
                          <chromosome><\t><start><\t><end>\t<region id>
                        example:
                          chr1 3172266 3172488 dmr_0
                          ...
                          chr19 61316546 61316778 dmr_513260
  -f DMR_OVERLAP_FRACTION, --DMR-overlap-fraction DMR_OVERLAP_FRACTION
                        Minimum fraction of DMR in overlapped interval required to
                        consider the DMR to be overlapping with region.
                        default: 1.0
  -e EXTENDING_DMR, --extending-DMR EXTENDING_DMR
                        Based to extend for DMRs to include flanking regions in 
                        calculation of the signal for epigenetic marks. It is ususally
                        helpful because DMRs are often nucleosome free and informative
                        informative epigenetic marks locate in flanking regions.
                        Start of DMRs will be subtracted by EXTENDING_DMR and the end
                        will be increase by EXTENDING_DMR in epigenetic mark signal
                        calculation. EXTENDING_DMR can be negative if you like
                        only consider center regions of DMRs. Be careful about choosing
                        value in that situation. The length of some DMRs may because negative
                        and if so, this program will stop and pump out an error
                        This option does not affect overlap between DMRs and query regions.
                        default: 0
  -p NUM_PROCESS, --num-process NUM_PROCESS
                        Number of threads/processes
  -n NUM_OUTPUT_FILE, --num-output-file NUM_OUTPUT_FILE
                        If greater than 1, split output files equally. 
  -g, --for-genome-wide-prediction
                        This option is used for generating the input files used for genome-wide
                        prediction. Please don't use this option when trying to generate inputs for
                        model training or prediction on given query regions.
                          1. When this option is used, this program will not try to overlap DMRs
                        with query regions. Instead, it will just output the values epigenetic
                        marks of query regions and DMRs.
                          2. When this option is not used, DMRs will be overlapped with query
                        regions and the values of epigenetic marks will be calculated for the DMR
                        in each overlap instance. If a DMR overlaps with multiple query regions,
                        its values of epigenetic marks will be calculated multiple times. In the
                        output file(s), the DMR id will be merged with the id of overlapping query
                        region to specify the overlapping relationship. For example:
                          dmr_23:query_reg_33
                          dmr_23:query_reg_34
                          dmr_23:query_reg_35
                          dmr_24:query_reg_36
                          dmr_24:query_reg_37
                        The overlapping information is required for training. The DMRs with no
                        overlapping query regions will be ignored.

Add option --help or -h to get help information
```


## Train enhancer model based on known enhancers
### REPTILE_train.R
```
Usage: ./REPTILE_train.R \
		 -i data_info_file \
		 -a query_region_epimark_file \
		 -d DMR_epimark_file \
		 -l label_file \
		 -s sample_for_training \
		 -o output_prefix

	Training enhancer model from annotated regions (known enhancers and known negative
	regions). "REPTILE_preprocess.py" can be used to prepare the input files.
	Please email Yupeng He (yupeng.he.bioinfo at gmail) for feedbacks, questions or bugs.

Options:
	-i DATA_INFO_FILE, --data-info-file=DATA_INFO_FILE
		Tab-separated file providing information about samples,
		epigenetic marks and paths to corresponding bigwig files.
		No duplicated mark name is allowed for each sample and
		no duplicated sample name is allowed for each mark. "samples"
		can be different cell/tissue types or different conditions.
		format:
		  <sample name><\t><mark name><\t><path to bigwig file>
		example (header is required):
		  sample mark bw_file
		  heart H3K4me1 bigwig/heart_H3K4me1.bw
		  ...
		  brain H3K27ac bigwig/brain_H3K27ac.bw
	    **Please use the same data info file as the one used in
	      data preprocessing (using "REPTILE_preprocess.py").


	-a ANNOTATED_REGION_EPIMARK_FILE, --annotated-region-epimark-file=ANNOTATED_REGION_EPIMARK_FILE
		Tab-separated file of epigenetic profiles of annotated regions
		format:
		  First line is a header indicating the content of each
		  column. First four columns are chromosome, start, end and
		  region id. Following columns are the scores/enrichment of
		  epigenetic marks for annotated regions.
		example:
		  chr start end id H3K4me1_H1 H3K4me1_H9 H3K4me1_IMR90 ...
		  chr1 3172000 3173000 reg_0 1.43 1.50 0.03 ...
		  ...
		  chr9 124412000 124413000 reg_302031 0.34 0.44 2.42 ...
	    **Highly recommend to to use "REPTILE_preprocess.py" script to generate
	      this file. Run "./REPTILE_preprocess.py -h" for more information.
	      Make sure the same data info file is used.


	-d DMR_EPIMARK_FILE, --DMR-epimark-file=DMR_EPIMARK_FILE
		Tab-separated file of epigenetic profiles of differentially
		methyated regions (DMRs). Format is same as that of 
		annotated_region_epimark_file. DMRs are used as high-resolution
		enhancer candidates to increase the resolution of training and
		prediction. The candidate loci can also come from assays like
		DNase-seq or ATAC-seq or analysis on motifs.
		format:
		  First line is a header indicating the content of each
		  column. First four columns are chromosome, start, end and
		  region id. Following columns are the scores/enrichment of
		  epigenetic marks for annotated regions.
		example:
		  chr start end id H3K4me1_H1 H3K4me1_H9 H3K4me1_IMR90 ...
		  chr1 3172266 3172488 dmr_0 1.43 1.50 0.03 ...
		  ...
		  chr19 61316546 61316778 dmr_513260 0.34 0.44 2.42 ...
	    **Highly recommend to to use "REPTILE_preprocess.py" script to generate
	      this file. Run "./REPTILE_preprocess.py -h" for more information.
	      Make sure the same data info file is used.


	-l LABEL_FILE, --label-file=LABEL_FILE
		Tab-separated file labeling what annotated regions are active in
		which sample(s).
		format:
		  The file has multiple columns. First column is region id.
		  Each following column corresponds to one sample and the values
		  indicate whether the region is active in the sample:
		    - 1: active,
		    - 0: no activity
		    - NA: unknown
		  First line is a header indicating the content of each column.
		  Name of first column is "id" and others are sample names.
		example:
		  id H1 H9 IMR90 ...
		  reg_0 1 1 0 ...
		  ...
		  reg_302031 0 0 1 ....


	-s SAMPLES_FOR_TRAINING, --samples-for-training=SAMPLES_FOR_TRAINING
		Samples in which the activities of annotated regions are used for
		training model.
		format:
		  Sample names separated by comma.
		example:
		  H1,H9,IMR90


	-r REF_SAMPLE, --reference-samples=REF_SAMPLE
		Samples used as reference to calculate intensity deviation
		format:
		  Sample names separated by comma.
		example:
		  E11_5_FB,E11_5_HT,E11_5_MB


	-o OUTPUT_PREFIX, --output-prefix=OUTPUT_PREFIX
		Prefix of the output file storing the model obtained from
		training. The output file is "<OUTPUT_PREFIX>.reptile".
		example: enhancer_model
		  The output file is enhancer_model.reptile


	-c CLASSIFIER_FAMILY, --classifier-family=CLASSIFIER_FAMILY
		Classifier family to use in the prediction model
		  default: RandomForest
		Classifiers available:
		 - RandomForest: random forest
		 - Logistic: logistic regression
		 - SVM: support vector machine 
		 - NaiveBayes: naive bayes model


	-x, --no-intensity-deviation
		If this option is used, REPTILE will not compute the intensity
		deviation feature, which captures the tissue-specificity of
		enhancers.


	-t NUM_TREES, --number-of-trees=NUM_TREES
		Number of trees to be constructed in random forest
		  classifier. Ignored when other classifiers are
		  used.
		  default: 2000


	-h, --help
		Show this help message and exit
```


## Calculating enhancer scores
### REPTILE_compute_score.R
```
Usage: ./REPTILE_compute_score.R \
		 -i data_info_file \
		 -m model_file \
		 -a query_region_epimark_file \
		 -d DMR_epimark_file \
		 -s sample_for_prediction \
		 -o output_prefix

	Predicting enhancer activities of query regions in target sample. This script will
	calculate enhancer confidence score for each query region and each DMR.
	Please email Yupeng He (yupeng.he.bioinfo at gmail) for feedbacks, questions or bugs.

Options:
	-i DATA_INFO_FILE, --data-info-file=DATA_INFO_FILE
		Tab-separated file providing information about samples,
		epigenetic marks and paths to corresponding bigwig files.
		No duplicated mark name is allowed for each sample and
		no duplicated sample name is allowed for each mark. "samples"
		can be different cell/tissue types or different conditions.
		format:
		  <sample name><\t><mark name><\t><path to bigwig file>
		example (header is required):
		  sample mark bw_file
		  heart H3K4me1 bigwig/heart_H3K4me1.bw
		  ...
		  brain H3K27ac bigwig/brain_H3K27ac.bw
	    **Please use the same data info file as the one used in
	      data preprocessing (generation of input file).


	-m MODEL_FILE, --model-file=MODEL_FILE
		Enhancer model learned from known enhancers and known negative
		regions. It is the output file from "REPTILE_train.R".
		example:
		  enhancer_model.reptile


	-s SAMPLE_FOR_PREDICTION, --sample-for-prediction=SAMPLE_FOR_PREDICTION
		Sample in which the activities of query regions are to be predicted.
		format:
		  Sampl name
		example:
		  E11_5_FB


	-r REF_SAMPLE, --reference-samples=REF_SAMPLE
		Samples used as reference to calculate intensity deviation
		format:
		  Sample names separated by comma.
		example:
		  E11_5_FB,E11_5_HT,E11_5_MB


	-a QUERY_REGION_EPIMARK_FILE, --query-region-epimark-file=QUERY_REGION_EPIMARK_FILE
		Tab-separated file of epigenetic profiles of regions to be predicted.
		If "-n" option is used and the number is greater than 1, this option
		should be set as the prefix of query region files. For example, 
		if this option is "mm10_genome.region_epimark" and "-n" is 8,
		REPTILE expects input files:
		  mm10_genome.region_epimark.0.tsv
		  mm10_genome.region_epimark.1.tsv
		  ...
		  mm10_genome.region_epimark.6.tsv
		  mm10_genome.region_epimark.7.tsv
		The format of the file(s) is:
		  First line is a header indicating the content of each
		  column. First four columns are chromosome, start, end and
		  region id. Following columns are the scores/enrichment of
		  epigenetic marks for annotated regions.
		example:
		  chr start end id H3K4me1_H1 H3K4me1_H9 H3K4me1_IMR90 ...
		  chr1 3172000 3173000 reg_0 1.43 1.50 0.03 ...
		  ...
		  chr9 124412000 124413000 reg_302031 0.34 0.44 2.42 ...
	    **Highly recommend to to use "REPTILE_preprocess.py" script to generate
	      this file. Run "./REPTILE_preprocess.py -h" for more information.
	      Make sure the same data info file is used.


	-d DMR_EPIMARK_FILE, --DMR-epimark-file=DMR_EPIMARK_FILE
		Tab-separated file of epigenetic profiles of differentially
		methyated regions (DMRs). Format is same as that of 
		annotated_region_epimark_file. DMRs are used as high-resolution
		enhancer candidates to increase the resolution of training and
		prediction. The candidate loci can also come from assays like
		DNase-seq or ATAC-seq or analysis on motifs.
		Same as "-a": If "-n" option is used and the number is
		greater than 1, this option should be set as the prefix of files.
		For example, if this option is "mm10_genome.DMR_epimark" and "-n" is 8
		REPTILE expects input files:
		  mm10_genome.DMR_epimark.0.tsv
		  mm10_genome.DMR_epimark.1.tsv
		  ...
		  mm10_genome.DMR_epimark.6.tsv
		  mm10_genome.DMR_epimark.7.tsv
		The format of the file(s) is:
		  First line is a header indicating the content of each
		  column. First four columns are chromosome, start, end and
		  region id. Following columns are the scores/enrichment of
		  epigenetic marks for annotated regions.
		example:
		  chr start end id H3K4me1_H1 H3K4me1_H9 H3K4me1_IMR90 ...
		  chr1 3172266 3172488 dmr_0 1.43 1.50 0.03 ...
		  ...
		  chr19 61316546 61316778 dmr_513260 0.34 0.44 2.42 ...
	    **Highly recommend to to use "REPTILE_preprocess.py" script to generate
	      this file. Run "./REPTILE_preprocess.py -h" for more information.
	      Make sure the same data info file is used.


	-o OUTPUT_PREFIX, --output-prefix=OUTPUT_PREFIX
		Prefix of the output file storing the predictions.
		The output file is "<OUTPUT_PREFIX>.DMR.bed"
		and "<OUTPUT_PREFIX>.R.bed" corresponding to the predictions
		from the classifier for DMRs on DMRs and classifier for query
		regions on query regions, respectively. If "-w" option is used,
		REPTILE will also generate "<OUTPUT_PREFIX>.D.bed" file as the
		combined prediction. This file will store the enhancer score of
		each query region,which is defined as the maximum of scores of
		the whole region and the overlapping DMRs
		example:
		  With "--output-prefix E11_5_FB_enhancer", two output files
		  are:
		  - E11_5_FB_enhancer.DMR.bed
		  - E11_5_FB_enhancer.R.bed
		  - E11_5_FB_enhancer.D.bed (if "-w" is used)
	    **output format:
	      BED format with 5th column as the enhancer confidence score.


	-c CLASSIFIER_FAMILY, --classifier-family=CLASSIFIER_FAMILY
		Classifier family to use in prediction model. "Logistic" is
		faster but slightly less accurate compared to "RandomForest".
		It is useful in generating quick and good results.
		  default: RandomForest
		Classifiers available:
		 - RandomForest: random forest
		 - Logistic: logistic regression
		 - SVM: support vector machine 
		 - NaiveBayes: naive bayes model


	-x, --no-intensity-deviation
		If this option is used, REPTILE will not compute the intensity
		deviation feature, which captures the tissue-specificity of
		enhancers.


	-w, --not-genome-wide-prediction
		This option is used to generate enhancer scores for given regions
		by combining the results from two classifiers (one for DMRs and.
		one for query regions). If this option is used, REPTILE will
		generate one additional file, "<OUTPUT_PREFIX>.D.bed", to show
		the combined scores of each query region. The combined score is
		defined as the maximum of scores of the whole region and the
		overlapping DMRs.
		  To use this option, the input files should be generated without
		the "-g" option during preprocessing step ("REPTILE_preprocess.py").
		That is, the ids in DMR_epimark_file ("-d") should contain the ids of
		overlapping query regions. For example, "dmr_304:reg_1750".


	-p NUM_PROCS, --number-of-processors=NUM_PROCS
		Number of processors to use. This option requires two R
		packages: "foreach" and "doMC".


	-n NUM_SPLITS, --number-of-splits=NUM_SPLITS
		REPTILE accepts splitted input files and this option
		specifies the the number of DMR epimark files (or
		query region epimark files after the split). For
		example, if "-a" is "mm10_genome.region_epimark" and
		"-d" is "mm10_genome.DMR_epimark and this option is 8,
		REPTILE will search for below input files:
		  mm10_genome.DMR_epimark.0.tsv
		  ...
		  mm10_genome.DMR_epimark.7.tsv
		and
		  mm10_genome.region_epimark.0.tsv
		  ...
		  mm10_genome.region_epimark.7.tsv

		If this option is 1, there is no split and REPTILE will
		assume "-a" and "-d" options give the exact names of
		of input files. Thus, it will search for two files:
		  mm10_genome.DMR_epimark
		and
		  mm10_genome.region_epimark


	-h, --help
		Show this help message and exit
```


## Generating enhancer calls in the genome
### REPTILE_call_enhancer.py
```
usage: REPTILE_call_enhancer.py [-h] [-d DMR_FILE] [-o OUTPUT_FILE]
                                [-p CUTOFF] [-s PEAK_SPACING]
                                query_region_file

Used to call enhancers based on REPTILE enhancer confidence scores 
of DMRs and sliding windows. Please email Yupeng He (yupeng.he.bioinfo at gmail) for 
feedbacks, questions or bugs.

positional arguments:
  query_region_file     bed file of query regions with enhancer scores such as the
                        *.R.bed files from REPTILE_compute_score.R.
                        format:
                          <chromosome><\t><start><\t><end>\t<region id>
                        example:
                          chr1 3172000 3173000 reg_0
                          ...
                          chr9 124412000 124413000 reg_302031

optional arguments:
  -h, --help            show this help message and exit
  -d DMR_FILE, --DMR-file DMR_FILE
                        bed file of DMRs with enhancer scores such as the
                        *.DMR.bed file from REPTILE_compute_score.R.
                        format:
                          <chromosome><\t><start><\t><end>\t<DMR id>
                        example:
                          chr1 5143000 5143500 dmr_0
                          ...
                          chr19 4412008 4412403 dmr_2052
  -o OUTPUT_FILE, --output-file OUTPUT_FILE
                        Name of output BED file of predicted enhancers and their scores.
  -p CUTOFF, --cutoff CUTOFF
                        Minimum enhancer confidence score for a sliding window
                        to be consider as a peak
                        default: 0.5
  -s PEAK_SPACING, --peak-spacing PEAK_SPACING
                        Minimum distance for two sliding windows
                        to be called as separate peaks

Add option --help or -h to get help information
```


## Evaluation the prediction results
### REPTILE_evaluate_prediction.R
```
Usage: ./REPTILE_evaluate_prediction.R \
		 -s sample_for_prediction \
		 -l label_file \
		 -p query_region_file_with_score

	Evaluating the prediction accuracy by calculating the Area Under Receiver Operating.
	Characteristic (AUROC) and the Area Under Precision-Recall curve (AUPR). It will also
	calculate the percent of true positives in the top 5, 10 and 20 predictions
	Please email Yupeng He (yupeng.he.bioinfo at gmail) for feedbacks, questions or bugs.

Options:
	-s SAMPLE_FOR_PREDICTION, --sample-for-prediction=SAMPLE_FOR_PREDICTION
		Sample where the enhancer confidence scores are generated for query
		regions.
		format:
		  Sampl name
		example:
		  E11_5_FB


	-l LABEL_FILE, --label-file=LABEL_FILE
		Tab-separated file labeling what query (annotated) regions are active
		active enhancers in which sample(s).
		format:
		  The file has multiple columns. First column is region id.
		  Each following column corresponds to one sample and the values
		  indicate whether the region is active in the sample:
		    - 1: active,
		    - 0: no activity
		    - NA: unknown
		  First line is a header indicating the content of each column.
		  Name of first column is "id" and others are sample names.
		example:
		  id E11_5_FB E11_5_MB E11_5_HB ...
		  reg_0 1 1 0 ...
		  ...
		  reg_302031 0 0 1 ....


	-p PREDICTION_RESULT_FILE, --prediction-result-file=PREDICTION_RESULT_FILE
		bed file of query regions with enhancer scores (5th column) such as the
		"*.R.bed" files and "*.D.bed" files from REPTILE_compute_score.R.
		format:
		  BED file with five columns, chromosome, start, end, region id
		  and score of enhancer activity.
		example:
		  chr1 3172000 3173000 reg_0 0.7 ...
		  ...
		  chr9 124412000 124413000 reg_302031 0.1 ...
	    **Output files from "REPTILE_compute_score.R" are in this format
	      and can be directly used as input for this script


	-q, --suppress-warnings
		Showing warnings is default. If this option is enabled, warnings will
		not be shown.

	-h, --help
		Show this help message and exit
```

