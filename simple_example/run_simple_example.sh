#!/usr/bin/env sh

num_procs=1
num_splits=12

path_to_reptile=../bin/
data_info_file=data/data_info_mESC_E11_5.tsv
dmr_file=data/DMR_CG_mESC_E11_5_ext150.bed

mkdir -p tmp/ results/
rm -rf tmp/* results/* data/

# 0 - Download data
echo -n "Downloading data: "
wget -q neomorph.salk.edu/yupeng/share/REPTILE_simple_example_data.tar
tar xf REPTILE_simple_example_data.tar
rm REPTILE_simple_example_data.tar
echo -e "Done!\n"

# 1 - Training
echo "Training enhancer model"
training_region=data/training_data/mESC_region_for_train.bed
training_label=data/training_data/mESC_region_for_train_label.tsv
## Preprocessing training data
echo -n " - Preprocessing training data: "
${path_to_reptile}/REPTILE_preprocess.py ${data_info_file} ${training_region} tmp/training_region -d ${dmr_file} -p $num_procs
echo "Done!"
## Training
echo -n " - Training: "
${path_to_reptile}/REPTILE_train.R -i ${data_info_file} -d tmp/training_region.DMR_with_epimark.tsv -a tmp/training_region.region_with_epimark.tsv -l ${training_label} -s mESC -o tmp/REPTILE_model
echo -e "Done!\n"

# 2 - Test: Predict enhancer activities of given regions (query regions)
echo "Making prediction on given regions"
test_region=data/test_data/vista_enhancer_mm10_mm.bed
test_label=data/test_data/vista_enhancer_state.tsv
## Preprocessing test data
echo -n " - Preprocessing test data: "
${path_to_reptile}/REPTILE_preprocess.py ${data_info_file} ${test_region} tmp/test_region -d ${dmr_file}
echo "Done!"
## Prediction & Evaluation
echo -n " - Generating enhancer score for:"
for sample in E11_5_FB E11_5_MB E11_5_HB E11_5_HT E11_5_LM E11_5_NT
do
    echo -n " ${sample}"
    ## Generating enhancer scores
    ${path_to_reptile}/REPTILE_compute_score.R -i ${data_info_file} -m tmp/REPTILE_model.reptile -d tmp/test_region.DMR_with_epimark.tsv -a tmp/test_region.region_with_epimark.tsv -s ${sample} -o results/${sample}_pred -w
done
echo -e "\nEvaluation results:"
echo -e "  Sample AUROC AUPR top5 top10 top20"
for sample in E11_5_FB E11_5_MB E11_5_HB E11_5_HT E11_5_LM E11_5_NT
do
    ## Evaluation
    echo -e "  ${sample}" $(${path_to_reptile}/REPTILE_evaluate_prediction.R -p results/${sample}_pred.D.bed -s ${sample} -l ${test_label})
done
echo -e "Prediction done!\n"

# 3 - Generate genome-wide prediction
## Divide genome into sliding windows
echo "Generating genome-wide prediction"
sliding_windows=tmp/mm10_w2kb_s100bp.bed
echo -n " - Making sliding windows: "
bedtools makewindows -w 2000 -s 100 -g data/mm10_chrLen.tsv |awk '{print $_"\tbin_"i++}' > ${sliding_windows}
echo "Done!"
## Preprocessing data for genome-wide prediction
echo -n " - Preprocessing data of DMRs and sliding bins across genome: "
${path_to_reptile}/REPTILE_preprocess.py ${data_info_file} ${sliding_windows} tmp/mm10_w2kb_s100bp -d ${dmr_file} -g -n $num_splits -p $num_procs
echo "Done!"
## Prediction
echo -n " - Generating genome-wide prediction for:"
for sample in E11_5_FB E11_5_MB E11_5_HB E11_5_HT E11_5_LM E11_5_NT
do
    echo -n " ${sample}"
    ## Generating enhancer scores
    ${path_to_reptile}/REPTILE_compute_score.R -i ${data_info_file} -m tmp/REPTILE_model.reptile -d tmp/mm10_w2kb_s100bp.DMR_with_epimark -a tmp/mm10_w2kb_s100bp.region_with_epimark -s ${sample} -o tmp/${sample}_pred -p $num_procs -n $num_splits
    
    ## Call putative enhancers based on the scores
    ${path_to_reptile}/REPTILE_call_enhancer.py -p 0.5 -d tmp/${sample}_pred.DMR.bed tmp/${sample}_pred.R.bed -o results/enhancer_${sample}.bed
done
echo -e "\nGenome-wide prediction done!\n"

