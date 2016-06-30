#!/usr/bin/env Rscript
######################################################################################
## REPTILE
## - Regulatory Element Prediction based on TIssue-specific Local Epigenetic marks
## by Yupeng He (yupeng.he.bioinfo at gmail)
## Salk Institute for Biological Studies
######################################################################################

##  Evaluate performance by comparing prediction and annotations/labels
##  Please email Yupeng He (yupeng.he.bioinfo at gmail) for feedback, question and bug. 


## option parsing
suppressPackageStartupMessages(library("REPTILE",verbose=FALSE))

option_parser = get_option_parser_evaluation()
options <- optparse::parse_args(option_parser)
##No input
if(length(commandArgs(TRUE)) == 0){
    optparse::print_help(option_parser)
    quit("no")
}
##End of option parsing
######################################################################################
## Input from command line options
##Samples in which the activities of annotated regions are used for training
sample_for_prediction = options$sample_for_prediction
##Label file
label_file = options$label_file
##Prediction result file
prediction_result_file = options$prediction_result_file
##Whether warnings should be shown
suppress_warnings = options$suppress_warnings

## Read labels of predicted regions
label_region = read_label(label_file,sample_for_prediction)

## Read prediction result
pred = read.table(prediction_result_file,header=F,stringsAsFactors=F)

## Initialization
AUROC = NA;
AUPR = NA;
ROC = NULL;
PR = NULL;

## Get the intersect of predictions and annotations
## AUROC and AUPC can only be computeted using the intersect
region_id_prediction_with_label = intersect(pred[,4],names(label_region))
if(length(region_id_prediction_with_label) == 0){
    options_printout <- paste0(
        "\n  sample_for_prediction = ",sample_for_prediction,
        "\n  label_file = ",label_file,
        "\n  prediction_result_file = ",prediction_result_file)
    if(!suppress_warnings){
        warning(paste0("No prediction has known label/annotation. Cannot ",
                       "evaluate accuracy.\n",
                       "Inputs are:", options_printout)
                )
    }
}else{
    ## Get the predictions without annotation to know whether they are true
    region_id_prediction_without_label = setdiff(pred[,4],names(label_region))
    if((!suppress_warnings) & length(region_id_prediction_without_label) > 0){
        warning(paste0(length(region_id_prediction_without_label),
                       " predictions have no label/annotation."))
    }
    suppressPackageStartupMessages(require(flux))
    
    ## Get predictions and corresponding annotations
    predictions <- pred[,5]
    names(predictions) = pred[,4]
    if(length(region_id_prediction_without_label) > 0){
        predictions = predictions[region_id_prediction_with_label]
    }
    annotations = label_region[region_id_prediction_with_label]
    
    ## Compute metrics
    num_positives = sum(annotations == 1)
    num_negatives = sum(annotations == 0)
    if(num_positives == 0 | num_negatives == 0){
        if(!suppress_warnings){
            options_printout <- paste0(
                "\n  sample_for_prediction = ",sample_for_prediction,
                "\n  label_file = ",label_file,
                "\n  prediction_result_file = ",prediction_result_file)            
            warning(paste0("Only one class of labels presents. Cannot ",
                           "evaluate accuracy.\n",
                           "Inputs are:", options_printout)
                    )
        }
        cat(paste(AUROC,AUPR,sep="\t"),"\n")
        quit()
    }

    s = sort(predictions,decreasing=T,index.return=T)
    for(i in 1:length(predictions)){
        cutoff = s$x[i]
        num_true_positives = sum(predictions >= cutoff & annotations == 1)
        num_false_positives = sum(predictions >= cutoff & annotations == 0)
        ROC <- cbind(ROC,c(
            ## False Positive Rate
            num_false_positives / num_negatives,
            ## True Positive Rate
            num_true_positives / num_positives
            ))        
        PR <- cbind(PR,c(
            ## Recall (True Positve Rate)
            num_true_positives / num_positives,
            ## Precision
            num_true_positives / sum(predictions >= cutoff)
            ))        
    }
    ROC = t(cbind(c(0,0),ROC,c(1,1)))
    PR = t(cbind(c(0,1),PR[,!is.nan(PR[1,])]))
    AUROC = auc(ROC[,1],ROC[,2])
    AUPR = auc(PR[,1],PR[,2])

}
cat(AUROC,AUPR,
    sum(annotations[s$ix[1:5]] == 1),
    sum(annotations[s$ix[1:10]] == 1),
    sum(annotations[s$ix[1:20]] == 1),
    fill=TRUE,
    sep="\t")

