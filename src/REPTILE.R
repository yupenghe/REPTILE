######################################################################################
## REPTILE
## - Regulatory Element Prediction based on TIssue-specific Local Epigenetic marks
## by Yupeng He (yupeng.he.bioinfo at gmail)
## Salk Institute for Biological Studies
## Jun. 6, 2016
######################################################################################

##  Functions required for other REPTILE R scripts
##  Please email Yupeng He (yupeng.he.bioinfo at gmail) for feedbacks, questions or bugs. 

read_epigenomic_data <- function(data_info,
                                 epimark_file,
                                 query_sample,
                                 ref_sample = NULL,
                                 incl_dev = T){
    ##
    ## Function to read epimark file and generate featuers for classifier
    ##

    ##Check if the query sample is in input data
    if(!(query_sample %in% unique(data_info$sample))){
        avail_samples <- paste0(paste0("    ",
                                       unique(data_info$sample)
                                       ),
                                sep="\n")
        stop(paste0(query_sample,
                    " is not in data_info file!\n",                   
                    "  Available samples are:\n"),
             avail_samples)        
    }
    ##Check if the reference sample(s) is in input data
    if(!is.null(ref_sample) &
       sum(ref_sample %in% unique(data_info$sample)) != length(ref_sample)
       ){
        avail_samples <- paste0(paste0("    ",
                                       unique(data_info$sample)
                                       ),
                                sep="\n")
        stop(paste0(paste0(ref_sample[!(ref_sample %in% unique(data_info$sample))],
                           collapse=","),
                    " is not in data_info file!\n",                   
                    "  Available samples are:\n"),
             avail_samples)        
    }
    
    ##Read epimark file
    epimark <- read.table(epimark_file,sep = "\t",header=T,stringsAsFactors=F)
    if(nrow(epimark) == 0){
        return(NULL)
    }       
    rownames(epimark) = epimark$id
    epimark = data.matrix(epimark[,5:ncol(epimark)])

    ##Check if data_info is compatible with header of epimark_file
    if( (nrow(data_info) != ncol(epimark)) ){
        stop(paste0("Sample information and the header of \"",
                    epimark_file,"\" are not compatible!\n",
                    "  Please check whether the same data_info file was used for generating \"",
                    epimark_file,"\""))        
    }
    else if(sum(paste0(data_info$mark,"_",data_info$sample) != colnames(epimark))>0){
        stop(paste0("Sample information and the header of \"",
                    epimark_file,"\" are not compatible!\n",
                    "  Please check whether the same data_info file was used for generating \"",
                    epimark_file,"\""))        
    }
    
    ##Calculate deviation as features
    if(incl_dev){
        epimark_dev <- calculate_epimark_deviation(data_info,
                                                   epimark,
                                                   query_sample,
                                                   ref_sample
                                                   )
    }
    
    ##Only include epigenetic marks in query sample
    epimark = epimark[,data_info$sample==query_sample]
    colnames(epimark) = data_info$mark[data_info$sample==query_sample]
    
    ##Add deviations as features
    if(incl_dev){
        epimark = cbind(epimark,epimark_dev)
    }
    
    return(epimark)
}

calculate_epimark_deviation <- function(data_info,x,query_sample,ref_sample=NULL){
    ##
    ## Function to calculate deviation of each epigenetic mark in
    ## query sample compared to the average of other samples
    ##
    
    res = NULL
    dev_names = NULL
    if(is.null(ref_sample)){
        for(query_mark in data_info$mark[data_info$sample==query_sample]){
            ind = (data_info$mark == query_mark)
            if(sum(ind) == 1){
                next
            }else if(sum(ind) == 2){
                d = x[,ind & (data_info$sample == query_sample)]
                d = d - x[,ind & (data_info$sample != query_sample)]
            }else{            
                ## calculation
                d = x[,ind & (data_info$sample == query_sample)]
                d = d - rowMeans(x[,ind & (data_info$sample != query_sample)])
            }
            dev_names = c(dev_names,paste0(query_mark,"_dev"))
            res = cbind(res,d)
        }
    }else{
        for(query_mark in data_info$mark[data_info$sample==query_sample]){
            ind = (data_info$mark == query_mark)
            if(sum(ind & data_info$sample==query_sample) == 0){
                next
            }
            
            if(sum(ind) == 1){
                next
            }else if(sum(ind & (data_info$sample %in% ref_sample)) == 1){
                d = x[,ind & (data_info$sample == query_sample)]
                d = d - x[,ind & (data_info$sample %in% ref_sample)]
            }else{            
                ## calculation
                d = x[,ind & (data_info$sample == query_sample)]
                d = d - rowMeans(x[,ind & (data_info$sample %in% ref_sample)])
            }
            dev_names = c(dev_names,paste0(query_mark,"_dev"))
            res = cbind(res,d)
        }
    }
    colnames(res) = dev_names    
    return(res)
}

read_label <- function(label_file,query_sample){
    ##
    ## Function to read the labels (enhancer activities) of
    ## query regions
    ##

    activity <- read.table(label_file,sep = "\t",
                           header=T,
                           stringsAsFactors=F)    
    if(sum(query_sample %in% colnames(activity)[-1]) != length(query_sample)){
        avail_samples <- paste0(paste0("    ",
                                       colnames(activity[,-1])
                                       ),
                                sep="\n")
        stop(paste0(paste0(query_sample,collapse=","),
                    " has no data in label file!\n",                   
                    "  Label file contains data of:\n"),
             avail_samples)        
    }
    res = activity[,query_sample]
    names(res) = activity[,1]
    return(res)
}

reptile_train_one_mode <- function(epimark,label,
                                   family,
                                   ntree,nodesize){
    ##
    ## Function to train one mode of REPTILE
    ##
    if(family == "Logistic"){## Logistric regression
        df <- data.frame(label=as.numeric(as.character(label)),
                         epimark)
        mod <- glm(label~.,data=df,family=binomial)
    }
    else{## Randome forest
        suppressPackageStartupMessages(requireNamespace("randomForest",
                                                        quietly = TRUE))
        mod <- randomForest::randomForest(epimark,label,
                                          nodesize=nodesize,
                                          ntree=ntree,importance=T)
    }
    return(mod)
}

reptile_train <- function(epimark_region,label_region,
                          epimark_DMR=NULL,label_DMR=NULL,
                          family="randomForest",
                          ntree=2000,nodesize=1){
    ##
    ## Function to train both modes of REPTILE
    ## - D (DMR-based mode)
    ## - R (Region-based mode)
    ##

    reptile = list("D"=NULL,"R"=NULL)
    ## R (Region-based)
    reptile$R <- reptile_train_one_mode(epimark_region,label_region,
                                        family,
                                        nodesize=nodesize,
                                        ntree=ntree)
    ## D (DMR-based)
    if(!is.null(epimark_DMR) & !is.null(label_DMR)){
        reptile$D <- reptile_train_one_mode(epimark_DMR,label_DMR,
                                            family,
                                            nodesize=nodesize,
                                            ntree=ntree)
    }
    return(reptile)
}

reptile_predict_one_mode <- function(reptile_classifier,
                                     epimark,
                                     family){
    ##
    ## Function to make prediction using one mode of REPTILE
    ##

    if(family == "Logistic"){
        ##Logistric regression
        pred = predict(reptile_classifier,data.frame(epimark))
        names(pred) = rownames(epimark)
    }
    else{
        suppressPackageStartupMessages(requireNamespace("randomForest",
                                                        quietly=T))
        ##Random Forest
        pred = predict(reptile_classifier,epimark,type="prob")[,2]
    }
    return(pred)

}

reptile_predict <- function(reptile_model,
                            epimark_region,
                            epimark_DMR=NULL,
                            family="randomForest"
                            ){
    ##
    ## Function to make prediction using both modes of REPTILE
    ## - D (DMR-based mode)
    ## - R (Region-based mode)
    ##

    ##Input format
    ##rownames of epimark_region follow the format
    ## - reg_id
    ##rownames of epimark_DMR follow the format
    ## - dmr_id:reg_id

    
    pred = list("D"=NULL,"R"=NULL,"DMR"=NULL)
    ## R
    pred$R = reptile_predict_one_mode(reptile_model$R,epimark_region,family)

    ## D
    if(!is.null(epimark_DMR)){
        pred$D = reptile_predict_one_mode(reptile_model$D,epimark_DMR,family)
        pred$DMR = pred$D

        ## For regions containing no DMRs, use prediction from R
        ind = ! ((names(pred$R) %in% sapply(strsplit(names(pred$D),":"),
                                            function(x) return(x[2]))))
        pred$D = c(pred$D,pred$R[ind])
        ##pred$D = c(pred$D,pred$R)
        names(pred$D) <- sapply(strsplit(names(pred$D),":"),
                                function(x) {
                                    if(length(x) == 1){
                                        return(x)}
                                    else{
                                        return(x[2])}
                                }
                                )
        ## aggregate
        pred$D = aggregate(pred$D,list(names(pred$D)),max)
        ## Reformat
        region_id = pred$D[,1]
        pred$D = pred$D[,2]
        names(pred$D) = region_id
    }else{
        pred$D = pred$R
    }
    return(pred)
}

reptile_predict_genome_wide <- function(reptile_model,
                                        epimark_region,
                                        epimark_DMR=NULL,
                                        family="randomForest"
                                        ){
    ##
    ## Function to make prediction using both modes of REPTILE:
    ## - D (DMR-based mode)
    ## - R (Region-based mode)
    ##

    ##Input format
    ##rownames of epimark_region follow the format
    ## - reg_id
    ##rownames of epimark_DMR follow the format
    ## - dmr_id
        
    pred = list("R"=NULL,"DMR"=NULL)
    ## R
    pred$R = reptile_predict_one_mode(reptile_model$R,epimark_region,family)

    ## D
    if(!is.null(epimark_DMR)){
        pred$DMR = reptile_predict_one_mode(reptile_model$D,epimark_DMR,family)
    }
    return(pred)
}

reptile_eval_prediction <- function(predictions, annotations){
    ##
    ## Function to evaluate the results of prediction
    ## - AUROC and AUPR will be calculated
    ##
    annotations = annotations[names(predictions)]
    num_positives = sum(annotations == 1)
    num_negatives = sum(annotations == 0)
    curve_ROC = NULL # ROC curve
    curve_PR = NULL # Precision-Recall curve
    s = sort(predictions,decreasing=T,index.return=T)
    for(i in 1:length(predictions)){
        cutoff = s$x[i]
        num_true_positives = sum(predictions >= cutoff & annotations == 1)
        num_false_positives = sum(predictions >= cutoff & annotations == 0)
        FPR = num_false_positives / num_negatives # False Positive Rate
        TPR = num_true_positives / num_positives # True Positive Rate
        PR = num_true_positives / sum(predictions >= cutoff) # Precision
        curve_ROC <- cbind(curve_ROC,c(FPR,TPR))
        curve_PR <- cbind(curve_PR,c(TPR,PR))
    }
    curve_ROC = t(cbind(c(0,0),curve_ROC,c(1,1)))
    curve_PR = t(cbind(c(0,1),curve_PR[,!is.nan(curve_PR[1,])]))
    suppressPackageStartupMessages(requireNamespace("flux",quietly = TRUE))
    AUROC = flux::auc(curve_ROC[,1],curve_ROC[,2])
    AUPR = flux::auc(curve_PR[,1],curve_PR[,2])
    return(list("AUROC"=AUROC, "AUPR"=AUPR))
}


get_option_parser_training <- function(){
    ##
    ## Function to generate option parser for training script
    ##
    suppressPackageStartupMessages(requireNamespace("optparse",
                                                    quietly = TRUE))
    option_list <- list(
        optparse::make_option(c("-i", "--data-info-file"),
                              type="character",
                              metavar="data_info_file",
                              dest="data_info_file",
                              help = paste0(
                                  "Tab-separated file providing information about samples,\n",
                                  "\t\tepigenetic marks and paths to corresponding bigwig files.\n",
                                  "\t\tNo duplicated mark name is allowed for each sample and\n",
                                  "\t\tno duplicated sample name is allowed for each mark. \"samples\"\n",
                                  "\t\tcan be different cell/tissue types or different conditions.\n",
                                  "\t\tformat:\n",
                                  "\t\t  <sample name><\\t><mark name><\\t><path to bigwig file>\n",
                                  "\t\texample (header is required):\n",
                                  "\t\t  sample mark bw_file\n",
                                  "\t\t  heart H3K4me1 bigwig/heart_H3K4me1.bw\n",
                                  "\t\t  ...\n",
                                  "\t\t  brain H3K27ac bigwig/brain_H3K27ac.bw\n",
                                  "\t    **Please use the same data info file as the one used in\n",
                                  "\t      data preprocessing (using \"REPTILE_preprocess.py\").\n")

                              ),
        
        optparse::make_option(c("-a", "--annotated-region-epimark-file"),
                              type="character",
                              metavar="annotated_region_epimark_file",
                              dest="annotated_region_epimark_file",
                              help=paste0(
                                  "Tab-separated file of epigenetic profiles of annotated regions\n",
                                  "\t\tformat:\n",
                                  "\t\t  First line is a header indicating the content of each\n",
                                  "\t\t  column. First four columns are chromosome, start, end and\n",
                                  "\t\t  region id. Following columns are the scores/enrichment of\n",
                                  "\t\t  epigenetic marks for annotated regions.\n",
                                  "\t\texample:\n",
                                  "\t\t  chr start end id H3K4me1_H1 H3K4me1_H9 H3K4me1_IMR90 ...\n",
                                  "\t\t  chr1 3172000 3173000 reg_0 1.43 1.50 0.03 ...\n",
                                  "\t\t  ...\n",
                                  "\t\t  chr9 124412000 124413000 reg_302031 0.34 0.44 2.42 ...\n",
                                  "\t    **Highly recommend to to use \"REPTILE_preprocess.py\" script to generate\n",
                                  "\t      this file. Run \"./REPTILE_preprocess.py -h\" for more information.\n",
                                  "\t      Make sure the same data info file is used.\n")
                              ),
        
        optparse::make_option(c("-d", "--DMR-epimark-file"),
                              type="character",
                              metavar="DMR_epimark_file",
                              dest="DMR_epimark_file",
                              help=paste0(
                                  "Tab-separated file of epigenetic profiles of differentially\n",
                                  "\t\tmethyated regions (DMRs). Format is same as that of \n",
                                  "\t\tannotated_region_epimark_file. DMRs are used as high-resolution\n",
                                  "\t\tenhancer candidates to increase the resolution of training and\n",
                                  "\t\tprediction. The candidate loci can also come from assays like\n",
                                  "\t\tDNase-seq or ATAC-seq or analysis on motifs.\n",
                                  "\t\tformat:\n",
                                  "\t\t  First line is a header indicating the content of each\n",
                                  "\t\t  column. First four columns are chromosome, start, end and\n",
                                  "\t\t  region id. Following columns are the scores/enrichment of\n",
                                  "\t\t  epigenetic marks for annotated regions.\n",
                                  "\t\texample:\n",
                                  "\t\t  chr start end id H3K4me1_H1 H3K4me1_H9 H3K4me1_IMR90 ...\n",
                                  "\t\t  chr1 3172266 3172488 dmr_0 1.43 1.50 0.03 ...\n",
                                  "\t\t  ...\n",
                                  "\t\t  chr19 61316546 61316778 dmr_513260 0.34 0.44 2.42 ...\n",
                                  "\t    **Highly recommend to to use \"REPTILE_preprocess.py\" script to generate\n",
                                  "\t      this file. Run \"./REPTILE_preprocess.py -h\" for more information.\n",
                                  "\t      Make sure the same data info file is used.\n")
                              ),
        
        optparse::make_option(c("-l", "--label-file"),
                              type="character",
                              metavar="label_file",
                              dest="label_file",
                              help=paste0(
                                  "Tab-separated file labeling what annotated regions are active in\n",
                                  "\t\twhich sample(s).\n",
                                  "\t\tformat:\n",
                                  "\t\t  The file has multiple columns. First column is region id.\n",
                                  "\t\t  Each following column corresponds to one sample and the values\n",
                                  "\t\t  indicate whether the region is active in the sample:\n",
                                  "\t\t    - 1: active,\n",
                                  "\t\t    - 0: no activity\n",
                                  "\t\t    - NA: unknown\n",
                                  "\t\t  First line is a header indicating the content of each column.\n",
                                  "\t\t  Name of first column is \"id\" and others are sample names.\n",
                                  "\t\texample:\n",
                                  "\t\t  id H1 H9 IMR90 ...\n",
                                  "\t\t  reg_0 1 1 0 ...\n",
                                  "\t\t  ...\n",
                                  "\t\t  reg_302031 0 0 1 ....\n")
                              ),
        
        optparse::make_option(c("-s", "--samples-for-training"),
                              type="character",
                              metavar="samples_for_training",
                              dest="samples_for_training",
                              help=paste0(
                                  "Samples in which the activities of annotated regions are used for\n",
                                  "\t\ttraining model.\n",
                                  "\t\tformat:\n",
                                  "\t\t  Sample names separated by comma.\n",
                                  "\t\texample:\n",
                                  "\t\t  H1,H9,IMR90\n")
                              ),

        optparse::make_option(c("-r", "--reference-samples"),
                              type="character",
                              metavar="ref_sample",
                              dest="ref_sample",
                              default=NULL,
                              help=paste0(
                                  "Samples used as reference to calculate intensity deviation\n",
                                  "\t\tformat:\n",
                                  "\t\t  Sample names separated by comma.\n",
                                  "\t\texample:\n",
                                  "\t\t  E11_5_FB,E11_5_HT,E11_5_MB\n")
                              ),

        optparse::make_option(c("-o", "--output-prefix"),
                              type="character",
                              metavar="output_prefix",
                              dest="output_prefix",
                              help=paste0(
                                  "Prefix of the output file storing the model obtained from\n",
                                  "\t\ttraining. The output file is \"<OUTPUT_PREFIX>.reptile\".\n",
                                  "\t\texample: enhancer_model\n",
                                  "\t\t  The output file is enhancer_model.reptile\n")
                              ),
        
        optparse::make_option(c("-c", "--classifier-family"),
                              type="character",
                              default="RandomForest",
                              metavar="classifier_family",
                              dest="classifier_family",
                              help=paste0(
                                  "Classifier family to use in the prediction model\n",
                                  "\t\t  default: RandomForest\n",
                                  "\t\tClassifiers available:\n",
                                  "\t\t - RandomForest: random forest\n",
                                  "\t\t - Logistic: logistic regression\n"
                                  )
                              ),

        optparse::make_option(c("-x", "--no-intensity-deviation"),
                              type='logical',
                              action = "store_false",
                              default=TRUE,
                              metavar="incl_dev",
                              dest="incl_dev",
                              help=paste0(
                                  "If this option is used, REPTILE will not compute the intensity\n",
                                  "\t\tdeviation feature, which captures the tissue-specificity of\n",
                                  "\t\tenhancers.\n")
                              ),

        optparse::make_option(c("-t", "--number-of-trees"),
                              type="double",
                              default=2000,
                              metavar="num_trees",
                              dest="num_trees",
                              help=paste0(
                                  "Number of trees to be constructed in random forest\n",
                                  "\t\t  classifier. Ignored when other classifiers are\n",
                                  "\t\t  used.\n",
                                  "\t\t  default: 2000\n")
                              )
        )
        
    description <- paste0("\tTraining enhancer model from annotated regions (known enhancers and known negative\n",
                          "\tregions). \"REPTILE_preprocess.py\" can be used to prepare the input files.\n",
                          "\tPlease email Yupeng He (yupeng.he.bioinfo at gmail) for feedbacks, questions or bugs.")

    usage <- paste0("Usage: ./REPTILE_train.R \\\n",
                    "\t\t -i data_info_file \\\n",
                    "\t\t -a query_region_epimark_file \\\n",
                    "\t\t -d DMR_epimark_file \\\n",
                    "\t\t -l label_file \\\n",
                    "\t\t -s sample_for_training \\\n",
                    "\t\t -o output_prefix\n"
                    )

    option_parser <- optparse::OptionParser(usage = usage,
                                            description = description,
                                            option_list=option_list)
    return(option_parser)
}

get_option_parser_compute_score <- function(){
    ##
    ## Function to generate option parser for REPTILE_compute_score.R
    ##

    suppressPackageStartupMessages(requireNamespace("optparse",
                                                    quietly = TRUE))

    option_list <- list(
        optparse::make_option(c("-i", "--data-info-file"),
                              type="character",
                              metavar="data_info_file",
                              dest="data_info_file",
                              help = paste0(
                                  "Tab-separated file providing information about samples,\n",
                                  "\t\tepigenetic marks and paths to corresponding bigwig files.\n",
                                  "\t\tNo duplicated mark name is allowed for each sample and\n",
                                  "\t\tno duplicated sample name is allowed for each mark. \"samples\"\n",
                                  "\t\tcan be different cell/tissue types or different conditions.\n",
                                  "\t\tformat:\n",
                                  "\t\t  <sample name><\\t><mark name><\\t><path to bigwig file>\n",
                                  "\t\texample (header is required):\n",
                                  "\t\t  sample mark bw_file\n",
                                  "\t\t  heart H3K4me1 bigwig/heart_H3K4me1.bw\n",
                                  "\t\t  ...\n",
                                  "\t\t  brain H3K27ac bigwig/brain_H3K27ac.bw\n",
                                  "\t    **Please use the same data info file as the one used in\n",
                                  "\t      data preprocessing (generation of input file).\n")
                              ),
        
        optparse::make_option(c("-m", "--model-file"),
                              type="character",
                              metavar="model_file",
                              dest="model_file",
                              help=paste0(
                                  "Enhancer model learned from known enhancers and known negative\n",
                                  "\t\tregions. It is the output file from \"REPTILE_train.R\".\n",
                                  "\t\texample:\n",
                                  "\t\t  enhancer_model.reptile\n")
                              ),
        
        optparse::make_option(c("-s", "--sample-for-prediction"),
                              type="character",
                              metavar="sample_for_prediction",
                              dest="sample_for_prediction",
                              help=paste0(
                                  "Sample in which the activities of query regions are to be predicted.\n",
                                  "\t\tformat:\n",
                                  "\t\t  Sampl name\n",
                                  "\t\texample:\n",
                                  "\t\t  E11_5_FB\n")
                              ),
        
        optparse::make_option(c("-r", "--reference-samples"),
                              type="character",
                              metavar="ref_sample",
                              dest="ref_sample",
                              default=NULL,
                              help=paste0(
                                  "Samples used as reference to calculate intensity deviation\n",
                                  "\t\tformat:\n",
                                  "\t\t  Sample names separated by comma.\n",
                                  "\t\texample:\n",
                                  "\t\t  E11_5_FB,E11_5_HT,E11_5_MB\n")
                              ),

        optparse::make_option(c("-a", "--query-region-epimark-file"),
                              type="character",
                              metavar="query_region_epimark_file",
                              dest="query_region_epimark_file",
                              help=paste0(
                                  "Tab-separated file of epigenetic profiles of regions to be predicted.\n",
                                  "\t\tIf \"-n\" option is used and the number is greater than 1, this option\n",
                                  "\t\tshould be set as the prefix of query region files. For example, \n",
                                  "\t\tif this option is \"mm10_genome.region_epimark\" and \"-n\" is 8,\n",
                                  "\t\tREPTILE expects input files:\n",
                                  "\t\t  mm10_genome.region_epimark.0.tsv\n",
                                  "\t\t  mm10_genome.region_epimark.1.tsv\n",
                                  "\t\t  ...\n",
                                  "\t\t  mm10_genome.region_epimark.6.tsv\n",
                                  "\t\t  mm10_genome.region_epimark.7.tsv\n",
                                  "\t\tThe format of the file(s) is:\n",
                                  "\t\t  First line is a header indicating the content of each\n",
                                  "\t\t  column. First four columns are chromosome, start, end and\n",
                                  "\t\t  region id. Following columns are the scores/enrichment of\n",
                                  "\t\t  epigenetic marks for annotated regions.\n",
                                  "\t\texample:\n",
                                  "\t\t  chr start end id H3K4me1_H1 H3K4me1_H9 H3K4me1_IMR90 ...\n",
                                  "\t\t  chr1 3172000 3173000 reg_0 1.43 1.50 0.03 ...\n",
                                  "\t\t  ...\n",
                                  "\t\t  chr9 124412000 124413000 reg_302031 0.34 0.44 2.42 ...\n",
                                  "\t    **Highly recommend to to use \"REPTILE_preprocess.py\" script to generate\n",
                                  "\t      this file. Run \"./REPTILE_preprocess.py -h\" for more information.\n",
                                  "\t      Make sure the same data info file is used.\n")
                              ),
        
        optparse::make_option(c("-d", "--DMR-epimark-file"),
                              type="character",
                              metavar="DMR_epimark_file",
                              dest="DMR_epimark_file",
                              help=paste0(
                                  "Tab-separated file of epigenetic profiles of differentially\n",
                                  "\t\tmethyated regions (DMRs). Format is same as that of \n",
                                  "\t\tannotated_region_epimark_file. DMRs are used as high-resolution\n",
                                  "\t\tenhancer candidates to increase the resolution of training and\n",
                                  "\t\tprediction. The candidate loci can also come from assays like\n",
                                  "\t\tDNase-seq or ATAC-seq or analysis on motifs.\n",
                                  "\t\tSame as \"-a\": If \"-n\" option is used and the number is\n",
                                  "\t\tgreater than 1, this option should be set as the prefix of files.\n",
                                  "\t\tFor example, if this option is \"mm10_genome.DMR_epimark\" and \"-n\" is 8\n",
                                  "\t\tREPTILE expects input files:\n",
                                  "\t\t  mm10_genome.DMR_epimark.0.tsv\n",
                                  "\t\t  mm10_genome.DMR_epimark.1.tsv\n",
                                  "\t\t  ...\n",
                                  "\t\t  mm10_genome.DMR_epimark.6.tsv\n",
                                  "\t\t  mm10_genome.DMR_epimark.7.tsv\n",
                                  "\t\tThe format of the file(s) is:\n",
                                  "\t\t  First line is a header indicating the content of each\n",
                                  "\t\t  column. First four columns are chromosome, start, end and\n",
                                  "\t\t  region id. Following columns are the scores/enrichment of\n",
                                  "\t\t  epigenetic marks for annotated regions.\n",
                                  "\t\texample:\n",
                                  "\t\t  chr start end id H3K4me1_H1 H3K4me1_H9 H3K4me1_IMR90 ...\n",
                                  "\t\t  chr1 3172266 3172488 dmr_0 1.43 1.50 0.03 ...\n",
                                  "\t\t  ...\n",
                                  "\t\t  chr19 61316546 61316778 dmr_513260 0.34 0.44 2.42 ...\n",
                                  "\t    **Highly recommend to to use \"REPTILE_preprocess.py\" script to generate\n",
                                  "\t      this file. Run \"./REPTILE_preprocess.py -h\" for more information.\n",
                                  "\t      Make sure the same data info file is used.\n")
                              ),
        
        optparse::make_option(c("-o", "--output-prefix"),
                              type="character",
                              metavar="output_prefix",
                              dest="output_prefix",
                              default="REPTILE_output",
                              help=paste0(
                                  "Prefix of the output file storing the predictions.\n",
                                  "\t\tThe output file is \"<OUTPUT_PREFIX>.DMR.bed\"\n",
                                  "\t\tand \"<OUTPUT_PREFIX>.R.bed\" corresponding to the predictions\n",
                                  "\t\tfrom the classifier for DMRs on DMRs and classifier for query\n",
                                  "\t\tregions on query regions, respectively. If \"-w\" option is used,\n",
                                  "\t\tREPTILE will also generate \"<OUTPUT_PREFIX>.D.bed\" file as the\n",
                                  "\t\tcombined prediction. This file will store the enhancer score of\n",
                                  "\t\teach query region,which is defined as the maximum of scores of\n",
                                  "\t\tthe whole region and the overlapping DMRs. Also, when \"-w\" is enabled,\n",
                                  "\t\t\"<OUTPUT_PREFIX>.DMR.bed\" will only report the scores of DMRs that\n",
                                  "\t\toverlap with any query regions.\n",
                                  "\t\texample:\n",
                                  "\t\t  With \"--output-prefix E11_5_FB_enhancer\", two output files\n",
                                  "\t\t  are:\n",
                                  "\t\t  - E11_5_FB_enhancer.DMR.bed\n",
                                  "\t\t  - E11_5_FB_enhancer.R.bed\n",
                                  "\t\t  - E11_5_FB_enhancer.D.bed (if \"-w\" is used)\n",
                                  "\t    **output format:\n",
                                  "\t      BED format with 5th column as the enhancer confidence score.\n"
                                  )
                              ),
        
        optparse::make_option(c("-c", "--classifier-family"),
                              type="character",
                              default="RandomForest",
                              metavar="classifier_family",
                              dest="classifier_family",
                              help=paste0(
                                  "Classifier family to use in prediction model. \"Logistic\" is\n",
                                  "\t\tfaster but slightly less accurate compared to \"RandomForest\".\n",
                                  "\t\tIt is useful in generating quick and good results.\n",
                                  "\t\t  default: RandomForest\n",
                                  "\t\tClassifiers available:\n",
                                  "\t\t - RandomForest: random forest\n",
                                  "\t\t - Logistic: logistic regression\n"
                                  )
                              ),                                                            
        
        optparse::make_option(c("-x", "--no-intensity-deviation"),
                              type='logical',
                              action = "store_false",
                              default=TRUE,
                              metavar="incl_dev",
                              dest="incl_dev",
                              help=paste0(
                                  "If this option is used, REPTILE will not compute the intensity\n",
                                  "\t\tdeviation feature, which captures the tissue-specificity of\n",
                                  "\t\tenhancers.\n")
                              ),

        optparse::make_option(c("-w", "--not-genome-wide-prediction"),
                              action = "store_true",
                              default = FALSE,
                              type="logical",
                              metavar="not_genome_wide",
                              dest="not_genome_wide",
                              help=paste0(
                                  "This option is used to generate enhancer scores for given regions\n",
                                  "\t\tby combining the results from two classifiers (one for DMRs and.\n",
                                  "\t\tone for query regions). If this option is used, REPTILE will\n",
                                  "\t\tgenerate one additional file, \"<OUTPUT_PREFIX>.D.bed\", to show\n",
                                  "\t\tthe combined scores of each query region. The combined score is\n",
                                  "\t\tdefined as the maximum of scores of the whole region and the\n",
                                  "\t\toverlapping DMRs.\n",
                                  "\t\t  To use this option, the input files should be generated without\n",
                                  "\t\tthe \"-g\" option during preprocessing step (\"REPTILE_preprocess.py\").\n",
                                  "\t\tThat is, the ids in DMR_epimark_file (\"-d\") should contain the ids of\n",
                                  "\t\toverlapping query regions. For example, \"dmr_304:reg_1750\".\n"
                                  )
                              ),

        optparse::make_option(c("-p", "--number-of-processors"),
                              type="integer",
                              default=1,
                              metavar="num_procs",
                              dest="num_procs",
                              help=paste0(
                                  "Number of processors to use. REPTILE supports multiprocessing\n"
                                  )
                              ),

        optparse::make_option(c("-n", "--number-of-splits"),
                              type="integer",
                              default=1,
                              metavar="num_splits",
                              dest="num_splits",
                              help=paste0(                                  
                                  "REPTILE accepts splitted input files and this option\n",
                                  "\t\tspecifies the the number of DMR epimark files (or\n",
                                  "\t\tquery region epimark files after the split). For\n",
                                  "\t\texample, if \"-a\" is \"mm10_genome.region_epimark\" and\n",
                                  "\t\t\"-d\" is \"mm10_genome.DMR_epimark and this option is 8,\n",
                                  "\t\tREPTILE will search for below input files:\n",
                                  "\t\t  mm10_genome.DMR_epimark.0.tsv\n",
                                  "\t\t  ...\n",
                                  "\t\t  mm10_genome.DMR_epimark.7.tsv\n",
                                  "\t\tand\n",
                                  "\t\t  mm10_genome.region_epimark.0.tsv\n",
                                  "\t\t  ...\n",
                                  "\t\t  mm10_genome.region_epimark.7.tsv\n\n",
                                  "\t\tIf this option is 1, there is no split and REPTILE will\n",
                                  "\t\tassume \"-a\" and \"-d\" options give the exact names of\n",
                                  "\t\tof input files. Thus, it will search for two files:\n",
                                  "\t\t  mm10_genome.DMR_epimark\n",
                                  "\t\tand\n",
                                  "\t\t  mm10_genome.region_epimark\n"
                                  )
                              )
        )

    description <- paste0("\tPredicting enhancer activities of query regions in target sample. This script will\n",
                          "\tcalculate enhancer confidence score for each query region and each DMR.\n",
                          "\tPlease email Yupeng He (yupeng.he.bioinfo at gmail) for feedbacks, questions or bugs.")

    usage <- paste0("Usage: ./REPTILE_compute_score.R \\\n",
                    "\t\t -i data_info_file \\\n",
                    "\t\t -m model_file \\\n",
                    "\t\t -a query_region_epimark_file \\\n",
                    "\t\t -d DMR_epimark_file \\\n",
                    "\t\t -s sample_for_prediction \\\n",
                    "\t\t -o output_prefix\n"
                    )

    option_parser <- optparse::OptionParser(usage = usage,
                                            description = description,
                                            option_list=option_list)
    return(option_parser)
}

get_option_parser_evaluation <- function(){
    ##
    ## Function to generate option parser for performance evaluation
    ##

    suppressPackageStartupMessages(requireNamespace("optparse",
                                                    quietly = TRUE))

    option_list <- list(
        optparse::make_option(c("-s", "--sample-for-prediction"),
                              type="character",
                              metavar="sample_for_prediction",
                              dest="sample_for_prediction",
                              help=paste0(
                                  "Sample where the enhancer confidence scores are generated for query\n",
                                  "\t\tregions.\n",
                                  "\t\tformat:\n",
                                  "\t\t  Sampl name\n",
                                  "\t\texample:\n",
                                  "\t\t  E11_5_FB\n")
                              ),

        optparse::make_option(c("-l", "--label-file"),
                              type="character",
                              metavar="label_file",
                              dest="label_file",
                              help=paste0(
                                  "Tab-separated file labeling what query (annotated) regions are active\n",
                                  "\t\tactive enhancers in which sample(s).\n",
                                  "\t\tformat:\n",
                                  "\t\t  The file has multiple columns. First column is region id.\n",
                                  "\t\t  Each following column corresponds to one sample and the values\n",
                                  "\t\t  indicate whether the region is active in the sample:\n",
                                  "\t\t    - 1: active,\n",
                                  "\t\t    - 0: no activity\n",
                                  "\t\t    - NA: unknown\n",
                                  "\t\t  First line is a header indicating the content of each column.\n",
                                  "\t\t  Name of first column is \"id\" and others are sample names.\n",
                                  "\t\texample:\n",
                                  "\t\t  id E11_5_FB E11_5_MB E11_5_HB ...\n",
                                  "\t\t  reg_0 1 1 0 ...\n",
                                  "\t\t  ...\n",
                                  "\t\t  reg_302031 0 0 1 ....\n")
                              ),
        
        optparse::make_option(c("-p", "--prediction-result-file"),
                              type="character",
                              metavar="prediction_result_file",
                              dest="prediction_result_file",
                              help=paste0(
                                  "bed file of query regions with enhancer scores (5th column) such as the\n",
                                  "\t\t\"*.R.bed\" files and \"*.D.bed\" files from REPTILE_compute_score.R.\n",
                                  "\t\tformat:\n",
                                  "\t\t  BED file with five columns, chromosome, start, end, region id\n",
                                  "\t\t  and score of enhancer activity.\n",
                                  "\t\texample:\n",
                                  "\t\t  chr1 3172000 3173000 reg_0 0.7 ...\n",
                                  "\t\t  ...\n",
                                  "\t\t  chr9 124412000 124413000 reg_302031 0.1 ...\n",
                                  "\t    **Output files from \"REPTILE_compute_score.R\" are in this format\n",
                                  "\t      and can be directly used as input for this script\n")
                              ),
        
        optparse::make_option(c("-q", "--suppress-warnings"),
                              action = "store_true",
                              default = FALSE,
                              type = "logical",
                              dest = "suppress_warnings",
                              help=paste0(
                                  "Showing warnings is default. If this option is enabled, warnings will\n",
                                  "\t\tnot be shown."
                                  )
                              )
        )

    description <- paste0("\tEvaluating the prediction accuracy by calculating the Area Under Receiver Operating.\n",
                          "\tCharacteristic (AUROC) and the Area Under Precision-Recall curve (AUPR). It will also\n",
                          "\tcalculate the percent of true positives in the top 5, 10 and 20 predictions\n",
                          "\tPlease email Yupeng He (yupeng.he.bioinfo at gmail) for feedbacks, questions or bugs.")
    
    usage <- paste0("Usage: ./REPTILE_evaluate_prediction.R \\\n",
                    "\t\t -s sample_for_prediction \\\n",
                    "\t\t -l label_file \\\n",
                    "\t\t -p query_region_file_with_score\n"
                    )
    
    option_parser <- optparse::OptionParser(usage = usage,
                                            description = description,
                                            option_list=option_list)
    return(option_parser)
}

