#!/usr/bin/env Rscript
######################################################################################
## REPTILE
## - Regulatory Element Prediction based on TIssue-specific Local Epigenetic marks
## by Yupeng He (yupeng.he.bioinfo at gmail)
## Salk Institute for Biological Studies
######################################################################################

##  Predicting enhancer activities of query regions in target sample
##  Please email Yupeng He (yupeng.he.bioinfo at gmail) for feedback, question and bug. 

## option parsing
suppressPackageStartupMessages(library("REPTILE",verbose=FALSE))

option_parser = get_option_parser_compute_score()
options <- optparse::parse_args(option_parser)
## No input
if(length(commandArgs(TRUE)) == 0){
    optparse::print_help(option_parser)
    quit("no")
}
##End of option parsing
######################################################################################
## Input from command line options
##Data information file
data_info_file = options$data_info_file
##Model file
model_file = options$model_file
##Sample in which prediction is made
query_sample = options$sample_for_prediction
##Reference sample(s)
ref_sample = options$ref_sample
if(!is.null(ref_sample)){
    ref_sample = unlist(strsplit(options$ref_sample,","))
}
##Epimark file of DMRs
DMR_epimark_file = options$DMR_epimark_file
##Epimark file of input regions
query_region_epimark_file = options$query_region_epimark_file
##Prefix of output files
output_prefix = options$output_prefix
##Classifier family to use
classifier_family = options$classifier_family
##Whether to compute deviation as additional features
incl_dev = options$incl_dev
##Number of precessors to use
num_procs = options$num_procs
##Number of splits
num_splits = options$num_splits
##Whether to combine score of query regions and the scores from overlapping DMRs
not_genome_wide = options$not_genome_wide

## Load enhancer model
load(model_file)

## Read data informaiton
data_info = read.table(data_info_file,header=T,stringsAsFactors=F)[,1:2]

## No split
if(num_splits == 1){
    ## Check if the files exist
    expected_files <- c(query_region_epimark_file)
    if(!is.null(DMR_epimark_file)){
        expected_files <- c(expected_files,DMR_epimark_file)
    }
    if_exist = file.exists(expected_files)
    if(sum(if_exist) != length(expected_files)){
        stop(paste0("Input files cannot be found! Please check the input options.\n\n",
                    "- query region epimark file is:\n",
                    "    ",query_region_epimark_file,"\n\n",
                    "- DMR epimark file is:\n",
                    "    ",DMR_epimark_file,"\n\n",
                    "- Below expected input file(s) cannot be found:\n",
                    paste0("    ",expected_files[!if_exist],"\n",collapse="")
                    ),
             "\n"
             )
    }
    remove(expected_files,if_exist)
    
    ## Begin
    epimark_region <- read_epigenomic_data(data_info,
                                           query_region_epimark_file,
                                           query_sample=query_sample,
                                           ref_sample=ref_sample,
                                           incl_dev=incl_dev
                                           )
    epimark_region = epimark_region[rowSums(is.na(epimark_region))==0,]
    
    epimark_DMR = NULL
    if(!is.null(DMR_epimark_file)){
        epimark_DMR <- read_epigenomic_data(data_info,
                                            DMR_epimark_file,
                                            query_sample=query_sample,
                                            ref_sample=ref_sample,
                                            incl_dev=incl_dev
                                            )
        if(!is.null(epimark_DMR)){
            epimark_DMR = epimark_DMR[rowSums(is.na(epimark_DMR))==0,]
        }
    }
    
    ## Make predictions using the two classifiers
    if(not_genome_wide){
        pred <- reptile_predict(reptile,epimark_region,epimark_DMR,
                                family=classifier_family)
    }else{
        pred <- reptile_predict_genome_wide(reptile,epimark_region,epimark_DMR,
                                            family=classifier_family)
    }
    remove(epimark_region,epimark_DMR)
    
    ## Output
    ## Query regions
    query_region <- read.table(query_region_epimark_file,
                               header=T,stringsAsFactors=F)
    rownames(query_region) = query_region$id
    query_region = query_region[,1:3]

    if(not_genome_wide & (!is.null(DMR_epimark_file))){
        ## *.D.bed file with combined scores
        write.table(
            cbind(query_region,rownames(query_region),
                  pred$D[rownames(query_region)]),
            paste0(output_prefix,".D.bed"),
            quote=F,row.names=F,col.names=F,sep="\t")
    }
    ## *.R.bed file with score based on the epimark of whole regions
    write.table(
        cbind(query_region,rownames(query_region),
              pred$R[rownames(query_region)]),
        paste0(output_prefix,".R.bed"),
        quote=F,row.names=F,col.names=F,sep="\t")    
    
    ## DMRs
    if(!is.null(DMR_epimark_file)){
        DMR = read.table(DMR_epimark_file,header=T,stringsAsFactors=F)
        if(not_genome_wide){
            ## Get unique DMR ids
            DMR_id_uniq = unique(sapply(strsplit(DMR$id,":"),function(x) return(x[1])))
            ## Get unique DMRs
            DMR = aggregate(DMR[,1:3],list(sapply(strsplit(DMR$id,":"),function(x) return(x[1]))),max)
            DMR_id = DMR[,1]
            DMR = DMR[,2:4]
            rownames(DMR) = DMR_id
            
            pred_DMR = pred$DMR
            remove(pred)
            names(pred_DMR) = sapply(strsplit(names(pred_DMR),":"),function(x) return(x[1]))
            
            write.table(
                cbind(DMR[DMR_id_uniq,],DMR_id_uniq,
                      pred_DMR[DMR_id_uniq]),
                paste0(output_prefix,".DMR.bed"),
                quote=F,row.names=F,col.names=F,sep="\t")
        }else{
            DMR_id = DMR[,4]
            DMR = DMR[,1:3]
            rownames(DMR) = DMR_id
            
            write.table(
                cbind(DMR,DMR_id,
                      pred$DMR[DMR_id]),
                paste0(output_prefix,".DMR.bed"),
                quote=F,row.names=F,col.names=F,sep="\t")
        }
    }
    
}else{
    ## Check if the files exist
    expected_files <- paste0(query_region_epimark_file,".",(1:num_splits)-1,".tsv")
    if(!is.null(DMR_epimark_file)){
        expected_files <- c(expected_files,
                            paste0(DMR_epimark_file,".",(1:num_splits)-1,".tsv"))
    }
    if_exist = file.exists(expected_files)
    if(sum(if_exist) != length(expected_files)){
        write(paste0("Error: Input files cannot be found! Please check the input options.\n",
                     "       Note that number of splits (\"-n\") is greater than 1. Now the\n",
                     "       \"-a\" and \"-d\" options are specifying the prefixes of input\n",
                     "       files!\n\n",
                     "- Prefix for query region epimark file is:\n",
                     "    ",query_region_epimark_file,"\n\n",
                     "- Prefix for DMR epimark file is:\n",
                     "    ",DMR_epimark_file,"\n\n",
                     "- Number of splits is:\n",
                     "    ",num_splits,"\n\n",
                     "- Below expected input file(s) cannot be found:\n",
                     paste0("    ",expected_files[!if_exist],"\n",collapse="")
                     ),            
              stderr()
              )
        stop("Input files cannot be found! Please check the input options.")
    }
    remove(expected_files,if_exist)

    ## Multiprocessing
    suppressPackageStartupMessages(library("doParallel",,quietly=T,verbose=F))   
    suppressPackageStartupMessages(library("foreach",quietly=T,verbose=F))

    cl <- makeCluster(num_procs)
    registerDoParallel(cl)
    foreach(s = 1:num_splits) %dopar% {
        ## Loading REPTILE library
        suppressPackageStartupMessages(library("REPTILE",verbose=FALSE))
        
        ## Read data
        epimark_region <- read_epigenomic_data(data_info,
                                               paste0(query_region_epimark_file,".",s-1,".tsv"),
                                               query_sample=query_sample,
                                               ref_sample=ref_sample,
                                               incl_dev=incl_dev
                                               )
        epimark_region = epimark_region[rowSums(is.na(epimark_region))==0,]

        epimark_DMR = NULL
        if(!is.null(DMR_epimark_file)){
            epimark_DMR <- read_epigenomic_data(data_info,
                                                paste0(DMR_epimark_file,".",s-1,".tsv"),
                                                query_sample=query_sample,
                                                ref_sample=ref_sample,
                                                incl_dev=incl_dev
                                                )
            if(!is.null(epimark_DMR)){
                epimark_DMR = epimark_DMR[rowSums(is.na(epimark_DMR))==0,]
            }
        }

        ## Make predictions
        if(not_genome_wide){
            pred <- reptile_predict(reptile,epimark_region,epimark_DMR,
                                                family=classifier_family)        
        }else{
            pred <- reptile_predict_genome_wide(reptile,epimark_region,epimark_DMR,
                                                family=classifier_family)
        }
        remove(epimark_region,epimark_DMR)

        ## Output
        ##Get coordinates of query regions
        query_region <- read.table(paste0(query_region_epimark_file,".",s-1,".tsv"),
                                   header=T,stringsAsFactors=F)

        rownames(query_region) = query_region$id
        query_region = query_region[,1:3]
        
        if(not_genome_wide & (!is.null(DMR_epimark_file))){
            write.table(
                cbind(query_region,rownames(query_region),
                      pred$D[rownames(query_region)]),
                paste0(output_prefix,".",s-1,".D.bed"),
                quote=F,row.names=F,col.names=F,sep="\t")
        }
        
        write.table(
            cbind(query_region,rownames(query_region),
                  pred$R[rownames(query_region)]),
            paste0(output_prefix,".",s-1,".R.bed"),
            quote=F,row.names=F,col.names=F,sep="\t")    
        
        ##Print out DMRs and their enhancer score
        if(!is.null(DMR_epimark_file)){
            DMR = read.table(
                paste0(DMR_epimark_file,".",s-1,".tsv"),
                header=T,stringsAsFactors=F)
            if(not_genome_wide){
                DMR_id_uniq = unique(sapply(strsplit(DMR$id,":"),function(x) return(x[1])))
                ## 
                DMR = aggregate(DMR[,1:3],list(sapply(strsplit(DMR$id,":"),function(x) return(x[1]))),max)
                DMR_id = DMR[,1]
                DMR = DMR[,2:4]
                rownames(DMR) = DMR_id
                
                pred_DMR = pred$DMR
                remove(pred)
                names(pred_DMR) = sapply(strsplit(names(pred_DMR),":"),function(x) return(x[1]))
                
                write.table(
                    cbind(DMR[DMR_id_uniq,],DMR_id_uniq,
                          pred_DMR[DMR_id_uniq]),
                    paste0(output_prefix,".",s-1,".DMR.bed"),
                    quote=F,row.names=F,col.names=F,sep="\t")
            }else{
                
                DMR_id = DMR[,4]
                DMR = DMR[,1:3]
                rownames(DMR) = DMR_id
                
                write.table(
                    cbind(DMR,DMR_id,
                          pred$DMR[DMR_id]),
                    paste0(output_prefix,".",s-1,".DMR.bed"),
                    quote=F,row.names=F,col.names=F,sep="\t")
            }
        }
    }
    stopCluster(cl)
    
    ## Combine output files
    ## DMR
    if(!is.null(DMR_epimark_file)){
        output_DMR_files = paste0(output_prefix,".",(1:num_splits)-1,".DMR.bed")
        cmd_combine <- paste0("cat ",paste(output_DMR_files,collapse=" "),
                              "> ",paste0(output_prefix,".DMR.bed"))
        system(cmd_combine)
        cmd_remove <- paste0("rm -r ",
                         paste(output_DMR_files,collapse=" "))
        system(cmd_remove)
    }

    ## D
    if(not_genome_wide & !is.null(DMR_epimark_file)){
        output_D_files = paste0(output_prefix,".",(1:num_splits)-1,".D.bed")
        cmd_combine <- paste0("cat ",paste(output_D_files,collapse=" "),
                              "> ",paste0(output_prefix,".D.bed"))
        system(cmd_combine)
        cmd_remove <- paste0("rm -r ",
                             paste(output_D_files,collapse=" "))
        system(cmd_remove)
    }
    ## R
    output_R_files = paste0(output_prefix,".",(1:num_splits)-1,".R.bed")
    cmd_combine <- paste0("cat ",paste(output_R_files,collapse=" "),
                          "> ",paste0(output_prefix,".R.bed"))
    system(cmd_combine)   
    cmd_remove <- paste0("rm -r ",
                         paste(output_R_files,collapse=" "))
    system(cmd_remove)
}
