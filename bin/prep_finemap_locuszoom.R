#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
credCutoff=0.01
#args = c("/home/sukmb465/Documents/Eike/Nextflow/finepipe/Results/GerNorItaSpa/2/cred/", "sss", 1, 2, "GerNorItaSpa")
#args = list("SCZ","CD")
if (length(args)!=4) {
  stop("You need: tag1 sss_or_cond n-causal-snps", call.=FALSE)
} 
library(data.table)
tag <- args[[4]]

###
output_dir <- paste0(args[[1]])
#chunk_tbl <- fread(paste0("chunk_tbl.tsv"))
###über alle chunks mitfilfe eines Tags
fmapforlocuszoom<-function(tag=tag,output_dir){
  #for (j in 1:dim(chunk_tbl)[1]) {
    #filename_sum_stats <- paste0("results_",args[[1]],"/",args[[1]],".z")
    #filename_cred <- paste0("results_",args[[1]],"/",args[[1]],".cred") 
    n_highest_given<-as.numeric(args[[3]])
    if(args[[2]]=="sss"){
      #count down to find the highest available n-causal if this is not ideal, reduce setting in pipeline starting command
      while (!file.exists(paste0(output_dir,tag,".cred",as.character(n_highest_given)))&n_highest_given>1) {
        n_highest_given <- n_highest_given - 1
      }
    }
    n_highest_given<-as.character(n_highest_given)
    if(file.exists(paste0(tag,".z"))&(file.exists(paste0(tag,".cred"))|file.exists(paste0(tag,".cred",n_highest_given)))){
      filename_sum_stats <- paste0(tag,".z")
      print(paste0(tag,".z"))
      # if(file.exists(paste0(tag,".cred"))){
      #   filename_cred <- paste0(tag,".cred")
      }
      if(file.exists(paste0(tag,".cred",n_highest_given))){
        filename_cred <- paste0(tag,".cred",n_highest_given)
        print(paste0(tag,".cred",n_highest_given))
      }else{filename_cred <- paste0(tag,".cred")}
    
      #filename_cred <- paste0(output_dir,"results_",tag,"/",tag,".cred")
      sum_stats <- fread(filename_sum_stats)
      print("sumstats read in")
      finemap_result <- fread(filename_cred, skip = "index")
      
      n_signals <- (ncol(finemap_result)-1)/2
      len_fmap <- nrow(finemap_result)[1]
      
      snp_vector<- vector(mode="character",length=0)
      pp_vector<- vector(mode="numeric",length=0)
      group_vector<- vector(mode="character",length=0)
      color_vector<- vector(mode="character",length=0)
      for (i in 1:n_signals) {
        colnames(finemap_result)[2*i] <- paste0("set",i)
        colnames(finemap_result)[2*i+1] <- paste0("pp_set",i)
        colnum1<-2*i
        colnum2<-2*i+1
        snp_vector<- c(snp_vector, as.vector(as.matrix(finemap_result[,..colnum1])))
        pp_vector<- c(pp_vector, as.vector(as.matrix(finemap_result[,..colnum2])))
        group_vector<-c(group_vector,rep(paste0("Set",i),len_fmap))
        color_vector<-c(color_vector,rep(switch(i,"red","blue","purple","green","orange"),len_fmap))
      }
      
      finemap_reformat <- data.table(snp=snp_vector, pp=pp_vector, group=group_vector, color=color_vector)
      
      if (exists("credCutoff")){
      	finemap_reformat <- finemap_reformat[pp>=credCutoff,] #if you take this out, the 95% credible set will contain some genes which have less than 1 % pp
      }
      finemap_reformat <- finemap_reformat[sum_stats[,.(rsid,chromosome,position)],on="snp==rsid",nomatch=NULL]
      finemap_reformat<-finemap_reformat[,.(snp,chromosome,position,pp,group,color)]
      colnames(finemap_reformat) <- c("snp","chr","pos","pp","group","color")
      
      filename_finemap <- paste0(tag,".cred.zoom")
      fwrite(finemap_reformat, file = filename_finemap, sep = "\t")
    #}
  }
#}
fmapforlocuszoom(tag = args[[4]], output_dir = output_dir)
#fmapforlocuszoom(chunk_tbl = chunk_tbl, tag = args[[2]], output_dir = output_dir)
