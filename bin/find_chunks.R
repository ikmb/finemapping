#!/usr/bin/env Rscript
###This script shall define the boundaries for the FINEMAP-Tool. It uses the given gene_start and gene_end columns from the .olap file.
###Currently, you get a +-1M frame around the gene start, which is approx. the TSS. It is taken care that long genes are fully covered.
###Then, genes are chunked to capture any shared QTLs that might exist. The boundaries are adapted to PLINK.
#expected input: Overlapdoc
#like Rscript find_chunks.R /path/to/overlap.olap path/to/sumstats1 path/to/sumstats2 tag1 tag2)
args = commandArgs(trailingOnly=TRUE)
#args = list("/home/user/Desktop/GBXMap/UC.info.0.4.sumstats_SCZ.info.0.6.maf.0.01.sumstats_div10.segbfs","sum_stats/UC.sumstats.info04.filtered.withheader",
#            "sum_stats/SCZ.info.0.6.maf.0.01.sumstats.filtered.withheader.fmap","UC","SCZ")
# test if there are five arguments: if not, return an error
#if (length(args)!=3) {
#  stop("You need: path/to/.olap path/to/sumstats1 tag1", call.=FALSE)
#} 
#output_dir <- paste0(args[[3]],"_","output/")
#if(!dir.exists(output_dir)){dir.create(output_dir)}else{stop("your desired output directory already exists. Rename your directory to avoid overwriting of your results.")}
library(data.table)
AB_overlap <- fread(args[[1]])
sumstats1<-fread(args[[2]])
ref<-args[[3]]
snps<-fread(args[[4]],header=FALSE)$V1#added to consider snp list (in conflict with UTMOST case -> either remove l76ff or change args)
print("read in all files")
#sumstats2<-fread(args[[3]])
#tag1<-args[[3]]
#tag2<-args[[5]]
if(colnames(AB_overlap)[1]=="chunk"){chunk_mode<-"PWGWAS"}else{chunk_mode<-"UTMOST"}
### organize the segbfs that it generates a table with the exact columns of the UTMOST mode
if(chunk_mode=="PWGWAS"){
  AB_overlap<-AB_overlap[,c("chunk","NSNP","chr","st","sp","PPA_3")]
  setnames(AB_overlap,"chunk","chunk_PW")
  AB_overlap[,chr:=as.numeric(gsub("chr","",chr))]
  setorder(AB_overlap,chr,st)
  AB_overlap<-AB_overlap[PPA_3>0.5]
  AB_overlap[,chunk:=1]
  chunk_num <- 1
  j <- 1
  ###following loop should see if there are chunks following directly on each other in the PWGWAS output, so they will be interpreted as one chunk
  ### is not useful when doing this pipeline for GWAS purposes, because numbering of loci is arbitrary then
  ###therefore, condition is deactivated
  for (i in 1:(dim(AB_overlap)[1])) {
    AB_overlap[j,"chunk"] <- chunk_num
    if(j<dim(AB_overlap)[1]){
      if (AB_overlap[j+1,chr] == AB_overlap[j,chr] && (AB_overlap[j+1,chunk_PW] == AB_overlap[j,chunk_PW] + 1) && FALSE){ #condition deactivated by insertion of &&F
        chunk_num <- chunk_num}
      else{chunk_num<-chunk_num+1}
      j<-j+1
    }
  }
  print(AB_overlap)
  fwrite(AB_overlap, file = paste0("chunk_chunk.csv"), sep = ",")
  print("Written chunk_chunk.csv")
  #group by chunk to 
  n_chunks <- max(AB_overlap[,chunk])
  chunk_tbl <- data.table(chunk=1:n_chunks, 
                          first_gene_start=AB_overlap[,min(st),by=chunk][,V1],
                          last_gene_end=AB_overlap[,max(sp),by=chunk][,V1],
                          CHR=unique(AB_overlap[,chr,by="chunk"])[,chr],
                          left_cis_boundary=AB_overlap[,min(st),by=chunk][,V1], 
                          right_cis_boundary=AB_overlap[,max(sp),by=chunk][,V1])
  chunk_tbl[,plink_left_bound:=left_cis_boundary%/%1000][,plink_right_bound:=(right_cis_boundary+999)%/%1000
                                                         ][,subset_left_bound:=plink_left_bound*1000
                                                           ][,subset_right_bound:=plink_right_bound*1000]
  n_genes=AB_overlap[, .N, by=chunk][,N]
  chunk_tbl<-cbind(chunk_tbl, n_genes)
  print(chunk_tbl)
  findleadSNP<-function(sumstats,chromosome,subset_left_bound,subset_right_bound,snps){
    DT<-sumstats[CHR==chromosome&BP>=subset_left_bound&BP<=subset_right_bound&SNP%in%snps,.(SNP,P)]#[P==min(P),SNP]
    leadSNP<-DT[P==min(P),SNP][1]
    print(leadSNP)
    return(leadSNP)
  }
  print("looking for the leadSNP...")
  chunk_tbl[,paste0("leadSNP_","1"#,ref
  ):=mapply(findleadSNP,chromosome=CHR,subset_left_bound=subset_left_bound,
                                                  subset_right_bound=subset_right_bound, MoreArgs = list(sumstats=sumstats1,snps=snps))]
#  chunk_tbl[,paste0("leadSNP_",args[[5]]):=mapply(findleadSNP,chromosome=CHR,subset_left_bound=subset_left_bound,
#                                                  subset_right_bound=subset_right_bound, MoreArgs = list(sumstats=sumstats2))]
  ###
  fwrite(chunk_tbl,file=paste0("chunk_tbl.csv"),sep = ",",scipen=100)
}

####
if(chunk_mode=="UTMOST"){
  setorder(AB_overlap, CHR, gene_start)
  ###1M from gene start which is approx TSS is eQTL window.
  tssRightBoundary<-function(gene_start,gene_end){
    if(gene_end>(gene_start+900000)){return (gene_end+100000)}else{return (gene_start+1000000)}
  }
  AB_overlap[,left_cis_boundary:=gene_start-1000000][,right_cis_boundary:=mapply(tssRightBoundary,gene_start,gene_end)]
  #chunk the SNPS if closer than 250000bp to the next one
  AB_overlap[,chunk:=1]
  chunk_num <- 1
  j <- 1
  for (i in 1:(dim(AB_overlap)[1])) {
    AB_overlap[j,"chunk"] <- chunk_num
    if(j<dim(AB_overlap)[1]){
      if (AB_overlap[j+1,CHR]==AB_overlap[j,CHR] && AB_overlap[j+1,left_cis_boundary] < AB_overlap[j,right_cis_boundary] + 250000){
        chunk_num <- chunk_num}
      else{chunk_num<-chunk_num+1}
      j<-j+1
    }
  }
  fwrite(AB_overlap, file = paste0("gene_chunk.csv"), sep = ",")
  print("written gene_chunk")
  #group by chunk to 
  n_chunks <- max(AB_overlap[,chunk])
  chunk_tbl <- data.table(chunk=1:n_chunks, 
                          first_gene_start=AB_overlap[,min(gene_start),by=chunk][,V1],
                          last_gene_end=AB_overlap[,max(gene_end),by=chunk][,V1],
                          CHR=unique(AB_overlap[,CHR,by="chunk"])[,CHR],
                          left_cis_boundary=AB_overlap[,min(left_cis_boundary),by=chunk][,V1], 
                          right_cis_boundary=AB_overlap[,max(right_cis_boundary),by=chunk][,V1])
  
  chunk_tbl[,plink_left_bound:=left_cis_boundary%/%1000][,plink_right_bound:=(right_cis_boundary+999)%/%1000
                                                         ][,subset_left_bound:=plink_left_bound*1000
                                                           ][,subset_right_bound:=plink_right_bound*1000]
  n_genes=AB_overlap[, .N, by=chunk][,N]
  chunk_tbl<-cbind(chunk_tbl, n_genes)
  ###find the lead SNP in each chunk which is at least 100kb from the side for each disease. 
  findleadSNP<-function(sumstats,chromosome,subset_left_bound,subset_right_bound){
    DT<-sumstats[CHR==chromosome&BP>=subset_left_bound&BP<=subset_right_bound,.(SNP,P)]#[P==min(P),SNP]
    leadSNP<-DT[P==min(P),SNP][1]
    print(leadSNP)
    return(leadSNP)
  }
  ###fallback if the code below with mapply should cause problems
  #snpvector1<-vector(mode = "character", length = dim(chunk_tbl)[1])
  #snpvector2<-vector(mode = "character", length = dim(chunk_tbl)[1])
  #for (i in 1:dim(chunk_tbl)[[1]]) {
  #  snpvector1[i]<-findleadSNP(sumstats1,chunk_tbl[i,CHR],chunk_tbl[i,subset_left_bound],chunk_tbl[i,subset_right_bound])
  #  snpvector2[i]<-findleadSNP(sumstats2,chunk_tbl[i,CHR],chunk_tbl[i,subset_left_bound],chunk_tbl[i,subset_right_bound])
  #}
  #chunk_tbl[,paste0("leadSNP_",args[[4]]):=snpvector1]
  #chunk_tbl[,paste0("leadSNP_",args[[5]]):=snpvector2]
  ###should work
  chunk_tbl[,paste0("leadSNP_",args[[4]]):=mapply(findleadSNP,chromosome=CHR,subset_left_bound=subset_left_bound,
                                                  subset_right_bound=subset_right_bound, MoreArgs = list(sumstats=sumstats1))]
  chunk_tbl[,paste0("leadSNP_",args[[5]]):=mapply(findleadSNP,chromosome=CHR,subset_left_bound=subset_left_bound,
                                                  subset_right_bound=subset_right_bound, MoreArgs = list(sumstats=sumstats2))]
  
  ###
  fwrite(chunk_tbl,file=paste0("chunk_tbl.csv"),sep = ",",scipen=100)
}
