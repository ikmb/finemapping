#!/usr/bin/env nextflow
// -*- mode:groovy -*-

//TODO:chunk size test
//TODO: check for chr3:45800039:A:G as it is right now and possible 1:752721:A:G input in .bim and snplist

def helpMessage() {
  log.info"""
  =================================================================
   IKMB | Finemap Pipeline | v${workflow.manifest.version}
  =================================================================
  Usage:
  The typical command for running the pipeline is as follows:
  nextflow run IKMB/Finepipe --locus '/path/to/locusfile.sample' 
  Mandatory arguments:
  --locus 		The path to the .csv-File containing following columns: chunk,NSNP,chr,st,sp,PPA_3

  --sumstats        Summary stats
                    Header must exactly be:
                    CHR	BP	SNP	A1	A2	P	OR	BETA	SE	N	CHISQ	Z	SOURCE  FRQ_A_A1	FRQ_U_A1	INFO
                    only CHR, BP, SNP, A1, A2, P, BETA, SE, FRQ_U_A1 are in use. 
                    The other columns can be specified arbitrarily.
                    CHR has to be like 1, not like chr1, and use 23 instead of X.
                    FRQ_U_A1 is the reference allele frequency of A1. 
  
  --reference       Shared file name of bim bed fam of a reference population. 

  --snps            A sorted file with all the SNPs in the loci 
                    (every other marker will be filtered out),
                    without the header of the reference bimbedfam like.  
                    Default is all rsids from reference bim file
  --nsum            Number of probands in reference population.

  --method          Method of Finemap (sss or cond)

  --nsignal         n of causal signals

  --locuszoomdb     location of locusZoom Db 
  --output          directory for result files

  

  Optonal arguments:
  --dprime        result in d-prime instead of r2
  -profile      	The nextflow execution profile to use, standard is medcluster ("local" or "standard")
  -resume         Continue the workflow
  """.stripIndent()
}

log.info "=========================================" 
log.info "F I N E M A P P I N G   P I P E L I N E"
log.info "IKMB pipeline version v${workflow.manifest.version}" 
log.info "Nextflow Version: 	$workflow.nextflow.version" 
log.info "=== Inputs =============================="
log.info "Boundries-File:		${params.locus}"
log.info "Summary stats:    ${params.sumstats}"
log.info "Reference population:   ${params.reference}"
log.info "SNPs in the loci: 		${params.snps}"
log.info "Number of probands: 		${params.nsum}"
log.info "Finemap method: 		${params.method}"
log.info "Number of signals: 		${params.nsignal}"
log.info "dprime instead of r2:    ${params.dprime}"
log.info "LocusZoom DB    ${params.locuszoomdb}"
log.info "Output directory: 		${params.output}"
log.info "=========================================="
log.info "Command Line:         	$workflow.commandLine"
if (workflow.containerEngine) {
	log.info "Container Engine: 	${workflow.containerEngine}"
}
log.info "=========================================" 

// Show help message
if (params.help){
	helpMessage()
	exit 0
}


process set_snps {
  scratch true
  input:
  output:
  file(thesnplist) into thesnplist_ch
  shell:
  thesnplist = "thesnplist.txt"
  """
  if [ -f "!{params.snps}" ]
  then
    sort -u !{params.snps}  > $thesnplist
  else  
    awk '{print \$2}' !{params.reference}.bim | sort -u > $thesnplist
  fi
  """

}

thesnplist_ch.into{ thesnplist_ch2;thesnplist_ch3 }

process find_chunks {
    scratch true
    label 'Rscript'
    input:
    file(thesnpslist) from thesnplist_ch2
    output:
    file(chunk_tbl) into chunk_tbl_ch
    
    script:
    
    chunk_chunk = "chunk_chunk.csv"
    chunk_tbl = "chunk_tbl.csv"
    """
    Rscript $baseDir/bin/find_chunks.R ${params.locus} ${params.sumstats} ${params.reference} ${thesnpslist}
    """ 
}

chunk_tbl_ch.splitCsv(header: true, sep:',').map{it ->[it.chunk,it.first_gene_start,it.last_gene_end,it.CHR,it.left_cis_boundary,it.right_cis_boundary,it.plink_left_bound,it.plink_right_bound,it.subset_left_bound,it.subset_right_bound,it.n_genes,it.leadSNP_1 ]}.into {chunks_ch2; chunks_ch3}


process sort_snplist {
  scratch true
  input:
  output:
  file(sorted) into snplist_sorted_ch
  shell:

  sorted = "reference_snplist_sorted.txt"
  """
  
  gawk '{print \$2}' !{params.reference}.bim | sort > reference_snplist_sorted.txt
  """
}

// -F \$'\t'
Channel.fromPath(params.sumstats, checkIfExists: true).into{ summarystats_ch; summarystats_ch2 }

process subset_sumstats {
    scratch true
    beforeScript 'ulimit -Ss unlimited'   
    input:
    tuple chunk,first_gene_start,last_gene_end,CHR,left_cis_boundary,right_cis_boundary,plink_left_bound,plink_right_bound,subset_left_bound,subset_right_bound,n_genes,leadSNP_1 from chunks_ch2
    each file(sorted) from snplist_sorted_ch
    each file(summarystats) from summarystats_ch
    each file(thesnpslist) from thesnplist_ch3
    output:
    tuple val(chunk),file(snplist) into key_snplist_ch
    tuple val(chunk),file(zfile) into key_zfile_ch
    tuple val(chunk), file(metal) into metal_ch
    shell:
    
    chromosome=CHR
    left_bound=subset_left_bound
    right_bound=subset_right_bound
  	filetag=chunk
    zraw=filetag+".txt"
    snplistraw=filetag+"_snplistraw"
    presnplist=filetag+"_presnplist"
    snplist=filetag+"_snplist"
    zfile=filetag+".z"
    metal=filetag+".metal"
    """
    #!/bin/bash
  	cat !{params.sumstats}| awk 'BEGIN {CONVFMT = "%.6e"; OFMT = "%.6e"; print "rsid chromosome position allele1 allele2 maf beta se"} NR!=1 {print \$3 " " \$1 " " \$2 " " \$4 " " \$5 " " \$14 " " \$8 " " \$9}' | gawk 'BEGIN {CONVFMT = "%.6e"; OFMT = "%.6e"; print "rsid chromosome position allele1 allele2 maf beta se"} (NR!=1 && \$2 == '$chromosome' && \$3 > '$left_bound' && \$3 < '$right_bound') {print }' | awk 'BEGIN {print "rsid chromosome position allele1 allele2 maf beta se"}{if (NR!=1 && \$6 <= 0.5) {print } ;if (NR!=1 && \$6 > 0.5) {mafcor = 1 - \$6; minusbeta = -1 * \$7; print \$1 " " \$2 " " \$3 " " \$5 " " \$4 " " mafcor " " minusbeta " " \$8}}' | awk 'BEGIN {print "rsid chromosome position allele1 allele2 maf beta se"}{if (NR!=1 && \$6 >= 0.0) {print }}' >> $zraw 
    cat $zraw | gawk 'NR!=1{print \$1}' | sort -u > $snplistraw
        #now save the real snplist as the comm of two snplists. Each one has to be sorted!!!
	  #CHANGED BY Mareike (to consider snplist): 
    comm -12 <(cat $snplistraw) <(cat $sorted) | sort -u > $presnplist
    comm -12 <(cat $presnplist) <(cat !{thesnpslist} | sort -u) > $snplist
	  #comm -12 <(cat $snplistraw) <(cat $sorted)|grep -F -f !{thesnpslist} > $snplist

        #now subset the zfile to be concordant with the snplist
	  head -n 1 $zraw > $zfile
	  grep -Fwf $snplist $zraw >> $zfile
	  echo -e 'MarkerName\tP-value' >> $metal

          grep -f $snplist !{params.sumstats} | awk -v OFS="\\t" '\$1=\$1' |cut -f3,6 >> $metal
	  #CHANGED BY Mareike (I guess the old version is better but failed for some reason) #join -1 1 -2 3 -o 1.1,2.6 -t \$'\t' <(sort $snplist) <(sort -k 3 !{params.sumstats}) >> $metal
    """
} 

key_zfile_ch.into{key_zfile_ch1;key_zfile_ch2}

key_zfile_ch2.join(key_snplist_ch).join(chunks_ch3).into{plink_input_ch;plink_input_ch2}

process plink {
    scratch true
    label 'plink'
    input:
    tuple val(chunk),file(zfile),file(snplist),first_gene_start,last_gene_end,CHR,left_cis_boundary,right_cis_boundary,plink_left_bound,plink_right_bound,subset_left_bound,subset_right_bound,n_genes,leadSNP_1 from plink_input_ch
    
    output:
    tuple val(chunk),file(ld) into plink_output
    shell:
    output=chunk
    ld=output+".ld"
    logfile=output+".log"
    """
    plink --bfile !{params.reference} --extract $snplist --r square 'spaces' --a1-allele $zfile 4 1 --from-kb $plink_left_bound --to-kb $plink_right_bound --chr $CHR --out "$output" --threads ${task.cpus} 
    """ 
}


key_zfile_ch1.join(plink_output).into{joinedzfileld_ch;joinedzfileld_ch2}

process master {
    scratch true
    input:
    tuple val(chunk),val(zfile), val(ldfile) from joinedzfileld_ch
    output:
    tuple val(chunk), file(output) into master_output_ch
    script:
    output=chunk+"_output.master"
    """
    echo "z;ld;snp;config;cred;log;n_samples" > $output
    echo "$zfile;$ldfile;${chunk}.snp;${chunk}.config;${chunk}.cred;${chunk}.log;${params.nsum}" >> $output
    """ 
}

Channel.fromPath(params.reference).map { file -> tuple(file.baseName, file) }.into{identifier_ch;identifier_ch2}

master_output_ch.set{master_output_ch1}

identifier_ch.combine(master_output_ch1).set{finemap_input_ch}

if(params.method == "sss"){
  process finemap_sss {
      scratch true
      label 'finemap'
      publishDir "${params.output}/${datasetID}/${chunk}/cred/", mode: 'copy'
      input:
      tuple datasetID, datasetFile, val(chunk), file(finemap_output) from finemap_input_ch
      output:
      file(error_finemap) optional true 
      file('*.cred*')
      file('*.cred') optional true
      tuple val(chunk), file(NAME) into finemap_output_ch
      shell:
      error_finemap="error_finemap.log"
      NAME=chunk+".cred"+params.nsignal
      
      '''
      NAME=!{chunk}".cred"!{params.nsignal}
      j=0
      while [ ! -f $NAME ] && [ $j -lt 10 ]
      do
        if [ "!{params.method}" == "sss" ]; then
          !{baseDir}/bin/finemap_v1.4_x86_64/finemap_v1.4_x86_64 --!{params.method} --in-files !{finemap_output} --log --n-causal-snps !{params.nsignal} 
          j=$(( $j + 1 ))
        else
          !{baseDir}/bin/finemap_v1.4_x86_64/finemap_v1.4_x86_64 --!{params.method} --in-files !{finemap_output} --log --n-causal-snps !{params.nsignal} --cond-pvalue 0.0001
          j=$(( $j + 1 ))
        fi
      done
      '''
  }
}else{
    process finemap_cond {
      scratch true
      label 'finemap'
      publishDir "${params.output}/${datasetID}/${chunk}/cred/", mode: 'copy'
      input:
      tuple datasetID, datasetFile, val(chunk), file(finemap_output) from finemap_input_ch
      output:
      file(error_finemap) optional true 
      file('*.cred*') optional true
      file('*.cred') optional true
      tuple val(chunk), file(NAME) into finemap_output_ch
      shell:
      error_finemap="error_finemap.log"
      NAME=chunk+".cred"
      
      '''
      NAME=!{chunk}".cred"!{params.nsignal}
      j=0
      while [ ! -f $NAME ] && [ $j -lt 10 ]
      do
        if [ "!{params.method}" == "sss" ]; then
          !{baseDir}/bin/finemap_v1.4_x86_64/finemap_v1.4_x86_64 --!{params.method} --in-files !{finemap_output} --log --n-causal-snps !{params.nsignal} 
          j=$(( $j + 1 ))
        else
          !{baseDir}/bin/finemap_v1.4_x86_64/finemap_v1.4_x86_64 --!{params.method} --in-files !{finemap_output} --log --n-causal-snps !{params.nsignal} --cond-pvalue 0.0001
          j=$(( $j + 1 ))
        fi
      done
      '''
  }
}


finemap_output_ch.combine(identifier_ch2).into{prep_finemap_locuszoom_ch;prep_finemap_locuszoom_ch2}

joinedzfileld_ch2.join(prep_finemap_locuszoom_ch).set{prep_finemap_locuszoom_input}

process prep_finemap_locuszoom {
    scratch true
    label 'Rscript'
    input:
    tuple val(chunk),file(zfile), val(ldfile), file(finemap),val(datasetID), path(datasetFile) from prep_finemap_locuszoom_input
    output:
    tuple val(chunk),file(credzoom) into prep_finemap_locuszoom_output_ch
    script:
    inputpath="${params.output}/${datasetID}/${chunk}/cred/"
    errorlog = "error_prep_finemap_locuszoom.log"
    credzoom = "${chunk}"+".cred.zoom"
    """
    Rscript $baseDir/bin/prep_finemap_locuszoom.R ${inputpath} ${params.method} ${params.nsignal} ${chunk} 2> $errorlog
    """ 
}

plink_input_ch2.join(prep_finemap_locuszoom_ch2).join(prep_finemap_locuszoom_output_ch).join(metal_ch).set{get_r2_from_single_SNP_input}

if(params.dprime){
  process get_dprime_from_single_SNP {
      scratch true
      label 'locuszoom'
      publishDir "${params.output}/${datasetID}/${chunk}/", mode: 'copy'
      input:
      each file(summarystats) from summarystats_ch2
      tuple val(chunk),file(zfile),file(snplist),first_gene_start,last_gene_end,CHR,left_cis_boundary,right_cis_boundary,plink_left_bound,plink_right_bound,subset_left_bound,subset_right_bound,n_genes,leadSNP_1, path(finemap) ,val(datasetID), path(datasetFile),file(credzoom),file(metal) from get_r2_from_single_SNP_input
      output:
      file('**/*')
      shell:
      rawld="raw${datasetID}.${chunk}.${leadSNP_1}.ld"
      chunkleadsnpld="${datasetID}.${chunk}.${leadSNP_1}.ld"
      outputdir="${chunk}_chr${CHR}_*"
      //outputdir="${chunk}_${CHR}_*"
      """
      plink --bfile ${params.reference} --extract ${snplist} --r2 inter-chr dprime --ld-snp ${leadSNP_1} --ld-window-r2 0 --a1-allele ${zfile} 4 1 --from-kb ${plink_left_bound} --to-kb ${plink_right_bound} --chr ${CHR} --out raw${datasetID}.${chunk}.${leadSNP_1} --threads ${task.cpus}
      cat $rawld | awk 'BEGIN{print "snp1 snp2 rsquare dprime"} NR!=1 {print \$6 " " \$3 " " \$7 " " \$8}' > $chunkleadsnpld
      /opt/locuszoom/locuszoom/bin/locuszoom --metal "${metal}" --refsnp "${leadSNP_1}" --chr ${CHR} --start ${subset_left_bound} --end ${subset_right_bound} --prefix ${chunk} --build hg38 --ld $chunkleadsnpld --db "${params.locuszoomdb}" fineMap="${credzoom}" showAnnot=T showRefsnpAnnot=T annotPch="24,24,25,25,22,22,8,7,21" ldCol="dprime" --no-date --ld-measure 'dprime'
    
      """ 
  }
}else{
  process get_r2_from_single_SNP {
      scratch true
      label 'locuszoom'
      publishDir "${params.output}/${datasetID}/${chunk}/", mode: 'copy'
      input:
      each file(summarystats) from summarystats_ch2
      tuple val(chunk),file(zfile),file(snplist),first_gene_start,last_gene_end,CHR,left_cis_boundary,right_cis_boundary,plink_left_bound,plink_right_bound,subset_left_bound,subset_right_bound,n_genes,leadSNP_1, path(finemap) ,val(datasetID), path(datasetFile),file(credzoom),file(metal) from get_r2_from_single_SNP_input
      output:
      file('**/*')
      shell:
      rawld="raw${datasetID}.${chunk}.${leadSNP_1}.ld"
      chunkleadsnpld="${datasetID}.${chunk}.${leadSNP_1}.ld"
      outputdir="${chunk}_chr${CHR}_*"
      //outputdir="${chunk}_${CHR}_*"
      """
      plink --bfile ${params.reference} --extract ${snplist} --r2 inter-chr --ld-snp ${leadSNP_1} --ld-window-r2 0 --a1-allele ${zfile} 4 1 --from-kb ${plink_left_bound} --to-kb ${plink_right_bound} --chr ${CHR} --out raw${datasetID}.${chunk}.${leadSNP_1} --threads ${task.cpus}
      cat $rawld | awk 'BEGIN{print "snp1 snp2 rsquare dprime"} NR!=1 {print \$6 " " \$3 " " \$7 " NA"}' > $chunkleadsnpld
      /opt/locuszoom/locuszoom/bin/locuszoom --metal "${metal}" --refsnp "${leadSNP_1}" --chr ${CHR} --start ${subset_left_bound} --end ${subset_right_bound} --prefix ${chunk} --gwas-cat whole-cat_significant-only --build hg38 --ld $chunkleadsnpld --db "${params.locuszoomdb}" fineMap="${credzoom}" showAnnot=T showRefsnpAnnot=T annotPch="24,24,25,25,22,22,8,7,21" --no-date 
    
      """ 
  }
}
workflow.onComplete { 
	log.info ( workflow.success ? "\nDone! Open the following directory for the outputs --> $params.output\n" : "Oops .. something went wrong" )
}
