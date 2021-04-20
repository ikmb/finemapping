# finemapping
A Nextflow framework for finemap

Installation:
Download this pipeline locally or on the medcluster.
Install finemap_v1.4_x86_64.tgz within the bin folder so finemap is located in bin/finemap_v1.4_x86_64/finemap_v1.4_x86_64.

Pipeline needs 6 input files:
reference: in bim, bed, fam (3 files with same basename)
locusfile: csv-file setting the boundaries of the finemap plots; columns: chunk,NSNP,chr,st,sp,PPA_3
snplist: file with 1 snp per row

Call pipeline with:
nextflow run main.nf 
-profile local/medcluster 
--boundaries /home/user/finepipe/example/locusfile.sample 
--snps /home/user/finepipe/example/snplist.sample 
--reference /home/user/finepipe/example/GerNorItaSpa.chr3 
--sumstats /home/user/finepipe/example/sumstats.sample 
--nsum 15743 
--nsignal 1 
--ld r2 
--method sss 
-resume  

If locally, download the locuszoom database and set it with:
--locuszoomdb /home/user/locuszoom/data/database/locuszoom_hg38.db
