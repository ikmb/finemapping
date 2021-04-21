# finemapping
## A Nextflow framework for finemap

### Installation:
Download this pipeline locally or on the medcluster.
Install finemap_v1.4_x86_64.tgz within the bin folder so finemap is located in bin/finemap_v1.4_x86_64/finemap_v1.4_x86_64.

Pipeline needs 6 input files:
**Reference**: in bim, bed, fam (3 files with same basename)
**Locus-file**: csv-file setting the boundaries of the finemap plots; columns: chunk,NSNP,chr,st,sp,PPA_3
**SNP-List**: file with 1 snp per row

### Usage:
Call pipeline with:
```bash
nextflow run main.nf -profile local/medcluster --locus /home/user/finepipe/example/locusfile.sample --snps /home/user/finepipe/example/snplist.sample --reference /home/user/finepipe/example/GerNorItaSpa.chr3 --sumstats /home/user/finepipe/example/sumstats.sample --nsum 15743 --nsignal 1 --ld r2 --method sss -resume  
```
#### Mandatory:
**-profile**    "local" or "medcluster"
**--locus** "path/to/locus.file"
**--snps**  "path/to/snp.list"
**reference**   "path/to/bimbedfam" (no file extension)
**--sumstats**  "path/to/sumstat.file"
**--nsum**  N of studysize
**--ld**    "r2" (not implemented yet: "dprime")
**--method**    "sss" or "cond"

If you run it locally, download the locuszoom database and set it with:
**--locuszoomdb** "/path/to/locuszoom_hg38.db"

#### Optional:
**-resume** Continue a run with cached processes


