# finemapping
## A Nextflow framework for finemap

### Installation:
Download this pipeline locally or on the medcluster.<br />
Install finemap_v1.4_x86_64.tgz within the bin folder so finemap is located in bin/finemap_v1.4_x86_64/finemap_v1.4_x86_64.<br />

Pipeline needs 6 input files:<br />
**Reference**: in bim, bed, fam (3 files with same basename)<br />
**Locus-file**: csv-file setting the boundaries of the finemap plots; columns: chunk,NSNP,chr,st,sp,PPA_3<br />
**SNP-List**: file with 1 snp per row<br />

### Usage:
Call pipeline with:<br />
```bash
nextflow run main.nf -profile standard --locus /home/user/finepipe/example/locusfile.sample --snps /home/user/finepipe/example/snplist.sample --reference /home/user/finepipe/example/GerNorItaSpa.chr3 --sumstats /home/user/finepipe/example/sumstats.sample --nsum 15743 --nsignal 1 --method sss -resume  
```
#### Mandatory:
**--locus** "path/to/locus.file"<br />
**--snps**  "path/to/snp.list"<br />
**--reference**   "path/to/bimbedfam" (no file extension)<br />
**--sumstats**  "path/to/sumstat.file"<br />
**--nsum**  N of studysize<br />
**--method**    "sss" or "cond"<br />
**--nsignal**   N of max signals<br />

If you run it locally, download the locuszoom database and set it with:<br />
**--locuszoomdb** "/path/to/locuszoom_hg38.db"<br />
and set profile to local:
**-profile**    "local" or "standard"<br />

#### Optional:
**--dprime**    sets ld method from the default r² to dprime<br />
**--output** "/path/to/output" if not set output is baseDir/Results<br />
**-resume** Continue a run with cached processes<br />

