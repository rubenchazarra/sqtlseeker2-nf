#!/bin/bash

module load nextflow/21.04.0
module load singularity/3.11.5
module load R/4.3.2


chr="21"
in_dir="/gpfs/projects/bsc83/Projects/scRNAseq/imestres/sQTLseeker/runs/data/mod/"
# out_dir="/gpfs/projects/bsc83/Projects/scRNAseq/imestres/sQTLseeker/runs/results_permuted_covs_pca/MONO/chr_${chr}/"
image="/gpfs/projects/bsc83/Projects/scRNAseq/imestres/sQTLseeker/runs/sqtlseeker2-nf.sif"
prog="/gpfs/projects/bsc83/Projects/scRNAseq/rchazarr/01_Projects/04_Single-cell-AS/02_Anquino2023/05_sQTL/03_sQTLseekeR-nf/sqtlseeker2-nf/sqtlseeker2.nf"

f_genotype="${in_dir}/genotype_chr${chr}_EUB.AFB.vcf.gz"
f_tcc="${in_dir}/MONO_tcc-exp-EUB.AFB.tsv.gz"
f_metadata="${in_dir}/MONO_metadata-EUB.AFB.tsv"
f_genes="${in_dir}/MONO_genes-EUB.AFB.bed"

mkdir -p $out_dir
cd $out_dir

nextflow run $prog \
        --genotype $f_genotype \
        --trexp $f_tcc \
        --metadata $f_metadata \
        --genes $f_genes \
        --dir $out_dir \
        --mode "permuted" \
        --covariates "true" \
        --svqtl "false" \
        --with-singularity $image \
        --resume

