# sQTLseekeR2 Modifications

This file is just for me to keep track of all the changes made in the original
pipeline.

## TODOs:

- [x] Make pipeline full custom
- [x] Change `store_true` to `store` in R scripts
- [x] Use the MOD version of `sqtls.p.R` in `sqtlseeker2.nf`
- [] Add info of new options 
- [x] Use the data in the `RUNS/` directory 
- [x] Set the correct default parameters

## Make Pipeline FULL custom

The parameters to add from the scripts `prepare_trexp.R`, `sqtlseeker.R`, 
`sqtlseeker.p.R` `sqtls.R` and `sqtls.p.R`.

> Note that the hard-coded values are taken from the pipeline not
> from the functions in the R package.

### prepare_trexp.R

To add:
    - `min_gene_exp` hard-coded to 1; set to 0.1
    - `min_proportion` hard-coded to 0.8; set to 0.4 (this proportion is per
    GENE, not per transcript)
    - `min_transcript_expr` hard-coded to 0.1; set to 0.1
    - `min_dispersion` hard-coded to 0.1; set to 0.01

### sqtlseeker.R

To add:
    - `asympt` hard-coded as a flag; set to true/false
    - `nb_perm_max` hard-coded to 1000; set to 1000
    - `min_nb_ext_scores` hard-coded to 1000; set to 1000
    - `nb_perm_max` hard-coded to 1e6; set to 1e6
    - `nb_perm_max_svqtl` hard-coded to 1e4; set to 1e4
    - `min_nb_ind_geno` hard-coded to 10; set to 10
    - `ld` hard-coded to 0; set to 0.1
    - `svqtl` hard-coded as a flag; set as a flag

### sqtlseeker.p.R

To add:
    - `min_nb_ext_scores` hard-coded to 100; set to 100
    - `nb_perm_min` hard-coded to 100; set to 100
    - `min_nb_ind_geno` hard-coded to 10; set to 10

### sqtls.R and sqtls.p.R

To add:
    - `fdr` hard-coded to 0.05; set to 0.05
    - `fdr_svqtl` hard-coded to 0.05; set to 0.05
    - `rm_svqtl` hard-coded as a flag; set to true/false (?)
    - `md_min` hard-coded to 0.05; set to 0
    - `type_fdr` hard-coded to BH: set to qvalue
    - `plot_pdf` hard-coded to NULL; set to NULL

