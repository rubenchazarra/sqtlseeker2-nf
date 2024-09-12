/*
 * Copyright (c) 2019, Centre for Genomic Regulation (CRG)
 *
 * Copyright (c) 2019, Diego Garrido-Martín
 * 
 * This file is part of 'sqtlseeker2-nf': 
 * sQTLseekeR2 in Nextflow, a pipeline for splicing QTL mapping
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

// SEP-2024 Edits to the pipeline made by Ruben Chazarra-Gil & Iris Mestres at
// Barcelona Supercomputing Center

// Define parameters

params.genotype = null
params.trexp = null
params.metadata = null
params.genes = null
params.dir = "result"
params.mode = "nominal"
params.win = 5000
params.covariates = false
params.kn = 10
params.kp = 1
params.fdr = 0.05
params.svqtl = false
params.ld = 0
params.min_md = 0.05
params.max_perm = 1000 
params.help = false

// Define full-custom additional parameters
// prepare_trexp

params.min_gene_exp = 1
params.min_proportion = 0.8
params.min_transcript_expr = 0.1
params.min_dispersion = 0.1

// sqtlseeker and sqtlseeker.p

params.asympt = true
params.min_nb_ext_scores_nominal = 1000
params.min_nb_ext_scores_permuted = 100
params.nb_perm_max_nominal = 1e6
params.nb_perm_max_svqtl = 1e4
params.nb_perm_min = 100
params.min_nb_ind_geno = 10

// sqtls and sqtls.p

params.fdr_svqtl = 0.05
params.type_fdr = "BH"
params.plot_pdf = null


// Print usage

if (params.help) {
  log.info ''
  log.info 'sqtlseeker2-nf ~ A pipeline for splicing QTL mapping'
  log.info '----------------------------------------------------'
  log.info 'Run sQTLseekeR2 on a set of data.'
  log.info ''
  log.info 'Usage: '
  log.info "    ${workflow.projectDir.baseName} [options]"
  log.info ''
  log.info 'Options:'
  log.info '--genotype GENOTYPE_FILE    the genotype file'
  log.info '--trexp EXPRESSION_FILE     the transcript expression file'
  log.info '--metadata METADATA_FILE    the metadata file'
  log.info '--genes GENES_FILE          the gene location file' 
  log.info "--dir DIRECTORY             the output directory (default: $params.dir)"
  log.info "--mode MODE                 the run mode: nominal or permuted (default: $params.mode)"
  log.info "--win WINDOW		the cis window in bp (default: $params.win)"
  log.info "--covariates COVARIATES     include covariates in the model (default: $params.covariates)"
  log.info "--fdr FDR                   false discovery rate level (default: $params.fdr)"
  log.info "--min_md MIN_MD             minimum effect size reported (default: $params.min_md)"
  log.info "--svqtl SVQTLS              test for svQTLs (default: $params.svqtl)"
  log.info ''
  log.info 'Additional parameters for mode = nominal:'
  log.info "--ld LD                     threshold for LD-based variant clustering (default: $params.ld, no clustering)"
  log.info "--kn KN                     number of genes per batch in nominal pass (default: $params.kn)"
  log.info ''
  log.info 'Additional parameters for mode = permuted:'
  log.info "--kp KP                     number of genes per batch in permuted pass (default: $params.kp)"
  log.info "--max_perm MAX_PERM         maximum number of permutations (default: $params.max_perm)"
  log.info ''
  exit 1
}


// Check mandatory options

if (!params.genotype) {
    exit 1, "Genotype file not specified."
} else if (!params.trexp){
    exit 1, "Transcript expression file not specified."
} else if (!params.metadata){
    exit 1, "Metadata file not specified."
} else if (!params.genes){
    exit 1, "Gene location file not specified."
}

 
// Print selected options

log.info ""
log.info "sqtlseeker2-nf ~ A pipeline for splicing QTL mapping"
log.info ""
log.info "General parameters"
log.info '------------------'
log.info "Genotype file                      : ${params.genotype}"
log.info "Transcript expression file         : ${params.trexp}"
log.info "Metadata file                      : ${params.metadata}"
log.info "Gene location file                 : ${params.genes}"
log.info "Output directory                   : ${params.dir}"
log.info "Run mode                           : ${params.mode}"
log.info "Cis window                         : ${params.win}"
log.info "Covariates                         : ${params.covariates}"
log.info "FDR level                          : ${params.fdr}"
log.info "Min. effect size                   : ${params.min_md}"
log.info "Test for svQTLs                    : ${params.svqtl}"
log.info ""

if(params.mode == "nominal"){
  log.info 'Additional parameters for mode = nominal'
  log.info '----------------------------------------'
  log.info "LD-based clustering threshold      : ${params.ld}"
  log.info "Genes/batch in nominal pass        : ${params.kn}"
  log.info ""
} else if(params.mode == "permuted"){
  log.info 'Additional parameters for mode = permuted'
  log.info '-----------------------------------------'
  log.info "Genes/batch in permuted pass       : ${params.kp}"
  log.info "Max. number of permutations        : ${params.max_perm}"
  log.info ""
}


// Create file objects given parameters

genotype_file = file(params.genotype)
trexp_file = file(params.trexp)
metadata_file = file(params.metadata)
genes_file = file(params.genes)


// Obtain the 'groups' list

def groups = []
myReader = metadata_file.newReader()
String line
while( line = myReader.readLine() ) {
    def (sampleId, indId, group, covariates) = line.tokenize('\t')
    if( group != 'group' ) {
    	groups += group
    }
}
myReader.close()
groups.unique()
String show = groups.join(", ")

println "groups: $show\n"


// Index genotype file

process index {

    publishDir "${params.dir}/00_idx_genotype", mode: 'copy'

    input:
    file genotype from genotype_file

    output:    
    set file("${genotype.baseName}.*bgz"), file("${genotype.baseName}.*bgz.tbi") into index_ch
    
    script:
    """
    index_geno.R --genotype_file $genotype 

    """
}

index_ch.into{index2nominal_ch; index2permuted_ch}


// Preprocess input data

process prepare {

    publishDir "${params.dir}/01_tx_exp/$group", mode: 'copy'

    tag { group }

    input:
    val group from Channel.from(groups)
    file te from trexp_file 
    file metadata from metadata_file
    file genes from genes_file

    output:
    set val(group), file('tre.df.RData') into tre_ch
    set val(group), file('genes.ss.bed') into genes_ch
    set val(group), file('covariates.df.RData') into cov_ch
    
    script:
    if (params.covariates == true)
    """
    prepare_trexp.R \
        --transcript_expr $te \
        --metadata $metadata \
        --gene_location $genes \
        --group "$group" \
        --covariates \
        --min_gene_expr ${params.min_gene_exp} \
        --min_proportion ${params.min_proportion} \
        --min_transcript_expr ${params.min_transcript_expr} \
        --min_dispersion ${params.min_dispersion} \
        --output_tre tre.df.RData \
        --output_gene genes.ss.bed \
        --output_cov covariates.df.RData

    """
    else
    """
    prepare_trexp.R \
        --transcript_expr $te \
        --metadata $metadata \
        --gene_location $genes \
        --group "$group" \
        --min_gene_expr ${params.min_gene_exp} \
        --min_proportion ${params.min_proportion} \
        --min_transcript_expr ${params.min_transcript_expr} \
        --min_dispersion ${params.min_dispersion} \
        --output_tre tre.df.RData \
        --output_gene genes.ss.bed \
        --output_cov covariates.df.RData

    """
}

tre_ch.into {tre2nominal_ch; tre2permuted_ch}
genes_ch.into {genes2nominal_ch; genes2permuted_ch} 
cov_ch.into {cov2nominal_ch; cov2permuted_ch}


// Run sQTLseekeR2 (nominal)

tre2nominal_ch.join(cov2nominal_ch).combine(genes2nominal_ch.splitText( by: params.kn, file: "nominal_in" ), by: 0).set{nominal_in_ch}

process nominal_test {

    publishDir "${params.dir}/02_sQTL_nominal/$group", mode: 'copy'

    tag {"$group, $chunk"}
 
    input:
    set val(group), file(tre_rdata), file(cov_rdata), file (chunk) from nominal_in_ch
    set file(indexed_geno), file(tbi) from index2nominal_ch 

    output:
    set val(group), file('nominal_out.*') into nominal_out_ch 
 
    script: 
    """
    res=\$(echo $chunk | sed 's/_in/_out/')
    sqtlseeker.R \
        --transcript_expr $tre_rdata \
        --indexed_geno $indexed_geno \
        --gene_location $chunk \
        --covariates $cov_rdata \
        --asympt ${params.asympt} \
        --min_nb_ext_scores ${params.min_nb_ext_scores_nominal} \
        --min_nb_ind_geno ${params.min_nb_ind_geno} \
        --nb_perm_max ${params.nb_perm_max_nominal} \
        --nb_perm_max_svqtl ${params.nb_perm_max_svqtl} \
        --ld ${params.ld} \
        --window ${params.win} \
        --svqtl ${params.svqtl} \
        --output_file \$res 

    """
}

nominal_out_ch.collectFile(sort: { it.name }).map() {[it.name, it]}.into{all_nominal_tests_ch1; all_nominal_tests_ch2}

process nominal_mtc {

    // publishDir "${params.dir}/groups/$group"
    publishDir "${params.dir}/02_sQTL_nominal/$group", mode: 'copy'

    tag { group }

    input: 
    set val(group), file('all-tests.nominal.tsv') from all_nominal_tests_ch1

    output:
    set val(group), file('all-tests.nominal.tsv'), file ("sqtls-${params.fdr}fdr.nominal.tsv") into nominal_end_ch
 
    script:
    """
    sqtls.R \
        --nominal all-tests.nominal.tsv \
        --fdr ${params.fdr} \
        --rm_svqtl ${params.svqtl} \
        --md_min ${params.min_md} \
        --type_fdr ${params.type_fdr} \
        --plot_pdf ${params.plot_pdf} \
        -output sqtls-${params.fdr}fdr.nominal.tsv

    """       
}


// Run sQTLseekeR2 (permuted)
 
if (params.mode == "permuted") {

    tre2permuted_ch.join(cov2permuted_ch).combine(genes2permuted_ch.splitText( by: params.kp, file: "permuted_in" ), by: 0).set{permuted_in_ch}

    process permuted_test {

        publishDir "${params.dir}/03_sQTL_permuted/$group", mode: 'copy'

        tag {"$group, $chunk"}

        input:
	    set val(group), file(tre_rdata), file(cov_rdata), file (chunk) from permuted_in_ch
        set file(indexed_geno), file(tbi) from index2permuted_ch

        output:
        set val(group), file('permuted_out.*') into permuted_out_ch

        script:
        """
        res=\$(echo $chunk | sed 's/_in/_out/')
        sqtlseeker.p.R -t $tre_rdata -i $indexed_geno \
            --transcript_expr $tre_rdata \
            --indexed_geno $indexed_geno \
            --gene_location $chunk \
            --covariates $cov_rdata \
            --min_nb_ext_scores ${params.min_nb_ext_scores_permuted} \
            --min_nb_ind_geno ${params.min_nb_ind_geno} \
            --nb_perm_max ${params.max_perm} \
            --nb_perm_min ${params.nb_perm_min} \
            --window ${params.win} \
            --output_file \$res 
        """
    }    

    permuted_out_ch.collectFile(sort: { it.name }).map() {[it.name, it]}.set{all_permuted_tests_ch}

    all_nominal_tests_ch2.join(all_permuted_tests_ch).set{all_tests_ch}

    process permuted_mtc {

        // publishDir "${params.dir}/groups/$group"
        publishDir "${params.dir}/03_sQTL_permuted/$group", mode: 'copy'
        
        tag { group }

        input:
        
        set val(group), file('all-tests.nominal.tsv'), file('all-tests.permuted.tsv') from all_tests_ch

        output:
        set val(group), file('all-tests.permuted.tsv'), file ("sqtls-${params.fdr}fdr.permuted.tsv") into permuted_end_ch
 
        script:
        """
        sqtls.p.MOD.R \
            --nominal all-tests.nominal.tsv \
            --permuted all-tests.permuted.tsv \
            --fdr ${params.fdr} \
            --type_fdr ${params.type_fdr} \
            --rm_svqtl ${params.svqtl} \
            --md_min ${params.min_md} \
            --output sqtls-${params.fdr}fdr.permuted.tsv
        """
    }
}

