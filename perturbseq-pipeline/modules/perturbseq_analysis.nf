// Author: Mohamed Kassam
// Date: 2026-07-19
process PERTURBSEQ_ANALYSIS {
 tag 'analysis'
 publishDir "${params.outdir}/analysis", mode:'copy'
 input: path h5_list; path samplesheet; path config
 output: path 'perturbseq_results/**', emit:results
 script:
 """
 mkdir -p perturbseq_results
 Rscript ${projectDir}/R/run_pipeline.R --samplesheet $samplesheet --config $config --outdir perturbseq_results
 """
}
