// Author: Mohamed Kassam
// Date: 2026-07-19
process CELLRANGER_COUNT {
 tag { sample.sample_id }
 publishDir "${params.outdir}/cellranger", mode:'symlink'
 input:
 path samplesheet
 path reads
 path transcriptome
 path guide_reference
 output:
 tuple val(sample), path("${sample.sample_id}/outs/filtered_feature_bc_matrix.h5"), emit:h5
 script:
 sample=[sample_id:'sample']
 """
 ${params.cellranger_path} count --id=${sample.sample_id} --transcriptome=$transcriptome --fastqs=. --sample=${sample.sample_id} --localcores=${task.cpus} --localmem=60 ${guide_reference ? "--feature-ref=$guide_reference" : ''}
 """
}
