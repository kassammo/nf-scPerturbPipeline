// Author: Mohamed Kassam
// Date: 2026-07-19
process FASTQC {
 tag { reads.simpleName }
 publishDir "${params.outdir}/fastqc", mode:'copy'
 input: path reads
 output: path '*_fastqc.{html,zip}', emit:reports
 script: "fastqc --threads ${task.cpus} $reads"
}
