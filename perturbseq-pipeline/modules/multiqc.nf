// Author: Mohamed Kassam
// Date: 2026-07-19
process MULTIQC {
 publishDir "${params.outdir}/multiqc", mode:'copy'
 input: path reports
 output: path 'multiqc_report.html'
 script: "multiqc ."
}
