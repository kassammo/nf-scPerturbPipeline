#!/usr/bin/env nextflow
// Author: Mohamed Kassam
// Date: 2026-07-19
nextflow.enable.dsl=2
include { VALIDATE_INPUTS } from './modules/validate_inputs'
include { FASTQC } from './modules/fastqc'
include { MULTIQC } from './modules/multiqc'
include { CELLRANGER_COUNT } from './modules/cellranger_count'
include { PERTURBSEQ_ANALYSIS } from './modules/perturbseq_analysis'

workflow {
  if (!params.samplesheet) error "--samplesheet is required"
  validated = VALIDATE_INPUTS(file(params.samplesheet), file(params.config))
  if (params.run_cellranger) {
    if (!params.transcriptome) error "--transcriptome is required with --run_cellranger"
    reads = Channel.fromPath(params.input ?: 'NO_FASTQS_FOUND', checkIfExists:true)
    FASTQC(reads)
    MULTIQC(FASTQC.out.reports.collect())
    CELLRANGER_COUNT(validated.out.samplesheet, reads.collect(), file(params.transcriptome), params.guide_reference ? file(params.guide_reference) : [])
    h5_ch = CELLRANGER_COUNT.out.h5
  } else {
    h5_ch = validated.out.h5_inputs
  }
  PERTURBSEQ_ANALYSIS(h5_ch.collect(), validated.out.samplesheet, validated.out.config)
}
