// Author: Mohamed Kassam
// Date: 2026-07-19
process VALIDATE_INPUTS {
 tag "validate"
 publishDir "${params.outdir}/validation", mode:'copy'
 input: path samplesheet; path config
 output:
 path 'validated_samplesheet.csv', emit:samplesheet
 path 'analysis.yaml', emit:config
 path 'h5_inputs.txt', emit:h5_inputs
 script:
 """
 python ${projectDir}/bin/validate_inputs.py --samplesheet $samplesheet --config $config --out validated_samplesheet.csv --h5-list h5_inputs.txt
 cp $config analysis.yaml
 """
}
