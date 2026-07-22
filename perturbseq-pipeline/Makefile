# Author: Mohamed Kassam
# Date: 2026-07-19
install:
	micromamba create -f environment.yml
validate:
	python bin/validate_inputs.py --samplesheet example/samplesheet.csv --config conf/analysis.yaml --out /tmp/validated.csv --h5-list /tmp/h5.txt
test:
	pytest -q
run:
	nextflow run main.nf -profile conda --samplesheet example/samplesheet.csv --config conf/analysis.yaml
