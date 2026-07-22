# Author: Mohamed Kassam
# Date: 2026-07-19
FROM mambaorg/micromamba:2.0.5
COPY --chown=$MAMBA_USER:$MAMBA_USER environment.yml /tmp/environment.yml
RUN micromamba install -y -n base -f /tmp/environment.yml && micromamba clean --all --yes
COPY . /opt/perturbseq
WORKDIR /opt/perturbseq
ENTRYPOINT ["/usr/local/bin/_entrypoint.sh"]
CMD ["nextflow", "run", "main.nf", "--help"]
