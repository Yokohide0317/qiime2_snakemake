#!/bin/bash

set -eux

CONDA_ROOT="${HOME}/miniconda3"
ENV_NAME="qiime2-2022.8"
NUM_CORE=12


source ${CONDA_ROOT}/etc/profile.d/conda.sh
conda ${ENV_NAME}

snakemake --cores ${NUM_CORE}

conda deactivate
