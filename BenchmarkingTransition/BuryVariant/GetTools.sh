#!/usr/bin/env bash

set -ex

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
OPT_DIR=${DIR}/opt

mkdir -p ${OPT_DIR}
pushd ${OPT_DIR}

PYTHON_DIR=${OPT_DIR}/miniconda3
CONDA=Miniconda3-latest-Linux-x86_64.sh

if [[ ! -d ${PYTHON_DIR} ]]; then
wget -q https://repo.continuum.io/miniconda/${CONDA}\
    && sh ${CONDA} -b -p ${PYTHON_DIR}\
    && ${PYTHON_DIR}/bin/python ${PYTHON_DIR}/bin/pip install pyvcf==0.6.8\
    && ${PYTHON_DIR}/bin/python ${PYTHON_DIR}/bin/pip install scipy==1.1.0\
    && ${PYTHON_DIR}/bin/python ${PYTHON_DIR}/bin/conda install --yes -c bioconda samtools=${samtools_version} htslib=1.9 bcftools=1.9 \
    && rm -f ${CONDA}
fi


