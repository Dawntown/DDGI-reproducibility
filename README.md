# DDGI Reproducibility

This repository contains figure reproduction and analysis scripts for [DDGI](https://github.com/Dawntown/DDGI).

The notebooks are organized by analysis section:

- `1_preprocessing/`: dataset summary and preprocessing overview
- `2_benchmarking/`: benchmarking analyses
- `3_tcells/`: T cell analysis
- `4_macrophages/`: macrophage analysis
- `5_crispr/`: CRISPR perturbation analyses
- `tutorial/`: example DDGI workflow
- `datasets/`: dataset preprocessing and data-generation notebooks

To run the notebooks and scripts, use a Python environment where DDGI and its dependencies are installed. The analysis notebooks expect preprocessed datasets and generated model outputs to be available at the relative paths used in the notebooks and configuration files.

For `5_crispr/5.6_run.sh`, set `DDGI_PYTHON` to the Python executable from the environment where DDGI is installed, for example:

```bash
DDGI_PYTHON=/path/to/ddgi/bin/python bash 5_crispr/5.6_run.sh 100
```
