# Geneset Enrichment Analysis

This repository contains a script [gsea.py], that runs a geneset enrichment analysis
for genesets described in terms of *M. tuberculosis H37Rv* locus tags. It uses the
[COMBAT TB NeoDb](https://neodb.sanbi.ac.za/browser/).

Requirements:

```txt
click
statsmodels
scipy
neo4j_python_driver
```

and Python > 3.6. To install these requirements with conda, run `conda create --name gsea --file environment.txt` where
`gsea` is the name of the conda environment you are creating (you can change it if you like). Then activate the environment with `conda activate gsea` before running the `gsea.py` script.
