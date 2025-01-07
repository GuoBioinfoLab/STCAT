# STCAT  <a href="https://www.python.org/"><img src="https://img.shields.io/badge/python-3.8+-brightgreen.svg" alt="Python Versions" width="80"></a>

STCAT is an automated T cell type annotation tool for scRNA-seq datasets. 
It based on a high-confidence T cell subtypes and states reference. 
The reference can be found in our TCellAtlas portal. 
STCAT can automatically annotate T cell subtypes and states for scRNA-seq data of different conditions and tissues.

# TCellAtlats Website
TCellAtlas contains 1,677,799 high-quality T cells of 339 samples from 38 10x Genomics projects across 37 conditions and 16 tissues. It also includes 47,215 high-quality T cells in 21 conditions and 8 tissues from 18 Smart-seq projects. TCellAtlas contains all 68 T cell subtypes/states, which makes it the most comprehensive T cell subtypes/states and T cell database with the largest number of cells.
Information of STCAT can be also found in our TCellAtlas portal. 
The database is accessible at [TCellAtlas](https://guolab.wchscu.cn/TCellAtlas/#/).
# Install STCAT
## Create environment

```
conda create -n STCAT python=3.9.16 pandas=2.2.0
conda activate STCAT
```
## Using pip
```console
pip install STCAT
```
# Usage 
STCAT expects a raw count matrix as input and can be implemented with only one line of code in Python. 
STCAT expects to use an Anndata object ( .h5ad file format ) as input, and at the same time, a raw count matrix ( reads or UMIs ) is required. 
The file input is in a cell-by-gene format ( cells as rows and genes as columns ). For more information, please see [anndata](https://anndata.readthedocs.io/en/latest/).
The barcode should be unique for each cell, with no duplicates.
As for the annotation result, STCAT will be automatically added to the common anndata format of scRNA-seq analysis for easy viewing.
```python
import scanpy as sc
import STCAT
adata = sc.read_h5ad(<file_path>)
results = STCAT.STCAT(adata)
```
## Example:
Here is an example for guidance, and the demo.h5ad file mentioned in the example can be found below.
[Tutorial](tutorial.ipynb)
### demo.h5ad file in Tutorial
[demo.h5ad](demo.h5ad.bz2)
