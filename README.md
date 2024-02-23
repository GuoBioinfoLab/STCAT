# TCAT

TCAT is an automated T cell type annotation tool for scRNA-seq datasets. It is based on CellTypist, a logistic regression classifier optimized by the stochastic gradient descent algorithm. We made some changes to CellTypist and developed this tool in conjunction with our marker. Models trained by CellTypist are used in our automated annotation process.[Python Versions](https://img.shields.io/badge/python-3.8+-brightgreen.svg)

# Install TCAT
### Using pip 
```console
pip install TCAT
```

# Usage 
TCAT expects to use an Anndata object (.h5ad file format) as input, and at the same time, a raw count matrix (reads or UMIs) is required. The file input is in a cell-by-gene format (cells as rows and genes as columns).
```python
import TCAT
results = TCAT.TCAT(<your_adata>)
```

# Citation
