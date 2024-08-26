# STCAT  <a href="https://www.python.org/"><img src="https://img.shields.io/badge/python-3.8+-brightgreen.svg" alt="Python Versions" width="80"></a>

STCAT is an automated T cell type annotation tool for scRNA-seq datasets. It is based on CellTypist, a logistic regression classifier optimized by the stochastic gradient descent algorithm. We made some changes to CellTypist and developed this tool in conjunction with our marker. Models trained by CellTypist are used in our automated annotation process.

# TcellAtlats website
Information of STCAT can be also found in our TcellAtlas portal. 
<a href="https://github.com/GuoBioinfoLab/STCAT"><img src="[https://img.shields.io/badge/TcellAtlas-blue](https://guolab.wchscu.cn/TCellAtlas/#/)" alt="TcellAtlas website" width="50"></a>

# Install STCAT
### Using pip
```console
pip install STCAT
```

# Usage 
STCAT expects to use an Anndata object ( .h5ad file format ) as input, and at the same time, a raw count matrix ( reads or UMIs ) is required. The file input is in a cell-by-gene format ( cells as rows and genes as columns ).
```python
import STCAT
results = STCAT.STCAT(<your_adata>)
```

# Citation


