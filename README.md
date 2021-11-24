# annoyance

[**annoy**](https://github.com/spotify/annoy) for the **an**notation of **ce**ll types. 

This is a simple API for classifying cell types as annotated in an **`adata.obs`** table. 

## Example usage

```python
import annoyance as oy
import anndata as a

# read AnnData
adata = a.read_h5ad("/path/to/data.h5ad")

# build and run classifier
oy.build(adata)

adata

...
```

## Installation

```python
pip install annoyance
```

### To install the development version

```python
git clone https://github.com/mvinyard/annoyance.git

cd ./annoyance/
pip install -e .
```
