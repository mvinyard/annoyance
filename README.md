# ![logo](/docs/images/annoyance.logo.svg)

### [**annoy**](https://github.com/spotify/annoy) for the **an**notation of **ce**ll types. 


This is a simple API for classifying cell types as annotated in an **`adata.obs`** table. An example notebook may be found [**here**](docs/notebooks/annoyance.example.pbmc3k.ipynb).

### Overview
```python
import annoyance
import anndata as a

# load data
adata = a.read_h5ad("./data/pbmc3k.processed.h5ad")
```

**Step 1**. Instantiate the model and split the data on categorical labels
```python
annoy = annoyance.SpotifyAnnoy()
annoy.split_data(
    adata, annot_col="Annotation", labels=["Monocyte", "Neutrophil"], on="X_pca"
)
```
```
label mask created:

CD4 T cells          1144
CD14+ Monocytes       480
B cells               342
CD8 T cells           316
NK cells              154
FCGR3A+ Monocytes     150
Dendritic cells        37
Megakaryocytes         15
Name: louvain, dtype: int64
```

**Step 2**. Build and evaluate the predictive model based on the data split

```python
annoy.prebuild()
```
```
                   precision    recall  f1-score   support

          B cells       1.00      1.00      1.00        86
  CD14+ Monocytes       0.94      0.99      0.97       120
      CD4 T cells       0.98      0.98      0.98       286
      CD8 T cells       0.91      0.91      0.91        79
  Dendritic cells       1.00      0.67      0.80         9
FCGR3A+ Monocytes       0.95      0.95      0.95        37
   Megakaryocytes       1.00      0.75      0.86         4
         NK cells       0.97      0.95      0.96        39

         accuracy                           0.97       660
        macro avg       0.97      0.90      0.93       660
     weighted avg       0.97      0.97      0.96       660
```

**Step 3**. Build the classifier using all of the data then save the `AnnoyIndex`. 
```python
annoy.build()
annoy.save()
```
```
Annoyance annoy index saved to:

	annoyance_index/annoyance.annoy_idx.50_dims.20_neighbors.10_trees.ann
	annoyance_index/annoyance.annoy_idx.50_dims.20_neighbors.10_trees.txt
```

### Installation

```BASH
pip install annoyance
```

#### To install the development version

```BASH
git clone https://github.com/mvinyard/annoyance.git

cd ./annoyance/
pip install -e .
```

### Notes

This project uses open-source code from [**spotify/annoy**](https://github.com/spotify/annoy). However, **this repo is in no way affiliated with Spotify**. 

Interested? Questions and discussion may be directed to [**Michael Vinyard**](https://github.com/mvinyard) at: [mvinyard@broadinstitute.org](mailto:mvinyard@broadinstitute.org). 
