# ![logo](/docs/images/annoyance.logo.svg)

### [**annoy**](https://github.com/spotify/annoy) for the **an**notation of **ce**ll types. 


This is a simple API for classifying cell types as annotated in an **`adata.obs`** table. An example notebook may be found [**here**](docs/notebooks/annoyance.example.pbmc3k.ipynb).

### Overview
```python
import annoyance
import anndata as a

# load data
adata = a.read_h5ad("../inputs/weinreb2020.h5ad")

# instantiate model
annoy = annoyance.SpotifyAnnoy()
```

**Step 1**. split data on categorical labels
```python
annoy.split_data(
    adata, annot_col="Annotation", labels=["Monocyte", "Neutrophil"], on="X_pca"
)
```

**Step 2**. Build and evaluate the predictive model based on the data split

```python
annoy.prebuild()
```

**Step 3**. Build the classifier using all of the data then save the `AnnoyIndex`. 
```python
annoy.build()
annoy.save()
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

This project uses open-source code from [**spotify/annoy**](https://github.com/spotify/annoy). However, **this repo is in no way affiliated with Spotify**. Questions and interests may be directed to [**Michael Vinyard**](https://github.com/mvinyard) at: [mvinyard@broadinstitute.org](mailto:mvinyard@broadinstitute.org). 
