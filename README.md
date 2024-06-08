# LeafDom
---

LeafDom is a command-line tool for preprocessing, performing PCA, clustering, and visualizing single-cell RNA sequencing (scRNA-seq) data stored in an `anndata` object. It supports flexible analysis by allowing users to choose different clustering methods (k-means and DBSCAN) and visualize the clusters in PCA plots. This project is inspired by the PCA clustering tools in `scanpy`.

We hope to implement more clustering methods in the near future such as `t-sne` and `UMAP`.

## Features

- **Preprocessing**:Standardize the scRNA-seq data.
- **PCA**: Perform Principal Component Analysis for dimensionality reduction.
- **Clustering**: Apply clustering algorithms (DBSCAN or k-means).
- **Visualization**: Visualize clusters in PCA plots and save the plots to specified output files.

## Requirements

- Python 3.x
- `numpy`
- `matplotlib`
- `anndata`
- `scikit-learn`

## Installation

1. **Clone the repository**:
   ```bash
   git clone https://github.com/Daraudom/LeafDom.git
   cd LeafDom
   ```

2. **Install the required Python packages**:
   ```bash
   pip install numpy matplotlib anndata scikit-learn
   ```
3. **Check to see if the tool is working**
   ```bash
   python leafdom.py --help
   ```
   Use `--help` to print out some helpful statements on how to operate Leafdom.

## Usage

### Basic Usage

```bash
python leafdom.py -i path_to_your_anndata.h5ad -o output_plot.png
```
### Command Line Arguments

#### Necessary Parameters

- `-i, --input`: Path to the input anndata file (h5ad format).
- `-o, --output_file`: Path to the output file for the PCA plot.

#### Optional Parameters

- `-m, --method`: Clustering method (`dbscan` or `kmeans`). Default is `kmeans`.
- `-k, --components`: Number of principal components to select for PCA. Default is `2`.
- `-e, --eps`: Maximum distance between two samples for DBSCAN. Default is `0.5`.
- `-s, --min_samples`: Number of samples in a neighborhood for a point to be considered as a core point in DBSCAN. Default is `10`.
- `-c, --clusters`: Number of clusters to form for k-means. Default is `3`.

### Clustering Methods

#### K-means

```bash
python leafdom.py -i path_to_your_anndata.h5ad -o output_plot.png -m kmeans -c <numer of clusters> -k <number of PCs>
```

- **Description**: K-means clustering partitions the data into `k` clusters by minimizing the variance within each cluster.
- **Advantages**: Simple and fast; works well when clusters are spherical and equally sized.
- **Disadvantages**: Requires the number of clusters (`k`) to be specified beforehand; sensitive to outliers; assumes clusters are spherical.

#### DBSCAN

```bash
python leafdom.py -i path_to_your_anndata.h5ad -o output_plot.png -m dbscan -eps <max distance between two samples> -s <min-samples in neighborhood> -k <number of PCs>
```

- **Description**: DBSCAN (Density-Based Spatial Clustering of Applications with Noise) clustering groups points that are closely packed together and marks as outliers points that lie alone in low-density regions.
- **Advantages**: Does not require the number of clusters to be specified; can find arbitrarily shaped clusters; robust to outliers.
- **Disadvantages**: Parameters `eps` (maximum distance) and `min_samples` (minimum number of points) need to be chosen carefully; may struggle with varying densities.

## Running the Tool

### Step-by-Step Instructions

1. **Prepare Your Data**:
   - **This is the important step!** Ensure your scRNA-seq data is stored in an `anndata` object and saved in h5ad format (this is the valid file extension of an anndata). The `anndata` object should have already been filtered with the quality control metrics. Additionally, it should have been normalized and logarithmically transformed. See `scanpy` tools, `sc.pp.normalize_per_cell()` and `sc.pp.log1p()` on how to do this. After this, ensure you have already identify the highly variable genes you're interested in. 
   
   The code blocks below is a summary on how to prepare your preprocessed `anndata` object. This is just an example and is not exhaustive of what you are limited to.

```python
# create a new python file
import anndata as ad
import scanpy as sc

# read your 10x cell ranger file (if only one dataset exists)
sc.read_10x_mtx(/path/to/10x_mtx_file, cache=True)
# OR
# read multiple 10x cell-ranger files
adatas = {} # prepare an empty dictionary
datasets = ["dataset 1", "dataset 2", "dataset 3"] # replace with the prefix of your datasets (e.g. GSM5114464_S7_D20)
DATADIR = "/path/to/directory_of_datasets" # save the directory of where your datasets is
for ds in datasets:
    print(f"Processing {ds}...") # log statement
    adatas[ds] = sc.read_10x_mtx(DATADIR, prefix=ds, cache=True)

combined = ad.concat(adatas, label="dataset") # combined the different adatas into one object using concat
combined.obs_names_make_unique()
```
```python
# Filter cells based on minimum counts and genes
sc.pp.filter_cells(combined, min_counts=...)
sc.pp.filter_cells(combined, min_genes=...)
# Filter genes based on minimum counts and cells
sc.pp.filter_genes(combined, min_counts=..)
sc.pp.filter_genes(combined, min_cells=..)

# Filtering based on QC metrics (e.g. Mitochondrial, Ribosomes,..)
combined.var["mt"] = combined.var_names.str.startswith("MT-") 
sc.pp.calculate_qc_metrics(
    combined, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True
)

# Extra filtering if needed...
adata_filt = combined[(combined.obs.pct_counts_mt < MT%) & ..., :] 
adata_filt
```

```python
# Data Normalization and Transformation
sc.pp.normalize_per_cell(adata_filt, counts_per_cell_after=1e4) # normalize to 10,000 reads/cell
sc.pp.log1p(adata_filt) # log transform

# Identify the Highly Variable Genes
sc.pp.highly_variable_genes(adata_filt, batch_key="dataset", n_top_genes=INTEGER)
final_adata = adata_filt[:, adata_filt.var["highly_variable"])]

# Save anndata in h5ad format
final_adata.write('final_adata.h5ad') # will automatically save to your current working diredctory
# you can move it to a different directory using `mv` in the terminal
```

2. **Run the tool from the command line** with the appropriate arguments.

```bash
python leafdom.py -i path_to_your_anndata.h5ad -o output_plot.png -m dbscan -k 2 -e 0.5 -s 10
```

## License

This project is licensed under the MIT License. See the [LICENSE](LICENSE) file for details.

---


functions include:
<!-- - preprocessing (normalization and log transformation)
- selecting highly variable genes -> new adata object with only the highest variable genes (select top n genes) -->
- performing PCA on the new data
- visualize the PCA plots using either k-means or DBscans
