# LeafDom
---

LeafDom is a command-line tool for preprocessing, performing PCA, clustering, and visualizing single-cell RNA sequencing (scRNA-seq) data stored in an `anndata` object. It supports flexible analysis by allowing users to choose different clustering methods and visualize the clusters in PCA plots.

## Features

- **Preprocessing**: Normalize, log-transform, and standardize scRNA-seq data.
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
   git clone https://github.com/yourusername/LeafDom.git
   cd LeafDom
   ```

2. **Install the required Python packages**:
   ```bash
   pip install numpy matplotlib anndata scikit-learn
   ```

## Usage

### Command Line Arguments

#### Necessary Parameters

- `-i, --input`: Path to the input anndata file (h5ad format).
- `-o, --output_file`: Path to the output file for the PCA plot.

#### Optional Parameters

- `-m, --method`: Clustering method (`dbscan` or `kmeans`). Default is `dbscan`.
- `-k, --components`: Number of principal components to select for PCA. Default is `2`.
- `-e, --eps`: Maximum distance between two samples for DBSCAN. Default is `0.5`.
- `-s, --min_samples`: Number of samples in a neighborhood for a point to be considered as a core point in DBSCAN. Default is `10`.
- `-c, --clusters`: Number of clusters to form for k-means. Default is `3`.

### Clustering Methods

#### K-means

- **Description**: K-means clustering partitions the data into `k` clusters by minimizing the variance within each cluster.
- **Advantages**: Simple and fast; works well when clusters are spherical and equally sized.
- **Disadvantages**: Requires the number of clusters (`k`) to be specified beforehand; sensitive to outliers; assumes clusters are spherical.

#### DBSCAN

- **Description**: DBSCAN (Density-Based Spatial Clustering of Applications with Noise) clustering groups points that are closely packed together and marks as outliers points that lie alone in low-density regions.
- **Advantages**: Does not require the number of clusters to be specified; can find arbitrarily shaped clusters; robust to outliers.
- **Disadvantages**: Parameters `eps` (maximum distance) and `min_samples` (minimum number of points) need to be chosen carefully; may struggle with varying densities.

## Running the Tool

### Step-by-Step Instructions

1. **Prepare Your Data**:
   - Ensure your scRNA-seq data is stored in an `anndata` object and saved in h5ad format.

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