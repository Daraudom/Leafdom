import argparse
import anndata as ad
from utils import preprocess_anndata, compute_pca, cluster_data, visualize_pca_clusters
from memory_profiler import profile

@profile
def main(input_path, output_file, method='kmeans', k=2, eps=0.5, min_samples=10, n_clusters=3):
    print("Reading the anndata...")
    adata = ad.read_h5ad(input_path)
    print("Processing the anndata...")
    data = preprocess_anndata(adata)
    print("Computing PCA...")
    pca_data = compute_pca(data, n_components=k)
    print("Determining Clusters...")
    cluster_labels = cluster_data(pca_data, method=method, eps=eps, min_samples=min_samples, n_clusters=n_clusters)
    print("Visualizing Plot...")
    visualize_pca_clusters(pca_data, cluster_labels, adata, output_file)
    print("Plot saved and ready to view...")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="LeafDom: A tool for scRNA-seq data preprocessing, PCA, clustering, and visualization.")
    parser.add_argument('-i', '--input', required=True, help="Path to the input anndata file (h5ad format).")
    parser.add_argument('-o', '--output_file', required=True, help="Path to the output file for the PCA plot.")
    parser.add_argument('-m', '--method', choices=['dbscan', 'kmeans'], default='kmeans', help="Clustering method (default: dbscan).")
    parser.add_argument('-k', '--components', type=int, default=2, help="Number of principal components to select for PCA (default: 2).")
    parser.add_argument('-e', '--eps', type=float, default=0.5, help="Maximum distance between two samples for DBSCAN (default: 0.5).")
    parser.add_argument('-s', '--min_samples', type=int, default=10, help="Number of samples in a neighborhood for a point to be considered as a core point in DBSCAN (default: 10).")
    parser.add_argument('-c', '--clusters', type=int, default=3, help="Number of clusters to form for k-means (default: 3).")
    
    args = parser.parse_args()
    main(args.input, args.output_file, args.method, args.components, args.eps, args.min_samples, args.clusters)
