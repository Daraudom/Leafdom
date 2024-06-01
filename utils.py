import numpy as np
import matplotlib.pyplot as plt
from sklearn.preprocessing import StandardScaler
from sklearn.cluster import DBSCAN, KMeans
from sklearn.decomposition import PCA

def preprocess_anndata(adata):
    """
    Preprocess scRNA-seq data stored in an anndata object.
    
    Parameters:
    adata (anndata.AnnData): Annotated data object containing the scRNA-seq data.
    
    Returns:
    np.ndarray: Standardized data ready for PCA.
    """
    data = adata.X

    # Normalizing the data per cell
    if not isinstance(data, np.ndarray):
        data = data.toarray()
    
    data = data / np.sum(data, axis=1, keepdims=True) * 1e4
    
    # Logarithmizing the data
    data = np.log1p(data)
    
    # Standardizing the data
    scaler = StandardScaler()
    scaled_data = scaler.fit_transform(data)
    return scaled_data

def compute_pca(data, n_components=2):
    """
    Compute PCA on the data.
    
    Parameters:
    data (np.ndarray): Standardized data.
    n_components (int): Number of principal components to select.
    
    Returns:
    np.ndarray: PCA-transformed data.
    """
    pca = PCA(n_components=n_components)
    pca_data = pca.fit_transform(data)
    return pca_data

def cluster_data(pca_data, method='dbscan', eps=0.5, min_samples=10, n_clusters=3):
    """
    Cluster the PCA-transformed data using the specified clustering method.
    
    Parameters:
    pca_data (np.ndarray): Data transformed by PCA.
    method (str): Clustering method ('dbscan' or 'kmeans').
    eps (float): The maximum distance between two samples for them to be considered as in the same neighborhood (for DBSCAN).
    min_samples (int): The number of samples in a neighborhood for a point to be considered as a core point (for DBSCAN).
    n_clusters (int): The number of clusters to form (for k-means).
    
    Returns:
    np.ndarray: Cluster labels for each data point.
    """
    if method == 'dbscan':
        dbscan = DBSCAN(eps=eps, min_samples=min_samples)
        cluster_labels = dbscan.fit_predict(pca_data)
    elif method == 'kmeans':
        kmeans = KMeans(n_clusters=n_clusters, random_state=42)
        cluster_labels = kmeans.fit_predict(pca_data)
    else:
        raise ValueError("Unsupported clustering method. Choose 'dbscan' or 'kmeans'.")
    return cluster_labels

def visualize_pca_clusters(pca_data, cluster_labels, output_file):
    """
    Visualize the clusters in the PCA-transformed data with different colors and save the plot to a file.
    
    Parameters:
    pca_data (np.ndarray): Data transformed by PCA.
    cluster_labels (np.ndarray): Cluster labels for each data point.
    output_file (str): Path to the output file for the PCA plot.
    
    Returns:
    None
    """
    plt.figure(figsize=(10, 7))
    scatter = plt.scatter(pca_data[:, 0], pca_data[:, 1], c=cluster_labels, cmap='viridis', marker='o', s=3)
    plt.xlabel('Principal Component 1')
    plt.ylabel('Principal Component 2')
    plt.title('PCA of scRNA-seq Data')
    plt.colorbar(scatter, label='Cluster')
    plt.savefig(output_file)
    plt.close()
