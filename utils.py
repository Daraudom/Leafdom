import numpy as np
import matplotlib.pyplot as plt
from sklearn.preprocessing import StandardScaler
import pandas as pd

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
    Compute PCA manually on the data.
    
    Parameters:
    data (np.ndarray): Standardized data.
    n_components (int): Number of principal components to select.
    
    Returns:
    np.ndarray: PCA-transformed data.
    """
    # Step 2: Compute the covariance matrix
    covariance_matrix = np.cov(data, rowvar=False)
    
    # Step 3: Perform eigendecomposition
    eigenvalues, eigenvectors = np.linalg.eigh(covariance_matrix)
    
    # Step 4: Sort eigenvalues and corresponding eigenvectors
    sorted_idx = np.argsort(eigenvalues)[::-1]
    sorted_eigenvalues = eigenvalues[sorted_idx]
    sorted_eigenvectors = eigenvectors[:, sorted_idx]
    
    # Step 5: Select the top n_components eigenvectors
    top_eigenvectors = sorted_eigenvectors[:, :n_components]
    
    # Step 6: Project the data onto the top eigenvectors
    pca_data = np.dot(data, top_eigenvectors)
    
    return pca_data

def kmeans(pca_data, n_clusters, max_iters=300, tol=1e-4):
    """
    Perform K-means clustering.
    
    Parameters:
    pca_data (np.ndarray): Data transformed by PCA.
    n_clusters (int): Number of clusters.
    max_iters (int): Maximum number of iterations.
    tol (float): Tolerance to declare convergence.
    
    Returns:
    np.ndarray: Cluster labels for each data point.
    np.ndarray: Centroids of the clusters.
    """
    np.random.seed(42)
    n_samples = pca_data.shape[0]
    centroids = pca_data[np.random.choice(n_samples, n_clusters, replace=False)]
    
    for _ in range(max_iters):
        distances = np.sqrt(((pca_data - centroids[:, np.newaxis])**2).sum(axis=2))
        labels = np.argmin(distances, axis=0)
        new_centroids = np.array([pca_data[labels == i].mean(axis=0) for i in range(n_clusters)])
        
        if np.all(np.abs(new_centroids - centroids) < tol):
            break
        centroids = new_centroids
        
    return labels, centroids

def dbscan(pca_data, eps, min_samples):
    """
    Perform DBSCAN clustering.
    
    Parameters:
    pca_data (np.ndarray): Data transformed by PCA.
    eps (float): The maximum distance between two samples for them to be considered as in the same neighborhood.
    min_samples (int): The number of samples in a neighborhood for a point to be considered as a core point.
    
    Returns:
    np.ndarray: Cluster labels for each data point.
    """
    n_samples = pca_data.shape[0]
    labels = -np.ones(n_samples, dtype=int)
    core_samples_mask = np.zeros_like(labels, dtype=bool)
    cluster_id = 0
    
    for i in range(n_samples):
        if labels[i] != -1:
            continue
        
        neighbors = np.where(np.linalg.norm(pca_data - pca_data[i], axis=1) < eps)[0]
        
        if len(neighbors) < min_samples:
            labels[i] = -1
        else:
            labels[i] = cluster_id
            core_samples_mask[neighbors] = True
            stack = list(neighbors)
            
            while stack:
                point = stack.pop()
                if labels[point] == -1:
                    labels[point] = cluster_id
                
                if labels[point] != cluster_id:
                    continue
                
                point_neighbors = np.where(np.linalg.norm(pca_data - pca_data[point], axis=1) < eps)[0]
                if len(point_neighbors) >= min_samples:
                    for pn in point_neighbors:
                        if labels[pn] == -1:
                            labels[pn] = cluster_id
                        if not core_samples_mask[pn]:
                            stack.append(pn)
                            core_samples_mask[pn] = True
            
            cluster_id += 1
    
    return labels

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
        cluster_labels = dbscan(pca_data, eps, min_samples)
    elif method == 'kmeans':
        cluster_labels, _ = kmeans(pca_data, n_clusters)
    else:
        raise ValueError("Unsupported clustering method. Choose 'dbscan' or 'kmeans'.")
    return cluster_labels

# def identify_marker_genes(adata, cluster_labels, top_n=3):
#     """
#     Identify marker genes for each cluster.
    
#     Parameters:
#     adata (anndata.AnnData): Annotated data object containing the scRNA-seq data.
#     cluster_labels (np.ndarray): Cluster labels for each data point.
#     top_n (int): Number of top marker genes to select for each cluster.
    
#     Returns:
#     dict: Dictionary with cluster labels as keys and list of marker genes as values.
#     """
#     markers = {}
#     unique_clusters = np.unique(cluster_labels)
    
#     var_names = np.array(adata.var_names)  # Convert to numpy array
    
#     for cluster in unique_clusters:
#         if cluster == -1:
#             continue  # Skip noise points
#         cluster_data = adata[cluster_labels == cluster].X
#         mean_expression = np.mean(cluster_data, axis=0)
#         top_genes_idx = np.argsort(mean_expression)[::-1][:top_n]
#         markers[cluster] = var_names[top_genes_idx].tolist()
    
#     return markers

def visualize_pca_clusters(pca_data, cluster_labels, adata, output_file):
    """
    Visualize the clusters in the PCA-transformed data with different colors and save the plot to a file.
    
    Parameters:
    pca_data (np.ndarray): Data transformed by PCA.
    cluster_labels (np.ndarray): Cluster labels for each data point.
    output_file (str): Path to the output file for the PCA plot.
    
    Returns:
    None
    """
    # Identify marker genes for each cluster
   # markers = identify_marker_genes(adata, cluster_labels, top_n=3)

    plt.figure(figsize=(10, 7))
    scatter = plt.scatter(pca_data[:, 0], pca_data[:, 1], c=cluster_labels, cmap='viridis', marker='o', s=3)
    plt.xlabel('Principal Component 1')
    plt.ylabel('Principal Component 2')
    plt.title('PCA of scRNA-seq Data')

    # Calculate centroids for each cluster and label them
    unique_clusters = np.unique(cluster_labels)
    for cluster in unique_clusters:
        if cluster == -1:
            # Skip labeling for noise points in DBSCAN (cluster label -1)
            continue
        centroid = pca_data[cluster_labels == cluster].mean(axis=0)

        # marker_genes = ", ".join(map(str, markers[cluster][:1]))  # added
        # plt.text(centroid[0], centroid[1], marker_genes, fontsize=8, weight='bold', color='white', ha='center', va='center', 
        #          bbox=dict(facecolor='black', alpha=0.5, edgecolor='none', boxstyle='round,pad=0.3'))

        plt.text(centroid[0], centroid[1], str(cluster + 1), fontsize=12, weight='bold', color='white', ha='center', va='center', 
                 bbox=dict(facecolor='black', alpha=0.5, edgecolor='none', boxstyle='round,pad=0.3'))
        
    plt.savefig(output_file)
    plt.close()
