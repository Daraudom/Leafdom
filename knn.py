import numpy as np
from scipy.spatial.distance import cdist
    
# Step 2: Compute the distance matrix
def compute_distance_matrix(data):
    return cdist(data, data, metric='euclidean')

# Step 3: Find the K nearest neighbors
def find_k_nearest_neighbors(distance_matrix, k):
    num_points = distance_matrix.shape[0]
    knn_indices = np.zeros((num_points, k), dtype=int)
    knn_distances = np.zeros((num_points, k))
    
    for i in range(num_points):
        sorted_indices = np.argsort(distance_matrix[i])
        knn_indices[i] = sorted_indices[1:k+1]  # Exclude the point itself (first index)
        knn_distances[i] = distance_matrix[i, knn_indices[i]]
    
    return knn_indices, knn_distances

# Step 4: Construct the KNN graph
def knn_graph(data, k):
    distance_matrix = compute_distance_matrix(data)
    knn_indices, knn_distances = find_k_nearest_neighbors(distance_matrix, k)
    return knn_indices, knn_distances

