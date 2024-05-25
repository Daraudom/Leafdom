import numpy as np
from scipy.spatial import cKDTree
from scipy.sparse.csgraph import laplacian
from scipy.linalg import eigh
from scipy.optimize import curve_fit

# Define the smooth approximation of Φ
def phi_smooth(x, a, b):
    return (1 + a * (x ** 2) ** b) ** -1

# Define the Ψ function
def psi(x, y, min_dist):
    distance_squared = np.linalg.norm(x - y) ** 2
    if distance_squared <= min_dist:
        return 1.0
    else:
        return np.exp(-(distance_squared - min_dist))

# fuzzy set
def LocalFuzzySimplicialSetForEachPoint(X, n):
    all_fs_sets_0 = []
    all_fs_sets_1 = []
   
    for x in X:
        fs_set_0, fs_set_1 = LocalFuzzySimplicialSet(X, x, n)
        all_fs_sets_0.append(fs_set_0)
        all_fs_sets_1.append(fs_set_1)
   
    return all_fs_sets_0, all_fs_sets_1

def LocalFuzzySimplicialSet(X, x, n: int):
    # Step 1: Find the k nearest neighbors and their distances
    knn_tree = cKDTree(X)
    _, knn_indices = knn_tree.query(x, k=n+1)  # Include x itself
    knn_indices = knn_indices[1:]  # Exclude x itself
    knn_dists = np.linalg.norm(X[knn_indices] - x, axis=1)
    rho = knn_dists[0]  # Distance to nearest neighbor
   
    # Step 2: Smooth approximation of knn distances
    sigma = SmoothKNNDist(knn_dists, n, rho)
   
    # Step 3: Initialize fuzzy simplicial set
    fs_set_0 = set(range(len(X)))
    fs_set_1 = {frozenset([x_idx, y_idx]): 0 for y_idx in knn_indices}
   
    # Step 4: Calculate weights for edges
    for y_idx in knn_indices:
        dx_y = np.maximum(0, np.linalg.norm(X[y_idx] - x) - rho) / sigma
        weight = np.exp(-dx_y)
        fs_set_1[frozenset([x_idx, y_idx])] = weight
   
    return fs_set_0, fs_set_1

def SmoothKNNDist(knn_dists, n, rho):
    # Binary search for sigma
    lower_bound = 0
    upper_bound = 1e5  # Some large value
    target_sum = np.log2(n)
   
    while upper_bound - lower_bound > 1e-5:
        sigma = (upper_bound + lower_bound) / 2
        sum_exp = np.sum(np.exp(-(knn_dists - rho) / sigma))
       
        if sum_exp > target_sum:
            lower_bound = sigma
        else:
            upper_bound = sigma
   
    return (upper_bound + lower_bound) / 2

def UMAP(X, n, d, min_dist, n_epochs):
    # Construct the relevant weighted graph
    all_fs_sets_0, all_fs_sets_1 = LocalFuzzySimplicialSetForEachPoint(X, n)
    top_rep = merge_fs_sets(all_fs_sets_1)  # Assuming merge function is defined

    # Perform optimization of the graph layout
    Y = SpectralEmbedding(top_rep, d)
    Y = optimize_embedding(top_rep, Y, min_dist, n_epochs, n_neg_samples=1)
   
    return Y

def SpectralEmbedding(all_fs_sets_1, d):
    # Merge all fuzzy simplicial sets into a single topological representation
    merged_top_rep = merge_fs_sets(all_fs_sets_1)

    # Convert topological representation to adjacency matrix
    A = top_rep_to_adjacency(merged_top_rep)

    # Compute degree matrix
    D = np.diag(np.sum(A, axis=1))

    # Compute Laplacian matrix
    L = laplacian(A, normed=False)

    # Eigen decomposition
    _, evec = eigh(L, eigvals=(0, d))  # Assumes eigenvalues are sorted in ascending order

    # Select embedding
    Y = evec[:, 1:d+1]  # Exclude the trivial constant eigenvector

    return Y

def optimize_embedding(top_rep, Y, min_dist, n_epochs, n_neg_samples):
    alpha = 1.0
   
    # Fit Φ from Ψ defined by min_dist
    X = np.random.rand(len(Y), 2)  # Sample X for fitting Φ from Ψ
    y = np.random.rand(len(Y))  # Sample y for fitting Φ from Ψ
    params, _ = curve_fit(phi_smooth, y, psi(X[0], X[1], min_dist))
    a_fit, b_fit = params

    for epoch in range(1, n_epochs + 1):
        for ([a, b], p) in top_rep:
            if np.random.random() <= p:
                # Sample simplex with probability p
                ya = Y[a]
                yb = Y[b]
               
                # Update ya using gradient of log(Φ)
                gradient_phi = gradient_log_phi(ya, yb, min_dist)
                ya += alpha * gradient_phi
               
                # Update yb using gradient of log(Φ)
                gradient_phi = gradient_log_phi(yb, ya, min_dist)
                yb += alpha * gradient_phi
               
                for _ in range(n_neg_samples):
                    c = np.random.choice(len(Y))
                    yc = Y[c]
                   
                    # Update ya using gradient of log(1 - Φ)
                    gradient_neg_phi = gradient_log_neg_phi(ya, yc, min_dist)
                    ya += alpha * gradient_neg_phi
               
        # Update learning rate
        alpha = 1.0 - epoch / n_epochs
   
    return Y

# Helper functions for merging, converting, and computing gradients
def merge_fs_sets(all_fs_sets_1):
    merged_top_rep = {}
    for fs_set in all_fs_sets_1:
        for edge, weight in fs_set.items():
            if edge not in merged_top_rep:
                merged_top_rep[edge] = weight
            else:
                merged_top_rep[edge] = max(merged_top_rep[edge], weight)
    return merged_top_rep

def top_rep_to_adjacency(top_rep):
    n_vertices = len(top_rep)
    A = np.zeros((n_vertices, n_vertices))
    for edge, weight in top_rep.items():
        u, v = edge
        A[u][v] = weight
        A[v][u] = weight  # Assuming undirected graph
    return A

def gradient_log_phi(x, y, min_dist):
    distance_squared = np.linalg.norm(x - y) ** 2
    if distance_squared <= min_dist:
        return np.zeros_like(x)
    else:
        return -2 * a_fit * b_fit * (1 + a_fit * distance_squared / 2) ** (-b_fit - 1) * (x - y)

def gradient_log_neg_phi(x, y, min_dist):
    distance_squared = np.linalg.norm(x - y) ** 2
    if distance_squared <= min_dist:
        return np.zeros_like(x)
    else:
        return -np.exp(-(distance_squared - min_dist)) * (x - y)

# Example usage
X = np.random.rand(100, 2)  # Sample dataset (100 points in 2D)
n = 5  # Number of nearest neighbors
d = 2  # Desired dimensionality of the embedding
min_dist = 0.1  # Minimum distance between embedded points
n_epochs = 50  # Number of optimization epochs
embedding = UMAP(X, n, d, min_dist, n_epochs)
print("Embedded points shape:", embedding.shape)
