# leafdom.py
import argparse
import numpy as np
from utils import UMAP

def main():
    parser = argparse.ArgumentParser(description='Run UMAP dimensionality reduction.')
    parser.add_argument('--input', type=str, required=True, help='Path to the input CSV file.')
    parser.add_argument('--output', type=str, required=True, help='Path to the output CSV file.')
    parser.add_argument('--n_neighbors', type=int, default=15, help='Number of nearest neighbors (default: 15).')
    parser.add_argument('--min_dist', type=float, default=0.1, help='Minimum distance between points in the low-dimensional space (default: 0.1).')
    parser.add_argument('--n_epochs', type=int, default=200, help='Number of training epochs (default: 200).')
    parser.add_argument('--n_components', type=int, default=2, help='Number of dimensions in the output (default: 2).')

    args = parser.parse_args()

    # Load the input data
    data = np.loadtxt(args.input, delimiter=',')

    # Run UMAP
    embedding = UMAP(data, n=args.n_neighbors, d=args.n_components, min_dist=args.min_dist, n_epochs=args.n_epochs)

    # Save the output
    np.savetxt(args.output, embedding, delimiter=',')

if __name__ == '__main__':
    main()
