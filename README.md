# Leafdom

Leafdom is a command-line tool that implements the Uniform Manifold Approximation and Projection (UMAP) algorithm from scratch. UMAP is a dimension reduction technique that can be used for visualization and preserving structure in high-dimensional data. Leafdom allows users to easily apply UMAP to their datasets via the command line.

## Features

- Implements UMAP from scratch
- Supports various distance metrics
- Allows easy integration into command-line workflows
- Provides options for customizing the UMAP parameters
- Outputs transformed data for visualization and further analysis

## Installation

Clone the repository:

```bash
git clone https://github.com/yourusername/leafdom.git
cd leafdom
```

Install the required dependencies:

```bash
pip install -r requirements.txt
```

## Usage

To use Leafdom, run the following command with your dataset:

```bash
python leafdom.py --input data.csv --output transformed_data.csv --n_neighbors 15 --min_dist 0.1 --metric euclidean
```

### Command-Line Options

- `--input`: Path to the input dataset (CSV file).
- `--output`: Path to the output transformed dataset (CSV file).
- `--n_neighbors`: The size of local neighborhood (in terms of number of neighboring sample points) used for manifold approximation. (default: 15)
- `--min_dist`: The effective minimum distance between embedded points. Smaller values will result in a more clustered/clumped embedding. (default: 0.1)
- `--metric`: The metric to use for measuring distance in the input space. (default: euclidean)

## Examples

### Basic Usage

Transform a dataset with default parameters:

```bash
python leafdom.py --input mydata.csv --output mytransformeddata.csv
```

### Custom Parameters

Transform a dataset with custom parameters for `n_neighbors` and `min_dist`:

```bash
python leafdom.py --input mydata.csv --output mytransformeddata.csv --n_neighbors 10 --min_dist 0.2
```

## Development

To contribute to Leafdom, follow these steps:

1. Fork the repository.
2. Create a new branch: `git checkout -b feature-branch`.
3. Make your changes and commit them: `git commit -m 'Add some feature'`.
4. Push to the branch: `git push origin feature-branch`.
5. Submit a pull request.

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Acknowledgements

- UMAP: Uniform Manifold Approximation and Projection by Leland McInnes, John Healy, and James Melville.

## Contact

For any questions or suggestions, please contact [your email address].

---