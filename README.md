# SparseClustering.jl

A Julia package specialized for hierarchical clustering on sparse data structures with spatial and temporal structure. The clustering algorithm in this package can keep track of this structure, and cluster only across the time dimension, while leaving unclustered different regions in space at the same time point. The clusters are represented with an efficient implementation of the Union-Find data structure to support fast clustering operations. Between that and its use of sorting the distance matrix, this package achieves computational complexity of $O(n k \log(n k) + n t)$, where $n$ is the number of ROIs, $k$ is the average number of nonzero elements per row in the distance matrix, and $t$ is the number of time points. This can be much faster than the typical $O(n^2)$ complexity of hierarchical clustering.

## Features:

- Hierarchical clustering for sparse matrices with constraints on merging based on overlap and height thresholds.
- Efficient Union-Find data structure tailored to support the package's clustering methods.

## Installation

```julia
using Pkg
Pkg.add("git@github.com:flavell-lab/SparseClustering.jl")
```

## Usage:

### 1. Hierarchical Clustering with Thresholds

The primary feature of `SparseClustering.jl` is its ability to cluster a distance matrix while considering specific overlap and height thresholds.

```julia
using SparseClustering

ds = ... # Your distance matrix
inv_map = ... # Your mapping of ROIs to time points
overlap_threshold = 0.5
height_threshold = 2.0

# Perform clustering
clusters = hclust_minimum_threshold_sparse(ds, inv_map, overlap_threshold, height_threshold)
```

In this method:

- `ds`: Represents a matrix of pairwise distances between ROIs.
- `inv_map`: A dictionary mapping ROIs to time points.
- `overlap_threshold`: Clusters won't be merged if the merged cluster would contain more than this fraction of its ROIs from the same time point.
- `height_threshold`: The maximum distance between two ROIs that can be merged. If any two ROIs have a distance larger than this value, the clustering process terminates, returning the current clusters.

### 2. UnionFind Data Structure

While the main feature of the package is the clustering algorithm, at its core is the Union-Find data structure, which can also be used directly:

```julia
# Create a new UnionFind structure with 5 elements
uf = UnionFind(5)

# Merge subsets containing elements 1 and 2
union!(uf, 1, 2)

# Find the representative of the element 1
root_of_1 = find(uf, 1)
```

The `UnionFind` structure provides two primary operations: `find`, which returns the representative element of a subset, and `union!`, which merges two subsets. It operates in amoritized time $O(\alpha(n))$, where $\alpha$ is the inverse Ackermann function, which is essentially constant for any practical value of $n$.

## License:

[MIT License](link-to-license)
