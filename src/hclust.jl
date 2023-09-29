"""
    hclust_minimum_threshold_sparse(ds::AbstractMatrix{T}, inv_map, overlap_threshold, height_threshold; use_sparse=true, pair_match=false) where T<:Real

Clusters a matrix of pairwise distances between ROIs. Avoids merginc clusters that would cause ROIs from the same time point to be clustered together.
Stops the clustering algorithm when the distances become too large and terminates, returning the current clusters.

# Arguments:
- `ds::AbstractMatrix{T}`: A matrix of pairwise distances between ROIs.
- `inv_map`: A dictionary mapping ROIs to time points.
- `overlap_threshold`: A fraction between 0 and 1. If a cluster has more than `overlap_threshold` of its ROIs from the same time point, it is not merged.
- `height_threshold`: The maximum distance between two ROIs that can be merged.
- `use_sparse`: A boolean indicating whether to use sparse matrices to store intermediate data.
- `pair_match`: A boolean indicating whether to only merge clusters of size at most 2 (turning it into a matching algorithm).

# Returns:
- `UnionFind`: A UnionFind object representing the current clusters.
"""
function hclust_minimum_threshold_sparse(ds::AbstractMatrix{T}, inv_map, overlap_threshold, height_threshold; use_sparse=true, pair_match=false) where T<:Real
    if use_sparse
        d = copy(sparse(ds))
    else
        d = Matrix(ds)
    end

    n = size(d, 1)
    
    sorted_pairs_list = sorted_pairs(d)
    
    max_frame = maximum([maximum([y for y in keys(x)]) for x in values(inv_map)])
    tree_frames = generate_tree_frames(inv_map, n, max_frame)

    curr_cluster_ids = UnionFind(n)

    if pair_match
        merged_nodes = zeros(Bool, n)
    end

    new_tree_frame = similar(tree_frames, size(tree_frames, 2))

    for pair in sorted_pairs_list
        i, j = pair
        
        # too high - finish clustering
        if d[i, j] > height_threshold
            break
        end

        clust_i = find(curr_cluster_ids, i)
        clust_j = find(curr_cluster_ids, j)

        # If the clusters i and j have already been merged, skip
        if clust_i == clust_j || (pair_match && (merged_nodes[i] || merged_nodes[j]))
            continue
        end
        
        # Check if merge is valid under ROI collision criterion
        new_tree_frame .= tree_frames[clust_i,:] .+ tree_frames[clust_j,:]
        overlaps = count(x -> x > 1, new_tree_frame)
        ratio = overlaps / count(>(0), new_tree_frame)
        
        if ratio > overlap_threshold
            continue
        end
        
        # Update clusters
        union!(curr_cluster_ids, clust_i, clust_j)
        
        # Handle blocked merges for i, j, and last_tree
        # (similar to what you had in your original code)

        tree_frames[clust_i, :] .= new_tree_frame    
        tree_frames[clust_j, :] .= new_tree_frame

        if pair_match
            merged_nodes[i] = true
            merged_nodes[j] = true
        end
    end

    return curr_cluster_ids
end
