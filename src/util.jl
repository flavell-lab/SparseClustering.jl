"""
    sort_distance_matrix(s::SparseMatrixCSC{T}) where T<:Real

Sort the non-zero entries of the upper triangular part of a sparse distance matrix 
and return the (i,j) indices of these entries sorted by their values.

# Arguments:
- `s::SparseMatrixCSC{T}`: A sparse matrix, typically representing a distance matrix.

# Returns:
- `Array{Tuple{Int, Int}}`: A list of (i,j) tuples sorted by the distance values 
   from the matrix `s`.

# Notes:
The function only considers the upper triangular part of the matrix 
since distance matrices are typically symmetric and we avoid considering duplicate entries.

# Example:
```julia
s = sparse([0.0 2.0; 2.0 0.0])
result = sort_distance_matrix(s)
"""
function sort_distance_matrix(s::SparseMatrixCSC{T}) where T<:Real
    pairs = Vector{Tuple{Int, Int, T}}()

    for col in 1:size(s, 2)
        for idx in nzrange(s, col)
            row = rowvals(s)[idx]
            value = s.nzval[idx]
            if row < col  # This ensures we only look at the upper triangle of the matrix
                push!(pairs, (row, col, value))
            end
        end
    end
    
    # Now, sort the pairs by the distance values
    sort!(pairs, by=x -> x[3])

    # Return only the (i,j) pairs, discarding the values
    return [(x[1], x[2]) for x in pairs]
end


"""
    generate_timepoint_map(inv_map, n, max_timept)

Generate a matrix representing time points where each ROI was detected at.

# Arguments:
- `inv_map`: A dictionary where each key is an ROI and each value is 
  a dictionary whose keys are time points where that ROI was found
- `n`: Number of ROIs.
- `max_timept`: The maximum number of time points in the dataset.

# Returns:
- `Matrix{Int32}`: A matrix of dimensions `n x max_timept`. The entry (i,j) 
  is 1 if the `i`-th tree node matches the `j`-th frame, and 0 otherwise.

# Example:
```julia
inv_map = Dict(1 => Set([2]), 2 => Set([1,3]))
n = 3
max_timept = 2
result = generate_timepoint_map(inv_map, n, max_frame)
```
"""
function generate_timepoint_map(inv_map, n, max_timept)
    timepoint_map = zeros(Int32, n, max_timept)
    for i=1:n
        matching_keys = collect(keys(inv_map[i]))
        if length(matching_keys) > 0
            for key in matching_keys
                timepoint_map[i, key] = 1
            end
        end
    end
    return timepoint_map
end
