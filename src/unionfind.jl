"""
    UnionFind(n::Int)

A data structure that keeps track of a partition of a set into disjoint subsets.
Supports two primary operations: `find`, which determines the representative 
element of a subset, and `union`, which merges two subsets into one.

# Fields:
- `parent::Vector{Int}`: An array that holds the representative element for each subset.
- `rank::Vector{Int}`: An array that holds the depth of trees. It's used to keep the tree flat 
  when two trees of equal rank are merged.

# Example:
```julia
uf = UnionFind(5)
union!(uf, 1, 2)
union!(uf, 3, 4)
println(find(uf, 1))  # Could be 1, 2, 3, or 4, depending on the internal structure.
```
"""
mutable struct UnionFind
    parent::Vector{Int}
    rank::Vector{Int}
    
    UnionFind(n::Int) = new(collect(1:n), zeros(Int, n))
end


"""
    find(uf::UnionFind, x::Int) -> Int

Returns the representative of the subset containing `x`.

# Arguments:
- `uf::UnionFind`: A `UnionFind` object.
- `x::Int`: An integer representing the element whose representative we want to find.

# Returns:
- `Int`: The representative of the subset containing `x`.

Uses path compression to ensure that every visited node points directly 
to the root in subsequent operations.
"""
function find(uf::UnionFind, x::Int)
    if uf.parent[x] != x
        uf.parent[x] = find(uf, uf.parent[x])  # Path compression
    end
    return uf.parent[x]
end


"""
    union!(uf::UnionFind, x::Int, y::Int)

Merges the subsets containing `x` and `y`. If `x` and `y` are already in the 
same subset, the function does nothing. Otherwise, it combines the two subsets.

# Arguments:
- `uf::UnionFind`: A `UnionFind` object.
- `x::Int` and `y::Int`: Integers representing elements we wish to union.

Uses the "union by rank" optimization to keep the tree as flat as possible:
1. If the ranks of the two subsets' roots are different, the root of the subset with 
   the smaller rank is attached to the root of the subset with the larger rank.
2. If the ranks are the same, one is arbitrarily attached to the other and 
   the rank of the resulting tree is increased by 1.
"""
function union!(uf::UnionFind, x::Int, y::Int)
    rootX = find(uf, x)
    rootY = find(uf, y)
    
    if rootX != rootY
        # Union by rank to ensure the tree remains flat
        if uf.rank[rootX] > uf.rank[rootY]
            uf.parent[rootY] = rootX
        else
            uf.parent[rootX] = rootY
            if uf.rank[rootX] == uf.rank[rootY]
                uf.rank[rootY] += 1
            end
        end
    end
end
