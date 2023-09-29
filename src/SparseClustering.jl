module SparseClustering
    using SparseArrays

    include("unionfind.jl")
    include("hclust.jl")
    include("util.jl")
    export
        UnionFind,
        find,
        union!,
        hclust_minimum_threshold_sparse
end
