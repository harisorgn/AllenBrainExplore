module AllenBrainExplore

using AllenSDK
using PythonCall
using DataFrames

using Muon

using OhMyThreads

using SparseArrays
using AxisArrays: AxisArray

using HTTP
using FileIO

using GLMakie

include("utils.jl")
include("connectivity.jl")
include("cell_types.jl")
include("genes.jl")

export plot_connection_density
export get_expression_matrix, get_gene_description

end