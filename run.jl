using AllenBrainExplore
using DataFrames
using Statistics
using CSV

datasets = ["MERFISH-C57BL6J-638850", "Zhuang-ABCA-1", "Zhuang-ABCA-2", "Zhuang-ABCA-3", "Zhuang-ABCA-4", "MERFISH-C57BL6J-638850-imputed"]

region = "LH"

df = get_expression_matrix(datasets[1], region)
GC.gc()

dfg = get_gene_description(datasets[1], region)
filter!(r -> occursin("receptor", r.name), dfg)

G = Matrix(df[:, lowercase.(dfg.gene_symbol)])

Gm = vec(mean(G; dims=1))

idxs = sortperm(Gm; rev=true)

dff = DataFrame([:gene => dfg.gene_symbol[idxs], :description => dfg.name[idxs], :expression_LH => Gm[idxs]])

CSV.write("receptor_genes_$(region).csv", dff)
