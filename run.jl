using AllenBrainExplore
using DataFrames
using Statistics
using CSV

datasets = ["MERFISH-C57BL6J-638850", "Zhuang-ABCA-1", "Zhuang-ABCA-2", "Zhuang-ABCA-3", "Zhuang-ABCA-4", "MERFISH-C57BL6J-638850-imputed"]

df = get_expression_matrix(datasets[1], "VTA")

GC.gc()

G = Matrix(df[:, Not(:obs)])

Gm = vec(mean(G; dims=1))

idxs = sortperm(Gm; rev=true)

genes = names(df[:, Not(:obs)])

dfg = DataFrame([:gene => genes[idxs], :expression_VTA => Gm[idxs]])

CSV.write("genes_VTA.csv", dfg)