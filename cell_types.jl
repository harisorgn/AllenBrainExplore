using PythonCall
using DataFrames

include("utils.jl")

abcpc = pyimport("abc_atlas_access.abc_atlas_cache.abc_project_cache")

pc = abcpc.AbcProjectCache
abc_cache = pc.from_cache_dir("./data/abc_atlas/")

datasets = ["Zhuang-ABCA-1", "Zhuang-ABCA-2", "Zhuang-ABCA-3", "Zhuang-ABCA-4"]

cell_metadata = [abc_cache.get_metadata_dataframe(directory=d, file_name="cell_metadata") for d in datasets]
df_cell = reduce(vcat, to_dataframe.(cell_metadata))

cluster_details = abc_cache.get_metadata_dataframe(
    directory="WMB-taxonomy",
    file_name="cluster_to_cluster_annotation_membership_pivoted",
    keep_default_na=false
)
df_cluster = to_dataframe(cluster_details)

ccf_coordinates = [abc_cache.get_metadata_dataframe(directory="$(d)-CCF", file_name="ccf_coordinates") for d in datasets]
df_ccf = reduce(vcat, to_dataframe.(ccf_coordinates))

parcellation_annotation = abc_cache.get_metadata_dataframe(
    directory="Allen-CCF-2020",
    file_name="parcellation_to_parcellation_term_membership_acronym"
)
df_parcel = to_dataframe(parcellation_annotation)

vta_idxs = df_parcel[df_parcel.structure .== "VTA", :parcellation_index]
vta_labels = df_ccf[vta_idxs, :cell_label]

cell_idxs = map(vta_labels) do l
    findfirst(x -> occursin(l, x), df_cell.cell_label)
end

cluster_ids = df_cell[cell_idxs, :cluster_alias]

df_cluster[cluster_ids, :]