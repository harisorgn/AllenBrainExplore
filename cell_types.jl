using PythonCall
using DataFrames

include("utils.jl")

abcpc = pyimport("abc_atlas_access.abc_atlas_cache.abc_project_cache")

pc = abcpc.AbcProjectCache
abc_cache = pc.from_cache_dir("./data/abc_atlas/")

datasets = ["Zhuang-ABCA-1", "Zhuang-ABCA-2", "Zhuang-ABCA-3", "Zhuang-ABCA-4"]

df = abc_cache.get_metadata_dataframe(directory=datasets[1], file_name="cell_metadata")

cluster_details = abc_cache.get_metadata_dataframe(
    directory="WMB-taxonomy",
    file_name="cluster_to_cluster_annotation_membership_pivoted",
    keep_default_na=false
)

cluster_details.set_index("cluster_alias", inplace=true)

ccf_coordinates = abc_cache.get_metadata_dataframe(directory="$(datasets[1])-CCF", file_name="ccf_coordinates")
#ccf_coordinates.set_index("cell_label", inplace=true)
#ccf_coordinates.rename(columns=Dict("x" => "x_ccf", "y" => "y_ccf", "z" => "z_ccf"), inplace=true)

df = to_dataframe(ccf_coordinates)

parcellation_annotation = abc_cache.get_metadata_dataframe(
    directory="Allen-CCF-2020",
    file_name="parcellation_to_parcellation_term_membership_acronym"
)
parcellation_annotation.set_index("parcellation_index", inplace=true)


