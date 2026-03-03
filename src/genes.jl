function get_expression_matrix(dataset_name, region)
    abcpc = pyimport("abc_atlas_access.abc_atlas_cache.abc_project_cache")
    pc = abcpc.AbcProjectCache
    abc_cache = pc.from_cache_dir("./data/abc_atlas/")

    dn = occursin("imputed", dataset_name) ? join(split(dataset_name, "-")[1:end-1], "-") : dataset_name

    cell_metadata = abc_cache.get_metadata_dataframe(directory=dn, file_name="cell_metadata")
    df_cell = to_dataframe(cell_metadata)

    gene = abc_cache.get_metadata_dataframe(directory=dataset_name, file_name="gene")
    df_gene = to_dataframe(gene)

    cluster_details = abc_cache.get_metadata_dataframe(
        directory="WMB-taxonomy",
        file_name="cluster_to_cluster_annotation_membership_pivoted",
        keep_default_na=false
    )
    df_cluster = to_dataframe(cluster_details)

    neuron_idxs = (.!occursin.("astro", lowercase.(df_cluster.class))) .& (.!occursin.("opc", lowercase.(df_cluster.class))) .& (.!occursin.("vascular", lowercase.(df_cluster.class))) .& (.!occursin.("immune", lowercase.(df_cluster.class)))
    neuron_cluster_alias = df_cluster.cluster_alias[neuron_idxs]
    filter!(r -> r.cluster_alias in neuron_cluster_alias, df_cell)

    ccf_coordinates = abc_cache.get_metadata_dataframe(directory="$(dn)-CCF", file_name="ccf_coordinates")
    df_ccf = to_dataframe(ccf_coordinates)

    parcellation_annotation = abc_cache.get_metadata_dataframe(
        directory="Allen-CCF-2020",
        file_name="parcellation_to_parcellation_term_membership_acronym"
    )
    df_parcel = to_dataframe(parcellation_annotation)

    parcel_idxs = df_parcel[df_parcel.structure .== region, :parcellation_index]
    if isempty(parcel_idxs)
        parcel_idxs = df_parcel[df_parcel.substructure .== region, :parcellation_index]
    end
    
    ccf_idxs = findall(idx -> idx in parcel_idxs, df_ccf.parcellation_index)
    region_labels = df_ccf[ccf_idxs, :cell_label]

    cell_idxs = tmap(region_labels) do l
        findfirst(x -> occursin(l, x), df_cell.cell_label)
    end
    filter!(!isnothing, cell_idxs)
    # explicitly cast Int to change from Vector{Union{Int, Nothing}} to Vector{Int} so that indexing adata doesn't error
    cell_idxs = Int.(cell_idxs) 

    @show length(cell_idxs)

    expr_mat_file_name = occursin("MERFISH", dataset_name) ? join(split(dataset_name, "-")[2:end], "-") : dataset_name

    file = abc_cache.get_file_path(directory=dataset_name, file_name="$(expr_mat_file_name)/log2")
    file = join([pyconvert(String, p) for p in file.parts], "/")

    adata = readh5ad(file)
    # HACK: Type conversion is necessary because the original type contains Int32 instead of Int(64)
    # and some type assertions in Tables fail when converting AnnData to DataFrame.
    #adata.X = convert(Adjoint{Float64, SparseMatrixCSC{Float64, Int}}, adata.X)

    df = DataFrame(adata[cell_idxs, :])
    rename!(df, vcat(["obs"], lowercase.(df_gene.gene_symbol));  makeunique=true)

    return df
end

function get_gene_description(dataset_name, region)
    abcpc = pyimport("abc_atlas_access.abc_atlas_cache.abc_project_cache")
    pc = abcpc.AbcProjectCache
    abc_cache = pc.from_cache_dir("./data/abc_atlas/")

    dn = occursin("imputed", dataset_name) ? join(split(dataset_name, "-")[1:end-1], "-") : dataset_name

    cell_metadata = abc_cache.get_metadata_dataframe(directory=dn, file_name="cell_metadata")
    df_cell = to_dataframe(cell_metadata)

    gene = abc_cache.get_metadata_dataframe(directory=dataset_name, file_name="gene")
    df_gene = to_dataframe(gene)

    return df_gene[!, [:gene_symbol, :name]]
end
