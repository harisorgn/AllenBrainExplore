function get_expression_matrix(dataset_name, region)
    abcpc = pyimport("abc_atlas_access.abc_atlas_cache.abc_project_cache")
    pc = abcpc.AbcProjectCache
    abc_cache = pc.from_cache_dir("./data/abc_atlas/")

    dn = occursin("imputed", dataset_name) ? join(split(dataset_name, "-")[1:end-1], "-") : dataset_name

    cell_metadata = abc_cache.get_metadata_dataframe(directory=dn, file_name="cell_metadata")
    df_cell = to_dataframe(cell_metadata)

    gene = abc_cache.get_metadata_dataframe(directory=dataset_name, file_name="gene")
    df_gene = to_dataframe(gene)

    #cluster_details = abc_cache.get_metadata_dataframe(
    #    directory="WMB-taxonomy",
    #    file_name="cluster_to_cluster_annotation_membership_pivoted",
    #    keep_default_na=false
    #)
    #df_cluster = to_dataframe(cluster_details)

    ccf_coordinates = abc_cache.get_metadata_dataframe(directory="$(dn)-CCF", file_name="ccf_coordinates")
    df_ccf = to_dataframe(ccf_coordinates)

    parcellation_annotation = abc_cache.get_metadata_dataframe(
        directory="Allen-CCF-2020",
        file_name="parcellation_to_parcellation_term_membership_acronym"
    )
    df_parcel = to_dataframe(parcellation_annotation)

    parcel_idxs = df_parcel[df_parcel.structure .== region, :parcellation_index]
    ccf_idxs = findall(idx -> idx in parcel_idxs, df_ccf.parcellation_index)
    region_labels = df_ccf[ccf_idxs, :cell_label]

    cell_idxs = tmap(region_labels) do l
        findfirst(x -> occursin(l, x), df_cell.cell_label)
    end
    @show length(cell_idxs)

    expr_mat_file_name = occursin("MERFISH", dataset_name) ? join(split(dataset_name, "-")[2:end], "-") : dataset_name

    file = abc_cache.get_file_path(directory=dataset_name, file_name="$(expr_mat_file_name)/log2")
    file = join([pyconvert(String, p) for p in file.parts], "/")

    adata = readh5ad(file)
    df = DataFrame(adata[cell_idxs, :])
    rename!(df, vcat(["obs"], lowercase.(df_gene.gene_symbol)))

    # HACK: Type conversion is necessary because the original type contains Int32 instead of Int(64)
    # and some type assertions in Tables fail.
    for gene_name in names(df)
        if !occursin("obs", gene_name)
            df[!, gene_name] = convert(SparseArrays.SparseVector{Float64, Int}, df[:, gene_name])
        end
    end
    
    return df
end

