get_data(X::AbstractArray) = X
get_data(X::AxisArray) = X.data

function download_projection_density(id, filename, resolution)
    download("http://api.brain-map.org/grid_data/download_file/$(id)??image=projection_density&resolution=$(resolution)", filename * ".nrrd")

    return get_data(pyconvert(AbstractArray{Float64}, load("$(filename).nrrd")))
end

function download_injection_density(id, filename, resolution)
    download("http://api.brain-map.org/grid_data/download_file/$(id)??image=injection_density&resolution=$(resolution)", filename * ".nrrd")

    return get_data(pyconvert(AbstractArray{Float64}, load("$(filename).nrrd")))
end

function get_structure_mask(rspc, id)
    annotation, meta = rspc.get_annotation_volume()
    rsp = rspc.get_reference_space()
    stmap = rsp.make_structure_mask([id])

    return pyconvert(AbstractArray{Bool}, stmap)
end

function connection_density!(gridpos, D, MASK; title="", hemisphere=2)
    X = copy(D)
    X[.!(MASK)] .= zero(eltype(X))

    coords_dens = findall(!iszero, X)
    vals_dens = X[coords_dens]

    points_dens = map(coords_dens) do i                                                                                                    
        Point3f(Tuple(i)...)                                                                                                  
    end

    coords_region = findall(!iszero, MASK)
    lu_x = extrema(getindex.(coords_region, 1))
    lu_y = extrema(getindex.(coords_region, 2))
    z = getindex.(coords_region, 3)
    if hemisphere == 1 # left hemisphere
        z_hem = z[1:div(length(z), 2)]
        z_ticklabels = ["Lateral", "Medial"]
    elseif  hemisphere == 2 # right hemisphere
        z_hem = z[div(length(z), 2):end]
        z_ticklabels = ["Medial", "Lateral"]
    else # both hemispheres
        z_hem = z
        z_ticklabels = ["Left", "Right"]
    end
    lu_z = extrema(z_hem)

    points_region = map(coords_region) do i                                                                                                    
        Point3f(Tuple(i)...)                                                                                                  
    end

    ax = Axis3(
            gridpos[1,1], 
            xlabel="",
            ylabel="",
            zlabel="",
            xticks=([lu_x...], ["Anterior", "Posterior"]), 
            yticks=([lu_y...], ["Dorsal", "Ventral"]),
            zticks=([lu_z...], z_ticklabels),
            limits=(lu_x..., lu_y..., lu_z...),
            title=title
    )

    scatter!(ax, points_dens; color=vals_dens)
    scatter!(ax, points_region; color=(:grey, 0.1))
end

function plot_connection_density(S, D, S_mask, D_mask; source="", destination="", hemisphere=[1])
    
    N_experiments = length(S)

    fig = Figure()

    for i in Base.OneTo(N_experiments)
        g_src = fig[1, i] = GridLayout()
        g_dst = fig[2, i] = GridLayout()
        
        S_exp = @views S[i]
        D_exp = @views D[i]

        connection_density!(g_src[1,1], S_exp, S_mask; hemisphere=hemisphere[i])
        connection_density!(g_dst[1,1], D_exp, D_mask; hemisphere=hemisphere[i])

        Label(fig[0, i], text = "Experiment $i", fontsize=20)
        colsize!(fig.layout, i, Relative(1/(N_experiments+1)))
    end

    Colorbar(fig[1:2, N_experiments+1]; limits=(0, 1), ticks=0:0.1:1, label="Density")

    Label(fig[1,0], text = "$source [Source]", rotation=pi/2, fontsize=20)
    Label(fig[2,0], text = "$destination [Destination]", rotation=pi/2, fontsize=20)

    rowsize!(fig.layout, 1, Relative(1/2))
    rowsize!(fig.layout, 2, Relative(1/2))
    
    #rowgap!(fig.layout, 1, Relative(0.001))
    #colgap!(fig.layout, 2, Relative(0.001))

    fig
end

function plot_connection_density(source::String, destination::String; resolution=50, N_experiments=1)
    mouse_conn = pyimport("allensdk.core.mouse_connectivity_cache")
    mcc = mouse_conn.MouseConnectivityCache()
    struct_tree = mcc.get_structure_tree()

    lh = pyconvert(Dict, struct_tree.get_structures_by_acronym([source])[0])
    lh_experiments = pyconvert(AbstractArray{Dict}, mcc.get_experiments(injection_structure_ids=[lh["id"]]))
    lh_experiments = DataFrame(lh_experiments)

    vta = pyconvert(Dict, struct_tree.get_structures_by_acronym([destination])[0])

    pandas_df = mcc.get_structure_unionizes(
        lh_experiments.id, 
        is_injection=false,
        structure_ids=[vta["id"]],
        include_descendants=true
    )
    df = to_dataframe(pandas_df)
    sort!(df, :projection_density, rev=true)

    rsp = pyimport("allensdk.core.reference_space_cache")
    rspc = rsp.ReferenceSpaceCache(resolution, "annotation/ccf_2017", manifest="manifest.json")

    S_lh = get_structure_mask(rspc, lh["id"])
    S_vta = get_structure_mask(rspc, vta["id"])

    X = [download_projection_density(exp_id, "proj_dens", resolution) for exp_id in df.experiment_id[1:N_experiments]]

    Y = [download_injection_density(exp_id, "inj_dens", resolution) for exp_id in df.experiment_id[1:N_experiments]]

    fig = plot_connection_density(X, Y, S_lh, S_vta; source="LHb", destination="VTA", hemisphere=df.hemisphere_id[1:N_experiments])

    return fig
end