using AllenSDK
using PythonCall
using DataFrames

using HTTP
using FileIO

using GLMakie

function download_projection_density(id, filename, resolution)
    download("http://api.brain-map.org/grid_data/download_file/$(id)??image=projection_density&resolution=$(resolution)", filename * ".nrrd")

    return pyconvert(AbstractArray{Float64}, load("$(filename).nrrd"))
end

function download_injection_density(id, filename, resolution)
    download("http://api.brain-map.org/grid_data/download_file/$(id)??image=injection_density&resolution=$(resolution)", filename * ".nrrd")

    return pyconvert(AbstractArray{Float64}, load("$(filename).nrrd"))
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

function plot_connection_density(S, D, S_mask, D_mask; source="", destination="", hemisphere=2)
    
    fig = Figure()
    ga = fig[1, 1] = GridLayout()
    gb = fig[1, 2] = GridLayout()
    gc = fig[1, 3] = GridLayout()

    connection_density!(ga, S, S_mask; hemisphere, title="$source [Source]")
    connection_density!(gb, D, D_mask; hemisphere, title="$destination [Destination]")

    Colorbar(gc[1,1]; limits=(0, 1), ticks=0:0.1:1, label="Density")

    fig
end

resolution = 50

mouse_conn = pyimport("allensdk.core.mouse_connectivity_cache")
mcc = mouse_conn.MouseConnectivityCache()
struct_tree = mcc.get_structure_tree()

lh = pyconvert(Dict, struct_tree.get_structures_by_acronym(["LH"])[0])
lh_experiments = pyconvert(AbstractArray{Dict}, mcc.get_experiments(injection_structure_ids=[lh["id"]]))
lh_experiments = DataFrame(lh_experiments)

vta = pyconvert(Dict, struct_tree.get_structures_by_acronym(["VTA"])[0])

pandas_df = mcc.get_structure_unionizes(
    lh_experiments.id, 
    is_injection=false,
    structure_ids=[vta["id"]],
    include_descendants=true
)

experiment_ids = pyconvert(Vector{Int}, pandas_df["experiment_id"])
hemisphere_ids = pyconvert(Vector{Int}, pandas_df["hemisphere_id"])
vals_dens = pyconvert(Vector{Float64}, pandas_df["projection_density"])

idx = argmax(vals_dens)

rsp = pyimport("allensdk.core.reference_space_cache")
rspc = rsp.ReferenceSpaceCache(resolution, "annotation/ccf_2017", manifest="manifest.json")

S_lh = get_structure_mask(rspc, lh["id"])
S_vta = get_structure_mask(rspc, vta["id"])

X = download_projection_density(experiment_id[idx], "proj_dens", resolution)

Y = download_injection_density(experiment_id[idx], "inj_dens", resolution)

plot_connection_density(X.data, Y.data, S_lh, S_vta; source="LHb", destination="VTA", hemisphere=hemisphere_ids[idx])
