module SCModel

using NCDatasets
using Statistics
using Interpolations

export run_SCAMPy
export obs_LES
export get_profile

function run_SCAMPy(u::Array{FT, 1},
                    u_names::Array{String, 1},
                    y_names::Union{Array{String, 1}, Array{Array{String,1},1}},
                    scm_dir::String,
                    ti::Union{FT, Array{FT,1}},
                    tf::Union{FT, Array{FT,1}},
                    ) where {FT<:AbstractFloat}
    
    exe_path = string(scm_dir, "call_SCAMPy.sh")
    sim_uuid  = u[1]
    for i in 2:length(u_names)
        sim_uuid = string(sim_uuid,u[i])
    end
    command = `bash $exe_path $u $u_names`
    run(command)

    # SCAMPy file descriptor
    sim_uuid = string(sim_uuid, ".txt")
    sim_dirs = readlines(sim_uuid)
    run(`rm $sim_uuid`)

    y_scm = zeros(0)
    for i in 1:length(sim_dirs)
        sim_dir = sim_dirs[i]
        if length(ti) > 1
            ti_ = ti[i]
            tf_ = tf[i]
        else
            ti_ = ti
            tf_ = tf
        end

        if typeof(y_names)==Array{Array{String,1},1}
            y_names_ = y_names[i]
        else
            y_names_ = y_names
        end
        append!(y_scm, get_profile(sim_dir, y_names_, ti = ti_, tf = tf_))
        run(`rm -r $sim_dir`)
    end

    for i in eachindex(y_scm)
        if isnan(y_scm[i])
            y_scm[i] = 1.0e4
        end
    end

    return y_scm
end

function obs_LES(y_names::Array{String, 1},
                    sim_dir::String,
                    ti::Float64,
                    tf::Float64;
                    z_scm::Union{Array{Float64, 1}, Nothing} = nothing,
                    ) where {FT<:AbstractFloat}
    
    y_names_les = get_les_names(y_names, sim_dir)
    y_highres = get_profile(sim_dir, y_names_les, ti = ti, tf = tf)
    y_tvar_highres = get_timevar_profile(sim_dir, y_names_les, ti = ti, tf = tf)
    if !isnothing(z_scm)
        y_ = zeros(0)
        y_tvar = zeros(0)
        z_les = get_profile(sim_dir, ["z_half"])
        num_outputs = Integer(length(y_highres)/length(z_les))
        for i in 1:num_outputs
            y_itp = interpolate( (z_les,), 
                y_highres[1 + length(z_les)*(i-1) : i*length(z_les)],
                Gridded(Linear()) )
            append!(y_, y_itp(z_scm))

            y_tvar_itp = interpolate( (z_les,), 
                y_tvar_highres[1 + length(z_les)*(i-1) : i*length(z_les)],
                Gridded(Linear()) )
            append!(y_tvar, y_tvar_itp(z_scm))
        end
    else
        y_ = y_highres
        y_tvar = y_tvar_highres
    end
    return y_, y_tvar
end

function get_profile(sim_dir::String,
                     var_name::Array{String,1};
                     ti::Float64=0.0,
                     tf::Float64=0.0,
                     getFullHeights=false)

    if length(var_name) == 1 && occursin("z_half", var_name[1])
        prof_vec = nc_fetch(sim_dir, "profiles", var_name[1])
    else
        t = nc_fetch(sim_dir, "timeseries", "t")
        dt = t[2]-t[1]
        ti_diff, ti_index = findmin( broadcast(abs, t.-ti) )
        tf_diff, tf_index = findmin( broadcast(abs, t.-tf) )
        
        prof_vec = zeros(0)
        # If simulation does not contain values for ti or tf, return high value
        if ti_diff > dt || tf_diff > dt
            for i in 1:length(var_name)
                var_ = nc_fetch(sim_dir, "profiles", var_name[i])
                append!(prof_vec, 1.0e4*ones(length(var_[:, 1])))
            end
        else
            for i in 1:length(var_name)
                var_ = nc_fetch(sim_dir, "profiles", var_name[i])
                # LES vertical fluxes are per volume, not mass
                if occursin("resolved_z_flux", var_name[i])
                    rho_half=nc_fetch(sim_dir, "reference", "rho0_half")
                    var_ = var_.*rho_half
                end
                append!(prof_vec, mean(var_[:, ti_index:tf_index], dims=2))
            end
        end
    end
    return prof_vec 
end

function get_timevar_profile(sim_dir::String,
                     var_name::Array{String,1};
                     ti::Float64=0.0,
                     tf::Float64=0.0,
                     getFullHeights=false)

    t = nc_fetch(sim_dir, "timeseries", "t")
    dt = t[2]-t[1]
    ti_diff, ti_index = findmin( broadcast(abs, t.-ti) )
    tf_diff, tf_index = findmin( broadcast(abs, t.-tf) )
    
    prof_vec = zeros(0)
    # If simulation does not contain values for ti or tf, return high value
    for i in 1:length(var_name)
        var_ = nc_fetch(sim_dir, "profiles", var_name[i])
        # LES vertical fluxes are per volume, not mass
        if occursin("resolved_z_flux", var_name[i])
            rho_half=nc_fetch(sim_dir, "reference", "rho0_half")
            var_ = var_.*rho_half
        end
        # append!(prof_vec, var(var_[:, ti_index:tf_index], dims=2) )
        append!(prof_vec, maximum(var(var_[:, ti_index:tf_index], dims=2))
            *ones(length(var_[:, 1]) ) )
    end

    return prof_vec
end

function get_les_names(scm_y_names::Array{String,1}, sim_dir::String)
    y_names = deepcopy(scm_y_names)
    if "thetal_mean" in y_names
        if occursin("GABLS",sim_dir) || occursin("Soares",sim_dir)
            y_names[findall(x->x=="thetal_mean", y_names)] .= "theta_mean"
        else
            y_names[findall(x->x=="thetal_mean", y_names)] .= "thetali_mean"
        end
    end
    if "total_flux_qt" in y_names
        y_names[findall(x->x=="total_flux_qt", y_names)] .= "resolved_z_flux_qt"
    end
    if "total_flux_h" in y_names && (occursin("GABLS",sim_dir) || occursin("Soares",sim_dir))
        y_names[findall(x->x=="total_flux_h", y_names)] .= "resolved_z_flux_theta"
    elseif "total_flux_h" in y_names
        y_names[findall(x->x=="total_flux_h", y_names)] .= "resolved_z_flux_thetali"
    end
    if "u_mean" in y_names
        y_names[findall(x->x=="u_mean", y_names)] .= "u_translational_mean"
    end
    if "v_mean" in y_names
        y_names[findall(x->x=="v_mean", y_names)] .= "v_translational_mean"
    end
    if "tke_mean" in y_names
        y_names[findall(x->x=="tke_mean", y_names)] .= "tke_nd_mean"
    end
    return y_names
end

function nc_fetch(dir, nc_group, var_name)
    find_prev_to_name(x) = occursin("Output", x)
    split_dir = split(dir, ".")
    sim_name = split_dir[findall(find_prev_to_name, split_dir)[1]+1]
    ds = NCDataset(string(dir, "/stats/Stats.", sim_name, ".nc"))
    ds_group = ds.group[nc_group]
    ds_var = deepcopy( Array(ds_group[var_name]) )
    close(ds)
    return Array(ds_var)
end

end #module