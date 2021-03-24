# reference in tree version of CalibrateEmulateSample
prepend!(LOAD_PATH, [joinpath(@__DIR__, "..", "..")])

include("read_mopex.jl")

function nash_cascade(inflow, out_flux, k)
    outflow = 0.0
    for i in eachindex(out_flux)
        out_flux[i] = ((1.0 - k) * out_flux[i]) + ((1.0 - k) * inflow)
        outflow = (k / (1.0 - k)) * out_flux[i]
        inflow = outflow
    end
    return outflow, out_flux
end


function hymod_simulate(forcings; C_max=1.0, B_exp=0.3, alpha=0.4, ks=0.001, kq=0.5, init_flow=false)
    nsteps = length(forcings.precip)
    res_fast = zeros(Float64, 3)
    res_slow = zeros(Float64, 1)
    if init_flow
        res_slow[1] = 2.3503 / (ks * 22.5)
    end

    # simulated fluxes
    channel_flow = zeros(Float64, nsteps)
    ground_flow = zeros(Float64, nsteps)
    surface_flow = zeros(Float64, nsteps)
    simulated_et = zeros(Float64, nsteps)

    # initial states
    soil_moist, ea = 0.0, 0.0
    for t in eachindex(forcings.precip)
        precip = forcings.precip[t]
        pet = forcings.pet[t]
        temp = forcings.temp[t]

        # update soil moisture state
        soil_moist = max(soil_moist + precip - pet - ea, 0.0)
    
        # ensure that soil moisture is not more than
        # C_max soil storage capacity
        precip_excess = max(soil_moist - C_max, 0.0)
        soil_moist -= precip_excess
    
        # effective precip routed though soil
        precip_eff = precip - precip_excess

        # update
        ea = min(min(soil_moist / C_max, 1.0) * pet, soil_moist)
        deff_ea = max(pet - ea, 0.0)
        pe = precip_eff * (1.0 - (1.0 - min(soil_moist / C_max, 1.0)) ^ B_exp)

        # extract deff evap from excess runoff
        precip_excess -= min(precip_excess, deff_ea)
        ea += min(precip_excess, deff_ea)

        # Route slow soil flow component with single nash cascade
        Us = (1.0 - alpha) * pe
        Qs, res_slow = nash_cascade(Us, res_slow, ks)

        # Route quick overland flow component with nash cascade
        Uq = precip_excess + (alpha * pe)
        Qq, res_fast = nash_cascade(Uq, res_fast, kq)
    
        simulated_et[t] = ea 
        ground_flow[t] = Qs
        surface_flow[t] = Qq
        channel_flow[t] = max(Qs + Qq, 0.0)
    end
    return (pet=simulated_et, ground_flow=ground_flow, surface_flow=surface_flow, channel_flow=channel_flow) 
end

function nse(est, actual)
  res = est .- actual
  return 1.0 - sum(res.^2) / sum((actual - mean(actual)).^2)
end

# calibration parameters priors
# cmax = uniform(1, 100-400) 
# alpha = Beta(2,3) || uniform(0.2, 0.99)
# Bexp = Gamma(2,6) || uniform(0.0, 2.0)
# ks = uniform(0.01, kq)
# kq = uniform(ks, 1.5)
#
function hymod_calibrate()
    # read in the mopex dataset
    v = read_mopex(url="file://" * joinpath(@__DIR__, "data", "03455000.dly"))
    res = hymod_simulate(v)
end

# plot results
#include("hydrograph.jl")
#Hydrograph.plot(v)
