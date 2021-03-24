module Hydrograph

import Plots
Plots.GRBackend()

function plot(v)
    flow_plot = Plots.plot(v.datetime, v.flow, title="flow rate (mm / day)", seriescolor=:black, label="")
    precip_plot = Plots.sticks(v.datetime, v.precip, title="precip (mm / day)", seriescolor=:black, label="")
    Plots.plot(flow_plot, precip_plot, layout = Plots.grid(2, 1, heights=[0.7, 0.3]))
end

end
