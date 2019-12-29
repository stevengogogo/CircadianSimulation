using Interact
using Plotsa
include("MODEL_PerGoldbeter1996.jl")
gr()
# Define adjustable function with predefined parametets
function simulation_adjust(PER; vd=0.95, tspan2=100.0)
    # The parameterset is from Ingalls
    #SET VALUE
    u0 = [2,0,2.0,2.0,2.0]
    p = (
        vs=0.76, vm=0.65, vd= vd,
        ks=0.38, k1=1.9, k2=1.3,
        V1=3.2, V2=1.58, V3=5.0, V4=2.5,
        K1=1.0, K2=1.0, K3=1.0, K4=1.0,
        KI=1.0, Km1=0.5, Kd=0.2, n=4.0
        )


    tspan = (0.0, tspan2)

    #SOLVE
    prob = ODEProblem(PER, u0, tspan, values(p))
    sol = solve(prob)
    return sol
end

# Interact setting

## TIME COURSE
PER = get_PER_MODEL()
labels = ["mRNA(m)","PER(p0)","PER-P(p1)","PER-PP(p2)","PER-Nucleus(pn)"]
ui = @manipulate for t_end in 100.0:10.:300.0,
                        vd in slider(0.1:0.01:6., value=0.95, label="Vd (Decay rate of PER(uM/h))")
    sol = simulation_adjust(PER; vd=vd, tspan2=t_end)
    plot(sol, label=labels)
    ylabel!("Concentration(uM)")
    xlabel!("t(hr)")
end


## Phase Plot

ui = @manipulate for t_end in 100.0:10.:300.0,
                        vd in slider(0.1:0.01:6., value=0.95, label="Vd (Decay rate of PER(uM/h))"),
                        v1 in slider(1:1:5, value=1, label="u1"),
                        v2 in slider(1:1:5, value=2, label="u2"),
                        v3 in slider(1:1:5, value=3, label="u3")
    sol = simulation_adjust(PER; vd=vd, tspan2=t_end)
    plot(sol, vars=(v1,v2,v3), title="Phase Plot (vd = $vd)", xlabel=labels[v1], ylabel=labels[v2], zlabel=labels[v3])
end
