using ModelingToolkit
using DifferentialEquations
#CHEMICAL REACTIONS
mm(v, k, sub) = v*sub/(k+sub)
hill(v,p,k,n) = v/(1+(p/k)^n)


function get_PER_MODEL()
    # labels = ["mRNA(m)" "PER(p0)" "PER-P(p1)" "PER-PP(p2)" "PER-Nucleus(pn)"]
    # From Goldbeter, 1996
    #PARAMETERS
    @parameters t
    @parameters vs vm vd
    @parameters ks k1 k2
    @parameters V1 V2 V3 V4
    @parameters K1 K2 K3 K4 KI Km1 Kd n
    @variables m(t) p0(t) p1(t) p2(t) pn(t)
    @derivatives D'~t

    eqs = [
        D(m)  ~ hill(vs, pn, KI, n) - mm(vm, Km1, m), #mRNA(m)
        D(p0) ~ ks*m - mm(V1, K1,p0) + mm(V2, K2, p1), #PER(p0)
        D(p1) ~ mm(V1, K1, p0) - mm(V2, K2,p1) - mm(V3, K3, p1) + mm(V4, K4, p2), #PER-P(p1)
        D(p2) ~ mm(V3, K3, p1) - mm(V4, K4, p2) - k1*p2 + k2*pn - mm(vd,Kd,p2), #PER-PP(p2)
        D(pn) ~ k1*p2-k2*pn #PER-Nucleus(pn)
    ]
    de = ODESystem(eqs)
    PER = ODEFunction(de,
                    [m,p0, p1, p2, pn],
                    [vs,vm,vd,ks,k1,k2,V1,V2,V3,V4,K1,K2,K3,K4,KI,Km1,Kd,n]
                    )

    return PER
end


function get_param(;vd=0.95, namedtuple=false)
    # From Goldbeter, 1996
    u0 = [2,0,2.0,2.0,2.0]
    p = (
        vs=0.76, vm=0.65, vd=vd,
        ks=0.38, k1=1.9, k2=1.3,
        V1=3.2, V2=1.58, V3=5.0, V4=2.5,
        K1=1.0, K2=1.0, K3=1.0, K4=1.0,
        KI=1.0, Km1=0.5, Kd=0.2, n=4.0
        )
    if namedtuple == false
        p = values(p)
    end
    return u0, p
end

function solve_TimeSeries(PER; vd=0.95)
    # The parameterset is from Ingalls
    #SET VALUE
    u0, p = get_param(vd=vd)
    tspan = (0.0, 100.0)

    #SOLVE
    prob = ODEProblem(PER, u0, tspan, p)
    sol = solve(prob)
    return sol
end

function solve_steady_state(PER;vd=0.95)
    u0, p = get_param(vd=vd)
    prob = SteadyStateProblem(PER, u0, p)
    cp = solve(prob, SSRootfind())
    return cp
end

function simulation(;vd=0.95)
    labels = ["mRNA(m)" "PER(p0)" "PER-P(p1)" "PER-PP(p2)" "PER-Nucleus(pn)"]
    PER= get_PER_MODEL()
    sol = simulation(PER; vd=vd)
    cp = solve_steady_state(PER; vd=vd)
    return labels, sol, cp
end


# MAIN
#simulation
function main()
    labels, sol, cp = simulation(vd=0.95)
    # Static Plotting
    plot_num = (1, 4, 5)

    # Phase Plot with critical point
    phase = plot(sol, vars=plot_num, title="Phase Plot (vd=$vd)", xlabel=labels[plot_num[1]], ylabel=labels[plot_num[2]], zlabel=labels[plot_num[3]])
    plot!(phase, [cp[plot_num[1]]], [cp[plot_num[2]]], [cp[plot_num[3]]],markershape=:circle, markersize=10, label="critical point")

    #Plot mRna, p2,pn
    time = plot(sol, title="PER Oscillator (vd=$vd)", vars=collect(plot_num) , labels= labels[ collect(plot_num)' ])
    ylabel!(time, "Concentration(uM)")
    xlabel!(time, "t(hr)")

    # savefig("img/phasePlot_criticalpoint_$vd")
    # savefig("img/timeSeries_three_$vd.png")
    display(phase)
    display(time)
end
######
if abspath(PROGRAM_FILE) == @__FILE__
    using Plots
    gr(dpi=400)
    main()
end
