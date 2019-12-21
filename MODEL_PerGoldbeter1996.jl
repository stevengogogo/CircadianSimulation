using ModelingToolkit
using DifferentialEquations
using Plots
plotly()

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
        D(m)  ~ hill(vs, pn, KI, n) - mm(vm, Km1, m),
        D(p0) ~ ks*m - mm(V1, K1,p0) + mm(V2, K2, p1),
        D(p1) ~ mm(V1, K1, p0) - mm(V2, K2,p1) - mm(V3, K3, p1) + mm(V4, K4, p2),
        D(p2) ~ mm(V3, K3, p1) - mm(V4, K4, p2) - k1*p2 + k2*pn - mm(vd,Kd,p2),
        D(pn) ~ k1*p2-k2*pn
    ]
    de = ODESystem(eqs)
    PER = ODEFunction(de,
                    [m,p0, p1, p2, pn],
                    [vs,vm,vd,ks,k1,k2,V1,V2,V3,V4,K1,K2,K3,K4,KI,Km1,Kd,n]
                    )

    return PER
end


function get_param(;vd=0.95)
    # From Goldbeter, 1996
    u0 = [2,0,2.0,2.0,2.0]
    p = (
        vs=0.76, vm=0.65, vd=vd,
        ks=0.38, k1=1.9, k2=1.3,
        V1=3.2, V2=1.58, V3=5.0, V4=2.5,
        K1=1.0, K2=1.0, K3=1.0, K4=1.0,
        KI=1.0, Km1=0.5, Kd=0.2, n=4.0
        )
    return u0, values(p)
end

function simulation(PER; vd=0.95)
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


# MAIN
#Simulation

vd = 1.9
PER= get_PER_MODEL()
sol = simulation(PER; vd=vd)
cp = solve_steady_state(PER; vd=vd)


# Plotting

labels = ["mRNA(m)" "PER(p0)" "PER-P(p1)" "PER-PP(p2)" "PER-Nucleus(pn)"]
plot_num = (1, 4, 5)
gr()

# Phase Plot with critical point
plot(sol, vars=plot_num, title="Phase Plot (vd=$vd)", xlabel=labels[plot_num[1]], ylabel=labels[plot_num[2]], zlabel=labels[plot_num[3]])
plot!([cp[plot_num[1]]], [cp[plot_num[2]]], [cp[plot_num[3]]],markershape=:circle, markersize=10, label="critical point")
#savefig("img/phasePlot_criticalpoint_$vd")



#Plot mRna, p2,pn
plot(sol, title="PER Oscillator (vd=$vd)", vars=collect(plot_num) , labels= labels[ collect(plot_num)' ])
ylabel!("Concentration(uM)")
xlabel!("t(hr)")
#savefig("img/timeSeries_three_$vd.png")


# Gif plot
animTimeSeries = @animate for vd=0.1:0.1:3
    sol = simulation(PER; vd=vd)
    plot(sol, title="PER Oscillator (vd=$vd)", vars=collect(plot_num) , labels= labels[ collect(plot_num)' ])
    ylabel!("Concentration(uM)")
    xlabel!("t(hr)")
end

animPhase = @animate for vd=0.1:0.1:3
    sol = simulation(PER; vd=vd)
    cp = solve_steady_state(PER; vd=vd)
    plot(sol, vars=plot_num, title="Phase Plot (vd=$vd)", xlabel=labels[plot_num[1]], ylabel=labels[plot_num[2]], zlabel=labels[plot_num[3]])
    plot!([cp[plot_num[1]]], [cp[plot_num[2]]], [cp[plot_num[3]]],markershape=:circle, markersize=10, label="critical point")
end

gif(animTimeSeries, "img/anim_timeSeries.gif", fps = 8)
gif(animPhase, "img/anim_phasePlot.gif", fps = 8)
