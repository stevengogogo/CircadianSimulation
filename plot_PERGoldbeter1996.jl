include("MODEL_PerGoldbeter1996.jl")
using Plots
gr(dpi=400)

# SETTING
function setting(;vd=0.95)
    #= Get parameters, labels and solutions =#
    _, p = get_param(vd=vd)
    labels, sol, cp = simulation(vd=vd)
    return p, labels, sol, cp
end


#= STATIC PLOTS =#

function PhasePlot_criticalpoint(;vd=vd, cp=cp, plot_num=plot_num, sol=sol, labels=labels)
    #= Phase Plot with critical point =#
    g = plot(sol, vars=plot_num, title="Phase Plot (vd=$vd)", xlabel=labels[plot_num[1]], ylabel=labels[plot_num[2]], zlabel=labels[plot_num[3]])

    plot!(g, [cp[plot_num[1]]], [cp[plot_num[2]]], [cp[plot_num[3]]],markershape=:circle, markersize=10, label="critical point")

    return g
end



function PlotTimeSeries(;vd=vd, plot_num=plot_num, sol=sol, labels=labels)
    #= Time-series Plot mRna, p2,pn =#
    g = plot(sol, title="PER Oscillator (vd=$vd)", vars=collect(plot_num) , labels= labels[ collect(plot_num)' ])
    ylabel!(g, "Concentration(uM)")
    xlabel!(g, "t(hr)")
    return g
end


#= ANIMATION =#
function ANIM_TimeSeries(;vd_span=0.1:0.1:3)
    #=Animate time series with multiple vd values =#
    animTimeSeries = @animate for vd=vd_span
        sol = simulation(PER; vd=vd)
        plot(sol, title="PER Oscillator (vd=$vd)", vars=collect(plot_num) , labels= labels[ collect(plot_num)' ])
        ylabel!("Concentration(uM)")
        xlabel!("t(hr)")
    end
    return animTimeSeries
end

function ANIM_Phase(;vd_span=0.1:0.1:3)
    #=Animate phase plot with multiple vd values =#
    animPhase = @animate for vd=vd_span
        sol = simulation(PER; vd=vd)
        cp = solve_steady_state(PER; vd=vd)
        plot(sol, vars=plot_num, title="Phase Plot (vd=$vd)", xlabel=labels[plot_num[1]], ylabel=labels[plot_num[2]], zlabel=labels[plot_num[3]])
        plot!([cp[plot_num[1]]], [cp[plot_num[2]]], [cp[plot_num[3]]],markershape=:circle, markersize=10, label="critical point")
    end
    return animPhase
end



### Tutorial of Negative feedback
function animTimeSeries_moving(;sol=sol, p=p, plot_num = (1, 5),t_tick=30)
    #= Plot time-series with moving labels=#

    animTime_movingCircle = @animate for i= 1:t_tick:length(sol.t)
        g = plot(sol, title="PER Oscillator (vd=$vd)", vars=collect(plot_num) , labels= labels[ collect(plot_num)' ])

        for k in plot_num
            plot!(g, [sol.t[i]], [sol.u[i][k]] ,markershape=:circle, markersize=10, labels= labels[k])
        end

        ylabel!("Concentration(uM)")
        xlabel!("t(hr)")
    end
    return animTime_movingCircle
end


#### With differentiation of mRNA
function animTimeSeries_movingD(;t_tick=30, labels=labels, sol=sol)
    #= Plot time-series with moving labels with With differentiation of mRNA=#
    plot_num = (1, 5)
    u, p = get_param(;namedtuple=true)
    d(pn,m;p=p)  =  hill(p.vs, pn, p.KI, p.n) - mm(p.vm, p.Km1, m)

    name = Dict("input"=> labels[plot_num[2]],"output"=>labels[plot_num[1]])

    animTime_movingCircleD = @animate for i= 1:t_tick:length(sol.t)
        g = plot(sol, title="PER Oscillator (vd=$vd)", vars=collect(plot_num) , labels= labels[ collect(plot_num)' ])


        dg = plot(sol.t, d.([i[plot_num[end]] for i in sol.u],
                            [i[plot_num[1]] for i in sol.u]), linewidth=3, labels=string("Rate of ", name["output"]),xlims = (0,sol.t[end]))

        for k in plot_num
            plot!(g, [sol.t[i]], [sol.u[i][k]] ,markershape=:circle, markersize=10, labels= labels[k])
        end


        hline!(dg, [0.0], labels=["rate=0"])
        plot!(dg, [sol.t[i]], [d(sol.u[i][plot_num[2]],sol.u[i][plot_num[1]])], markershape=:circle, markersize=10, labels=name["output"])


        ylabel!(g, "Concentration(uM)")
        xlabel!(g, "t(hr)")
        ylabel!(dg, "rate(uM/hr)")
        xlabel!(dg, "t(hr)")
        plot(g, dg, layout=(2,1))
    end
    return animTime_movingCircleD
end


# draw differential of RNA: without mRNA time series
function animTimeSeries_movingD_without_mRNA(;t_tick=30, labels=labels, sol=sol)
    #= Plot time-series with moving labels with With differentiation of mRNA, and remove mRNA time-series for clearance=#
    plot_num = (1,)
    name = Dict("input"=> 5,"output"=>1)
    u, p = get_param(;namedtuple=true)

    # rate of mRNA production
    d(pn,m;p=p)  =  hill(p.vs, pn, p.KI, p.n) - mm(p.vm, p.Km1, m)


    animTime_movingCircleD_mRNA_removed = @animate for i= 1:t_tick:length(sol.t)
        g = plot(sol, title="PER Oscillator (vd=$vd)", vars=collect([name["input"]]) , labels= labels[ collect([name["input"]])' ], color="red")


        dg = plot(sol.t, d.([i[name["input"]] for i in sol.u],
                            [i[name["output"]] for i in sol.u]), linewidth=3, labels=string("Rate of ", labels[name["output"]]), xlims = (0,sol.t[end]))

        for k in [name["input"]]
            plot!(g, [sol.t[i]], [sol.u[i][k]] ,markershape=:circle, markersize=10, labels= labels[k])
        end


        hline!(dg, [0.0], labels=["rate=0"]) #mark dv/dt = 0
        plot!(dg, [sol.t[i]], [d(sol.u[i][name["input"]],sol.u[i][name["output"]])], markershape=:circle, markersize=10, labels=labels[name["output"]])


        ylabel!(g, "Concentration(uM)")
        xlabel!(g, "t(hr)")
        ylabel!(dg, "rate(uM/hr)")
        xlabel!(dg, "t(hr)")

        plot(g, dg, layout=(2,1))
    end
    return animTime_movingCircleD_mRNA_removed
end


function ANIM_phaseplot(;vd = vd, t_tick=30, labels=labels, sol=sol, plot_num=plot_num)
    #= Display what phase plot is by animation =#
    function PlotTimeSeries_nolegend(;vd=vd, plot_num=plot_num, sol=sol, labels=labels)
        #= Time-series Plot mRna, p2,pn =#
        g = plot(sol, title="PER Oscillator (vd=$vd)", vars=collect(plot_num) , labels= labels[ collect(plot_num)' ], legend=nothing)
        ylabel!(g, "Concentration(uM)")
        xlabel!(g, "t(hr)")
        return g
    end
    animatephase = @animate for i= 1:t_tick:length(sol.t)
        g_ph = PhasePlot_criticalpoint(;vd=vd, plot_num=plot_num, sol=sol, labels=labels)
        g_time = PlotTimeSeries_nolegend()

        for k in plot_num
            plot!(g_time, [sol.t[i]], [sol.u[i][k]] ,markershape=:circle, markersize=10, labels= labels[k])
        end

        plot!(g_ph, [sol.u[i][plot_num[1]]], [sol.u[i][plot_num[2]]],[sol.u[i][plot_num[3]]] ,markershape=:circle, markersize=10, labels= labels[collect(plot_num)'])

        plot(g_time, g_ph, layout=(1,2))
    end
    return animatephase
end



## MAIN

function main()
    p, labels, sol, cp=  setting(vd=0.95)
    plot_num = (1, 4, 5) # Plot index

    # Static Plotting
    phaseplot_criticalpoint = PhasePlot_criticalpoint()
    savefig("img/phasePlot_criticalpoint_$vd")

    plottimeseries = PlotTimeSeries()
    savefig("img/timeSeries_three_$vd.png")

    # Animation
    animTimeSeries = ANIM_TimeSeries()
    gif(animTimeSeries, "img/anim_timeSeries.gif", fps = 8)


    animPhase = ANIM_Phase()
    gif(animPhase, "img/anim_phasePlot.gif", fps = 8)

    animTime_movingCircle = animTimeSeries_moving(plot_num = (1,5),t_tick=1)
    gif(animTime_movingCircle, "img/animTime_movingCircle.gif", fps = 8)

    animTime_movingCircle_p012 = animTimeSeries_moving(plot_num = (2,3,4),t_tick=1)
    gif(animTime_movingCircle_p012, "img/animTime_movingCircle_p012.gif",fps=8)

    animTime_movingCircleD = animTimeSeries_movingD()
    gif(animTime_movingCircleD, "img/animTime_movingCircleD.gif", fps = 8)

    animTime_movingCircleD_mRNA_removed = animTimeSeries_movingD_without_mRNA()
    gif(animTime_movingCircleD_without_mRNA, "img/animTime_movingCircleD_without_mRNA.gif", fps = 8)

    animatephase = ANIM_phaseplot(t_tick=1)
    gif(animatephase,"img/animatephase.gif", fps = 8 )
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
