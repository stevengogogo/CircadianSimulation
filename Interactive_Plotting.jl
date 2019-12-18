using Interact
using Plots
include("MODEL_PerGoldbeter1996.jl")

# Define adjustable function with predefined parametets
function simulation_Adjust(PER;vs=0.76, vm=0.65, vd=0.95,ks=0.38, k1=1.9, k2=1.3,V1=3.2, V2=1.58, V3=5.0, V4=2.5,K1=1.0, K2=1.0, K3=1.0, K4=1.0,KI=1.0, Km1=0.5, Kd=0.2, n=4.0, u01=2.0, u02=2.0, u03=2.0, u04=2.0, tspan1=0.0, tspan2=100.0)



    tspan = (tspan1, tspan2)
    u0 = [u01, u02, u03, u04]
    p = (vs, vm, vd,
        ks, k1, k2,
        V1, V2, V3, V4,
        K1, K2, K3, K4,
        KI, Km1, Kd, n)

    prob = ODEProblem(PER, u0, tspan, p)
    sol= solve(prob)
    return sol
end

# Interact setting
ui = button()
display(ui)

PER, labels = get_PER_MODEL()

ui = @manipulate for vd in 0.2:0.01:2.0
    vd = Float64(vd)
    sol = simulation_Adjust(PER; vd=vd)
    plot(sol)
end
