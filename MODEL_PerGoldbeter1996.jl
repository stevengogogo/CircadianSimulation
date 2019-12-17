using ModelingToolkit
using DifferentialEquations
using Plots

@parameters t
@parameters vs vm vd
@parameters ks k1 k2
@parameters V1 V2 V3 V4
@parameters K1 K2 K3 K4 KI Km1 Kd n
@variables m(t) p0(t) p1(t) p2(t) pn(t)
@derivatives D'~t

mm(v, k, sub) = v*sub/(k+sub)
hill(v,p,k,n) = v/(1+(p/k)^n)


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

u0 = [2,0,2.0,2.0,2.0]


p = ( 
    vs=0.76, vm=0.65, vd=0.95,
    ks=0.38, k1=1.9, k2=1.3,
    V1=3.2, V2=1.58, V3=5.0, V4=2.5,
    K1=1.0, K2=1.0, K3=1.0, K4=1.0,
    KI=1.0, Km1=0.5, Kd=0.2, n=4.0
    )

tspan = (0.0, 100.0)

prob = ODEProblem(PER, u0, tspan, values(p))
sol = solve(prob)

plot(sol, vars=(1,2,3), lw=1)
