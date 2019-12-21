using Interact
@parameters t
@parameters vs vm vd
@parameters ks k1 k2
@parameters V1 V2 V3 V4
@parameters K1 K2 K3 K4 KI Km1 Kd n
@parameters a, ω, ϕ, b
@variables m(t) p0(t) p1(t) p2(t) pn(t)
@derivatives D'~t


mm(v, k, sub) = v*sub/(k+sub)
hill(v,p,k,n) = v/(1+(p/k)^n)

SIN(a, ω,t,ϕ,b) = a*sin(ω*t + ϕ) + a +b


eqs = [
    D(m)  ~ hill(vs, pn, KI, n) - mm(vm, Km1, m),
    D(p0) ~ ks*m - mm(V1, K1,p0) + mm(V2, K2, p1),
    D(p1) ~ mm(V1, K1, p0) - mm(V2, K2,p1) - mm(V3, K3, p1) + mm(V4, K4, p2),
    D(p2) ~ mm(V3, K3, p1) - mm(V4, K4, p2) - k1*p2 + k2*pn - mm(SIN(a,ω,t,ϕ,b) ,Kd,p2),
    D(pn) ~ k1*p2-k2*pn
]

de = ODESystem(eqs)

PER = ODEFunction(de,
                [m,p0, p1, p2, pn],
                [b, a, ω, ϕ, vs,vm,ks,k1,k2,V1,V2,V3,V4,K1,K2,K3,K4,KI,Km1,Kd,n]
                )



u0 = [2,0,2.0,2.0,2.0]
p = (
    b= 0.95, a=1., ω=2. *pi/24, ϕ=0.,
    vs=0.76, vm=0.65,
    ks=0.38, k1=1.9, k2=1.3,
    V1=3.2, V2=1.58, V3=5.0, V4=2.5,
    K1=1.0, K2=1.0, K3=1.0, K4=1.0,
    KI=1.0, Km1=0.5, Kd=0.2, n=4.0
    )

tspan = (0.0, 100.0)

#SOLVE
prob = ODEProblem(PER, u0, tspan, values(p))
sol = solve(prob)


ui = @manipulate for amp in slider(0.:0.01:1., value=0.0, label="amplitude"),
                     omega in slider((2. *pi/48.):0.01:(2. *pi/12.), value=2. *pi/24.0, label="angular frequency"),
                     ba in slider(0.6:0.01:2.9, value=0.95, label="vd_baseline"),
                     phi in slider(0.:0.1:(2.0 * pi), value=0.0, label="phase")
     p = (
         b= ba, a= amp, ω=omega, ϕ=phi,
         vs=0.76, vm=0.65,
         ks=0.38, k1=1.9, k2=1.3,
         V1=3.2, V2=1.58, V3=5.0, V4=2.5,
         K1=1.0, K2=1.0, K3=1.0, K4=1.0,
         KI=1.0, Km1=0.5, Kd=0.2, n=4.0
         )
    prob = ODEProblem(PER, u0, tspan, values(p))
    sol = solve(prob)
    plot(sol, vars=[1])
    plot!(sol.t , [SIN(amp, omega, i, phi, ba) for i in sol.t ] )
end
