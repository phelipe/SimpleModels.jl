using Plots
using StaticArrays
using DifferentialEquations
pyplot()

function simulation(kp::T, kv::T, Ts::Y, t0::Y, tend::Y, r1::Y , r2::Y) where {T<:AbstractMatrix, Y<: AbstractFloat}

    function myrobot(du, u, p, t)
        m = SVector{2}([23.902, 1.285])
        l= SVector{2}([.45, .45])
        r= SVector{2}([0.091, 0.048])
        I= SVector{2}([1.266, 0.093])
        g = 9.81
        θ = SVector{2}(u[1:2])
        dθ = SVector{2}(u[3:4])
        const cos1 = cos(θ[1])
        const cos2 = cos(θ[2])
        const sin1 = sin(θ[1])
        const sin2 = sin(θ[2])
        const cos12 = cos(θ[1]+θ[2])

        #Matriz de inércia
        m11 = m[1] * r[1]^2 + m[2]* (l[1]^2 + r[2]^2 + 2*l[1] * r[2] * cos2)+ I[1] + I[2]
        m12 = m[2] * (r[2]^2 + l[1] *r[2] * cos2) + I[2]
        m22 = m[2] * r[2]^2 + I[2]
        M = SMatrix{2,2}([m11 m12; m12 m22])

        #Matriz C
        h = -m[2] * l[1] * r[2] * sin2
        c11 = h * dθ[2]
        c12 = h * (dθ[1] + dθ[2])
        c21 = -h * dθ[1]
        c22 = 0
        C = SMatrix{2,2}([c11 c12; c21 c22])

        #Matriz gravidade
        g2 = m[2] * r[2] * g * cos12
        g1 = (m[1] * r[1] + m[2] * l[1]) * g * cos1 + g2
        G = SVector{2}([g1; g2])

        xr = SVector{2}([r1, r2])
        e = xr - θ
        tau = kp*e - kv*dθ
        du[1:2] = dθ
        du[3:4] = inv(M)*(tau - C * dθ - G)
    end

    tspan = (t0,tend)
    prob2 = ODEProblem(myrobot,[0., 0., 0., 0.], tspan)
    sol2 = solve(prob2, Tsit5(), saveat = Ts,force_dtmin=true, maxiters = 1e10)

    organize(2,sol2)
end

Ts = 0.08 # Intervalo entre leituras da saída
tend = 20.0 # tempo de simulação
t0 = 0.0 # instante inicial
r1 =  0.6#0.6#1.6# referência junta 1
r2 =  2.0#1.8#2.0# referência junta 2


kp_vec = [3500.,60.]
kv_vec = [200.,35.]
kp = SMatrix{2,2}(diagm(kp_vec))
kv = SMatrix{2,2}(diagm(kv_vec))
x, v, t, a, ta, j, tj = simulation(kp, kv, Ts, t0, tend, r1, r2)

p1 = plot(t,x[1], label = "desejado 1")
plot!([r1],seriestype= :hline, label = "referência")
p2 = plot(t,x[2], label = "desejado 2")
plot!([r2],seriestype= :hline, label = "referência")

plot(p1,p2)
