push!(LOAD_PATH, expanduser("/home/phelipe/Documentos"))
using SimpleModels
using Trajectory
using Plots
using StaticArrays
using DifferentialEquations
pyplot()


#robot = Dof2([23.902, 1.285], [.45, .45], [0.091, 0.048], [1.266, 0.093])
#model = dinamic(robot)

Ts = 0.08
tend = 2.0
t0 = 0.0
x1,v1,a1,j1 = minimumjerkf(0.,0.,0.,t0,1.05,0.,0.,tend)
x2,v2,a2,j2 = minimumjerkf(0.,0.,0.,t0,1.57,0.,0.,tend)
#kp = SMatrix{2,2}(diagm([3300.,400.]));kv = SMatrix{2,2}(diagm([200.,20.]))
kp = SMatrix{2,2}(diagm([7000.,430.]));kv = SMatrix{2,2}(diagm([450.,30.]))


function myrobot(du, u, p, t)
    m=[23.902, 1.285]
    l= [.45, .45]
    r= [0.091, 0.048]
    I= [1.266, 0.093]
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

    xr = [x1(t), x2(t)]
    vr = [v1(t), v2(t)]
    e = xr - θ
    de = vr - dθ
    tau = kp*e + kv*de

    du[1:2] = dθ
    du[3:4] = inv(M)*(tau - C * dθ - G)
end

tspan = (t0,tend)
prob2 = ODEProblem(myrobot,[0.,0.,0.,0.],tspan)
sol2 = solve(prob2,Tsit5(),saveat = Ts)

x,v,t,a,ta,j,tj = organize(2,sol2)

xd1 = map(x->x1(x),t)
xd2 = map(x->x2(x),t)
vd1 = map(x->v1(x),t)
vd2 = map(x->v2(x),t)
ad1 = map(x->a1(x),ta)
ad2 = map(x->a2(x),ta)
jd1 = map(x->j1(x),tj)
jd2 = map(x->j2(x),tj)

function plotx()
    p1 = plot(t,[xd1,x[1]], xlabel="tempo",ylabel="posição")
    p2 = plot(t,[xd2,x[2]], xlabel="tempo",ylabel="posição")
    plot(p1,p2, tittle="Posições")
end

function plotex()
    p1 = plot(t,xd1 - x[1], xlabel="tempo",ylabel="erro posição")
    p2 = plot(t,xd2 - x[2], xlabel="tempo",ylabel="erro posição")
    plot(p1,p2, tittle="Posições")
end

function plotv()
    p1 = plot(t,[vd1,v[1]], xlabel="tempo",ylabel="velocidade")
    p2 = plot(t,[vd2,v[2]], xlabel="tempo",ylabel="velocidade")
    plot(p1,p2, tittle="Velocidades")
end

function plotev()
    p1 = plot(t,vd1 - v[1], xlabel="tempo",ylabel="erro velocidade")
    p2 = plot(t,vd2 - v[2], xlabel="tempo",ylabel="erro velocidade")
    plot(p1,p2, tittle="Velocidades")
end

function plota()
    p1 = plot(t,[ad1,a[1]], xlabel="tempo",ylabel="aceleração")
    p2 = plot(t,[ad2,a[2]], xlabel="tempo",ylabel="aceleração")
    plot(p1,p2, tittle="Acelerações")
end

function plotea()
    p1 = plot(t,ad1 - a[1], xlabel="tempo",ylabel="erro aceleração")
    p2 = plot(t,ad2 - a[2], xlabel="tempo",ylabel="erro aceleração")
    plot(p1,p2, tittle="Acelerações")
end

function plotj()
    p1 = plot(t,[jd1,j[1]], xlabel="tempo",ylabel="tranco")
    p2 = plot(t,[jd2,j[2]], xlabel="tempo",ylabel="tranco")
    plot(p1,p2, tittle="Tranco")
end

function plotej()
    p1 = plot(t,jd1 - j[1] , xlabel="tempo",ylabel="erro tranco")
    p2 = plot(t,jd2 - j[2], xlabel="tempo",ylabel="erro tranco")
    plot(p1,p2, tittle="Tranco")
end

function plotallx()
    plot(plotx(),plotex())
end

function plotallv()
    plot(plotv(),plotev())
end

function plotalla()
    plot(plota(),plotea())
end

function plotallj()
    plot(plotj(),plotej())
end


#gui()
#plotallv()
