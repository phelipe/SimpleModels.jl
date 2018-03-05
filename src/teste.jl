push!(LOAD_PATH, expanduser("/home/phelipe/Documentos/programas/"))
using SimpleModels
using Trajectory
using Plots
using StaticArrays
using DifferentialEquations
pyplot()

type InputRobot{T} <: DEDataVector{T}
    x::Array{T,1}
    tau::Array{T,1}
end

robot = Dof2([23.902, 1.285], [.45, .45], [0.091, 0.048], [1.266, 0.093])
model = dinamic(robot)

Ts = 0.005
tend = 2.0
t0 = 0.0
x_0 = [0.,0.]
v_0 = [0.,0.]
x1,v1,a1,j1 = minimumjerkf(0.,0.,0.,0.,1.05,0.,0.,tend)
x2,v2,a2,j2 = minimumjerkf(0.,0.,0.,0.,1.57,0.,0.,tend)
kp = SMatrix{2,2}(diagm([7300.,600.]))
kv = SMatrix{2,2}(diagm([700.,20.]))

function controlador(integrator)
    data = user_cache(integrator)
    data = data[length(data)]
    #q = data.x[1:2]
    #dotq = data.x[3:4]
    xr = [x1(integrator.t), x2(integrator.t)]
    vr = [v1(integrator.t), v2(integrator.t)]
    #e = xr - q
    #de = vr - dotq
    #torque = (kp*e + kv*de)
    #println(integrator.sol)
    for c in user_cache(integrator)
        q = c.x[1:2]
        dotq = c.x[3:4]
        e = xr - q
        de = vr - dotq
        torque = (kp*e + kv*de)
        c.tau = torque
    end
end

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
    kp = SMatrix{2,2}(diagm([7300.,600.]))
    kv = SMatrix{2,2}(diagm([700.,20.]))
    tau = kp*e + kv*de

    du[1:2] = dθ
    du[3:4] = inv(M)*(tau - C * dθ - G)
end


cbs = PeriodicCallback(controlador,Ts)
tspan = (t0,tend)
start = InputRobot(vcat(x_0,v_0), zeros(2))
#prob = ODEProblem(model,start,tspan)
#prob = ODEProblem(myrobot2,start,tspan)

prob2 = ODEProblem(myrobot,[0.,0.,0.,0.],tspan)
sol2 = solve(prob2,Tsit5(),saveat = 0.08)
#sol = solve(prob,Tsit5(),saveat = 0.2,callback = cbs,force_dtmin=true)
#sol = solve(prob,Tsit5(),callback = cbs)#,force_dtmin=true)#,reltol=1e-3,abstol=1e-6,,saveat = Ts/4
#(x,v,t,a,ta,j,tj) = organize(robot,sol)
(x,v,t,a,ta,j,tj) = organize(2,sol2)
p1 = plot(j2,tj)
p2 = plot(tj,j[2])
plot(p1,p2)

#(x,v,t,a,ta,j,tj) = organize(2,sol2)
# plot(x1,t)
# plot(t,x[1])
#plot(p1,p2)
