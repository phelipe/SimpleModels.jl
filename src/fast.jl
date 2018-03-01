push!(LOAD_PATH, expanduser("/home/phelipe/Documentos"))
using SimpleModels
using DifferentialEquations

robot = Dof2([23.902, 1.285], [.45, .45], [0.091, 0.048], [1.266, 0.093])
model = dinamic(robot)

type InputRobot{T} <: DEDataVector{T}
    x::Array{T,1}
    tau::Array{T,1}
end

function controlador(integrator)
    for c in user_cache(integrator)
        #q = c.x[1:2]
        #dotq = c.x[3:4]
        #e = r - q[2]
        #println(integrator.t)
        #println(e)
        c.tau = [0.,0.]#[(kp*e - kv*dotq[2]), 0]
    end
end

Ts = 0.05
tend = 5.
t0 = 0.
x_0 = [0.,0.]
v_0 = [0.,0.]

cbs = PeriodicCallback(controlador,Ts)
tspan = (t0,tend)
start = InputRobot(vcat(x_0,v_0), zeros(2))
prob = ODEProblem(model,start,tspan)
sol = solve(prob,Tsit5(),callback = cbs,saveat = Ts/2)
#out_x = map(x -> x[1:2],sol.u) #posição
#out_dx = map(x -> x[3:end],sol.u)
