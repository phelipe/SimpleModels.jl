# Arquivo contendo comparação cinemática, neste arquivo farei a geração de trajetória no domínio cartesiano e no domínio das juntas para comparação.
push!(LOAD_PATH, expanduser("/home/phelipe/Documentos"))
using Trajectory
using Plots

m=[23.902, 1.285]
l= [.45, .45]
r= [0.091, 0.048]
I= [1.266, 0.093]
g = 9.81
t0 = 0.
tend = 3.
Ts = 0.01

function cinematica_direta(θ::T) where T<:AbstractArray
    x = l[1] * sin(θ[1]) + l[2] * sin(θ[1] + θ[2])
    y = -(l[1] * cos(θ[1]) + l[2] * cos(θ[1] + θ[2]))
    [x,y]
end

function cinematica_inversa(cord::Vector{T}) where T<:AbstractFloat
    x = cord[1]
    y = cord[2]
    t = (x^2 + y^2 - l[1]^2 - l[2]^2) / (2 * l[1] * l[2])
    θ2 = atan2( (1 - t^2), t)
    k1 = l[1] + l[2] * cos(θ2)
    k2 = l[2] * sin(θ2)
    γ = atan2(k2,k1)
    θ1 = atan2(x,-y) - γ
    [θ1, θ2]
end


position_start = cinematica_direta([0.,0.])
position_end = cinematica_direta([deg2rad(90),deg2rad(90)])

#teta1 = deg2rad(45)
#teta2 = deg2rad(0)
#println("Valores na cinemática direta $(cord = cinematica_direta([teta1,teta2])) metros")
#println("Valores na cinemática inversa $(rad2deg.(cinematica_inversa(cord))) deg")

#mínimo jerk cartesiano
x1,v1,a1,j1 = minimumjerkf(position_start[1],0.,0.,t0,position_end[1],0.,0.,tend)
x2,v2,a2,j2 = minimumjerkf(position_start[2],0.,0.,t0,position_end[2],0.,0.,tend)

time = collect(t0:Ts:tend)
x1data = map(x-> x1(x),time)
x2data = map(x-> x1(x),time)
