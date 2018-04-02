# Arquivo contendo comparação cinemática, neste arquivo farei a geração de trajetória no domínio cartesiano e no domínio das juntas para comparação.
push!(LOAD_PATH, expanduser("/home/phelipe/Documentos"))
using Trajectory
using Plots
pyplot()

m=[23.902, 1.285]
l= [.45, .45]
r= [0.091, 0.048]
I= [1.266, 0.093]
g = 9.81
t0 = 0.
tend = 4.
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

θ_start = [0.,0.]
#θ_end = [deg2rad(10),deg2rad(90)]
θ_end = [1.8, 2.0]
position_start = cinematica_direta(θ_start)
position_end = cinematica_direta(θ_end)

#teta1 = deg2rad(45)
#teta2 = deg2rad(0)
#println("Valores na cinemática direta $(cord = cinematica_direta([teta1,teta2])) metros")
#println("Valores na cinemática inversa $(rad2deg.(cinematica_inversa(cord))) deg")

#mínimo jerk cartesiano
x1,v1,a1,j1 = minimumjerkf(position_start[1],0.,0.,t0,position_end[1],0.,0.,tend)
x2,v2,a2,j2 = minimumjerkf(position_start[2],0.,0.,t0,position_end[2],0.,0.,tend)

time = collect(t0:Ts:tend)
x1data = map(a-> x1(a),time)
x2data = map(a-> x2(a),time)
function plotcartesiano()
    p1 = plot(time,x1data,xlabel="tempo",ylabel="posição x")
    p2 = plot(time,x2data,xlabel="tempo",ylabel="posição y")
    plot(p1,p2, tittle="Cartesiano")
end

θ = map((x,y)-> cinematica_inversa([x,y]),x1data,x2data)
θ1 = map(x->x[1],θ)
θ2 = map(x->x[2],θ)
vθ1 = diff(θ1)./Ts
vθ2 = diff(θ2)./Ts
aθ1 = diff(vθ1)./Ts
aθ2 = diff(vθ2)./Ts
jθ1 = diff(aθ1)./Ts
jθ2 = diff(aθ2)./Ts

function plotjuntas()
    p1 = plot(time,θ1,xlabel="tempo",ylabel="posição θ1")
    p2 = plot(time,θ2,xlabel="tempo",ylabel="posição θ2")
    plot(p1,p2, tittle="Jerk cartesiano nas juntas")
end

function plotvjuntas()
    p1 = plot(time[1:length(vθ1)],vθ1,xlabel="tempo",ylabel="velocidade θ1")
    p2 = plot(time[1:length(vθ1)],vθ2,xlabel="tempo",ylabel="velocidade θ2")
    plot(p1,p2, tittle="Jerk cartesiano nas juntas")
end

function plotajuntas()
    p1 = plot(time[1:length(aθ1)],aθ1,xlabel="tempo",ylabel="aceleração θ1")
    p2 = plot(time[1:length(aθ1)],aθ2,xlabel="tempo",ylabel="aceleração θ2")
    plot(p1,p2, tittle="Jerk cartesiano nas juntas")
end

function plotjjuntas()
    p1 = plot(time[1:length(jθ1)],jθ1,xlabel="tempo",ylabel="tranco θ1")
    p2 = plot(time[1:length(jθ2)],jθ2,xlabel="tempo",ylabel="tranco θ2")
    plot(p1,p2, tittle="Jerk cartesiano nas juntas")
end


# domínio das juntas
xj1,vj1,aj1,jj1 = minimumjerkf(θ_start[1],0.,0.,t0,θ_end[1],0.,0.,tend)
xj2,vj2,aj2,jj2 = minimumjerkf(θ_start[2],0.,0.,t0,θ_end[2],0.,0.,tend)

θ_junta1 = map(a-> xj1(a),time)
θ_junta2 = map(a-> xj2(a),time)
v_junta1 = map(a-> vj1(a),time)
v_junta2 = map(a-> vj2(a),time)
a_junta1 = map(a-> aj1(a),time)
a_junta2 = map(a-> aj2(a),time)
j_junta1 = map(a-> jj1(a),time)
j_junta2 = map(a-> jj2(a),time)

cart = map((x,y)-> cinematica_direta([x,y]),θ_junta1,θ_junta2)
cart_x = map(x->x[1],cart)
cart_y = map(x->x[2],cart)

function plotjerks()
    p1 = plot(time,[jθ1,j_junta1],xlabel="tempo",ylabel="tranco θ1")
    p2 = plot(time,[jθ2,j_junta2],xlabel="tempo",ylabel="tranco θ2")
    plot(p1,p2)
end

function plotpositions()
    p1 = plot(time,[θ1,θ_junta1],xlabel="tempo",ylabel="posição θ1")
    p2 = plot(time,[θ2,θ_junta2],xlabel="tempo",ylabel="posição θ2")
    plot(p1,p2)
end

function plotvelocitys()
    p1 = plot(time,[vθ1,v_junta1],xlabel="tempo",ylabel="velocidade θ1")
    p2 = plot(time,[vθ2,v_junta2],xlabel="tempo",ylabel="velocidade θ2")
    plot(p1,p2)
end

function plotacelerations()
    p1 = plot(time,[aθ1,a_junta1],xlabel="tempo",ylabel="aceleração θ1")
    p2 = plot(time,[aθ2,a_junta2],xlabel="tempo",ylabel="aceleração θ2")
    plot(p1,p2)
end

function plotxy()
    p1 = plot(x1data,x2data)
    p2 = plot(cart_x,cart_y)
    p = plot(p1,p2)
end

#Essas funções mostram o que acontece nas juntas quando o jerk é feito no domínio cartesiano
#plotjuntas()
#plotvjuntas()
#plotajuntas()
#plotjjuntas()

# Essas funções comparam as ações nas juntas
plotxy()
plotjerks()
