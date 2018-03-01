struct Dof2 <:Robot
    m::Vector{AbstractFloat}
    l::Vector{AbstractFloat}
    r::Vector{AbstractFloat}
    I::Vector{AbstractFloat}
    g::AbstractFloat
    dof::Int64  #::AbstractFloat

    function Dof2(m::Vector{T}, l::Vector{T}, r::Vector{T}, I::Vector{T}, g = 9.81) where T<:AbstractFloat
        if (length(m) != 2 || length(l) != 2 || length(r) != 2 || length(I) != 2)
            error("The length must be 2")
        else
            return new(m,l,r,I,g,2)
        end
    end
end

function dinamic(robot::Dof2)
    out = function(du, u, p, t)
        θ = u[1:2]
        dθ = u[3:4]
        cos1 = cos(θ[1])
        cos2 = cos(θ[2])
        sin1 = sin(θ[1])
        sin2 = sin(θ[2])
        cos12 = cos(θ[1]+θ[2])

        #Matriz de inércia
        m11 = robot.m[1] * robot.r[1]^2 + robot.m[2]* (robot.l[1]^2 + robot.r[2]^2 + 2*robot.l[1] * robot.r[2] * cos2)+ robot.I[1] + robot.I[2]
        m12 = robot.m[2] * (robot.r[2]^2 + robot.l[1] *robot.r[2] * cos2) + robot.I[2]
        m22 = robot.m[2] * robot.r[2]^2 + robot.I[2]
        M = [m11 m12; m12 m22]

        #Matriz C
        h = -robot.m[2] * robot.l[1] * robot.r[2] * sin2
        c11 = h * dθ[2]
        c12 = h * (dθ[1] + dθ[2])
        c21 = -h * dθ[1]
        c22 = 0
        C = [c11 c12; c21 c22]

        #Matriz gravidade
        g2 = robot.m[2] * robot.r[2] * robot.g * cos12
        g1 = (robot.m[1] * robot.r[1] + robot.m[2] * robot.l[1]) * robot.g * cos1 + g2
        G = [g1; g2]

        du[1:2] = dθ
        du[3:4] = inv(M)*(u.tau - C * dθ - G)
    end

end

function organize(model::Dof2, data)
    out_x = map(x -> x[1:model.dof],data.u)
    out_dx = map(x -> x[(model.dof+1):end],data.u)
    #velocidade
    #Aqui estou fazendo uma aproximação da aceleração e do jerk
    #out_d2x = diff(out_dx)/Ts #aceleração
    #out_d3x = diff(out_d2x)/Ts #jerk
    θ = []
    for i=1:model.dof
        push!(θ, map(x -> x[i],out_x))
    end
    ω = []
    for i=1:model.dof
        push!(ω, map(x -> x[i],out_dx))
    end
    α = []
    for i=1:model.dof
        push!(α, diff(ω[i])./diff(data.t)[1])
    end
    ta = collect(0:(data.t[end]/(length(α[1])-1)):data.t[end])
    J = []
    for i=1:model.dof
        push!(J, diff(α[i])./diff(ta)[1])
    end
    tj = collect(0:(data.t[end]/(length(J[1])-1)):data.t[end])

    θ, ω, data.t, α, ta, J, tj
end
