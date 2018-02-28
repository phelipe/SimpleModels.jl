struct Dof2 <:Robot
    m::Vector{AbstractFloat}
    l::Vector{AbstractFloat}
    r::Vector{AbstractFloat}
    I::Vector{AbstractFloat}
    g::AbstractFloat
    const dof = 2 #::AbstractFloat

    function Dof2(m::Vector{T}, l::Vector{T}, r::Vector{T}, I::Vector{T}, g = 9.81) where T<:AbstractFloat
        if (length(m) != dof || length(l) != dof || length(r) != dof || length(I) != dof)
            error("The length must be 2")
        else
            return new(m,l,r,I,g)
        end
    end

    function dinamic(du, u, p, t)
        θ = u[1:2]
        dθ = u[3:4]
        cos1 = cos(θ[1])
        cos2 = cos(θ[2])
        sin1 = sin(θ[1])
        sin2 = sin(θ[2])
        cos12 = cos(θ[1]+θ[2])

        #Matriz de inércia
        m11 = m[1] * r[1]^2 + m[2]* (l[1]^2 + r[2]^2 + 2*l[1] * r[2] * cos2)+ I[1] + I[2]
        m12 = m[2] * (r[2]^2 + l[1] *r[2] * cos2) +I[2]
        m22 = m[2] * r[2]^2 + I[2]
        M = [m11 m12; m12 m22]

        #Matriz C
        h = -m[2] * l[1] * r[2] * sin2
        c11 = h * dθ[2]
        c12 = h * (dθ[1] + dθ[2])
        c21 = -h * dθ[1]

        c22 = 0
        C = [c11 c12; c21 c22]

        #Matriz gravidade
        g2 = m[2] * r[2] * g * cos12
        g1 = (m[1] * r[1] + m[2] * l[1]) * g * cos1 + g[2]
        G = [g1; g2]

        du[1:2] = dθ
        du[3:4] = inv(M)*(u.tau - C * dθ - G)
    end

end
