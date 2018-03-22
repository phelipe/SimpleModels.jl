#Arquivo que contém o código para a busca dos melhores ganhos utilizando evolução diferencial
include("search-model.jl")

# Parâmetros da simulação
Ts = 0.08 # Intervalo entre leituras da saída
tend = 2.0 # tempo de simulação
t0 = 0.0 # instante inicial
r1 = 1.2 # referência junta 1
r2 = 0.6 # referência junta 2


function custo(x::Vector{Float64})

    kp = SMatrix{2,2}(diagm([x[1], x[2]]))
    kv = SMatrix{2,2}(diagm([x[3], x[4]]))
    x, v, t, a, ta, j, tj = simulation(kp, kv, Ts, t0, tend, r1, r2)
    erro1 = -(x[1] - r1)
    erro2 = -(x[2] - r2)
    out = sum(abs.(erro1)) + sum(abs.(erro2)) + sum(abs.(j[1])) + sum(abs.(j[2]))
    out
end


# N  = número de parâmetros do vetor de entrada
# λ  = descendência
# μ  = número de pais
# ρ  = número de pais envolvidos na procriação de uma prole

result, fitness, cnt = es(custo, N;
            initPopulation = generatePositions,
            mutation = mutationwrapper(muts),
μ = 15, ρ = 1, λ = P)
