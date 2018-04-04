#Arquivo que contém o código para a busca dos melhores ganhos utilizando evolução diferencial
include("search-model.jl")
using Evolutionary
# Parâmetros da simulação
Ts = 0.05 # Intervalo entre leituras da saída
tend = 2.0 # tempo final para estabilização
t0 = 0.0 # instante inicial
r1 =  0.6#1.6# referência junta 1
r2 =  0.8#2.0# referência junta 2


function custo(x::Vector{Float64})
    kp = SMatrix{2,2}(diagm([x[1], x[2]]))
    kv = SMatrix{2,2}(diagm([x[3], x[4]]))
    x, v, t, a, ta, j, tj = simulation(kp, kv, Ts, t0, tend, r1, r2)
    erro1 = -(x[1] - r1)
    erro2 = -(x[2] - r2)
    size = length(erro1)
    # essa função custo é um somátorio do módulo das normas dos erros mais um somatório da norma dos trancos
    #out = sum(abs.(erro1))/maximum(erro1) + sum(abs.(erro2))/maximum(erro2) + sum(abs.(j[1]))/maximum(j[1]) + sum(abs.(j[2]))/maximum(j[2])

    # essa função custo leva em conta o erro inicial e final e o jerk, procurando assim minimizar a diferença entre a posição desejado e o erro inicial e minimizar o erro final bem como o jerk
    erro_init_1 = sum(abs.(x[1][1:floor(Integer,size/3)]))/(maximum(abs.(x[1][1:floor(Integer,size/3)])))
    erro_init_2 = sum(abs.(x[2][1:floor(Integer,size/3)]))/(maximum(abs.(x[2][1:floor(Integer,size/3)])))


    erro_end_1 = sum(abs.(erro1[floor(Integer,size/3):end]))/maximum(abs.(erro1))
    erro_end_2 = sum(abs.(erro2[floor(Integer,size/3):end]))/maximum(abs.(erro2))
    # erro_end = sum(abs.(erro1[floor(Integer,size/3):end]))/r1 + sum(abs.(erro2[floor(Integer,size/3):end]))/r2

    jerk_1 = sum(abs.(j[1]))/maximum(abs.(j[1]))
    jerk_2 = sum(abs.(j[2]))/maximum(abs.(j[2]))

    erro_init = erro_init_1 + erro_init_2
    erro_end = erro_end_1 + erro_end_2
    jerk = jerk_1 + jerk_2
    println("$(erro_init) + $(erro_end) + $(jerk)")
    out =  erro_init + erro_end + jerk

    # essa função é a mesma da anterio porém considerando metade metado como valores de erro inicial e final
    #out = sum(abs.(erro1[1:floor(Integer,size/2)] - r1))/(abs(maximum(erro1)- r1)) + sum(abs.(erro2[1:floor(Integer,size/2)] - r2))/(abs(maximum(erro2)- r2)) +  sum(abs.(erro1[2*floor(Integer,size/2):end]))/maximum(erro1) + sum(abs.(erro2[2*floor(Integer,size/2):end]))/maximum(erro2) + sum(abs.(j[1]))/maximum(j[1]) + sum(abs.(j[2]))/maximum(j[2])

    out
end

#### Para a evolução diferencial
# N  : número de parâmetros do vetor de entrada
# λ  : descendência
# μ  : número de pais
# ρ  : número de pais envolvidos na procriação de uma prole

#N = 4
#P = 100
#result, fitness, cnt = cmaes(custo, N; μ = 15,  λ = P, iterations = 1)

#### Para algoritmo genético
# N : número de parâmetros do vetor de entrada
# ɛ : Quantidade de indivíduos da geração atual que terão sobrevivência garantida na próxima geração

N = 4
println("iniciando")
result, fitness, cnt = ga(custo, N; initPopulation = (n -> rand(n).*[1000., 100., 1000., 100.]), populationSize = 50, ɛ = 0.1, selection = sus, crossover = intermediate(0.25), mutation = domainrange(fill(0.5,N)), iterations = 50)

t_end_new = tend + 3.
kp = SMatrix{2,2}(diagm(result[1:2]))
kv = SMatrix{2,2}(diagm(result[3:4]))

x, v, t, a, ta, j, tj = simulation(kp, kv, Ts, t0, t_end_new, r1, r2);

p1 = plot(t,x[1], label = "desejado 1")
p1= plot!([r1],seriestype= :hline, label = "referência");
p2 = plot(t,x[2], label = "desejado 2")
p2 = plot!([r2],seriestype= :hline, label = "referência");
println("------------------------------------------------------------")
println("O resultado foi $(result) com $(cnt) iterações e $(fitness)")
println("posição 1 final: $(x[1][end]), posição 2 final $(x[2][end])")
println("erro 1 final: $(rad2deg(x[1][end] - r1)) graus, posição 2 final $(rad2deg(x[2][end] - r2)) graus")
println("erro 1% = $(((r1 - x[1][end])/r1)*100)%, erro 2% = $(((r2 - x[2][end])/r2)*100)%")
plot(p1,p2)

# BUG: Encontre um problema, quando ocorre um erro na integração eu mando o vetor incompleto para o lado de fora e ele usa esses valores no cálculo do fitness. Devo corrigir isso!!!
