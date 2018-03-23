#Arquivo que contém o código para a busca dos melhores ganhos utilizando evolução diferencial
include("search-model.jl")
using Evolutionary
# Parâmetros da simulação
Ts = 0.08 # Intervalo entre leituras da saída
tend = 4.0 # tempo de simulação
t0 = 0.0 # instante inicial
r1 =  1.8# referência junta 1
r2 =  2.0# referência junta 2


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
    #
    out = sum(abs.(erro1[1:floor(Integer,size/3)] - r1))/(abs(maximum(erro1)- r1)) + sum(abs.(erro2[1:floor(Integer,size/3)] - r2))/(abs(maximum(erro2)- r2)) +  sum(abs.(erro1[2*floor(Integer,size/3):end]))/maximum(erro1) + sum(abs.(erro2[2*floor(Integer,size/3):end]))/maximum(erro2) + sum(abs.(j[1]))/maximum(j[1]) + sum(abs.(j[2]))/maximum(j[2])
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
tic()
result, fitness, cnt = ga(custo, N; initPopulation = (n -> rand(n).*[1000., 100., 100., 10.]), populationSize = 100, ɛ = 0.1, selection = sus, crossover = intermediate(0.25), mutation = domainrange(fill(0.5,N)), iterations = 1000)
toc()

kp = SMatrix{2,2}(diagm(result[1:2]))
kv = SMatrix{2,2}(diagm(result[3:4]))
println("O resultado foi $(result) com $(cnt) iterações e $(fitness)")
x, v, t, a, ta, j, tj = simulation(kp, kv, Ts, t0, tend, r1, r2);
plot(x[1])


# NOTE: O jerk tem um valor muito maior que o erro de posição então na hora de minimizar, reduzir o valor do jerk se torna muito mais importante que reduzir o erro. Para resolver isso fiz uma normalização do jerk e do erro e já melhorou o resultado, poréma apresentava sobressinal e isso foi resolvido com uma melhoria na equação base.
# NOTE: Agora devo cuidar da comparação com o PID tradicional para verificar isso.
# NOTE: Comparar com o pid do minimo jerk
