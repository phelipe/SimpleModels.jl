#Arquivo que contém o código para a busca dos melhores ganhos utilizando evolução diferencial
include("search-model.jl")

# Parâmetros da simulação
Ts = 0.08 # Intervalo entre leituras da saída
tend = 2.0 # tempo de simulação
t0 = 0.0 # instante inicial
r1 = 1.2 # referência junta 1
r2 = 0.6 # referência junta 2




# NOTE: este código ficará dentro da parte do otimizador para obter os valores de posição e jerk que serão utilizados na função objetivo.
# NOTE: Os ganhos devem ser todos maiores ou iguais a zero
kp = SMatrix{2,2}(diagm([1000., 30.]))
kv = SMatrix{2,2}(diagm([250., 3.]))
x, v, t, a, ta, j, tj = simulation(kp, kv, Ts, t0, tend, r1, r2)
erro1 = -(x[1] - r1)
erro2 = -(x[2] - r2)

############
