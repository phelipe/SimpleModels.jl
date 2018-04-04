# Arquivo contendo simulação do robô 2D com PID
include("search-model.jl")
# Parâmetros da simulação
Ts = 0.08 # Intervalo entre leituras da saída
tend =  4.0 # tempo de simulação
t0 = 0.0 # instante inicial
r1 =  .8# referência junta 1
r2 =  1.6# referência junta 2

kp_pid = SMatrix{2,2}(diagm([00., 400.]))
kv_pid = SMatrix{2,2}(diagm([00., 00.]))
x_pid, v_pid, t_pid, a_pid, ta_pid, j_pid, tj_pid = simulation(kp, kv, Ts, t0, tend, r1, r2)
erro1 = -(x_pid[1] - r1)
erro2 = -(x_pid[2] - r2)
plot(t_pid,x_pid)
