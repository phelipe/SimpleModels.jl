# Arquivo contendo simulação do robô 2D com PID
include("search-model.jl")
# Parâmetros da simulação
Ts = 0.05 # Intervalo entre leituras da saída
tend =  2.0 # tempo de simulação
t0 = 0.0 # instante inicial
r1 =  .8# referência junta 1
r2 =  1.6# referência junta 2

kp_pid = SMatrix{2,2}(diagm([2800., 80.]))
kv_pid = SMatrix{2,2}(diagm([315., 15.]))
x_pid, v_pid, t_pid, a_pid, ta_pid, j_pid, tj_pid = simulation(kp_pid, kv_pid, Ts, t0, tend, r1, r2)
erro1 = -(x_pid[1] - r1)
erro2 = -(x_pid[2] - r2)
plot(t_pid,x_pid)

p1 = plot(t_pid,x_pid[1], label = "pid 1")
p1= plot!([r1],seriestype= :hline, label = "referência");
p2 = plot(t_pid,x_pid[2], label = "pid 2")
p2 = plot!([r2],seriestype= :hline, label = "referência");


println("------------------------------------------------------------")
println("posição 1 final: $(x_pid[1][end]), posição 2 final $(x_pid[2][end])")
println("erro 1 final: $(rad2deg(x_pid[1][end] - r1)) graus, posição 2 final $(rad2deg(x_pid[2][end] - r2)) graus")
println("erro 1% = $(((r1 - x_pid[1][end])/r1)*100)%, erro 2% = $(((r2 - x_pid[2][end])/r2)*100)%")
println("Total jerk 1 = $(sum(abs.(j_pid[1]))), Total jerk 2 = $(sum(abs.(j_pid[2])))")
println("Max jerk 1 = $(maximum(abs.(j_pid[1]))), Max jerk 2 = $(maximum(abs.(j_pid[2])))")
plot(p1,p2)
