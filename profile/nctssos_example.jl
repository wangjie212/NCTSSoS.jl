using BenchmarkTools
using MosekTools, NCTSSoS
using Profile 

n = 10 
@ncpolyvar x[1:n]
f = 0.0
for i = 1:n
    jset = max(1, i - 5):min(n, i + 1)
    jset = setdiff(jset, i)
    f += (2x[i] + 5 * x[i]^3 + 1)^2
    f -= sum([4x[i] * x[j] + 10x[i]^3 * x[j] + 2x[j] + 4x[i] * x[j]^2 + 10x[i]^3 * x[j]^2 + 2x[j]^2 for j in jset])
    f += sum([x[j] * x[k] + 2x[j]^2 * x[k] + x[j]^2 * x[k]^2 for j in jset for k in jset])
end

cons = typeof(f)[]
for i = 1:n
    push!(cons, 1 - x[i]^2)
    push!(cons, x[i] - 1 / 3)
end
pop =  PolyOpt(f; constraints = cons);
solver_config = SolverConfig(optimizer=Mosek.Optimizer; order=3, cs_algo=MF(), ts_algo=AsIsElimination())

result = cs_nctssos(pop, solver_config)
# Problem
#   Name                   :                 
#   Objective sense        : maximize        
#   Type                   : CONIC (conic optimization problem)
#   Constraints            : 56743           
#   Affine conic cons.     : 0               
#   Disjunctive cons.      : 0               
#   Cones                  : 0               
#   Scalar variables       : 1               
#   Matrix variables       : 1068 (scalarized: 137680)
#   Integer variables      : 0      

result.objective

Profile.clear()
@profile cs_nctssos(pop, solver_config)
Profile.print(mincount=100,format=:flat)
