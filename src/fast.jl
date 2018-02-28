workspace()
push!(LOAD_PATH, expanduser("/home/phelipe/Documentos"))
using SimpleModels
model = Dof2([1., 1.], [2., 1.], [1.0, 0.5], [0.3, 0.5],0.)
