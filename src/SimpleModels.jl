module SimpleModels

#=using DifferentialEquations
type InputRobot{T} <: DEDataVector{T}
    x::Array{T,1}
    tau::Array{T,1}
end
=#
export Dof2

abstract type Robot end

include("1dofrobot.jl")
include("2dofrobot.jl")


end
