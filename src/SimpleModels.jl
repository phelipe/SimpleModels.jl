module SimpleModels

import StaticArrays
import StaticArrays: SVector, SMatrix

export Dof2
export dinamic
export organize
abstract type Robot end

include("1dofrobot.jl")
include("2dofrobot.jl")

end
