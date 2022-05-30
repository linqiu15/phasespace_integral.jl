module PhaseSpaceIntegral

export Proc, refresh_proc
export ps12, ps13_dalitz, ps13_kumar, ps14_kumar, ps15_kumar, ps23_byckling, ps24_byckling


import QuadGK: quadgk, gauss
import Cubature: hcubature as ccubature
import HCubature: hcubature
import Cuba: vegas, suave, divonne, cuhre
import StaticArrays: @SVector
import LinearAlgebra: det

include("others/quadgauss.jl")

λ(x, y, z) = x^2 + y^2 + z^2 - 2 * x * y - 2 * x * z - 2 * y * z
diff(x) = abs(x[1] - x[2])
partition(x, v) = v[1] * (1 - x) + v[2] * x

global xxx, www = gauss(10, 0, 1)
function setquadgauss(N::Int)
    global xxx, www = gauss(N, 0, 1)
    println("Gauss Quadrature Points:$N !")
    nothing
end

mutable struct Proc
    init_m::Vector{Float64}
    final_m::Vector{Float64}

    function Proc(init_v, final_v)
        if length(init_v) == 1 && sum(init_v) < sum(final_v)
            error("Error: the initial threshold should be above the final one!")
        end
        new(init_v, final_v)
    end
end

function refresh_proc(proc::Proc, init_v, final_v)
    if length(init_v) == 1 && sum(init_v) < sum(final_v)
        error("Error: the initial threshold should be above the final one!")
    end
    proc.init_m = init_v
    proc.final_m = final_v
    nothing
end

integral_methods = ["quadgk", "quadgauss", "ccubature", "hcubature", "vegas", "suave", "divonne", "cuhre"]

include("quadratic.jl")
include("1n.jl")
include("2n.jl")


end