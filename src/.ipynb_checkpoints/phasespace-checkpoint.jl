module PhaseSpaceIntegral

export Proc, refresh_proc
export ps12, ps13_dalitz, ps13_kumar, ps14_kumar


using QuadGK: quadgk, gauss
using HCubature: hcubature
using Cuba: divonne, cuhre
using StaticArrays: @SVector

include("others/quadgauss.jl")

λ(x, y, z) = x^2 + y^2 + z^2 - 2 * x * y - 2 * x * z - 2 * y * z
diff(x) = abs(x[1] - x[2])
partition(x, v) = v[1] * (1 - x) + v[2] * x

global xxx, www = gauss(10, 0, 1)
function setquadgauss(N::Int)
    global xxx, www = gauss(N, 0, 1)
    nothing
end

mutable struct Proc
    init_m::Vector{Float64}
    final_m::Vector{Float64}

    function Proc(init_v, final_v)
        if sum(init_v) < sum(final_v)
            error("Error: the initial threshold should be above the final one!")
        end
        new(init_v, final_v)
    end
end

function refresh_proc(proc::Proc, init_v, final_v)
    if sum(init_v) < sum(final_v)
        error("Error: the initial threshold should be above the final one!")
    end
    proc.init_m = init_v
    proc.final_m = final_v
    nothing
end

integral_methods = ["quadgk", "quadgauss", "hcubature", "divonne", "cuhre"]


#----------------------------------1n--------------------------------------------

##---------------------------------12--------------------------------------------
function ps12(ampsq::Function, proc::Proc)
    M = proc.init_m[1]
    m1, m2 = proc.final_m
    return ampsq() * sqrt(λ(M^2, m1^2, m2^2)) / (8 * π * M^2)
end

##---------------------------------13(dalitz)--------------------------------------------

function s23f(init_v::T, final_v::T, m12::Float64) where {T<:Vector{Float64}}
    M = init_v[1]
    m1, m2, m3 = final_v
    E2star, E3star = (m12^2 - m1^2 + m2^2) / (2 * m12), (M^2 - m12^2 - m3^2) / (2 * m12)
    return @SVector [(E2star + E3star)^2 - abs2(sqrt(abs(E2star^2 - m2^2)) + sqrt(abs(E3star^2 - m3^2))), (E2star + E3star)^2 - abs2(sqrt(abs(E2star^2 - m2^2)) - sqrt(abs(E3star^2 - m3^2)))]
end

function ps13_dalitz(ampsq::Function, proc::Proc; method="hcubature", xxx=xxx, www=www)
    M = proc.init_m[1]
    m1, m2, m3 = proc.final_m
    res = 0.0

    s12(x) = (m1 + m2)^2 * (1 - x) + (M - m3)^2 * x
    function s23(x, y)
        s23_min, s23_max = s23f(proc.init_m, proc.final_m, sqrt(s12(x)))
        return s23_min * (1 - y) + s23_max * y
    end


    if method == "quadgk" || !(method in integral_methods)
        error("Error: no such an integration method!")
    elseif method == "quadgauss"
        res = quadgauss(x -> quadgauss(y -> ampsq(s12(x), s23(x, y)) * diff(s23f(proc.init_m, proc.final_m, sqrt(s12(x)))), xxx, www), xxx, www)
    elseif method == "hcubature"
        res = hcubature(x -> ampsq(s12(x[1]), s23(x...)) * diff(s23f(proc.init_m, proc.final_m, sqrt(s12(x[1])))), [0, 0], [1, 1])[1]
    elseif method == "divonne"
        res = divonne((x, f) -> f[1] = ampsq(s12(x[1]), s23(x...)) * diff(s23f(proc.init_m, proc.final_m, sqrt(s12(x[1])))))[1]
    elseif method == "cuhre"
        res = cuhre((x, f) -> f[1] = ampsq(s12(x[1]), s23(x...)) * diff(s23f(proc.init_m, proc.final_m, sqrt(s12(x[1])))))[1]
    end

    return 1 / ((2π)^3 * 16 * M^2) * abs((M - m3)^2 - (m1 + m2)^2) * res
end

## --------------------------13(kumar)----------------------------------------------------------
function kumar_u1_13(init_v::T, final_v::T, s1::Float64) where {T<:Vector{Float64}}
    M = init_v[1]
    m1, m2, m3 = final_v
    u1_min = M^2 + m2^2 - (s1 + m2^2 - m3^2) * (M^2 + s1 - m1^2) / (2 * s1) - sqrt(abs(λ(s1, m2^2, m3^2) * λ(M^2, s1, m1^2))) / (2 * s1)
    u1_max = M^2 + m2^2 - (s1 + m2^2 - m3^2) * (M^2 + s1 - m1^2) / (2 * s1) + sqrt(abs(λ(s1, m2^2, m3^2) * λ(M^2, s1, m1^2))) / (2 * s1)
    return @SVector [u1_min, u1_max]
end

function ps13_kumar(ampsq::Function, proc::Proc; method="hcubature", xxx=xxx, www=www)
    M = proc.init_m[1]
    m1, m2, m3 = proc.final_m
    res = 0.0

    s1(x::Float64) = (m2 + m3)^2 * (1 - x) + (M - m1)^2 * x
    function u1(x::Float64, y::Float64)
        u1_min, u1_max = kumar_u1_13(proc.init_m, proc.final_m, s1(x))
        return u1_min * (1 - y) + u1_max * y
    end

    if method == "quadgk" || !(method in integral_methods)
        error("Error: no such an integration method!")
    elseif method == "quadgauss"
        res = quadgauss(x -> quadgauss(y -> ampsq(s1(x), u1(x, y)) * diff(kumar_u1_13(proc.init_m, proc.final_m, s1(x))), xxx, www), xxx, www)
    elseif method == "hcubature"
        res = hcubature(x -> ampsq(s1(x[1]), u1(x[1], x[2])) * diff(kumar_u1_13(proc.init_m, proc.final_m, s1(x[1]))), [0, 0], [1, 1])[1]
    elseif method == "divonne"
        res = divonne((x, f) -> f[1] = ampsq(s1(x[1]), u1(x[1], x[2])) * diff(kumar_u1_13(proc.init_m, proc.final_m, s1(x[1]))))[1]
    elseif method == "cuhre"
        res = cuhre((x, f) -> f[1] = ampsq(s1(x[1]), u1(x[1], x[2])) * diff(kumar_u1_13(proc.init_m, proc.final_m, s1(x[1]))))[1]
    end

    return 1 / (2^7 * π^3 * M^2) * ((M - m1)^2 - (m2 + m3)^2) * res
end



#------------------------------14(kumar)----------------------------------------------------
# invariants: s1,u1,s2,u2,t2
kumar_s2p(M, m1, m2, s1, u1, s2) = s2 + M^2 + (m1^2 + m2^2) - (s1 + u1)
kumar_ξ2(M, m1, m2, u1, s2, s2p) = ((M^2 + s2p - s2) * (M^2 + m2^2 - u1) - 2 * M^2 * (s2p + m2^2 - m1^2)) / sqrt(abs(λ(M^2, s2, s2p) * λ(M^2, u1, m2^2)))
kumar_η2(M, m3, m4, s2, u2, s2p) = (2 * M^2 * (s2 + m3^2 - m4^2) - (M^2 + m3^2 - u2) * (M^2 + s2 - s2p)) / sqrt(abs(λ(M^2, s2, s2p) * λ(M^2, m3^2, u2)))
kumar_ω2(M, m2, m3, u1, u2, t2) = (2 * M^2 * (u1 + m3^2 - t2) - (M^2 + m3^2 - u2) * (M^2 + u1 - m2^2)) / sqrt(abs(λ(M^2, u1, m2^2) * λ(M^2, m3^2, u2)))
kumar_ζ2(ξ2, η2, ω2) = (ω2 - ξ2 * η2) / sqrt(abs((1 - ξ2^2) * (1 - η2^2)))

function kumar_u1_14(init_v::T, final_v::T, s1::U, s2::U) where {T<:Vector{Float64},U<:Float64}
    M = init_v[1]
    m1, m2, m3, m4 = final_v
    u1_min = M^2 + m2^2 - (s1 + m2^2 - s2) * (M^2 + s1 - m1^2) / (2 * s1) - sqrt(abs(λ(s1, m2^2, s2) * λ(M^2, s1, m1^2))) / (2 * s1)
    u1_max = M^2 + m2^2 - (s1 + m2^2 - s2) * (M^2 + s1 - m1^2) / (2 * s1) + sqrt(abs(λ(s1, m2^2, s2) * λ(M^2, s1, m1^2))) / (2 * s1)

    return @SVector [u1_min, u1_max]
end

function kumar_u2_14(init_v::T, final_v::T, s2::U, s2p::U) where {T<:Vector{Float64},U<:Float64}
    M = init_v[1]
    m1, m2, m3, m4 = final_v
    u2_min = M^2 + m3^2 - (s2 + m3^2 - m4^2) * (M^2 + s2 - s2p) / (2 * s2) - sqrt(abs(λ(s2, m3^2, m4^2) * λ(M^2, s2, s2p))) / (2 * s2)
    u2_max = M^2 + m3^2 - (s2 + m3^2 - m4^2) * (M^2 + s2 - s2p) / (2 * s2) + sqrt(abs(λ(s2, m3^2, m4^2) * λ(M^2, s2, s2p))) / (2 * s2)

    return @SVector [u2_min, u2_max]
end

function kumar_t2_14(init_v::T, final_v::T, u1::U, u2::U, ξ2::U, η2::U) where {T<:Vector{Float64},U<:Float64}
    M = init_v[1]
    m1, m2, m3, m4 = final_v
    t2_min = u1 + m3^2 - (M^2 + m3^2 - u2) * (M^2 + u1 - m2^2) / (2 * M^2) + sqrt(abs(λ(M^2, m3^2, u2) * λ(M^2, u1, m2^2))) / (2 * M^2) * (-ξ2 * η2 - sqrt(abs((1 - ξ2^2) * (1 - η2^2))))
    t2_max = u1 + m3^2 - (M^2 + m3^2 - u2) * (M^2 + u1 - m2^2) / (2 * M^2) + sqrt(abs(λ(M^2, m3^2, u2) * λ(M^2, u1, m2^2))) / (2 * M^2) * (-ξ2 * η2 + sqrt(abs((1 - ξ2^2) * (1 - η2^2))))
    return @SVector [t2_min, t2_max]
end


function kumar_integrand_14(ampsq::Function, init_v::T, final_v::T, x) where {T<:Vector{Float64}}
    M = init_v[1]
    m1, m2, m3, m4 = final_v
    x1, x2, x3, x4, x5 = x

    s1 = (m2 + m3 + m4)^2 * (1 - x1) + (M - m1)^2 * x1
    s2 = (m3 + m4)^2 * (1 - x2) + (sqrt(s1) - m2)^2 * x2
    u1 = partition(x3, kumar_u1_14(init_v, final_v, s1, s2))
    s2p = kumar_s2p(M, m1, m2, s1, u1, s2)
    u2 = partition(x4, kumar_u2_14(init_v, final_v, s2, s2p))
    ξ2 = kumar_ξ2(M, m1, m2, u1, s2, s2p)
    η2 = kumar_η2(M, m3, m4, s2, u2, s2p)
    t2 = partition(x5, kumar_t2_14(init_v, final_v, u1, u2, ξ2, η2))
    ω2 = kumar_ω2(M, m2, m3, u1, u2, t2)
    ζ2 = kumar_ζ2(ξ2, η2, ω2)

    res = (((M - m1)^2 - (m2 + m3 + m4)^2) * ((sqrt(s1) - m2)^2 - (m3 + m4)^2) * diff(kumar_u1_14(init_v, final_v, s1, s2)) * diff(kumar_u2_14(init_v, final_v, s2, s2p)) * diff(kumar_t2_14(init_v, final_v, u1, u2, ξ2, η2)) * ampsq(s1, s2, u1, u2, t2) / sqrt(abs(λ(M^2, s2, s2p) * λ(M^2, m2^2, u1) * (1 - ξ2^2) * λ(M^2, m3^2, u2) * (1 - η2^2) * (1 - ζ2^2))))

    if isnan(res) || isinf(res)
        error("Error: the integrand evaluates to NaN or Inf!")
    end

    return res
end

function ps14_kumar(ampsq::Function, proc::Proc)
    res = cuhre((x, f) -> f[1] = kumar_integrand_14(ampsq, proc.init_m, proc.final_m, x), 5)[1]
    return 1 / (2^10 * π^6) * res
end


end