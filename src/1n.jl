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

    function integrand(x, y)
        s12 = (m1 + m2)^2 * (1 - x) + (M - m3)^2 * x
        s23_min, s23_max = s23f(proc.init_m, proc.final_m, sqrt(s12))
        s23 = s23_min * (1 - y) + s23_max * y
        return ampsq(s12, s23) * abs(s23_max - s23_min)
    end

    if method == "quadgk" || !(method in integral_methods)
        error("Error: no such an integration method!")
    elseif method == "quadgauss"
        res = quadgauss(x -> quadgauss(y -> integrand(x, y), xxx, www), xxx, www)
    elseif method == "ccubature"
        res = ccubature(x -> integrand(x...), [0, 0], [1, 1])[1]
    elseif method == "hcubature"
        res = hcubature(x -> integrand(x...), [0, 0], [1, 1])[1]
    elseif method == "vegas"
        res = vegas((x, f) -> f[1] = integrand(x...), 2)[1]
    elseif method == "suave"
        res = suave((x, f) -> f[1] = integrand(x...), 2)[1]
    elseif method == "divonne"
        res = divonne((x, f) -> f[1] = integrand(x...))[1]
    elseif method == "cuhre"
        res = cuhre((x, f) -> f[1] = integrand(x...))[1]
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

    function integrand(x, y)
        s1 = (m2 + m3)^2 * (1 - x) + (M - m1)^2 * x
        u1_min, u1_max = kumar_u1_13(proc.init_m, proc.final_m, s1)
        u1 = u1_min * (1 - y) + u1_max * y
        return ampsq(s1, u1) * abs(u1_max - u1_min)
    end

    if method == "quadgk" || !(method in integral_methods)
        error("Error: no such an integration method!")
    elseif method == "quadgauss"
        res = quadgauss(x -> quadgauss(y -> integrand(x, y), xxx, www), xxx, www)
    elseif method == "ccubature"
        res = ccubature(x -> integrand(x...), [0, 0], [1, 1])[1]
    elseif method == "hcubature"
        res = hcubature(x -> integrand(x...), [0, 0], [1, 1])[1]
    elseif method == "vegas"
        res = vegas((x, f) -> f[1] = integrand(x...), 2)[1]
    elseif method == "suave"
        res = suave((x, f) -> f[1] = integrand(x...), 2)[1]
    elseif method == "divonne"
        res = divonne((x, f) -> f[1] = integrand(x...))[1]
    elseif method == "cuhre"
        res = cuhre((x, f) -> f[1] = integrand(x...))[1]
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
    u1_span = kumar_u1_14(init_v, final_v, s1, s2)
    u1 = partition(x3, u1_span)
    s2p = kumar_s2p(M, m1, m2, s1, u1, s2)
    u2_span = kumar_u2_14(init_v, final_v, s2, s2p)
    u2 = partition(x4, u2_span)
    ξ2 = kumar_ξ2(M, m1, m2, u1, s2, s2p)
    η2 = kumar_η2(M, m3, m4, s2, u2, s2p)
    t2_span = kumar_t2_14(init_v, final_v, u1, u2, ξ2, η2)
    t2 = partition(x5, t2_span)
    ω2 = kumar_ω2(M, m2, m3, u1, u2, t2)
    ζ2 = kumar_ζ2(ξ2, η2, ω2)

    res = (((M - m1)^2 - (m2 + m3 + m4)^2) * ((sqrt(s1) - m2)^2 - (m3 + m4)^2) * diff(u1_span) * diff(u2_span) * diff(t2_span) * ampsq(s1, s2, u1, u2, t2) / sqrt(abs(λ(M^2, s2, s2p) * λ(M^2, m2^2, u1) * (1 - ξ2^2) * λ(M^2, m3^2, u2) * (1 - η2^2) * (1 - ζ2^2))))

    if isnan(res) || isinf(res)
        error("Error: the integrand evaluates to NaN or Inf!")
    end

    return res
end

function ps14_kumar(ampsq::Function, proc::Proc; method="cuhre")
    # if (method in ["quadgk"]) || !(method in integral_methods)
    #     error("Error: no such an integration method!")
    # elseif method == "quadgauss"
    #     res = quadgauss(x1 -> quadgauss(x2 -> quadgauss(x3 -> quadgauss(x4 -> quadgauss(x5 -> kumar_integrand_14(ampsq, proc.init_m, proc.final_m, [x1, x2, x3, x4, x5]), xxx, www), xxx, www), xxx, www), xxx, www), xxx, www)
    if (method in ["quadgk", "quadgauss"]) || !(method in integral_methods)
        error("Error: no such an integration method!")
    elseif method == "ccubature"
        res = ccubature(x -> kumar_integrand_14(ampsq, proc.init_m, proc.final_m, x), zeros(5), ones(5))[1]
    elseif method == "hcubature"
        res = hcubature(x -> kumar_integrand_14(ampsq, proc.init_m, proc.final_m, x), zeros(5), ones(5))[1]
    elseif method == "vegas"
        res = vegas((x, f) -> f[1] = kumar_integrand_14(ampsq, proc.init_m, proc.final_m, x), 5)[1]
    elseif method == "suave"
        res = suave((x, f) -> f[1] = kumar_integrand_14(ampsq, proc.init_m, proc.final_m, x), 5)[1]
    elseif method == "divonne"
        res = divonne((x, f) -> f[1] = kumar_integrand_14(ampsq, proc.init_m, proc.final_m, x), 5)[1]
    elseif method == "cuhre"
        res = cuhre((x, f) -> f[1] = kumar_integrand_14(ampsq, proc.init_m, proc.final_m, x), 5)[1]
    end

    return 1 / (2^10 * π^6) * res
end

#------------------------------15(kumar)----------------------------------------------------
# invariants: s1,s2,s3,u1,u2,u3,t2,t3;s34,s0,s4,u0,t1 
kumar15_s2p(m, m1, m2, s1, s2, u1) = s2 + m^2 + m1^2 + m2^2 - s1 - u1
kumar15_s3p(m, m1, m2, m3, s1, s3, u1, u2) = s3 + 2 * m^2 + m1^2 + m2^2 + m3^2 - s1 - u1 - u2
kumar15_t1p(m2) = m2^2
kumar15_t2p(m, m2, m3, u1, u2, t2) = t2 + m^2 + m2^2 + m3^2 - u1 - u2
kumar15_t3p(m, m2, m3, m4, u1, u2, u3, t3) = t3 + 2 * m^2 + m2^2 + m3^2 + m4^2 - u1 - u2 - u3
# adding "abs" to avoid the possible negative sqrt-root due to overflow of machine
kumar15_ξ2(m, m1, s2, t1, s2p, t1p) = ((m^2 + s2p - s2) * (m^2 + t1p - t1) - 2 * m^2 * (s2p + t1p - m1^2)) / sqrt(abs(λ(m^2, s2, s2p) * λ(m^2, t1, t1p)))
kumar15_ξ3(m, m1, s3, t2, s3p, t2p) = ((m^2 + s3p - s3) * (m^2 + t2p - t2) - 2 * m^2 * (s3p + t2p - m1^2)) / sqrt(abs(λ(m^2, s3, s3p) * λ(m^2, t2, t2p)))
kumar15_η2(m, m3, s2, s3, u2, s2p) = (2 * m^2 * (s2 + m3^2 - s3) - (m^2 + m3^2 - u2) * (m^2 + s2 - s2p)) / sqrt(abs(λ(m^2, s2, s2p) * λ(m^2, m3^2, u2)))
kumar15_η3(m, m4, m5, s3, u3, s3p) = (2 * m^2 * (s3 + m4^2 - m5^2) - (m^2 + m4^2 - u3) * (m^2 + s3 - s3p)) / sqrt(abs(λ(m^2, s3, s3p) * λ(m^2, m4^2, u3)))
kumar15_ω2(m, m3, u2, t1, t2, t1p) = (2 * m^2 * (t1 + m3^2 - t2) - (m^2 + m3^2 - u2) * (m^2 + t1 - t1p)) / sqrt(abs(λ(m^2, t1, t1p) * λ(m^2, m3^2, u2)))
kumar15_ω3(m, m4, u3, t2, t3, t2p) = (2 * m^2 * (t2 + m4^2 - t3) - (m^2 + m4^2 - u3) * (m^2 + t2 - t2p)) / sqrt(abs(λ(m^2, t2, t2p) * λ(m^2, m4^2, u3)))
kumar15_ζ2(ξ2, η2, ω2) = (ω2 - ξ2 * η2) / sqrt(abs((1 - ξ2^2) * (1 - η2^2)))
kumar15_ζ3(ξ3, η3, ω3) = (ω3 - ξ3 * η3) / sqrt(abs((1 - ξ3^2) * (1 - η3^2)))

function kumar_s1_15(init_v::T, final_v::T) where {T<:Vector{Float64}}
    m = init_v[1]
    m1, m2, m3, m4, m5 = final_v
    return @SVector [(m2 + m3 + m4 + m5)^2, (m - m1)^2]
end

function kumar_s2_15(init_v::T, final_v::T, s1::U) where {T<:Vector{Float64},U<:Float64}
    m = init_v[1]
    m1, m2, m3, m4, m5 = final_v
    return @SVector [(m3 + m4 + m5)^2, (sqrt(s1) - m2)^2]
end
function kumar_s3_15(init_v::T, final_v::T, s2::U) where {T<:Vector{Float64},U<:Float64}
    m = init_v[1]
    m1, m2, m3, m4, m5 = final_v
    return @SVector [(m4 + m5)^2, (sqrt(s2) - m3)^2]
end
function kumar_u1_15(init_v::T, final_v::T, s1::U, s2::U, s1p::U) where {T<:Vector{Float64},U<:Float64}
    m = init_v[1]
    m1, m2, m3, m4, m5 = final_v
    u1_min = m^2 + m2^2 - (s1 + m2^2 - s2) * (m^2 + s1 - s1p) / (2 * s1) - sqrt(abs(λ(s1, m2^2, s2) * λ(m^2, s1, s1p))) / (2 * s1)
    u1_max = m^2 + m2^2 - (s1 + m2^2 - s2) * (m^2 + s1 - s1p) / (2 * s1) + sqrt(abs(λ(s1, m2^2, s2) * λ(m^2, s1, s1p))) / (2 * s1)

    return @SVector [u1_min, u1_max]
end
function kumar_u2_15(init_v::T, final_v::T, s2::U, s3::U, s2p::U) where {T<:Vector{Float64},U<:Float64}
    m = init_v[1]
    m1, m2, m3, m4, m5 = final_v
    u2_min = m^2 + m3^2 - (s2 + m3^2 - s3) * (m^2 + s2 - s2p) / (2 * s2) - sqrt(abs(λ(s2, m3^2, s3) * λ(m^2, s2, s2p))) / (2 * s2)
    u2_max = m^2 + m3^2 - (s2 + m3^2 - s3) * (m^2 + s2 - s2p) / (2 * s2) + sqrt(abs(λ(s2, m3^2, s3) * λ(m^2, s2, s2p))) / (2 * s2)

    return @SVector [u2_min, u2_max]
end
function kumar_u3_15(init_v::T, final_v::T, s3::U, s3p::U) where {T<:Vector{Float64},U<:Float64}
    m = init_v[1]
    m1, m2, m3, m4, m5 = final_v
    u3_min = m^2 + m4^2 - (s3 + m4^2 - m5^2) * (m^2 + s3 - s3p) / (2 * s3) - sqrt(abs(λ(s3, m4^2, m5^2) * λ(m^2, s3, s3p))) / (2 * s3)
    u3_max = m^2 + m4^2 - (s3 + m4^2 - m5^2) * (m^2 + s3 - s3p) / (2 * s3) + sqrt(abs(λ(s3, m4^2, m5^2) * λ(m^2, s3, s3p))) / (2 * s3)

    return @SVector [u3_min, u3_max]
end

function kumar_t2_15(init_v::T, final_v::T, t1::U, u2::U, t1p::U, ξ2::U, η2::U) where {T<:Vector{Float64},U<:Float64}
    m = init_v[1]
    m1, m2, m3, m4, m5 = final_v
    t2_min = t1 + m3^2 - (m^2 + m3^2 - u2) * (m^2 + t1 - t1p) / (2 * m^2) + sqrt(abs(λ(m^2, m3^2, u2) * λ(m^2, t1, t1p))) / (2 * m^2) * (-ξ2 * η2 - sqrt(abs((1 - ξ2^2) * (1 - η2^2))))
    t2_max = t1 + m3^2 - (m^2 + m3^2 - u2) * (m^2 + t1 - t1p) / (2 * m^2) + sqrt(abs(λ(m^2, m3^2, u2) * λ(m^2, t1, t1p))) / (2 * m^2) * (-ξ2 * η2 + sqrt(abs((1 - ξ2^2) * (1 - η2^2))))
    return @SVector [t2_min, t2_max]
end

function kumar_t3_15(init_v::T, final_v::T, t2::U, u3::U, t2p::U, ξ3::U, η3::U) where {T<:Vector{Float64},U<:Float64}
    m = init_v[1]
    m1, m2, m3, m4, m5 = final_v
    t3_min = t2 + m4^2 - (m^2 + m4^2 - u3) * (m^2 + t2 - t2p) / (2 * m^2) + sqrt(abs(λ(m^2, m4^2, u3) * λ(m^2, t2, t2p))) / (2 * m^2) * (-ξ3 * η3 - sqrt(abs((1 - ξ3^2) * (1 - η3^2))))
    t3_max = t2 + m4^2 - (m^2 + m4^2 - u3) * (m^2 + t2 - t2p) / (2 * m^2) + sqrt(abs(λ(m^2, m4^2, u3) * λ(m^2, t2, t2p))) / (2 * m^2) * (-ξ3 * η3 + sqrt(abs((1 - ξ3^2) * (1 - η3^2))))
    return @SVector [t3_min, t3_max]
end

function kumar_integrand_15(ampsq::Function, init_v::T, final_v::T, x) where {T<:Vector{Float64}}
    m = init_v[1]
    m1, m2, m3, m4, m5 = final_v
    x1, x2, x3, x4, x5, x6, x7, x8 = x
    s1_span = kumar_s1_15(init_v, final_v)
    s1 = partition(x1, s1_span)
    s2_span = kumar_s2_15(init_v, final_v, s1)
    s2 = partition(x2, s2_span)
    s3_span = kumar_s3_15(init_v, final_v, s2)
    s3 = partition(x3, s3_span)
    s1p = m1^2
    u1_span = kumar_u1_15(init_v, final_v, s1, s2, s1p)
    u1 = partition(x4, u1_span)
    s2p = kumar15_s2p(m, m1, m2, s1, s2, u1)
    u2_span = kumar_u2_15(init_v, final_v, s2, s3, s2p)
    u2 = partition(x5, u2_span)
    s3p = kumar15_s3p(m, m1, m2, m3, s1, s3, u1, u2)
    u3_span = kumar_u3_15(init_v, final_v, s3, s3p)
    u3 = partition(x6, u3_span)
    t1 = u1
    t1p = kumar15_t1p(m2)
    ξ2 = kumar15_ξ2(m, m1, s2, t1, s2p, t1p)
    η2 = kumar15_η2(m, m3, s2, s3, u2, s2p)
    t2_span = kumar_t2_15(init_v, final_v, t1, u2, t1p, ξ2, η2)
    t2 = partition(x7, t2_span)
    t2p = kumar15_t2p(m, m2, m3, u1, u2, t2)
    ξ3 = kumar15_ξ3(m, m1, s3, t2, s3p, t2p)
    η3 = kumar15_η3(m, m4, m5, s3, u3, s3p)
    t3_span = kumar_t3_15(init_v, final_v, t2, u3, t2p, ξ3, η3)
    t3 = partition(x8, t3_span)
    ω2 = kumar15_ω2(m, m3, u2, t1, t2, t1p)
    ω3 = kumar15_ω3(m, m4, u3, t2, t3, t2p)
    ζ2 = kumar15_ζ2(ξ2, η2, ω2)
    ζ3 = kumar15_ζ3(ξ3, η3, ω3)

    s34_span = kumar15_quadratic(m, m1, m2, m3, m4, m5, s1, s2, s3, u1, u2, u3, t2, t3)


    res = (diff(s1_span) * diff(s2_span) * diff(s3_span) * diff(u1_span) * diff(u2_span) * diff(u3_span) * diff(t2_span) * diff(t3_span) / sqrt(abs(λ(m^2, s2, s2p) * λ(m^2, m3^2, u2) * λ(m^2, s3, s3p) * λ(m^2, m4^2, u3) * λ(m^2, t1, t1p) * λ(m^2, t2, t2p) * (1 - ξ2^2) * (1 - η2^2) * (1 - ζ2^2) * (1 - ξ3^2) * (1 - η3^2) * (1 - ζ3^2))) * (ampsq(s1, s2, s3, u1, u2, u3, t2, t3, s34_span[1]) + ampsq(s1, s2, s3, u1, u2, u3, t2, t3, s34_span[2])) / 2)

    if isnan(res) || isinf(res)
        error("Error: the integrand evaluates to NaN or Inf!")
    end

    return res * m^2
end

function ps15_kumar(ampsq::Function, proc::Proc; method="cuhre")
    if (method in ["quadgk", "quadgauss"]) || !(method in integral_methods)
        error("Error: no such an integration method!")
    elseif method == "ccubature"
        res = ccubature(x -> kumar_integrand_15(ampsq, proc.init_m, proc.final_m, x), zeros(8), ones(8))[1]
    elseif method == "hcubature"
        res = hcubature(x -> kumar_integrand_15(ampsq, proc.init_m, proc.final_m, x), zeros(8), ones(8))[1]
    elseif method == "vegas"
        res = vegas((x, f) -> f[1] = kumar_integrand_15(ampsq, proc.init_m, proc.final_m, x), 8)[1]
    elseif method == "suave"
        res = suave((x, f) -> f[1] = kumar_integrand_15(ampsq, proc.init_m, proc.final_m, x), 8)[1]
    elseif method == "divonne"
        res = divonne((x, f) -> f[1] = kumar_integrand_15(ampsq, proc.init_m, proc.final_m, x), 8)[1]
    elseif method == "cuhre"
        res = cuhre((x, f) -> f[1] = kumar_integrand_15(ampsq, proc.init_m, proc.final_m, x), 8)[1]
    end

    return 1 / (2^13 * π^9) * res
end