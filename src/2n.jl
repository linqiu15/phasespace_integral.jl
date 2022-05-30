Gram(x, y, z, u, v, w) = x^2 * y + x * y^2 + z^2 * u + z * u^2 + v * w^2 + v^2 * w + x * z * w + x * u * v + y * z * v + y * u * w - x * y * (z + u + v + w) - z * u * (x + y + v + w) - v * w * (x + y + z + u)

Δ4(masq, mi2sq, mi3sq, hsi1, hsi2, hsi3, ti1, ti2, ti3, si2) = 1 / 16 * det([2*masq masq+ti1-hsi1 masq+ti2-hsi2 masq+ti3-hsi3; masq+ti1-hsi1 2*ti1 ti1+ti2-mi2sq ti1+ti3-si2; masq+ti2-hsi2 ti1+ti2-mi2sq 2*ti2 ti2+ti3-mi3sq; masq+ti3-hsi3 ti1+ti3-si2 ti2+ti3-mi3sq 2*ti3])

V(hsi1, hsi2, hsi3, ti1, ti2, ti3, mi2sq, mi3sq, masq) = -1 / 8 * det([2*hsi2 masq+hsi2-ti2 hsi1+hsi2-mi2sq; masq+hsi2-ti2 2*masq masq+hsi1-ti1; hsi3+hsi2-mi3sq masq+hsi3-ti3 0])


function byckling_ti1_limits(hsi1, hsi2, ti2, mi2sq, masq)
    ti1_min = hsi1 + masq + (-hsi2 + ti2 - masq) * (hsi1 + hsi2 - mi2sq) / (2 * hsi2) - sqrt(abs(λ(hsi2, ti2, masq) * λ(hsi2, hsi1, mi2sq))) / (2 * hsi2)
    ti1_max = hsi1 + masq + (-hsi2 + ti2 - masq) * (hsi1 + hsi2 - mi2sq) / (2 * hsi2) + sqrt(abs(λ(hsi2, ti2, masq) * λ(hsi2, hsi1, mi2sq))) / (2 * hsi2)
    return @SVector [ti1_min, ti1_max]
end

function byckling_si2_limits(hsi1, hsi2, hsi3, ti1, ti2, ti3, mi2sq, mi3sq, masq)
    si2_min = hsi1 + hsi3 + 2 / λ(hsi2, ti2, masq) * (4 * V(hsi1, hsi2, hsi3, ti1, ti2, ti3, mi2sq, mi3sq, masq) - sqrt(abs(Gram(ti2, hsi3, hsi2, ti3, mi3sq, masq) * Gram(ti1, hsi2, hsi1, ti2, mi2sq, masq))))
    si2_max = hsi1 + hsi3 + 2 / λ(hsi2, ti2, masq) * (4 * V(hsi1, hsi2, hsi3, ti1, ti2, ti3, mi2sq, mi3sq, masq) + sqrt(abs(Gram(ti2, hsi3, hsi2, ti3, mi3sq, masq) * Gram(ti1, hsi2, hsi1, ti2, mi2sq, masq))))
    return @SVector [si2_min, si2_max]
end

function byckling_hsi2_limits(hsi1, hsi3, mi2, mi3)
    return @SVector [(sqrt(hsi1) + mi2)^2, (sqrt(hsi3) - mi3)^2]
end


##---------------------------------23(Byckling)--------------------------------------------
function byckling_integrandR3_23(ampsq::Function, init_v::T, final_v::T, s::U, x) where {T<:Vector{Float64},U<:Float64}
    ma, mb = init_v
    m1, m2, m3 = final_v
    x1, x2, x3, x4 = x
    hs2_span = byckling_hsi2_limits(m1^2, s, m2, m3)
    hs2 = partition(x1, hs2_span)
    t2_span = byckling_ti1_limits(hs2, s, mb^2, m3^2, ma^2)
    t2 = partition(x2, t2_span)
    t1_span = byckling_ti1_limits(m1^2, hs2, t2, m2^2, ma^2)
    t1 = partition(x3, t1_span)
    s2_span = byckling_si2_limits(m1^2, hs2, s, t1, t2, mb^2, m2^2, m3^2, ma^2)
    s2 = partition(x4, s2_span)

    res = diff(hs2_span) * diff(t2_span) * diff(t1_span) * diff(s2_span) / sqrt(abs(-Δ4(ma^2, m2^2, m3^2, m1^2, hs2, s, t1, t2, mb^2, s2))) * (2 * π) * ampsq(hs2, t2, t1, s2)

    if isnan(res) || isinf(res)
        error("Error: the integrand evaluates to NaN or Inf!")
    end

    return res / (32 * (2 * π)^5 * sqrt(λ(s, mb^2, ma^2)))
end

function ps23_byckling(ampsq::Function, proc::Proc, s::Float64; method="cuhre")
    if (method in ["quadgk", "quadgauss"]) || !(method in integral_methods)
        error("Error: no such an integration method!")
    elseif method == "ccubature"
        res = ccubature(x -> byckling_integrandR3_23(ampsq, proc.init_m, proc.final_m, s, x), zeros(4), ones(4))[1]
    elseif method == "hcubature"
        res = hcubature(x -> byckling_integrandR3_23(ampsq, proc.init_m, proc.final_m, s, x), zeros(4), ones(4))[1]
    elseif method == "vegas"
        res = vegas((x, f) -> f[1] = byckling_integrandR3_23(ampsq, proc.init_m, proc.final_m, s, x), 4)[1]
    elseif method == "suave"
        res = suave((x, f) -> f[1] = byckling_integrandR3_23(ampsq, proc.init_m, proc.final_m, s, x), 4)[1]
    elseif method == "divonne"
        res = divonne((x, f) -> f[1] = byckling_integrandR3_23(ampsq, proc.init_m, proc.final_m, s, x), 4)[1]
    elseif method == "cuhre"
        res = cuhre((x, f) -> f[1] = byckling_integrandR3_23(ampsq, proc.init_m, proc.final_m, s, x), 4)[1]
    end
    return res
end


##---------------------------------24(Byckling)--------------------------------------------
function byckling_integrandR4_24(ampsq::Function, init_v::T, final_v::T, s::U, x) where {T<:Vector{Float64},U<:Float64}
    ma, mb = init_v
    m1, m2, m3, m4 = final_v
    x1, x2, x3, x4, x5, x6, x7 = x
    hs3_span = byckling_hsi2_limits((m1 + m2)^2, s, m3, m4)
    hs3 = partition(x1, hs3_span)
    t3_span = byckling_ti1_limits(hs3, s, mb^2, m4^2, ma^2)
    t3 = partition(x2, t3_span)
    hs2_span = byckling_hsi2_limits(m1^2, hs3, m2, m3)
    hs2 = partition(x3, hs2_span)
    t2_span = byckling_ti1_limits(hs2, hs3, t3, m3^2, ma^2)
    t2 = partition(x4, t2_span)
    s3_span = byckling_si2_limits(hs2, hs3, s, t2, t3, mb^2, m3^2, m4^2, ma^2)
    s3 = partition(x5, s3_span)
    t1_span = byckling_ti1_limits(m1^2, hs2, t2, m2^2, ma^2)
    t1 = partition(x6, t1_span)
    s2_span = byckling_si2_limits(m1^2, hs2, hs3, t1, t2, t3, m2^2, m3^2, ma^2)
    s2 = partition(x7, s2_span)

    s234_span = byckling24_quadratic(s, ma, mb, m1, m2, m3, m4, hs3, t3, hs2, t2, s3, t1, s2)

    res = diff(hs3_span) * diff(t3_span) * diff(hs2_span) * diff(t2_span) * diff(s3_span) * diff(t1_span) * diff(s2_span) / sqrt(abs(Δ4(ma^2, m3^2, m4^2, hs2, hs3, s, t2, t3, mb^2, s3) * Δ4(ma^2, m2^2, m3^2, m1^2, hs2, hs3, t1, t2, t3, s2))) * (2 * π) * (ampsq(hs3, t3, hs2, t2, s3, t1, s2, s234_span[1]) + ampsq(hs3, t3, hs2, t2, s3, t1, s2, s234_span[2])) / 2

    if isnan(res) || isinf(res)
        error("Error: the integrand evaluates to NaN or Inf!")
    end

    return res / (4 * 8^2 * (2 * π)^8 * sqrt(λ(s, mb^2, ma^2)))
end

function ps24_byckling(ampsq::Function, proc::Proc, s::Float64; method="cuhre")
    if (method in ["quadgk", "quadgauss"]) || !(method in integral_methods)
        error("Error: no such an integration method!")
    elseif method == "ccubature"
        res = ccubature(x -> byckling_integrandR4_24(ampsq, proc.init_m, proc.final_m, s, x), zeros(7), ones(7))[1]
    elseif method == "hcubature"
        res = hcubature(x -> byckling_integrandR4_24(ampsq, proc.init_m, proc.final_m, s, x), zeros(7), ones(7))[1]
    elseif method == "vegas"
        res = vegas((x, f) -> f[1] = byckling_integrandR4_24(ampsq, proc.init_m, proc.final_m, s, x), 7)[1]
    elseif method == "suave"
        res = suave((x, f) -> f[1] = byckling_integrandR4_24(ampsq, proc.init_m, proc.final_m, s, x), 7)[1]
    elseif method == "divonne"
        res = divonne((x, f) -> f[1] = byckling_integrandR4_24(ampsq, proc.init_m, proc.final_m, s, x), 7)[1]
    elseif method == "cuhre"
        res = cuhre((x, f) -> f[1] = byckling_integrandR4_24(ampsq, proc.init_m, proc.final_m, s, x), 7)[1]
    end
    return res
end