module random_walk

using HCubature

"""
    triple_integral()

三重積分
∫_{-π}^{π} ∫_{-π}^{π} ∫_{-π}^{π} 1/(3 - cos(x) - cos(y) - cos(z)) dx dy dz
を数値的に計算して返します。
"""
function triple_integral()
    f(x) = 1 / (3 - cos(x[1]) - cos(x[2]) - cos(x[3]))
    val_0, err_0 = hcubature(f, [-π, -π, -π], [π, π, π]; rtol=1e-8)
    val = (3 / (2 * π)^3) * val_0
    err = (3 / (2 * π)^3) * err_0
    return val, err
end

using QuadGK

function triple_integral_quadgk()
    f(x, y, z) = 1 / (3 - cos(x) - cos(y) - cos(z))
    # zについて積分
    inner(z, x, y) = f(x, y, z)
    # yについて積分
    middle(y, x) = quadgk(z -> inner(z, x, y), -π, π; rtol=1e-8)[1]
    # xについて積分
    val_0, err_0 = quadgk(x -> quadgk(y -> middle(y, x), -π, π)[1], -π, π; rtol=1e-8)
    val = (3 / (2 * π)^3) * val_0
    err = (3 / (2 * π)^3) * err_0
    return val, err
end

using Cuba

using Cuba

function test_pathological()
    function integrand(x)
        θ = 2π .* x .- π  # [0,1]^3 → [-π, π]^3 に変換
        denom = 3 - sum(cos.(θ))
        if abs(denom) < 1e-12 || !isfinite(denom)
            return 0.0
        else
            return (2π)^3 / denom  # ヤコビアン補正も忘れずに！
        end
    end

    cuba_func(x, f) = (f[1] = integrand(x) * 3 / (2π)^3; return 0)

    result = cuhre(cuba_func, 3, 1; rtol=1e-6)

    println("Computed: ", result.integral[1], " ± ", result.error[1])
end

test_pathological()



function test_integral()
    f(x) = (x[1]) + (x[2])^2
    val_0, err_0 = hcubature(f, [-π, -π], [π, π]; rtol=1e-8)
    val = 1 * val_0
    err = 1 * err_0
    return val, err
end

function p_3()
    val, err = triple_integral()
    value = 1 - 1 / val
    return value
end

end # module random_walk