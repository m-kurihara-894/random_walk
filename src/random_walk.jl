module random_walk

using Cuba

function u3_cuba()
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
    val = result.integral[1]
    err = result.error[1]
    return val, err
end

function p_3()
    val, err = u3_cuba()
    value = 1 - 1 / val
    return value
end

end