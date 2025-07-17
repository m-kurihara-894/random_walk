module random_walk
# 多重積分のパッケージとしてHCubature, Cubature, Cubaなどが知られている。
# ここではCubaを使う。
using Cuba

# (3 / (2 π)^3 ) * ∫_-π^π ∫_-π^π ∫_-π^π 1 / (3 - cos(x) - cos(y) - cos(z)) dx dy dzを計算したい
# まずは、∫_-π^π ∫_-π^π ∫_-π^π 1 / (3 - cos(x) - cos(y) - cos(z)) dx dy dzを計算したい
function u3_cuba()
    function integrand(x)
        # Cubaの積分区間は[0,1]^3なので、変数変換をして積分区間を対応させたい
        # [0,1]^3 → [-π, π]^3 に変換
        θ = 2π .* x .- π
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