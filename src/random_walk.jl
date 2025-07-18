module random_walk
# 多重積分のパッケージとしてHCubature, Cubature, Cubaなどが知られている。
# ここではCubaを使う。
using Cuba

# u(3) ≡ (3 / (2 π)^3 ) * ∫_-π^π ∫_-π^π ∫_-π^π 1 / (3 - cos(x) - cos(y) - cos(z)) dx dy dzを計算したい
# まずは、∫_[-π, π]^3 1 / (3 - cos(x) - cos(y) - cos(z)) dx dy dzを計算したい
function u3_cuba()
    function integrand(x)
        # Cubaの積分区間は[0,1]^3なので、変数変換をして積分区間を対応させたい
        # [0,1]^3 → [-π, π]^3 に変換
        θ = 2π .* x .- π
        # 原点で発散があるので、原点近傍では0にする
        denom = 3 - sum(cos.(θ))
        if abs(denom) < 1e-12 || !isfinite(denom)
            return 0.0
        else
            return (2π)^3 / denom  # ヤコビアン補正は(2π)^3
        end
    end
    # u(3)は(3 / (2π)^3) * integrand(x)を積分したもの
    cuba_func(x, f) = (f[1] = integrand(x) * 3 / (2π)^3; return 0)

    result = cuhre(cuba_func, 3, 1; rtol=1e-6)
    val = result.integral[1]
    err = result.error[1]
    return val, err
end

# p(3) ≡ 1 - 1 / u(3)を計算する
function p3()
    val, err = u3_cuba()
    value = 1 - 1 / val
    return value
end

# -----------------------------------------------

# 1次元ランダムウォークと熱拡散方程式との対応を見たい。
using Random, Distributions, Plots

export simulate_rw, plot_rw_vs_diffusion_save

function simulate_rw(n_particles::Int, n_steps::Int; Δx::Float64 = 1.0)
    steps = rand([-1, 1], n_particles, n_steps)
    positions = sum(steps, dims=2) .* Δx
    return vec(positions)
end

"""
    plot_rw_vs_diffusion_save(positions::Vector{Float64}, n_steps::Int; Δx=1.0, Δt=1.0, filename="plot.png")

ランダムウォークのヒストグラムと熱拡散方程式の解析解の分布を重ねて描画し、
表示はせずにファイルへ保存します。
"""
function plot_rw_vs_diffusion_save(positions::Vector{Float64}, n_steps::Int; Δx=1.0, Δt=1.0, filename="plot.png")
    D = Δx^2 / (2Δt)
    σ² = 2D * n_steps * Δt
    x = range(minimum(positions)-10, maximum(positions)+10, length=1000)
    pdf_gauss = pdf.(Normal(0, sqrt(σ²)), x)

    plt = histogram(positions, bins=100, normalize=true, label="Random Walk", xlabel="x", ylabel="Probability")
    plot!(plt, x, pdf_gauss, lw=3, label="Diffusion Equation")

    savefig(plt, filename)  # ファイルに保存
    return filename
end

end # module