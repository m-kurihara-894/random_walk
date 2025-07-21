module random_walk
# 多重積分のパッケージとしてHCubature, Cubature, Cubaなどが知られている。
# ここではCubaを使う。
using Cuba

# u(3) ≡ (3 / (2 π)^3 ) * ∫_-π^π ∫_-π^π ∫_-π^π 1 / (3 - cos(x) - cos(y) - cos(z)) dx dy dz を計算したい
# まずは、 ∫_[-π, π]^3 1 / (3 - cos(x) - cos(y) - cos(z)) dx dy dz を計算したい
function u3_cuba()
    function integrand(x)
        # Cubaの積分区間は [0,1]^3 なので、変数変換をして積分区間を対応させたい
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
    # u(3) は (3 / (2π)^3) * integrand(x) を積分したもの
    cuba_func(x, f) = (f[1] = integrand(x) * 3 / (2π)^3; return 0)

    result = cuhre(cuba_func, 3, 1; rtol=1e-6)
    val = result.integral[1]
    err = result.error[1]
    return val, err
end

# p(3) ≡ 1 - 1 / u(3) を計算する
function p3()
    val, err = u3_cuba()
    value = 1 - 1 / val
    return value
end

# -----------------------------------------------
# -----------------------------------------------

# 1次元ランダムウォークと熱拡散方程式との対応を見たい。
using Random, Distributions, Plots

export simulate_rw, plot_rw_vs_diffusion_save

"""
n_particles 個の粒子が n_steps ステップランダムウォークする。
各ステップで右に進む確率を p_right とする（デフォルトは対称 0.5）。
"""
function simulate_rw(n_particles::Int, n_steps::Int; Δx::Float64 = 1.0, p_right::Float64 = 0.5)
    @assert 0.0 ≤ p_right ≤ 1.0 "p_right must be between 0.0 and 1.0"

    # Bernoulli分布に従って 1（右）か -1（左）を選ぶ（つまり、確率 p で右へ、確率 1 - p で左に行く）
    bern = rand(Bernoulli(p_right), n_particles, n_steps)
    steps = 2 .* bern .- 1  # 1 → +1, 0 → -1 に変換
    positions = sum(steps, dims=2) .* Δx
    return vec(positions), p_right
end

"""
ランダムウォークのヒストグラムと熱拡散方程式の解析解の分布を重ねて描画し、ファイルへ保存。
"""
function plot_rw_vs_diffusion_save(positions::Vector{Float64}, n_steps::Int, p_right::Float64; Δx=1.0, Δt=1.0, filename="plot.png")
    # 移流拡散方程式（初期条件 u (x, 0) = δ (x) の場合、基本解は u (x, t) = 1 / √(4 π D t) exp (- (x - v t)^2 / (4 D t)) である）
    D = Δx^2 / (2Δt)
    σ² = 2D * n_steps * Δt
    v = (2 * p_right - 1) * (Δx / Δt)
    x = range(minimum(positions)-10, maximum(positions)+10, length=1000)
    pdf_gauss = pdf.(Normal(v * n_steps * Δt, sqrt(σ²)), x)

    plt = histogram(positions, bins=100, normalize=true, label="Random Walk", xlabel="x", ylabel="Probability")
    plot!(plt, x, pdf_gauss, lw=3, label="Diffusion Equation")

    savefig(plt, filename)
    return filename
end

end  # module
