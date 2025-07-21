# random_walk
- 協力者：[k-pine-su](https://github.com/k-pine-su/randomwalk)

## 概要
- [3次元対称単純ランダムウォークの再帰確率 $` p (3) `$ の数値計算](https://github.com/m-kurihara-894/random_walk/tree/main?tab=readme-ov-file#p3%E3%81%AE%E8%A7%A3%E8%AA%AC)
  - 理論式の導出（[1次元](https://github.com/m-kurihara-894/random_walk/tree/main?tab=readme-ov-file#1%E6%AC%A1%E5%85%83%E5%AF%BE%E7%A7%B0%E5%8D%98%E7%B4%94%E3%83%A9%E3%83%B3%E3%83%80%E3%83%A0%E3%82%A6%E3%82%A9%E3%83%BC%E3%82%AF%E3%81%A8fourier%E5%A4%89%E6%8F%9B)、[3次元](https://github.com/m-kurihara-894/random_walk/tree/main?tab=readme-ov-file#3%E6%AC%A1%E5%85%83%E5%AF%BE%E7%A7%B0%E5%8D%98%E7%B4%94%E3%83%A9%E3%83%B3%E3%83%80%E3%83%A0%E3%82%A6%E3%82%A9%E3%83%BC%E3%82%AF%E3%81%A8fourier%E5%A4%89%E6%8F%9B)）
  - [数値計算のアイディア](https://github.com/m-kurihara-894/random_walk/tree/main?tab=readme-ov-file#cuba%E3%82%92%E7%94%A8%E3%81%84%E3%81%9F-u-3-%E3%81%AE%E6%95%B0%E5%80%A4%E8%A7%A3%E3%81%A8%E3%82%A2%E3%82%A4%E3%83%87%E3%82%A3%E3%82%A2)
  - [実行結果](https://github.com/m-kurihara-894/random_walk/tree/main?tab=readme-ov-file#p-3-%E3%81%AE%E6%95%B0%E5%80%A4%E8%A7%A3)
- [1次元単純ランダムウォークのシミュレーションと熱拡散方程式の解との比較](https://github.com/m-kurihara-894/random_walk/tree/main?tab=readme-ov-file#simulate_rw-plot_rw_vs_diffusion_save%E3%81%AE%E4%BD%BF%E3%81%84%E6%96%B9)
  - [理論式の導出](https://github.com/m-kurihara-894/random_walk/tree/main?tab=readme-ov-file#1%E6%AC%A1%E5%85%83%E5%8D%98%E7%B4%94%E3%83%A9%E3%83%B3%E3%83%80%E3%83%A0%E3%82%A6%E3%82%A9%E3%83%BC%E3%82%AF%E3%81%A8%E7%86%B1%E6%8B%A1%E6%95%A3%E6%96%B9%E7%A8%8B%E5%BC%8F)
  - [実行結果](https://github.com/m-kurihara-894/random_walk/tree/main?tab=readme-ov-file#%E4%BD%BF%E3%81%84%E6%96%B9)

をjuliaで実装した。

---

---

## p3の解説
p3は

- 理論的な側面
    - 資料1：https://www2.math.kyushu-u.ac.jp/~hara/lectures/02/pr-grad-all.pdf
- 数値結果
    - 資料2：https://mathworld.wolfram.com/PolyasRandomWalkConstants.html

に基づき、3次元対称単純ランダムウォークの再帰確率（有限回数で再び原点に戻る確率） $` p (3) `$ を数値的に求めることが目的のプログラムである。

---

### 1次元対称単純ランダムウォークとFourier変換
1次元対称単純ランダムウォークにおいて、点 $` x `$ が次に点 $` y `$ に飛ぶ確率 $` p_{x, y} `$ は
```math
p_{x, y} = \frac{1}{2} \delta_{x, x - 1} + \frac{1}{2} \delta_{x, x + 1}
```
と表せる（ $` \delta_{a, b} `$ はクロネッカーのデルタ）。

これを関数 $` f (x) `$ で表すと
```math
f (x) = 
\begin{cases}
    \frac{1}{2} & x \in \pm 1 \\
    0 & \mathrm{otherwise}
  \end{cases}
```
と表せる。

ここで、格子点 $` x \in \mathbb{Z} `$ で定義された関数 $` f (x) `$ と $` k \in [- \pi, \pi] `$ に対するFourier変換（離散時間Fourier変換：DTFTという）：
```math
\hat{f} (k) \equiv \sum_{x \in \mathbb{Z}} e^{- i k x} f (x), \quad f (x) \equiv \int_{- \pi}^\pi \hat{f} (k)\,e^{i k x}\,\frac{d k}{2 \pi}
```
における畳み込み：
```math
(f \ast g) (x) \equiv \sum_{y \in \mathbb{Z}} f (x - y) g (y)
```
を用いると、 $` f \ast g `$ のFourier変換 $` \widehat{(f \ast g)} `$ は
```math
\widehat{(f \ast g)} (k) = \hat{f} (k) \hat{g} (k)
```
を満たす。

原点から出発したランダムウォークが時刻 $` \ell `$ に点 $` z `$ にいる確率 $` P_\ell (z) `$ は $` \ell `$ 個の $` f (x) = p_{0, x} `$ の畳み込み：
```math
P_\ell (z) = (f \ast f \ast \dots \ast f) (z) = \int_{- \pi}^\pi (\hat{f} (k))^\ell\,e^{i k z}\,\frac{d k}{2 \pi}
```
で表せて（畳み込みは今の場合においては、原点から $` \ell `$ 回の操作で点 $` z `$ に到着できる全ての可能性を足し上げていると思えばよい）、特に $` z = 0 `$ （原点に戻る）の場合は
```math
P_\ell (0) = \int_{- \pi}^\pi (\hat{f} (k))^\ell\,\frac{d k}{2 \pi}
```
である。

 $` f (x) `$ のFourier変換は
```math
\hat{f} (k) = \frac{1}{2} e^{- i k \cdot 1} + \frac{1}{2} e^{- i k \cdot (- 1)} = \cos k
```
より
```math
\sum_{\ell = 0}^\infty P_\ell (0) = \int_{- \pi}^\pi \frac{1}{1 - \hat{f} (k)}\,\frac{d k}{2 \pi}  = \int_{- \pi}^\pi \frac{1}{1 - \cos k}\,\frac{d k}{2 \pi} = \infty
```
である。

1次元ランダムウォークの再帰確率 $` p (1) `$ は
```math
p (1) = 1 - \frac{1}{\sum_{\ell = 0}^\infty P_\ell (0)} = 1
```
である（最初の等号の導出は資料1の3.2.3を参照）。

---

### 3次元対称単純ランダムウォークとFourier変換
1次元の場合と同様に3次元対称単純ランダムウォークの再帰確率 $` p (3) `$ を考える。

点 $` \vec{x} = (x_1, x_2, x_3) `$ が次に点 $` \vec{y} = (y_1, y_2, y_3) `$ に飛ぶ確率 $` p_{\vec{x}, \vec{y}} `$ は
```math
p_{\vec{x}, \vec{y}} = \frac{1}{6} \delta_{(x_1, x_2, x_3), (x_1 - 1, x_2, x_3)} + \frac{1}{6} \delta_{(x_1, x_2, x_3), (x_1 + 1, x_2, x_3)} + \frac{1}{6} \delta_{(x_1, x_2, x_3), (x_1, x_2 - 1, x_3)} + \frac{1}{6} \delta_{(x_1, x_2, x_3), (x_1, x_2 + 1, x_3)} + \frac{1}{6} \delta_{(x_1, x_2, x_3), (x_1, x_2, x_3 - 1)} + \frac{1}{6} \delta_{(x_1, x_2, x_3), (x_1, x_2, x_3 + 1)}
```
と表せる（ $` \delta_{\vec{a}, \vec{b}} `$ はクロネッカーのデルタ）。

これを関数 $` f (\vec{x}) `$ で表すと
```math
f (\vec{x}) = 
\begin{cases}
    \frac{1}{6} & \vec{x} \in (\pm 1, 0, 0), (0, \pm 1, 0), (0, 0, \pm 1) \\
    0 & \mathrm{otherwise}
  \end{cases}
```
と表せる。

格子点 $` \vec{x} \in \mathbb{Z} \times \mathbb{Z} \times \mathbb{Z} `$ で定義された関数 $` f (\vec{x}) `$ と $` \vec{k} \in [- \pi, \pi] \times [- \pi, \pi] \times [- \pi, \pi] `$ に対するFourier変換（ $` \vec{k} = (k_1, k_2, k_3) `$ とした）：
```math
\hat{f} (\vec{k}) \equiv \sum_{\vec{x} \in \mathbb{Z} \times \mathbb{Z} \times \mathbb{Z}} e^{- i \vec{k} \vec{x}} f (\vec{x}), \quad f (\vec{x}) \equiv \int_{- \pi}^\pi \int_{- \pi}^\pi \int_{- \pi}^\pi \hat{f} (\vec{k})\,e^{i \vec{k} \cdot \vec{x}}\,\frac{d^3 \vec{k}}{(2 \pi)^3}
```
における畳み込み：
```math
(f \ast g) (\vec{x}) \equiv \sum_{\vec{y} \in \mathbb{Z} \times \mathbb{Z} \times \mathbb{Z}} f (\vec{x} - \vec{y}) g (\vec{y})
```
を用いると、
原点から出発したランダムウォークが時刻 $` \ell `$ に点 $` \vec{z} = (z_!, z_2, z_3) `$ にいる確率 $` P_\ell (\vec{z}) `$ は $` \ell `$ 個の $` f (\vec{x}) = p_{(0, 0, 0), \vec{x}} `$ の畳み込み：
```math
P_\ell (\vec{z}) = (f \ast f \ast \dots \ast f) (\vec{z}) = \int_{- \pi}^\pi \int_{- \pi}^\pi \int_{- \pi}^\pi (\hat{f} (\vec{k}))^\ell\,e^{i \vec{k} \cdot \vec{z}}\,\frac{d^3 \vec{k}}{(2 \pi)^3}
```
で表せて、特に $` \vec{z} = (0, 0, 0) `$ の場合は
```math
P_\ell ((0, 0, 0)) = \int_{- \pi}^\pi \int_{- \pi}^\pi \int_{- \pi}^\pi (\hat{f} (\vec{k}))^\ell\,\frac{d^3 \vec{k}}{(2 \pi)^3}
```
である。

 $` f (\vec{x}) `$ のFourier変換は
```math
\hat{f} (\vec{k}) = \frac{1}{6} e^{- i (k_1, k_2, k_3) \cdot (- 1, 0, 0)} + \frac{1}{6} e^{- i (k_1, k_2, k_3) \cdot (1, 0, 0)} + \frac{1}{6} e^{- i (k_1, k_2, k_3) \cdot (0, - 1, 0)} + \frac{1}{6} e^{- i (k_1, k_2, k_3) \cdot (0, 1, 0)} + \frac{1}{6} e^{- i (k_1, k_2, k_3) \cdot (0, 0, - 1)} + \frac{1}{6} e^{- i (k_1, k_2, k_3) \cdot (0, 0, 1)} = \frac{1}{3} \cos k_1 + \frac{1}{3} \cos k_2 + \frac{1}{3} \cos k_3
```
より
```math
u (3) \equiv \sum_{\ell = 0}^\infty P_\ell ((0, 0, 0)) = \int_{- \pi}^\pi \int_{- \pi}^\pi \int_{- \pi}^\pi \frac{1}{1 - \hat{f} (\vec{k})}\,\frac{d^3 \vec{k}}{(2 \pi)^3}  = \int_{- \pi}^\pi \int_{- \pi}^\pi \int_{- \pi}^\pi \frac{3}{3 - \cos k_1 - \cos k_2 - \cos k_3}\,\frac{d^3 \vec{k}}{(2 \pi)^3}
```
である。

3次元ランダムウォークの再帰確率 $` p (3) `$ は
```math
p (3) = 1 - \frac{1}{u (3)}
```
より、 $` u (3) `$ を求めれば、再帰確率 $` p (3) `$ が求められる。

---

### Cubaを用いた $` u (3) `$ の数値解とアイディア
Cubaを用いて $` u (3) `$ の値を求めるが、非積分関数が原点において発散してしまう。

しかしながら、今回は原点近傍の値を捨てても積分結果に及ぼす影響が小さいことが分かる。

原点を中心とする半径 $` \rho `$ 以下の球 $` B_\rho `$ を考える。

積分範囲について
```math
[- \pi, \pi] \times [- \pi, \pi] \times [- \pi, \pi] = B_\rho \cup ([- \pi, \pi] \times [- \pi, \pi] \times [- \pi, \pi] \setminus B_\rho)
```
と分けて、 $` B_\rho `$ において非積分関数を積分してみると
```math
\int_{B_\rho} \frac{3}{3 - \cos k_1 - \cos k_2 - \cos k_3}\,\frac{d^3 \vec{k}}{(2 \pi)^3} \leq \int_{B_\rho} \frac{6}{k_1^2 + k_2^2 + k_3^2}\,\frac{d^3 \vec{k}}{(2 \pi)^3} = \frac{3 \rho}{\pi^2}
```
となるため、十分に小さい $` \rho `$ を取れば、原点近傍における積分の寄与は無視できる。

そのため、p3においては、場合分けで発散する点の値を0としている。

また、Cubaは積分区間が $` [0, 1] \times [0, 1] \times [0, 1] `$ となっているため、適当な変数変換を行い、後でJacobianを乗じている。

---

###  $` p (3) `$ の数値解
p3を実行する：
```julia
(@v1.11) pkg> activate .

julia> import random_walk

julia> random_walk.u3_cuba()
(1.5163840716236674, 1.5121509978460541e-6) # (val, err)

julia> random_walk.p3()
0.3405364651916645
```
 $` p (3) \approx 0.34 `$ となり、再帰確率が1ではないことが分かる。

---

---

## simulate_rw, plot_rw_vs_diffusion_saveの解説
plot_rw_vs_diffusion_saveはsimulate_rwでシミュレーションした1次元ランダムウォークを熱拡散方程式の基本解と比較することが目的のプログラムである。

---

### 1次元単純ランダムウォークと熱拡散方程式
1ステップ後（ $` \Delta t `$ 後）に右へ $` \Delta x `$ 進む確率が $` p `$ 、左へ $` \Delta x `$ 進む確率が $` 1 - p `$ である1次元単純ランダムウォーク：
```math
u (x, t + \Delta t) = p\,u (x + \Delta x, t) + (1 - p)\,u (x -  \Delta x, t)
```
を考える。

ここで、 $` \Delta t `$ に関しては1次、$` \Delta x `$ に関しては2次までTaylor展開すると
```math
\frac{\partial u}{\partial t} + (1 - 2p) \frac{\Delta x}{\Delta t} \frac{\partial u}{\partial x} = \frac{(\Delta x)^2}{2 \Delta t} \frac{\partial^2 u}{\partial x^2}
```
と表せ、移流速度 $` v \equiv (1 - 2p) \frac{\Delta x}{\Delta t} `$ 、拡散係数 $`D \equiv \frac{(\Delta x)^2}{2 \Delta t} `$ を導入すると、熱拡散方程式：
```math
\frac{\partial u}{\partial t} + v \frac{\partial u}{\partial x} = D \frac{\partial^2 u}{\partial x^2}
```
を得る。
---

### 使い方
```julia
(@v1.11) pkg> activate .

julia> import random_walk

import random_walk

using random_walk

positions, p_right = simulate_rw(10000, 1000; p_right=0.55)  # 右に偏る確率が0.55のランダムウォーク

file = random_walk.plot_rw_vs_diffusion_save(positions, 1000, p_right; filename="rw_vs_diff.png") # 熱拡散方程式の解との比較を表したグラフ

println("保存しました：", file)
```
でランダムウォークと熱拡散方程式の比較ができる。