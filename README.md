# random_walk
## src/p_3.jlの解説
src/p_3.jlは

- 理論的な側面
    - 資料1：https://www2.math.kyushu-u.ac.jp/~hara/lectures/02/pr-grad-all.pdf
- 数値結果
    - 資料2：https://mathworld.wolfram.com/PolyasRandomWalkConstants.html

に基づき、3次元対称単純ランダムウォークの再帰確率（有限回数で再び原点に戻る確率）を数値的に求めることが目的である。

---

1. 1次元対称単純ランダムウォークにおいて、点 $` x `$ が次に点 $` y `$ に飛ぶ確率 $` p_{x, y} `$ は
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

2. 1次元の場合と同様に3次元対称単純ランダムウォークの再帰確率 $` p (3) `$ を考える。

点 $` \vec{x} = (x_1, x_2, x_3) `$ が次に点 $` \vec{y} = (y_1, y_2, y_3) `$ に飛ぶ確率 $` p_{\vec{x}, \vec{y}} `$ は
```math
p_{\vec{x}, \vec{y}} = \frac{1}{6} \delta_{(x_1, x_2, x_3), (x_1 - 1, x_2, x_3)} + \frac{1}{6} \delta_{(x_1, x_2, x_3), (x_1 + 1, x_2, x_3)} + \frac{1}{6} \delta_{(x_1, x_2, x_3), (x_1, x_2 - 1, x_3)} + \frac{1}{6} \delta_{(x_1, x_2, x_3), (x_1, x_2 + 1, x_3)} + \frac{1}{6} \delta_{(x_1, x_2, x_3), (x_1, x_2, x_3 - 1)} + \frac{1}{6} \delta_{(x_1, x_2, x_3), (x_1, x_2, x_3 + 1)}
```
と表せる（ $` \delta_{\vec{a}, \vec{b}} `$ はクロネッカーのデルタ）。

これを関数 $` f (\vec{x}) `$ で表すと
```math
f (\vec{x}) = 
\begin{cases}
    \frac{1}{6} & \vec{x} \in {\pm 1, 0, 0}, {0, \pm 1, 0}, {0, 0, \pm 1} \\
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
P_\ell (\vec{z}) = (f \ast f \ast \dots \ast f) (\vec{z}) = \int_{- \pi}^\pi \int_{- \pi}^\pi \int_{- \pi}^\pi (\hat{f} (\vec{k}))^\ell\,e^{i \vec{k} \cdot \vec{z}}\,\frac{d \vec{k}}{(2 \pi)^3}
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
u (3) \equiv \sum_{\ell = 0}^\infty P_\ell ((0, 0, 0)) = \int_{- \pi}^\pi \int_{- \pi}^\pi \int_{- \pi}^\pi \frac{1}{1 - \hat{f} (\vec{k})}\,\frac{d \vec{k}}{(2 \pi)^3}  = \int_{- \pi}^\pi \int_{- \pi}^\pi \int_{- \pi}^\pi \frac{3}{3 - \cos k_1 - \cos k_2 - \cos k_3}\,\frac{d \vec{k}}{(2 \pi)^3}
```
である。

3次元ランダムウォークの再帰確率 $` p (3) `$ は
```math
p (3) = 1 - \frac{1}{u (3)}
```
より、 $` u (3) `$ を求めれば、再帰確率 $` p (3) `$ が求められる。

---

3. Cubaを用いて $` u (3) `$ の値を求める。

そもそも $` u (3) `$ は


<div style="border-left: 4px solid #007acc; padding-left: 10px; margin: 10px 0;">
これは色付き縦棒を使ったブロックです。<br>
ここに注意点や補足などを書けます。
</div>

<div style="background: #e0f7fa; border-left: 4px solid #00acc1; padding: 10px; margin: 10px 0;">
📘 補足：このセクションは関連情報を含みます。
</div>
