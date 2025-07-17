# random_walk
## src/p_3.jlの解説
src/p_3.jlの内容は

- 理論的な側面
    - https://www.math.kyoto-u.ac.jp/~ichiro/lectures/07bpr.pdf
    - https://www2.math.kyushu-u.ac.jp/~hara/lectures/02/pr-grad-all.pdf
- 数値結果
    - https://mathworld.wolfram.com/PolyasRandomWalkConstants.html

に基づき、3次元対称単純ランダムウォークの再帰確率（有限回数で再び原点に戻る確率）を数値的に求めることが目的である。

大雑把な流れを説明すると、

1. 以下の3つを満たす確率変数$X_n$を考える：
   - $\mathbb{Z} \times \mathbb{Z} \times \mathbb{Z}$に値を取る
   - 独立かつ同分布（i.i.d. = independent idencally distributedという）
   - $P (X_n = (\pm 1, 0, 0)) = P (X_n = (0, \pm 1, 0)) = P (X_n = (0, 0, \pm 1)) = \frac{1}{6}$

   このとき、確率過程$S_n (= \sum_{k = 1}^n X_k)$（単純ランダムウォークという）に対し
   ```math
   P(\exist n \in \mathbb{Z}_+ \mathrm{s.t.} S_n = 0)
   ```
   を再起確率という。
