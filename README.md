# random_walk
## src/p_3.jlの解説
src/p_3.jlの内容は https://www2.math.kyushu-u.ac.jp/~hara/lectures/02/pr-grad-all.pdf に基づく。

大雑把な流れを説明すると、

1. 3次元対称単純ランダムウォークの再帰確率（有限回数で再び原点に戻る確率）を求めるために、
確率変数$X_n$は
- $\mathbb{Z} \times \mathbb{Z} \times \mathbb{Z}$に値を取る
- 独立
- 同分布
- $$
  P (X_n = (\pm 1, 0, 0)) = P (X_n = (0, \pm 1, 0)) = P (X_n = (0, 0, \pm 1)) = \frac{1}{6}
  $$
を満たす。