# 内容

本コードは宇宙マイクロ波背景放射(CMB)の温度ゆらぎデータからミンコフスキー汎関数を計算するコードです。

ミンコフスキー汎関数は構造のトポロジーや形状を定量化する統計量であり、パワースペクトルなどの2点統計量では引き出すことができない高次相関(非ガウス性)の情報を引き出すことができます。
CMBゆらぎの統計はガウス統計によく従いますが、宇宙初期に起こったとされる急激な膨張期「インフレーション」の物理機構によっては、検出可能な非ガウス性が生まれる可能性があります。
ミンコフスキー汎関数は、わずかな非ガウス性のシグナルを捉えることができる統計量の1つです。

CMBのような2次元の場では、ミンコフスキー汎関数には表面積、曲率、オイラー数の3種類があり、これらを密度のしきい値の関数としてプロットします。

<img width="737" alt="スクリーンショット 2021-10-13 17 25 49" src="https://user-images.githubusercontent.com/86592645/137096609-e099b62b-a07b-47d0-bb37-7256ea0f4d8c.png">

上のグラフはCMBマップのシミュレーションから測定したミンコフスキー汎関数で、左から右の順に表面積(V0)、曲率(V1)、オイラー数(V2)を表しています。
一番上のパネルはガウス場におけるミンコフスキー汎関数の理論表式(赤線)とシミュレーションの測定結果の平均値(点)を比較したものです。
両者は非常によく一致しているように見えますが、シミュレーションには0.５%程度(f_NL=50)のわずかな非ガウス性が含まれています。
ミンコフスキー汎関数の非ガウス成分を取り出したものが下のパネルです。
非ガウス性が非常に小さい場合に適用できる摂動論な理論表式が導出されていて(Matusbara 2003)、シミュレーションの結果と誤差の範囲内で非常によく一致しています。

本コードでは、ガウス場の理論表式に加え、1次の摂動論の表式も合わせて出力できます。

なお、本コードは2次元の球面上の場として表されるものであればCMB温度ゆらぎ以外の場(例えば、CMBの偏光場や弱い重力レンズ場)に応用できます。

# References
- Primordial Non-Gaussianity and Analytical Formula for Minkowski Functionals of the Cosmic Microwave Background and Large-scale Structure  
Chiaki Hikage, Eiichiro Komatsu, Takahiko Matsubara  
Astrophys. J., Vol.653 No.1 (2006), pp.11-26

- Limits on Primordial Non-Gaussianity from Minkowski Functionals of the WMAP Temperature Anisotropies  
Chiaki Hikage, Takahiko Matsubara, Peter Coles, Michele Liguori, Frode K. Hansen, Sabino Matarrese  
Mon. Not. Roy. Astron. Soc., Vol.389 Issue 3 (2008), pp.1439-1446
