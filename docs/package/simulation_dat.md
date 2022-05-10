# Simulation.dat
## 目次
1. [変数](#anchor1)
2. [関数](#anchor2)

## 前提条件
    from KudoCop.simulationdat import SimulationDat
    sdat = SimulationDat()
<a id="anchor1"></a>
# 変数

## len(sdat.atoms)

 全原子数


## atoms
 inputやdump.posの列情報についてpandasのDataFrameで与える。

覚えておくべきcolumnsについて

 - "id" : 原子のidについての一次元pandas.Seriesを格納
 - "type" : 原子のtypeについての一次元pandas.Seriesを格納
 - "mask" : 原子のマスク変数についての一次元pandas.Seriesを格納
 - "pos" : 原子の座標についての二次元pandas.Seriesを格納
 - "velo" : 原子の速度についての二次元pandas.Seriesを格納
 - "force" : 原子にかかる力のについての二次元pandas.Seriesを格納

注意
laich内での原子のidは1から始まり、sdat内では0から始まります。

つまり　laich内のid - 1 = SimulationDatの原子のindex　です。

例　laich内のidが7の原子の情報を取りたいとき
```
>>> sdat.atoms.iloc[6]
type     5.000
mask     0.000
x       66.745
y       46.264
z       11.566
Name: 6, dtype: float64
```

## cell

 セルの最大値についての配列 [x, y, z]


## atom_symbol_to_type
para.rdから読み込まれる。
キーを元素記号(str)、値を原子のタイプ(int)とした辞書
例）もしatom_symbol_to_type={"Fe":1, "Cu":2}ならば
sdat.atom_symbol_to_type["Fe"]は1となる

## atom_type_to_symbol
para.rdから読み込まれる。
キーを原子のタイプ(int)、値を元素記号(str)とした辞書

## atom_type_set
input.rd, dumppos, xyzから読み込まれる。
系内のatomのタイプ(int)が入ったset。

## atom_type_to_mass
キーを原子タイプ(int)、値を原子量(float)とした辞書。

## mpigrid
mpi分割についての配列(x,y,zの3方向分)
## ompgrid
omp分割についての配列(x,y,zの3方向分)
## cutoff
laich計算におけるカットオフ距離
## margin
laich計算におけるmargin
## total_step
laich計算におけるtotal_step
## file_step
laich計算におけるfile_step
## flagdecomp
Decompの有無。bool型
## flagconnect
結合情報記述の有無。自ら設定する必要あり。(デフォルトはFalse)  
使用例)sdat.flagconnect = Trueにしてからtrimming_particles()を実行することで，connect_listもtrimmingにより残った原子ものだけが残る．  
一般にconnect_listを使用した処理は重くなりがちだが，この機能を用いれば一部の原子，分子のみを対象に処理できるようになる．

## fix_info
inputにおけるfixについての配列
## move_info
inputにおけるmoveについての配列
## press_info
inputにおけるpressについての配列
## sumforce_info
inputにおけるsumforceについての配列
## thermofree_info
inputにおけるthermofreeについての配列

<a id="anchor2"></a>
# 関数
注意点 : すべての関数は返り値を持っていない。

## wrap_prticles()
原子の座標をセルの範囲内に収める。

## shift_particles()
第一引数(shift) : シフトさせたい数値。数値が一つならすべての座標がその数値分ずれる。要素が3つの配列ならば、x, y, zの座標が対応する配列の要素分ずれる。デフォルトではセルの最小値が0になるような相対移動が行われる。

## min_length()
第一引数(length) : 二つの原子間のベクトルを示す配列[x, y , z]。もしくはそのベクトルを行とした二次元配列。

周期的境界条件を考慮したベクトルにしてくれる。

## get_connect_list(cut_off)
第一引数(cut_off) : カットオフ距離

cut_offを上回るボンドオーダーを示す原子ペアについてconnect_listに格納する。

## read_para(input_filename)
para.rdを読み込む。
dumppos, dumpbond, input.rd, config.rd, xyzファイルよりも先に読み込んでください。

## read_input(input_filename)
input.rdを読み込む。

## read_config(input_filename)
config.rdを読み込む。

## read_xyz(input_filename)
xyzファイルを読み込む。

## read_dumppos(input_filename)
dumpposファイルを読み込む。

## read_dumpbond(input_filename)
dumpbondファイルを読み込む。


## to_input(output_filename)
input.rdを出力する。
第一引数(output_filename) : 保存するファイル名

## to_dumppos(output_filename, time_step=None, out_columns=None)
dumpposを出力する。
第一引数(output_filename) : 保存するファイル名

## to_xyz(output_filename, out_columns=None, structure_name=None)
xyzファイルを出力する。
第一引数(output_filename) : 保存するファイル名


<!-- 
## strain_cell
セルと原子位置をひずませる関数  
第一引数(strain_cell) : ひずみの大きさ(x,y,z) -->


<!-- ## concate_particles
新たな原子を追加する関数  
第一引数(concate_particles) : sdat.particlesに追加したい情報を有した辞書。

sdat.particles の末尾に第一引数の要素を追加する。第一引数のキーがsdat.particlesに存在しない場合、sdat.particlesに値をゼロ配列として該当キーを生成する(下の例だとmask変数になる)。ちなみに、sdat.particlesのキーが第一引数に存在しない場合は対応してません。
例）
sdat.particles = {"pos" : [[0,0,0],[1,1,1]]}の場合。
sdat.concate_particles({"pos":[[2,2,2],[3,3,3]], "mask" : [1,1]})
とすると、
sdat.particles : {"pos":[[0,0,0],[1,1,1],[2,2,2],[3,3,3]], "mask" : [0,0,1,1]}
となる。 -->

<!-- ## replicate_particles
第一引数(replicate_particles) : 等倍したいx,y,z方向へのリスト．

周期境界方向にセルを等倍し，粒子をコピーする．
この関数を呼び出した際，関数内では[self.shift_particles()](#shift_particles)が呼び出され，セルの最小値が0になるように相対移動が行われる．

*例) x,y,z方向にそれぞれセルを2倍にしてコピーする*
```py
sdat.replicate_particles([2,2,2]) -->
<!-- ```

## trimming_particles
flagの原子のみを抽出する関数  
第一引数(flag) : True か False を格納した配列。(一次元)
sdat.particlesの中から、第一引数のTrueのインデックスだけ取り出す。

第二引数(reindex): True→IDを振りなおす． False→IDを振りなおさない(トリミングした原子のIDは空行になる) <br>
                   普通に使うときはTrueでよい．どのIDを持つ原子が切り取られたかを確認したいときなどにFalseを使うといい.

例）
モデルをsinカーブに切り取りたかったら
flag = sdat.particles["pos"][0] < np.sin(...)とか？
いろいろできそうですね～

<!-- ## set_decomp
decompされたセル長を設定する関数 -->