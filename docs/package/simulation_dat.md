# Simulation.dat
## 目次
1. [変数](#anchor1)
2. [関数](#anchor2)

## 前提条件
    from KudoCop.simulationdat import SimulationDat
    sdat = SimulationDat()
<a id="anchor1"></a>
# 変数



## atoms
 inputやdump.posの列情報についてpandasのDataFrameで与える。

覚えておくべきcolumnsについて

 - "id" : 原子のidについての一次元pandas.Seriesを格納
 - "type" : 原子のtypeについての一次元pandas.Seriesを格納
 - "mask" : 原子のマスク変数についての一次元pandas.Seriesを格納
 - "x" : 原子のx座標についての一次元pandas.Seriesを格納
 - "y" : 原子のy座標についての一次元pandas.Seriesを格納
 - "z" : 原子のz座標についての一次元pandas.Seriesを格納
 - "vx" : 原子のx方向の速度についての一次元pandas.Seriesを格納
 - "vy" : 原子のy方向の速度についての一次元pandas.Seriesを格納
 - "vz" : 原子のz方向の速度についての一次元pandas.Seriesを格納
 - "fx" : 原子にかかるx方向の力についての一次元pandas.Seriesを格納
 - "fy" : 原子にかかるy方向の力についての一次元pandas.Seriesを格納
 - "fz" : 原子にかかるz方向の力についての一次元pandas.Seriesを格納

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

## get_total_atoms()

 全原子数を返す

## wrap_prticles()
原子の座標をセルの範囲内に収める。


## get_connect_list(cut_off)
第一引数(cut_off) : カットオフ距離

cut_offを上回るボンドオーダーを示す原子ペアについてconnect_listに格納し、connect_listを返す

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

## to_dumppos(output_filename, time_step, out_columns)
dumpposを出力する。
第一引数(output_filename) : 保存するファイル名
time_step : タイムステップ、デフォルトは0
out_columns : 出力する列を指定

## to_xyz(output_filename, out_columns, structure_name)
xyzファイルを出力する。
第一引数(output_filename) : 保存するファイル名
out_columns : 出力する列を指定
structure_name : 構造の名前

## replicate_particles([x, y, z])
第一引数(replicate_particles) : 等倍したいx,y,z方向へのリスト．

周期境界方向にセルを等倍し，粒子をコピーする．
```
*例) x,y,z方向にそれぞれセルを2倍にしてコピーする*
sdat.replicate_particles([2,2,2])
```
## concat_atoms(outer_sdat)
sdatにouter_sdatを結合する
```
*例) sdatにouter_sdatを結合する*
sdat.concat_atoms(outer_sdat)
```

## count_bonds(cut_off)
結合の強さがcut_off以上の結合数をカウントする。
```
>>>sdat.count_bonds(0.3)
{(2, 2): 0, (2, 3): 3967, (2, 4): 0, (2, 5): 36, (3, 3): 0, (3, 4): 0, (3, 5): 791, (4, 4): 3344, (4, 5): 35466, (5, 5): 8090}]
>>> sdat.atom_type_to_symbol
{1: 'C', 2: 'H', 3: 'O', 4: 'N', 5: 'Si'}
```
この場合原子2と原子3の結合数が3967個ある
## count_mols(self, cut_off, lower_mol_limit=1, upper_mol_limit=10)
cut_off以上の結合を結合とみなし、分子の数を数える。
分子内の原子の数がlower_mol_limit以上upper_mol_limit以下の分子のみがカウントされる。
```
>>> sdat.count_mols(0.3, upper_mol_limit=100)
{(0, 0, 8, 4, 1, 1): 1, (0, 0, 0, 0, 1, 0): 4, (0, 0, 2, 1, 0, 0): 1225, (0, 0, 3, 1, 0, 0): 1}
>>> sdat.atom_type_to_symbol
{1: 'C', 2: 'H', 3: 'O', 4: 'N', 5: 'Si'}
```
この場合H2O分子が1225個ある
## get_atom_idx_from_mol(cuf_off, target_mol)
target_molの個数を数える。target_molにはタプルを渡す。
```
例水分子の個数を数える。
>>> sdat.get_atom_idx_from_mol(0.3, target_mol=(0, 0, 2, 1, 0, 0))[:10]
[[28313, 28315, 28314], [28325, 28327, 28326], [28334, 28336, 28335], [28340, 28342, 28341], [28343, 28345, 28344], [28349, 28351, 28350], [28355, 28357, 28356], [28358, 28360, 28359], [28361, 28363, 28362], [28364, 28366, 28365]]
```



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