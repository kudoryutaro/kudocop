# Simulationdat
# 目次
1. [変数](#anchor1)
2. [関数](#anchor2)

## import
```python
from KudoCop import SimulationDat
sdat = SimulationDat()
# (推奨)初めにpara.rdを読み込む
sdat.import_para('para.rd')
# ファイルの読み込み
sdat.import_input('input.rd')
```
<a id="anchor1"></a>
# 変数



## atoms
 inputやdump.posの列情報についてpandasのDataFrameで与える。

覚えておくべきcolumnsについて

 - "id" : 原子のid
- "type" : 原子のtype
 - "mask" : 原子のマスク変数
 - "x", "y", "z" : 原子のx, y, z座標
 - "vx", "vy", "vz" : 原子のx, y, z方向の速度
 - "fx", "fy", "fz" : 原子にかかるx, y, z方向の力
 - "q" : 原子の持つ電荷
```python
>>> sdat.atoms
       type       x       y       z     vx     vy     vz      q  mask
0         5  36.205  47.491   0.093  0.000  0.000  0.000  0.585     4
1         5  86.205  47.491   0.093  0.000  0.000  0.000  0.585     4
2         5  47.122  33.792  43.391 -0.001 -0.002 -0.005  0.902     5
3         5  47.133  33.647  74.308  0.000  0.000 -0.000  0.879     2
4         5  97.133  33.647  74.308  0.000  0.000 -0.000  0.878     2
...     ...     ...     ...     ...    ...    ...    ...    ...   ...
34338     2  72.366  22.848  43.915  0.008  0.003 -0.025  0.414     5
34339     3  71.804  23.420  43.409 -0.001  0.003  0.013 -0.757     5
34340     2  21.110  48.754  41.645 -0.014  0.018 -0.022  0.337     5
34341     2  22.073  48.246  42.493  0.007  0.000  0.013  0.352     5
34342     3  21.942  48.335  41.588  0.001  0.006  0.002 -0.807     5

[34343 rows x 9 columns]
```
注意
laich内での原子のidは1から始まり、SimulationDat内では0から始まります。
つまり laich内のid - 1 = SimulationDatの原子のindexです。

例　laich内のidが7の原子の情報を取りたいとき
```python
>>> sdat.atoms.iloc[6]
type     5.000
mask     0.000
x       66.745
y       46.264
z       11.566
Name: 6, dtype: float64
```

atomsにあるcolumnを取り出す時はatomsを省略して書くことが出来る。
```python
>>> sdat[['x', 'y', 'z']]
            x       y       z
0      16.849  20.389  20.174
1      32.722  10.533   4.074
2      23.587  11.436  13.742
3       3.904  28.283   0.572
4      48.649  31.810  28.123
...       ...     ...     ...
37354  96.081  16.863  25.520
37355  95.477  17.429  25.990
37356  81.027   4.444  43.712
37357  81.443   3.474  39.168
37358  82.514   3.136  40.906
```
ただし、代入は出来ない


## cell
```python
>>> sdat.cell
[100.0, 50.0, 110.0]
```
 セルの最大値についての配列 [x, y, z]


## atom_symbol_to_type
```python
>>> sdat.atom_symbol_to_type
{'C': 1, 'H': 2, 'O': 3, 'N': 4, 'S': 5, 'Si': 6, 'Na': 7, 'F': 8, 'P': 9}
```
para.rdから読み込まれる。
キーを元素記号(str)、値を原子のタイプ(int)とした辞書
例）もしatom_symbol_to_type={"Fe":1, "Cu":2}ならば
sdat.atom_symbol_to_type["Fe"]は1となる

## atom_type_to_symbol
```python
>>> sdat.atom_type_to_symbol
{1: 'C', 2: 'H', 3: 'O', 4: 'N', 5: 'S', 6: 'Si', 7: 'Na', 8: 'F', 9: 'P'}
```
para.rdから読み込まれる。
キーを原子のタイプ(int)、値を元素記号(str)とした辞書


## atom_type_to_mass
```python
>>> sdat.atom_type_to_mass
{1: 12.0107, 2: 1.00794, 3: 15.9994, 4: 14.0067, 5: 32.065, 6: 28.0855, 7: 22.98976928, 8: 18.9984032, 9: 30.973762}
```
キーを原子タイプ(int)、値を原子量(float)とした辞書。

## bondorder_connect_list
dumpbondで読み込んだそのままの隣接リスト
```python
sdat.bondorder_connect_list[atom_idx]
[  703  4910 10212 12492 20328 25777 22654]
```
## bondorder_list
dumpbondで読み込んだボンドオーダーが入っているリスト
```python
>>>sdat.bondorder_list[atom_idx]
[[0.004 0.    0.    0.004]
 [0.002 0.    0.    0.002]
 [0.593 0.    0.    0.593]
 [0.729 0.003 0.    0.733]
 [0.91  0.402 0.    1.312]
 [0.017 0.    0.    0.017]
 [0.122 0.    0.    0.122]]
```
左からσ、π、ππ、結合の合計
bondorder_connect_list[atom_idx]の行とbondorder_list[atom_idx]の列が対応している

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
## contidion関数
原子の条件を設定する関数。関数にコンディション関数を渡すことで限られた条件を満たす原子のみを考慮することができる。
作り方
第一引数をsdatとし、返す値はpd.seriesまたはnp.arrayのboolの配列を返すようにする。Trueのものは条件を満たしていて、Falseのものは条件を満たしていない。例えば、原子のz座標が15以下のもののみで結合の数を数えたいときは、
```python
>>> def condition(sdat):
...     return sdat['z'] <= 15
>>> print(sdat.count_bonds(0.5, bond_type='dumpbond', condition=condition))
{'H-H': 0, 'H-O': 0, 'H-N': 0, 'H-Si': 0, 'O-O': 0, 'O-N': 0, 'O-Si': 0, 'N-N': 0, 'N-Si': 9422, 'Si-Si': 167}
```
## 結合の条件
cut_off:カットオフ
bond_type:ボンドの作り方。'dumppos'とすると、それぞれの原子からcut_off以下の原子を結合しているとみなす。'dumpbond'とすると、bondorderがcut_off以上の原子を結合しているとみなす

## get_total_atoms()
```python
>>> sdat.get_total_atoms()
37359
```
 全原子数を返す

## wrap_prticles()
原子の座標をセルの範囲内に収める。


## get_connect_list(cut_off, bond_type)
第一引数(cut_off) : カットオフ距離
第二引数(bond_type) : ボンドの作り方。'dumppos'とすると、それぞれの原子からcut_off以下の原子を結合しているとみなす。'dumpbond'とすると、bondorderがcut_off以上の原子を結合しているとみなす
```python
>>> connect_list = sdat.get_connect_list(0.5, bond_type='dumpbond')
>>> connect_list[3]
[1452, 6471]
```
この場合atomのindexが3の原子はatomのindexが1452, 6471の２つと結合している

## import_para(input_filename)
para.rdを読み込む。
dumppos, dumpbond, input.rd, config.rd, xyzファイルよりも先に読み込んでください。

## import_input(input_filename)
input.rdを読み込む。

## import_config(input_filename)
config.rdを読み込む。

## import_xyz(input_filename)
xyzファイルを読み込む。

## import_dumppos(input_filename)
dumpposファイルを読み込む。

## import_dumpbond(input_filename)
dumpbondファイルを読み込む。


## export_input(output_filename)
input.rdを出力する。
第一引数(output_filename) : 保存するファイル名

## export_dumppos(output_filename, time_step, out_columns)
dumpposを出力する。
第一引数(output_filename) : 保存するファイル名
time_step : タイムステップ、デフォルトは0
out_columns : 出力する列を指定

## export_xyz(output_filename, out_columns, structure_name)
xyzファイルを出力する。
第一引数(output_filename) : 保存するファイル名
out_columns : 出力する列を指定
structure_name : 構造の名前

## replicate_particles([x, y, z])
第一引数(replicate_particles) : 等倍したいx,y,z方向へのリスト．

周期境界方向にセルを等倍し，粒子をコピーする．
```python
#例) x,y,z方向にそれぞれセルを2倍にしてコピーする
sdat.replicate_particles([2,2,2])
```
## concat_atoms(outer_sdat)
sdatにouter_sdatを結合する
```python
#例) sdatにouter_sdatを結合する
sdat.concat_atoms(outer_sdat)
```

## delete_atoms(condition, reindex)
条件に当てはまる原子を削除する

condition : function or np.array
    削除したい原子の条件を指定する関数またはnp.array
reindex : bool
    reindex == Trueの時は原子のid(idx)は新しく割り振られる
    reindex == Falseの時は原子のid(idx)は新しく割り振られず、削除前のものと変わらない

## count_bonds(cut_off, bond_type)
結合数をカウントする関数。
cut_off:カットオフ
bond_type:'dumppos' or 'dumpbond'
'dumppos'の場合はcut_off以下の原子を結合しているとみなす。
'dumpbond'の場合はcut_off以上のbond orderの結合を結合しているとみなす。
```python
>>> sdat.count_bonds(2, bond_type='dumppos')
{'H-H': 4796, 'H-O': 8659, 'H-N': 465, 'H-Si': 0, 'O-O': 0, 'O-N': 228, 'O-Si': 375, 'N-N': 60, 'N-Si': 43673, 'Si-Si': 0}
>>> sdat.count_bonds(0.5, bond_type='dumpbond')
{'H-H': 0, 'H-O': 5649, 'H-N': 350, 'H-Si': 0, 'O-O': 0, 'O-N': 169, 'O-Si': 385, 'N-N': 7, 'N-Si': 44848, 'Si-Si': 881}
```
## count_mols(self, cut_off, lower_mol_limit=1, upper_mol_limit=10)
分子数を数える関数。
結合数をカウントする関数。
cut_off:カットオフ
bond_type:'dumppos' or 'dumpbond'
'dumppos'の場合はcut_off以下の原子を結合しているとみなす。
'dumpbond'の場合はcut_off以上のbond orderの結合を結合しているとみなす。
分子内の原子の数がlower_mol_limit以上upper_mol_limit以下の分子のみがカウントされる。
```python
>>> sdat.count_mols(0.5, bond_type='dumpbond')
{'N1': 2, 'H5O2N1': 1, 'H4O1N1': 1, 'H2O1': 2395, 'H3O1': 101, 'H1': 1}
```

## get_atom_idx_from_mol(cuf_off, target_mol)
target_molに対応する原子ののindexを取り出す。target_molにはタプルを渡す。
```python
例水分子の個数を数える。
#          ↓ゼロ,C, H, O, N, S, Si
target_mol=(0, 0, 2, 1, 0, 0)
>>> sdat.get_atom_idx_from_mol(0.3, target_mol=(0, 0, 2, 1, 0, 0))[:10]
[[28313, 28315, 28314], [28325, 28327, 28326], [28334, 28336, 28335], [28340, 28342, 28341], [28343, 28345, 28344], [28349, 28351, 28350], [28355, 28357, 28356], [28358, 28360, 28359], [28361, 28363, 28362], [28364, 28366, 28365]]
```

## count_triplets(cut_off, bond_type, condition)
３体間の原子をカウントする関数
例えば、水分子一つにH-O-Hは1個含まれている。
メタン分子一つにはcombination(4,2) = 6より、H-C-Hは6個含まれている。
```python
>>> print(sdat.count_bonds(0.5, bond_type='dumpbond', condition=condition))
{'H-H': 0, 'H-O': 0, 'H-N': 0, 'H-Si': 0, 'O-O': 0, 'O-N': 0, 'O-Si': 0, 'N-N': 0, 'N-Si': 9422, 'Si-Si': 167}
>>> sdat.count_triplets(0.5, bond_type='dumpbond')
{'H-H-H': 0, 'H-H-O': 0, 'H-H-N': 0, 'H-H-Si': 0, 'H-O-H': 2784, 'H-O-O': 0, 'H-O-N': 251, 'H-O-Si': 321, 'H-N-H': 19, 'H-N-O': 84, 'H-N-N': 1, 'H-N-Si': 579, 'H-Si-H': 0, 'H-Si-O': 0, 'H-Si-N': 0, 'H-Si-Si': 0, 'O-H-H': 0, 'O-H-O': 0, 'O-H-N': 0, 'O-H-Si': 0, 'O-O-H': 0, 'O-O-O': 0, 'O-O-N': 0, 'O-O-Si': 0, 'O-N-H': 84, 'O-N-O': 14, 'O-N-N': 0, 'O-N-Si': 226, 'O-Si-H': 0, 'O-Si-O': 57, 'O-Si-N': 922, 'O-Si-Si': 67, 'N-H-H': 0, 'N-H-O': 0, 'N-H-N': 0, 'N-H-Si': 0, 'N-O-H': 251, 'N-O-O': 0, 'N-O-N': 13, 'N-O-Si': 3, 'N-N-H': 1, 'N-N-O': 0, 'N-N-N': 0, 'N-N-Si': 19, 'N-Si-H': 0, 'N-Si-O': 922, 'N-Si-N': 62049, 'N-Si-Si': 5811, 'Si-H-H': 0, 'Si-H-O': 0, 'Si-H-N': 0, 'Si-H-Si': 0, 'Si-O-H': 321, 'Si-O-O': 0, 'Si-O-N': 3, 'Si-O-Si': 34, 'Si-N-H': 579, 'Si-N-O': 226, 'Si-N-N': 19, 'Si-N-Si': 41681, 'Si-Si-H': 0, 'Si-Si-O': 67, 'Si-Si-N': 5811, 'Si-Si-Si': 62}
```

## count_atom_types(res_type, condition)
原子のtype別に原子の数をカウントする関数
restype='series'の場合はseriesを返す
restype='dict'の場合はdictを返す
```python
>>> sdat.count_atom_types(res_type='series')
N     16211
Si    12148
H      6000
O      3000
Name: type, dtype: int64
>>> sdat.count_atom_types(res_type='dict')
{'N': 16211, 'Si': 12148, 'H': 6000, 'O': 3000}
```
## density(x_max, x_min, y_max, y_min, z_max, x_min)
セル内の密度を計算する関数
x_max, x_min, y_max, y_min, z_max, x_minを指定しなければセル全体についての密度を計算する
```python
>>>sdat.density()
2.79432 
```
## get_coordination_number(cut_off=0.5,　bond_type='dumpbond', condition=None)
ある原子からみて平均して何個の原子が結合しているかを調べる
```python
>>> sdat.get_coordination_number(0.5, bond_type='dumpbond')
{'H-H': 0.0, 'H-O': 0.9415, 'H-N': 0.058333333333333334, 'H-Si': 0.0, 'O-H': 1.883, 'O-O': 0.0, 'O-N': 0.05633333333333333, 'O-Si': 0.12833333333333333, 'N-H': 0.021590278206156315, 'N-O': 0.010425020048115477, 'N-N': 0.00043180556412312625, 'N-Si': 2.7665165628277095, 'Si-H': 0.0, 'Si-O': 0.03169245966414225, 'Si-N': 3.6918011195258478, 'Si-Si': 0.07252222588080343}
```
この場合Oからみて平均して1.8個のHと結合している
この場合Hからみて平均して0.94個のOと結合している

## get_coordination_number(cut_off=0.5, target_atom_type=-1, bond_type='dumpbond', condition=None)
target_atom_typeの原子の配位数の分布を調べる

```python
>>> sdat.get_coordination_number_distribution(0.5, bond_type='dumpbond', target_atom_type=6)
{0: 0, 1: 10, 2: 149, 3: 1959, 4: 9342, 5: 686, 6: 2, 7: 0, 8: 0}
```
この場合、target_atom_typeが6の原子はSiであり、4配位となっているSiは9342個あることがわかる。

## count_terminal(cut_off=0.5, bond_type='dumpbond', condition=None)
H終端とOH終端の数をカウントする。
```python
>>> sdat.count_terminal(0.5, bond_type='dumpbond')
{'H-N': 350, 'H-Si': 0, 'Si-O-H': 313, 'N-O-H': 55}
```

## get_bond_angle(cut_off=0.5, bond_type='dumpbond', ondition=None) 
原子の角度の配列が入ったdictを返す関数
```python
>>> angles = sdat.get_bond_angle(0.5, bond_type='dumpbond')
>>> angles['H-O-H'][:10]
[102.37058492347465, 102.54134180886581, 104.93414922331057, 100.02985070789435, 106.36254007212739, 106.91552652630523, 101.24480868980716, 109.46491665905091, 105.13736475941324, 102.8663440653305]
plt.hist(angles['Si-O-Si'],bins=100)
plt.show()
```
H-O-H間の角度が約104°になっていることがわかる。
```python
>>> angles['Si-O-Si'][:10]
[124.48262528281757, 145.64000243424962, 91.99056007784137, 165.9609035258014, 147.9001538447931, 118.66324945244097, 95.66016925923351, 118.3555365076579, 128.83645357178645, 95.706759951125]
plt.hist(angles['Si-O-Si'],bins=100)
plt.show()
```
Si-O-Siの角度は95~170とばらつきがある。


## get_connected_atoms_from_atom_idx(cut_off=0.5, bond_type='dumpbond', start_atom_idx=0, res_type='list')
start_atom_idxと連結な原子すべてを探索して原子のindexを返す関数
res_type='list'のときはlistに原子のindexを入れて返す
```python
>>> sdat.get_connected_atoms_from_atom_idx(0.5, bond_type='dumpbond', start_atom_idx=30000, res_type='list')
[30000, 30001, 30002]
```
この場合この分子は水で、原子のindexが30000の原子は30001と30002と連結である。

res_type='bool'のときはnp.arrayで連結している原子はすべてTrue、連結ではない原子はすべてFalseとなるようなnp.arrayを返す
```python
>>> sdat.get_connected_atoms_from_atom_idx(0.5, bond_type='dumpbond', start_atom_idx=30000, res_type='bool')
array([False, False, False, ..., False, False, False])
```


## add_mol_q_to_atoms(cut_off=0.5, bond_type='dumpbond')
原子ごとの電荷のsumを計算し、sdat.atomsにmol_qというカラムで追加する
```python
>>> sdat.add_mol_q_to_atoms(0.5, bond_type='dumpbond')
>>> sdat.atoms
       type       x       y       z     vx     vy     vz      q  mask   mol_q
0         6  16.849  20.389  20.174  0.001 -0.004  0.000  0.527     5 -10.157
1         6  32.722  10.533   4.074  0.000 -0.000  0.000  0.521     4 -10.157
2         6  23.587  11.436  13.742  0.001 -0.004  0.000  0.509     0 -10.157
3         6   3.904  28.283   0.572  0.000  0.000 -0.000  0.116     4 -10.157
4         6  48.649  31.810  28.123  0.001 -0.004  0.000  0.501     5 -10.157
...     ...     ...     ...     ...    ...    ...    ...    ...   ...     ...
37354     2  96.081  16.863  25.520 -0.003  0.012  0.010  0.382     5  -0.005
37355     3  95.477  17.429  25.990 -0.003  0.004  0.000 -0.783     5  -0.005
37356     2  81.027   4.444  43.712  0.009 -0.032  0.019  0.266     5  -9.300
37357     2  81.443   3.474  39.168 -0.001 -0.016  0.016  0.398     5   0.025
37358     3  82.514   3.136  40.906 -0.007  0.002  0.004 -0.596     5  -9.300

[37359 rows x 10 columns]
```