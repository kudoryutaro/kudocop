# Simulationdat
## 目次
1. [変数](#anchor1)
2. [関数](#anchor2)

## Import
```python
from KudoCop import SimulationDats
sdats = SimulationDats(para_file_name, dir_name, step_nums, skip_nums)
```
para_file_name :　パラメータファイルの名前、デフォルトはpara.rd
dir_name : dump.posとdump.bondファイルが入っているディレクトリのパス、
与えなかった場合はカレントディレクトリになる
step_nums : 読み込むステップ、0 ~ 50000で2000ステップごとに読み込む場合はstep_nums=range(0, 50000, 2000)
skip_num:何個つづ飛ばして読み込むか
skip_num=1の場合は全部読み込み、skip_num=100の場合は100個飛ばしで読み込む
skip_numを使うことで読み込み時間を削減できる

<a id="anchor1"></a>
# 変数

## step_nums
読み込んだステップ数が入っているリスト
```python
>>> sdats.step_nums
[0, 1200000, 2400000, 3600000, 4800000, 6000000, 7200000, 8400000, 9600000, 10800000, 12000000, 13200000]
```
## cell
 セルの最大値についての配列 [x, y, z]
```python
>>> sdat.cell
[100.0, 50.0, 110.0]
```
 ## atoms
 それぞれのステップの原子の情報の入ったDataFrameを含むリスト
 ```python
 >>>sdat.atoms[step_idx]
        type       x       y       z     vx     vy     vz      q  mask
0         5  36.205  47.491   0.093  0.000  0.000  0.000  0.585     4
1         5  86.205  47.491   0.093  0.000  0.000  0.000  0.585     4
2         5  47.128  33.330  43.737 -0.005  0.003  0.005  0.836     5
3         5  47.133  33.647  73.775  0.000  0.000 -0.000  0.881     2
4         5  97.133  33.647  73.775  0.000  0.000 -0.000  0.878     2
...     ...     ...     ...     ...    ...    ...    ...    ...   ...
34338     2  82.687  21.250  33.712 -0.004  0.001 -0.005  0.413     5
34339     3  82.107  20.832  34.307 -0.000  0.003  0.001 -0.767     5
34340     2  20.423  43.491  40.693 -0.005  0.026 -0.003  0.346     5
34341     2  20.588  43.798  42.087 -0.000 -0.004  0.005  0.335     5
34342     3  20.981  43.183  41.445  0.002  0.003 -0.001 -0.837     5

[34343 rows x 9 columns]
 ```

## atom_symbol_to_type
para.rdから読み込まれる。
キーを元素記号(str)、値を原子のタイプ(int)とした辞書
例）もしatom_symbol_to_type={"Fe":1, "Cu":2}ならば
sdat.atom_symbol_to_type["Fe"]は1となる
```python
>>> sdats.atom_symbol_to_type
{'C': 1, 'H': 2, 'O': 3, 'N': 4, 'S': 5, 'Si': 6, 'Na': 7, 'F': 8, 'P': 9}
```
## atom_type_to_symbol
para.rdから読み込まれる。
キーを原子のタイプ(int)、値を元素記号(str)とした辞書
```python
>>> sdats.atom_type_to_symbol
{1: 'C', 2: 'H', 3: 'O', 4: 'N', 5: 'S', 6: 'Si', 7: 'Na', 8: 'F', 9: 'P'}
```

## atom_type_to_mass
キーを原子タイプ(int)、値を原子量(float)とした辞書。
```python
>>> sdats.atom_type_to_mass
{1: 12.0107, 2: 1.00794, 3: 15.9994, 4: 14.0067, 5: 32.065, 6: 28.0855, 7: 22.98976928, 8: 18.9984032, 9: 30.973762}
```
## bondorder_lists
dumpbondで読み込んだボンドオーダーが入っているリスト
```python
>>>sdat.bondorder_lists[step_idx][atom_idx]
[[0.004 0.    0.    0.004]
 [0.002 0.    0.    0.002]
 [0.593 0.    0.    0.593]
 [0.729 0.003 0.    0.733]
 [0.91  0.402 0.    1.312]
 [0.017 0.    0.    0.017]
 [0.122 0.    0.    0.122]]
```
左からσ、π、ππ、結合の合計
## bondorder_connect_lists
dumpbondで読み込んだそのままの隣接リスト

```python
>>>sdat.bondorder_connect_list[step_idx][atom_idx]
[  703  4910 10212 12492 20328 25777 22654]
```
<a id="anchor2"></a>
# 関数

## contidion関数
原子の条件を設定する関数。関数にコンディション関数を渡すことで限られた条件を満たす原子のみを考慮することができる。
作り方
第一引数をsdats, 第二引数をskip_idxとし、返す値はpd.seriesまたはnp.arrayのboolの配列を返すようにする。Trueのものは条件を満たしていて、Falseのものは条件を満たしていない。例えば、原子のz座標が15以下のもののみで結合の数を数えたいときは、
```python
>>> def condition(sdats, skip_idx):
...        return sdats.atoms[skip_idx]['z'] <= 15
>>>sdats.count_bonds(0.5, bond_type='dumpbond', contidion=condition)
          H-H  H-O  H-N  H-Si  O-O  O-N  O-Si  N-N  N-Si  Si-Si
0           0    0    0     0    0    0     0    0  9422    167
1200000     0    0    0     0    0    0     0    0  9413    169
2400000     0    0    0     0    0    0     0    0  9408    168
3600000     0    0    0     0    0    0     0    0  9412    170
4800000     0    0    0     0    0    0     0    1  9427    169
6000000     0    0    0     0    0    0     0    0  9427    165
7200000     0    0    0     0    0    0     0    0  9418    169
8400000     0    0    0     0    0    0     0    0  9429    169
9600000     0    0    0     0    0    0     0    0  9416    169
10800000    0    0    0     0    0    0     0    0  9432    170
12000000    0    0    0     0    0    0     0    0  9407    168
13200000    0    0    0     0    0    0     0    1  9439    172
```
## get_total_atoms()
 全原子数を返す

## get_connect_lists(cut_off)
結合の強さの合計がcut_off以上のもので隣接リストを作り返す
```python
>>>sdats.get_connect_lists(0.3)[step_idx][atom_idx]
[10213, 12493, 20329]
```

## count_bonds(cut_off, bond_type, condition) -> pd.DataFrame
結合の強さの合計がcut_off以上のもので結合数をカウントする
```python
>>>sdats.count_bonds(0.5)
      H-H   H-O  H-N  H-Si  O-O  O-N  O-Si   N-N   N-Si  Si-Si
0       0  3969    0    37    0    0   770  3359  35465   8071
1000    0  3968    0    37    0    0   760  3359  35464   8078
2000    0  3968    0    36    0    0   768  3359  35464   8064
3000    0  3967    0    36    0    0   773  3359  35464   8097
4000    0  3969    0    36    0    0   765  3359  35465   8067
5000    0  3968    0    36    0    0   769  3359  35466   8076
6000    0  3968    0    36    0    0   760  3359  35465   8077
7000    0  3968    0    36    0    0   765  3359  35463   8096
8000    0  3969    0    36    0    0   771  3359  35464   8116
9000    0  3969    0    36    0    0   772  3359  35464   8075
```

### count_mols(cut_off, lower_mol_limit=1, upper_mol_limit=10, rename_columns=True)-> pd.DataFrame 
結合の強さの合計がcut_off以上のものを結合しているとみなし、
分子の数を数える
分子内の原子の数はlower_mol_limit以上upper_mol_limit以下となる
rename_colummsをFalseにするとカラム名はタプルのままになる
```python
>>> sdats.count_mols(0.5, bond_type='dumpbond', upper_mol_limit=5)
          N1  H1  H1O1  H2N1  H2O1  H3N1  H3O1  H3O1N1  H4N1
0          2   1     0     0  2395     0   101       0     0
1200000    2   1     0     2  1957    42    18       2     8
2400000    2   0     7     2  1444   118     9       6    11
3600000    2   0     7     2   996   244     3      12    17
4800000    2   0    11     2   678   308     0      19    23
6000000    2   1    10     3   420   413     0       5    22
7200000    2   1     9     2   268   437     0      13    16
8400000    2   0     5     0   145   509     0       6    16
9600000    2   1     4     0    75   477     0       6    21
10800000   2   0     7     1    31   435     0       7    11
12000000   2   0     2     0    25   392     0       2    12
13200000   2   1     0     1    14   356     0       3    10
```


## export_dumpposes(output_folder, out_columns=None)
dumpposを出力する。
第一引数(output_folder) : 保存するフォルダー名
out_columns : 出力する列を指定ないときは自動で決定される

## export_dp_system(dp_system_dir:Path, set_dir_name:str, export_properties=['coord', 'box', 'energy', 'force'])
NNP4laich, DeePMD-kitで読み込むことの出来るdp形式のデータセットを作成する.
dp_system_dir : Path
      systemのディレクトリのパス
set_dir_name : str
      system内のsetの名前
```
作成されるデータセットの形式
      dp_system_dir
            ├── type_map.raw         # 原子のsymbolとtypeの対応関係
            ├── type.raw             # 原子のtype
            ├── set_dir_name.0
            │   ├── coord.npy        # 原子の座標
            │   ├── force.npy        # それぞれの原子にかかる力
            │   ├── energy.npy       # それぞれのフレームのpotential energy
            │   └── box.npy          # それぞれのフレームのセルの大きさ
            ├── set_dir_name.1
            │   ├── coord.npy        # 原子の座標
            │   ├── force.npy        # それぞれの原子にかかる力
            │   ├── energy.npy       # それぞれのフレームのpotential energy
            │   └── box.npy          # それぞれのフレームのセルの大きさ
            ├....
```