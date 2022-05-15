# Simulationdat
## 目次
1. [変数](#anchor1)
2. [関数](#anchor2)

## 前提条件
    from KudoCop.SimulationDats import SimulationDats
    sdats = SimulationDats(para_file_name, dir_name, step_nums,
             import_dumpposes_flag, import_dumpbonds_flag)
para_file_name :　パラメータファイルの名前、デフォルトはpara.rd
dir_name : dump.posとdump.bondファイルが入っているディレクトリのパス、
与えなかった場合はカレントディレクトリになる
step_nums : 読み込むステップ、0~50000で2000ステップごとに読み込む場合はstep_nums=range(0, 50000, 2000)
与えなかった場合はディレクトリ内のすべてのdump.pos, bump.bondが読み込まれる
import_dumpposes_flag : dumpposを読み込むか、デフォルトはTrue
import_dumpbondes_flag : dumpbondを読み込むか、デフォルトはTrue

<a id="anchor1"></a>
# 変数

## step_nums
読み込んだステップ数が入っているリスト

## cell
 セルの最大値についての配列 [x, y, z]

 ## atoms
 それぞれのステップの原子の情報の入ったDataFrameを含むリスト
 ```
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

## atom_type_to_symbol
para.rdから読み込まれる。
キーを原子のタイプ(int)、値を元素記号(str)とした辞書

## atom_type_set
input.rd, dumppos, xyzから読み込まれる。
系内のみのatomのタイプ(int)が入ったset。

## atom_type_to_mass
キーを原子タイプ(int)、値を原子量(float)とした辞書。

## bondorder_lists
dumpbondで読み込んだボンドオーダーが入っているリスト
```
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

```
>>>sdat.bondorder_connect_list[step_idx][atom_idx]
[  703  4910 10212 12492 20328 25777 22654]
```
<a id="anchor2"></a>
# 関数
## get_total_atoms()
 全原子数を返す

## get_connect_lists(cut_off)
結合の強さの合計がcut_off以上のもので隣接リストを作り返す
```
>>>sdats.get_connect_lists(0.3)[step_idx][atom_idx]
[10213, 12493, 20329]
```

## count_bonds(cut_off) -> pd.DataFrame
結合の強さの合計がcut_off以上のもので結合数をカウントする
```
>>>sdats.count_bonds(0.3)
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
```
>>> sdats.count_mols(0.3)
       N1  N2  H2O1  H3O1  H5O2
0       4  15  1246     1   0.0
1000    4  15  1256     1   0.0
2000    4  15  1248     1   0.0
3000    4  15  1243     1   0.0
4000    4  15  1251     1   0.0
...    ..  ..   ...   ...   ...
95000   4  15  1181     2   0.0
96000   4  15  1185     2   0.0
97000   4  15  1181     2   0.0
98000   4  15  1183     2   0.0
99000   4  15  1178     2   0.0
```


## export_dumpposes(self, output_folder, out_columns=None)
dumpposを出力する。
第一引数(output_folder) : 保存するフォルダー名
out_columns : 出力する列を指定ないときは自動で決定される