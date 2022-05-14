# Simulationdat
## 目次
1. [変数](#anchor1)
2. [関数](#anchor2)

## 前提条件
    from KudoCop.SimulationDats import SimulationDats
    sdats = SimulationDats(, para_file_name, dir_name=None, step_nums=None)
para_file_name :　パラメータファイルの名前
dir_name : dump.posとdump.bondファイルが入っているディレクトリのパス、
与えなかった場合はカレントディレクトリになる
step_nums : 読み込むステップをリストで渡す
与えなかった場合はディレクトリ内のすべてのdump.pos, bump.bondが読み込まれる
<a id="anchor1"></a>
# 変数

## step_nums
読み込んだステップ数が入っているリスト

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



<a id="anchor2"></a>
# 関数

## get_total_atoms()
 全原子数を返す

## export_dumppos(output_filename, time_step, out_columns)
dumpposを出力する。
第一引数(output_filename) : 保存するファイル名
time_step : タイムステップ、デフォルトは0
out_columns : 出力する列を指定