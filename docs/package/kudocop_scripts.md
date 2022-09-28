# KudoCop scriptとは
ターミナル上で簡単に分析できるスクリプト集

# Add Path and Authorize
cshrcでパスを追加
```sh
setenv PATH :$HOME/kudocop/package/KudoCop_scripts:$PATH
```
pythonファイルを実行できるようにする
```sh
chmod u+x $HOME/kudocop/package/KudoCop_scripts/*
```

# How to use
例えば、分子数をカウントするときはcount_mols.pyを使う。
python3 count_mols.pyとせずに直接count_mols.pyだけで実行できる。
```sh
count_mols.py -s 500 | less
```
```sh
bond_type: dumpbond, cut_off: 0.5
condition: None <= x <= None
condition: None <= y <= None
condition: None <= z <= None
          N1  H1  H1O1  H2N1  H2O1  H3N1  H3O1  H3O1N1  H4N1  H4O1N1  H4O1N2  \
0          2   1     0     0  2395     0   101       0     0       1       0
2000000    2   0     5     0  1631    78     4       7    12       5       0
4000000    2   2    17     2   905   278     3       9    17      22       0
6000000    2   1    10     3   420   413     0       5    22      20       0
8000000    2   1     6     1   177   458     0       5    24      19       0
10000000   2   0     9     1    55   473     0       6     5       9       1
12000000   2   0     2     0    25   392     0       2    12       6       0

          H4O2N1  H4O2N1Si1  H4O4Si1  H4O4N1Si1  H5O1N2  H5O2N1  H5O2N2Si1  \
0              0          0        0          0       0       1          0
2000000        1          0        0          0       0       0          0
4000000        3          1        1          0       2       0          0
6000000        6          0        0          1       3       0          0
8000000        0          0        0          0       1       0          1
10000000       3          0        0          0       0       1          0
12000000       1          0        0          0       0       0          0

          H5O3N1  H5O3N1Si1  H6N3Si1  H6O1N2  H6O2N1  H6O2N2  H7O1N2
0              0          0        0       0       0       0       0
2000000        0          2        0       0       1       1       0
4000000        1          2        0       5       1       0       0
6000000        0          1        0       0       0       0       1
8000000        0          0        0       2       0       1       1
10000000       0          0        1       0       0       2       0
12000000       0          0        0       0       0       0       0
~
(END)
```

コマンドの使い方を忘れたときは--helpオプションで使い方を見ることが出来る
```sh
count_mols.py --help
```
```sh
usage: count_mols.py [-h] [-s SKIP_NUM] [-f FILE_NAME] [-c CUT_OFF] [-b {dumppos,dumpbond}]
                     [-p PARA_FILE_NAME] [--x_min X_MIN] [--y_min Y_MIN] [--z_min Z_MIN] [--x_max X_MAX]
                     [--y_max Y_MAX] [--z_max Z_MAX]

分子数をカウントする.

optional arguments:
  -h, --help            show this help message and exit
  -s SKIP_NUM, --skip_num SKIP_NUM
                        何個おきにファイルを読み込むか
  -f FILE_NAME, --file_name FILE_NAME
                        読み込むファイル名.指定しない場合はdump.pos.0から全て読み込む
  -c CUT_OFF, --cut_off CUT_OFF
                        カットオフ
  -b {dumppos,dumpbond}, --bond_type {dumppos,dumpbond}
                        結合の種類, dumppos:原子間の距離から結合を作る, dumpbond:bond orderから結合を作る
  -p PARA_FILE_NAME, --para_file_name PARA_FILE_NAME
                        para.rdのファイル名
  --x_min X_MIN         min x
  --y_min Y_MIN         min y
  --z_min Z_MIN         min z
  --x_max X_MAX         max x
  --y_max Y_MAX         max y
  --z_max Z_MAX         max z
```



# Scripts

## count_atom_types
原子数をtype別にカウントする.
```sh
usage: count_atom_types.py [-h] [-s SKIP_NUM] [-f FILE_NAME] [-p PARA_FILE_NAME] [--x_min X_MIN]
                           [--y_min Y_MIN] [--z_min Z_MIN] [--x_max X_MAX] [--y_max Y_MAX] [--z_max Z_MAX]

原子数をtype別にカウントする.

optional arguments:
  -h, --help            show this help message and exit
  -s SKIP_NUM, --skip_num SKIP_NUM
                        何個おきにファイルを読み込むか
  -f FILE_NAME, --file_name FILE_NAME
                        読み込むファイル名.指定しない場合はdump.pos.0から全て読み込む
  -p PARA_FILE_NAME, --para_file_name PARA_FILE_NAME
                        para.rdのファイル名
  --x_min X_MIN         min x
  --y_min Y_MIN         min y
  --z_min Z_MIN         min z
  --x_max X_MAX         max x
  --y_max Y_MAX         max y
  --z_max Z_MAX         max z
```
```sh
count_atom_types.py -f input.rd
```
```sh
condition: None <= x <= None
condition: None <= y <= None
condition: None <= z <= None
N     16208
Si    12148
H      5992
O      2997
Name: type, dtype: int64
```

## convert_file_type.py
input.rd, xyz, car file, dumpposファイルを相互に変換することが出来る
ファイルの種類はファイルの名前から自動で推測される。
```sh
usage: convert_file_type.py [-h] [-f FILE_NAME] [-o OUTPUT_FILE_NAME] [-p PARA_FILE_NAME]

ファイルを変換する

optional arguments:
  -h, --help            show this help message and exit
  -f FILE_NAME, --file_name FILE_NAME
                        変換前のファイル名
  -o OUTPUT_FILE_NAME, --output_file_name OUTPUT_FILE_NAME
                        変換後のファイル名
  -p PARA_FILE_NAME, --para_file_name PARA_FILE_NAME
                        para.rdのファイル名
```

## count_mols
分子数をカウントする

```sh
usage: count_mols.py [-h] [-s SKIP_NUM] [-f FILE_NAME] [-c CUT_OFF] [-b {dumppos,dumpbond}]
                     [-p PARA_FILE_NAME] [--x_min X_MIN] [--y_min Y_MIN] [--z_min Z_MIN] [--x_max X_MAX]
                     [--y_max Y_MAX] [--z_max Z_MAX]

分子数をカウントする.

optional arguments:
  -h, --help            show this help message and exit
  -s SKIP_NUM, --skip_num SKIP_NUM
                        何個おきにファイルを読み込むか
  -f FILE_NAME, --file_name FILE_NAME
                        読み込むファイル名.指定しない場合はdump.pos.0から全て読み込む
  -c CUT_OFF, --cut_off CUT_OFF
                        カットオフ
  -b {dumppos,dumpbond}, --bond_type {dumppos,dumpbond}
                        結合の種類, dumppos:原子間の距離から結合を作る, dumpbond:bond orderから結合を作る
  -p PARA_FILE_NAME, --para_file_name PARA_FILE_NAME
                        para.rdのファイル名
  --x_min X_MIN         min x
  --y_min Y_MIN         min y
  --z_min Z_MIN         min z
  --x_max X_MAX         max x
  --y_max Y_MAX         max y
  --z_max Z_MAX         max z
```

```sh
count_mols.py -s 3000
```
```sh
bond_type: dumpbond, cut_off: 0.5
condition: None <= x <= None
condition: None <= y <= None
condition: None <= z <= None
          N1  H1  H1O1  H2O1  H3N1  H3O1  H3O1N1  H4N1  H4O1N1  H4O2N1  H5O2N1
0          2   1     0  2395     0   101       0     0       1       0       1
12000000   2   0     2    25   392     0       2    12       6       1       0
```

## count_bonds
結合数をカウントする

```sh
usage: count_bonds.py [-h] [-s SKIP_NUM] [-f FILE_NAME] [-c CUT_OFF] [-b {dumppos,dumpbond}]
                      [-p PARA_FILE_NAME] [--x_min X_MIN] [--y_min Y_MIN] [--z_min Z_MIN] [--x_max X_MAX]
                      [--y_max Y_MAX] [--z_max Z_MAX]

結合数をカウントする.

optional arguments:
  -h, --help            show this help message and exit
  -s SKIP_NUM, --skip_num SKIP_NUM
                        何個おきにファイルを読み込むか
  -f FILE_NAME, --file_name FILE_NAME
                        読み込むファイル名.指定しない場合はdump.pos.0から全て読み込む
  -c CUT_OFF, --cut_off CUT_OFF
                        カットオフ
  -b {dumppos,dumpbond}, --bond_type {dumppos,dumpbond}
                        結合の種類, dumppos:原子間の距離から結合を作る, dumpbond:bond orderから結合を作る
  -p PARA_FILE_NAME, --para_file_name PARA_FILE_NAME
                        para.rdのファイル名
  --x_min X_MIN         min x
  --y_min Y_MIN         min y
  --z_min Z_MIN         min z
  --x_max X_MAX         max x
  --y_max Y_MAX         max y
  --z_max Z_MAX         max z
```

```sh
count_bonds.py -s 3000
```
```sh
bond_type: dumpbond, cut_off: 0.5
condition: None <= x <= None
condition: None <= y <= None
condition: None <= z <= None
          H-H   H-O   H-N  H-Si  O-O  O-N  O-Si  N-N   N-Si  Si-Si
0           0  5649   350     0    0  169   385    7  44848    881
12000000    0   349  5652    19    0  173  5551    2  40407    724
```


## count_triplets
3体間(X-Y-Z)の個数をカウントする
```sh
usage: count_triplets.py [-h] [-s SKIP_NUM] [-f FILE_NAME] [-c CUT_OFF] [-b {dumppos,dumpbond}]
                         [-p PARA_FILE_NAME] [--x_min X_MIN] [--y_min Y_MIN] [--z_min Z_MIN]
                         [--x_max X_MAX] [--y_max Y_MAX] [--z_max Z_MAX]

X-Y-Zの数をカウントする.

optional arguments:
  -h, --help            show this help message and exit
  -s SKIP_NUM, --skip_num SKIP_NUM
                        何個おきにファイルを読み込むか
  -f FILE_NAME, --file_name FILE_NAME
                        読み込むファイル名.指定しない場合はdump.pos.0から全て読み込む
  -c CUT_OFF, --cut_off CUT_OFF
                        カットオフ
  -b {dumppos,dumpbond}, --bond_type {dumppos,dumpbond}
                        結合の種類, dumppos:原子間の距離から結合を作る, dumpbond:bond orderから結合を作る
  -p PARA_FILE_NAME, --para_file_name PARA_FILE_NAME
                        para.rdのファイル名
  --x_min X_MIN         min x
  --y_min Y_MIN         min y
  --z_min Z_MIN         min z
  --x_max X_MAX         max x
  --y_max Y_MAX         max y
  --z_max Z_MAX         max z
```
```sh
count_triplets.py -s 3000
```
```sh
bond_type: dumpbond, cut_off: 0.5
condition: None <= x <= None
condition: None <= y <= None
condition: None <= z <= None
          H-H-H  H-H-O  H-H-N  H-H-Si  H-O-H  H-O-O  H-O-N  H-O-Si  H-N-H  \
0             0      0      0       0   2784      0    251     321     19
12000000      0      0      0       0     26      0     54     260   3031

          H-N-O  H-N-N  H-N-Si  H-Si-H  H-Si-O  H-Si-N  H-Si-Si  O-H-H  O-H-O  \
0            84      1     579       0       0       0        0      0      0
12000000    259      0    5818       0       5      52        1      0      0

          O-H-N  O-H-Si  O-O-H  O-O-O  O-O-N  O-O-Si  O-N-H  O-N-O  O-N-N  \
0             0       0      0      0      0       0     84     14      0
12000000      1       3      0      0      0       0    259      7      0

          O-N-Si  O-Si-H  O-Si-O  O-Si-N  O-Si-Si  N-H-H  N-H-O  N-H-N  \
0            226       0      57     922       67      0      0      0
12000000     159       5    5042    6666      115      0      1      0

          N-H-Si  N-O-H  N-O-O  N-O-N  N-O-Si  N-N-H  N-N-O  N-N-N  N-N-Si  \
0              0    251      0     13       3      1      0      0      19
12000000      16     54      0     18     135      0      0      0       8

          N-Si-H  N-Si-O  N-Si-N  N-Si-Si  Si-H-H  Si-H-O  Si-H-N  Si-H-Si  \
0              0     922   62049     5811       0       0       0        0
12000000      52    6666   53600     4783       0       3      16        0

          Si-O-H  Si-O-O  Si-O-N  Si-O-Si  Si-N-H  Si-N-O  Si-N-N  Si-N-Si  \
0            321       0       3       34     579     226      19    41681
12000000     260       0     135     2695    5818     159       8    35387

          Si-Si-H  Si-Si-O  Si-Si-N  Si-Si-Si
0               0       67     5811        62
12000000        1      115     4783        42
```

## count_terminal
H終端またはOH終端の数をカウントする
```sh
usage: count_terminal.py [-h] [-s SKIP_NUM] [-f FILE_NAME] [-c CUT_OFF] [-b {dumppos,dumpbond}]
                         [-p PARA_FILE_NAME] [--x_min X_MIN] [--y_min Y_MIN] [--z_min Z_MIN]
                         [--x_max X_MAX] [--y_max Y_MAX] [--z_max Z_MAX]

H終端とOH終端の数をカウントする.

optional arguments:
  -h, --help            show this help message and exit
  -s SKIP_NUM, --skip_num SKIP_NUM
                        何個おきにファイルを読み込むか
  -f FILE_NAME, --file_name FILE_NAME
                        読み込むファイル名.指定しない場合はdump.pos.0から全て読み込む
  -c CUT_OFF, --cut_off CUT_OFF
                        カットオフ
  -b {dumppos,dumpbond}, --bond_type {dumppos,dumpbond}
                        結合の種類, dumppos:原子間の距離から結合を作る, dumpbond:bond orderから結合を作る
  -p PARA_FILE_NAME, --para_file_name PARA_FILE_NAME
                        para.rdのファイル名
  --x_min X_MIN         min x
  --y_min Y_MIN         min y
  --z_min Z_MIN         min z
  --x_max X_MAX         max x
  --y_max Y_MAX         max y
  --z_max Z_MAX         max z
```
```sh
count_terminal.py -s 3000
```
```sh
bond_type: dumpbond, cut_off: 0.5
condition: None <= x <= None
condition: None <= y <= None
condition: None <= z <= None
           H-N  H-Si  Si-O-H  N-O-H
0          350     0     313     55
12000000  5652    19     255     23
```

## coor
対象となる原子の配位数をカウントする
```sh
usage: coor.py [-h] [-s SKIP_NUM] [-f FILE_NAME] [-t TARGET_ATOM_TYPE] [-c CUT_OFF]
               [-b {dumppos,dumpbond}] [-p PARA_FILE_NAME] [--x_min X_MIN] [--y_min Y_MIN] [--z_min Z_MIN]
               [--x_max X_MAX] [--y_max Y_MAX] [--z_max Z_MAX]

target_atom_typeで指定した原子の配位数をカウントする.

optional arguments:
  -h, --help            show this help message and exit
  -s SKIP_NUM, --skip_num SKIP_NUM
                        何個おきにファイルを読み込むか
  -f FILE_NAME, --file_name FILE_NAME
                        読み込むファイル名.指定しない場合はdump.pos.0から全て読み込む
  -t TARGET_ATOM_TYPE, --target_atom_type TARGET_ATOM_TYPE
                        対象となる原子のtype
  -c CUT_OFF, --cut_off CUT_OFF
                        カットオフ
  -b {dumppos,dumpbond}, --bond_type {dumppos,dumpbond}
                        結合の種類, dumppos:原子間の距離から結合を作る, dumpbond:bond orderから結合を作る
  -p PARA_FILE_NAME, --para_file_name PARA_FILE_NAME
                        para.rdのファイル名
  --x_min X_MIN         min x
  --y_min Y_MIN         min y
  --z_min Z_MIN         min z
  --x_max X_MAX         max x
  --y_max Y_MAX         max y
  --z_max Z_MAX         max z
```
```sh
coor.py -s 3000 -t 6
```
target_atom_typeが6の原子はSiで、4配位となっているSiが多い
```sh
bond_type: dumpbond, cut_off: 0.5
condition: None <= x <= None
condition: None <= y <= None
condition: None <= z <= None
          0   1    2     3     4    5   6   7   8   9   10  11  12
0          0  10  149  1959  9342  686   2   0   0   0   0   0   0
4000000    0  10  151  1650  9615  720   2   0   0   0   0   0   0
8000000    0  10  148  1587  9698  702   3   0   0   0   0   0   0
12000000   0  10  148  1578  9679  729   4   0   0   0   0   0   0
```

## density
密度を計算する
```sh
usage: density.py [-h] [-f FILE_NAME] [-p PARA_FILE_NAME] [--x_min X_MIN] [--y_min Y_MIN] [--z_min Z_MIN]
                  [--x_max X_MAX] [--y_max Y_MAX] [--z_max Z_MAX]

密度を計算する

optional arguments:
  -h, --help            show this help message and exit
  -f FILE_NAME, --file_name FILE_NAME
                        読み込むファイル名.指定しない場合はdump.pos.0から全て読み込む
  -p PARA_FILE_NAME, --para_file_name PARA_FILE_NAME
                        para.rdのファイル名
  --x_min X_MIN         min x
  --y_min Y_MIN         min y
  --z_min Z_MIN         min z
  --x_max X_MAX         max x
  --y_max Y_MAX         max y
  --z_max Z_MAX         max z
```
```sh
 density.py -f input.rd --x_min 10 --x_max 50 --y_min 10 --y_max 20
```
```sh
condition: 10.0 <= x <= 50.0
condition: 10.0 <= y <= 20.0
condition: 0 <= z <= 110.0
density (g/cm^3): 1.8906165604139749
```
