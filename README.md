# KudoCop
Kudo's COworking python Package

# What you can do with KudoCop
- 分子動力学シミュレーションの結果を読み込み, 書き込み (対応ファイル : [Data Formats](./docs/package/data_formats.md))
- 分子動力学シミュレーションの結果を分析
- DMol3を用いた第一原理計算のラッパー
- Deep Learning Potentialを学習させるためのデータセットを作成


# Requirements
```
numpy>=1.20.3
pandas>=1.3.4
tqdm>=4.62.3
cython=0.29.24
ase>=3.22.1
```
# Install
インストール方法
```sh
git clone git@github.com:kudoryutaro/kudocop.git
cd ./kudocop/package
python3 ./KudoCop/setup.py build_ext --inplace
```

# Add Path
```sh
#bashの場合 -> $HOME/.bashrcに下記を追記
export PYTHONPATH="$PYTHONPATH:$HOME/kudocop/package"
#cshellの場合 -> $HOME/.cshrcに下記を追記
# すでにPYTHONPATHがある場合はすでにあるPYTHONPATHに追加
setenv PYTHONPATH ${PYTHONPATH}:$HOME/kudocop/package
# PYTHONPATHが無い場合
setenv PYTHONPATH $HOME/kudocop/package
```
# How to update
mainにpushしないでください。
各自ブランチを作りpull requestをしてください。
Docstringを必ず書いてください。

# Setting of VSCode
補完が出るようになります
```
例 (nfshome15/rkudo を変えてください)
"python.analysis.extraPaths": [
    "/nfshome15/rkudo/kudocop/package"
],

```

# Documents of KudoCop
- [Tutorial of sdat](./docs/tutorial/tutorial_sdat.ipynb)
- [Data Formats](./docs/package/data_formats.md)
- [SimulationDat(sdat)](./docs/package/SimulationDat.md)
- [SimulationDats(sdats)](./docs/package/SimulationDats.md)
- [KudoCop Scripts](./docs/package/kudocop_scripts.md)
- [Structure of KudoCop](./docs/package/kudocop_structure.md)
