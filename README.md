# KudoCop
Kudo's COworking python Package

分子動力学シミュレーションの結果を読み込み, 書き込み, 分析などができるライブラリです。

# Requirements
```
numpy>=1.20.3
pandas>=1.3.4
tqdm>=4.62.3
cython=0.29.24
```
# Install
インストール方法
```sh
git clone git@github.com:kudoryutaro/kudocop.git
cd ./kudocop/package/KudoCop
python3 setup.py build_ext --inplace
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
# 更新方法
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
- [SimulationDat(sdat)](./docs/package/SimulationDat.md)
- [SimulationDats(sdats)](./docs/package/SimulationDats.md)
- [Structure of KudoCop](./docs/package/kudocop_structure.md)
