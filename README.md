# KudoCop
Kudo's COworking python Package

分子動力学シミュレーションの結果を読み込み書き込み分析ができるライブラリです。

# Install
インストール方法
```sh
git clone git@github.com:kudoryutaro/kudocop.git
cd ./KudoCop/package/KudoCop
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

# Data structure of KudoCop
- [SimulationDat(sdat)](./docs/package/SimulationDat.md)
- [SimulationDats(sdats)](./docs/package/SimulationDats.md)

