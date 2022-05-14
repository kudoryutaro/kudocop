# KudoCop

# KudoCopとMolCopの違い

- I/OがKudoCopの方が8倍程度速い
- よりシンプルに書ける
- dumppos, dumpbondなどをステップごとの操作を一度にできる



# Install & Update
インストール方法
```sh
git clone git@github.com:kudoryutaro/kudocop.git
```
更新方法
```sh
git push origin main
```

# How to use KudoCop
```sh
#bashの場合 -> $HOME/.bashrcに下記を追記
export PYTHONPATH="$PYTHONPATH:$HOME/kudocop/package"
#cshellの場合 -> $HOME/.cshrcに下記を追記
setenv PYTHONPATH "$HOME/kudocop/package/:$PYTHONPATH"
```


# Data structure of KudoCop
- [SimulationDat(sdat)](./docs/package/SimulationDat.md)
- [SimulationDats(sdats)](./docs/package/SimulationDats.md)

