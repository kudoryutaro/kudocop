import numpy as np
import pandas as pd
import time


class ImportSumforce():
    def __init__(self):
        pass

    def import_sumforce(self, cluster_idx=None, ifn=None) -> None:
        """
        sumforceを読み込む関数
        読み込んだらsdat.sumforceにデータが入る
        cluster_idxで指定した場合はsdat.sumforce[cluster_idx]にデータフレームが入る
        ifnで指定した場合はsdat.sumforce[ifn]にデータフレームが入る

        Parameters
        ----------
        cluster_idx : int
            粒子団のインデックス(1-index)

            input.rdに以下が書かれている場合
            ------------
            ##sumforce 2
            1 2
            3
            ------------
            sumforceを測定する粒子団は二つあり、一つ目はmaskが1と2の粒子団
            二つ目はmaskが3の粒子団
            maskが1と2の粒子団のparticle_cluster_idxが1であり、
            maskが3の粒子団のparticle_cluster_idxが2である

        Returns
        -------
        None
        """
        if self.sumforce is None:
            self.sumforce = {}
        if cluster_idx is not None:
            df_sumforce = pd.read_csv(
                f'sumforce_{cluster_idx}.dat', sep='\s+', index_col=0)
        else:
            df_sumforce = pd.read_csv(ifn, sep='\s+')

        if ifn is not None:
            self.sumforce[ifn] = df_sumforce
        else:
            self.sumforce[cluster_idx] = df_sumforce
