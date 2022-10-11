import os
from pathlib import Path
import subprocess


class Packmol():
    def __init__():
        pass

    def packmol(self, sdat_list:list, pack_numbers_list:list, tolerance=2.0, 
                    packmol_tmp_dir='./packmol_tmp',xyz_condition=None, seed=-1):
        """packmolで原子を詰める関数
        Parameters
        ----------
            sdat_list : list of SimulationDat
                詰めるsdatのlist
            pack_numbers_list : list
                詰めるsdatの個数のリスト
            tolerance : float
                最小の原子間距離, 原子を詰める時に原子間がtolerance以上になるように詰める
            xyz_condition : list
                詰める原子の座標の条件
            packmol_tmp_dir : Path or str
                packmolを動かすときの一時的なディレクトリ
            seed : int
                シード値
                seed = -1のときはseedは時間で決定される
                
        Example
        -------
        xyz_condition = [
          # [xmin, ymin, zmin, xmax, ymax, zmax] で指定する
            [2, 2, 2, 8, 18, 28],  # sdat2の条件
            [2, 2, 2, 8, 10, 10]   # sdat3の条件
        ]
        sdat1.packmol(sdat_list=[sdat2, sdat3], pack_numbers_list=[5, 8], xyz_condition=xyz_condition)
        とすると, sdat1にsdat2が5個, sdat3が8個詰められる
        sdat2は2 <= x <= 8 かつ 2 <= y <= 18 かつ 2 <= z <= 28 の位置のみに詰められる
        sdat3は2 <= x <= 8 かつ 2 <= y <= 10 かつ 2 <= z <= 10 の位置のみに詰められる

        Note
        ----
            周期境界条件の場合は端まで詰めると計算が回らなくなるので注意.
            境界には間を開けることを推奨

            この関数を使うには、$ packmol のみでコマンドが使えるようにパスを通しておく必要がある
        """
        if xyz_condition is not None:
            assert len(xyz_condition) == len(sdat_list), 'sdat_listとxyz_conditionの長さが合いません'

        packmol_tmp_dir = Path(packmol_tmp_dir)
        os.makedirs(packmol_tmp_dir, exist_ok=True)

        # make config of packmol
        with open(packmol_tmp_dir / 'packmol_mixture_comment.inp', 'w') as f:
            f.write(f'tolerance {tolerance}\n')
            f.write(f'filetype xyz\n')
            f.write(f'seed {seed}\n')
            f.write(f'output packmol_mixture_result.xyz\n')
            f.write(f'\n')

            if self.atoms is not None and self.get_total_atoms() >= 1:
                f.write(f'structure this_sdat.xyz\n')
                f.write(f'\tnumber 1\n')
                f.write(f'\tfixed 0. 0. 0. 0. 0. 0. \n')
                f.write(f'end structure\n')
                f.write(f'\n')

            for sdat_idx in range(len(sdat_list)):
                f.write(f'structure sdat_idx_{sdat_idx}.xyz\n')
                f.write(f'\tnumber {pack_numbers_list[sdat_idx]}\n')
                if xyz_condition is not None:
                    f.write(f'\tinside box {xyz_condition[sdat_idx][0]} {xyz_condition[sdat_idx][1]} {xyz_condition[sdat_idx][2]} {xyz_condition[sdat_idx][3]} {xyz_condition[sdat_idx][4]} {xyz_condition[sdat_idx][5]}\n')
                f.write(f'end structure\n')
                f.write(f'\n')
        

                
        # export xyz file
        if self.atoms is not None and self.get_total_atoms() >= 1:
            self.export_xyz(packmol_tmp_dir / 'this_sdat.xyz', structure_name=f'this_sdat')

        for sdat_idx in range(len(sdat_list)):
            sdat_list[sdat_idx].export_xyz(packmol_tmp_dir / f'sdat_idx_{sdat_idx}.xyz',
                                            structure_name=f'sdat_idx_{sdat_idx}')

        # excecute packmol
        cmd = f"packmol < {'packmol_mixture_comment.inp'}"
        print(cmd)
        p = subprocess.Popen(cmd, cwd=packmol_tmp_dir, shell=True)
        p.wait()

        # import result of packmol
        self.import_xyz(packmol_tmp_dir / 'packmol_mixture_result.xyz')
