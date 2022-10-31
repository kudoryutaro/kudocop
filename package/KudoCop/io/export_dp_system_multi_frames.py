import numpy as np
import pandas as pd
import numpy as np
from pathlib import Path
import os

class ExportDPSystemMultiFrames():
    def __init__(self):
        pass

    def export_dp_system_atom_type(self, dp_system_dir:Path):
        """DeePMDのtype.rawを作成する
        Parameter
        ---------
        dp_system_dir : Path
            type.rawを作成するディレクトリのパス

        Note
        ----
            type.rawは0から始まる, sdat内の原子のtypeは1から始まる.
        """
        dp_system_dir = Path(dp_system_dir)
        atom_type = (self.atoms[0]['type'] - 1).astype('str').to_list()
        atom_type = map(lambda s: s + '\n', atom_type)
        with open(dp_system_dir / 'type.raw', 'w') as f:
            f.writelines(atom_type)

    def export_dp_system_atom_type_map(self, dp_system_dir:Path):
        """DeePMDのtype_map.rawを作成する
        Parameter
        ---------
        dp_system_dir : Path
            type.rawを作成するディレクトリのパス
        Note
        ----
            type.rawは0から始まる, sdat内の原子のtypeは1から始まる.
        """
        dp_system_dir = Path(dp_system_dir)
        atom_type_map = []
        for atom_type in range(1, len(self.atom_type_to_symbol) + 1):
            atom_type_map.append(self.atom_type_to_symbol[atom_type] + '\n')
        with open(dp_system_dir / 'type_map.raw', 'w') as f:
            f.writelines(atom_type_map)


    def export_dp_system_set_coord(self, dp_system_set_dir:Path):
        """DeePMDのcoord.npyを作成する
        Parameter
        ---------
        dp_system_set_dir : Path
            coord.npyを作成するディレクトリのパス
        Note
        ----
            coord.npyのshape : (Nframes, Natoms*3)
        """
        dp_system_set_dir = Path(dp_system_set_dir)
        coord = []
        for step_idx in range(len(self.step_nums)):
            coord.append(self.atoms[step_idx][['x', 'y', 'z']].values.ravel())
        coord = np.array(coord)
        np.save(dp_system_set_dir / 'coord', coord)


    def export_dp_system_set_box(self, dp_system_set_dir:Path):
        """DeePMDのbox.npyを作成する
        Parameter
        ---------
        dp_system_set_dir : Path
            box.npyを作成するディレクトリのパス

        Note
        ----
            boxのshape : (Nframes, 9)
                [[XX XY XZ YX YY YZ ZX ZY ZZ]]
        """
        dp_system_set_dir = Path(dp_system_set_dir)
        box = np.array([[0. for _ in range(9)] for _ in range(len(self.step_nums))])
        for step_idx in range(len(self.step_nums)):
            box[step_idx][0] = self.cell[0]
            box[step_idx][4] = self.cell[1]
            box[step_idx][8] = self.cell[2]
        np.save(dp_system_set_dir / 'box', box)

    def export_dp_system_set_energy(self, dp_system_set_dir:Path):
        """DeePMDのenergy.npyを作成する
        Parameter
        ---------
        dp_system_set_dir : Path
            energy.npyを作成するディレクトリのパス
        Note
        ----
            energy.npyのshape : (Nframes, )
        """
        dp_system_set_dir = Path(dp_system_set_dir)
        energy = self.potential_energy
        energy = np.array(energy)
        np.save(dp_system_set_dir / 'energy', energy)

    def export_dp_system_set_force(self, dp_system_set_dir:Path):
        """DeePMDのforce.npyを作成する
        Parameter
        ---------
        dp_system_set_dir : Path
            force.npyを作成するディレクトリのパス
        
        Note
        ----
        force : DFTの結果のforce
                shape : (Nframes, Natoms*3)

        """
        dp_system_set_dir = Path(dp_system_set_dir)
        if 'fx' not in self.atoms[0].columns:
            self.add_force_to_atoms()
        force = []
        for step_idx in range(len(self.step_nums)):
            force.append(self.atoms[step_idx][['fx', 'fy', 'fz']].values.ravel())
        np.save(dp_system_set_dir / 'force', np.array(force))

    def export_dp_system(self, dp_system_dir:Path, set_dir_name:str,
                        export_properties=['coord', 'box', 'energy', 'force']):
        """DeePMDで読み込むことの出来る形式のデータセットを作成する.
        Parameters
        ----------
            dp_system_dir : Path
                systemのディレクトリのパス
            set_dir_name : str
                system内のsetの名前
        Note
        ----
            作成されるデータセットの形式
                dp_system_dir
                    ├── dp_set_name
                    │   ├── coord.npy        # 原子の座標
                    │   ├── force.npy        # それぞれの原子にかかる力
                    │   ├── energy.npy       # それぞれのフレームのpotential energy
                    │   └── box.npy          # それぞれのフレームのセルの大きさ
                    ├── type_map.raw         # 原子のsymbolとtypeの対応関係
                    └── type.raw             # 原子のtype
        """
        assert len(self.atoms) != 0, 'len(sdats.atoms)が0です'
        dp_system_dir = Path(dp_system_dir)

        os.makedirs(dp_system_dir, exist_ok=True)
        dp_system_set_dir = dp_system_dir / set_dir_name
        os.makedirs(dp_system_set_dir, exist_ok=True)

        self.export_dp_system_atom_type(dp_system_dir=dp_system_dir)
        self.export_dp_system_atom_type_map(dp_system_dir=dp_system_dir)

        if 'coord' in export_properties:
            self.export_dp_system_set_coord(dp_system_set_dir=dp_system_set_dir)
        if 'box' in export_properties:
            self.export_dp_system_set_box(dp_system_set_dir=dp_system_set_dir)
        if 'energy' in export_properties:
            self.export_dp_system_set_energy(dp_system_set_dir=dp_system_set_dir)
        if 'force' in export_properties:
            self.export_dp_system_set_force(dp_system_set_dir=dp_system_set_dir)
        




