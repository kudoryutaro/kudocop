import numpy as np
import pandas as pd
import numpy as np
from pathlib import Path
import os

class ExportDPSystem():
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
        atom_type = self.atoms['type'].astype('str').to_list()
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
            self.atoms内にある原子のタイプのみ出力される
        """
        dp_system_dir = Path(dp_system_dir)
        atom_type_map = []
        for atom_type in range(1, len(self.atom_type_to_symbol) + 1):
            if atom_type not in self.get_atom_type_set():
                continue
            atom_type_map.append(self.atom_type_to_symbol[atom_type] + '\n')
        with open(dp_system_dir / 'type_map.raw', 'w') as f:
            f.writelines(atom_type_map)


    def export_dp_system_set_coord(self, dp_system_set_dir:Path):
        """DeePMDのcoord.npyを作成する
        Parameter
        ---------
        dp_system_set_dir : Path
            coord.npyを作成するディレクトリのパス
        """
        dp_system_set_dir = Path(dp_system_set_dir)
        coord = self.atoms[['x', 'y', 'z']].values.reshape(1, -1)
        np.save(dp_system_set_dir / 'coord', coord)


    def export_dp_system_set_box(self, dp_system_set_dir:Path):
        """DeePMDのbox.npyを作成する
        Parameter
        ---------
        dp_system_set_dir : Path
            box.npyを作成するディレクトリのパス

        Note:
            boxのshape : (1, 9)
                [[XX XY XZ YX YY YZ ZX ZY ZZ]]
        """
        dp_system_set_dir = Path(dp_system_set_dir)
        box = np.array([[0. for _ in range(9)]])
        box[0][0] = self.cell[0]
        box[0][4] = self.cell[1]
        box[0][8] = self.cell[2]
        np.save(dp_system_set_dir / 'box', box)

    def export_dp_system_set_energy(self, dp_system_set_dir:Path):
        """DeePMDのenergy.npyを作成する
        Parameter
        ---------
        dp_system_set_dir : Path
            energy.npyを作成するディレクトリのパス
        """
        assert hasattr(self, 'atoms_calc'), 'dmol3_calcで先に第一原理計算してください'
        dp_system_set_dir = Path(dp_system_set_dir)
        dp_system_set_dir = Path(dp_system_set_dir)
        energy = self.atoms_calc.get_total_energy()
        np.save(dp_system_set_dir / 'energy', np.array([energy]))

    def export_dp_system_set_force(self, dp_system_set_dir:Path):
        """DeePMDのforce.npyを作成する
        Parameter
        ---------
        dp_system_set_dir : Path
            force.npyを作成するディレクトリのパス
        
        Note
        ----
        force : DFTの結果のforce
                shape : (Natoms, 3)

        """
        assert hasattr(self, 'atoms_calc'), 'dmol3_calcで先に第一原理計算してください'
        dp_system_set_dir = Path(dp_system_set_dir)
        force = self.atoms_calc.get_forces()
        force = force.reshape(1, -1)
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
        
        """
        assert hasattr(self, 'atoms_calc'), 'dmol3_calcで先に第一原理計算してください'
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
        




