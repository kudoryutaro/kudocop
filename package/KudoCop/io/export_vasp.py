import numpy as np
import pandas as pd
import pathlib
import subprocess
from typing import List
import os
import time

class ExportVasp():
    def __init__(self):
        pass

    def export_vasp_poscar(
            self, 
            ofn: str, 
            comment: str="", 
            scaling_factor: float=1.0,
        ) -> None:
        """vaspのincarファイルを作成する
        注意:この関数はsdat.atomsの原子をtypeごとにinplaceに並べ替えます.

        Parameters
        ----------
            ofn: str
                出力先のパス
            comment: str
                POSCARの1行目に書かれるコメント
            scaling_factor: float
                VASPを参照してください. 基本1.0でok
        """
        for dim in range(3):
            if self.cell[dim] == 0:
                print(f'warning : cell[{dim}] is 0')

            if self.cell[dim] is None:
                print(f'warning : cell[{dim}] is not defined')
                print(f'warning : cell[{dim}] has been initialized to 0')
                self.cell[dim] = 0.0

        header_line = []
        header_line.append(f"{comment.strip()}\n")
        header_line.append(f"{scaling_factor}\n")
        header_line.append(f"{self.cell[0]:.10f} 0.0000000000 0.0000000000\n")
        header_line.append(f"0.0000000000 {self.cell[1]:.10f} 0.0000000000\n")
        header_line.append(f"0.0000000000 0.0000000000 {self.cell[2]:.10f}\n")

        atom_symbol_to_atom_counter = self.count_atom_types(res_type='dict')
        atom_symbols = []
        atom_types_counter = []
        for atom_type in range(1, len(self.atom_type_to_symbol) + 1):
            atom_symbol = self.atom_type_to_symbol[atom_type]
            if atom_symbol in atom_symbol_to_atom_counter and atom_symbol_to_atom_counter[atom_symbol] != 0:
                atom_types_counter.append(atom_symbol_to_atom_counter[atom_symbol])
                atom_symbols.append(atom_symbol)

        header_line.append(" ".join(atom_symbols) + "\n")
        header_line.append(" ".join(map(str, atom_types_counter)) + "\n")
        header_line.append("Cartesian\n")

        self.atoms['symbol'] = self.atoms['type'].replace(self.atom_type_to_symbol)
        self.atoms = self.atoms.sort_values('type').reset_index(drop=True)
        with open(ofn, 'w') as ofp:
            ofp.writelines(header_line)

        self.atoms.to_csv(ofn, columns=['x', 'y', 'z'], mode='a', header=False,
                          sep=' ', float_format='%.10f', index=False)



    def export_vasp_incar(
            self, 
            ofn: str,
            config: dict,  
        ):
        """vaspのINCARファイルを作成する.
        Parameters
        ----------
            ofn: str
                出力先のファイルパス
            config: dict
                vaspの設定
                config['key'] = 'value'とすると、
                INCARでは、
                key = value
                となる
        """
        config_list = []
        for config_key, config_value in config.items():
            config_list.append(f"{config_key} = {config_value}\n")
        with open(ofn, "w") as ofp:
            ofp.writelines(config_list)

    def export_vasp_potcar(
            self, 
            ofn: str,
            potcar_root: str,  
        ):
        """vaspのINCARファイルを作成する.
        Parameters
        ----------
            ofn: str
                出力先のファイルパス
            potcar_root: str
                potcarが入っているフォルダのパス
        """
        potcar_root = pathlib.Path(potcar_root)
        make_potcar_command_list = ["cat"]

        atom_symbol_to_atom_counter = self.count_atom_types(res_type='dict')
        for atom_type in range(1, len(self.atom_type_to_symbol) + 1):
            atom_symbol = self.atom_type_to_symbol[atom_type]
            if atom_symbol in atom_symbol_to_atom_counter and atom_symbol_to_atom_counter[atom_symbol] != 0:
                potcar_path = potcar_root / atom_symbol / "POTCAR"
                make_potcar_command_list.append(f"{potcar_path.resolve()}")
        make_potcar_command_list.append(f" > {ofn}")
        make_potcar_command = " ".join(make_potcar_command_list)
        subprocess.run(make_potcar_command, shell=True)

    def export_vasp_kpoints(
            self, 
            ofn: str,
            comment: str="",
            kx: int=1,
            ky: int=1,
            kz: int=1,
        ):
        """vaspのINCARファイルを作成する.
        Parameters
        ----------
            ofn: str
                出力先のファイルパス
            kx: int
            ky: int
            kz: int
                x, y, z方向のK点
        """

        output = []
        output.append(comment.strip() + "\n")
        output.append("0\n")
        output.append("Monkhorst\n")
        output.append(f"{kx} {ky} {kz}\n")
        output.append("0 0 0\n")

        with open(ofn, "w") as ofp:
            ofp.writelines(output)



    def export_vasp_iconst(
            self, 
            ofn: str,
            config: List[str],  
        ):
        """vaspのICONSTファイルを作成する.
        Parameters
        ----------
            ofn: str
                出力先のファイルパス
            config: List[str]
                vaspの設定
                config['key'] = 'value'とすると、
                INCARでは、
                key value
                となる
        Examples
        --------
            VASPでNPT計算をする時で、セルの角度を固定したい時、
            config = ['LA 1 2 0',
                      'LA 1 3 0',
                      'LA 2 3 0'
                    ]
            とすると、ICONSTには
            ```
            LA 1 2 0
            LA 1 3 0
            LA 2 3 0
            ```
            と出力される
        """
        for line_idx in range(len(config)):
            config[line_idx] = config[line_idx].rstrip() + '\n'
        
        with open(ofn, "w") as ofp:
            ofp.writelines(config)