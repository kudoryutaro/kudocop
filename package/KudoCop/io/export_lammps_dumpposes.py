import numpy as np
import pandas as pd
from tqdm import trange

class ExportLammpsDumpposes():
    def __init__(self):
        pass

    def export_lammps_dumpposes(self, ofn: str, out_columns=None) -> None:
        """lammps形式のdumpposを出力する
        """
        if out_columns is None:
            out_columns = ['type', 'x', 'y', 'z']

        with open(ofn, 'w') as f:
            f.write('')

        for step_idx in trange(len(self.step_nums), desc='[exporting lammps dumpposes]'):
            header = []
            header.append(f'ITEM: TIMESTEP\n')
            header.append(f'{self.step_nums[step_idx]}\n')
            header.append(f'ITEM: NUMBER OF ATOMS\n')
            header.append(f'{self.get_total_atoms()}\n')
            header.append(f'ITEM: BOX BOUNDS xy xz yz pp pp pp\n')
            header.append(f'0.0000000000000000e+00 {self.cell[0]:.16e} 0.0000000000000000e+00\n')
            header.append(f'0.0000000000000000e+00 {self.cell[1]:.16e} 0.0000000000000000e+00\n')
            header.append(f'0.0000000000000000e+00 {self.cell[2]:.16e} 0.0000000000000000e+00\n')
            header.append(f'ITEM: ATOMS id {" ".join(out_columns)}\n')

            with open(ofn, 'a') as f:
                f.writelines(header)
            
            # 1-index
            self.atoms[step_idx].index += 1
            self.atoms[step_idx].to_csv(ofn, columns=out_columns, sep=' ', header=None, mode='a')
            # 0-index
            self.atoms[step_idx].index -= 1