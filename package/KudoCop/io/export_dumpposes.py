import os
import numpy as np
import pandas as pd
from tqdm import tqdm


class ExportDumpposes():
    def __init__():
        pass

    def export_dumpposes(self, output_folder: str, out_columns=None) -> None:
        try:
            os.makedirs(output_folder)
        except FileExistsError:
            pass

        if out_columns is None:
            out_columns = ['type', 'mask', 'x', 'y', 'z']

        if 'mask' not in self.atoms[0]:
            print('warning : mask is not defined')
            print('warning : mask has been initialized to 0')
            for step_idx, step_num in enumerate(self.step_nums):
                self.atoms[step_idx]['mask'] = np.zeros(self.get_total_atoms(),
                                                        dtype=int)

        for dim in range(3):
            if self.cell[dim] == 0:
                print(f'warning : cell[{dim}] is 0')
            if self.cell[dim] is None:
                print(f'warning : cell[{dim}] is not defined')
                print(f'warning : cell[{dim}] has been initialized to 0')
                self.cell[dim] = 0.0

        header_line = [
            "ITEM: TIMESTEP\n",
            "\n",
            "ITEM: NUMBER OF ATOMS\n",
            f"{self.get_total_atoms()}\n",
            "ITEM: BOX BOUNDS pp pp pp\n",
            f"{0.0} {self.cell[0]}\n",
            f"{0.0} {self.cell[1]}\n",
            f"{0.0} {self.cell[2]}\n",
            " ".join(["ITEM: ATOMS"] + ['id'] + out_columns) + "\n"
        ]

        for step_idx, step_num in enumerate(tqdm(self.step_nums, desc='[exporting dumppos]')):
            ofn = f'{output_folder}/dump.pos.{step_num}'
            header_line[1] = f"{step_num}\n"
            with open(ofn, 'w') as ofp:
                ofp.writelines(header_line)

            # 1-indexed
            self.atoms[step_idx].index = self.atoms[step_idx].index + 1
            self.atoms[step_idx].to_csv(ofn, columns=out_columns, mode='a', header=False,
                                        sep=' ', float_format='%.6f')
            # 0-indexed
            self.atoms[step_idx].index = self.atoms[step_idx].index - 1
