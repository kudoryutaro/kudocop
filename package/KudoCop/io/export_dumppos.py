import numpy as np
import pandas as pd


class ExportDumppos():
    def __init__(self):
        pass

    def export_dumppos(self, ofn: str, time_step=None, out_columns=None) -> None:
        if out_columns is None:
            out_columns = ['type', 'mask', 'x', 'y', 'z']

        if 'mask' not in self.atoms:
            print('warning : mask is not defined')
            print('warning : mask has been initialized to 0')
            self.atoms['mask'] = np.zeros(self.get_total_atoms(),
                                          dtype=int)
        if time_step is None:
            print('warning : time_step is not defined')
            print('warning : time_step has been initialized to 0')
            time_step = 0

        for dim in range(3):
            if self.cell[dim] == 0:
                print(f'warning : cell[{dim}] is 0')

            if self.cell[dim] is None:
                print(f'warning : cell[{dim}] is not defined')
                print(f'warning : cell[{dim}] has been initialized to 0')
                self.cell[dim] = 0.0

        header_line = []
        header_line.append("ITEM: TIMESTEP\n")
        header_line.append(f"{time_step}\n")
        header_line.append("ITEM: NUMBER OF ATOMS\n")
        header_line.append(f"{self.get_total_atoms()}\n")
        header_line.append("ITEM: BOX BOUNDS pp pp pp\n")
        for dim in range(3):
            header_line.append(f"{0.0} {self.cell[dim]}\n")

        header_line.append(
            " ".join(["ITEM: ATOMS"] + ['id'] + out_columns) + '\n')

        with open(ofn, 'w') as ofp:
            ofp.writelines(header_line)

        # 1-indexed
        self.atoms.index = self.atoms.index + 1
        self.atoms.to_csv(ofn, columns=out_columns, mode='a', header=False,
                          sep='\t', float_format='%.6f')
        # 0-indexed
        self.atoms.index = self.atoms.index - 1
