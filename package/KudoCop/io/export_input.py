import numpy as np
import pandas as pd


class ExportInput():
    def __init__(self):
        pass

    def export_input(self, ofn: str) -> None:
        for dim in range(3):
            if self.cell[dim] == 0:
                print(f'warning : cell[{dim}] is 0')

            if self.cell[dim] is None:
                print(f'warning : cell[{dim}] is not defined')
                print(f'warning : cell[{dim}] has been initialized to 0')
                self.cell[dim] = 0.0

        # header_line
        header_line = []
        header_line.append(f"#cellx {0.0}  {self.cell[0]}\n")
        header_line.append(f"#celly {0.0}  {self.cell[1]}\n")
        header_line.append(f"#cellz {0.0}  {self.cell[2]}\n")

        header_line.append("\n")

        atom_type_set = set(self.atoms['type'])

        header_line.append(f"#masses {len(atom_type_set)}\n")
        for atom_type in atom_type_set:
            header_line.append(
                f"{atom_type} {self.atom_type_to_mass[atom_type]}\n")

        header_line.append("\n")

        for fix in self.fix_info:
            header_line.append(fix + '\n')

        for move in self.move_info:
            header_line.append(move + '\n')

        for press in self.press_info:
            header_line.append(press + '\n')

        for sumforce in self.sumforce_info:
            header_line.append(sumforce + '\n')

        for thermofree in self.thermofree_info:
            header_line.append(thermofree + '\n')

        for wall in self.wall_info:
            header_line.append(wall + '\n')

        header_line.append("")

        header_line.append(f"#atoms {self.get_total_atoms()}\n")

        with open(ofn, 'w') as ofs:
            ofs.writelines(header_line)

        # id, type, mask, pos, velo
        out_columns = ['type', 'mask', 'x', 'y', 'z']

        if 'mask' not in self.atoms:
            print('warning : mask is not defined.')
            print('warning : mask has been initialized to 0.')
            self.atoms['mask'] = np.zeros(self.get_total_atoms(),
                                          dtype=int)
        if 'vx' in self.atoms and 'vy' in self.atoms and 'vz' in self.atoms:
            out_columns.append('vx')
            out_columns.append('vy')
            out_columns.append('vz')

        # 1-indexed
        self.atoms.index = self.atoms.index + 1
        self.atoms.to_csv(ofn, columns=out_columns, mode='a', header=False,
                          sep='\t', float_format='%.6f')
        # 0-indexed
        self.atoms.index = self.atoms.index - 1

