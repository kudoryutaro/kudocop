import numpy as np
import pandas as pd


class ExportLammpsInput():
    def __init__(self):
        pass

    def export_lammps_input(self, ofn: str) -> None:
        for dim in range(3):
            if self.cell[dim] == 0:
                print(f'warning : cell[{dim}] is 0')

            if self.cell[dim] is None:
                print(f'warning : cell[{dim}] is not defined')
                print(f'warning : cell[{dim}] has been initialized to 0')
                self.cell[dim] = 0.0

        # header_line
        header_line = []
        header_line.append(f"# LAMMPS data\n")
        header_line.append(f"{self.get_total_atoms()} atoms\n")
        header_line.append(f"{max([atom_type for atom_type in self.atom_type_to_mass.keys()])} atom types\n")
        header_line.append(f"{0.0} {self.cell[0]} xlo xhi\n")
        header_line.append(f"{0.0} {self.cell[1]} ylo yhi\n")
        header_line.append(f"{0.0} {self.cell[2]} zlo zhi\n")
        header_line.append(f"0.0 0.0 0.0 xy xz yz\n")
        header_line.append(f"\n")
        header_line.append(f"Atoms\n")
        header_line.append(f"\n")
        

        with open(ofn, 'w') as ofs:
            ofs.writelines(header_line)

        out_columns = ['type', 'x', 'y', 'z']

        # if 'mask' not in self.atoms:
        #     print('warning : mask is not defined.')
        #     print('warning : mask has been initialized to 0.')
        #     self.atoms['mask'] = np.zeros(self.get_total_atoms(),
                                        #   dtype=int)
        # if 'vx' in self.atoms and 'vy' in self.atoms and 'vz' in self.atoms:
        #     out_columns.append('vx')
        #     out_columns.append('vy')
        #     out_columns.append('vz')
        self.atoms[['x', 'y', 'z']] = self.atoms[['x', 'y', 'z']].round(6)

        self.atoms['x'] = np.where( 
            self.atoms['x'] == 0, 0.001, self.atoms['x'])
        self.atoms['y'] = np.where(
            self.atoms['y'] == 0, 0.001, self.atoms['y'])
        self.atoms['z'] = np.where(
            self.atoms['z'] == 0, 0.001, self.atoms['z'])

        # 1-indexed
        self.atoms.index = self.atoms.index + 1

        self.atoms[out_columns].to_csv(
            ofn,
            sep=' ',
            mode='a',
            header=None,
        )

        # 0-indexed
        self.atoms.index = self.atoms.index - 1
