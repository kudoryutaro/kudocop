import numpy as np
import pandas as pd


class ImportDumppos():
    def __init__(self):
        pass

    def import_dumppos(self, ifn: str) -> None:
        current_row = 0
        with open(ifn, 'r') as ifp:
            while True:
                current_row += 1
                spline = ifp.readline().split()
                if len(spline) == 0:
                    continue
                if spline[0] == "ITEM:" and spline[1] == "BOX":
                    for dim in range(3):
                        spline = ifp.readline().split()
                        self.cell[dim] = float(spline[1])
                    current_row += 3
                    continue
                if spline[0] == 'ITEM:' and spline[1] == 'ATOMS':
                    columns = spline[3:]
                    break

        skip_rows = current_row
        self.atoms = pd.read_csv(
            ifn, skiprows=skip_rows, sep=' ', names=columns)
        self.atoms[['type', 'mask']] = self.atoms[['type', 'mask']].astype(int)
        self.atoms.index = self.atoms.index - 1
        self.atoms.sort_index(inplace=True)
        self.atom_type_set |= set(self.atoms['type'])
