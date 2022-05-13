import numpy as np
import pandas as pd


class ReadDumppos():
    def __init__(self):
        pass

    def read_file(self, sdat, ifn: str) -> None:
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
                        sdat.cell[dim] = float(spline[1])
                    current_row += 3
                    continue
                if spline[0] == 'ITEM:' and spline[1] == 'ATOMS':
                    columns = spline[3:]
                    break

        skip_rows = current_row
        sdat.atoms = pd.read_csv(
            ifn, skiprows=skip_rows, sep=' ', names=columns)
        sdat.atoms[['type', 'mask']] = sdat.atoms[['type', 'mask']].astype(int)
        sdat.atoms.index = sdat.atoms.index - 1
        sdat.atoms.sort_index(inplace=True)
        sdat.atom_type_set |= set(sdat.atoms['type'])
