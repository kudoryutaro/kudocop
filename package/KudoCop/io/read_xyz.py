import numpy as np
import pandas as pd
import re
import sys


class ReadXyz():
    def __init__(self):
        pass

    def read_file(self, sdat, ifn: str) -> None:
        with open(ifn, 'r') as ifp:
            lines = ifp.readlines()

        total_particle = int(lines[0])
        lattice_value = re.search('Lattice="(.*?)"', lines[1])

        if lattice_value is not None:
            cellsize = lattice_value.group(1).split()
            cellx = float(cellsize[0])
            celly = float(cellsize[4])
            cellz = float(cellsize[8])
            sdat.cell = [cellx, celly, cellz]

        splines = np.array([l.split() for l in lines[2:2+total_particle]])

        atom_data = dict()
        try:
            atom_data['type'] = splines[:, 0].astype(int)
        except:
            atom_symbols = splines[:, 0].astype(str)
            if sdat.atom_symbol_to_type is None:
                print('error : atom_tymbol_to_type is not defined')
                print('Read para first')
                sys.exit(-1)
            atom_data['type'] = np.array(
                [sdat.atom_symbol_to_type[atom_symbol] for atom_symbol in atom_symbols])
        # 0-indexed
        index = np.arange(total_particle)
        atom_data['x'] = splines[:, 1].astype(float)
        atom_data['y'] = splines[:, 2].astype(float)
        atom_data['z'] = splines[:, 3].astype(float)

        sdat.atoms = pd.DataFrame(data=atom_data, index=index)
        sdat.atoms.sort_index(inplace=True)
        sdat.atom_type_set |= set(sdat.atoms['type'])
