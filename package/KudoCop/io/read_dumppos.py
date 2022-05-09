import numpy as np
import pandas as pd


class ReadDumppos():
    def __init__(self):
        pass

    def read_file(self, sdat, ifn: str) -> None:
        with open(ifn, 'r') as ifp:
            lines = ifp.readlines()

        for ind, line in enumerate(lines):
            spline = line.split()
            if len(spline) == 0:
                continue
            if spline[0] == "ITEM:" and spline[1] == "NUMBER":
                total_particle = int(lines[ind+1])

            elif spline[0] == "ITEM:" and spline[1] == "BOX":
                for dim in range(3):
                    input_box = lines[ind+dim+1].split()
                    sdat.cell[dim] = float(input_box[1])

            elif spline[0] == "ITEM:" and spline[1] == "ATOMS":
                spline = spline[2:]
                atom_data = dict()
                splines = np.array([l.split()
                                   for l in lines[ind+1:ind+1+total_particle]])
                for c in range(len(spline)):
                    if spline[c] == 'id':
                        # 0-indexed
                        index = splines[:, c].astype(int) - 1
                    elif spline[c] in ['type', 'mask']:
                        atom_data[spline[c]] = splines[:, c].astype(int)
                    else:
                        atom_data[spline[c]] = splines[:, c].astype(float)

                sdat.atoms = pd.DataFrame(data=atom_data, index=index)
                sdat.atoms.sort_index(inplace=True)

        sdat.atom_type_set |= set(sdat.atoms['type'])
