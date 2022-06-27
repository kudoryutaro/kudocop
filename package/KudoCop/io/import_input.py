import numpy as np
import pandas as pd
import sys


class ImportInput():
    def __init__(self):
        pass

    def import_input(self, ifn: str):
        if self.atom_symbol_to_type is None:
            print('error : atom_symbol_to_type is not defined')
            print('Import para first')
            sys.exit(-1)

        with open(ifn, 'r') as ifp:
            lines = ifp.readlines()

        for idx, line in enumerate(lines):
            spline = line.split()
            if len(spline) == 0:
                continue
            if spline[0] == "#cellx":
                self.cell[0] = float(spline[2])

            if spline[0] == "#celly":
                self.cell[1] = float(spline[2])

            if spline[0] == "#cellz":
                self.cell[2] = float(spline[2])

            if spline[0] == "#masses":
                elem_num = int(spline[1])
                for _line in lines[idx+1:idx+1+elem_num]:
                    _spline = _line.split()
                    self.atom_type_to_mass[int(_spline[0])] = float(_spline[1])

            if spline[0] == "#fix":
                self.fix_info.append(line.rstrip())

            if spline[0] == "#move":
                self.move_info.append(line.rstrip())

            if spline[0] == "#press":
                self.press_info.append(line.rstrip())

            if spline[0] == "#thermofree":
                self.thermofree_info.append(line.rstrip())

            if spline[0] == "#wall":
                self.wall_info.append(line.rstrip())

            if spline[0] == "#sumforce":
                self.sumforce_info.append(line.rstrip())
                sumforce_num = int(spline[1])
                for _line in lines[idx+1:idx+1+sumforce_num]:
                    self.sumforce_info.append(_line.rstrip())

            atom_data = dict()
            if spline[0] == "#atoms":
                splines = np.array([l.split()
                                   for l in lines[idx+1:idx+1+int(spline[1])]])
                # 0-indexed
                index = splines[:, 0].astype(int) - 1
                atom_data['type'] = splines[:, 1].astype(int)
                atom_data['mask'] = splines[:, 2].astype(int)
                atom_data['x'] = splines[:, 3].astype(float)
                atom_data['y'] = splines[:, 4].astype(float)
                atom_data['z'] = splines[:, 5].astype(float)
                try:
                    atom_data['vx'] = splines[:, 6].astype(float)
                    atom_data['vy'] = splines[:, 7].astype(float)
                    atom_data['vz'] = splines[:, 8].astype(float)
                except:
                    pass

                self.atoms = pd.DataFrame(data=atom_data, index=index)

            if spline[0] == "#connect":
                read_l = int(spline[1])
                splines = [l.split()
                           for l in lines[idx+1:idx+1+int(spline[1])]]
                self.connect_list = [[int(v) - 1 for v in l[2:]]
                                     for l in splines]
