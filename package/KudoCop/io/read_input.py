import numpy as np
import pandas as pd


class ReadInput():
    def __init__(self):
        pass

    def read_file(self, sdat, ifn: str):
        with open(ifn, 'r') as ifp:
            lines = ifp.readlines()

        for idx, line in enumerate(lines):
            spline = line.split()
            if len(spline) == 0:
                continue
            if spline[0] == "#cellx":
                sdat.cell[0] = float(spline[2])

            if spline[0] == "#celly":
                sdat.cell[1] = float(spline[2])

            if spline[0] == "#cellz":
                sdat.cell[2] = float(spline[2])

            if spline[0] == "#masses":
                elem_num = int(spline[1])
                for _line in lines[idx+1:idx+1+elem_num]:
                    _spline = _line.split()
                    sdat.atom_type_to_mass[int(_spline[0])] = float(_spline[1])
                    sdat.atom_type_set.add(int(_spline[0]))

            if spline[0] == "#fix":
                sdat.fix_info.append(line.rstrip())

            if spline[0] == "#move":
                sdat.move_info.append(line.rstrip())

            if spline[0] == "#press":
                sdat.press_info.append(line.rstrip())

            if spline[0] == "#thermostatfree":
                sdat.thermofree_info.append(line.rstrip())

            if spline[0] == "#wall":
                sdat.wall_info.append(line.rstrip())

            if spline[0] == "#sumforce":
                sdat.sumforce_info.append(line.rstrip())
                sumforce_num = int(spline[1])
                for _line in lines[idx+1:idx+1+sumforce_num]:
                    sdat.sumforce_info.append(_line.rstrip())

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

                sdat.atoms = pd.DataFrame(data=atom_data, index=index)

            if spline[0] == "#connect":
                read_l = int(spline[1])
                splines = [l.split()
                           for l in lines[idx+1:idx+1+int(spline[1])]]
                sdat.connect_list = [[int(v) - 1 for v in l[2:]]
                                     for l in splines]
