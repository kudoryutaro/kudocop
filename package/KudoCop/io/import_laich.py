#!/usr/bin/env python3

import numpy as np

from . import import_abstract

class ImportLaich(import_abstract.ImportFile):
    def __init__(self, ifn : str):
        super().__init__()
        self.ifn = ifn

    def import_file(self,sdat):
        with open(self.ifn, 'r') as ifp:
            lines = ifp.readlines()

        for ind, line in enumerate(lines):
            spline = line.split()
            if len(spline) == 0:
                continue
            if spline[0] == "#cellx":
                sdat.dcell[0] = float(spline[1])
                sdat.cell[0] = float(spline[2])
                sdat.newcell[0] = float(spline[2]) - float(spline[1])

            if spline[0] == "#celly":
                sdat.dcell[1] = float(spline[1])
                sdat.cell[1] = float(spline[2])
                sdat.newcell[1] = float(spline[2]) - float(spline[1])

            if spline[0] == "#cellz":
                sdat.dcell[2] = float(spline[1])
                sdat.cell[2] = float(spline[2])
                sdat.newcell[2] = float(spline[2]) - float(spline[1])

            if spline[0] == "#masses":
                elem_num = int(spline[1])
                for _line in lines[ind+1:ind+1+elem_num]:
                    _spline = _line.split()
                    sdat.type_to_mass[int(_spline[0])] = float(_spline[1])
                    sdat.type_set.add(int(_spline[0]))

            if spline[0] == "#fix":
                sdat.fix_info.append(line.replace("\n",""))

            if spline[0] == "#move":
                sdat.move_info.append(line.replace("\n", ""))

            if spline[0] == "#press":
                sdat.press_info.append(line.replace("\n",""))

            if spline[0] == "#thermostatfree":
                sdat.thermofree_info.append(line.replace("\n",""))

            if spline[0] == "#wall":
                sdat.wall_info.append(line.replace("\n",""))

            if spline[0] == "#sumforce":
                sdat.sumforce_info.append(line.replace("\n", ""))
                sumforce_num = int(spline[1])
                for _line in lines[ind+1:ind+1+sumforce_num]:
                    sdat.sumforce_info.append(_line.replace("\n", ""))


            if spline[0] == "#atoms":
                sdat.total_particle = int(spline[1])
                splines = np.array( [ l.split() for l in lines[ind+1:ind+1+int(spline[1])] ] )

                part_id = splines[:, 0:1].astype(int)
                part_type = splines[:, 1:2].astype(int)
                part_mask = splines[:, 2:3].astype(int)
                part_pos = splines[:, 3:6].astype(float)

                sdat.particles['id'] = part_id[:, 0]
                sdat.particles['type'] = part_type[:, 0]
                sdat.particles['mask'] = part_mask[:, 0]
                sdat.particles['pos'] = part_pos[:, :]

                if np.shape(splines)[1] > 6:
                    part_velocity = splines[:, 6:9].astype(float)
                    sdat.particles['velo'] = part_velocity[:, :]

            if spline[0] == "#connect":
                read_l =  int(spline[1])
                splines = [l.split() for l in lines[ind+1:ind+1+int(spline[1])] ]
                sdat.connect_list = [ [int(v) - 1 for v in l[2:]] for l in splines ]
