#!/usr/bin/env python3

import numpy as np

from . import import_abstract
from .. import simulationdat as sdat

class ImportLaich(import_abstract.ImportFile):
    def __init__(self, ifn : str):
        super().__init__()
        self.ifn = ifn

    def import_file(self, data : sdat.SimulationDat):
        with open(self.ifn, 'r') as ifp:
            lines = ifp.readlines()

        for ind, line in enumerate(lines):
            spline = line.split()
            if len(spline) == 0:
                continue
            if spline[0] == "#cellx":
                data.dcell[0] = float(spline[1])
                data.cell[0] = float(spline[2])
                data.newcell[0] = float(spline[2]) - float(spline[1])

            if spline[0] == "#celly":
                data.dcell[1] = float(spline[1])
                data.cell[1] = float(spline[2])
                data.newcell[1] = float(spline[2]) - float(spline[1])

            if spline[0] == "#cellz":
                data.dcell[2] = float(spline[1])
                data.cell[2] = float(spline[2])
                data.newcell[2] = float(spline[2]) - float(spline[1])

            if spline[0] == "#masses":
                elem_num = int(spline[1])
                for _line in lines[ind+1:ind+1+elem_num]:
                    _spline = _line.split()
                    data.type_to_mass[int(_spline[0])] = float(_spline[1])
                    data.type_set.add(int(_spline[0]))

            if spline[0] == "#fix":
                data.fix_info.append(line.replace("\n",""))

            if spline[0] == "#move":
                data.move_info.append(line.replace("\n", ""))

            if spline[0] == "#press":
                data.press_info.append(line.replace("\n",""))

            if spline[0] == "#thermostatfree":
                data.thermofree_info.append(line.replace("\n",""))

            if spline[0] == "#wall":
                data.wall_info.append(line.replace("\n",""))

            if spline[0] == "#sumforce":
                data.sumforce_info.append(line.replace("\n", ""))
                sumforce_num = int(spline[1])
                for _line in lines[ind+1:ind+1+sumforce_num]:
                    data.sumforce_info.append(_line.replace("\n", ""))


            if spline[0] == "#atoms":
                data.total_particle = int(spline[1])
                splines = np.array( [ l.split() for l in lines[ind+1:ind+1+int(spline[1])] ] )

                part_id = splines[:, 0:1].astype(int)
                part_type = splines[:, 1:2].astype(int)
                part_mask = splines[:, 2:3].astype(int)
                part_pos = splines[:, 3:6].astype(float)

                data.particles['id'] = part_id[:, 0]
                data.particles['type'] = part_type[:, 0]
                data.particles['mask'] = part_mask[:, 0]
                data.particles['pos'] = part_pos[:, :]

                if np.shape(splines)[1] > 6:
                    part_velocity = splines[:, 6:9].astype(float)
                    data.particles['velo'] = part_velocity[:, :]

            if spline[0] == "#connect":
                read_l =  int(spline[1])
                splines = [l.split() for l in lines[ind+1:ind+1+int(spline[1])] ]
                data.connect_list = [ [int(v) - 1 for v in l[2:]] for l in splines ]
