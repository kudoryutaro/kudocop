#!/usr/bin/env python3

import numpy as np

from . import output_abstract
from .. import simulationdat as sdat

class OutputLaich(output_abstract.OutputFile):
    def __init__(self, ofn : str):
        super().__init__()
        self.ofn = ofn

    def output_file(self, data : sdat.SimulationDat(), out_columns, time_step):

        header_line = []
        header_line.append(f"#cellx {0.0}  {data.newcell[0]}")
        header_line.append(f"#celly {0.0}  {data.newcell[1]}")
        header_line.append(f"#cellz {0.0}  {data.newcell[2]}")
        header_line.append("")

        #if not len(data.type_set):
        data.type_set = set(data.particles['type'])

        header_line.append(f"#masses {len(data.type_set)}")
        for _type in data.type_set:
            header_line.append(f"{_type} {data.type_to_mass[_type]}")

        header_line.append("")

        for fix in data.fix_info:
            header_line.append(fix)

        for move in data.move_info:
            header_line.append(move)

        for press in data.press_info:
            header_line.append(press)

        for sumforce in data.sumforce_info:
            header_line.append(sumforce)

        for thermofree in data.thermofree_info:
            header_line.append(thermofree)

        for wall in data.wall_info:
            header_line.append(wall)

        header_line.append("")


        header_line.append(f"#atoms {data.total_particle}")

        body = [''] * data.total_particle
        out_columns = ['id', 'type', 'mask', 'pos']

        if 'mask' not in data.particles:
            data.particles['mask'] = np.zeros(data.total_particle,
                                              dtype=int)
        if 'velo' in data.particles:
            out_columns.append('velo')

        output_list = [ np.round(data.particles[colum], decimals=6).astype(str) if data.particles[colum].ndim > 1 else
                       np.round(data.particles[colum], decimals=6).reshape(data.total_particle,1).astype(str) for colum in out_columns ]

        body_data = np.hstack(output_list)

        np.savetxt(self.ofn, body_data, fmt="%s", header="\n".join(header_line), comments='', delimiter='    ')


        if data.flagconnect:
            with open(self.ofn, 'a') as ofs:
                print(f"#connect {len(data.connect_list)}", file=ofs)
                for ind, connect in enumerate(data.connect_list):
                    tmp_list = [str(ind+1), str(len(connect))]
                    for target in connect:
                        tmp_list.append(str(target+1))
                    print(" ".join(tmp_list), file=ofs)


