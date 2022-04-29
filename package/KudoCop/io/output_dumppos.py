#!/usr/bin/env python3

import numpy as np

from . import output_abstract
from .. import simulationdat as sdat

class OutputDumpPos(output_abstract.OutputFile):
    def __init__(self, ofn : str):
        super().__init__()
        self.ofn = ofn

    def output_file(self, data : sdat, out_columns, time_step):

        header_line = []
        header_line.append("ITEM: TIMESTEP")
        header_line.append(f"{time_step}")
        header_line.append("ITEM: NUMBER OF ATOMS")
        header_line.append(f"{data.total_particle}")
        header_line.append("ITEM: BOX BOUNDS pp pp pp")
        for dim in range(3):
            header_line.append(f"{data.dcell[dim]} {data.cell[dim]}")

        dump_columns = self.detect_clumns(out_columns,data.colum2suffix)
        header_line.append(" ".join(["ITEM: ATOMS"] + dump_columns))

        body = [''] * data.total_particle

        output_list = [ data.particles[colum].astype(str) if data.particles[colum].ndim > 1 else
                       data.particles[colum].reshape(data.total_particle,1).astype(str) for colum in out_columns ]

        body_data = np.hstack(output_list)

        np.savetxt(self.ofn, body_data, fmt="%s", header = "\n".join(header_line), comments='')

    @staticmethod
    def detect_clumns(columns_list, colum2suffix):

        colum_suffix = []

        ind = 0
        for colum in columns_list:

            if colum in colum2suffix:
                for value in colum2suffix[colum]:
                    colum_suffix.append(value)
            else:
                colum_suffix.append(colum)

        return colum_suffix
