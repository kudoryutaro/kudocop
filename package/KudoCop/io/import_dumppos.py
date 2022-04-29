#!/usr/bin/env python3

import numpy as np

from . import import_abstract

class ImportDumpPos(import_abstract.ImportFile):
    def __init__(self, ifn : str):
        super().__init__()
        self.ifn = ifn

    def import_file(self, data):
        with open(self.ifn, 'r') as ifp:
            lines = ifp.readlines()

        for ind, line in enumerate(lines):
            spline = line.split()
            if len(spline) == 0:
                continue
            if spline[0] == "ITEM:" and spline[1] == "NUMBER":
                data.total_particle = int(lines[ind+1])

            elif spline[0] == "ITEM:" and spline[1] == "BOX":
                for dim in range(3):
                    input_box = lines[ind+dim+1].split()
                    data.dcell[dim] = float(input_box[0])
                    data.cell[dim] = float(input_box[1])
                    data.newcell[dim] = float(input_box[1]) - float(input_box[0])

            elif spline[0] == "ITEM:" and spline[1] == "ATOMS":

                columns = []
                colum_suffix_num = self.detect_clumns(spline[2:])
                splines = np.array([ l.split() for l in lines[ind+1:ind+1+data.total_particle] ])

                for colum, suffix_id, num in colum_suffix_num:
                    _dtype = self.detect_dtype(splines[0, suffix_id])
                    partial = splines[:, suffix_id:suffix_id+num].astype(_dtype)
                    if num > 1:
                        data.particles[colum] = partial[:,:]
                    else:
                        data.particles[colum] = partial[:,0]
                    columns.append(colum)
                data.sort_particles_by_id(columns)


    @staticmethod
    def detect_clumns(columns_list):

        colum_suffix = []
        colum2suffix = {'x':('pos', 3), 'vx':('velo', 3), 'fx':('force', 3)}

        ind = 0
        while ind < len(columns_list):
            colum = columns_list[ind]
            if colum in colum2suffix:
                values = colum2suffix[colum]
                colum_suffix.append([values[0], ind, values[1]])
                ind = ind + values[1]
            else:
                colum_suffix.append([colum, ind, 1])
                ind = ind + 1
        return colum_suffix

    @staticmethod
    def detect_dtype(da: str):
        try:
            int(da)
            return int
        except ValueError:
            try:
                float(da)
                return float
            except ValueError:
                return str
