#!/usr/bin/env python3

import numpy as np

from . import import_abstract

class ImportConfig(import_abstract.ImportFile):
    def __init__(self, ifn : str):
        super().__init__()
        self.ifn = ifn

    def import_file(self, data):
        with open(self.ifn, 'r') as ifp:
            lines = ifp.readlines()

        for line in lines:
            spline = line.split()
            if len(spline) == 0:
                continue
            target = spline[0]
            value = spline[1]
            if target == "MPIGridX":
                data.mpigrid[0] = int(value)
            elif target == "MPIGridY":
                data.mpigrid[1] = int(value)
            elif target == "MPIGridZ":
                data.mpigrid[2] = int(value)
            elif target == "OMPGridX":
                data.ompgrid[0] = int(value)
            elif target == "OMPGridY":
                data.ompgrid[1] = int(value)
            elif target == "OMPGridZ":
                data.ompgrid[2] = int(value)
            elif target == "CUTOFF":
                data.cutoff = float(value)
            elif target == "MARGIN":
                data.margin = float(value)
            elif target == "DecompX":
                data.flagdecomp = True
                data.decompx = [float(_) for _ in spline[1:]]
            elif target == "DecompY":
                data.flagdecomp = True
                data.decompy = [float(_) for _ in spline[1:]]
            elif target == "DecompZ":
                data.flagdecomp = True
                data.decompz = [float(_) for _ in spline[1:]]
            elif target == "TotalStep":
                data.total_step = int(value)
            elif target == "FileStep":
                data.file_step = int(value)

        if (np.all(np.array(data.newcell) > 0.0 )):
            data.set_decomp()
            self.check_grid(data)

    def check_grid(self, data):
        totalgrid = 1
        for i in range(3):
            totalgrid *= data.mpigrid[i]
            totalgrid *= data.ompgrid[i]
        if totalgrid == 0:
            print("# Error : totoal grid = 0")
            exit()

        cut_off = data.cutoff + data.margin
        for i in range(3):
            grid = data.mpigrid[i] * data.ompgrid[i]
            localcell = data.newcell[i] / grid
            if localcell < cut_off:
                print("# Error : local cell < CUTOFF+MARGIN (%d)" % i)
                exit()

        if(data.flagdecomp):
            for i in range(3):
                decomp = []
                if i == 0:
                    decomp = data.decompx
                elif i == 1:
                    decomp = data.decompy
                elif i == 2:
                    decomp = data.decompz
                totaldecomp = 0.0
                cell = data.newcell[i]
                for ind, val in enumerate(decomp):
                    totaldecomp += val
                    localcell = val * cell
                    if localcell < cut_off:
                        print("# Error : decomp < CUTOFF+MARGIN (%d in %d)" % (ind ,i))
                        exit()
                if not totaldecomp == 1.0:
                    print("# Error : total decomp != 1.0 (%d)" % i)
                    exit()
