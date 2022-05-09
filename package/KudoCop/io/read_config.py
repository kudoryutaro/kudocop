import numpy as np


class ReadConfig():
    def __init__(self):
        pass

    def read_file(self, sdat, ifn: str) -> None:
        with open(ifn, 'r') as ifp:
            lines = ifp.readlines()

        for line in lines:
            spline = line.split()
            if len(spline) == 0:
                continue
            target = spline[0]
            value = spline[1]
            if target == "MPIGridX":
                sdat.mpigrid[0] = int(value)
            elif target == "MPIGridY":
                sdat.mpigrid[1] = int(value)
            elif target == "MPIGridZ":
                sdat.mpigrid[2] = int(value)
            elif target == "OMPGridX":
                sdat.ompgrid[0] = int(value)
            elif target == "OMPGridY":
                sdat.ompgrid[1] = int(value)
            elif target == "OMPGridZ":
                sdat.ompgrid[2] = int(value)
            elif target == "CUTOFF":
                sdat.cutoff = float(value)
            elif target == "MARGIN":
                sdat.margin = float(value)
            elif target == "DecompX":
                sdat.flagdecomp = True
                sdat.decompx = [float(_) for _ in spline[1:]]
            elif target == "DecompY":
                sdat.flagdecomp = True
                sdat.decompy = [float(_) for _ in spline[1:]]
            elif target == "DecompZ":
                sdat.flagdecomp = True
                sdat.decompz = [float(_) for _ in spline[1:]]
            elif target == "TotalStep":
                sdat.total_step = int(value)
            elif target == "FileStep":
                sdat.file_step = int(value)

        sdat.set_decomp()
        self.check_grid(sdat)

    def check_grid(self, sdat):
        totalgrid = 1
        for i in range(3):
            totalgrid *= sdat.mpigrid[i]
            totalgrid *= sdat.ompgrid[i]
        if totalgrid == 0:
            print("# Error : totoal grid = 0")
            exit()

        cut_off = sdat.cutoff + sdat.margin
        for i in range(3):
            grid = sdat.mpigrid[i] * sdat.ompgrid[i]
            localcell = sdat.cell[i] / grid
            if localcell < cut_off:
                print("# Error : local cell < CUTOFF+MARGIN (%d)" % i)
                exit()

        if(sdat.flagdecomp):
            for i in range(3):
                decomp = []
                if i == 0:
                    decomp = sdat.decompx
                elif i == 1:
                    decomp = sdat.decompy
                elif i == 2:
                    decomp = sdat.decompz
                totaldecomp = 0.0
                cell = sdat.cell[i]
                for ind, val in enumerate(decomp):
                    totaldecomp += val
                    localcell = val * cell
                    if localcell < cut_off:
                        print("# Error : decomp < CUTOFF+MARGIN (%d in %d)" %
                              (ind, i))
                        exit()
                if not totaldecomp == 1.0:
                    print("# Error : total decomp != 1.0 (%d)" % i)
                    exit()
