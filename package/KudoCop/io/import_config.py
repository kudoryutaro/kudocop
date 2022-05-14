import numpy as np


class ImportConfig():
    def __init__(self):
        pass

    def import_config(self, ifn: str) -> None:
        with open(ifn, 'r') as ifp:
            lines = ifp.readlines()

        for line in lines:
            spline = line.split()
            if len(spline) == 0:
                continue
            target = spline[0]
            value = spline[1]
            if target == "MPIGridX":
                self.mpigrid[0] = int(value)
            elif target == "MPIGridY":
                self.mpigrid[1] = int(value)
            elif target == "MPIGridZ":
                self.mpigrid[2] = int(value)
            elif target == "OMPGridX":
                self.ompgrid[0] = int(value)
            elif target == "OMPGridY":
                self.ompgrid[1] = int(value)
            elif target == "OMPGridZ":
                self.ompgrid[2] = int(value)
            elif target == "CUTOFF":
                self.cutoff = float(value)
            elif target == "MARGIN":
                self.margin = float(value)
            elif target == "DecompX":
                self.flagdecomp = True
                self.decompx = [float(_) for _ in spline[1:]]
            elif target == "DecompY":
                self.flagdecomp = True
                self.decompy = [float(_) for _ in spline[1:]]
            elif target == "DecompZ":
                self.flagdecomp = True
                self.decompz = [float(_) for _ in spline[1:]]
            elif target == "TotalStep":
                self.total_step = int(value)
            elif target == "FileStep":
                self.file_step = int(value)

