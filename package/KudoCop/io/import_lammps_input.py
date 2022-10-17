import numpy as np
import pandas as pd
import sys


class ImportLammpsInput():
    def __init__(self):
        pass

    def import_lammps_input(self, ifn: str):

        with open(ifn, 'r') as ifp:
            lines = ifp.readlines()

        for idx, line in enumerate(lines):
            spline = line.split()
            if len(spline) == 0:
                continue

            if line.startswith('#'):
                continue
            
            # cell - x direction
            if len(spline) >= 4 and spline[2] == 'xlo' and spline[3] == 'xhi':
                self.cell[0] = float(spline[1])
            
            # cell - y direction
            if len(spline) >= 4 and spline[2] == 'ylo' and spline[3] == 'yhi':
                self.cell[1] = float(spline[1])
            
            # cell - z direction
            if len(spline) >= 4 and spline[2] == 'zlo' and spline[3] == 'zhi':
                self.cell[2] = float(spline[1])
            
            if spline[0] == 'Atoms':
                self.atoms = pd.read_csv(
                    ifn,
                    skiprows=idx + 1,
                    names=['type', 'x', 'y', 'z'],
                    dtype={'type':'int', 'x':'float', 'y':'float', 'z':'float'},
                    index_col=0,
                    sep='\s+',
                )
                self.atoms.index -= 1

