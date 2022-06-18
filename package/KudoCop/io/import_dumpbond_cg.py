import numpy as np
import pandas as pd


class ImportDumpbondCG():
    def __init__(self):
        pass

    def import_dumpbond_cg(self, ifn: str) -> None:
        with open(ifn, 'r') as ifp:
            lines = ifp.readlines()
        total_atoms = len(lines) - 1
        self.connect_list_from_dumpbond_gc = [[] for _ in range(total_atoms)]
        for line in lines[1:]:
            spline = list(map(int, line.split()))
            atom_idx = spline[0] - 1
            c_list = spline[2:]
            for i in range(len(c_list)):
                c_list[i] -= 1

            self.connect_list_from_dumpbond_gc[atom_idx] = c_list
