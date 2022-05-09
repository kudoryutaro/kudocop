import numpy as np
import pandas as pd


class ReadDumpbond():
    def __init__(self):
        pass

    def read_file(self, sdat, ifn: str) -> None:
        with open(ifn, 'r') as ifp:
            lines = ifp.readlines()

        read_flag = False
        for ind, line in enumerate(lines):
            spline = line.split()
            if len(spline) == 0:
                continue

            if spline[0] == "ITEM:" and spline[1] == "NUMBER":
                read_particle_num = int(lines[ind+1])
                sdat.bondorder_list = [[] for _ in range(read_particle_num)]

            elif spline[0] == "Atom":
                # 0-indexed
                atom_idx = int(spline[1]) - 1
                read_flag = True

            elif read_flag:
                sdat.bondorder_list[atom_idx].append(
                    [int(s) - 1 if s.isdigit() else float(s) for s in spline])
        