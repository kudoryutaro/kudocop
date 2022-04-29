#!/usr/bin/env python3

import numpy as np

from . import import_abstract

class ImportDumpBond(import_abstract.ImportFile):
    def __init__(self, ifn : str):
        super().__init__()
        self.ifn = ifn

    def import_file(self, data):
        with open(self.ifn, 'r') as ifp:
            lines = ifp.readlines()

        read_flag = False
        for ind, line in enumerate(lines):
            spline = line.split()
            if len(spline) == 0:
                continue

            if spline[0] == "ITEM:" and spline[1] == "NUMBER":
                read_particle_num = int(lines[ind+1])
                data.bondorder_list = [ [] for _ in range(read_particle_num) ]

            elif spline[0] == "Atom":
                atom_id = int(spline[1])
                read_flag = True

            elif read_flag:
                data.bondorder_list[atom_id-1].append([int(s) - 1 if s.isdigit() else float(s) for s in spline])

