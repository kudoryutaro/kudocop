#!/usr/bin/env python3

import os
import re

import numpy as np

from . import import_abstract

class ImportXYZ(import_abstract.ImportFile):
    def __init__(self, ifn : str):
        super().__init__()
        self.ifn = ifn

    def import_file(self, data):

        with open(self.ifn, 'r') as ifp:
            lines = ifp.readlines()

        data.total_particle = int(lines[0])

        lattice_value = re.search('Lattice="(.*?)"', lines[1])
        if lattice_value is not None:
            cellsize = lattice_value.group(1).split()
            cellx = float(cellsize[0])
            celly = float(cellsize[4])
            cellz = float(cellsize[8])
            data.set_cellsize([[0,cellx], [0,celly], [0,cellz]])

        splines = np.array( [ l.split() for l in lines[2:2+data.total_particle] ] )

        _dtype = self.detect_dtype(splines[0, 0])
        part_elem = splines[:, 0].astype(_dtype)
        if _dtype == int:
            data.particles['type'] = part_elem

        else:
            data.particles['element'] = part_elem
            if data.elem_to_type:
                part_type = np.array([data.elem_to_type[key] for key in part_elem], int)
                data.particles['type'] = part_type
            else:
                print("elem_to_type is not defined")
                exit()

        part_pos = splines[:, 1:4].astype(float)

        data.particles['id'] = np.arange(1,data.total_particle+1)
        data.particles['pos'] = part_pos

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
