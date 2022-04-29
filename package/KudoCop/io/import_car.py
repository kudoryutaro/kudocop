#!/usr/bin/env python3

import os

import numpy as np

from . import import_abstract
from .. import simulationdat as sdat
from ..atomic_info.atomic_weight import get_weight

class ImportCar(import_abstract.ImportFile):
    def __init__(self, ifn : str):
        super().__init__()
        self.ifn = ifn

    def import_file(self, data : sdat):
        with open(self.ifn, 'r') as ifp:
            lines = ifp.readlines()

        body_start = -1
        body_end = -1

        for ind, line in enumerate(lines):
            spline = line.split()

            if not spline:
                pass
            elif spline[0] == "end":
                if lines[ind+1].split()[0] == "end":
                    body_end = ind + 1
                    break
            elif spline[0] == "PBC":
                data.cell[:] = list(map(float, spline[1:4]))
                data.newcell[:] = list(map(float, spline[1:4]))
                body_start = ind + 1

        splines = np.array( [l.split() for l in lines[body_start:body_end] if len(l) > 4] )

        data.total_particle = len(splines)

        part_id = np.arange(1,data.total_particle+1)
        part_element = splines[:, 7].astype(str)
        part_pos = splines[:, 1:4].astype(float)

        data.particles['id'] = part_id
        data.particles['element'] = part_element
        data.particles['pos'] = part_pos[:, :]

        if data.elem_to_type:
            part_type = np.array([data.elem_to_type[key] for key in part_element], int)
            data.particles['type'] = part_type
        else:
            print("elem_to_type is not defined")
            exit()

        data.wrap_particles()

        for elem in set(part_element):
            data.type_to_mass[data.elem_to_type[elem]] = get_weight(elem)
            data.type_set.add(data.elem_to_type[elem])

        ifn_mdf = os.path.splitext(self.ifn)[0] + ".mdf"
        if os.path.exists("./"+ifn_mdf):

            elemtag_to_id = {elemtag:ind for ind, elemtag
                             in enumerate(splines[:, 0].astype(str)) }

            with open(ifn_mdf, 'r') as ifp:
                lines = ifp.readlines()

            for ind, line in enumerate(lines):
                spline = line.split()
                if len(spline) <= 0:
                    continue
                if spline[0] == "@molecule":
                    body_start = ind + 2
                    break
            splines = [l.split() for l
                                 in lines[body_start:body_start+data.total_particle] ]
            data.connect_list = [ [ elemtag_to_id[elemtag.partition("%")[0].partition("/")[0]] for elemtag in s[12:] ]
                                    for s in splines]
