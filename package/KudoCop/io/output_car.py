#!/usr/bin/env python3

import time

import numpy as np

from . import output_abstract
from .. import simulationdat as sdat

class OutputCar(output_abstract.OutputFile):
    def __init__(self, ofn : str):
        super().__init__()
        self.ofn = ofn

    def output_file(self, data : sdat.SimulationDat, out_columns, time_step):

        header_line = [
            "!BIOSYM archive 3\n"
            "PBC=ON\n",
            "Materials Studio Generated CAR File\n",
            "!DATE %s\n" % (time.strftime(
                '%c', time.localtime(time.time()))),
            "PBC  %8.4f  %8.4f  %8.4f   90.0000   90.0000   90.0000 (P1)\n"
                % (data.newcell[0], data.newcell[1], data.newcell[2])
        ]

        fmt = "{ec:5s}  {p[0]:13.9f}  {p[1]:13.9f}  {p[2]:13.9f} XXXX 1      xx      {el:3s} 0.000\n"

        _type = data.particles["type"]
        elem = np.array([data.type_to_elem[key] for key in _type])
        pos = data.particles["pos"][:,:]

        dictionary = {elem:0 for elem in set(elem)}

        elem_cnt = []
        for el in elem:
            dictionary[el] += 1
            elem_cnt.append(el+str(dictionary[el]))

        body_line = [fmt.format(ec=ec, el=e, p=p)
                     for ec, e, p in zip(elem_cnt, elem, pos)]

        footer_line = ['end\n', 'end\n']

        with open(self.ofn, 'w') as ofp:
            ofp.write(''.join(header_line) + ''.join(body_line) + ''.join(footer_line))


