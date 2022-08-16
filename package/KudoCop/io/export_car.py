import numpy as np
import pandas as pd
import time
import csv


class ExportCar():
    def __init__(self):
        pass

    def export_car(self, ofn: str) -> None:
        # header_line = [
        #     "!BIOSYM archive 3\n"
        #     "PBC=ON\n",
        #     "Materials Studio Generated CAR File\n",
        #     "!DATE %s\n" % (time.strftime(
        #         '%c', time.localtime(time.time()))),
        #     "PBC  %8.4f  %8.4f  %8.4f   90.0000   90.0000   90.0000 (P1)\n"
        #     % (0, self.cell[1], self.cell[2])
        # ]
        header_line = [
            '!BIOSYM archive 3\n',
        'PBC=OFF\n'
        'Materials Studio Generated CAR File\n'
        '!DATE Thu Aug  4 15:49:49 2022\n'
        ]
        fmt = "{ec:5s}  {p[0]:13.9f}  {p[1]:13.9f}  {p[2]:13.9f} XXXX 1      xx      {el:3s} 0.000\n"

        _type = self.atoms["type"].values
        elem = np.array([self.atom_type_to_symbol[key] for key in _type])
        pos = self.atoms[['x','y','z']].values

        dictionary = {elem:0 for elem in set(elem)}

        elem_cnt = []
        for el in elem:
            dictionary[el] += 1
            elem_cnt.append(el+str(dictionary[el]))

        body_line = [fmt.format(ec=ec, el=e, p=p)
                     for ec, e, p in zip(elem_cnt, elem, pos)]

        footer_line = ['end\n', 'end\n']

        with open(ofn, 'w') as ofp:
            ofp.write(''.join(header_line) + ''.join(body_line) + ''.join(footer_line))