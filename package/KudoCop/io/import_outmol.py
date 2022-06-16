import numpy as np
import pandas as pd
import sys


class ImportOutmol():
    def __init__(self):
        pass

    def import_outmol(self, ifn: str):
        '''
        構造最適化によって最適化された最後の構造を読み込む
        '''
        with open(ifn, 'r') as ifp:
            lines = ifp.readlines()

        for index, line in enumerate(lines):
            spline = line.split()
            if len(spline) == 0:
                continue
            if len(spline) >= 2 and spline[0] == 'Final' and spline[1] == 'Coordinates':
                start_index = index + 3
            if len(spline) >= 2 and spline[0] == '$cell' and spline[1] == 'vectors':
                self.cell[0] = round(float(lines[index + 1].split()
                                     [0]) / 1.889726124993590, 2)
                self.cell[1] = round(float(lines[index + 2].split()
                                     [1]) / 1.889726124993590, 2)
                self.cell[2] = round(float(lines[index + 3].split()
                                     [2]) / 1.889726124993590, 2)

        index = start_index
        data = []

        while True:
            spline = lines[index].split()
            try:
                int(spline[0])
                data.append(spline)
            except:
                break
            index += 1

        df_outmol = pd.DataFrame(data, columns=['id', 'type', 'x', 'y', 'z'])
        df_outmol.drop('id', axis=1, inplace=True)
        df_outmol[['x', 'y', 'z']] = df_outmol[['x', 'y', 'z']].astype(float)

        df_outmol['type'] = df_outmol['type'].replace(self.atom_symbol_to_type)

        self.atoms = df_outmol.copy()
