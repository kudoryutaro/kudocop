import numpy as np
import pandas as pd


class ExportXyz():
    def __init__(self):
        pass

    def export_xyz(self, ofn: str, out_columns=None, structure_name=None) -> None:
        if out_columns is None:
            out_columns = ['type', 'x', 'y', 'z']
        if structure_name is None:
            print('warning : structure_name is not defined')
            print('warning : structure_name has been initialized to \'structure\'')
            structure_name = 'structure'

        # header_line
        header_line = []
        header_line.append(str(self.get_total_atoms()) + '\n')
        header_line.append(structure_name + '\n')

        with open(ofn, 'w') as ofp:
            ofp.writelines(header_line)
        self.atoms.to_csv(ofn, columns=out_columns, sep='\t',
                          mode='a', header=False, index=False, float_format='%.6f')
