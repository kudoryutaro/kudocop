import numpy as np
import pandas as pd


class ImportDumpbondCG():
    def __init__(self):
        pass

    def import_dumpbond_cg(self, ifn: str, max_coordination_num=6) -> None:
        df = pd.read_csv(ifn, skiprows=1, sep=' ',
                         header=None, names=range(max_coordination_num))
        df.drop([0, 2], axis=1, inplace=True)
        df = df.fillna(-1).astype(int) - 1
        self.connect_list_from_dumpbond_cg = [
            None for _ in range(len(df.index))]

        values = df.values

        for row in values:
            atom_idx = row[0]
            for col_i, next_atom_idx in enumerate(row):
                if next_atom_idx < 0:
                    break
            self.connect_list_from_dumpbond_cg[atom_idx] = row[1:col_i]
