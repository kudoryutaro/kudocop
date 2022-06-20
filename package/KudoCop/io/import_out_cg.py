import numpy as np
import pandas as pd
import sys


class ImportOutCg():
    def __init__(self):
        pass

    def import_out_cg(self, ifn: str):
        current_row = 0
        with open(ifn, 'r') as ifp:
            while True:
                current_row += 1
                spline = ifp.readline().split()
                if len(spline) == 0:
                    continue
                if len(spline) >= 2 and spline[0] == '#' and spline[1] == 'STEP':
                    columns = spline[1:]
                    skiprows = current_row
                    break
        columns.append('sec')
        df_out_cg = pd.read_csv(
            ifn, skiprows=skiprows, sep='\s+', header=None, names=columns, index_col=False, low_memory=False)

        df_out_cg.drop('sec', axis=1, inplace=True)

        break_points = df_out_cg.loc[df_out_cg['STEP'] == 'break', [
            'TEMP', 'KineticE', 'PotentialE']].values.tolist()

        self.break_points = [None for _ in range(len(break_points))]

        for i in range(len(self.break_points)):
            self.break_points[i] = [int(break_points[i][0]), int(
                break_points[i][1]), float(break_points[i][2])]
        df_out_cg = df_out_cg.loc[~(
            (df_out_cg['STEP'] == 'break') | (df_out_cg['STEP'] == '#')), :]

        df_out_cg['STEP'] = df_out_cg['STEP'].astype(int)
        df_out_cg[['TEMP', 'KineticE', 'PotentialE', 'TotalE', 'WTime']] = df_out_cg[[
            'TEMP', 'KineticE', 'PotentialE', 'TotalE', 'WTime']].astype(float)
        # df_out_cg.reset_index(inplace=True)
        df_out_cg.set_index('STEP', inplace=True)
        self.out_cg = df_out_cg
