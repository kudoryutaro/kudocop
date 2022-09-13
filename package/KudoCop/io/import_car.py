import numpy as np
import pandas as pd
import time


class ImportCar():
    def __init__(self):
        pass

    def import_car(self, ifn: str) -> None:
        df = pd.read_csv(ifn, names=['symbol+id', 'x', 'y', 'z', 'XXXX', '1', 'xx', 'symbol', '0.000'],
                         usecols=['x', 'y', 'z', 'symbol'],
                         skiprows=4,  sep='\s+')
        df['type'] = df['symbol'].map(self.atom_symbol_to_type)
        df.drop('symbol', axis=1, inplace=True)
        df.dropna(inplace=True)
        df = df[['type', 'x', 'y', 'z']]
        df['type'] = df['type'].astype(int)
        self.atoms = df
