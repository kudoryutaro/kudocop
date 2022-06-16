import numpy as np
import pandas as pd
from tqdm import tqdm


class ImportEnergy():
    def __init__(self):
        pass

    def import_energy(self, ifn: str) -> None:
        cols = [
            'bond',
            'lone_pair',
            'overcoordinate',
            'undercoordinate',
            'valence_angle',
            'penalty',
            'three_body_conjugation',
            'torsion_angle',
            'four_body_conjugation',
            'hydrogen_bond',
            'coulomb',
            'van_der_waals',
            'not implemented'
        ]
        self.energy = pd.read_csv(
            ifn, skiprows=1, index_col=0, names=cols, sep=' ')
