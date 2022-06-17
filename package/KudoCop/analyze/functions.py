import numpy as np
import pandas as pd


def get_target_atoms(sdat, target_atoms=None, condition=None, step_idx=None):
    if target_atoms is not None:
        return target_atoms

    if condition is None:
        return np.array([True] * sdat.get_total_atoms())

    if step_idx is None:
        return condition(sdat)
    else:
        return condition(sdat, step_idx)
